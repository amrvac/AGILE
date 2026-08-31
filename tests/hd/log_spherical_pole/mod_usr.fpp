!> Polar-axis tests on a spherical mesh: three runs from one build.
!>
!> Every test directory rebuilds AGILE from scratch, so the cost of the test
!> suite is dominated by compilation rather than by the runs. Cases that agree
!> on the *compile-time* parameters - `phys`, `geometry`, `specialboundary` and
!> `refine_usr`, the ones make/config_schema.toml turns into fypp defines - can
!> share a single build and differ only in their par file. All the spherical
!> hd pole cases do agree on those, so they live here together and pick their
!> initial condition through the `setup` entry of the &usr_list namelist:
!>
!>   uflow.par       setup = 'uniform'  a uniform Cartesian flow, single level
!>   uflow_amr.par   setup = 'uniform'  the same, with half the domain in phi
!>                                      refined so blocks meet across the axis
!>                                      at different levels
!>   blast_amr.par   setup = 'blast'    a blast wave that starts off the axis
!>                                      and expands across it, Lohner-refined
!>
!> Both setups are the same uniform Cartesian state - constant density,
!> constant pressure, constant Cartesian velocity, which is an exact steady
!> solution of the Euler equations - and 'blast' simply adds an over-pressured
!> sphere to it. Written in spherical components the velocity is a non-trivial
!> function of theta and phi, so keeping the uniform state uniform exercises
!> every curvilinear flux and every geometric source term.
!>
!> What all three have in common is that the domain runs onto the polar axis at
!> both ends, where the ghost cells come from the block half a revolution away
!> in phi rather than from a boundary condition. Note that the volume averages
!> in the log are almost blind to whether that copy is right - a cell at the
!> axis has vanishing volume, the face on the axis has zero area, and the TVD
!> limiter clips what is left - so for the uniform setups the real test is
!> check_pole_ghosts below, which compares the ghost cells themselves against
!> the exact solution. The blast is the exception: a strong shock crossing the
!> axis leaves the state genuinely different on either side, and its log *is*
!> sensitive to the pole copy.
!>
!> This is the logarithmically stretched variant of tests/hd/spherical_pole:
!> geometry is 'logSpherical', so the radial mesh is uniform in ln(r) while
!> theta still runs through both poles. It is the only directory in the tree
!> where a stretched radius meets the polar axis, which is the combination it
!> exists to cover - the two are independent (LOG_RADIUS touches only the
!> radial coordinate, the pole only theta and phi), but nothing else tests
!> them together.
!>
!> The physical domains are those of tests/hd/spherical_pole unchanged, so the
!> two directories are directly comparable; only the radial *spacing* differs.
!> movie.par spans r in [0.4, 2.8], a factor of seven in cell width from the
!> inner edge to the outer.
!>
!> check_pole_ghosts additionally calls check_log_grid, which asserts the mesh
!> against its closed form cell by cell at it = 0.
!>
!> One consequence shows up in check_blast_fits below: read_par_files stores
!> the radial domain bounds as ln(r), so comparing them against a physical
!> radius needs an exp.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> which initial condition to run: 'uniform' or 'blast'. Read from the
  !> &usr_list namelist of the par file; see params_read_user.
  character(len=20) :: setup = 'uniform'

  !> ambient state, shared by both setups
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
  !> the uniform Cartesian background velocity. Set per setup in usr_init: the
  !> blast wants one pointing mostly along +z, so that the hot spot is carried
  !> onto the axis rather than only expanding into it.
  double precision :: v0(3) = [1.0d0, 0.5d0, -0.3d0]

  !> 'blast' only: over-pressure, radius and centre of the initial hot spot.
  !> theta_c is small but non-zero, so the spot sits near the north pole
  !> without straddling it - check_blast_fits verifies that.
  double precision :: pblast = 1.0d1
  double precision :: rblast = 0.25d0
  double precision :: rc     = 1.0d0
  double precision :: thetac = 0.10d0 * dpi
  double precision :: phic   = 0.50d0 * dpi

  !> the device-side form of `setup`: specialbound_usr and analytic_state run
  !> on the GPU, where comparing a character string is awkward
  logical :: is_blast = .false.
  !$acc declare copyin(rho0, p0, v0, pblast, rblast, rc, thetac, phic, is_blast)

contains

  subroutine usr_init()
    use mod_comm_lib, only: mpistop

    ! the coordinate system comes from `geometry` in &meshlist

    call params_read_user(par_files)

    select case (trim(setup))
    case ('uniform')
       is_blast = .false.
    case ('blast')
       is_blast = .true.
       ! mostly along +z, so the spot is advected onto the north pole
       v0 = [0.1d0, -0.05d0, 0.35d0]
    case default
       call mpistop("&usr_list setup must be 'uniform' or 'blast'")
    end select
    !$acc update device(is_blast, v0)

    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_process_grid  => check_pole_ghosts

    call phys_activate()

  end subroutine usr_init

  !> Read the &usr_list namelist from every par file given on the command
  !> line. par_files is filled by read_arguments, which runs before usr_init.
  subroutine params_read_user(files)
    use mod_global_parameters, only: unitpar
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ setup

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read_user

  !> Analytic state at x, in conserved variables.
  !>
  !> Valid at any coordinates, the ghost cells beyond the axis included, whose
  !> stored theta is negative. That is not a coincidence, and it is what makes
  !> check_pole_ghosts possible: the Cartesian point built from (r, -theta,
  !> phi) is the very point at (r, theta, phi + pi) that the pole copy fetches
  !> from, and writing a uniform Cartesian vector in the ghost cell's own basis
  !> reproduces the sign flips the copy applies. So this one routine is the
  !> initial condition, the boundary condition and the reference the pole ghost
  !> cells are checked against.
  pure subroutine analytic_state(x_loc, wpt)
    !$acc routine seq
    double precision, intent(in)  :: x_loc(1:ndim)
    double precision, intent(out) :: wpt(1:nw)
    ! .. local ..
    double precision :: sint, cost, sinp, cosp
    double precision :: xcart(3), xc(3), v(3), pres

    sint = sin(x_loc(2)); cost = cos(x_loc(2))
    sinp = sin(x_loc(3)); cosp = cos(x_loc(3))

    pres = p0
    if (is_blast) then
       ! the spot is spherical in physical space, so place it by Cartesian
       ! distance
       xcart(1) = x_loc(1) * sint * cosp
       xcart(2) = x_loc(1) * sint * sinp
       xcart(3) = x_loc(1) * cost
       xc(1)    = rc * sin(thetac) * cos(phic)
       xc(2)    = rc * sin(thetac) * sin(phic)
       xc(3)    = rc * cos(thetac)
       if (sum((xcart - xc)**2) <= rblast**2) pres = pblast
    end if

    ! spherical components of the uniform Cartesian velocity
    v(1) =  v0(1)*sint*cosp + v0(2)*sint*sinp + v0(3)*cost
    v(2) =  v0(1)*cost*cosp + v0(2)*cost*sinp - v0(3)*sint
    v(3) = -v0(1)*sinp      + v0(2)*cosp

    wpt(1:nw)      = 0.0d0
    wpt(iw_rho)    = rho0
    wpt(iw_mom(1)) = rho0 * v(1)
    wpt(iw_mom(2)) = rho0 * v(2)
    wpt(iw_mom(3)) = rho0 * v(3)
    wpt(iw_e)      = pres * inv_gamma_1 + 0.5d0 * rho0 * sum(v(1:3)**2)

  end subroutine analytic_state

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    ! .. local ..
    double precision                :: wpt(1:nw), x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    if (is_blast) call check_blast_fits

    ! analytic_state returns conserved variables directly, so there is no
    ! phys_to_conserved call here
    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call analytic_state(x_loc, wpt)
             w(ix1,ix2,ix3,1:nw) = wpt(1:nw)
          end do
       end do
    end do

  end subroutine initonegrid_usr

  !> Abort unless the hot spot starts wholly inside the radial domain and
  !> wholly off the axis.
  !>
  !> Reaching the axis is not an error here - the whole point of the blast is
  !> that its shell crosses it - but the *initial* spot has to stay clear of
  !> it, so that the state at t = 0 is smooth across the axis and
  !> check_pole_ghosts can compare against an unambiguous analytic value. The
  !> shell is free to cross it a few steps later.
  subroutine check_blast_fits
    use mod_comm_lib, only: mpistop
    ! .. local ..
    double precision :: axis_distance

    ! geometry is 'logSpherical', so read_par_files has replaced the radial
    ! domain bounds by their logarithms - undo that to compare against a
    ! physical radius.
    if (rc - rblast < exp(xprobmin1) .or. rc + rblast > exp(xprobmax1)) &
       call mpistop("log_spherical_pole blast: hot spot sticks out of the &
          &radial domain")

    ! perpendicular distance from the centre of the spot to the polar axis
    axis_distance = rc * sin(thetac)
    if (axis_distance <= rblast) &
       call mpistop("log_spherical_pole blast: hot spot already straddles the &
          &polar axis at t = 0; move it out in theta")

  end subroutine check_blast_fits

  !> Keep the exact solution in the ghost cells at the radial boundaries, so
  !> that only the interior discretisation is under test. Used by the uniform
  !> setups; the blast leaves r open with 'cont' instead.
  subroutine specialbound_usr(qt, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB, w, x)
    !$acc routine vector
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3
    integer, intent(in)             :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    integer, intent(in)             :: iB
    double precision, intent(in)    :: qt
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    ! .. local ..
    double precision                :: wpt(1:nw), x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    ! Vector on the innermost loop, deliberately not collapse(3): nvfortran's
    ! OpenACC miscompiles a collapsed vector loop that calls another !$acc
    ! routine in its body, reached through this !$acc routine vector from the
    ! gang loops in fill_boundary_before_gc / fill_boundary_after_gc. It bit
    ! the hd curvilinear pole cases via analytic_state; the fix is to drop
    ! the collapse or inline the callee. See CLAUDE.md - very likely the
    ! same defect as issue #154.
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          !$acc loop vector private(wpt, x_loc)
          do ix1 = ixOmin1, ixOmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call analytic_state(x_loc, wpt)
             w(ix1,ix2,ix3,1:nw) = wpt(1:nw)
          end do
       end do
    end do

  end subroutine specialbound_usr

  !> Refine the half of the domain with phi < pi, and only that half.
  !>
  !> Used by uflow_amr.par, which sets refine_usr = .true. and
  !> refine_criterion = 1. The point is the level jump this leaves along
  !> phi = 0 and phi = pi: a block touching the axis inside the wedge has its
  !> pi-periodic partner outside it, one level coarser, so the pole copy has to
  !> go through the restriction and prolongation paths rather than the
  !> same-level one.
  !>
  !> Called by name rather than through a procedure pointer: forcedrefine_grid
  !> reaches it from inside an OpenACC kernel, so it is a compile-time
  !> dependency, which is why refine_usr is .true. in agile.par even though
  !> only one of the three par files switches it on at run time.
  subroutine usr_refine_grid(igrid,level,ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
     ixGmax2,ixGmax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,qt,w,x,refine,&
     coarsen)
    !$acc routine vector
    use mod_global_parameters
    integer, intent(in)          :: igrid, level, ixGmin1,ixGmin2,ixGmin3,&
       ixGmax1,ixGmax2,ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: qt, w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw), x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    if (x(ixmin1,ixmin2,ixmin3,3) < dpi) then
       refine  =  1
       coarsen = -1
    else
       refine  = -1
       coarsen =  1
    end if

  end subroutine usr_refine_grid

  !> Assert that the ghost cells across the axis carry the exact solution.
  !>
  !> Only meaningful at it = 0, where the interior is still exactly the
  !> analytic state, so the pole copy of it has to reproduce that state in the
  !> ghost cells to round-off. Silence means the check passed.
  !>
  !> For the 'uniform' setup the transverse loops run over the whole block
  !> (ixI), not just its interior mesh (ixO), so the edge and corner cells of
  !> the pole layer are covered too - in particular where the axis meets the
  !> radial physical boundary, filled by bc_phys rather than the pole copy.
  !> 'blast' keeps the face-only range: near its discontinuous spot surface a
  !> ghost cell holds a limited cell average that differs from the point value.
  !>
  !> How far it can be pushed across a *level jump* depends on the setup,
  !> because there the ghost has been restricted or prolonged on the way and
  !> comparing it against the analytic value at a point is only meaningful
  !> where the field is smooth on the scale of a coarse cell:
  !>
  !>  * 'uniform' is smooth everywhere, so both branches are checked -
  !>    round-off for a same-level pole neighbour and loose across a level
  !>    jump, where the measured error is 3.4e-2 against the 0.18 a wrong sign
  !>    produces.
  !>
  !>  * 'blast' has a discontinuous spot surface, and a coarse cell straddling
  !>    it holds the 2:1 average of hot and cold gas, which differs from the
  !>    analytic value at its centre by a large fraction of the jump however
  !>    correct the pole copy is. So it checks same-level neighbours only, and
  !>    covers the restricting and prolonging paths through its log - which it
  !>    can afford to do because its log is the one that is sensitive to them.
  subroutine check_pole_ghosts(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    ! .. local ..
    double precision :: wpt(1:nw), x_loc(1:ndim), err, tol
    integer          :: ix1, ix2, ix3, iside, i2, jxmin2, jxmax2
    integer          :: tx1lo, tx1hi, tx3lo, tx3hi

    if (it /= 0) return

    ! This directory is the one place where a logarithmically stretched radius
    ! meets the polar axis, so assert both: the mesh itself, and the pole copy
    ! running on it.
    call check_log_grid(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)

    ! process() runs before the solution is pulled back for output, so the
    ! host copy of w is stale unless we fetch this block ourselves
    !$acc update host(ps(igrid)%w)

    ! 'uniform' is smooth, so check the whole ghost layer including its edges
    ! and corners; 'blast' keeps the interior-transverse face only (see the
    ! docstring).
    if (is_blast) then
       tx1lo = ixOmin1; tx1hi = ixOmax1; tx3lo = ixOmin3; tx3hi = ixOmax3
    else
       tx1lo = ixImin1; tx1hi = ixImax1; tx3lo = ixImin3; tx3hi = ixImax3
    end if

    do iside = 1, 2
       i2 = 2*iside - 3
       if (neighbor_pole(0,i2,0,igrid) == 0) cycle
       if (neighbor_type(0,i2,0,igrid) == neighbor_sibling) then
          tol = 1.0d-10
       else
          if (is_blast) cycle
          tol = 1.0d-1
       end if
       if (iside == 1) then
          jxmin2 = ixImin2;   jxmax2 = ixOmin2-1
       else
          jxmin2 = ixOmax2+1; jxmax2 = ixImax2
       end if
       err = 0.0d0
       do ix3 = tx3lo, tx3hi
          do ix2 = jxmin2, jxmax2
             do ix1 = tx1lo, tx1hi
                x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
                call analytic_state(x_loc, wpt)
                err = max(err, maxval(abs(w(ix1,ix2,ix3,&
                   1:nwflux) - wpt(1:nwflux))))
             end do
          end do
       end do
       if (err > tol) then
          write(*,*) 'pole ghost cells deviate from the exact solution by',&
             err, ' tolerance', tol
          call mpistop('pole ghost-cell check failed')
       end if
    end do

  end subroutine check_pole_ghosts


  !> Assert the logarithmically stretched mesh against its closed form, at
  !> it = 0.
  !>
  !> A wrong metric on a stretched grid is close to invisible in the volume
  !> averages the log reports - the same reason the polar-axis cases check
  !> ghost cells directly - so the mesh is checked cell by cell here instead.
  !> Four things, all to round-off:
  !>
  !>   1. consecutive radial positions are in exact geometric progression,
  !>      x(i+1)/x(i) = q = exp(dxi). That is the defining property of a
  !>      logarithmic grid, and it holds for the volume *barycentres* stored in
  !>      x as well as for the faces, because both faces of every cell scale by
  !>      the same q. It is also the property that fails loudly if the
  !>      exponential is ever dropped or applied twice.
  !>   2. ds(1) is the physical width of the cell, and ds(1)/x(1) is therefore
  !>      the same constant everywhere - constant relative resolution.
  !>   3. each stored position is the volume barycentre of its own two faces.
  !>      The faces are rebuilt here from rnode, and the barycentre from the
  !>      raw quotient of quartic (spherical) or cubic (cylindrical)
  !>      differences - deliberately a different expression from the
  !>      cancellation-free one fill_geometry_device evaluates, so that the two
  !>      agreeing means something.
  !>   4. dvolume is the exact integral over those same faces, with the polar
  !>      factor written as cos(theta_L) - cos(theta_R) rather than the
  !>      2 sin(theta_c) sin(dtheta/2) the kernel uses.
  subroutine check_log_grid(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    ! .. local ..
    double precision :: d1, d2, d3, xlo1, xlo2, q, fL, fR, rbar, dvol
    double precision :: tL, tR, err
    double precision, parameter :: tol = 1.0d-11
    integer          :: ix1, ix2, ix3

    if (it /= 0) return
    ! process() syncs the positions for us, but not the metrics
    !$acc update host(ps(igrid)%ds, ps(igrid)%dvolume)

    d1 = rnode(rpdx1_,igrid); d2 = rnode(rpdx2_,igrid); d3 = rnode(rpdx3_,igrid)
    xlo1 = rnode(rpxmin1_,igrid); xlo2 = rnode(rpxmin2_,igrid)
    q = exp(d1)

    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

             ! the two radial faces, from the block corner in ln(r)
             fL = exp(xlo1 + dble(ix1-nghostcells-1)*d1)
             fR = exp(xlo1 + dble(ix1-nghostcells)*d1)
             tL = xlo2 + dble(ix2-nghostcells-1)*d2
             tR = tL + d2

             ! (1) geometric progression of the stored positions
             if (ix1 < ixOmax1) then
                err = abs(x(ix1+1,ix2,ix3,1)/x(ix1,ix2,ix3,1) - q)/q
                call assert_close(err, tol, 'radial positions are not in &
                   &geometric progression', x(ix1,ix2,ix3,1))
             end if

             ! (2) constant relative resolution
             err = abs(ps(igrid)%ds(ix1,ix2,ix3,1) - (fR-fL))/(fR-fL)
             call assert_close(err, tol, 'ds(1) is not the physical cell &
                &width', x(ix1,ix2,ix3,1))

             ! (3) the stored position is the volume barycentre
             rbar = 0.75d0*(fR**4 - fL**4)/(fR**3 - fL**3)
             err = abs(x(ix1,ix2,ix3,1) - rbar)/rbar
             call assert_close(err, tol, 'x(1) is not the volume barycentre &
                &of its faces', x(ix1,ix2,ix3,1))

             ! (4) the cell volume is the exact integral over those faces
             dvol = (fR**3 - fL**3)/3.0d0 * (cos(tL) - cos(tR)) * d3
             err = abs(ps(igrid)%dvolume(ix1,ix2,ix3) - dvol)/dvol
             call assert_close(err, tol, 'dvolume is not the exact cell &
                &volume', x(ix1,ix2,ix3,1))

          end do
       end do
    end do

  end subroutine check_log_grid

  !> Abort with a located message if err exceeds tol.
  subroutine assert_close(err, tol, what, r)
    use mod_global_parameters, only: mype, unitterm
    use mod_comm_lib, only: mpistop
    double precision, intent(in) :: err, tol, r
    character(len=*), intent(in) :: what

    if (err <= tol) return
    if (mype == 0) write(unitterm,*) 'check_log_grid: ', trim(what), &
       ' at r =', r, ' relative error', err
    call mpistop('check_log_grid: the logarithmic mesh is wrong')

  end subroutine assert_close

end module mod_usr
