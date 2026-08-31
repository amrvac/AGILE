!> Polar-axis tests on a logarithmically stretched cylindrical mesh: three
!> runs from one build.
!>
!> tests/hd/cylindrical_pole again, on geometry = 'logCylindrical' with a
!> non-zero log_r0, and organised the same way. Every test directory rebuilds AGILE from scratch, so the suite's
!> cost is dominated by compilation rather than by the runs; cases that agree
!> on the *compile-time* parameters - `phys`, `geometry`, `specialboundary` and
!> `refine_usr`, the ones make/config_schema.toml turns into fypp defines - can
!> share one build and differ only in their par file. These do, so they pick
!> their initial condition through the `setup` entry of the &usr_list namelist:
!>
!>   uflow.par       setup = 'uniform'  a uniform Cartesian flow, single level
!>   uflow_amr.par   setup = 'uniform'  the same, with half the domain in phi
!>                                      refined so blocks meet across the axis
!>                                      at different levels
!>   blast_amr.par   setup = 'blast'    a blast wave that starts off the axis
!>                                      and expands across it, Lohner-refined
!>
!> Both setups are the same uniform Cartesian state - constant density,
!> constant pressure, constant Cartesian velocity, an exact steady solution of
!> the Euler equations - and 'blast' adds an over-pressured sphere
!> to it. Written in cylindrical (r, z, phi) components the radial and
!> azimuthal velocity are non-trivial functions of phi (only the axial
!> component stays literally constant), so keeping the uniform state uniform
!> exercises every curvilinear flux and every geometric source term.
!>
!> All three run onto the cylindrical axis at r = 0, where the ghost cells come
!> from the block half a revolution away in phi rather than from a boundary
!> condition. The volume averages in the log are almost blind to whether that
!> copy is right - a cell at the axis has vanishing volume, the face on the
!> axis has zero area, and the TVD limiter clips what is left - so for the
!> uniform setups the real test is check_pole_ghosts below. The blast is the
!> exception: a strong shock crossing the axis leaves the state genuinely
!> different on either side, and its log *is* sensitive to the pole copy.
!>
!> These cases are also what pin down the sign of the radial momentum at the
!> cylindrical axis, where upstream MPI-AMRVAC and this fork disagree - see
!> "The polar axis" in CLAUDE.md.
!>
!> What this directory adds over tests/hd/cylindrical_pole is the *offset*
!> logarithmic radial map, s = ln(1 + r/r0), which is the only radial
!> stretching that reaches r = 0: the plain s = ln(r) of the other log_*
!> directories has no value there. So this is the one place in the tree where
!> a stretched radius meets the cylindrical axis, and check_pole_ghosts here
!> calls two extra checks at it = 0:
!>
!>   check_log_grid   the stretched mesh against its closed form, cell by
!>                    cell, as in tests/hd/log_cylindrical - but rewritten for
!>                    the offset map
!>   check_pole_mesh  the mesh in the ghost layer beyond the axis against the
!>                    exact mirror of the mesh inside it
!>
!> check_pole_mesh is the assertion this directory exists for. The pole copy
!> in getbc fills ghost cell j from interior cell j of the block half a
!> revolution away, so the two have to occupy mirrored volumes or the copied
!> cell average belongs to a cell of a different size. That holds because the
!> map is *odd* - s = 0 is r = 0, and r_of_s carries the sign through - and it
!> would not hold for the obvious alternative r = exp(s) - r0, whose ghost
!> cells shrink outward while the cells they mirror grow. With that map the
!> first ghost cell here would be misplaced by about 5% of a cell width, which
!> check_pole_ghosts' 1e-10 comparison against the analytic state would catch
!> and the volume averages in the log would not.
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
  !> blast wants one pointing mostly along -y, which is across the axis from
  !> the hot spot at phi = pi/2, so the spot is carried onto the axis rather
  !> than only expanding into it.
  double precision :: v0(3) = [1.0d0, 0.5d0, -0.3d0]

  !> 'blast' only: over-pressure, radius and centre of the initial hot spot, in
  !> cylindrical (r, z, phi). r_c is small but larger than the spot radius, so
  !> the spot sits near the axis without straddling it - check_blast_fits
  !> verifies that.
  double precision :: pblast = 1.0d1
  double precision :: rblast = 0.2d0
  double precision :: rc     = 0.3d0
  double precision :: zc     = 0.5d0
  double precision :: phic   = 0.50d0 * dpi

  !> the device-side form of `setup`: specialbound_usr and analytic_state run
  !> on the GPU, where comparing a character string is awkward
  logical :: is_blast = .false.
  !$acc declare copyin(rho0, p0, v0, pblast, rblast, rc, zc, phic, is_blast)

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
       ! mostly along -y, so the spot is advected across the axis
       v0 = [0.05d0, -0.35d0, 0.1d0]
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
  !> stored r is negative. That is not a coincidence, and it is what makes
  !> check_pole_ghosts possible: the Cartesian point built from (-r, z, phi) is
  !> the very point at (r, z, phi + pi) that the pole copy fetches from, and
  !> writing a uniform Cartesian vector in the ghost cell's own basis
  !> reproduces the sign flips the copy applies. So this one routine is the
  !> initial condition, the boundary condition and the reference the pole ghost
  !> cells are checked against.
  pure subroutine analytic_state(x_loc, wpt)
    !$acc routine seq
    double precision, intent(in)  :: x_loc(1:ndim)
    double precision, intent(out) :: wpt(1:nw)
    ! .. local ..
    double precision :: sinp, cosp
    double precision :: xcart(3), xc(3), v(3), pres

    sinp = sin(x_loc(3)); cosp = cos(x_loc(3))

    pres = p0
    if (is_blast) then
       ! the spot is spherical in physical space, so place it by Cartesian
       ! distance
       xcart(1) = x_loc(1) * cosp
       xcart(2) = x_loc(1) * sinp
       xcart(3) = x_loc(2)
       xc(1)    = rc * cos(phic)
       xc(2)    = rc * sin(phic)
       xc(3)    = zc
       if (sum((xcart - xc)**2) <= rblast**2) pres = pblast
    end if

    ! cylindrical (r, z, phi) components of the uniform Cartesian velocity
    v(1) =  v0(1)*cosp + v0(2)*sinp
    v(2) =  v0(3)
    v(3) = -v0(1)*sinp + v0(2)*cosp

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

  !> Abort unless the hot spot starts wholly inside the domain and wholly off
  !> the axis.
  !>
  !> Reaching the axis is not an error here - the whole point of the blast is
  !> that its shell crosses it - but the *initial* spot has to stay clear of
  !> it, so that the state at t = 0 is smooth across the axis and
  !> check_pole_ghosts can compare against an unambiguous analytic value. The
  !> shell is free to cross it a few steps later.
  subroutine check_blast_fits
    use mod_comm_lib, only: mpistop
    use mod_geometry, only: r_of_s

    ! in cylindrical the perpendicular distance to the axis is just r_c
    if (rc <= rblast) &
       call mpistop("cylindrical_pole blast: hot spot already straddles the &
          &axis at t = 0; move it out in r")

    ! xprobmax1 holds the *logical* radial bound in a LOG_RADIUS build, so it
    ! has to go back through the map before it can be compared with a radius
    if (rc + rblast > r_of_s(xprobmax1)) &
       call mpistop("cylindrical_pole blast: hot spot sticks out in r")

    if (zc - rblast < xprobmin2 .or. zc + rblast > xprobmax2) &
       call mpistop("cylindrical_pole blast: hot spot sticks out in z")

  end subroutine check_blast_fits

  !> Keep the exact solution in the ghost cells at the boundaries that are not
  !> the axis, so that only the interior discretisation is under test. Used by
  !> the uniform setups; the blast leaves them open with 'cont' instead.
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
  !> the pole layer are covered too - in particular where the axis meets the z
  !> physical boundary, filled by bc_phys rather than the pole copy. 'blast'
  !> keeps the face-only range: near its discontinuous spot surface a ghost
  !> cell holds a limited cell average that differs from the point value.
  !>
  !> How far it can be pushed across a *level jump* depends on the setup,
  !> because there the ghost has been restricted or prolonged on the way and
  !> comparing it against the analytic value at a point is only meaningful
  !> where the field is smooth on the scale of a coarse cell: 'uniform' is
  !> smooth everywhere, so both branches are checked, while 'blast' has a
  !> discontinuous spot surface, where a coarse cell straddling it holds the
  !> 2:1 average of hot and cold gas however correct the pole copy is. The
  !> blast covers the restricting and prolonging paths through its log instead,
  !> which it can afford to do because its log is sensitive to them.
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
    integer          :: ix1, ix2, ix3, tx2lo, tx2hi, tx3lo, tx3hi

    if (it /= 0) return

    ! the stretched mesh itself, on every block of the domain
    call check_log_grid(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)

    ! the cylindrical axis is the r = 0 face only
    if (neighbor_pole(-1,0,0,igrid) == 0) return
    if (neighbor_type(-1,0,0,igrid) == neighbor_sibling) then
       ! the geometry the pole copy assumes, on the blocks that do the copy
       call check_pole_mesh(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax2,ixOmax3,x)
       tol = 1.0d-10
    else
       if (is_blast) return
       tol = 1.0d-1
    end if
    ! process() runs before the solution is pulled back for output, so the
    ! host copy of w is stale unless we fetch this block ourselves
    !$acc update host(ps(igrid)%w)

    ! 'uniform' is smooth, so check the whole ghost layer including its edges
    ! and corners (where the axis meets the z physical boundary); 'blast'
    ! keeps the interior-transverse face only (see the docstring).
    if (is_blast) then
       tx2lo = ixOmin2; tx2hi = ixOmax2; tx3lo = ixOmin3; tx3hi = ixOmax3
    else
       tx2lo = ixImin2; tx2hi = ixImax2; tx3lo = ixImin3; tx3hi = ixImax3
    end if

    err = 0.0d0
    do ix3 = tx3lo, tx3hi
       do ix2 = tx2lo, tx2hi
          do ix1 = ixImin1, ixOmin1-1
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

  end subroutine check_pole_ghosts

  !> The offset logarithmic radial map, r as a function of the logical
  !> coordinate s, for s >= 0.
  !>
  !> Deliberately a different expression from the one fill_geometry_device
  !> evaluates: 2*r0*sinh(s/2)*exp(s/2) is exactly r0*(exp(s) - 1), but reaches
  !> it by a different rounding path, so the two agreeing is evidence rather
  !> than a tautology. Only the physical domain is covered, where s >= 0; the
  !> odd branch the real map takes beyond the axis is checked separately, by
  !> check_pole_mesh, and against the mesh itself rather than against a
  !> formula.
  pure function face_r(s) result(r)
    use mod_global_parameters, only: log_r0
    double precision, intent(in) :: s
    double precision             :: r

    r = 2.0d0*log_r0*sinh(0.5d0*s)*exp(0.5d0*s)

  end function face_r

  !> Assert that the stretched mesh is what its closed form says it is.
  !>
  !> A mesh that is subtly wrong barely moves a volume average - the same
  !> reason this directory asserts on ghost cells rather than on the log - so
  !> the mesh is compared against its closed form cell by cell at it = 0.
  subroutine check_log_grid(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    ! .. local ..
    double precision :: d1, d2, d3, xlo1, q, fL, fR, rbar, dvol, err
    double precision, parameter :: tol = 1.0d-11
    integer          :: ix1, ix2, ix3

    if (it /= 0) return
    ! process() syncs the positions for us, but not the metrics
    !$acc update host(ps(igrid)%ds, ps(igrid)%dvolume)

    d1 = rnode(rpdx1_,igrid); d2 = rnode(rpdx2_,igrid); d3 = rnode(rpdx3_,igrid)
    xlo1 = rnode(rpxmin1_,igrid)
    q = exp(d1)

    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

             ! the two radial faces, from the block corner in the logical
             ! coordinate s = ln(1 + r/r0)
             fL = face_r(xlo1 + dble(ix1-nghostcells-1)*d1)
             fR = face_r(xlo1 + dble(ix1-nghostcells)*d1)

             ! (1) the cell width grows by exactly exp(ds) per cell, which is
             ! what "uniform in ln(r + r0)" means and is the one statement
             ! about the mesh that does not go through face_r at all
             if (ix1 < ixOmax1) then
                err = abs(ps(igrid)%ds(ix1+1,ix2,ix3,1) &
                     /ps(igrid)%ds(ix1,ix2,ix3,1) - q)/q
                call assert_close(err, tol, 'radial cell widths are not in &
                   &geometric progression', x(ix1,ix2,ix3,1))
             end if

             ! (2) ds(1) is the physical cell width
             err = abs(ps(igrid)%ds(ix1,ix2,ix3,1) - (fR-fL))/(fR-fL)
             call assert_close(err, tol, 'ds(1) is not the physical cell &
                &width', x(ix1,ix2,ix3,1))

             ! (3) the stored position is the volume barycentre of its own two
             ! faces - as the raw quotient, not the cancellation-free rewrite
             ! the geometry kernel uses
             rbar = (2.0d0/3.0d0)*(fR**3 - fL**3)/(fR**2 - fL**2)
             err = abs(x(ix1,ix2,ix3,1) - rbar)/rbar
             call assert_close(err, tol, 'x(1) is not the volume barycentre &
                &of its faces', x(ix1,ix2,ix3,1))

             ! (4) the cell volume is the exact integral over those faces
             dvol = 0.5d0*(fR**2 - fL**2) * d2 * d3
             err = abs(ps(igrid)%dvolume(ix1,ix2,ix3) - dvol)/dvol
             call assert_close(err, tol, 'dvolume is not the exact cell &
                &volume', x(ix1,ix2,ix3,1))

          end do
       end do
    end do

  end subroutine check_log_grid

  !> Assert that the mesh beyond the axis is the exact mirror of the mesh
  !> inside it.
  !>
  !> This is the property the pole copy in getbc rests on, and the reason the
  !> radial map is odd rather than the naive r = exp(s) - r0: getbc fills ghost
  !> cell j from interior cell j of the block half a revolution away, so unless
  !> the two occupy mirrored volumes the copied cell average belongs to a cell
  !> of a different size, and the ghost is misplaced by an O(ds) fraction of a
  !> cell width however correct the copy itself is.
  !>
  !> The mirror is exact in exact arithmetic, and nearly so in practice.
  !> xprobmin1 is ln(1 + 0/r0) = 0 exactly, rnode's corner for the block on
  !> the axis is therefore 0 exactly, so the logical positions of a ghost cell
  !> and the interior cell it mirrors are exact negatives of one another;
  !> abs() and sign() in the map then carry that through, and the barycentre
  !> and volume expressions are even or odd in the pair of faces. Under a
  !> compiler that respects IEEE semantics (gfortran at -O3) the two sides
  !> agree bit-for-bit; ifx at -O3 defaults to -fp-model=fast and moves the
  !> closing bits, so assert_exact compares to 1e-11 rather than exactly -
  !> still nine orders of magnitude tighter than the error the naive map
  !> would introduce.
  subroutine check_pole_mesh(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax2,ixOmax3,x)
    use mod_global_parameters
    integer, intent(in)          :: igrid,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! .. local ..
    integer :: ix1, ix2, ix3, jx1

    if (it /= 0) return
    !$acc update host(ps(igrid)%ds, ps(igrid)%dvolume)

    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = 1, nghostcells
             ! the ghost cell ix1 layers below the axis, and the interior cell
             ! it is the mirror image of
             jx1 = ixOmin1 - ix1
             call assert_exact(x(jx1,ix2,ix3,1), -x(ixOmin1+ix1-1,ix2,ix3,1), &
                'a ghost cell beyond the axis is not the mirror of the &
                &interior cell it is filled from')
             call assert_exact(ps(igrid)%ds(jx1,ix2,ix3,1), &
                ps(igrid)%ds(ixOmin1+ix1-1,ix2,ix3,1), &
                'a ghost cell beyond the axis has a different width from the &
                &interior cell it is filled from')
             call assert_exact(ps(igrid)%dvolume(jx1,ix2,ix3), &
                ps(igrid)%dvolume(ixOmin1+ix1-1,ix2,ix3), &
                'a ghost cell beyond the axis has a different volume from the &
                &interior cell it is filled from')
          end do
       end do
    end do

  end subroutine check_pole_mesh

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

  !> Assert that got and want agree to a tight relative tolerance.
  !>
  !> The mesh mirror across the axis is exact in the source expressions (see
  !> check_pole_mesh above), so under a compiler that keeps to IEEE semantics
  !> got == want to the last bit. ifx at -O3 defaults to -fp-model=fast, which
  !> contracts and reassociates and uses a fast exp, and then the closing bits
  !> of the mirror move. The tolerance is 1e-11, matching check_log_grid's own
  !> tolerance for the same kind of geometry arithmetic in this file; it still
  !> catches the O(ds) misplacement the naive r = exp(s) - r0 map produces -
  !> several percent of a cell width - by nine orders of magnitude.
  subroutine assert_exact(got, want, what)
    use mod_global_parameters, only: mype, unitterm
    use mod_comm_lib, only: mpistop
    double precision, intent(in) :: got, want
    character(len=*), intent(in) :: what
    double precision, parameter  :: tol = 1.0d-11

    if (abs(got - want) <= tol * max(abs(want), tiny(1.0d0))) return
    if (mype == 0) write(unitterm,*) 'check_pole_mesh: ', trim(what), &
       ': got', got, ' expected', want
    call mpistop('check_pole_mesh: the mesh is not mirrored across the axis')

  end subroutine assert_exact

end module mod_usr
