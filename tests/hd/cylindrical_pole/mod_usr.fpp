!> Polar-axis tests on a cylindrical mesh: three runs from one build.
!>
!> The cylindrical counterpart of tests/hd/spherical_pole, and organised the
!> same way. Every test directory rebuilds AGILE from scratch, so the suite's
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

    ! in cylindrical the perpendicular distance to the axis is just r_c
    if (rc <= rblast) &
       call mpistop("cylindrical_pole blast: hot spot already straddles the &
          &axis at t = 0; move it out in r")

    if (rc + rblast > xprobmax1) &
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

    ! !$acc loop vector on the innermost loop only, deliberately not
    ! collapse(3): nvfortran's OpenACC miscompiles a collapsed vector loop
    ! inside this !$acc routine vector (reached from the gang loop in
    ! fill_boundary_before_gc / fill_boundary_after_gc) when the boundary
    ! region is corner-shaped - a ghost slab only nghostcells wide - giving
    ! wrong un-collapsed indices, so the ghost cells where the physical
    ! boundary meets the polar axis come out rotated. See CLAUDE.md.
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
    ! the cylindrical axis is the r = 0 face only
    if (neighbor_pole(-1,0,0,igrid) == 0) return
    if (neighbor_type(-1,0,0,igrid) == neighbor_sibling) then
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

end module mod_usr
