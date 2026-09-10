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

    if (rc - rblast < xprobmin1 .or. rc + rblast > xprobmax1) &
       call mpistop("spherical_pole blast: hot spot sticks out of the radial &
          &domain")

    ! perpendicular distance from the centre of the spot to the polar axis
    axis_distance = rc * sin(thetac)
    if (axis_distance <= rblast) &
       call mpistop("spherical_pole blast: hot spot already straddles the &
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
    ! wpt: compile-time nw_phys, NOT the runtime nw - see the note on the loop
    double precision                :: wpt(1:nw_phys), x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    ! nvfortran's OpenACC miscompiles a collapsed vector loop whose body has a
    ! call AND a runtime-sized private automatic array, reached through this
    ! !$acc routine vector from the gang loops in fill_boundary_before_gc /
    ! fill_boundary_after_gc: it gets the array's per-lane stride wrong and the
    ! ghost layer comes back index-rotated. Sizing wpt with the compile-time
    ! nw_phys is the minimal fix (dropping the collapse also works but still
    ! reorders the short vector loop). It bit the hd curvilinear pole cases via
    ! analytic_state. See CLAUDE.md ("Bug-hunting notes") and issue #154.
    !$acc loop collapse(3) vector private(wpt, x_loc)
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
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
  !> The transverse loops run over the whole block (ixI), not just its
  !> interior mesh (ixO), so the ghost cells at the *edges* and *corners* of
  !> the pole layer are covered too - in particular the cells where the axis
  !> meets the radial physical boundary, which a face-only check would miss.
  !> Only the interior-transverse *face* is held to round-off (the pole copy
  !> of an exactly-analytic interior); the edge/corner cells the widening
  !> adds are always checked loose (1e-1), because a corner cell can be
  !> filled by prolongation from a coarser neighbour in another direction,
  !> where an O(dx^2) difference from the point analytic value is expected.
  !> That loose bound still catches a wrong sign at the axis (O(0.1-1)).
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
    double precision :: wpt(1:nw), x_loc(1:ndim), err, errc, tol
    integer          :: ix1, ix2, ix3, iside, i2, jxmin2, jxmax2
    logical          :: face

    if (it /= 0) return
    ! process() runs before the solution is pulled back for output, so the
    ! host copy of w is stale unless we fetch this block ourselves
    !$acc update host(ps(igrid)%w)

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
       ! err: the interior-transverse face, checked at `tol` (round-off for a
       ! same-level pole copy). errc: the edge/corner cells the widening
       ! added, always checked loose - they can be filled by prolongation
       ! from a coarser neighbour in another direction, where an O(dx^2)
       ! difference from the point analytic value is expected, not a bug.
       err = 0.0d0; errc = 0.0d0
       do ix3 = ixImin3, ixImax3
          do ix2 = jxmin2, jxmax2
             do ix1 = ixImin1, ixImax1
                face = ix1>=ixOmin1 .and. ix1<=ixOmax1 .and. &
                       ix3>=ixOmin3 .and. ix3<=ixOmax3
                x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
                call analytic_state(x_loc, wpt)
                if (face) then
                   err  = max(err,  maxval(abs(w(ix1,ix2,ix3,1:nwflux) - wpt(1:nwflux))))
                else
                   errc = max(errc, maxval(abs(w(ix1,ix2,ix3,1:nwflux) - wpt(1:nwflux))))
                end if
             end do
          end do
       end do
       if (err > tol .or. errc > 1.0d-1) then
          write(*,*) 'pole ghost cells deviate from the exact solution by',&
             max(err, errc), ' tolerances', tol, 1.0d-1
          call mpistop('pole ghost-cell check failed')
       end if
    end do

  end subroutine check_pole_ghosts

end module mod_usr
