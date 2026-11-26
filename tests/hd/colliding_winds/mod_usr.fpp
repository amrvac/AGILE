#:mute
#:include "mod_gpu_directives.fpp"
#:endmute

!> Colliding stellar winds in a binary on a circular orbit.
!>
!> Two stars orbit their common centre of mass in the x-y plane, each blowing
!> an isotropic, supersonic wind. The winds are imposed by the internal
!> boundary hook usr_internal_bc: every cell within r_in(i) of star i is held
!> at the wind state of that star at the current time, so the source regions
!> follow the orbit. The winds collide in a bow-shaped interaction region
!> that the orbital motion winds into a spiral.
!>
!> Geometry (code units, separation a):
!>   star 1 (mass m1) at -a*q/(1+q) * e(phi),  star 2 at +a/(1+q) * e(phi),
!>   q = m2/m1, e(phi) = (cos phi, sin phi, 0), phi = 2*pi*t/period + phase0.
!>
!> Wind of star i, at distance r from it (lab frame):
!>   rho = mdot/(4 pi r^2 vw),   v = v_orb + vw * r-hat,   p = rho T_w,
!>   T_w = (vw/mach)^2/gamma     (p/rho, i.e. sound speed vw/mach)
!> an isothermal wind of Mach number `mach`. Inside the source region r is
!> floored at r_in/2 so the centre stays finite.
!>
!> Temperature floor: usr_internal_bc also keeps p/rho >= t_floor in every
!> cell, not just inside the stars (default: the colder wind's T_w; set
!> temp_floor > 0 to override). Real
!> hot-star winds are held near 1e4 K by photoionisation; without the floor
!> an adiabatic Mach-20 wind cools by r^(-4/3), its internal energy becomes a
!> tiny difference of large numbers, and the run produced negative pressures
!> (NaN sound speeds that then drop out of the CFL condition) within one time
!> unit. The shocked gas between the winds is far above the floor and is not
!> affected.
!>
!> The floor sits in the internal boundary rather than in a local source term
!> (addsource_usr) so that the finite-volume kernel stays free of the
!> SOURCE_LOCAL loop. getbc runs it after every stage, before the ghost
!> exchange, the next stage and setdt read the state. What it does not touch
!> are ghost cells (prolonged or physical-boundary ones); the Riemann solver's
!> face states there are protected by hd_small_pressure (fix_prim_state).
!>
!> The initial condition fills the whole domain with the wind whose ram
!> pressure dominates locally, so the run starts with the contact discontinuity
!> roughly in place instead of having to sweep out an ambient medium.
!>
!> Tracer 1 is the mass fraction of wind 1 (conservative: rho * fraction).
!>
!> Refinement: cells within refine_radius*r_in of either star are forced to the
!> finest level, so the source regions stay resolved as they move; elsewhere
!> the Lohner estimator refines the shocks and the contact, capped at
!> level lmax_far.
module mod_usr
  use mod_agile
  use mod_physics

  implicit none

  ! --- orbit ---
  double precision :: sep     = 1.0d0    !< binary separation a
  double precision :: q_mass  = 0.5d0    !< mass ratio m2/m1
  double precision :: period  = 12.0d0   !< orbital period
  double precision :: phase0  = 0.0d0    !< initial orbital phase [rad]
  ! --- winds, (1) and (2) ---
  double precision :: mdot(2) = [1.0d0, 0.25d0]  !< mass-loss rates
  double precision :: vw(2)   = [1.0d0, 1.0d0]   !< terminal wind speeds
  double precision :: mach(2) = [20.0d0, 20.0d0] !< wind Mach number at r_in
  double precision :: r_in(2) = [0.1d0, 0.1d0]   !< source-region radii
  ! --- refinement ---
  double precision :: refine_radius = 3.0d0  !< in units of r_in
  integer          :: lmax_far      = 3      !< max level away from the stars
  ! --- thermodynamics ---
  double precision :: temp_floor = -1.0d0 !< p/rho floor; <= 0: min of T_w
  double precision :: t_floor             !< the floor in use

  ${GPU_DECLARE_COPYIN('sep, q_mass, period, phase0, mdot, vw, mach, r_in, refine_radius, lmax_far, t_floor')}$

contains

  subroutine usr_init()
    use mod_comm_lib, only: mpistop

    call params_read_user(par_files)
    if (any(r_in >= 0.5d0*sep)) &
       call mpistop('colliding_winds: r_in must be below half the separation')
    if (any(mdot <= 0.0d0) .or. any(vw <= 0.0d0) .or. any(mach <= 0.0d0)) &
       call mpistop('colliding_winds: mdot, vw and mach must be positive')
    ${GPU_UPDATE_DEVICE('sep, q_mass, period, phase0, mdot, vw, mach, r_in, refine_radius, lmax_far')}$

    usr_init_one_grid => initonegrid_usr

    call phys_activate()

    ! hd_gamma is known once phys_activate has read &hd_list
    if (temp_floor > 0.0d0) then
       t_floor = temp_floor
    else
       t_floor = minval((vw/mach)**2) / hd_gamma
    end if
    ${GPU_UPDATE_DEVICE('t_floor')}$

  end subroutine usr_init

  subroutine params_read_user(files)
    use mod_global_parameters, only: unitpar
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ sep, q_mass, period, phase0, mdot, vw, mach, r_in, &
         refine_radius, lmax_far, temp_floor

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read_user

  !> Positions and velocities of both stars at time qt.
  pure subroutine star_orbits(qt, xs, vs)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: qt
    double precision, intent(out) :: xs(3,2), vs(3,2)
    double precision :: omega, phi, e(3)

    omega = 2.0d0*dpi/period
    phi   = omega*qt + phase0
    e     = [cos(phi), sin(phi), 0.0d0]
    xs(:,1) = -sep*q_mass/(1.0d0 + q_mass) * e
    xs(:,2) =  sep        /(1.0d0 + q_mass) * e
    ! v = omega z-hat x x
    vs(:,1) = omega * [-xs(2,1), xs(1,1), 0.0d0]
    vs(:,2) = omega * [-xs(2,2), xs(1,2), 0.0d0]

  end subroutine star_orbits

  !> Conservative state of the wind of star istar at offset dx from it (r is
  !> |dx|, floored at r_in/2), moving with the star's velocity vstar.
  pure subroutine wind_state(istar, dx, vstar, w)
    ${GPU_ROUTINE_SEQ()}$
    integer, intent(in)           :: istar
    double precision, intent(in)  :: dx(3), vstar(3)
    double precision, intent(out) :: w(1:nw_phys)
    double precision :: r, reff, rho, p, v(3)

    r      = sqrt(sum(dx**2))
    reff   = max(r, 0.5d0*r_in(istar))
    rho    = mdot(istar) / (4.0d0*dpi*reff**2*vw(istar))
    p      = rho * (vw(istar)/mach(istar))**2 / hd_gamma

    v = vstar
    if (r > 0.0d0) v = v + vw(istar) * dx / r

    w = 0.0d0
    w(iw_rho)      = rho
    w(iw_mom(1:3)) = rho * v
    w(iw_e)        = p/(hd_gamma - 1.0d0) + 0.5d0*rho*sum(v**2)
    if (istar == 1) w(tracer(1)) = rho

  end subroutine wind_state

  !> Internal boundary: hold each star's source region at its wind state,
  !> and floor the temperature everywhere else. Called by name from
  !> getintbc, on the device, for every interior cell.
  subroutine usr_internal_bc(level, qt, x, w)
    ${GPU_ROUTINE_SEQ()}$
    integer, intent(in)             :: level
    double precision, intent(in)    :: qt
    double precision, intent(in)    :: x(1:ndim)
    double precision, intent(inout) :: w(1:nw_phys)
    double precision :: xs(3,2), vs(3,2), dx(3), ekin, eint_min
    integer :: istar

    call star_orbits(qt, xs, vs)
    do istar = 1, 2
       dx = x(1:3) - xs(:,istar)
       if (sum(dx**2) <= r_in(istar)**2) then
          call wind_state(istar, dx, vs(:,istar), w)
          return
       end if
    end do

    ! temperature floor
    ekin     = 0.5d0 * sum(w(iw_mom(1:3))**2) / w(iw_rho)
    eint_min = w(iw_rho) * t_floor / (hd_gamma - 1.0d0)
    if (w(iw_e) - ekin < eint_min) w(iw_e) = ekin + eint_min

  end subroutine usr_internal_bc

  !> Fill the domain with whichever wind has the larger ram pressure
  !> rho v^2 ~ mdot vw / r^2 at that point.
  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision :: xs(3,2), vs(3,2), dx1(3), dx2(3), ram1, ram2
    double precision :: wpt(1:nw_phys)
    integer :: ix1, ix2, ix3

    call star_orbits(0.0d0, xs, vs)
    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             dx1 = x(ix1,ix2,ix3,1:3) - xs(:,1)
             dx2 = x(ix1,ix2,ix3,1:3) - xs(:,2)
             ram1 = mdot(1)*vw(1) / max(sum(dx1**2), tiny(1.0d0))
             ram2 = mdot(2)*vw(2) / max(sum(dx2**2), tiny(1.0d0))
             if (ram1 >= ram2) then
                call wind_state(1, dx1, vs(:,1), wpt)
             else
                call wind_state(2, dx2, vs(:,2), wpt)
             end if
             w(ix1,ix2,ix3,1:nw_phys) = wpt
          end do
       end do
    end do

  end subroutine initonegrid_usr

  !> Force the finest level around both stars, cap Lohner at lmax_far elsewhere.
  subroutine usr_refine_grid(level,qt,w,x,refineflag,coarsenflag,norefineflag,nocoarsenflag)
    use mod_global_parameters
    ${GPU_ROUTINE_SEQ()}$
    integer, intent(in)             :: level
    double precision, intent(in)    :: qt
    double precision, intent(in)    :: x(1:3)
    double precision, intent(in)    :: w(1:nw_phys)
    logical, intent(inout) :: refineflag, coarsenflag, norefineflag, nocoarsenflag
    double precision :: xs(3,2), vs(3,2)
    integer :: istar
    logical :: near

    call star_orbits(qt, xs, vs)
    near = .false.
    do istar = 1, 2
       if (sum((x - xs(:,istar))**2) <= (refine_radius*r_in(istar))**2) near = .true.
    end do

    if (near) then
       refineflag    = .true.
       nocoarsenflag = .true.
       coarsenflag   = .false.
       norefineflag  = .false.
    else
       ! leave the decision to Lohner, but never above lmax_far
       coarsenflag   = .false.
       norefineflag  = (level >= lmax_far)
    end if

  end subroutine usr_refine_grid

end module mod_usr
