module mod_usr
  use mod_amrvac
  use mod_physics
  implicit none

  ! --- User parameters (code units) ---
  double precision :: ca, mach, chi, rc, a_int, eps_rho
  double precision :: x1c, x2c, x3c
  !$acc declare create(ca, mach)

  ! --- Turbulence ---
  integer, parameter :: nmodes = 4096
  double precision :: vturb_rms
  double precision, allocatable :: kvec(:,:), phase(:,:), amp(:)
  logical :: turb_initialised = .false.
  double precision, parameter :: pi = 3.1415926535897932384626433832795d0

contains

  subroutine usr_init()
    use mod_global_parameters, only: unit_length, unit_temperature, unit_numberdensity

    ! ------------------------------------------------------------
    ! Set physical unit scales necessary for cooling module
    ! Girichidis+21 baseline:
    !   rc = 50 pc, nw = 1e-3 cm^-3, Tw = 1e6 K, vw = 100 km/s
    !  => code length = 50 pc, code number density = 1e-3 cm^-3, code temperature = 1e6 K
    ! ------------------------------------------------------------
    unit_length        = 50.d0 * 3.0856776d18   ! cm
    unit_numberdensity = 1.d-3                  ! cm^-3
    unit_temperature   = 1.d6                   ! K

    call set_coordinate_system("Cartesian_3D")

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr   ! impose uniform inflow at x_min

    call phys_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr()
    use mod_global_parameters, only: xprobmin1, xprobmax1, xprobmin2, xprobmax2, xprobmin3, xprobmax3
    use mod_physics, only: hd_gamma

    ! ---- Geometry ----
    rc   = 1.0d0                    ! cloud radius
    chi  = 140.0d0                  ! cloud/wind density contrast
    a_int = 3.d0 * 20.d0 / 1280.d0  ! cloud boundary width

    ca = sqrt(hd_gamma)
    mach = (100.d0 * 1.d5 / unit_velocity) / ca    ! so that vwind = mach*ca = 100 km s^-1

    ! Cloud centre (upstream so it has room to develop a tail)
    x1c = xprobmin1 + 0.25d0*(xprobmax1 - xprobmin1)
    x2c = 0.5d0*(xprobmin2 + xprobmax2)
    x3c = 0.5d0*(xprobmin3 + xprobmax3)

    eps_rho = 0.20d0  ! density gradient across cloud

    !$acc update device(ca, mach)

    call init_turbulence_modes()

    if (mype == 0) call print_timescales()

  end subroutine initglobaldata_usr


  subroutine initonegrid_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
                             ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters
    use mod_physics, only: hd_gamma
    implicit none

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
    integer, intent(in)             :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    integer :: i1, i2, i3
    double precision :: rho_w, p0, vwind
    double precision :: dx1, dx2, dx3
    double precision :: r2, r, S
    double precision :: sdiag, clamp_sdiag
    double precision :: rho_loc
    double precision, parameter :: inv_sqrt2 = 0.7071067811865475244d0
    integer :: n
    double precision :: dvx, dvy, dvz, arg

    rho_w = 1.0d0
    vwind = mach * ca
    p0    = rho_w * ca**2 / hd_gamma   ! so that ambient cs = ca

    do i3 = ixOmin3, ixOmax3
      do i2 = ixOmin2, ixOmax2
        do i1 = ixOmin1, ixOmax1

          ! distance to cloud centre
          dx1 = x(i1,i2,i3,1) - x1c
          dx2 = x(i1,i2,i3,2) - x2c
          dx3 = x(i1,i2,i3,3) - x3c

          r2 = dx1*dx1 + dx2*dx2 + dx3*dx3
          r  = sqrt(r2)

          ! smooth indicator (inside ~1, outside ~0)
          S = 0.5d0 * (1.0d0 - tanh((r - rc)/a_int))

          ! diagonal stratification
          sdiag = (( (x(i1,i2,i3,1) - x1c) + (x(i1,i2,i3,2) - x2c) ) * inv_sqrt2) / rc

          ! clamp to [-1, 1]
          if (sdiag < -1.0d0) then
            clamp_sdiag = -1.0d0
          else if (sdiag > 1.0d0) then
            clamp_sdiag =  1.0d0
          else
            clamp_sdiag = sdiag
          end if

          ! density: ambient + cloud * (1 + eps*sdiag)
          rho_loc = rho_w + (rho_w*(chi - 1.0d0)*S) * (1.0d0 + eps_rho*clamp_sdiag)
          w(i1,i2,i3, rho_) = rho_loc

          ! --- turbulent velocity
          dvx = 0.0d0; dvy = 0.0d0; dvz = 0.0d0
          if (S > 1.0d-12 .and. turb_initialised) then
            do n = 1, nmodes
              arg = kvec(1,n)*x(i1,i2,i3,1) + kvec(2,n)*x(i1,i2,i3,2) + kvec(3,n)*x(i1,i2,i3,3)
              dvx = dvx + amp(n) * sin(arg + phase(1,n))
              dvy = dvy + amp(n) * sin(arg + phase(2,n))
              dvz = dvz + amp(n) * sin(arg + phase(3,n))
            end do
          end if

          ! velocity: wind outside, cloud gets turbulence
          w(i1,i2,i3, mom(1)) = vwind * (1.0d0 - S) + S * dvx
          w(i1,i2,i3, mom(2)) = S * dvy
          w(i1,i2,i3, mom(3)) = S * dvz

          ! pressure (primitive slot e_ before conversion)
          w(i1,i2,i3, e_) = p0

        end do
      end do
    end do

    call phys_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
                          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)

  end subroutine initonegrid_usr


  subroutine init_turbulence_modes()
    use mod_global_parameters, only: unit_velocity
    integer :: n, i, seedsize, m
    integer, allocatable :: seed(:)
    double precision :: u1,u2,u3, phi
    double precision :: kmag, k0
    double precision :: khx,khy,khz, norm
    integer, parameter :: mmin = 4
    integer            :: mmax

    if (turb_initialised) return

    if (.not. allocated(kvec))  allocate(kvec(3,nmodes))
    if (.not. allocated(phase)) allocate(phase(3,nmodes))
    if (.not. allocated(amp))   allocate(amp(nmodes))

    ! 1 km/s in code units
    vturb_rms = (1.0d5) / unit_velocity
    mmax = 64

    ! Fixed deterministic seed
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    do i = 1, seedsize
      seed(i) = 123456 + 37*i
    end do
    call random_seed(put=seed)
    deallocate(seed)

    k0 = dble(mmin) * pi / rc   ! reference for scaling

    do n = 1, nmodes
      ! pick integer m in [mmin, mmax]
      call random_number(u1)
      m = mmin + int( (mmax - mmin + 1) * u1 )
      kmag = dble(m) * pi / rc

      ! random direction for k-hat
      call random_number(u1); call random_number(u2)
      ! u1 -> cos(theta) uniform in [-1,1], u2 -> phi uniform in [0,2pi)
      khz = 2.0d0*u1 - 1.0d0
      phi = 2.0d0*pi*u2
      norm = sqrt(max(1.0d-300, 1.0d0 - khz*khz))
      khx = norm*cos(phi)
      khy = norm*sin(phi)

      kvec(1,n) = kmag * khx
      kvec(2,n) = kmag * khy
      kvec(3,n) = kmag * khz

      ! random phases in [0,2pi)
      call random_number(u1); phase(1,n) = 2.0d0*pi*u1
      call random_number(u2); phase(2,n) = 2.0d0*pi*u2
      call random_number(u3); phase(3,n) = 2.0d0*pi*u3

      ! Burgers: P(k) ~ k^-4  =>  amp ~ k^-2
      amp(n) = (kmag/k0)**(-2.0d0)
    end do

    ! Normalise amplitudes
    amp(:) = amp(:) * (vturb_rms / sqrt(1.5d0*sum(amp(:)**2)))

    turb_initialised = .true.
  end subroutine init_turbulence_modes


  subroutine specialbound_usr(qt, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
                              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB, w, x)
    !$acc routine vector
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
    integer, intent(in)             :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    integer, intent(in)             :: iB
    double precision, intent(in)    :: qt
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    double precision :: rho_w, p0, vwind, inv_gamma_m1
    integer          :: ix1, ix2, ix3

    rho_w = 1.0d0
    vwind = mach * ca
    p0    = rho_w * ca**2 / hd_gamma

    select case(iB)

    ! --- Fixed inflow on left x boundary (id=1) ---
    case(1)

      inv_gamma_m1 = 1.0d0/(hd_gamma - 1.0d0)

      !$acc loop collapse(3) vector
      do ix3 = ixOmin3, ixOmax3
        do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

            ! primitives
            w(ix1,ix2,ix3, rho_)   = rho_w
            w(ix1,ix2,ix3, mom(1)) = vwind
            w(ix1,ix2,ix3, mom(2)) = 0.0d0
            w(ix1,ix2,ix3, mom(3)) = 0.0d0
            w(ix1,ix2,ix3, e_)     = p0

            ! inline to_conservative() for performance:
            w(ix1,ix2,ix3, e_) = w(ix1,ix2,ix3, e_) * inv_gamma_m1 + 0.5_dp * w(ix1,ix2,ix3, rho_) * &
                                 sum(w(ix1,ix2,ix3,iw_mom(1:ndim))**2)

            w(ix1,ix2,ix3, iw_mom(1)) = w(ix1,ix2,ix3, iw_rho) * w(ix1,ix2,ix3, iw_mom(1))
            w(ix1,ix2,ix3, iw_mom(2)) = w(ix1,ix2,ix3, iw_rho) * w(ix1,ix2,ix3, iw_mom(2))
            w(ix1,ix2,ix3, iw_mom(3)) = w(ix1,ix2,ix3, iw_rho) * w(ix1,ix2,ix3, iw_mom(3))

          end do
        end do
      end do

    end select

  end subroutine specialbound_usr


  subroutine print_timescales()
    use mod_global_parameters, only: unit_length, unit_time, xprobmin1, xprobmax1
    use mod_physics,          only: hd_gamma
    implicit none

    logical, save :: printed = .false.
    double precision :: vwind, t_adv, t_cc, t_s_w, t_s_c, t_box
    double precision :: myr_in_s, t0_Myr
    double precision :: cs_w, cs_c
    double precision :: mach_w, mach_c
    double precision :: Lx
    double precision :: vwind_kms, csw_kms, csc_kms
    double precision :: pc_in_cm, unit_length_pc

    if (printed) return
    printed = .true.

    vwind = mach * ca

    cs_w  = ca
    cs_c  = ca / sqrt(chi)

    mach_w = vwind / cs_w
    mach_c = vwind / cs_c

    t_adv = rc / vwind
    t_cc  = sqrt(chi) * rc / vwind
    t_s_w = rc / cs_w
    t_s_c = rc / cs_c
    Lx    = xprobmax1 - xprobmin1
    t_box = Lx / vwind

    myr_in_s = 1.d6 * 365.25d0 * 24.d0 * 3600.d0
    t0_Myr   = unit_time / myr_in_s

    vwind_kms = vwind * unit_velocity / 1.d5
    csw_kms   = cs_w  * unit_velocity / 1.d5
    csc_kms   = cs_c  * unit_velocity / 1.d5

    pc_in_cm       = 3.0856776d18
    unit_length_pc = unit_length / pc_in_cm

    if (mype == 0) then
      print *, '================================================='
      print *, ' Cloud-crushing diagnostics (startup)'
      print *, '-------------------------------------------------'
      print *, ' Flow:'
      print *, '   v_w        =', vwind, ' (code)  =', vwind_kms, ' km/s'
      print *, '   cs_w       =', cs_w,  ' (code)  =', csw_kms,  ' km/s'
      print *, '   cs_c       =', cs_c,  ' (code)  =', csc_kms,  ' km/s'
      print *, '   Mach_w     = v_w/cs_w =', mach_w
      print *, '   Mach_c     = v_w/cs_c =', mach_c
      print *, '-------------------------------------------------'
      print *, ' Timescales (code | Myr):'
      print *, '   t_adv  =', t_adv, ' | ', t_adv*t0_Myr, ' Myr'
      print *, '   t_cc   =', t_cc,  ' | ', t_cc *t0_Myr, ' Myr'
      print *, '   t_s,w  =', t_s_w, ' | ', t_s_w*t0_Myr, ' Myr'
      print *, '   t_s,c  =', t_s_c, ' | ', t_s_c*t0_Myr, ' Myr'
      print *, '   t_box  =', t_box, ' | ', t_box*t0_Myr, ' Myr'
      print *, '-------------------------------------------------'
      print *, ' Unit info:'
      print *, '   unit_length      =', unit_length, ' cm (= ', unit_length_pc, ' pc)'
      print *, '   unit_time        =', unit_time,   ' s  (= ', t0_Myr, ' Myr)'
      print *, '   unit_temperature =', unit_temperature, ' K'
      print *, '================================================='
    end if
  end subroutine print_timescales


end module mod_usr

