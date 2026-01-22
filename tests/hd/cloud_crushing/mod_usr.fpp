module mod_usr
  use mod_amrvac
  use mod_physics
  implicit none

  ! --- User parameters (code units) ---
  double precision :: ca, mach, chi, rc, a_int, eps_rho
  double precision :: x1c, x2c, x3c

  !$acc declare create(ca, mach)

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
    rc   = 1.0d0         ! cloud radius
    chi  = 140.0d0       ! cloud/wind density contrast
    a_int = 0.05d0 * rc  ! cloud boundary width

    ca = sqrt(hd_gamma)
    mach = (100.d0 * 1.d5 / unit_velocity) / ca    ! so that vwind = mach*ca = 100 km s^-1

    ! Cloud centre (upstream so it has room to develop a tail)
    x1c = xprobmin1 + 0.25d0*(xprobmax1 - xprobmin1)
    x2c = 0.5d0*(xprobmin2 + xprobmax2)
    x3c = 0.5d0*(xprobmin3 + xprobmax3)

    eps_rho = 0.05d0  ! 5% density gradient across cloud

    !$acc update device(ca, mach)

    if (mype == 0) call print_timescales()

  end subroutine initglobaldata_usr


  subroutine initonegrid_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
                             ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
    integer, intent(in)             :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    double precision :: rho_w, p0, vwind
    double precision :: r2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: r (ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: S (ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, parameter :: inv_sqrt2 = 0.7071067811865475244d0

    rho_w = 1.0d0
    vwind = mach * ca
    p0    = rho_w * ca**2 / hd_gamma   ! so that ambient cs = ca

    ! Distance to cloud centre
    r2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
         (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) - x1c)**2 + &
         (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) - x2c)**2 + &
         (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) - x3c)**2

    r(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = sqrt(r2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))

    ! Smooth indicator: ~1 inside, 0 outside
    S(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
      0.5d0 * (1.0d0 - tanh((r(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - rc)/a_int))

    ! Density: ambient + cloud excess * (1 + eps * clamp(sdiag))
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, rho_) = rho_w + &
      (rho_w * (chi - 1.d0) * S(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) * &
      (1.d0 + eps_rho * max(-1.d0, min(1.d0, &
        (((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) - x1c) + &
            (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) - x3c)) * inv_sqrt2) / rc )) )

    ! Velocity: wind outside, cloud initially at rest
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(1)) = & 
     vwind * (1.0d0 - S(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(2)) = 0.0d0
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(3)) = 0.0d0

    ! Pressure (primitive slot e_ before conversion)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, e_) = p0

    call phys_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
                           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)

  end subroutine initonegrid_usr


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

