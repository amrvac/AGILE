module mod_usr
  use mod_amrvac
  use mod_physics
  implicit none

  ! --- User parameters (code units) ---
  double precision :: ca, mach, chi, rc, a_int
  double precision :: x1c, x2c, x3c
  !$acc declare create(ca, mach, chi, rc, a_int, x1c, x2c, x3c)

contains

  subroutine usr_init()

    call set_coordinate_system("Cartesian_3D")

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr   ! impose uniform inflow at x_min

    call phys_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr()
    use mod_global_parameters, only: xprobmin1, xprobmax1, xprobmin2, xprobmax2, xprobmin3, xprobmax3

    ! Ambient sound speed and wind Mach number
    ca   = 1.0d0
    mach = 2.0d0

    ! Cloud / wind density contrast
    chi  = 100.0d0

    ! Cloud radius
    rc   = 1.0d0

    ! Interface thickness (small but nonzero)
    a_int = 0.05d0 * rc

    ! Cloud centre (a bit upstream so it has room to develop a tail)
    x1c = xprobmin1 + 0.25d0*(xprobmax1 - xprobmin1)
    x2c = 0.5d0*(xprobmin2 + xprobmax2)
    x3c = 0.5d0*(xprobmin3 + xprobmax3)

    !$acc update device(ca, mach, chi, rc, a_int, x1c, x2c, x3c)

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

    ! Density: ambient + overdense cloud
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, rho_) = &
     rho_w * (1.0d0 + (chi - 1.0d0) * S(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))

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

end module mod_usr

