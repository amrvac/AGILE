module mod_usr
  use mod_amrvac
  use mod_physics
  implicit none

  ! Zhou et al. (2025, ApJ 978:72), Sect. 2.3: oblique ring conduction test
  double precision, parameter :: ring_rmin    = 0.5d0
  double precision, parameter :: ring_rmax    = 0.7d0
  double precision, parameter :: ring_slab_hw = 0.2d0   ! half-width in ring-frame z
  double precision, parameter :: T_hot        = 12.0d0  ! pressure = T since rho=1
  double precision, parameter :: T_cold       = 10.0d0
  ! Constant parallel conductivity kappa*T^(5/2) = 0.01 (Zhou et al. 2025, Sect. 2.3)
  double precision, parameter :: ring_kappa0  = 0.01d0

contains

  subroutine usr_init()
    use mod_global_parameters
    implicit none

    ! kappa_0 * T^(5/2) = ring_kappa0 = const
    hypertc_kappa0      = ring_kappa0
    hypertc_const_kappa = .true.

    nwauxio = 1

    call set_coordinate_system("Cartesian_3D")
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => extra_var_output
    usr_add_aux_names => extra_var_names_output
    call phys_activate()
  end subroutine usr_init

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    use mod_global_parameters
    use mod_physics
    implicit none
    integer, intent(in)          :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    double precision :: r(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3)
    double precision :: theta(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3)

    ! --- Step 1: cylindrical coords ---
    r(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = &
       dsqrt(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)**2 &
           + x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)**2)

    ! theta = arctan(y/x)
    theta(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = &
       atan2(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2), &
             x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1))

    ! --- Step 2: thermodynamic state ---
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)   = 1.0d0
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(1)) = 0.0d0
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(2)) = 0.0d0
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(3)) = 0.0d0

    ! Hot arc on the left side of the ring, finite slab in z
    where (abs(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,3)) < ring_slab_hw &
           .and. r(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) > ring_rmin &
           .and. r(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) < ring_rmax &
           .and. theta(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) > 11.0d0/12.0d0*dpi &
           .and. theta(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) < 13.0d0/12.0d0*dpi)
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_) = T_hot
    elsewhere
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_) = T_cold
    end where

    ! --- Step 3: azimuthal unit B field: B = (-y, x, 0)/r ---
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(1)) = &
       -x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2) &
       / max(r(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3), 1.0d-15)
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(2)) = &
        x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1) &
       / max(r(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3), 1.0d-15)
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(3)) = 0.0d0

    ! --- Step 4: zero-initialise auxiliary scalars ---
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,psi_) = 0.0d0
#:if defined('HYPERTC') or defined('HYPERTC_ANISO')
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,q_)     = 0.0d0
#:endif
#:if defined('HYPERTC_ANISO')
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,qperp_) = 0.0d0
#:endif

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
  end subroutine initonegrid_usr

  subroutine extra_var_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
    use mod_global_parameters
    implicit none
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw+nwauxio)
    double precision             :: normconv(0:nw+nwauxio)

    ! T = p/rho = (gamma-1)*(e_tot - 0.5*|m|^2/rho - 0.5*B^2) / rho  (conserved input)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1) = &
       (mhd_gamma - 1.0d0) * ( &
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_) &
        - 0.5d0 * ( w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))**2 &
                  + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))**2 &
                  + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))**2 ) &
          / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) &
        - 0.5d0 * ( w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))**2 &
                  + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))**2 &
                  + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))**2 ) &
       ) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
  end subroutine extra_var_output

  subroutine extra_var_names_output(varnames)
    use mod_global_parameters
    implicit none
    character(len=*) :: varnames
    varnames = 'T'
  end subroutine extra_var_names_output

end module mod_usr
