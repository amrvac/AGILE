module mod_usr
  use mod_amrvac
  use mod_physics
  implicit none

  ! Hotplate conduction test from Navarro et al (2022, 2016)

  double precision, parameter :: hp_p0     = 0.1d0    ! uniform initial pressure
  double precision, parameter :: rho_hot   = 0.01d0   ! strip: T = p0/rho_hot = 10
  double precision, parameter :: rho_bg    = 0.1d0    ! ambient: T = p0/rho_bg = 1
  double precision, parameter :: hp_xhalf  = 0.1d0    ! hot-strip half-width in x
  double precision, parameter :: hp_theta  = dpi / 4.0d0  ! B tilt from y-axis (45 deg)
  double precision, parameter :: hp_kappa0 = 0.01d0   ! kappa_par * T^(5/2) = const

contains

  subroutine usr_init()
    use mod_global_parameters
    implicit none

    hypertc_kappa0      = hp_kappa0
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
    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3, &
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2, &
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, &
       ixGmin3:ixGmax3,1:nw)

    call set_hp_prim(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3, &
                     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3, &
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
  end subroutine initonegrid_usr

  ! Pin the y-minimum boundary to the hotplate profile
  ! v=0 always at the plate, so e_tot = p/(gamma-1) + 0.5*|B|^2
  subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB,w,x)
    !$acc routine vector
    use mod_global_parameters
    implicit none
    double precision, intent(in)    :: qt
    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3, &
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2, &
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, &
       ixGmin3:ixGmax3,1:nw)
    ! B = (sin theta, cos theta, 0), |B|^2 = 1
    double precision :: e_plate

    select case (iB)
    case (3)  ! y-min: hotplate
      e_plate = hp_p0 / (mhd_gamma - 1.0d0) + 0.5d0

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)   = rho_bg
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1)) = 0.0d0
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2)) = 0.0d0
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3)) = 0.0d0
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)     = e_plate

      where (abs(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)) < hp_xhalf)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) = rho_hot
      end where

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1)) = dsin(hp_theta)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2)) = dcos(hp_theta)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3)) = 0.0d0
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_)   = 0.0d0
#:if defined('HYPERTC') or defined('HYPERTC_ANISO')
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,q_)     = 0.0d0
#:endif
#:if defined('HYPERTC_ANISO')
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,qperp_) = 0.0d0
#:endif
    end select
  end subroutine specialbound_usr

  ! Set ICs
  subroutine set_hp_prim(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3, &
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    use mod_global_parameters
    implicit none
    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3, &
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2, &
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, &
       ixGmin3:ixGmax3,1:nw)

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)   = rho_bg
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(1)) = 0.0d0
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(2)) = 0.0d0
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(3)) = 0.0d0
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_)     = hp_p0

    where (abs(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)) < hp_xhalf &
           .and. x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2) <= xprobmin2)
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_) = rho_hot
    end where

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(1)) = dsin(hp_theta)
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(2)) = dcos(hp_theta)
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(3)) = 0.0d0
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,psi_)   = 0.0d0
#:if defined('HYPERTC') or defined('HYPERTC_ANISO')
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,q_)     = 0.0d0
#:endif
#:if defined('HYPERTC_ANISO')
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,qperp_) = 0.0d0
#:endif
  end subroutine set_hp_prim

  ! Auxiliary output: temperature T = p/rho recovered from conserved variables
  subroutine extra_var_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
    use mod_global_parameters
    implicit none
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, &
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2, &
       ixImin3:ixImax3,1:ndim)
    double precision              :: w(ixImin1:ixImax1,ixImin2:ixImax2, &
       ixImin3:ixImax3,nw+nwauxio)
    double precision              :: normconv(0:nw+nwauxio)

    ! T = (gamma-1)(e - |m|^2/(2*rho) - B^2/2) / rho
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
