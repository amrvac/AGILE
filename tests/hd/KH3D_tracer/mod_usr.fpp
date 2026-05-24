module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision :: rho_C = 1d1
  double precision :: rho_H = 1d0
  double precision :: kx    = 4d0*acos(-1d0)
  double precision :: kz    = 32d0*atan(1d0)
  double precision :: pint  = 2.5d0
  double precision :: dlin  = 2.5d-2
  double precision :: sig   = 1.25d-3

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_3D")
    usr_init_one_grid => initonegrid_usr
    call phys_activate()
  end subroutine usr_init

  subroutine initonegrid_usr(&
      ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3,&
      ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3,&
      w, x)
    integer, intent(in)             ::&
      ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3,&
      ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3
    double precision, intent(in)    :: x(&
      ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(&
      ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw)
    double precision :: vextra

    select case(iprob)
    case(1)
      vextra = 0d0
    case(2)
      vextra = 1d1
    case(3)
      vextra = 1d2
    case default
      error stop "Unknown iprob"
    endselect

    associate(&
      w_ => w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,:),&
      x_ => x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,:))

      where(x_(:,:,:,2) < 0.25d0 .or. x_(:,:,:,2) > 0.75d0)
        w_(:,:,:,rho_) = rho_H
      elsewhere
        w_(:,:,:,rho_) = rho_C
      endwhere
      w_(:,:,:,p_) = pint
      where((x_(:,:,:,2) <= 0.25d0-dlin) .or. (x_(:,:,:,2) >= 0.75d0+dlin))
        w_(:,:,:,mom(1)) = vextra-0.5d0
      endwhere
      where((x_(:,:,:,2) >= 0.25d0+dlin) .and. (x_(:,:,:,2) <= 0.75d0-dlin))
        w_(:,:,:,mom(1)) = vextra+0.5d0
      endwhere
      where((x_(:,:,:,2) > 0.25d0-dlin) .and. (x_(:,:,:,2) < 0.25d0+dlin))
        w_(:,:,:,mom(1)) = (x_(:,:,:,2)-(0.25d0-dlin))*&
          ((vextra+0.5d0)-(vextra-0.5d0))/(2.0d0*dlin)+(vextra-0.5d0)
      endwhere
      where((x_(:,:,:,2) > 0.75d0-dlin) .and. (x_(:,:,:,2) < 0.75d0+dlin))
        w_(:,:,:,mom(1)) = (x_(:,:,:,2)-(0.75d0-dlin))*&
          ((vextra-0.5d0)-(vextra+0.5d0))/(2.0d0*dlin)+(vextra+0.5d0)
      endwhere
      w_(:,:,:,mom(2)) = 0.01d0*dsin(kx*x_(:,:,:,1))*&
        (dexp(-0.5d0*(x_(:,:,:,2)-0.25d0)**2/sig)+&
         dexp(-0.5d0*(x_(:,:,:,2)-0.75d0)**2/sig))
      w_(:,:,:,mom(3)) = 0.1d0*dsin(kz*x_(:,:,:,3))*&
        (dexp(-0.5d0*(x_(:,:,:,2)-0.25d0)**2/sig)+&
         dexp(-0.5d0*(x_(:,:,:,2)-0.75d0)**2/sig))
      w_(:,:,:,tracer(1)) = 0d0
! NOTE: Tracers are always in conservative form.
      where((x_(:,:,:,2) >= 0.25d0) .and. (x_(:,:,:,2) <= 0.75d0))
        w_(:,:,:,tracer(1)) = 1d0*w_(:,:,:,rho_)
      endwhere

    end associate

    call phys_to_conserved(&
      ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3,&
      ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3,&
      w, x)
  end subroutine initonegrid_usr

end module mod_usr
