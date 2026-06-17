module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

contains

  subroutine usr_init()

    call set_coordinate_system("Cartesian_3D")

    usr_init_one_grid => initonegrid_usr

    call phys_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,rho_)=1.0d0
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(1))=0.5d0
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(2))=0.4d0
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(3))=0.6d0
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,p_)=1.0d0

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

end module mod_usr
