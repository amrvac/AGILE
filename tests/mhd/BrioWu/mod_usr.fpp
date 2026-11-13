!> Brio & Wu (1988) MHD shock tube, run along x on a thin 3D box.
!> gamma = 2 is set in &mhd_list. The left and right states meet at
!> x = 0.5; at t = 0.1 the solution holds, from left to right, a fast
!> rarefaction, a slow compound wave, the contact, a slow shock and a fast
!> rarefaction, which puts every wave family of the flux scheme under test.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

contains

  subroutine usr_init()

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

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,:) = zero

    where (x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1) < 0.5d0)
       w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)   = 1.0d0
       w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_)     = 1.0d0
       w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(2)) = 1.0d0
    elsewhere
       w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)   = 0.125d0
       w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_)     = 0.1d0
       w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(2)) = -1.0d0
    end where
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(1)) = 0.75d0

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

end module mod_usr
