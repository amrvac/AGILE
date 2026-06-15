module mod_usr

  use mod_amrvac
  use mod_physics
  implicit none

! User parameters.
  double precision :: rho_1,rho_2
  double precision :: p_1,p_2
  double precision :: angle_frac,x_offset
  double precision :: mx_post,my_post,epost,epre
  double precision :: rho_3,xpos_3,zpos_3,rad_3

!$acc declare create(rho_1,rho_2,p_1,p_2,angle_frac,x_offset,mx_post,my_post,epost,epre)
!$acc declare create(rho_3,xpos_3,zpos_3,rad_3)

contains

  subroutine usr_init()
    use mod_global_parameters
    implicit none
    call set_coordinate_system('Cartesian_3D')
    usr_set_parameters => set_parameters_usr
    usr_init_one_grid  => init_one_grid_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    call phys_activate()
  end subroutine usr_init

  subroutine print_params()
    implicit none
    character(len=*), parameter :: fmt = '(A20,1X,SP,ES12.2)'

    print *, '------------------ PARAMETERS ------------------'
    print fmt,  'rho_1',           rho_1
    print fmt,  'rho_2',           rho_2
    print fmt,  'rho_3',           rho_3
    print fmt,  'rad_3',           rad_3
    print fmt,  'xpos_3',          xpos_3
    print fmt,  'zpos_3',          zpos_3
    print fmt,  'p_1',             p_1
    print fmt,  'p_2',             p_2
    print fmt,  'hd_gamma',        hd_gamma
    print fmt,  'offset',          x_offset
    print fmt,  'angle in radians',dpi*angle_frac
    print fmt,  'mx_post',         mx_post
    print fmt,  'my_post',         my_post
    print fmt,  'epre',            epre
    print fmt,  'epost',           epost
    print *, '------------------------------------------------'
  end subroutine

  subroutine set_parameters_usr()
    use mod_global_parameters
    use mod_physics
    implicit none

    hd_gamma=1.4d0
    angle_frac=1.0d0/3.0d0
    x_offset=1.0d0/6.0d0
    xpos_3=7.0d0
    zpos_3=0.5d0
    rad_3=0.25d0
    rho_3=10.d0
    rho_1=hd_gamma
    rho_2=8.0d0
    p_1=1.0d0
    p_2=1.165d2
    mx_post=rho_2*8.25d0*dsin(dpi*angle_frac)
    my_post=-rho_2*8.25d0*dcos(dpi*angle_frac)
    epre=p_1/(hd_gamma-1.0d0)
    epost=p_2/(hd_gamma-1.0d0)+(mx_post**2+my_post**2)/(2.0d0*rho_2)
!$acc update device(hd_gamma,angle_frac,x_offset,rho_1,rho_2,p_1,p_2,mx_post,my_post,epre,epost)
!$acc update device(rad_3,rho_3,xpos_3,zpos_3)

    if (mype == 0) call print_params()
  end subroutine set_parameters_usr

  subroutine init_one_grid_usr(&
      ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3,&
      ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3,&
      w, x)
    use mod_global_parameters
    use mod_physics
    implicit none
    integer, intent(in) :: ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3
    integer, intent(in) :: ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3
    double precision, intent(inout),&
      dimension(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, 1:nw) :: w
    double precision, intent(in),&
      dimension(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, 1:ndim) :: x

    associate(&
      w_ => w(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, :),&
      x_ => x(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, :))

     where(x_(:,:,:,1)>(x_(:,:,:,2)/dtan(dpi*angle_frac)+x_offset))
       w_(:,:,:,iw_rho)    = rho_1
       w_(:,:,:,iw_e)      = epre
       w_(:,:,:,iw_mom(1)) = 0.0d0
       w_(:,:,:,iw_mom(2)) = 0.0d0
       w_(:,:,:,iw_mom(3)) = 0.0d0
     elsewhere
       w_(:,:,:,iw_rho)    = rho_2
       w_(:,:,:,iw_e)      = epost
       w_(:,:,:,iw_mom(1)) = mx_post
       w_(:,:,:,iw_mom(2)) = my_post
       w_(:,:,:,iw_mom(3)) = 0.0d0
     endwhere

     where(((x_(:,:,:,1)-xpos_3)**2+(x_(:,:,:,2))**2+(x_(:,:,:,3)-zpos_3)**2)<rad_3)
       w_(:,:,:,iw_rho)    = rho_3
     endwhere
    end associate

  end subroutine init_one_grid_usr

! NOTE: This must be named exactly `specialbound_usr` in AGILE since during
!   compilation it deals with this function.
! NOTE: This is ran on the GPU.
  subroutine specialbound_usr(&
      qt,&
      ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3,&
      ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3,&
      iB, w, x)
!$acc routine vector
    ! use mod_physics, only: to_conservative, to_primitive
    use mod_global_parameters
    use mod_physics
    implicit none
    integer, intent(in) :: ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3
    integer, intent(in) :: ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3
    integer, intent(in) :: iB
    double precision, intent(in) :: qt
    double precision, intent(inout),&
      dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:nw_phys) :: w
    double precision, intent(in),&
      dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:ndim) :: x
    integer :: ix1, ix2, ix3

    select case(iB)
    case(1)
      ! implementation of fixed postshock state at left boundary
!$acc loop collapse(3) vector
      do ix3 = ixOmin3, ixOmax3
      do ix2 = ixOmin2, ixOmax2
      do ix1 = ixOmin1, ixOmax1
        w(ix1,ix2,ix3,iw_rho)    = rho_2
        w(ix1,ix2,ix3,iw_mom(1)) = mx_post
        w(ix1,ix2,ix3,iw_mom(2)) = my_post
        w(ix1,ix2,ix3,iw_mom(3)) = 0.0d0
        w(ix1,ix2,ix3,iw_e)      = epost
      end do
      end do
      end do
    case(3)
      ! implementation of bottom boundary: fixed before x<1/6, solid wall x>=1/6
!$acc loop collapse(3) vector
      do ix3 = ixOmin3, ixOmax3
      do ix2 = ixOmin2, ixOmax2
      do ix1 = ixOmin1, ixOmax1
        if(x(ix1,ix2,ix3,1)<=x_offset)then
           w(ix1,ix2,ix3,iw_rho)    = rho_2
           w(ix1,ix2,ix3,iw_mom(1)) = mx_post
           w(ix1,ix2,ix3,iw_mom(2)) = my_post
           w(ix1,ix2,ix3,iw_mom(3)) = 0.0d0
           w(ix1,ix2,ix3,iw_e)      = epost
        else
           w(ix1,ix2,ix3,iw_rho)    =  w(ix1,2*ixOmax2-ix2+1,ix3,iw_rho)
           w(ix1,ix2,ix3,iw_mom(1)) =  w(ix1,2*ixOmax2-ix2+1,ix3,iw_mom(1))
           w(ix1,ix2,ix3,iw_mom(2)) = -w(ix1,2*ixOmax2-ix2+1,ix3,iw_mom(2))
           w(ix1,ix2,ix3,iw_mom(3)) =  w(ix1,2*ixOmax2-ix2+1,ix3,iw_mom(3))
           w(ix1,ix2,ix3,iw_e)      =  w(ix1,2*ixOmax2-ix2+1,ix3,iw_e)
        endif
      end do
      end do
      end do
    case(4)
      ! implementation of top: pre and post shock states, time dependent
!$acc loop collapse(3) vector
      do ix3 = ixOmin3, ixOmax3
      do ix2 = ixOmin2, ixOmax2
      do ix1 = ixOmin1, ixOmax1
        if(x(ix1,ix2,ix3,1)>(1.0d1*qt/dsin(dpi*angle_frac)+x(ix1,ix2,ix3,2)/dtan(dpi*angle_frac)+x_offset))then
           w(ix1,ix2,ix3,iw_rho)    = rho_1
           w(ix1,ix2,ix3,iw_mom(1)) = 0.0d0
           w(ix1,ix2,ix3,iw_mom(2)) = 0.0d0
           w(ix1,ix2,ix3,iw_mom(3)) = 0.0d0
           w(ix1,ix2,ix3,iw_e)      = epre
        else
           w(ix1,ix2,ix3,iw_rho)    = rho_2
           w(ix1,ix2,ix3,iw_mom(1)) = mx_post
           w(ix1,ix2,ix3,iw_mom(2)) = my_post
           w(ix1,ix2,ix3,iw_mom(3)) = 0.0d0
           w(ix1,ix2,ix3,iw_e)      = epost
        endif
      end do
      end do
      end do
    end select

  end subroutine specialbound_usr


   subroutine specialvar_output( &
      ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3,&
      ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3,&
      w, x, normconv)
   ! this subroutine can be used in convert, to add auxiliary variables to the
   ! converted output file, for further analysis using tecplot, paraview, ....
   ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
   !
   ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
   ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_physics
    implicit none
    integer, intent(in) :: ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3
    integer, intent(in) :: ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3
    double precision,&
      dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:nw+nwauxio) :: w
    double precision, intent(in),&
      dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:ndim) :: x
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: gradrho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision                   :: drho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision                   :: kk,kk0,grhomax,kk1
    integer                            :: idims

    gradrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    do idims=1,ndim
        call gradient(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw_rho),&
                    ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3,&
                    ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3,&
                    idims,drho)
       gradrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= &
          gradrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)+drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2.0d0
    enddo
    gradrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= &
        dsqrt(gradrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    kk=5.0d0
    kk0=0.01d0
    kk1=1.0d0
    grhomax=20000.0d0

    ! putting the schlierplot of density in nwauxio=1
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1)= &
       dexp(-kk*(gradrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
   end subroutine specialvar_output

   subroutine specialvarnames_output(varnames)

   ! newly added variables
   character(len=*) :: varnames
   varnames='schlierho'

   end subroutine specialvarnames_output


end module mod_usr
                                    
