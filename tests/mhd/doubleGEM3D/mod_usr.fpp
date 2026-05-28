module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision, parameter :: pi = 4d0*atan(1d0)

  double precision :: sheetl, rhorat, BB0, llx, lly, psi0bot, psi0top
  double precision :: llz,xmid,ymid,ysh1,ysh2,fkx,fky

  double precision, dimension(3) :: L
  double precision, dimension(2:3) :: v_per
  double precision, dimension(2:3) :: f_til
  double precision :: zmid,sig_z
  integer :: mode_root
  integer :: n_modes
  integer :: seed

! Precomputed random numbers.
! (4, n_modes, n_modes, 2:3)
  double precision, dimension(:,:,:,:), allocatable :: rand_1

!$acc declare create(sheetl,rhorat,BB0,llx,lly,psi0bot,psi0top)
!$acc declare create(llz,xmid,ymid,ysh1,ysh2,fkx,fky)
!$acc declare create(zmid, sig_z, mode_root, n_modes)
!$acc declare create(L, v_per, f_til)
!$acc declare create(rand_1)

contains

  subroutine usr_init()
    use mod_global_parameters
    implicit none
    call set_coordinate_system("Cartesian_3D")
    call params_read_usr(par_files)
    usr_set_parameters => set_parameters_usr
    usr_init_one_grid => initonegrid_usr

    call phys_activate()

  end subroutine usr_init

  subroutine params_read_usr(files)
    implicit none
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ seed, v_per, sig_z,&
      f_til, mode_root, n_modes, psi0bot, psi0top
    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status='old')
      read(unitpar, usr_list, end=111)
111   close(unitpar) 
    end do

!$acc update device(v_per,sig_z,f_til,mode_root,n_modes)
!$acc update device(psi0bot,psi0top)

  end subroutine params_read_usr

  subroutine set_parameters_usr()
    use mod_global_parameters
    use mod_physics
    use mod_random
    implicit none
    integer :: i1, i2, i3, i4
    integer, parameter :: i8 = selected_int_kind(18) ! as in mod_random
    type(rng_t) :: rng

    sheetl=1.0_dp
    rhorat=0.1_dp
    BB0=1.0_dp
    llx=xprobmax1-xprobmin1
    lly=xprobmax2-xprobmin2
    llz=xprobmax3-xprobmin3
    xmid=xprobmin1+0.5d0*llx
    ymid=xprobmin2+0.5d0*lly
    zmid=xprobmin3+0.5d0*llz
    ysh1=xprobmin2+0.25d0*lly
    ysh2=xprobmin2+0.75d0*lly
    fkx=2.0_dp*pi/llx
    fky=2.0_dp*pi/lly

    L = [llx,lly,llz]

    allocate(rand_1(4,n_modes,n_modes,2:3))

    if (mype == 0) then
      call rng%set_seed([int(seed, i8), 123456789_i8])
      do i4 = 2, 3
      do i3 = 1, n_modes
      do i2 = 1, n_modes
      do i1 = 1, 4
        rand_1(i1,i2,i3,i4) = rng%normal()
      end do
      end do
      end do
      end do
    end if
    if (npe > 0) then
      call MPI_BCAST(rand_1, size(rand_1), MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)
    end if

!$acc update device(sheetl,rhorat,BB0,llx,lly,llz,xmid,ymid,zmid,ysh1,ysh2,fkx,fky)
!$acc update device(L)
!$acc update device(rand_1)

    if (mype == 0) call print_params()
  end subroutine set_parameters_usr

  subroutine print_params()
    implicit none
    character(len=*), parameter :: fmt = '(A20,1X,SP,ES12.2)'
    character(len=*), parameter :: fmti = '(A20,1X,I12)'

    write(*,*) '=================DOUBLE GEM================'
    write(*,*) 'mhd_eta   =',mhd_eta
    write(*,*) 'sheetl    =',sheetl
    write(*,*) 'rhorat    =',rhorat
    write(*,*) 'BB0       =',BB0
    write(*,*) 'llx       =',llx
    write(*,*) 'lly       =',lly
    write(*,*) 'psi0bot   =',psi0bot
    write(*,*) 'psi0top   =',psi0top
    write(*,*) '==========================================='
    print *, '------------------ PARAMETERS ------------------'
    print fmti, 'seed',         seed
    print fmti, 'mode_root',    mode_root
    print fmti, 'n_modes',      n_modes
    print fmt,  'v_per(2)',     v_per(2)
    print fmt,  'v_per(3)',     v_per(3)
    print fmt,  'zmid',         zmid
    print fmt,  'sig_z',        sig_z
    print fmt,  'f_til(2)',     f_til(2)
    print fmt,  'f_til(3)',     f_til(3)
    print fmt,  'L(1)',         L(1)
    print fmt,  'L(2)',         L(2)
    print fmt,  'L(3)',         L(3)
    print *, '------------------------------------------------'
    ! NOTE: See if these are roughly 0 and 1, verifies random number generation
    print fmt, 'normal mean', sum(rand_1)/(8*n_modes**2)
    print fmt, 'normal var ', sum(rand_1**2)/(8*n_modes**2)-(sum(rand_1)/(8*n_modes**2))**2
    print *, '------------------------------------------------'
  end subroutine print_params

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    integer :: ix, iy, nx, ny

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(1)) = 0.0_dp
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(2)) = 0.0_dp
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(3)) = 0.0_dp

    associate(&
      w_ => w(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, :),&
      x_ => x(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, :))
      w_(:,:,:,iw_mom(2)) = 0
      w_(:,:,:,iw_mom(3)) = 0
      do iy = 1, n_modes-1
      do ix = 1, n_modes-1
        nx = mode_root+ix
        ny = mode_root+iy
        w_(:,:,:,iw_mom(2)) = w_(:,:,:,iw_mom(2))+&
          (rand_1(1,ix,iy,2)*sin(2*pi*x_(:,:,:,1)/L(1)*nx)+&
           rand_1(2,ix,iy,2)*cos(2*pi*x_(:,:,:,1)/L(1)*nx))*&
          (rand_1(3,ix,iy,2)*sin(2*pi*x_(:,:,:,2)/L(2)*ny)+&
           rand_1(4,ix,iy,2)*cos(2*pi*x_(:,:,:,2)/L(2)*ny))
      end do
      end do
      do iy = 1, n_modes-1
      do ix = 1, n_modes-1
        nx = mode_root+ix
        ny = mode_root+iy
        w_(:,:,:,iw_mom(3)) = w_(:,:,:,iw_mom(3))+&
          (rand_1(1,ix,iy,3)*sin(2*pi*x_(:,:,:,1)/L(1)*nx)+&
           rand_1(2,ix,iy,3)*cos(2*pi*x_(:,:,:,1)/L(1)*nx))*&
          (rand_1(3,ix,iy,3)*sin(2*pi*x_(:,:,:,2)/L(2)*ny)+&
           rand_1(4,ix,iy,3)*cos(2*pi*x_(:,:,:,2)/L(2)*ny))
      end do
      end do
      w_(:,:,:,iw_mom(2)) = w_(:,:,:,iw_mom(2))*f_til(2)/(n_modes-1)
      w_(:,:,:,iw_mom(3)) = w_(:,:,:,iw_mom(3))*f_til(3)/(n_modes-1)

      w_(:,:,:,iw_mom(2)) = v_per(2)*(&
        sin(2*pi*x_(:,:,:,1)/L(1)*mode_root)/2+&
        sin(2*pi*x_(:,:,:,2)/L(2)*mode_root)/2+&
        w_(:,:,:,iw_mom(2)))*&
        exp(-((x_(:,:,:,3)-zmid)/sig_z)**2)
      w_(:,:,:,iw_mom(3)) = v_per(3)*(&
        sin(2*pi*x_(:,:,:,1)/L(1)*mode_root)/2+&
        sin(2*pi*x_(:,:,:,2)/L(2)*mode_root)/2+&
        w_(:,:,:,iw_mom(3)))*&
        exp(-((x_(:,:,:,3)-zmid)/sig_z)**2)

    end associate

   
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(1))= &
       BB0*(-1.0_dp+dtanh((x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh1)/sheetl) &
                   +dtanh((ysh2-x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))/sheetl)) &
       -psi0bot*fky*dcos(fkx*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid)) &
                  *(dsin(fky*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh1)) &
                         +two*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh1)* &
                          dcos(fky*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh1))) &
                    *dexp(-fkx*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid)**2 &
                          -fky*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh1)**2) &
       +psi0top*fky*dcos(fkx*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid))              &
                   *(dsin(fky*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh2)) &
                         +two*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh2)* &
                          dcos(fky*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh2))) &
                    *dexp(-fkx*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid)**2 &
                          -fky*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh2)**2)

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(2))= &
       +psi0bot*fkx*dcos(fky*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh1))              &
                  *(dsin(fkx*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid)) &
                         +two*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid)* &
                         dcos(fkx*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid))) &
                    *dexp(-fkx*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid)**2 &
                          -fky*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh1)**2) &
       -psi0top*fkx*dcos(fky*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh2))              &
                  *(dsin(fkx*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid))&
                        +two*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid)* &
                         dcos(fkx*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid))) &
                    *dexp(-fkx*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)-xmid)**2 &
                          -fky*(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh2)**2)

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(3)) = zero

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_) = (0.5_dp*BB0**2) &
      *(rhorat+1.0_dp/(dcosh((x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh1)/sheetl)**2)  &
              +1.0_dp/(dcosh((x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh2)/sheetl)**2))
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_) = &
       (rhorat+1.0_dp/(dcosh((x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh1)/sheetl)**2)  &
              +1.0_dp/(dcosh((x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)-ysh2)/sheetl)**2))

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

end module mod_usr
