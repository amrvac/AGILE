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
    usr_print_log     => log_usr

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


  subroutine log_usr()
    use mod_timing, only: itTimeLast, timeLast
    use mod_forest, only: nleafs_active, nleafs_level
    use mod_functions_bfield, only: get_current
    use mod_global_parameters

    logical, save       :: opened = .false.
    integer             :: iigrid, igrid, amode, istatus(MPI_STATUS_SIZE)
    integer, parameter  :: my_unit = 20
    integer             :: level, i, idirmin
    character(len=80)   :: filename
    character(len=40)   :: fmt_string
    character(len=1024) :: line

    double precision :: dvolume, domain_volume
    double precision :: local_ME, local_IE, local_heating
    double precision :: local_maxJ, local_maxv
    double precision :: local_sums(3), global_sums(3)
    double precision :: local_maxs(2), global_maxs(2)
    double precision :: current&
      (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,7-2*ndir:3)

    integer          :: nx1, nx2, nx3, nc, ncells, dit
    double precision :: time, dtTimeLast, cellupdatesPerSecond
    double precision :: wctPerCodeTime, timeToFinish

    local_ME      = 0d0
    local_IE      = 0d0
    local_heating = 0d0
    local_maxJ    = 0d0
    local_maxv    = 0d0
    domain_volume = (xprobmax1-xprobmin1)*&
                    (xprobmax2-xprobmin2)*&
                    (xprobmax3-xprobmin3)

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       block => ps(igrid)
       dvolume = ps(igrid)%dvolume(ixMlo1, ixMlo2, ixMlo3)
       call get_current(ps(igrid)%w, ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3, idirmin, current)
       call compute_block_stats(ps(igrid)%w, current, dvolume,&
          local_ME, local_IE, local_heating, local_maxJ, local_maxv)
    end do

    local_sums = (/ local_ME, local_IE, local_heating /)
    local_maxs = (/ local_maxJ, local_maxv /)
    call MPI_ALLREDUCE(local_sums, global_sums, size(local_sums),&
       MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(local_maxs, global_maxs, size(local_maxs),&
       MPI_DOUBLE_PRECISION, MPI_MAX, icomm, ierrmpi)

    if (mype == 0) then

! average cell updates / rank / second
       nx1 = ixMhi1-ixMlo1+1
       nx2 = ixMhi2-ixMlo2+1
       nx3 = ixMhi3-ixMlo3+1
       nc     = nx1*nx2*nx3 ! per block
       ncells = nc*nleafs_active
       time       = MPI_WTIME()
       dit        = it-itTimeLast
       dtTimeLast = time-timeLast
       itTimeLast = it
       timeLast   = time
       cellupdatesPerSecond = dble(ncells)*dble(nstep)*dble(dit)/&
          (max(dtTimeLast, epsilon(1.0d0))*dble(npe))

! time to finish in hours
       wctPerCodeTime = dtTimeLast / max(dble(dit) * dt, epsilon(1.0d0))
       timeToFinish   = (time_max - global_time) * wctPerCodeTime / 3600.0d0

       filename = trim(base_filename)//".log"

       if (.not. opened) then
          if (restart_from_file == undefined) then
             open(unit=my_unit, file=trim(filename), status='replace')
             close(my_unit, status='delete')
          end if

          amode = ior(MPI_MODE_CREATE, MPI_MODE_WRONLY)
          amode = ior(amode, MPI_MODE_APPEND)
          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, MPI_INFO_NULL,&
             log_fh, ierrmpi)
          opened = .true.

          if (restart_from_file == undefined .or. reset_time) then
             line = 'it global_time dt'
             do level = 1, refine_max_level
                i = len_trim(line)+2
                write(line(i:), '(a,i0)') 'n', level
             end do
             line = trim(line)//' ME IE maxJ maxv heating'//&
                " 'cell updates/s/rank' 'time to finish [hrs]'"
             call MPI_FILE_WRITE(log_fh, trim(line)//new_line('a'),&
                len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
          end if
       end if

       write(line, '(i8,2ES13.4)') it, global_time, dt
       i = len_trim(line)+2
       write(fmt_string, '(a,i0,a)') '(', refine_max_level, 'i10)'
       write(line(i:), fmt_string) nleafs_level(1:refine_max_level)
       i = len_trim(line)+2
       write(line(i:), '(7ES13.4)')&
          global_sums(1), global_sums(2),& ! ME, IE
          global_maxs(1), global_maxs(2),& ! maxJ, maxv
          global_sums(3)/domain_volume,&   ! heating
          cellupdatesPerSecond, timeToFinish

       call MPI_FILE_WRITE(log_fh, trim(line)//new_line('a'),&
          len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

  contains

    subroutine compute_block_stats(w, current, dvolume,&
        bME, bIE, bheating, bmaxJ, bmaxv)
      double precision, intent(in)    :: w&
        (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)
      double precision, intent(in)    :: current&
        (ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,7-2*ndir:3)
      double precision, intent(in)    :: dvolume
      double precision, intent(inout) :: bME, bIE, bheating, bmaxJ, bmaxv
      integer :: i1, i2, i3
      double precision :: v2, B2, pth, J2

      do i3 = ixMlo3, ixMhi3
      do i2 = ixMlo2, ixMhi2
      do i1 = ixMlo1, ixMhi1
         associate(rho => w(i1,i2,i3,rho_))
         v2 = sum(w(i1,i2,i3,mom(:))**2)/rho**2
         B2 = sum(w(i1,i2,i3,mag(:))**2)
         pth   = (mhd_gamma-1d0)*(w(i1,i2,i3,e_)-0.5d0*rho*v2-0.5d0*B2)
         J2 = sum(current(i1,i2,i3,:)**2)
         bME    = bME+0.5d0*B2*dvolume
         bIE    = bIE+pth/(mhd_gamma-1d0)*dvolume
         bheating = bheating+mhd_eta*J2*dvolume
         bmaxJ = max(bmaxJ, sqrt(J2))
         bmaxv = max(bmaxv, sqrt(v2))
         end associate
      end do
      end do
      end do
    end subroutine compute_block_stats

  end subroutine log_usr

end module mod_usr
