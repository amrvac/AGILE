module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision :: v0, rho0, p0, T0, pbeta0, b0, mach0
  integer :: netheat_

contains

  subroutine usr_init()

    unit_length        = 1.d9 ! in cm (1 Mm)
    unit_temperature   = 1.d6 ! in K (1Mk)
    unit_numberdensity = 1.d9 ! in cm^-3

    call usr_params_read(par_files)

    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output
    usr_init_vector_potential => initvecpot_usr
    usr_source => special_source
    usr_modify_output => set_output_vars
    usr_print_log => special_log

    call set_coordinate_system("Cartesian_3D")

    call phys_activate()

    ! Extra variable stored in dat files
    netheat_ = var_set_extravar("netheat", "netheat")

    if (mype == 0) then
      write(*,'(A)')        ' ========================================'
      write(*,'(A)')        '  Orszag-Tang MHD vortex (3D) with cooling'
      write(*,'(A,ES12.5)') '    rho0  = ', rho0
      write(*,'(A,ES12.5)') '    p0    = ', p0
      write(*,'(A,ES12.5)') '    T0    = ', T0
      write(*,'(A,ES12.5)') '    b0    = ', b0
      write(*,'(A,ES12.5)') '    v0    = ', v0
      write(*,'(A,ES12.5)') '    mach0 = ', mach0
      write(*,'(A,ES12.5)') '    pbeta0= ', pbeta0
      write(*,'(A)')        ' ========================================'
    end if

  end subroutine usr_init

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ pbeta0, T0, rho0, mach0

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

    if(mach0<smalldouble) call mpistop("Wrong input values: mach0 must be positive")
    if(rho0<smalldouble) call mpistop("Wrong input values: rho0 must be positive")
    if(T0<smalldouble) call mpistop("Wrong input values: T0 must be positive")
    p0 = rho0*T0
    if(pbeta0>smalldouble)then
        b0=dsqrt(2.0d0*p0/pbeta0)
    else
        b0=zero
    endif
    v0 = mach0*dsqrt(mhd_gamma*p0/rho0)

  end subroutine usr_params_read


  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_) = rho0

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(1)) = &
       -v0*sin(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(2)) = &
        v0*sin(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1))

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(3)) = zero

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_) = p0

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(1)) = &
       -b0 * sin(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(2)) = &
        b0 * sin(4.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1))
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(3)) = zero

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  ! Initialise the vectorpotential on the corners - used with staggered grid
  subroutine initvecpot_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, xC, A, idir)
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idir
    double precision, intent(in)       :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)      :: A(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if (idir == 3) then
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = &
         b0 * half / dpi * (half * cos(4.0_dp * dpi * &
         xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)) + &
         cos(2.0_dp * dpi * xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)))
    else
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = zero
    end if

  end subroutine initvecpot_usr


  subroutine special_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,&
     qtC,wCT,qt,w,x)
    #:if defined('COOLING')
    use mod_radiative_cooling, only: rc_fl, getvar_cooling_exact
    #:endif
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       iwmin,iwmax
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    #:if defined('COOLING')
    double precision :: winit(nw_phys), coolrate
    integer :: i1, i2, i3

    if(mhd_radiative_cooling)then
       ! Subtract background cooling: uniform state at rest with p0
       winit(iw_rho) = rho0
       winit(iw_mom(1)) = 0.0d0
       winit(iw_mom(2)) = 0.0d0
       winit(iw_mom(3)) = 0.0d0
       winit(iw_e) = p0/(mhd_gamma-1.0d0)
       winit(iw_mag(1)) = 0.0d0
       winit(iw_mag(2)) = 0.0d0
       winit(iw_mag(3)) = 0.0d0
       do i3=ixOmin3,ixOmax3
       do i2=ixOmin2,ixOmax2
       do i1=ixOmin1,ixOmax1
          call getvar_cooling_exact(qdt, winit, winit, x(i1,i2,i3,:), coolrate, rc_fl)
          w(i1,i2,i3,e_) = w(i1,i2,i3,e_) + qdt*coolrate
       end do
       end do
       end do
    endif
    #:endif

  end subroutine special_source


  !> Auxiliary variables for VTU conversion only (autoconvert=T).
  !> These are computed on-the-fly during conversion and then discarded. 
  !> For datfiles, we need to use var_set_extravar + usr_modify_output instead.
  subroutine specialvar_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
    use mod_functions_bfield
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw+nwauxio)
    double precision             :: normconv(0:nw+nwauxio)

    double precision :: wlocal(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: i1, i2, i3

    wlocal(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw) = &
       w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    ! Convert to primitive for temperature and invbeta
    call phys_to_primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wlocal,x)

    ! nw+1: Temperature = p/rho
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1) = &
       wlocal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_) / &
       wlocal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)

    ! nw+2: 1/beta = B^2 / (2*p)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+2) = &
       (wlocal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))**2 + &
        wlocal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(2))**2 + &
        wlocal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3))**2) / &
       (2.0d0*wlocal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_))

    ! nw+3: div B
    call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+3) = &
       divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine specialvar_output


  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames = 'Temp invbeta divb'
  end subroutine specialvarnames_output

  !> Populate extra variables before dat output.
  !> Note: w may be stale on host when this is called (data lives on GPU).
  !> We sync the physics variables from device first.
  subroutine set_output_vars(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    #:if defined('COOLING')
    use mod_radiative_cooling, only: rc_fl, getvar_cooling_exact
    #:endif
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: coolrate, dt_use
    double precision, save :: dt_prev = 0.0d0
    integer :: i1, i2, i3

    ! usr_modify_output is called before saveamrfile syncs device->host, so w on host is stale
    !$acc update host(block%w)

    ! At the final output dt can be ~0 (time_max reached exactly) which causes division by zero in getvar_cooling_exact
    if (dt > smalldouble) then
       dt_use = dt
       dt_prev = dt
    else
       dt_use = dt_prev
    end if

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,netheat_) = zero
    #:if defined('COOLING')
    if(mhd_radiative_cooling) then
       do i3=ixOmin3,ixOmax3
       do i2=ixOmin2,ixOmax2
       do i1=ixOmin1,ixOmax1
          call getvar_cooling_exact(dt_use, w(i1,i2,i3,1:nw_phys), &
             w(i1,i2,i3,1:nw_phys), x(i1,i2,i3,:), coolrate, rc_fl)
          w(i1,i2,i3,netheat_) = coolrate
       end do
       end do
       end do
    endif
    #:endif

    ! Push netheat back to device so saveamrfile's device->host sync picks it up
    !$acc update device(block%w(:,:,:,netheat_))

  end subroutine set_output_vars


  subroutine special_log()
    use mod_forest, only: nleafs_active
    use mod_global_parameters

    logical, save        :: opened = .false.
    integer              :: iigrid, igrid, amode, istatus(MPI_STATUS_SIZE)
    integer, parameter   :: my_unit = 20
    character(len=80)    :: filename
    character(len=1024)  :: line

    double precision :: local_min_T, local_max_T, local_min_rho, local_max_rho
    double precision :: local_min_p
    double precision :: global_mins(3), global_maxs(2), local_mins(3), local_maxs(2)

    ! Compute local min/max over all blocks on this processor
    local_min_T   =  huge(1.0d0)
    local_max_T   = -huge(1.0d0)
    local_min_rho =  huge(1.0d0)
    local_max_rho = -huge(1.0d0)
    local_min_p   =  huge(1.0d0)

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       block => ps(igrid)

       call compute_block_minmax(ps(igrid)%w, &
          local_min_T, local_max_T, local_min_rho, local_max_rho, &
          local_min_p)
    end do

    ! MPI reduce
    local_mins = (/ local_min_T, local_min_rho, local_min_p /)
    call MPI_ALLREDUCE(local_mins, global_mins, 3, MPI_DOUBLE_PRECISION, &
       MPI_MIN, icomm, ierrmpi)

    local_maxs = (/ local_max_T, local_max_rho /)
    call MPI_ALLREDUCE(local_maxs, global_maxs, 2, MPI_DOUBLE_PRECISION, &
       MPI_MAX, icomm, ierrmpi)

    if (mype == 0) then
       filename = trim(base_filename) // "_special.log"

       if (.not. opened) then
          if (restart_from_file == undefined) then
             open(unit=my_unit, file=trim(filename), status='replace')
             close(my_unit, status='delete')
          end if

          amode = ior(MPI_MODE_CREATE, MPI_MODE_WRONLY)
          amode = ior(amode, MPI_MODE_APPEND)
          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, MPI_INFO_NULL, &
             log_fh, ierrmpi)
          opened = .true.

          if (restart_from_file == undefined .or. reset_time) then
             line = 'it global_time dt min_T max_T min_rho max_rho min_p nleafs'
             call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'), &
                len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
          end if
       end if

       write(line, '(i8,2ES16.8,5ES16.8,i10)') &
          it, global_time, dt, &
          global_mins(1), global_maxs(1), &
          global_mins(2), global_maxs(2), &
          global_mins(3), nleafs_active

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'), &
          len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

  contains

    subroutine compute_block_minmax(w, bmin_T, bmax_T, bmin_rho, bmax_rho, bmin_p)
      double precision, intent(in)    :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,1:nw)
      double precision, intent(inout) :: bmin_T, bmax_T, bmin_rho, bmax_rho, bmin_p

      integer :: i1, i2, i3
      double precision :: rho, vmag2, Bmag2, pth, Temp

      do i3 = ixMlo3, ixMhi3
      do i2 = ixMlo2, ixMhi2
      do i1 = ixMlo1, ixMhi1
         rho   = w(i1,i2,i3,rho_)
         vmag2 = (w(i1,i2,i3,mom(1))**2 + w(i1,i2,i3,mom(2))**2 + &
                  w(i1,i2,i3,mom(3))**2) / rho**2
         Bmag2 = w(i1,i2,i3,mag(1))**2 + w(i1,i2,i3,mag(2))**2 + &
                 w(i1,i2,i3,mag(3))**2
         pth   = (mhd_gamma - 1.0d0) * (w(i1,i2,i3,e_) - &
                  0.5d0*rho*vmag2 - 0.5d0*Bmag2)
         Temp  = pth / rho

         bmin_T   = min(bmin_T,   Temp)
         bmax_T   = max(bmax_T,   Temp)
         bmin_rho = min(bmin_rho, rho)
         bmax_rho = max(bmax_rho, rho)
         bmin_p   = min(bmin_p,   pth)
      end do
      end do
      end do
    end subroutine compute_block_minmax

  end subroutine special_log


end module mod_usr
