! This module is for imposing randomized heating in 3D setup
! Users can define the parameters "trelax,tramp,ntimes,vlim^D,periods,variation,nwaves" in the mod_usr.t file.

!periods=300,variation=75  The time durations in seconds of each heating episode is typically 300s +/- 75s
!ntimes = 100              100 heating episodes, e.g., 500 min. Enough for the simulation
!nwaves = 1000             1000 waves for each episode
!si=5./6.                  spectral index
!vlim^D could be chosen as twice of the maximal spatial resolution in the simulation

module mod_random_heating


  use mod_global_parameters
  use mod_physics_vars

  implicit none

  public :: rh_source, generateTV

  double precision, allocatable :: va1(:), va2(:), va3(:)
  double precision, allocatable :: dtimearr(:),timearray(:)
  !$acc declare create(va1,va2,va3,dtimearr,timearray)
  !$omp declare target(va1,va2,va3,dtimearr,timearray)

  double precision :: periods=300.d0
  double precision :: variation=75.d0
  integer :: ntimes=100
  integer :: nwaves=1000
  double precision :: si=5.d0/6.d0
  integer :: vlim1=100
  integer :: vlim2=100
  integer :: vlim3=100

  double precision :: trelax=10.d0
  double precision :: tramp=5.d0
  !$acc declare create(vlim1,vlim2,vlim3)
  !$acc declare create(trelax,tramp,ntimes)
  !$omp declare target(vlim1,vlim2,vlim3)
  !$omp declare target(trelax,tramp,ntimes)


  contains

    !> Read this module's parameters from a file
  subroutine rh_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /rh_list/ periods, variation, ntimes, nwaves, si, &
                        vlim1, vlim2, vlim3, trelax, tramp

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, rh_list, end=111)
111    close(unitpar)
    end do

  end subroutine rh_read_params

  subroutine rh_init()
    use mod_global_parameters

    call rh_read_params(par_files)
    !$acc update device(trelax,tramp,ntimes,vlim1,vlim2,vlim3)
    !$omp target update to(trelax,tramp,ntimes,vlim1,vlim2,vlim3)

  end subroutine rh_init

  subroutine rh_source(qt,lQgrid,x)
    use mod_global_parameters
    !$acc routine seq
    !$omp declare target

    double precision, intent(in)    :: qt
    double precision, intent(inout) :: lQgrid
    double precision, intent(in)    :: x(1:ndim)

    double precision                :: tshift,tr,stt1,stt2,stt3
    double precision                :: tl1,tl2,tl3,tl4,tl5,tl6
    double precision                :: pulse1,pulse2,pulse3,pulse4,pulse5,pulse6
    integer                         :: iip1,iip2,iip3,i

    tshift = qt - trelax

    if(qt .lt. trelax) then
      tr=zero
      lQgrid=zero
      return
    elseif(qt .lt. trelax+tramp) then
      tr=(qt-trelax)/tramp
    else
      tr=one
    endif

    do i = 4, ntimes
      if(qt .lt. trelax + timearray(i)) then
          tl1 = dexp(-(tshift - timearray(i-3))**2/(0.5d0*dtimearr(i-3))**2)
          tl2 = dexp(-(tshift - timearray(i-2))**2/(0.5d0*dtimearr(i-2))**2)
          tl3 = dexp(-(tshift - timearray(i-1))**2/(0.5d0*dtimearr(i-1))**2)
          tl4 = dexp(-(tshift - timearray(i)  )**2/(0.5d0*dtimearr(i)  )**2)
          tl5 = dexp(-(tshift - timearray(i+1))**2/(0.5d0*dtimearr(i+1))**2)
          tl6 = dexp(-(tshift - timearray(i+2))**2/(0.5d0*dtimearr(i+2))**2)
          ! remapped on the unit cube
          stt1=(xprobmax1-xprobmin1)/(dble(vlim1)-1.d0)
          stt2=(xprobmax2-xprobmin2)/(dble(vlim2)-1.d0)
          stt3=(xprobmax3-xprobmin3)/(dble(vlim3)-1.d0)
          iip1=floor((x(1)-xprobmin1)/stt1)+1
          iip2=floor((x(2)-xprobmin2)/stt2)+1
          iip3=floor((x(3)-xprobmin3)/stt3)+1
          iip1=max(1,iip1)
          iip2=max(1,iip2)
          iip3=max(1,iip3)
          pulse1=va1(iip1+(i-4)*vlim1)*va2(iip2+(i-4)*vlim2)*va3(iip3+(i-4)*vlim3)
          pulse2=va1(iip1+(i-3)*vlim1)*va2(iip2+(i-3)*vlim2)*va3(iip3+(i-3)*vlim3)
          pulse3=va1(iip1+(i-2)*vlim1)*va2(iip2+(i-2)*vlim2)*va3(iip3+(i-2)*vlim3)
          pulse4=va1(iip1+(i-1)*vlim1)*va2(iip2+(i-1)*vlim2)*va3(iip3+(i-1)*vlim3)
          pulse5=va1(iip1+(i)*vlim1)*va2(iip2+(i)*vlim2)*va3(iip3+(i)*vlim3)
          pulse6=va1(iip1+(i+1)*vlim1)*va2(iip2+(i+1)*vlim2)*va3(iip3+(i+1)*vlim3)
          lQgrid = tl1 * pulse1 +  tl2 *pulse2  &
                        +tl3 * pulse3 +  tl4 *pulse4  &
                        +tl5 * pulse5 +  tl6 *pulse6
          exit
      endif
    enddo

    lQgrid=tr*lQgrid

  end subroutine rh_source

  subroutine generateTV()

    integer                      :: i,vx1,vx2,vx3,j
    ! generate random T, typically 300s +/- 75s
    call rh_init()

    allocate(dtimearr(ntimes))
    allocate(timearray(ntimes))

    if (mype==0) then
      call randomT(timearray)
      dtimearr = timearray
      do i = 2, ntimes
        dtimearr(i) = timearray(i) - timearray(i-1)
      enddo
    endif
    call MPI_BARRIER(icomm,ierrmpi)
    if (npe>1) then
      call MPI_BCAST(timearray,ntimes,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(dtimearr,ntimes,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif
    !$acc update device(timearray,dtimearr)
    !$omp target update to(timearray,dtimearr)

    ! spatial distribution

    vx1=ntimes*vlim1
    vx2=ntimes*vlim2
    vx3=ntimes*vlim3
    allocate(va1(vx1))
    allocate(va2(vx2))
    allocate(va3(vx3))
    if (mype==0) then
       call randomV(1,va1)
       call randomV(2,va2)
       call randomV(3,va3)
    endif
    call MPI_BARRIER(icomm,ierrmpi)
    if (npe>1) then
      call MPI_BCAST(va1,vx1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(va2,vx2,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(va3,vx3,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif
    !$acc update device(va1,va2,va3)
    !$omp target update to(va1,va2,va3)

  contains

    subroutine randomT(tb)
        double precision, intent(inout) :: tb(:)

        character(len=100)  :: filename
        double precision                :: tt1,tt,mm
        integer                         :: i

        tt = 0.d0
        do i = 1, ntimes
        call random_number(mm)
        tt1 = periods + variation * 2. * (mm - 0.5)
        tt = tt + tt1
        tb(i) = tt / unit_time
        end do

        write(filename,"(a)") "randomtimes.dat"
        open(unit=21,file=filename,form='formatted',status='replace')
        write(21,'(es12.4)') tb
        close(21)

    end subroutine randomT

    subroutine randomV(ival,va)
        integer, intent(in)             :: ival
        double precision, intent(inout) :: va(:)

        character(len=100)  :: filename
        character(len=10)   :: xdirection

        double precision, allocatable   :: vn(:,:),vm(:)
        double precision                :: lambda,ampl,lambda_min,lambda_max,rms,phase
        integer                         :: i, j, kk, vlim

        select case (ival)
        case (1)
            vlim = vlim1
        case (2)
            vlim = vlim2
        case (3)
            vlim = vlim3
        end select

        allocate(vn(nwaves, vlim))
        allocate(vm(vlim))
        lambda_min = 1.d0/dble(vlim)
        lambda_max = 1.d0
        do kk = 1, ntimes
            do i = 1, nwaves
            lambda = (dble(i-1)/dble(nwaves-1))*(lambda_max-lambda_min)+lambda_min
            ampl = lambda**si
            call random_number(phase)
            ! phase to vary between -2pi-->+2pi
            phase=4.0d0*dpi*(phase-0.5d0)
            do j = 1, vlim
                vn(i,j) = ampl*dsin(2.d0*dpi*(dble(j-1)/dble(vlim-1))/lambda+phase)
            end do
            end do
            rms = 0.d0
            vm = 0.d0
            do j = 1, vlim
            ! sum all the waves in grid point j
            do i = 1, nwaves
                vm(j) = vm(j) + vn(i, j)
            end do
            rms = rms + vm(j)**2
            end do
            vm = vm / dsqrt(rms/vlim)    ! normalization
            do j = 1, vlim
            va(vlim*(kk-1)+j) = vm(j)**2
            end do
        end do

        write(xdirection,'(I2.2)') ival
        write(filename,"(a)") "randompositions"//TRIM(xdirection)//".dat"
        open(unit=22,file=filename,form='formatted',status='replace')
        write(22,'(es12.4)') va
        close(22)

        deallocate(vm)
        deallocate(vn)

    end subroutine randomV
  end subroutine generateTV


end module mod_random_heating
