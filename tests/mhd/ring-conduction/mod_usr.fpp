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
  double precision, parameter :: ring_Busr   = 1.0d-5  ! B = ring_Busr/r (line-current field, eq. 63)

  double precision :: kperp_integral_time = 0.0d0

contains

  subroutine usr_init()
    use mod_global_parameters
    implicit none

    nwauxio = 2

    call set_coordinate_system("Cartesian_3D")
    usr_init_one_grid  => initonegrid_usr
    usr_aux_output     => extra_var_output
    usr_add_aux_names  => extra_var_names_output
    usr_write_analysis => write_error_csv
    usr_process_global => compute_kperp
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
    double precision :: B(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3)

    ! --- Step 1: cylindrical coords ---
    r(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = &
       dsqrt(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)**2 &
           + x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)**2)

    where (x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1) > 0.0d0)
      theta(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = &
         atan(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2) &
              / x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1))
    elsewhere (x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1) < 0.0d0)
      theta(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = &
         dpi - atan(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2) &
                    / abs(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)))
    elsewhere
      theta(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = 0.0d0
    end where

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

    ! --- Step 3: azimuthal B field from straight-line current: |B| = ring_Busr/r ---
    B(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = &
       ring_Busr / max(r(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3), 1.0d-15)
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(1)) = &
       B(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) &
       * dcos(theta(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) + 0.5d0*dpi)
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(2)) = &
       B(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) &
       * dsin(theta(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) + 0.5d0*dpi)
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

    double precision :: r(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: T(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    ! T = p/rho = (gamma-1)*(e_tot - 0.5*|m|^2/rho - 0.5*B^2) / rho  (conserved input)
    T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
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

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1) = &
       T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    r(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       dsqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2 &
           + x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2)

    where (r(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) > ring_rmin &
     .and. r(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) < ring_rmax)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+2) = &
         dabs(T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - 61.0d0/6.0d0)
    elsewhere
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+2) = &
         dabs(T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - T_cold)
    end where
  end subroutine extra_var_output

  subroutine extra_var_names_output(varnames)
    use mod_global_parameters
    implicit none
    character(len=*) :: varnames
    varnames = 'T L'
  end subroutine extra_var_names_output

  subroutine write_error_csv()
    use mod_global_parameters
    implicit none
    call spatial_integral
  end subroutine write_error_csv

  subroutine spatial_integral
    use mod_global_parameters
    implicit none

    double precision :: dvolume(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision :: integral_ipe(5), integral_w(5)
    integer          :: iigrid, igrid, ni, i
    character(len=100)  :: filename
    character(len=1024) :: line, datastr
    logical :: alive

    ni = 5
    integral_ipe    = 0.0d0
    integral_w      = 0.0d0
    integral_ipe(5) = bigdouble

    do iigrid = 1, igridstail
      igrid  = igrids(iigrid)
      block => ps(igrid)
      dxlevel(1) = rnode(rpdx1_,igrid)
      dxlevel(2) = rnode(rpdx2_,igrid)
      dxlevel(3) = rnode(rpdx3_,igrid)
      if (slab) then
        dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) = &
           rnode(rpdx1_,igrid) * rnode(rpdx2_,igrid) * rnode(rpdx3_,igrid)
      else
        dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) = &
           block%dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)
      end if

      integral_ipe(1) = integral_ipe(1) + &
         integral_grid(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3, &
                       ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3, &
                       ps(igrid)%w,ps(igrid)%x,dvolume,1)
      integral_ipe(2) = integral_ipe(2) + &
         integral_grid(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3, &
                       ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3, &
                       ps(igrid)%w,ps(igrid)%x,dvolume,2)
      integral_ipe(3) = max(integral_ipe(3), &
         integral_grid(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3, &
                       ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3, &
                       ps(igrid)%w,ps(igrid)%x,dvolume,3))
      integral_ipe(4) = max(integral_ipe(4), &
         integral_grid(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3, &
                       ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3, &
                       ps(igrid)%w,ps(igrid)%x,dvolume,4))
      integral_ipe(5) = min(integral_ipe(5), &
         integral_grid(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3, &
                       ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3, &
                       ps(igrid)%w,ps(igrid)%x,dvolume,5))
    end do

    call MPI_ALLREDUCE(integral_ipe(1:2),integral_w(1:2),2, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(integral_ipe(3:4),integral_w(3:4),2, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)
    call MPI_ALLREDUCE(integral_ipe(5),integral_w(5),1, &
                       MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)

    ! normalise by domain volume
    integral_w(1:2) = integral_w(1:2) / &
       ((xprobmax1-xprobmin1)*(xprobmax2-xprobmin2)*(xprobmax3-xprobmin3))
    integral_w(2) = dsqrt(integral_w(2))

    if (mype == 0) then
      write(filename,'(a,a)') trim(base_filename), 'Lerrors.csv'
      inquire(file=filename, exist=alive)
      if (alive) then
        open(unit=21,file=filename,form='formatted',status='old',position='append')
      else
        open(unit=21,file=filename,form='formatted',status='new')
        write(21,'(a)') 'time, L1 error, L2 error, L infinity error, Tmax, Tmin'
      end if
      write(datastr,'(es11.4,2a)') global_time, ', '
      line = datastr
      do i = 1, ni-1
        write(datastr,'(es12.5,2a)') integral_w(i), ', '
        line = trim(line)//trim(datastr)
      end do
      write(datastr,'(es12.5)') integral_w(ni)
      line = trim(line)//trim(datastr)
      write(21,'(a)') trim(line)
      close(21)
    end if

  end subroutine spatial_integral

  double precision function integral_grid(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,dvolume,intval)
    use mod_global_parameters
    implicit none
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
    integer, intent(in) :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    integer, intent(in) :: intval
    double precision, intent(in) :: &
       x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in) :: &
       dvolume(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(in) :: &
       w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    double precision :: r(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: T(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: L(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: ix1, ix2, ix3

    T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
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

    r(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       dsqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2 &
           + x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2)

    where (r(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) > ring_rmin &
     .and. r(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) < ring_rmax)
      L(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
         dabs(T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - 61.0d0/6.0d0)
    elsewhere
      L(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
         dabs(T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - T_cold)
    end where

    integral_grid = 0.0d0
    select case(intval)
    case(1)
      do ix3 = ixOmin3, ixOmax3
        do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
            integral_grid = integral_grid + L(ix1,ix2,ix3) * dvolume(ix1,ix2,ix3)
          end do
        end do
      end do
    case(2)
      do ix3 = ixOmin3, ixOmax3
        do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
            integral_grid = integral_grid + L(ix1,ix2,ix3)**2 * dvolume(ix1,ix2,ix3)
          end do
        end do
      end do
    case(3)
      integral_grid = maxval(L(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    case(4)
      integral_grid = maxval(T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    case(5)
      integral_grid = minval(T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    case default
      call mpistop("integral_grid: intval not defined")
    end select

  end function integral_grid

  subroutine compute_kperp(iit, qt)
    use mod_global_parameters
    implicit none
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt

    double precision :: T(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision :: lapT(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)
    double precision :: integrals_ipe(2), integrals_g(2)
    double precision :: dv, r_here, kperp
    integer          :: iigrid, igrid, ix1, ix2, ix3
    character(len=100) :: filename
    logical :: alive

    integrals_ipe = 0.0d0

    do iigrid = 1, igridstail
      igrid  = igrids(iigrid)
      block => ps(igrid)
      dxlevel(1) = rnode(rpdx1_,igrid)
      dxlevel(2) = rnode(rpdx2_,igrid)
      dxlevel(3) = rnode(rpdx3_,igrid)
      dv = dxlevel(1) * dxlevel(2) * dxlevel(3)

      ! T from conserved variables over full ghost range (stencil needs ±1 neighbour)
      T(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) = &
         (mhd_gamma - 1.0d0) * ( &
            ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,p_) &
          - 0.5d0*(ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,mom(1))**2 &
                 + ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,mom(2))**2 &
                 + ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,mom(3))**2) &
            / ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,rho_) &
          - 0.5d0*(ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,mag(1))**2 &
                 + ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,mag(2))**2 &
                 + ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,mag(3))**2) &
         ) / ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,rho_)

      ! 2D centred Laplacian (x, y only — no z-variation in this test)
      lapT(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) = &
         ( T(ixMlo1+1:ixMhi1+1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) &
         - 2.0d0*T(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) &
         + T(ixMlo1-1:ixMhi1-1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) ) / dxlevel(1)**2 &
       + ( T(ixMlo1:ixMhi1,ixMlo2+1:ixMhi2+1,ixMlo3:ixMhi3) &
         - 2.0d0*T(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3) &
         + T(ixMlo1:ixMhi1,ixMlo2-1:ixMhi2-1,ixMlo3:ixMhi3) ) / dxlevel(2)**2

      do ix3 = ixMlo3, ixMhi3
        do ix2 = ixMlo2, ixMhi2
          do ix1 = ixMlo1, ixMhi1
            r_here = dsqrt(ps(igrid)%x(ix1,ix2,ix3,1)**2 &
                         + ps(igrid)%x(ix1,ix2,ix3,2)**2)
            if (r_here > ring_rmin .and. r_here < ring_rmax) then
              integrals_ipe(1) = integrals_ipe(1) + lapT(ix1,ix2,ix3) * dv
              integrals_ipe(2) = integrals_ipe(2) &
                               + (T(ix1,ix2,ix3) - 61.0d0/6.0d0) * dv
            end if
          end do
        end do
      end do
    end do

    call MPI_ALLREDUCE(integrals_ipe, integrals_g, 2, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, icomm, ierrmpi)

    kperp_integral_time = kperp_integral_time + dt * integrals_g(1)

    if (mype == 0 .and. &
        global_time >= tsavelast(fileanalysis_) + dtsave(fileanalysis_) - smalldouble) then
      if (dabs(kperp_integral_time) > 0.0d0) then
        kperp = integrals_g(2) / kperp_integral_time
      else
        kperp = 0.0d0
      end if
      write(filename,'(a,a)') trim(base_filename), 'kperp.csv'
      inquire(file=filename, exist=alive)
      if (alive) then
        open(unit=22,file=filename,form='formatted',status='old',position='append')
      else
        open(unit=22,file=filename,form='formatted',status='new')
        write(22,'(a)') 'time, k_perp, k_perp/k_par'
      end if
      ! k_par = tc_kappa_par = 0.01 from agile.par
      write(22,'(es11.4,a,es12.5,a,es12.5)') &
         global_time, ', ', kperp, ', ', kperp / 0.01d0
      close(22)
    end if

  end subroutine compute_kperp

end module mod_usr
