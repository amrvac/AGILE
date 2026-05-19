module mod_usr

  use mod_amrvac
  use mod_physics
  implicit none

  double precision, parameter :: pi = 4d0*atan(1d0)

! User parameters.
! NOTE: These are in CGS units.
  integer :: seed
  double precision :: xi
  double precision :: chi
  double precision :: mach_H
  double precision :: rho_C
  double precision :: L_x
  double precision :: v_sh
  double precision :: v_off
  double precision, dimension(2:3) :: v_per
  double precision :: z0
  double precision :: d_z
  double precision :: sig_z
  double precision, dimension(2:3) :: f_til
  integer :: mode_root
  integer :: n_modes
  double precision, dimension(2) :: refine_z

! Precomputed derived quantities.
! NOTE: These are in CGS units.
  double precision, dimension(3) :: L
  double precision :: rho_H
  double precision :: P
  double precision :: v_C
  double precision :: v_H
  double precision :: t_sh
  double precision :: c_H
  double precision :: t_cool
  double precision :: T_H
  double precision :: T_C
  double precision :: T_mix
  double precision :: rho_mix

! Precomputed random numbers.
! (4, n_modes, n_modes, 2:3)
  double precision, dimension(:,:,:,:), allocatable :: rand_1

!$acc declare create(P, rho_H, rho_C, z0, d_z, sig_z, v_H, v_C, v_sh, mode_root, n_modes)
!$acc declare create(L, v_per, f_til, refine_z)
!$acc declare create(rand_1)

contains

  ! subroutine seed_rng(seed)
  !   implicit none
  !   integer, intent(in) :: seed
  !   integer :: n, i
  !   integer, allocatable :: seed_array(:)
  !
  !   call random_seed(size=n)
  !   allocate(seed_array(n))
  !   do i = 1, n
  !     seed_array(i) = seed+37*i
  !   end do
  !   call random_seed(put=seed_array)
  !   deallocate(seed_array)
  ! end subroutine seed_rng
  !
  ! function randn() result(z)
  !   implicit none
  !   double precision :: z
  !   double precision, save :: spare
  !   logical, save :: has_spare = .false.
  !   double precision :: u, v, s, fac
  !
  !   if (has_spare) then
  !     z = spare
  !     has_spare = .false.
  !     return
  !   end if
  !   do
  !     call random_number(u)
  !     call random_number(v)
  !     u = 2d0*u-1d0
  !     v = 2d0*v-1d0
  !     s = u*u + v*v
  !     if (s > 0d0 .and. s < 1d0) exit
  !   end do
  !   fac    = sqrt(-2d0*log(s)/s)
  !   z      = u*fac
  !   spare  = v*fac
  !   has_spare = .true.
  ! end function randn

  subroutine params_read_usr(files)
    implicit none
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/&
      seed, xi, chi, mach_H, rho_C, L_x, v_sh, v_off, v_per, z0, d_z, sig_z,&
      f_til, mode_root, n_modes, refine_z
    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status='old')
      read(unitpar, usr_list, end=111)
111   close(unitpar)
    end do
  end subroutine params_read_usr

  subroutine usr_init()
    use mod_global_parameters
    implicit none
    call set_coordinate_system('Cartesian_3D')
    call params_read_usr(par_files)
    unit_length      = 2.03d+16 ! L_x   CGS
    unit_density     = 1.00d-20 ! rho_C CGS
    unit_temperature = 2.00d+4  ! T_C   CGS
    usr_set_parameters => set_parameters_usr
    usr_init_one_grid  => init_one_grid_usr
    ! usr_special_bc     => specialbound_usr
    call phys_activate()
    if (mype == 0) call print_units()
  end subroutine usr_init

  subroutine print_params()
    implicit none
    character(len=*), parameter :: fmt = '(A20,1X,SP,ES12.2)'
    character(len=*), parameter :: fmti = '(A20,1X,I12)'

    print *, '------------------ PARAMETERS ------------------'
    print fmt,  'xi',           xi
    print fmt,  'chi',          chi
    print fmt,  'mach_H',       mach_H
    print fmt,  'rho_C',        rho_C
    print fmt,  'L_x',          L_x
    print fmt,  't_sh',         t_sh
    print fmt,  'v_off',        v_off
    print fmt,  'v_per(2)',     v_per(2)
    print fmt,  'v_per(3)',     v_per(3)
    print fmt,  'z0',           z0
    print fmt,  'd_z',          d_z
    print fmt,  'sig_z',        sig_z
    print fmt,  'f_til(2)',     f_til(2)
    print fmt,  'f_til(3)',     f_til(3)
    print fmti, 'mode_root',    mode_root
    print fmti, 'n_modes',      n_modes
    print fmt,  'L(1)',         L(1)
    print fmt,  'L(2)',         L(2)
    print fmt,  'L(3)',         L(3)
    print fmt,  'rho_H',        rho_H
    print fmt,  'P',            P
    print fmt,  'v_C',          v_C
    print fmt,  'v_H',          v_H
    print fmt,  'v_sh',         v_sh
    print fmt,  'c_H',          c_H
    print fmt,  't_cool',       t_cool
    print fmt,  'T_H',          T_H
    print fmt,  'T_C',          T_C
    print fmt,  'T_mix',        T_mix
    print fmt,  'rho_mix',      rho_mix
    print fmt,  'refine_z(1)',  refine_z(1)
    print fmt,  'refine_z(2)',  refine_z(2)
    print *, '------------------------------------------------'
! NOTE: See if these are roughly 0 and 1, verifies if random number generation
!   works as expected.
    print fmt, 'normal mean', sum(rand_1)/(8*n_modes**2)
    print fmt, 'normal var ', sum(rand_1**2)/(8*n_modes**2)-(sum(rand_1)/(8*n_modes**2))**2
  end subroutine

  subroutine print_units()
    use mod_global_parameters
    implicit none
    character(len=*), parameter :: fmt = '(A20,1X,SP,ES12.2)'

    print *, '-------------------- UNITS ---------------------'
    print fmt, 'unit_length',        unit_length
    print fmt, 'unit_time',          unit_time
    print fmt, 'unit_velocity',      unit_velocity
    print fmt, 'unit_density',       unit_density
    print fmt, 'unit_numberdensity', unit_numberdensity
    print fmt, 'unit_pressure',      unit_pressure
    print fmt, 'unit_temperature',   unit_temperature
    print *, '------------------------------------------------'
  end subroutine

  subroutine set_parameters_usr()
    use mod_global_parameters
    use mod_physics
    use mod_random
    implicit none
    integer :: i1, i2, i3, i4
    integer, parameter :: i8 = selected_int_kind(18) ! as in mod_random
    double precision :: a, b, mu
    type(rng_t) :: rng

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

    L = [xprobmax1-xprobmin1,&
         xprobmax2-xprobmin2,&
         xprobmax3-xprobmin3]
    L = L*L_x
    z0 = z0*L_x
    d_z = d_z*L_x
    sig_z = sig_z*L_x

    a = 1d0+4d0*He_abundance
    b = 2d0+3d0*He_abundance
    mu = a/b
    t_sh = L(1)/v_sh
    v_C = -v_sh/2+v_off
    v_H = +v_sh/2+v_off
    c_H = v_sh/mach_H
    t_cool = t_sh/xi
    rho_H = rho_C/chi
    rho_mix = sqrt(rho_C*rho_H)
    P = c_H**2*rho_H/hd_gamma
    T_H = (mu*mp_cgs)/(hd_gamma*kb_cgs)*c_H**2
    T_C = T_H/chi
    T_mix = sqrt(T_C*T_H)

!$acc update device(P, rho_H, rho_C, z0, d_z, sig_z, v_H, v_C, v_sh, mode_root, n_modes)
!$acc update device(L, v_per, f_til, refine_z)
!$acc update device(rand_1)

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
    integer :: ix, iy, nx, ny

    associate(&
      w_ => w(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, :),&
      x_ => x(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, :))

! STRATEGY:
!   We are given empty w_, and x_ in primitive form, and w_ should be returned
!   in conserved form.
!   1. Compute w_ in CGS units. Take care, x_ is in code units.
!   2. Convert w_ to code units.
!   3. Convert w_ to conserved.

      w_(:,:,:,iw_rho)    = rho_H+(rho_C-rho_H)/2*&
        (1-tanh((x_(:,:,:,3)*unit_length-z0)/d_z))
      w_(:,:,:,iw_e)      = P
      w_(:,:,:,iw_mom(1)) = v_H-v_sh/2*&
        (1-tanh((x_(:,:,:,3)*unit_length-z0)/d_z))

! We first compute the the noise terms in the y- and z-velocity buffers.
      w_(:,:,:,iw_mom(2)) = 0
      w_(:,:,:,iw_mom(3)) = 0
      do iy = 1, n_modes-1
      do ix = 1, n_modes-1
        nx = mode_root+ix
        ny = mode_root+iy
        w_(:,:,:,iw_mom(2)) = w_(:,:,:,iw_mom(2))+&
          (rand_1(1,ix,iy,2)*sin(2*pi*x_(:,:,:,1)*unit_length/L(1)*nx)+&
           rand_1(2,ix,iy,2)*cos(2*pi*x_(:,:,:,1)*unit_length/L(1)*nx))*&
          (rand_1(3,ix,iy,2)*sin(2*pi*x_(:,:,:,2)*unit_length/L(2)*ny)+&
           rand_1(4,ix,iy,2)*cos(2*pi*x_(:,:,:,2)*unit_length/L(2)*ny))
      end do
      end do
      do iy = 1, n_modes-1
      do ix = 1, n_modes-1
        nx = mode_root+ix
        ny = mode_root+iy
        w_(:,:,:,iw_mom(3)) = w_(:,:,:,iw_mom(3))+&
          (rand_1(1,ix,iy,3)*sin(2*pi*x_(:,:,:,1)*unit_length/L(1)*nx)+&
           rand_1(2,ix,iy,3)*cos(2*pi*x_(:,:,:,1)*unit_length/L(1)*nx))*&
          (rand_1(3,ix,iy,3)*sin(2*pi*x_(:,:,:,2)*unit_length/L(2)*ny)+&
           rand_1(4,ix,iy,3)*cos(2*pi*x_(:,:,:,2)*unit_length/L(2)*ny))
      end do
      end do
      w_(:,:,:,iw_mom(2)) = w_(:,:,:,iw_mom(2))*f_til(2)/(n_modes-1)
      w_(:,:,:,iw_mom(3)) = w_(:,:,:,iw_mom(3))*f_til(3)/(n_modes-1)

! Then the velocities themselves.
      w_(:,:,:,iw_mom(2)) = v_per(2)*(&
        sin(2*pi*x_(:,:,:,1)*unit_length/L(1)*mode_root)/2+&
        sin(2*pi*x_(:,:,:,2)*unit_length/L(2)*mode_root)/2+&
        w_(:,:,:,iw_mom(2)))*&
        exp(-((x_(:,:,:,3)*unit_length-z0)/sig_z)**2)
      w_(:,:,:,iw_mom(3)) = v_per(3)*(&
        sin(2*pi*x_(:,:,:,1)*unit_length/L(1)*mode_root)/2+&
        sin(2*pi*x_(:,:,:,2)*unit_length/L(2)*mode_root)/2+&
        w_(:,:,:,iw_mom(3)))*&
        exp(-((x_(:,:,:,3)*unit_length-z0)/sig_z)**2)

      w_(:,:,:,iw_rho)    = w_(:,:,:,iw_rho)/unit_density
      w_(:,:,:,iw_e)      = w_(:,:,:,iw_e)/unit_pressure
      w_(:,:,:,iw_mom(1)) = w_(:,:,:,iw_mom(1))/unit_velocity
      w_(:,:,:,iw_mom(2)) = w_(:,:,:,iw_mom(2))/unit_velocity
      w_(:,:,:,iw_mom(3)) = w_(:,:,:,iw_mom(3))/unit_velocity

    end associate

    call phys_to_conserved(&
      ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3,&
      ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3,&
      w, x)
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
      !dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:nw) :: w
    double precision, intent(in),&
      dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:ndim) :: x
    integer :: ix1, ix2, ix3
    double precision :: inv_gamma_m1

    inv_gamma_m1 = 1d0/(hd_gamma-1d0)

! NOTE: Zero-gradient for vy and vz, i.e. copy the closest interior cell into
!   the halo. There are no corners.
! NOTE: Due to some quirks with primitive and conservative conversion at the
!   moment, we apply the boundary condition in conservative form.
    select case(iB)
    case(5)
!$acc loop collapse(3) vector
      do ix3 = ixOmin3, ixOmax3
      do ix2 = ixOmin2, ixOmax2
      do ix1 = ixOmin1, ixOmax1
        w(ix1,ix2,ix3,iw_rho) = rho_C/unit_density
        w(ix1,ix2,ix3,iw_mom(1)) = v_C/unit_velocity*w(ix1,ix2,ix3,iw_rho)
        w(ix1,ix2,ix3,iw_mom(2)) = w(ix1,ix2,ixOmax3+1,iw_mom(2))
        w(ix1,ix2,ix3,iw_mom(3)) = w(ix1,ix2,ixOmax3+1,iw_mom(3))
        w(ix1,ix2,ix3,iw_e) = P/unit_pressure*inv_gamma_m1+&
          5d-1*sum(w(ix1,ix2,ix3,iw_mom(1:ndim))**2)/w(ix1,ix2,ix3,iw_rho)
      end do
      end do
      end do
    case(6)
!$acc loop collapse(3) vector
      do ix3 = ixOmin3, ixOmax3
      do ix2 = ixOmin2, ixOmax2
      do ix1 = ixOmin1, ixOmax1
        w(ix1,ix2,ix3,iw_rho) = rho_H/unit_density
        w(ix1,ix2,ix3,iw_mom(1)) = v_H/unit_velocity*w(ix1,ix2,ix3,iw_rho)
        w(ix1,ix2,ix3,iw_mom(2)) = w(ix1,ix2,ixOmin3-1,iw_mom(2))
        w(ix1,ix2,ix3,iw_mom(3)) = w(ix1,ix2,ixOmin3-1,iw_mom(3))
        w(ix1,ix2,ix3,iw_e) = P/unit_pressure*inv_gamma_m1+&
          5d-1*sum(w(ix1,ix2,ix3,iw_mom(1:ndim))**2)/w(ix1,ix2,ix3,iw_rho)
      end do
      end do
      end do
    end select

  end subroutine specialbound_usr

  subroutine usr_refine_grid(&
    igrid, level,&
    ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3,&
    ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3,&
    qt, w, x, refine, coarsen)
#ifdef _OPENACC
! NOTE: The Cray compiler fails when trying to inline this routine, for now
!   disable inlining for Cray.
    !dir$ inlinenever usr_refine_grid
#endif
    !$acc routine seq
    use mod_global_parameters
    implicit none
    integer, intent(in) :: igrid, level
    integer, intent(in) :: ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3
    integer, intent(in) :: ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3
    double precision, intent(in) :: qt
    double precision, intent(in),&
      dimension(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, 1:nw) :: w
    double precision, intent(in),&
      dimension(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, 1:ndim) :: x
    integer, intent(inout) :: refine, coarsen

    associate(&
      w_ => w(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, :),&
      x_ => x(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, :))

      if (qt == 0d0) then
        if (any(x_(:,:,:,3) > refine_z(1)) .and.&
            any(x_(:,:,:,3) < refine_z(2))) then
          coarsen = -1
          refine = 1
        end if
      end if

    end associate

  end subroutine usr_refine_grid

end module mod_usr
