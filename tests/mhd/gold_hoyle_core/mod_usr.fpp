#:mute
#:include "../../../src/mod_gpu_directives.fpp"
#:endmute
#:if GEOM != 'cylindrical' or not (defined('B0SPLIT') and defined('COOLING') and defined('HYPERTC_ANISO') and defined('SOURCE_USR') and defined('REFINE_USR'))
  #:stop "gold_hoyle_core: cylindrical mesh with mhd_b0_split, cooling, anisotropic conduction, mhd_source_usr, refine_usr"
#:endif

!> Gold-Hoyle loop with a Gaussian cool core on the cylindrical mesh x = (r, z, phi):
!> the AGILE counterpart of legolas_runs/gold_hoyle_core, thermally unstable.
!>
!>   g(r)   = exp(-(r - r_core)^2 / (2 w_core^2))
!>   p0(r)  = rho_bg t_bg + p_excess g       T0(r) = t_bg + (t_core - t_bg) g       rho0 = p0 / T0
!>   B0_phi = alpha r / (1 + alpha^2 r^2)    B0_z  = sqrt(1 / (1 + alpha^2 r^2)^2 - 2 p_excess g)
!>
!> p_excess = 0 is the force-free loop at uniform pressure; otherwise the dip of B0_z confines the core.
!> lundquist = .true. replaces the field by the linear force-free one, B0_phi = J1(alpha r), B0_z = J0(alpha r)
!> (curl B0 = alpha B0; use with p_excess = 0).
!> Code units as in Legolas: B_c = L = v_A = 1, p = n k_B T (He_abundance = H_ion_fr = 0).
!> Rosner cooling, Spitzer conduction along B, lambda_perp times Braginskii across it,
!> and a heating per unit mass that cancels the radiative losses and the perpendicular
!> conduction of the equilibrium.
!> footpoint_length > 0 adds the conductive loss kappa_par T / footpoint_length^2 per unit
!> volume, also cancelled in the equilibrium.
!> The field is B0; the evolved B1 starts at zero.
module mod_usr
  use mod_amrvac
  use mod_physics
  use mod_gimli

  implicit none

  double precision :: rho_cgs = 1.6727d-15
  double precision :: bc_gauss = 10.0d0
  double precision :: l_cgs = 1.0d10
  double precision :: alpha = 20.0d0
  double precision :: rho_bg = 0.6993006993d0
  double precision :: t_bg = 2.255573870044d-3
  double precision :: t_core = 1.577324384646d-3
  double precision :: r_core = 0.6d0
  double precision :: w_core = 0.4d0
  double precision :: p_excess = 0.0d0
  logical :: lundquist = .false.
  double precision :: lambda_perp = 1.0d0
  logical :: heating = .true.
  double precision :: footpoint_length = 0.0d0
  !> Legolas eigenmode: added to the initial state if add_mode, and the template of mode_amp
  character(len=std_len) :: ldat_file = ''
  logical :: add_mode = .true.
  !> isobaric density ring: rho0 (1 + pulse_amp exp(-(r-pulse_r)^2/2 pulse_wr^2 - (z-pulse_z)^2/2 pulse_wz^2))
  double precision :: pulse_amp = 0.0d0
  double precision :: pulse_r = 0.92d0
  double precision :: pulse_wr = 0.025d0
  double precision :: pulse_z = 3.141592653589793d0
  double precision :: pulse_wz = 0.6d0
  !> isobaric broadband noise of Legolas (gold_hoyle_core/smod_user_defined.f08, p3 cosines, realisation p4):
  !> rho0 (1 + noise_amp noise(r) cos(noise_m phi + noise_k z)), switched off within 0.1 of the axis and the wall
  double precision :: noise_amp = 0.0d0
  integer :: noise_cosines = 200
  integer :: noise_seed = 1
  integer :: noise_m = 0
  double precision :: noise_k = 1.0d0
  !> the shell refine_rmin < r < refine_rmax is refined to refine_max_level
  double precision :: refine_rmin = 0.75d0
  double precision :: refine_rmax = 2.0d0

  double precision :: p0 = 0.001d0
  double precision :: sink_coeff = 0.0d0
  !> div(kappa_perp grad T0) / rho0, tabulated at r = (i-1) heat_dr
  integer, parameter :: n_heat = 32768
  double precision :: heat_table(n_heat) = 0.0d0
  double precision :: heat_dr = 1.0d0
  ${GPU_DECLARE_COPYIN('alpha, p0, p_excess, lundquist, t_bg, t_core, r_core, w_core, heating, refine_rmin, refine_rmax, sink_coeff')}$
  ${GPU_DECLARE_CREATE('heat_table, heat_dr')}$

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_constants, only: mp_cgs, kB_cgs
    double precision, parameter :: coulomb_log = 22.0d0

    call params_read_user(par_files)
    if (ldat_file /= '') call read_legolas_data(ldat_file, 21)

    unit_length      = l_cgs
    unit_density     = rho_cgs
    unit_temperature = bc_gauss**2 * mp_cgs / (4.0d0 * dpi * rho_cgs * kB_cgs)
    p0 = rho_bg * t_bg
    ${GPU_UPDATE_DEVICE('alpha, p0, p_excess, lundquist, t_bg, t_core, r_core, w_core, heating, refine_rmin, refine_rmax, sink_coeff')}$

    nwauxio = 3
    usr_init_one_grid    => initonegrid_usr
    usr_special_bc       => specialbound_usr
    usr_before_main_loop => print_setup
    usr_print_log        => write_log
    usr_aux_output       => lorentz_output
    usr_add_aux_names    => lorentz_names

    call phys_activate()

    ! Legolas, cgs: kappa_par = 1.8e-5 T^2.5 / ln L, kappa_perp = lambda 8.2e-13 ln L^2 kappa_par n^2 / (B^2 T^3)
    tc_kappa0_par = 1.8d-5 / coulomb_log * unit_temperature**3.5d0 &
       / (unit_length * unit_density * unit_velocity**3)
    tc_kappa_perp_factor = lambda_perp * 8.2d-13 * coulomb_log**2 * (htc_Cchi &
       * unit_numberdensity / (unit_magneticfield * unit_temperature**1.5d0))**2
    ${GPU_UPDATE_DEVICE('tc_kappa0_par, tc_kappa_perp_factor')}$
    if (footpoint_length > 0.0d0) sink_coeff = tc_kappa0_par / footpoint_length**2
    ${GPU_UPDATE_DEVICE('sink_coeff')}$

  end subroutine usr_init

  subroutine params_read_user(files)
    use mod_global_parameters, only: unitpar
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ rho_cgs, bc_gauss, l_cgs, alpha, rho_bg, t_bg, t_core, r_core, w_core, p_excess, lundquist, &
       lambda_perp, heating, footpoint_length, ldat_file, add_mode, pulse_amp, pulse_r, pulse_wr, pulse_z, pulse_wz, &
       refine_rmin, refine_rmax, noise_amp, noise_cosines, noise_seed, noise_m, noise_k

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read_user

  pure double precision function core_profile(r) result(g)
#ifdef _OPENMP
    ${GPU_ROUTINE_VECTOR()}$
#else
    ${GPU_ROUTINE_SEQ()}$
#endif
    double precision, intent(in) :: r

    g = exp(-0.5d0 * ((abs(r) - r_core) / w_core)**2)

  end function core_profile

  pure double precision function equilibrium_pressure(r) result(p)
#ifdef _OPENMP
    ${GPU_ROUTINE_VECTOR()}$
#else
    ${GPU_ROUTINE_SEQ()}$
#endif
    double precision, intent(in) :: r

    p = p0 + p_excess * core_profile(r)

  end function equilibrium_pressure

  pure double precision function equilibrium_density(r) result(rho)
#ifdef _OPENMP
    ${GPU_ROUTINE_VECTOR()}$
#else
    ${GPU_ROUTINE_SEQ()}$
#endif
    double precision, intent(in) :: r

    rho = equilibrium_pressure(r) / (t_bg + (t_core - t_bg) * core_profile(r))

  end function equilibrium_density

  !> r kappa_perp dT0/dr of the equilibrium
  pure double precision function equilibrium_heat_flux(r) result(flux)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in) :: r
    double precision :: rho, temp, x(1:ndim), b(1:ndir)

    rho  = equilibrium_density(r)
    temp = equilibrium_pressure(r) / rho
    x    = (/ r, 0.0d0, 0.0d0 /)
    call usr_b0_field(x, b)
    flux = -r * tc_kappa_perp_factor * tc_kappa0_par * temp**2.5d0 &
       / (1.0d0 + (htc_Cchi * temp**1.5d0 / rho)**2 * sum(b(1:3)**2)) &
       * (t_core - t_bg) * (r - r_core) / w_core**2 * core_profile(r)

  end function equilibrium_heat_flux

  !> -div(r kappa_perp dT0/dr) / (r rho0): the perpendicular-conduction heating of the equilibrium
  pure double precision function equilibrium_heating(r) result(h)
    double precision, intent(in) :: r
    double precision, parameter :: flux_step = 1.0d-4

    h = (equilibrium_heat_flux(r + flux_step) - equilibrium_heat_flux(r - flux_step)) &
       / (2.0d0 * flux_step * max(r, flux_step) * equilibrium_density(r))

  end function equilibrium_heating

  pure double precision function noise(r) result(f)
    double precision, intent(in) :: r
    integer :: j

    f = 0.0d0
    do j = 1, noise_cosines
       f = f + cos(j * dpi * r + 2.0d0 * dpi * modulo(sin(12.9898d0 * j + 78.233d0 * noise_seed) * 43758.5453d0, 1.0d0))
    end do
    f = f / sqrt(0.5d0 * noise_cosines) * 0.25d0 * (1.0d0 + tanh((r - 0.1d0) / 0.02d0)) * (1.0d0 - tanh((r - 0.9d0) / 0.02d0))

  end function noise

  pure subroutine usr_equilibrium(x, rho, p)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: rho, p

    p   = equilibrium_pressure(x(1))
    rho = equilibrium_density(x(1))

  end subroutine usr_equilibrium

  !> (B_r, B_z, B_phi) of the rope
  pure subroutine usr_b0_field(x, b)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: b(3)
    double precision :: gh

    b(1) = 0.0d0
    if (lundquist) then
       b(2) = bessel_series(0, alpha * x(1))
       b(3) = bessel_series(1, alpha * x(1))
       return
    end if
    gh   = 1.0d0 / (1.0d0 + (alpha * x(1))**2)
    b(2) = sqrt(gh**2 - 2.0d0 * p_excess * core_profile(x(1)))
    b(3) = alpha * x(1) * gh

  end subroutine usr_b0_field

  !> J_order(x) by its power series; exact to round-off for |x| < 3
  pure double precision function bessel_series(order, x) result(f)
    ${GPU_ROUTINE_SEQ()}$
    integer, intent(in)          :: order
    double precision, intent(in) :: x
    double precision :: term
    integer :: k

    term = merge(1.0d0, 0.5d0 * x, order == 0)
    f = term
    do k = 1, 25
       term = -term * (0.5d0 * x)**2 / (k * (k + order))
       f = f + term
    end do

  end function bessel_series

  !> curl B0 in the mesh's component order (r, z, phi), which is left-handed:
  !> j_z = -(r B_phi)' / r, j_phi = B_z'
  pure subroutine usr_j0_field(x, j)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: j(3)
    double precision :: gh, b(3)

    call usr_b0_field(x, b)
    j(1) = 0.0d0
    if (lundquist) then
       j(2) = -alpha * b(2)
       j(3) = -alpha * b(3)
       return
    end if
    gh = 1.0d0 / (1.0d0 + (alpha * x(1))**2)
    j(2) = -2.0d0 * alpha * gh**2
    j(3) = (-2.0d0 * alpha**2 * x(1) * gh**3 &
       + p_excess * (x(1) - r_core) / w_core**2 * core_profile(x(1))) / b(2)

  end subroutine usr_j0_field

  subroutine usr_refine_grid(level,qt,w,x,refineflag,coarsenflag,norefineflag,nocoarsenflag)
    ${GPU_ROUTINE_SEQ()}$
    integer, intent(in)             :: level
    double precision, intent(in)    :: qt
    double precision, intent(in)    :: x(1:3)
    double precision, intent(in)    :: w(1:nw)
    logical, intent(inout)          :: refineflag, coarsenflag, norefineflag, nocoarsenflag

    if (x(1) > refine_rmin .and. x(1) < refine_rmax) then
       refineflag = .true.
       norefineflag = .false.
       coarsenflag = .false.
       nocoarsenflag = .true.
    end if

  end subroutine usr_refine_grid

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision                :: r, z
    integer                         :: ix1, ix2, ix3, iw

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw) = 0.0d0
    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             r = x(ix1,ix2,ix3,1)
             z = x(ix1,ix2,ix3,2)
             w(ix1,ix2,ix3,rho_) = equilibrium_density(r) * (1.0d0 + pulse_amp * exp( &
                -0.5d0 * ((r - pulse_r) / pulse_wr)**2 - 0.5d0 * ((z - pulse_z) / pulse_wz)**2))
             if (noise_amp > 0.0d0) w(ix1,ix2,ix3,rho_) = w(ix1,ix2,ix3,rho_) + equilibrium_density(r) * noise_amp &
                * noise(r) * cos(noise_m * x(ix1,ix2,ix3,3) + noise_k * z)
             w(ix1,ix2,ix3,p_)   = equilibrium_pressure(r)
          end do
       end do
    end do

    if (ldat_file /= '' .and. add_mode) then
       do iw = 1, nw
          call add_perturbation_to_w_array(ixGmin1,ixGmin2,ixGmin3,&
             ixGmax1,ixGmax2,ixGmax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,iw,x)
       end do
    end if

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  !> the wall: the equilibrium at rest, B1 = 0, no heat flux
  subroutine specialbound_usr(qt, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB, w, x)
    use mod_global_parameters
    ${GPU_ROUTINE_VECTOR()}$
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3
    integer, intent(in)             :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    integer, intent(in)             :: iB
    double precision, intent(in)    :: qt
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer                         :: ix1, ix2, ix3

#ifndef _OPENMP
    ${GPU_LOOP_VECTOR("collapse(3)")}$
#endif
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
#:if defined('EQUISPLIT')
             w(ix1,ix2,ix3,iw_rho)    = 0.0d0
             w(ix1,ix2,ix3,iw_e)      = 0.0d0
#:else
             w(ix1,ix2,ix3,iw_rho)    = equilibrium_density(x(ix1,ix2,ix3,1))
             w(ix1,ix2,ix3,iw_e)      = equilibrium_pressure(x(ix1,ix2,ix3,1)) * inv_gamma_1
#:endif
             w(ix1,ix2,ix3,iw_mom(1)) = 0.0d0
             w(ix1,ix2,ix3,iw_mom(2)) = 0.0d0
             w(ix1,ix2,ix3,iw_mom(3)) = 0.0d0
             w(ix1,ix2,ix3,iw_mag(1)) = 0.0d0
             w(ix1,ix2,ix3,iw_mag(2)) = 0.0d0
             w(ix1,ix2,ix3,iw_mag(3)) = 0.0d0
             w(ix1,ix2,ix3,psi_)      = 0.0d0
             w(ix1,ix2,ix3,q_(1))     = 0.0d0
             w(ix1,ix2,ix3,q_(2))     = 0.0d0
             w(ix1,ix2,ix3,q_(3))     = 0.0d0
          end do
       end do
    end do

  end subroutine specialbound_usr

  !> heating per unit mass: the perpendicular conduction of the equilibrium, and its radiative loss at x, from the
  !> cooling's own integrator so the two cancel to round-off at t = 0
  subroutine addsource_usr(qdt, qt, wCT, wCTprim, wnew, x, split)
    use mod_radiative_cooling, only: rc_fl, getvar_cooling_exact
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)    :: qdt, qt, wCT(nw_phys), wCTprim(nw_phys)
    double precision, intent(inout) :: wnew(nw_phys)
    double precision, intent(in)    :: x(1:ndim)
    logical, intent(in)             :: split
    double precision :: coolrate, w_eq(nw_phys), temp, temp_eq, s
    integer :: i

    if (.not. heating) return
    w_eq(:)      = 0.0d0
    w_eq(iw_rho) = equilibrium_density(x(1))
    w_eq(iw_e)   = equilibrium_pressure(x(1)) * inv_gamma_1
    call getvar_cooling_exact(qdt, w_eq, w_eq, x, coolrate, rc_fl)
    wnew(iw_e) = wnew(iw_e) + qdt * coolrate * wCT(iw_rho) / w_eq(iw_rho)
    if (t_core /= t_bg) then
       s = abs(x(1)) / heat_dr
       i = min(int(s) + 1, n_heat - 1)
       s = s - (i - 1)
       wnew(iw_e) = wnew(iw_e) - qdt * ((1.0d0 - s) * heat_table(i) + s * heat_table(i + 1)) * wCT(iw_rho)
    end if
    temp    = max(wCTprim(iw_e) / wCTprim(iw_rho), 0.0d0)
    temp_eq = equilibrium_pressure(x(1)) / w_eq(iw_rho)
    wnew(iw_e) = wnew(iw_e) + qdt * sink_coeff * (temp_eq**3.5d0 * wCT(iw_rho) / w_eq(iw_rho) - temp**3.5d0)

  end subroutine addsource_usr

  subroutine print_setup()
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    use mod_radiative_cooling, only: rc_fl, findL
    double precision :: r(3), rho, temp, bmag, lam, kpar, kperp, b(3), x_loc(1:ndim)
    integer          :: i

    if (He_abundance /= 0.0d0 .or. H_ion_fr /= 0.0d0) call mpistop(&
       "gold_hoyle_core: the unit choice assumes He_abundance = 0 and H_ion_fr = 0")
    heat_dr = 1.25d0 * xprobmax1 / (n_heat - 1)
    if (t_core /= t_bg) then
       do i = 1, n_heat
          heat_table(i) = equilibrium_heating((i - 1) * heat_dr)
       end do
    end if
    ${GPU_UPDATE_DEVICE('heat_table, heat_dr')}$

    if (mype /= 0) return

    write(*,'(A)')          ' Gold-Hoyle loop with a cool core'
    write(*,'(A,ES12.5,A)') '   L                = ', unit_length, ' cm'
    write(*,'(A,ES12.5,A)') '   unit density     = ', unit_density, ' g cm^-3'
    write(*,'(A,ES12.5,A)') '   B_c              = ', unit_magneticfield, ' G'
    write(*,'(A,ES12.5,A)') '   unit temperature = ', unit_temperature, ' K'
    write(*,'(A,ES12.5,A)') '   tau_A            = ', unit_time, ' s'
    write(*,'(A,3ES12.5)')  '   alpha, p0, p_ex  = ', alpha, p0, p_excess
    write(*,'(A,L2)')       '   lundquist field  = ', lundquist
    write(*,'(A,ES12.5)')   '   lambda_perp      = ', lambda_perp
    write(*,'(A,L2)')       '   heating balance  = ', heating
    write(*,'(A,ES12.5)')   '   footpoint length = ', footpoint_length
    write(*,'(A)')          '      r       T [K]       rho0        kappa_par' // &
       '   kappa_perp  Lambda      tau_cool    beta'
    r = (/ 0.0d0, r_core, 1.0d0 /)
    do i = 1, 3
       rho   = equilibrium_density(r(i))
       temp  = equilibrium_pressure(r(i)) / rho
       x_loc = (/ r(i), 0.0d0, 0.0d0 /)
       call usr_b0_field(x_loc, b)
       bmag  = sqrt(sum(b(1:3)**2))
       call findL(temp, lam, rc_fl)
       kpar  = tc_kappa0_par * temp**2.5d0
       kperp = tc_kappa_perp_factor * kpar / (1.0d0 + (htc_Cchi * bmag * temp**1.5d0 / rho)**2)
       write(*,'(3X,F6.3,7ES12.4)') r(i), temp * unit_temperature, rho, kpar, kperp, &
          lam, rho * temp * inv_gamma_1 / max(rho**2 * lam, smalldouble), 2.0d0 * rho * temp / bmag**2
    end do

  end subroutine print_setup

  !> |j x B|, |j| |B| and div B of the total field on the interior of one block
  subroutine lorentz(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,jxb,jb,divb)
    use mod_global_parameters
    use mod_functions_bfield, only: get_current, get_divb
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out) :: jxb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out) :: jb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out) :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: j1(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)
    double precision :: j(3), b(3), f(3), b0(3), j0(3)
    integer          :: idirmin, ix1, ix2, ix3

    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,j1)
    call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)

    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
             call usr_b0_field(x(ix1,ix2,ix3,1:ndim), b0)
             call usr_j0_field(x(ix1,ix2,ix3,1:ndim), j0)
             j(1:3) = j1(ix1,ix2,ix3,1:3) + j0(1:3)
             b(1:3) = w(ix1,ix2,ix3,iw_mag(1:3)) + b0(1:3)
             f(1) = j(2)*b(3) - j(3)*b(2)
             f(2) = j(3)*b(1) - j(1)*b(3)
             f(3) = j(1)*b(2) - j(2)*b(1)
             jxb(ix1,ix2,ix3) = sqrt(sum(f(1:3)**2))
             jb(ix1,ix2,ix3)  = sqrt(sum(j(1:3)**2)) * sqrt(sum(b(1:3)**2))
          end do
       end do
    end do

  end subroutine lorentz

  !> it t dt Ekin Emag vmax dB_L2 drho_max sigmaJ jxb_max divb_mean divb_max Tmin Tmax drho_L2 mode_amp
  !>
  !> drho_max, drho_L2: maximum and rms of rho/rho0 - 1; sigmaJ: current-weighted sine of
  !> the angle between j and B; mode_amp: projection of rho - rho0 on the Legolas mode's density.
  subroutine write_log()
    use mod_global_parameters
    integer, parameter :: nsum = 11, nmax = 6
    double precision   :: jxb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision   :: jb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision   :: divb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision   :: sums_ipe(nsum), sums(nsum), maxs_ipe(nmax), maxs(nmax)
    double precision   :: dv, rho, rho_eq, v2, temp, b1(3), b0(3), x_loc(1:ndim), mode(8)
    integer            :: iigrid, igrid, ix1, ix2, ix3
    integer, parameter :: my_unit = 21
    logical, save      :: opened = .false.

    sums_ipe = 0.0d0
    maxs_ipe = 0.0d0
    maxs_ipe(5:6) = -bigdouble

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       block => ps(igrid)
       dxlevel(1) = rnode(rpdx1_,igrid)
       dxlevel(2) = rnode(rpdx2_,igrid)
       dxlevel(3) = rnode(rpdx3_,igrid)

       call lorentz(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid)%w,ps(igrid)%x,jxb,jb,divb)

       do ix3 = ixMlo3, ixMhi3
          do ix2 = ixMlo2, ixMhi2
             do ix1 = ixMlo1, ixMhi1
                dv = block%dvolume(ix1,ix2,ix3)
                x_loc(1:ndim) = ps(igrid)%x(ix1,ix2,ix3,1:ndim)
                call usr_b0_field(x_loc, b0)
                rho_eq  = equilibrium_density(x_loc(1))
                rho     = ps(igrid)%w(ix1,ix2,ix3,iw_rho)
                v2      = sum(ps(igrid)%w(ix1,ix2,ix3,iw_mom(1:3))**2) / rho**2
                b1(1:3) = ps(igrid)%w(ix1,ix2,ix3,iw_mag(1:3))
                temp    = gamma_1 * (ps(igrid)%w(ix1,ix2,ix3,iw_e) - 0.5d0 * rho * v2 &
                          - 0.5d0 * sum(b1(1:3)**2)) / rho

                sums_ipe(1) = sums_ipe(1) + 0.5d0 * rho * v2 * dv
                sums_ipe(2) = sums_ipe(2) + 0.5d0 * sum((b1(1:3) + b0(1:3))**2) * dv
                sums_ipe(3) = sums_ipe(3) + sum(b1(1:3)**2) * dv
                sums_ipe(4) = sums_ipe(4) + sum(b0(1:3)**2) * dv
                sums_ipe(5) = sums_ipe(5) + jxb(ix1,ix2,ix3) * dv
                sums_ipe(6) = sums_ipe(6) + jb(ix1,ix2,ix3) * dv
                sums_ipe(7) = sums_ipe(7) + abs(divb(ix1,ix2,ix3)) * dv
                sums_ipe(8) = sums_ipe(8) + dv
                sums_ipe(9) = sums_ipe(9) + (rho / rho_eq - 1.0d0)**2 * dv
                if (ldat_file /= '') then
                   mode = legolas_perturbation(x_loc)
                   sums_ipe(10) = sums_ipe(10) + (rho - rho_eq) * mode(1) * dv
                   sums_ipe(11) = sums_ipe(11) + mode(1)**2 * dv
                end if

                maxs_ipe(1) = max(maxs_ipe(1), sqrt(v2))
                maxs_ipe(2) = max(maxs_ipe(2), abs(rho / rho_eq - 1.0d0))
                maxs_ipe(3) = max(maxs_ipe(3), jxb(ix1,ix2,ix3))
                maxs_ipe(4) = max(maxs_ipe(4), abs(divb(ix1,ix2,ix3)))
                maxs_ipe(5) = max(maxs_ipe(5), -temp)
                maxs_ipe(6) = max(maxs_ipe(6), temp)
             end do
          end do
       end do
    end do

    call MPI_ALLREDUCE(sums_ipe, sums, nsum, MPI_DOUBLE_PRECISION, MPI_SUM,&
       icomm, ierrmpi)
    call MPI_ALLREDUCE(maxs_ipe, maxs, nmax, MPI_DOUBLE_PRECISION, MPI_MAX,&
       icomm, ierrmpi)
    if (mype /= 0) return

    if (.not. opened) then
       if (restart_from_file == undefined) then
          open(unit=my_unit, file=trim(base_filename) // ".log", status='replace')
          write(my_unit,'(A)') '# it t dt Ekin Emag vmax dB_L2 drho_max ' // &
             'sigmaJ jxb_max divb_mean divb_max Tmin Tmax drho_L2 mode_amp'
       else
          open(unit=my_unit, file=trim(base_filename) // ".log", status='old', position='append')
       end if
       opened = .true.
    end if
    write(my_unit,'(I8,15ES16.8)') it, global_time, dt, sums(1), sums(2), maxs(1),&
       sqrt(sums(3) / sums(4)), maxs(2), sums(5) / max(sums(6), smalldouble), maxs(3),&
       sums(7) / sums(8), maxs(4), -maxs(5), maxs(6), sqrt(sums(9) / sums(8)), &
       sums(10) / max(sums(11), smalldouble)
    flush(my_unit)

  end subroutine write_log

  !> div B, |j x B| and |j x B| / (|j| |B|)
  subroutine lorentz_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw+nwauxio)
    double precision             :: normconv(0:nw+nwauxio)
    double precision :: jxb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: jb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call lorentz(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw),x,jxb,jb,divb)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1) = &
       divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+2) = &
       jxb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+3) = &
       jxb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) / &
       max(jb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3), smalldouble)

  end subroutine lorentz_output

  subroutine lorentz_names(varnames)
    character(len=*) :: varnames
    varnames = 'divb jxb sinjb'
  end subroutine lorentz_names

end module mod_usr
