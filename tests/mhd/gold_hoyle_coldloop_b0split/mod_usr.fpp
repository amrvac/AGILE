#:mute
#:include "../../../src/mod_gpu_directives.fpp"
#:endmute

!> Gold-Hoyle "cold loop": the force-free flux rope of tests/mhd/gold_hoyle
!> in physical units, with optically thin radiative cooling, anisotropic
!> (hyperbolic) thermal conduction and a heating that balances the losses of
!> the initial state. A non-adiabatic stability test. With B0SPLIT the rope
!> is B0 and the evolved B1 starts at zero (tests/mhd/gold_hoyle_coldloop_b0split).
!>
!> Equilibrium (Gold & Hoyle 1960), with r in units of L:
!>
!>   B_z = B_c / (1 + alpha^2 r^2)      B_phi = B_c alpha r / (1 + alpha^2 r^2)
!>
!> uniform density rho0 = m_p n0 and uniform pressure p0 = beta_c B_c^2 / 2
!> (Gaussian units, B^2/8pi magnetic pressure). Force-free, so the thermal
!> state can be anything; with cooling and conduction on, what is under test
!> is whether the thermal (non-adiabatic) modes stay damped.
!>
!> Physical inputs (&usr_list): n0 [cm^-3], B_c [G], L [cm], beta_c, alpha.
!> Code units are chosen so that rho0 = B_c = L = 1 and the Alfven speed
!> v_A = B_c / sqrt(4 pi rho0) is 1, hence tau_A = L / v_A = 1 as well:
!>
!>   unit_length = L,  unit_numberdensity = n0,
!>   unit_temperature = B_c^2 / (8 pi n0 k_B)   (needs He_abundance = 0)
!>
!> The MHD module's gas law is p = 2 n k_B T (protons + electrons, R = 1 in
!> code units), so for the given p0 and rho0 the code's temperature is
!> T0 = p0 m_p / (2 k_B rho0) - half of what a single-species gas (p = n k T)
!> would call it. print_setup reports both, and the cooling and conduction
!> coefficients are evaluated at the code's T0.
!>
!> Heating: a rate per unit mass, fixed once and for all so that it cancels
!> the radiative losses of the initial uniform state exactly,
!> H = rho/rho0 * rho0^2 Lambda(T0). Evaluated with the cooling module's own
!> exact-integration routine on the background state, so the cancellation at
!> t = 0 is to round-off whatever the time step.
!>
!> Log (typefilelog = 'special'): the equilibrium columns of
!> tests/mhd/gold_hoyle plus Tmin and Tmax in code units.
module mod_usr
  use mod_amrvac
  use mod_physics
  use mod_gimli

  implicit none

  ! ---- physical inputs, &usr_list ----
  !> number density [cm^-3]
  double precision :: n0_cgs = 1.0d9
  !> axial field on the axis [G]
  double precision :: bc_gauss = 10.0d0
  !> length unit L [cm]; the wall is at r = L
  double precision :: l_cgs = 1.0d10
  !> plasma beta on the axis, 2 p0 / B_c^2
  double precision :: beta_c = 0.002d0
  !> twist per unit length alpha [1/L]; the core radius is 1/alpha
  double precision :: alpha = 20.0d0
  !> switch the background heating off (then the loop just cools)
  logical :: heating = .true.
  !> Gaussian hot column at (x, y) = (spot_r, 0), all z: T = T0 (1 + spot_amp exp(-d^2/spot_sigma^2))
  double precision :: spot_amp = 0.0d0
  double precision :: spot_r = 0.15d0
  double precision :: spot_sigma = 0.05d0
  !> Legolas eigenmode from GIMLI added to the initial state
  character(len=std_len) :: ldat_file = ''

  ! ---- code-unit state ----
  double precision :: rho0 = 1.0d0
  double precision :: b0 = 1.0d0
  double precision :: qtwist = 20.0d0
  double precision :: p0 = 0.001d0
  !> conserved background state the heating is balanced against
  double precision :: w_bg(20) = 0.0d0
  ${GPU_DECLARE_COPYIN('rho0, b0, qtwist, p0, heating, w_bg, spot_amp, spot_r, spot_sigma')}$

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_constants, only: mp_cgs, kB_cgs

    call params_read_user(par_files)
    if (ldat_file /= '') call read_legolas_data(ldat_file, 21)

    ! code units: rho0 = B_c = L = v_A = tau_A = 1, see the module header
    unit_length        = l_cgs
    unit_numberdensity = n0_cgs
    unit_temperature   = bc_gauss**2 / (8.0d0 * dpi * n0_cgs * kB_cgs)

    rho0   = 1.0d0
    b0     = 1.0d0
    qtwist = alpha
    p0     = 0.5d0 * beta_c * b0**2
    ${GPU_UPDATE_DEVICE('rho0, b0, qtwist, p0, heating, spot_amp, spot_r, spot_sigma')}$

    nwauxio = 3

    usr_init_one_grid    => initonegrid_usr
    usr_special_bc       => specialbound_usr
    usr_before_main_loop => print_setup
    usr_print_log        => gh_print_log
    usr_aux_output       => gh_aux_output
    usr_add_aux_names    => gh_aux_names

    call phys_activate()

    ! the background state, once the variable indices exist. B is left out:
    ! the cooling only needs rho and the thermal energy
    w_bg(:)    = 0.0d0
    w_bg(iw_rho) = rho0
    w_bg(iw_e)   = p0 * inv_gamma_1
    ${GPU_UPDATE_DEVICE('w_bg')}$

  end subroutine usr_init

  subroutine params_read_user(files)
    use mod_global_parameters, only: unitpar
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ n0_cgs, bc_gauss, l_cgs, beta_c, alpha, heating, spot_amp, spot_r, spot_sigma, &
       ldat_file

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read_user

  !> The Gold-Hoyle field at x, in the components of the mesh: (B_x, B_y,
  !> B_z) on the Cartesian mesh, (B_r, B_z, B_phi) on the cylindrical one.
  pure subroutine gold_hoyle_field(x, b)
#ifdef _OPENMP
    ${GPU_ROUTINE_VECTOR()}$
#else
    ${GPU_ROUTINE_SEQ()}$
#endif
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: b(1:3)
    double precision              :: f

#:if GEOM == 'cylindrical'
    f = b0 / (1.0d0 + (qtwist * x(1))**2)
    b(1) = 0.0d0
    b(2) = f
    b(3) = f * qtwist * x(1)
#:else
    f = b0 / (1.0d0 + qtwist**2 * (x(1)**2 + x(2)**2))
    b(1) = -f * qtwist * x(2)
    b(2) =  f * qtwist * x(1)
    b(3) =  f
#:endif

  end subroutine gold_hoyle_field

  pure double precision function radius(x)
    double precision, intent(in) :: x(1:ndim)
#:if GEOM == 'cylindrical'
    radius = abs(x(1))
#:else
    radius = sqrt(x(1)**2 + x(2)**2)
#:endif
  end function radius

  !> squared distance across the rope from the spot at (r, phi) = (spot_r, 0)
  pure double precision function spot_distance2(x)
    double precision, intent(in) :: x(1:ndim)
#:if GEOM == 'cylindrical'
    spot_distance2 = x(1)**2 + spot_r**2 - 2.0d0 * x(1) * spot_r * cos(x(3))
#:else
    spot_distance2 = (x(1) - spot_r)**2 + x(2)**2
#:endif
  end function spot_distance2

  !> curl of the Gold-Hoyle field: alpha B with alpha = 2q/(1+q^2 r^2)
  pure subroutine gold_hoyle_current(x, j)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: j(1:3)
    double precision              :: b(1:3), r2

    call gold_hoyle_field(x, b)
#:if GEOM == 'cylindrical'
    r2 = x(1)**2
#:else
    r2 = x(1)**2 + x(2)**2
#:endif
    j(1:3) = 2.0d0 * qtwist / (1.0d0 + qtwist**2 * r2) * b(1:3)

  end subroutine gold_hoyle_current

  pure subroutine usr_b0_field(x, b0f)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: b0f(3)
    call gold_hoyle_field(x, b0f)
  end subroutine usr_b0_field

  pure subroutine usr_j0_field(x, j0f)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: j0f(3)
    call gold_hoyle_current(x, j0f)
  end subroutine usr_j0_field

  !> the evolved field at t = 0: the rope, or nothing on top of B0
  pure subroutine evolved_field(x, b)
#ifdef _OPENMP
    ${GPU_ROUTINE_VECTOR()}$
#else
    ${GPU_ROUTINE_SEQ()}$
#endif
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: b(1:3)
#:if defined('B0SPLIT')
    b(1:3) = 0.0d0
#:else
    call gold_hoyle_field(x, b)
#:endif
  end subroutine evolved_field

  !> B0 and J0 where the code stores B1, zero otherwise
  pure subroutine background(x, b, j)
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: b(1:3), j(1:3)
#:if defined('B0SPLIT')
    call gold_hoyle_field(x, b)
    call gold_hoyle_current(x, j)
#:else
    b(1:3) = 0.0d0
    j(1:3) = 0.0d0
#:endif
  end subroutine background

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    ! .. local ..
    double precision                :: b(1:3)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3, iw

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw) = 0.0d0
    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call evolved_field(x_loc, b)
             w(ix1,ix2,ix3,rho_)   = rho0
             w(ix1,ix2,ix3,p_)     = p0 * (1.0d0 + spot_amp &
                * exp(-spot_distance2(x_loc) / spot_sigma**2))
             w(ix1,ix2,ix3,mag(1)) = b(1)
             w(ix1,ix2,ix3,mag(2)) = b(2)
             w(ix1,ix2,ix3,mag(3)) = b(3)
          end do
       end do
    end do

    if (ldat_file /= '') then
       do iw = 1, nw
          call add_perturbation_to_w_array(ixGmin1,ixGmin2,ixGmin3,&
             ixGmax1,ixGmax2,ixGmax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,iw,x)
       end do
    end if

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  !> The wall: the exact initial state held in the ghost cells, at rest,
  !> at T0 and with no heat flux. psi and the heat-flux vector take the
  !> boundary type of rho, so they have to be set here as well.
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
    ! .. local ..
    double precision                :: b(1:3)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

#ifndef _OPENMP
    ${GPU_LOOP_VECTOR("collapse(3) private(b, x_loc)")}$
#endif
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call evolved_field(x_loc, b)

             w(ix1,ix2,ix3,iw_rho)    = rho0
             w(ix1,ix2,ix3,iw_mom(1)) = 0.0d0
             w(ix1,ix2,ix3,iw_mom(2)) = 0.0d0
             w(ix1,ix2,ix3,iw_mom(3)) = 0.0d0
             w(ix1,ix2,ix3,iw_mag(1)) = b(1)
             w(ix1,ix2,ix3,iw_mag(2)) = b(2)
             w(ix1,ix2,ix3,iw_mag(3)) = b(3)
             w(ix1,ix2,ix3,psi_)      = 0.0d0
             w(ix1,ix2,ix3,iw_e)      = p0 * inv_gamma_1 + 0.5d0 * sum(b(1:3)**2)
#:if defined('HYPERTC_ANISO')
             w(ix1,ix2,ix3,q_(1))     = 0.0d0
             w(ix1,ix2,ix3,q_(2))     = 0.0d0
             w(ix1,ix2,ix3,q_(3))     = 0.0d0
#:elif defined('HYPERTC')
             w(ix1,ix2,ix3,q_)        = 0.0d0
#:endif

          end do
       end do
    end do

  end subroutine specialbound_usr

#:if defined('SOURCE_USR')
  !> Heating per unit mass that balances the radiative losses of the initial
  !> state: rho/rho0 times the loss rate of the uniform background, computed
  !> with the same exact-integration routine the cooling itself uses, so the
  !> two cancel to round-off at t = 0 for any time step. Must be called
  !> addsource_usr; runs on the device.
  subroutine addsource_usr(qdt, qt, wCT, wCTprim, wnew, x, split)
#:if defined('COOLING')
    use mod_radiative_cooling, only: rc_fl, getvar_cooling_exact
#:endif
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)    :: qdt, qt, wCT(nw_phys), wCTprim(nw_phys)
    double precision, intent(inout) :: wnew(nw_phys)
    double precision, intent(in)    :: x(1:ndim)
    logical, intent(in)             :: split
    ! .. local ..
    double precision :: coolrate

#:if defined('COOLING')
    if (heating) then
       call getvar_cooling_exact(qdt, w_bg(1:nw_phys), w_bg(1:nw_phys), x,&
          coolrate, rc_fl)
       wnew(iw_e) = wnew(iw_e) + qdt * coolrate * wCT(iw_rho) / rho0
    end if
#:endif

  end subroutine addsource_usr
#:endif

  !> The scales of the run, once the par file and the physics are set up.
  subroutine print_setup()
    use mod_global_parameters
    use mod_constants, only: mp_cgs, kB_cgs
    use mod_comm_lib, only: mpistop
#:if defined('COOLING')
    use mod_radiative_cooling, only: rc_fl, findL
#:endif
    double precision :: t0, lz, twist, lam, tcool, tcond

    if (He_abundance /= 0.0d0) call mpistop(&
       "gold_hoyle_coldloop: the unit choice assumes He_abundance = 0")

    t0 = p0 / rho0
#:if GEOM == 'cylindrical'
    lz = xprobmax2 - xprobmin2
#:else
    lz = xprobmax3 - xprobmin3
#:endif
    twist = qtwist * lz

    if (mype == 0) then
       write(*,'(A)')         ' ========================================'
       write(*,'(A)')         '  Gold-Hoyle cold loop'
       write(*,'(A,ES12.5,A)') '    L                    = ', unit_length, ' cm'
       write(*,'(A,ES12.5,A)') '    n0                   = ', unit_numberdensity, ' cm^-3'
       write(*,'(A,ES12.5,A)') '    rho0                 = ', unit_density, ' g cm^-3'
       write(*,'(A,ES12.5,A)') '    B_c                  = ', bc_gauss, ' G'
       write(*,'(A,ES12.5,A)') '    p0                   = ', p0 * unit_pressure, ' dyn cm^-2'
       write(*,'(A,ES12.5)')   '    beta_c               = ', beta_c
       write(*,'(A,ES12.5,A)') '    T0 (code gas law)    = ', t0 * unit_temperature, ' K'
       write(*,'(A,ES12.5,A)') '    T0 (p = n k T)       = ', 2.0d0 * t0 * unit_temperature, ' K'
       write(*,'(A,ES12.5,A)') '    v_A                  = ', unit_velocity, ' cm s^-1'
       write(*,'(A,ES12.5,A)') '    tau_A = L / v_A      = ', unit_time, ' s'
       write(*,'(A,ES12.5)')   '    alpha (twist/L)      = ', qtwist
       write(*,'(A,ES12.5)')   '    core radius 1/alpha  = ', 1.0d0 / qtwist
       write(*,'(A,ES12.5)')   '    rope length L_z / L  = ', lz
       write(*,'(A,ES12.5)')   '    total twist / 2pi    = ', twist / (2.0d0 * dpi)
       write(*,'(A,ES12.5)')   '    time_max / tau_A     = ', time_max
       write(*,'(A,2ES12.5)')  '    code: p0, T0 (rho0 = B0 = 1) = ', p0, t0
#:if defined('COOLING')
       call findL(t0, lam, rc_fl)
       tcool = p0 * inv_gamma_1 / (rho0**2 * lam)
       write(*,'(A,ES12.5)')   '    rho0^2 Lambda(T0)    = ', rho0**2 * lam
       write(*,'(A,ES12.5,A)') '    tau_cool / tau_A     = ', tcool, &
          '  (thermal energy / radiative loss)'
       write(*,'(A,L2)')       '    heating balance      = ', heating
#:endif
#:if defined('HYPERTC')
       tcond = rho0 * inv_gamma_1 / (qtwist**2 * tc_kappa0_par * t0**2.5d0)
       write(*,'(A,ES12.5)')   '    kappa0_par (code)    = ', tc_kappa0_par
       write(*,'(A,ES12.5,A)') '    tau_cond / tau_A     = ', tcond, &
          '  (Spitzer, across one core radius)'
#:endif
#:if defined('B0SPLIT')
       write(*,'(A)')          '    B0 split: B1 = 0 at t = 0'
#:endif
       if (twist > 2.0d0 * dpi) then
          write(*,'(A)')      '  note: total twist exceeds 2 pi. The unbounded'
          write(*,'(A)')      '  uniform-twist field is ideal-MHD stable (Voslamber'
          write(*,'(A)')      '  & Callebaut 1962); with a wall the m=1 stability'
          write(*,'(A)')      '  has to be checked case by case.'
       end if
       write(*,'(A)')         ' ========================================'
    end if

  end subroutine print_setup

  !> Per-cell force-free diagnostics of the conserved state w of one block:
  !> |j x B|, |j| |B| and div B on the interior.
  subroutine gh_lorentz(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
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
    ! .. local ..
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)
    double precision :: j(1:3), b(1:3), f(1:3), b0f(1:3), j0f(1:3), x_loc(1:ndim)
    integer          :: idirmin, ix1, ix2, ix3

    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)

    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call background(x_loc, b0f, j0f)
             j(1:3) = current(ix1,ix2,ix3,1:3) + j0f(1:3)
             b(1:3) = w(ix1,ix2,ix3,iw_mag(1:3)) + b0f(1:3)
             f(1) = j(2)*b(3) - j(3)*b(2)
             f(2) = j(3)*b(1) - j(1)*b(3)
             f(3) = j(1)*b(2) - j(2)*b(1)
             jxb(ix1,ix2,ix3) = sqrt(sum(f(1:3)**2))
             jb(ix1,ix2,ix3)  = sqrt(sum(j(1:3)**2)) * sqrt(sum(b(1:3)**2))
          end do
       end do
    end do

  end subroutine gh_lorentz

  !> The log. One line per log time:
  !>
  !>   it t dt Ekin Emag vmax dB_L2 drho_max sigmaJ jxb_max divb_mean divb_max Tmin Tmax sig_r2
  !>
  !> sig_r2 is the radial variance of the temperature excess about spot_r; its
  !> growth rate is twice the cross-field diffusivity.
  !>
  !> as in tests/mhd/gold_hoyle, plus the extreme temperatures in code units
  !> (multiply by unit_temperature for K). saveamrfile has refreshed the
  !> host copy of w and filled the ghost cells before calling this.
  subroutine gh_print_log()
    use mod_global_parameters
    ! .. local ..
    integer, parameter :: nsum = 10, nmax = 5
    double precision   :: jxb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision   :: jb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision   :: divb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision   :: sums_ipe(nsum), sums(nsum), maxs_ipe(nmax), maxs(nmax)
    double precision   :: dv, rho, v2, b(1:3), b0f(1:3), j0f(1:3), bex(1:3), x_loc(1:ndim), temp
    double precision   :: ekin, emag, vmax, db_l2, drho_max, sigmaj, jxb_max
    double precision   :: divb_mean, divb_max, tmin, tmax, sig_r2, dtemp, rc
    integer            :: iigrid, igrid, ix1, ix2, ix3
    integer, parameter :: my_unit = 21
    logical, save      :: opened = .false.
    character(len=80)  :: filename

    sums_ipe = 0.0d0
    maxs_ipe = 0.0d0
    maxs_ipe(5) = -bigdouble    ! carries -Tmin, so a MAX reduction gives the min

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       block => ps(igrid)
       dxlevel(1) = rnode(rpdx1_,igrid)
       dxlevel(2) = rnode(rpdx2_,igrid)
       dxlevel(3) = rnode(rpdx3_,igrid)

       call gh_lorentz(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid)%w,ps(igrid)%x,jxb,jb,divb)

       do ix3 = ixMlo3, ixMhi3
          do ix2 = ixMlo2, ixMhi2
             do ix1 = ixMlo1, ixMhi1
                if (slab_uniform) then
                   dv = dxlevel(1) * dxlevel(2) * dxlevel(3)
                else
                   dv = block%dvolume(ix1,ix2,ix3)
                end if

                x_loc(1:ndim) = ps(igrid)%x(ix1,ix2,ix3,1:ndim)
                call background(x_loc, b0f, j0f)
                call gold_hoyle_field(x_loc, bex)
                rho    = ps(igrid)%w(ix1,ix2,ix3,iw_rho)
                v2     = sum(ps(igrid)%w(ix1,ix2,ix3,iw_mom(1:3))**2) / rho**2
                temp   = gamma_1 * (ps(igrid)%w(ix1,ix2,ix3,iw_e) - 0.5d0 * rho * v2 &
                         - 0.5d0 * sum(ps(igrid)%w(ix1,ix2,ix3,iw_mag(1:3))**2)) / rho
                b(1:3) = ps(igrid)%w(ix1,ix2,ix3,iw_mag(1:3)) + b0f(1:3)

                sums_ipe(1) = sums_ipe(1) + 0.5d0 * rho * v2 * dv
                sums_ipe(2) = sums_ipe(2) + 0.5d0 * sum(b(1:3)**2) * dv
                sums_ipe(3) = sums_ipe(3) + sum((b(1:3) - bex(1:3))**2) * dv
                sums_ipe(4) = sums_ipe(4) + sum(bex(1:3)**2) * dv
                sums_ipe(5) = sums_ipe(5) + jxb(ix1,ix2,ix3) * dv
                sums_ipe(6) = sums_ipe(6) + jb(ix1,ix2,ix3) * dv
                sums_ipe(7) = sums_ipe(7) + abs(divb(ix1,ix2,ix3)) * dv
                sums_ipe(8) = sums_ipe(8) + dv
                dtemp = max(temp - p0 / rho0, 0.0d0)
                rc    = radius(x_loc)
                sums_ipe(9)  = sums_ipe(9)  + dtemp * dv
                sums_ipe(10) = sums_ipe(10) + dtemp * (rc - spot_r)**2 * dv

                maxs_ipe(1) = max(maxs_ipe(1), sqrt(v2))
                maxs_ipe(2) = max(maxs_ipe(2), abs(rho - rho0) / rho0)
                maxs_ipe(3) = max(maxs_ipe(3), jxb(ix1,ix2,ix3))
                maxs_ipe(4) = max(maxs_ipe(4), abs(divb(ix1,ix2,ix3)))
                maxs_ipe(5) = max(maxs_ipe(5), -temp)
             end do
          end do
       end do
    end do

    call MPI_ALLREDUCE(sums_ipe, sums, nsum, MPI_DOUBLE_PRECISION, MPI_SUM,&
       icomm, ierrmpi)
    call MPI_ALLREDUCE(maxs_ipe, maxs, nmax, MPI_DOUBLE_PRECISION, MPI_MAX,&
       icomm, ierrmpi)
    call max_temperature(tmax)

    ekin      = sums(1)
    emag      = sums(2)
    db_l2     = sqrt(sums(3) / sums(4))
    sigmaj    = sums(5) / max(sums(6), smalldouble)
    divb_mean = sums(7) / sums(8)
    vmax      = maxs(1)
    drho_max  = maxs(2)
    jxb_max   = maxs(3)
    divb_max  = maxs(4)
    tmin      = -maxs(5)
    sig_r2    = sums(10) / max(sums(9), smalldouble)

    if (mype == 0) then
       if (.not. opened) then
          filename = trim(base_filename) // ".log"
          if (restart_from_file == undefined) then
             open(unit=my_unit, file=filename, status='replace')
             write(my_unit,'(A)') '# it t dt Ekin Emag vmax dB_L2 drho_max ' // &
                'sigmaJ jxb_max divb_mean divb_max Tmin Tmax sig_r2'
          else
             open(unit=my_unit, file=filename, status='old', position='append')
          end if
          opened = .true.
       end if
       write(my_unit,'(I8,14ES16.8)') it, global_time, dt, ekin, emag, vmax,&
          db_l2, drho_max, sigmaj, jxb_max, divb_mean, divb_max, tmin, tmax, sig_r2
       flush(my_unit)
    end if

  end subroutine gh_print_log

  !> the largest temperature over the domain, in code units
  subroutine max_temperature(tmax)
    use mod_global_parameters
    double precision, intent(out) :: tmax
    double precision :: tmax_ipe, rho, v2, b2, temp
    integer          :: iigrid, igrid, ix1, ix2, ix3

    tmax_ipe = -bigdouble
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       do ix3 = ixMlo3, ixMhi3
          do ix2 = ixMlo2, ixMhi2
             do ix1 = ixMlo1, ixMhi1
                rho  = ps(igrid)%w(ix1,ix2,ix3,iw_rho)
                v2   = sum(ps(igrid)%w(ix1,ix2,ix3,iw_mom(1:3))**2) / rho**2
                b2   = sum(ps(igrid)%w(ix1,ix2,ix3,iw_mag(1:3))**2)
                temp = gamma_1 * (ps(igrid)%w(ix1,ix2,ix3,iw_e) &
                       - 0.5d0 * rho * v2 - 0.5d0 * b2) / rho
                tmax_ipe = max(tmax_ipe, temp)
             end do
          end do
       end do
    end do
    call MPI_ALLREDUCE(tmax_ipe, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX,&
       icomm, ierrmpi)

  end subroutine max_temperature

  !> div B, |j x B| and |j x B| / (|j| |B|) for the converted output
  subroutine gh_aux_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw+nwauxio)
    double precision             :: normconv(0:nw+nwauxio)
    ! .. local ..
    double precision :: jxb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: jb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call gh_lorentz(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw),x,jxb,jb,divb)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1) = &
       divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+2) = &
       jxb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+3) = &
       jxb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) / &
       max(jb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3), smalldouble)

  end subroutine gh_aux_output

  subroutine gh_aux_names(varnames)
    character(len=*) :: varnames
    varnames = 'divb jxb sinjb'
  end subroutine gh_aux_names

end module mod_usr
