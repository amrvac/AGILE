#:if PHYS == 'mhd'

#:if defined('N_TRACER')
#:set N_TRACER_ = N_TRACER
#:else
#:set N_TRACER_ = 0
#:endif

#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter, public              :: nw_phys=2+2*ndim+1+${N_TRACER_}$ &
#:if defined('HYPERTC_ANISO')
    + 2
#:elif defined('HYPERTC')
    + 1
#:else
    + 0
#:endif

  integer, parameter, public              :: nw_flux=2+2*ndim+1+${N_TRACER_}$ &
#:if defined('HYPERTC_ANISO')
    + 2
#:elif defined('HYPERTC')
    + 1
#:else
    + 0
#:endif

  !> Whether an energy equation is used
  logical, public                         :: mhd_energy = .true.
  !$acc declare copyin(mhd_energy)

  !> Index of the density (in the w array)
  integer, public                         :: rho_
  !$acc declare create(rho_)

  !> Indices of the momentum density
  integer, allocatable, public            :: mom(:)
  !$acc declare create(mom)

#:if defined('N_TRACER')
  !> Indices of the tracers
  integer, public                         :: tracer(${N_TRACER_}$)
  !$acc declare create(tracer)
#:endif

  !> Index of the energy density (-1 if not present)
  integer, public                         :: e_
  !$acc declare create(e_)

  !> Indices of the magnetic field
  integer, allocatable, public            :: mag(:)
  !$acc declare create(mag)

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public                         :: p_
  !$acc declare create(p_)

  !> Indices of the GLM psi
  integer, public :: psi_
  !$acc declare create(psi_)

  !> Number of tracer species
  integer, public                         :: mhd_n_tracer = 0
  !$acc declare copyin(mhd_n_tracer)

  !> The adiabatic index
  double precision, public                :: mhd_gamma = 5.d0/3.0d0
  !$acc declare copyin(mhd_gamma)

  !> The adiabatic index minus 1
  double precision, public                :: mhd_gamma_m1
  !$acc declare copyin(mhd_gamma_m1)

  !> Helium abundance over Hydrogen
  double precision, public  :: He_abundance=0.1d0
  !$acc declare copyin(He_abundance)
  !> Ionization fraction of H
  !> H_ion_fr = H+/(H+ + H)
  double precision, public  :: H_ion_fr=1d0
  !$acc declare copyin(H_ion_fr)
  !> Ionization fraction of He
  !> He_ion_fr = (He2+ + He+)/(He2+ + He+ + He)
  double precision, public  :: He_ion_fr=1d0
  !$acc declare copyin(He_ion_fr)
  !> Ratio of number He2+ / number He+ + He2+
  !> He_ion_fr2 = He2+/(He2+ + He+)
  double precision, public  :: He_ion_fr2=1d0
  !$acc declare copyin(He_ion_fr2)
  ! used for eq of state when it is not defined by units,
  ! the units do not contain terms related to ionization fraction
  ! and it is p = RR * rho * T
  double precision, public  :: RR=1d0
  !$acc declare copyin(RR)
  !> Index of field-aligned heat flux scalar q_; for HYPERTC_ANISO: qperp_ = q_+1
#:if defined('HYPERTC')
  integer, public                         :: q_ = -1
  !$acc declare create(q_)
#:endif
#:if defined('HYPERTC_ANISO')
  integer, public                         :: qperp_ = -1
  !$acc declare create(qperp_)
#:endif

#:if defined('HYPERTC') or defined('HYPERTC_ANISO')
  !> sig_par = tc_kappa0_par * Te^2.5; if <= 0 and tc_kappa_par <= 0, set from Spitzer in phys_units()
  double precision, public                :: tc_kappa0_par = -1.0d0
  !$acc declare copyin(tc_kappa0_par)
  !> sig_par = tc_kappa_par (constant, no T^2.5); takes precedence over tc_kappa0_par when > 0
  double precision, public                :: tc_kappa_par = -1.0d0
  !$acc declare copyin(tc_kappa_par)
#:endif
#:if defined('HYPERTC_ANISO')
  !> Magnetisation chi prefactor: chi = htc_Cchi * B * Te^1.5 / n
  double precision, public                :: htc_Cchi = 0.0d0
  !$acc declare copyin(htc_Cchi)
  !> sig_perp = tc_kappa0_perp * Te^2.5; if <= 0 and tc_kappa_perp <= 0, uses Braginskii sig_par/(1+chi^2)
  double precision, public                :: tc_kappa0_perp = -1.0d0
  !$acc declare copyin(tc_kappa0_perp)
  !> sig_perp = tc_kappa_perp (constant, no T^2.5); takes precedence over tc_kappa0_perp when > 0
  double precision, public                :: tc_kappa_perp = -1.0d0
  !$acc declare copyin(tc_kappa_perp)
#:endif

  !> Switch for hyperbolic thermal conduction
  logical, public                         :: mhd_hyperbolic_thermal_conduction = .false.
  !$acc declare copyin(mhd_hyperbolic_thermal_conduction)

  !> Freeze rho/v/B and evolve only the energy equation (compile-time MHD_ENERGY_ONLY flag)
  logical, public                         :: mhd_energy_only = .false.

  !> Compile-time selector for anisotropic HTC
  logical, public                         :: mhd_hyperbolic_thermal_conduction_anisotropic = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: mhd_glm_alpha = 0.5d0
  !$acc declare copyin(mhd_glm_alpha)

  !> Whether to use gravity
  logical, public                         :: mhd_gravity = .false.
  !$acc declare copyin(mhd_gravity)

  !> The resistivity
  double precision, public                :: mhd_eta = 0.0d0
  !$acc declare copyin(mhd_eta)

  !> switch for adding resistive terms
  logical, public                         :: mhd_resistivity = .false.
  !$acc declare copyin(mhd_resistivity)

  !> Whether plasma is partially ionized
  logical, public                         :: mhd_partial_ionization = .false.
  !$acc declare copyin(mhd_partial_ionization)

  !> switch for radiative cooling
  logical, public                         :: mhd_radiative_cooling = .false.
  !$acc declare copyin(mhd_radiative_cooling)

  !> Whether particles module is added
  logical, public                         :: mhd_particles = .false.
  !$acc declare copyin(mhd_particles)

  !> switch for source user
  logical, public                         :: mhd_source_usr = .false.
  !$acc declare copyin(mhd_source_usr)

#:enddef

#:def read_params()
    !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mhd_list/ mhd_energy, mhd_gamma, mhd_glm_alpha, mhd_gravity,&
      mhd_n_tracer, mhd_radiative_cooling, He_abundance, mhd_eta, mhd_source_usr, &
      mhd_resistivity, mhd_hyperbolic_thermal_conduction, &
      mhd_hyperbolic_thermal_conduction_anisotropic, &
#:if defined('HYPERTC')
      tc_kappa_par, &
#:endif
#:if defined('HYPERTC_ANISO')
      tc_kappa_perp, &
#:endif
      mhd_energy_only

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mhd_list, end=111)
111    close(unitpar)
    end do

#ifdef _OPENACC
    !$acc update device(mhd_energy, &
    !$acc&     mhd_gamma, mhd_glm_alpha, &
    !$acc&     mhd_gravity, mhd_n_tracer, mhd_radiative_cooling, &
    !$acc&     He_abundance, mhd_eta, mhd_source_usr, mhd_resistivity)
#:if defined('HYPERTC')
    !$acc update device(tc_kappa_par)
#:endif
#:if defined('HYPERTC_ANISO')
    !$acc update device(tc_kappa_perp)
#:endif
#endif

  end subroutine read_params
#:enddef

#:def phys_activate() 
  subroutine phys_activate()
    call phys_init()
  end subroutine phys_activate
#:enddef

#:def phys_units()
  subroutine phys_units()
    use mod_global_parameters
    double precision :: mp, kB
    double precision :: a,b

    !> here no SI_UNIT used by default, to be implemented
    mp = mp_cgs
    kB = kB_cgs
    !> eq_state_units by default, to be implemented
    a = 1.d0+4.d0*He_abundance
    b = 2.d0+3.d0*He_abundance

    if(unit_density/=1.d0 .or. unit_numberdensity/=1.d0) then
      if(unit_density/=1.d0) then
        unit_numberdensity=unit_density/(a*mp)
      else if(unit_numberdensity/=1.d0) then
        unit_density=a*mp*unit_numberdensity
      end if
      if(unit_temperature/=1.d0) then
        unit_pressure=b*unit_numberdensity*kB*unit_temperature
        unit_velocity=dsqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_pressure/=1.d0) then
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        unit_velocity=dsqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_velocity/=1.d0) then
        unit_pressure=unit_density*unit_velocity**2
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=1.d0) then
        unit_velocity=unit_length/unit_time
        unit_pressure=unit_density*unit_velocity**2
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    else if(unit_temperature/=1.d0) then
      ! units of temperature and velocity are dependent
      if(unit_pressure/=1.d0) then
        unit_numberdensity=unit_pressure/(b*unit_temperature*kB)
        unit_density=a*mp*unit_numberdensity
        unit_velocity=dsqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      end if
    else if(unit_pressure/=1.d0) then
      if(unit_velocity/=1.d0) then
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=0.d0) then
        unit_velocity=unit_length/unit_time
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    end if
    unit_mass=unit_density*unit_length**3

    !$acc update device(unit_density, unit_numberdensity, unit_temperature, unit_pressure, unit_velocity, unit_length, unit_time, unit_mass)

#:if defined('HYPERTC') or defined('HYPERTC_ANISO')
    if (tc_kappa0_par <= 0.0d0 .and. tc_kappa_par <= 0.0d0) &
      tc_kappa0_par = 8.0d-7 * unit_temperature**3.5d0 &
                    / (unit_length * unit_density * unit_velocity**3.0d0)
    !$acc update device(tc_kappa0_par)
    !$acc update device(tc_kappa_par)
#:endif
#:if defined('HYPERTC_ANISO')
    htc_Cchi = 0.823d0 * (4.753567596681522d6 / 20.0d0) &
             * unit_magneticfield * unit_temperature**1.5d0 / unit_numberdensity
    !$acc update device(htc_Cchi)
    !$acc update device(tc_kappa0_perp)
    !$acc update device(tc_kappa_perp)
#:endif
  end subroutine phys_units
#:enddef
  
#:def phys_init()
    !> Initialize the module
  subroutine phys_init()
    use mod_global_parameters
#:if defined('COOLING')
    use mod_radiative_cooling, only: rc_fl, radiative_cooling_init_params, radiative_cooling_init
#:endif

    integer      :: idir, idum

    call read_params(par_files)
    call phys_units()

    phys_energy  = mhd_energy
    phys_total_energy  = mhd_energy
    phys_internal_e = .false.
    phys_gamma = mhd_gamma
    phys_partial_ionization=mhd_partial_ionization
    need_global_cmax=.true.
    mhd_gamma_m1=mhd_gamma-1.0_dp
 !$acc update device(physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization,need_global_cmax,mhd_gamma_m1)

    use_particles = mhd_particles

    ! Determine flux variables
    rho_ = var_set_rho()
    !$acc update device(rho_)

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)
    !$acc update device(mom)

    ! Set index of energy variable
    if (mhd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if
    !$acc update device(e_,p_)

    ! Set index for heat flux variable(s)
#:if defined('HYPERTC_ANISO')
    q_     = var_set_q(need_bc=.false.)
    !$acc update device(q_)
    qperp_ = var_set_qperp(need_bc=.false.)
    !$acc update device(qperp_)
#:elif defined('HYPERTC')
    q_ = var_set_q(need_bc=.false.)
    !$acc update device(q_)
#:endif

    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)
    !$acc update device(mag)

    psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    !$acc update device(psi_)

    !> GLM MHD uses split source addition in psi:
    any_source_split = .true.
    !$acc update device(any_source_split)

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .true.

    ! Register tracer fields
#:if defined('N_TRACER')
    #:for i in range(1, N_TRACER_+1)
        tracer(${i}$) = var_set_fluxvar("trc", "trp", ${i}$, need_bc=.false.)
    #:endfor
    !$acc update device(tracer)
#:endif

    ! set number of variables which need update ghostcells
    nwgc=nwflux
    !$acc update device(nwgc)

    ! Define custom flux types:
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw_flux))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw_flux])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    ! BnormLF fix:
    flux_type(:,psi_)=flux_tvdlf
    do idir=1,ndir
       flux_type(idir,mag(idir))=flux_tvdlf
    end do
    !$acc update device(flux_type)

    
! use cycle, needs to be dealt with:    
!    ! Initialize particles module
!    if (mhd_particles) then
!       call particles_init()
!       phys_req_diagonal = .true.
!    end if

#:if defined('COOLING')
    call radiative_cooling_init_params(phys_gamma,He_abundance)
    call radiative_cooling_init(rc_fl)
    !$acc update device(rc_fl)
    !$acc enter data copyin(rc_fl%tcool,rc_fl%Lcool, rc_fl%Yc)
#:endif

  end subroutine phys_init
#:enddef

#:def phys_get_dt()
  subroutine phys_get_dt(w, x, dx, dtnew)
  !$acc routine seq
#:if defined('GRAVITY')
  use mod_usr, only: gravity_field
#:endif    
#:if defined('RESISTIVE')
  use mod_global_parameters, only: dtdiffpar
#:endif    
    real(dp), intent(in)   :: w(nw_phys), x(1:ndim), dx(1:ndim)
    real(dp), intent(out)  :: dtnew
    ! .. local ..
    integer                :: idim
    real(dp)               :: field

    dtnew = huge(1.0d0)
    
#:if defined('GRAVITY')
    do idim = 1, ndim
       field = gravity_field(w, x, idim)
       field = max( abs(field), epsilon(1.0d0) )
       dtnew = min( dtnew, 1_dp / sqrt( field/dx(idim) ) )
    end do
#:endif    
    
#:if defined('RESISTIVE')
    do idim = 1,ndim
       dtnew=min(dtnew, dtdiffpar*dx(idim)**2/mhd_eta)
    enddo
#:endif    

  end subroutine phys_get_dt
#:enddef  

#:def addsource_local()
subroutine addsource_local(qdt, dtfactor, qtC, wCT, wCTprim, qt, wnew, x, dr, &
    qsourcesplit)
  !$acc routine seq
#:if defined('SOURCE_USR')
  use mod_usr, only: addsource_usr
#:endif
#:if defined('GRAVITY')
  use mod_usr, only: gravity_field
#:endif
#:if defined('COOLING')
  use mod_radiative_cooling, only: rc_fl, radiative_cooling_add_source
#:endif

  use mod_global_parameters, only:cmax_global
  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCT(nw_phys), wCTprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim), dr(ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit
  ! .. local ..
  integer                  :: idim
  real(dp)                 :: field

  if (.not. qsourcesplit) then 
     !---------------------------------
     ! unsplit sources
     !---------------------------------

#:if defined('GRAVITY')
     do idim = 1, ndim
        field = gravity_field(wCT, x, idim)
        wnew(iw_mom(idim)) = wnew(iw_mom(idim)) + qdt * field * wCT(iw_rho)
        wnew(iw_e)         = wnew(iw_e) + qdt * field * wCT(iw_mom(idim))
     end do
#:endif

#:if defined('COOLING')
     call radiative_cooling_add_source(qdt,wCT,wCTprim,wnew,x)
#:endif

#:if defined('SOURCE_USR')
     call addsource_usr(qdt, qt, wCT, wCTprim, wnew, x, .false.)
#:endif

     
  else
     !---------------------------------
     ! split sources     
     !---------------------------------
     
     wnew(psi_)=wnew(psi_)*dexp(-qdt*cmax_global*mhd_glm_alpha/minval(dr))
     
  end if

end subroutine addsource_local
#:enddef

#:def addsource_compact()
subroutine addsource_compact(qdt, dtfactor, qtC, wCTprim1, wCTprim2, wCTprim3, qt, wnew, x, dx, &
     qsourcesplit)
  !$acc routine seq
#:if defined('HYPERTC')
  use mod_global_parameters, only: dt, cmax_global, ndir, smalldouble
#:endif

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCTprim1(nw_phys,3),wCTprim2(nw_phys,3),wCTprim3(nw_phys,3)
  real(dp), intent(in)     :: x(1:ndim), dx(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit
  ! .. local ..
  real(dp)                 :: laplb_cd2
  real(dp)                 :: Jdir1,Jdir2,Jdir3
  integer                  :: idir
#:if defined('HYPERTC')
  real(dp)   :: Te_c, rho_c, e_c, pth_c, sig_par
  real(dp)   :: gradT(3), bgradT
  real(dp)   :: Bmag2, Bmag, Bmag2_safe, tau_par
  real(dp)   :: q_sat, f_sat
  integer    :: k_tc
#:endif
#:if defined('HYPERTC_ANISO')
  real(dp)   :: sig_perp, chi, gradTperp_mag, tau_perp
#:endif

  if (.not. qsourcesplit) then 
     !---------------------------------
     ! unsplit sources
     !---------------------------------

#:if defined('RESISTIVE')
     ! using the compact stencil formulation and adopting constant eta

     ! > eta*(laplacian of B)_idir          added to B_idir
     ! > B_idir*[eta*(laplacian of B)_idir] added to e_
     do idir = 1,ndim
        laplb_cd2 =  &
             ( &
             wCTprim1(iw_mag(1)-1+idir,3) &
             - 2*wCTprim1(iw_mag(1)-1+idir,2) &
             + wCTprim1(iw_mag(1)-1+idir,1) &
             ) &
             / dx(1)**2 + &
             ( &
             wCTprim2(iw_mag(1)-1+idir,3) &
             - 2*wCTprim2(iw_mag(1)-1+idir,2) &
             + wCTprim2(iw_mag(1)-1+idir,1) &
             ) &
             / dx(2)**2 + &
             ( &
             wCTprim3(iw_mag(1)-1+idir,3) &
             - 2*wCTprim3(iw_mag(1)-1+idir,2) &
             + wCTprim3(iw_mag(1)-1+idir,1) &
             ) &
             / dx(3)**2 

        wnew(iw_mag(1)-1+idir) = wnew(iw_mag(1)-1+idir) + qdt*mhd_eta*laplb_cd2
        wnew(iw_e) = wnew(iw_e) + qdt * mhd_eta * laplb_cd2 * wCTprim1(iw_mag(1)-1+idir,2)
     enddo

     ! > eta* J**2 added to e
     Jdir1 = (wCTprim2(iw_mag(3),3) - wCTprim2(iw_mag(3),1)) &
          / 2.0_dp/dx(2) &
          - (wCTprim3(iw_mag(2),3) - wCTprim3(iw_mag(2),1)) &
          / 2.0_dp/dx(3)
     Jdir2 = (wCTprim3(iw_mag(1),3) - wCTprim3(iw_mag(1),1)) &
          /2.0_dp/dx(3) &
          - (wCTprim1(iw_mag(3),3) - wCTprim1(iw_mag(3),1)) &
          / 2.0_dp/dx(1)
     Jdir3 = (wCTprim1(iw_mag(2),3) - wCTprim1(iw_mag(2),1)) &
          /2.0_dp/dx(1) &
          - (wCTprim2(iw_mag(1),3) - wCTprim2(iw_mag(1),1)) &
          / 2.0_dp/dx(2)
     wnew(iw_e) = wnew(iw_e) + qdt*mhd_eta*(Jdir1**2+Jdir2**2+Jdir3**2)
#:endif

#:if defined('HYPERTC')
    Te_c  = wCTprim1(iw_e,2) / wCTprim1(iw_rho,2)
    rho_c = wCTprim1(iw_rho,2)
    pth_c = wCTprim1(iw_e,2)

    gradT(1) = (wCTprim1(iw_e,3)/wCTprim1(iw_rho,3) - wCTprim1(iw_e,1)/wCTprim1(iw_rho,1)) / (2.d0*dx(1))
    gradT(2) = (wCTprim2(iw_e,3)/wCTprim2(iw_rho,3) - wCTprim2(iw_e,1)/wCTprim2(iw_rho,1)) / (2.d0*dx(2))
    gradT(3) = (wCTprim3(iw_e,3)/wCTprim3(iw_rho,3) - wCTprim3(iw_e,1)/wCTprim3(iw_rho,1)) / (2.d0*dx(3))

    Bmag2 = 0.0d0
    do k_tc = 1, ndir
      Bmag2 = Bmag2 + wCTprim1(iw_mag(k_tc),2)**2
    end do
    Bmag       = sqrt(Bmag2)
    Bmag2_safe = max(Bmag2, smalldouble**2)

    bgradT = 0.0d0
    do k_tc = 1, ndir
      bgradT = bgradT + wCTprim1(iw_mag(k_tc),2) * gradT(k_tc)
    end do
    bgradT = bgradT * Bmag / Bmag2_safe

    if (tc_kappa_par > 0.0d0) then
      sig_par = tc_kappa_par
    else
      sig_par = tc_kappa0_par * Te_c**2.5d0
    end if

    e_c = pth_c / (mhd_gamma - 1.0d0) + 0.5d0 * Bmag2

    ! free-streaming limit: q_sat = 1.5 rho c_s^3, c_s = sqrt(p/rho)
    q_sat = 1.5d0 * rho_c * (pth_c / rho_c)**1.5d0
    f_sat = 1.0d0 / (1.0d0 + abs(sig_par * bgradT) / q_sat)
    tau_par = max(4.d0*dt, f_sat*sig_par*Te_c*(mhd_gamma-1.0d0) / (e_c*cmax_global**2))
    wnew(q_) = wnew(q_) - qdt*(f_sat*sig_par*bgradT + wCTprim1(q_,2))/tau_par

#:if defined('HYPERTC_ANISO')
    gradTperp_mag = sqrt(max(gradT(1)**2 + gradT(2)**2 + gradT(3)**2 &
                           - bgradT**2, 0.0d0))
    if (tc_kappa_perp > 0.0d0) then
      sig_perp = tc_kappa_perp
    else if (tc_kappa0_perp > 0.0d0) then
      sig_perp = tc_kappa0_perp * Te_c**2.5d0
    else
      chi      = htc_Cchi * Bmag * Te_c**1.5d0 / rho_c
      sig_perp = sig_par / (1.0d0 + chi**2)
    end if

    tau_perp = max(4.d0*dt, f_sat*sig_perp*Te_c*(mhd_gamma-1.0d0) / (e_c*cmax_global**2))
    wnew(qperp_) = wnew(qperp_) - qdt*(f_sat*sig_perp*gradTperp_mag + wCTprim1(qperp_,2))/tau_perp
#:endif
#:endif

  else
     !---------------------------------
     ! split sources
     !---------------------------------

     ! Not yet implemented

  end if

end subroutine addsource_compact
#:enddef

#:def to_primitive()
  pure subroutine to_primitive(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(nw_phys)

    u(iw_mom(1))=u(iw_mom(1))/u(iw_rho)
    u(iw_mom(2))=u(iw_mom(2))/u(iw_rho)
    u(iw_mom(3))=u(iw_mom(3))/u(iw_rho)
    u(iw_e)=mhd_gamma_m1*(u(iw_e)-0.5_dp*&
      (u(iw_rho)*(u(iw_mom(1))**2+u(iw_mom(2))**2+u(iw_mom(3))**2)+&
       u(iw_mag(1))**2+u(iw_mag(2))**2+u(iw_mag(3))**2))

  end subroutine to_primitive
#:enddef

#:def to_conservative()  
  pure subroutine to_conservative(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(nw_phys)

    ! Compute energy from pressure and kinetic energy
    u(iw_e)=u(iw_e)/mhd_gamma_m1+0.5_dp*&
      (u(iw_rho)*(u(iw_mom(1))**2+u(iw_mom(2))**2+u(iw_mom(3))**2)+&
       u(iw_mag(1))**2+u(iw_mag(2))**2+u(iw_mag(3))**2)
    ! Compute momentum from density and velocity components
    u(iw_mom(1))=u(iw_rho)*u(iw_mom(1))
    u(iw_mom(2))=u(iw_rho)*u(iw_mom(2))
    u(iw_mom(3))=u(iw_rho)*u(iw_mom(3))

  end subroutine to_conservative
#:enddef

#:def get_flux()
  subroutine get_flux(u, xC, flux_dim, flux)
#:if defined('HYPERTC')
    use mod_global_parameters, only: cmax_global, smalldouble
#:endif
    !$acc routine seq
    real(dp), intent(in)  :: u(nw_phys)
    real(dp), intent(in)  :: xC(1:ndim)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(nw_flux)
    real(dp)              :: ptotal
#:if defined('HYPERTC_ANISO')
    real(dp)              :: Bmag2_tc, b_fd
#:elif defined('HYPERTC')
    real(dp)              :: Bmag_tc
#:endif

    ! Hyperbolic TC field geometry
#:if defined('HYPERTC_ANISO')
    Bmag2_tc = u(iw_mag(1))**2 + u(iw_mag(2))**2 + u(iw_mag(3))**2
    b_fd     = u(iw_mag(flux_dim)) * sqrt(Bmag2_tc) / max(Bmag2_tc, smalldouble**2)
#:elif defined('HYPERTC')
    Bmag_tc = sqrt(max(u(iw_mag(1))**2 + u(iw_mag(2))**2 + u(iw_mag(3))**2, smalldouble**2))
#:endif

#:if defined('MHD_ENERGY_ONLY')
    flux = 0.0_dp
#:if defined('HYPERTC_ANISO')
    flux(iw_e) = u(iw_q)*b_fd + u(iw_qperp)*sqrt(max(1.0_dp - b_fd**2, 0.0_dp))
#:elif defined('HYPERTC')
    flux(iw_e) = u(iw_q) * u(iw_mag(flux_dim)) / Bmag_tc
#:endif
    return
#:endif

    ! Density flux
    flux(iw_rho)=u(iw_rho)*u(iw_mom(flux_dim))

    ! Momentum flux with pressure term
    flux(iw_mom(1))=u(iw_rho)*u(iw_mom(1))*u(iw_mom(flux_dim))-&
      u(iw_mag(flux_dim))*u(iw_mag(1))
    flux(iw_mom(2))=u(iw_rho)*u(iw_mom(2))*u(iw_mom(flux_dim))-&
      u(iw_mag(flux_dim))*u(iw_mag(2))
    flux(iw_mom(3))=u(iw_rho)*u(iw_mom(3))*u(iw_mom(flux_dim))-&
      u(iw_mag(flux_dim))*u(iw_mag(3))
    ptotal=u(iw_e)+0.5_dp*(u(iw_mag(1))**2+u(iw_mag(2))**2+u(iw_mag(3))**2)
    flux(iw_mom(flux_dim))=flux(iw_mom(flux_dim))+ptotal

    ! Energy flux
    flux(iw_e)=u(iw_mom(flux_dim))*(u(iw_e)/mhd_gamma_m1+0.5_dp*&
      u(iw_rho)*(u(iw_mom(1))**2+u(iw_mom(2))**2+u(iw_mom(3))**2)+&
      2.0_dp*ptotal-u(iw_e))-u(iw_mag(flux_dim))*&
      (u(iw_mag(1))*u(iw_mom(1))+u(iw_mag(2))*u(iw_mom(2))+u(iw_mag(3))*u(iw_mom(3)))

    ! Magnetic flux
    flux(iw_mag(1))=u(iw_mom(flux_dim))*u(iw_mag(1))-u(iw_mag(flux_dim))*u(iw_mom(1))
    flux(iw_mag(2))=u(iw_mom(flux_dim))*u(iw_mag(2))-u(iw_mag(flux_dim))*u(iw_mom(2))
    flux(iw_mag(3))=u(iw_mom(flux_dim))*u(iw_mag(3))-u(iw_mag(flux_dim))*u(iw_mom(3))

    ! GLM psi flux
    flux(iw_mag(flux_dim))=u(psi_)
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
    flux(psi_)=cmax_global**2*u(iw_mag(flux_dim))

    ! Tracer flux. Note that tracers stay conservative.
#:if defined('N_TRACER')
  #:for i in range(1, N_TRACER_+1)
      flux(tracer(${i}$)) = u(tracer(${i}$)) * u(iw_mom(flux_dim))
  #:endfor
#:endif

    ! Hyperbolic TC fluxes
#:if defined('HYPERTC_ANISO')
    flux(iw_e)     = flux(iw_e) + u(iw_q)*b_fd + u(iw_qperp)*sqrt(max(1.0_dp - b_fd**2, 0.0_dp))
    flux(iw_q)     = 0.0_dp
    flux(iw_qperp) = 0.0_dp
#:elif defined('HYPERTC')
    flux(iw_e) = flux(iw_e) + u(iw_q) * u(iw_mag(flux_dim)) / Bmag_tc
    flux(iw_q) = 0.0_dp
#:endif

  end subroutine get_flux
#:enddef

#:def get_cmax()
!> Returns maximum local signal speed |v_n| + c_f (fast magnetosonic) from primitive state u in direction flux_dim;
!> used in LLF/TVDLF flux estimation.
pure real(dp) function get_cmax(u, x, flux_dim) result(wC)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim

  real(dp) :: inv_rho, b2, cfast2

  inv_rho=1.0_dp/u(iw_rho)
  wC=mhd_gamma*u(iw_e)*inv_rho
  cfast2=(u(iw_mag(1))**2+u(iw_mag(2))**2+u(iw_mag(3))**2)*inv_rho+wC
  b2=cfast2**2-4.0_dp*wC*u(iw_mag(flux_dim))**2*inv_rho
  if(b2<0.d0) b2=0.d0
  wC=sqrt(0.5_dp*(cfast2+sqrt(b2)))+abs(u(iw_mom(flux_dim)))

end function get_cmax
#:enddef  

#:def get_rho()
  pure real(dp) function get_rho(w, x) result(rho)
    !$acc routine seq
    real(dp), intent(in)  :: w(nw_phys)
    real(dp), intent(in)  :: x(1:ndim)

    rho = w(iw_rho)
  end function get_rho
#:enddef

#:def get_pthermal()
pure real(dp) function get_pthermal(w, x) result(pth)
  !$acc routine seq
  real(dp), intent(in)  :: w(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)

  pth = (phys_gamma-1.0_dp)*(w(iw_e)-0.5_dp*sum(w(iw_mom(:))**2)/w(iw_rho) &
       -0.5_dp*(w(iw_mag(1))**2+w(iw_mag(2))**2+w(iw_mag(3))**2))
end function get_pthermal
#:enddef

#:def get_Rfactor()
pure real(dp) function get_Rfactor() result(Rfactor)
  !$acc routine seq
  Rfactor = 1.0d0
end function get_Rfactor
#:enddef

#:def estimate_speeds_minmax()
!> Davis (1988) min/max wave speed estimates wL = min(v_n - c_f), wR = max(v_n + c_f) (fast magnetosonic) over left/right states;
!> used in HLL flux estimation.
subroutine estimate_speeds_minmax(uL, uR, xC, flux_dim, wL, wR)
  !$acc routine seq
  real(dp), intent(in)  :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: wL, wR

  real(dp) :: inv_rho, cs2, cA2, cAn2, sum2, discriminant, cfL, cfR

  ! Left state
  inv_rho = 1.0_dp / uL(iw_rho)
  cs2 = mhd_gamma * uL(iw_e) * inv_rho
  cA2 = (uL(iw_mag(1))**2 + uL(iw_mag(2))**2 + uL(iw_mag(3))**2) * inv_rho
  cAn2 = uL(iw_mag(flux_dim))**2 * inv_rho
  sum2 = cs2 + cA2
  discriminant = max(sum2**2 - 4.0_dp*cs2*cAn2, 0.0_dp)
  cfL = sqrt(0.5_dp * (sum2 + sqrt(discriminant)))

  ! Right state
  inv_rho = 1.0_dp / uR(iw_rho)
  cs2 = mhd_gamma * uR(iw_e) * inv_rho
  cA2 = (uR(iw_mag(1))**2 + uR(iw_mag(2))**2 + uR(iw_mag(3))**2) * inv_rho
  cAn2 = uR(iw_mag(flux_dim))**2 * inv_rho
  sum2 = cs2 + cA2
  discriminant = max(sum2**2 - 4.0_dp*cs2*cAn2, 0.0_dp)
  cfR = sqrt(0.5_dp * (sum2 + sqrt(discriminant)))
  
  wL = min(uL(iw_mom(flux_dim)) - cfL, uR(iw_mom(flux_dim)) - cfR)
  wR = max(uL(iw_mom(flux_dim)) + cfL, uR(iw_mom(flux_dim)) + cfR)

end subroutine estimate_speeds_minmax
#:enddef


#:endif
