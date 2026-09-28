#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute

#:if PHYS == 'srmhd'

#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter, public              :: nw_phys=2+2*ndim+3
  integer, parameter, public              :: nw_flux=2+2*ndim+1
  
  !> Whether synge eos is used
  logical, public                         :: srmhd_eos = .false.
  !$acc declare copyin(srmhd_eos)

  !> Index of the density (in the w array) as primitive or conserved
  integer, public                         :: rho_
  integer, public                         :: d_
  !$acc declare create(rho_,d_)

  !> Indices of the momentum density
  integer, allocatable, public            :: mom(:)
  !$acc declare create(mom)

  !> Indices of the magnetic field
  integer, allocatable, public            :: mag(:)
  !$acc declare create(mag)


  !> Index of the energy density
  integer, public                         :: tau_
  !$acc declare create(tau_)

  !> Index of the gas pressure should equal tau_
  integer, public                         :: p_
  !$acc declare create(p_)


  !> Index of the auxiliary variable mu (1/(lfac*h))
  integer, public            :: mu_
  !$acc declare create(mu_)



  ! may add this for con2prim speedups
  ! !> Index of the upper limit on mu field
  ! integer, allocatable, public            :: mu_plus_
  ! !$acc declare create(mu_plus_)


  !> Index of the Lorentz factor
  integer, public     :: lfac_
  !$acc declare create(lfac_)

  !> Index of GLM psi
  integer, public :: psi_
  ${GPU_DECLARE_CREATE('psi_')}$



!   !> Number of tracer species
!   integer, public                         :: srmhd_n_tracer = 0
!   !$acc declare copyin(srmhd_n_tracer)

  !> The adiabatic index
  double precision, public                :: srmhd_gamma = 5.d0/3.0d0
  !$acc declare copyin(srmhd_gamma)

  ! !> derived values from adiabatic index 
  ! double precision, public                :: gamma_1,inv_gamma_1,gamma_to_gamma_1
  ! !$acc declare copyin(gamma_1,inv_gamma_1,gamma_to_gamma_1)

  !> Helium abundance over Hydrogen
  double precision, public  :: He_abundance=0.1d0
  !$acc declare copyin(He_abundance)

  !> Whether particles module is added
  logical, public                         :: srmhd_particles = .false.
  !$acc declare copyin(srmhd_particles)

  !> switch for source user
  logical, public                         :: srmhd_source_usr = .false.
  !$acc declare copyin(srmhd_source_usr)



#:enddef


#:def read_params()
    !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srmhd_list/ srmhd_eos,srmhd_gamma,srmhd_n_tracer, &
      He_abundance, srmhd_source_usr, &
      srmhd_small_pressure, srmhd_small_density

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, srmhd_list, end=111)
111    close(unitpar)
    end do

    if (srmhd_small_pressure < 0.0d0) call mpistop(&
       "srmhd_small_pressure should be positive.")
    if (srmhd_small_density < 0.0d0) call mpistop(&
       "srmhd_small_density should be positive.")

#ifdef _OPENACC
    !$acc update device(srmhd_eos, &
    !$acc&     srmhd_gamma, srmhd_n_tracer, &
    !$acc&     He_abundance, srmhd_source_usr)
    !$acc update device(srmhd_small_pressure, srmhd_small_density)
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
    double precision :: mp,kB

    !> here no SI_UNIT used by default, to be implemented
    mp = mp_cgs
    kB = kB_cgs
    unit_velocity=const_c

    ! we assume user sets: unit_numberdensity, unit_length, He_abundance
    ! then together with light speed c, all units fixed
    unit_density=(1.0d0+4.0d0*He_abundance)*mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.0d0+3.0d0*He_abundance)*unit_numberdensity*kB)
    unit_time=unit_length/unit_velocity
    unit_mass=unit_density*unit_length**3

    !$acc update device(unit_density, unit_numberdensity, unit_temperature, unit_pressure, unit_velocity, unit_length, unit_time, unit_mass)
  end subroutine phys_units
#:enddef


#:def phys_init()
    !> Initialize the module
  subroutine phys_init()
    use mod_global_parameters
    integer      :: idir

    call read_params(par_files)
    call phys_units()

    phys_energy  = .true.
    phys_total_energy  = .true.
    phys_gamma = srmhd_gamma

    ! gamma_1=srmhd_gamma-1.0d0
    ! inv_gamma_1=1.0d0/gamma_1
    ! gamma_to_gamma_1=srmhd_gamma/gamma_1
    ! !$acc update device(gamma_1,inv_gamma_1,gamma_to_gamma_1)

    phys_internal_e=.false.
    phys_partial_ionization=.false.
    need_global_cmax=.false.

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

 !$acc update device(physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization,need_global_cmax,phys_req_diagonal)

    use_particles = srmhd_particles

    ! Determine flux variables
    rho_ = var_set_rho()
    d_=rho_
    !$acc update device(rho_,d_)

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)
    !$acc update device(mom)

    ! Set index of energy variable
    tau_ = var_set_energy()
    p_ = tau_
    !$acc update device(tau_,p_)

    ! set b field indices
    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)
    !$acc update device(mag)

!     ! Register tracer fields
! #:if defined('N_TRACER')
!     #:for i in range(1, N_TRACER_+1)
!         tracer(${i}$) = var_set_fluxvar("trc", "trp", ${i}$, need_bc=.false.)
!     #:endfor
!     !$acc update device(tracer)
! #:endif

    ! Set index for auxiliary variables
    ! MUST be after the possible tracers (which have fluxes)
    mu_  = var_set_auxvar('mu','mu')
    lfac_= var_set_auxvar('lfac','lfac')
    psi_ = var_set_auxvar('psi', 'psi')
    !$acc update device(mu_,lfac_)

    ! set number of variables which need update ghostcells
    nwgc=nwflux+nwaux
    !$acc update device(nwgc)

    ! Define custom flux types:
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw_flux))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw_flux])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if
    !$acc update device(flux_type)

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1
    iw_vector(2) = mag(1) - 1
    !$acc update device(nvector, iw_vector)

! use cycle, needs to be dealt with:
!    ! Initialize particles module
!    if (srmhd_particles) then
!       call particles_init()
!       phys_req_diagonal = .true.
!    end if

  end subroutine phys_init
#:enddef

#:def phys_get_dt()
  subroutine phys_get_dt(w, x, dx, dtnew)
  !$acc routine seq
    real(dp), intent(in)   :: w(nw_phys), x(1:ndim), dx(1:ndim)
    real(dp), intent(out)  :: dtnew

    dtnew = huge(1.0d0)
    
#:if defined('SOURCE_USR')
    ! TODO: user-set time step limit, also in other physics modules!!!
#:endif    

  end subroutine phys_get_dt
#:enddef 

#:def to_primitive()
  pure subroutine to_primitive(u)
    !$acc routine seq
    use srmhd_con2prim
    real(dp), intent(inout) :: u(nw_phys)

    real(dp) :: p, d, x, tau, lfac
    real(dp) :: mu
    real(dp) :: vel(1:ndim)
    real(dp) :: s_sqr, r_bar_sqr, q_bar, b_sqr, s_dot_b, eps

    s_sqr = u(iw_mom(1))**2 + u(iw_mom(2))**2 + u(iw_mom(3))**2

    b_sqr = u(iw_mag(1))**2 + u(iw_mag(2))**2 + u(iw_mag(3))**2

    s_dot_b = u(iw_mom(1)) * u(iw_mag(1)) + u(iw_mom(2)) * u(iw_mag(2))& 
             +u(iw_mom(3)) * u(iw_mag(3))

    d=u(iw_rho)

    tau=u(iw_tau)

    call con2prim(mu, u(iw_rho), u(iw_tau), s_sqr, b_sqr, s_dot_b)

    x=1/(1 + mu*b_sqr/d) !(26)

    vel(1)= mu*x *(u(iw_mom(1))/d + mu*s_dot_b*u(iw_mag(1))/d**2)

    vel(2)= mu*x *(u(iw_mom(2))/d + mu*s_dot_b*u(iw_mag(2))/d**2)

    vel(3)= mu*x *(u(iw_mom(3))/d + mu*s_dot_b*u(iw_mag(3))/d**2)

    lfac = 1/sqrt(1-(vel(1)**2 + vel(2)**2 + vel(3)**2))

    r_bar_sqr= s_sqr*x**2/d**2 + mu * x * (1+x) * s_dot_b**2/d**3 !(38)

    q_bar=  tau/d - 0.5_dp*b_sqr/d - 0.5_dp* mu**2 * x**2 * (s_sqr*b_sqr/d**3-s_dot_b**2/d**3) !(39)

    eps=lfac*(q_bar - mu*r_bar_sqr) + sum(vi**2)*lfac**2/(1 + lfac) !(42)

    u(iw_rho) = d/lfac  

    u(iw_tau) = ideal_eos(u(iw_rho), eps)

    u(iw_mom(1)) = vel(1)*lfac

    u(iw_mom(2)) = vel(2)*lfac

    u(iw_mom(3)) = vel(3)*lfac

    u(lfac_) = lfac
    

  end subroutine to_primitive
#:enddef


#:def to_conservative
  pure subroutine to_conservative(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(nw_phys)

    ! real(dp)    ::  tau, d
    real(dp)    ::  v_sqr, b_sqr, v_dot_b

    ! d=u(iw_rho)*u(lfac_)

    v_sqr = (u(iw_mom(1))**2 + u(iw_mom(2))**2 + u(iw_mom(3))**2)/u(lfac_)**2

    b_sqr= u(iw_mag(1))**2 + u(iw_mag(2))**2 + u(iw_mag(3))**2

    v_dot_b = (u(iw_mag(1))*u(iw_mom(1)) + u(iw_mag(2))*u(iw_mom(2)) + u(iw_mag(3))*u(iw_mom(3)))/u(lfac_)

    u(iw_mom(1)) = (u(iw_rho)/u(mu_) + b_sqr/u(lfac_))*u(iw_mom(1)) - v_dot_b*u(iw_mag(1))

    u(iw_mom(2)) = (u(iw_rho)/u(mu_) + b_sqr/u(lfac_))*u(iw_mom(2)) - v_dot_b*u(iw_mag(2))

    u(iw_mom(3)) = (u(iw_rho)/u(mu_) + b_sqr/u(lfac_))*u(iw_mom(3)) - v_dot_b*u(iw_mag(3))

    u(iw_tau) = u(iw_rho) / u(mu_) * u(lfac_) - u(iw_tau) + 0.5_dp*b_sqr * (1 + v_sqr)&
              - 0.5_dp*v_dot_b**2 - u(iw_rho)*u(lfac_) 

    u(iw_rho)= u(iw_rho)*u(lfac_)
  
  end subroutine to_conservative
#:enddef

!input primitive
#:def get_flux()
  subroutine get_flux(u, xC, flux_dim, flux)
    use mod_global_parameters, only: cmax_global
    !$acc routine seq
    real(dp), intent(in)  :: u(nw_phys)
    real(dp), intent(in)  :: xC(1:ndim)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(nw_flux)

    real(dp) :: vel(1:ndim), si


    vel(1) = u(iw_mom(1)) / u(lfac_)
    vel(2) = u(iw_mom(2)) / u(lfac_)
    vel(3) = u(iw_mom(3)) / u(lfac_) 

    v_sqr = vel(1)**2 + vel(2)**2 + vel(3)**2 

    b_sqr= u(iw_mag(1))**2 + u(iw_mag(2))**2 + u(iw_mag(3))**2

    v_dot_b = (u(iw_mag(1))*u(iw_mom(1)) + u(iw_mag(2))*u(iw_mom(2)) + u(iw_mag(3))*u(iw_mom(3)))/u(lfac_)


    si = (u(iw_rho) / u(mu_) * u(lfac_) + b_sqr)*vel(flux_dim) - v_dot_b*u(iw_mag(flux_dim))

    ! density flux
    flux(iw_rho) = u(iw_rho) * u(lfac_) * vel(flux_dim)

    ! momentum flux
    flux(iw_mom(1)) = si*vel(1) - u(iw_mag(1)) * u(iw_mag(flux_dim)) / u(lfac_)**2 - v_dot_b* vel(flux_dim)* u(iw_mag(1))

    flux(iw_mom(2)) = si*vel(2) - u(iw_mag(2)) * u(iw_mag(flux_dim)) / u(lfac_)**2 - v_dot_b* vel(flux_dim)* u(iw_mag(2))

    flux(iw_mom(3)) = si*vel(1) - u(iw_mag(3)) * u(iw_mag(flux_dim)) / u(lfac_)**2 - v_dot_b* vel(flux_dim)* u(iw_mag(3))

    flux(iw_mom(flux_dim)) = flux(iw_mom(flux_dim)) + u(iw_tau) + 0.5_dp*(b_sqr*(1 + v_sqr) - v_dot_b**2)

    ! energy flux
    flux(iw_tau) = si - vel(flux_dim)*u(iw_rho)*u(lfac_)

    ! magnetic flux
    flux(iw_mag(1)) = vel(flux_dim) * u(iw_mag(1)) - vel(1) * u(iw_mag(flux_dim))

    flux(iw_mag(2)) = vel(flux_dim) * u(iw_mag(2)) - vel(2) * u(iw_mag(flux_dim))

    flux(iw_mag(3)) = vel(flux_dim) * u(iw_mag(3)) - vel(3) * u(iw_mag(flux_dim))

    !GLM psi flux
    flux(iw_mag(flux_dim))=u(psi_)
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
    flux(psi_)=cmax_global**2*u(iw_mag(flux_dim))


    end subroutine get_flux

  #:enddef


#:def get_cmax()
!> Returns maximum local signal speed from primitive state u in direction flux_dim;
!> used in LLF/TVDLF flux estimation.
pure function get_cmax(u, x, flux_dim) result(wC)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim

  real(dp) :: wC

  wC=1

end function get_cmax
#:enddef  

