#:if PHYS == 'srhd'

#:if defined('N_TRACER')
#:set N_TRACER_ = N_TRACER
#:else
#:set N_TRACER_ = 0
#:endif

#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter, public              :: nw_phys=2+ndim+2+${N_TRACER_}$
  integer, parameter, public              :: nw_flux=2+ndim+${N_TRACER_}$
  
  !> Whether synge eos is used
  logical, public                         :: srhd_eos = .false.
  !$acc declare copyin(srhd_eos)

  !> Index of the density (in the w array) as primitive or conserved
  integer, public                         :: rho_
  integer, public                         :: d_
  !$acc declare create(rho_,d_)

  !> Indices of the momentum density
  integer, allocatable, public            :: mom(:)
  !$acc declare create(mom)

#:if defined('N_TRACER')
  !> Indices of the tracers
  integer, public                         :: tracer(${N_TRACER_}$)
  !$acc declare create(tracer)
#:endif

  !> Index of the energy density
  integer, public                         :: e_
  !$acc declare create(e_)

  !> Index of the gas pressure should equal e_
  integer, public                         :: p_
  !$acc declare create(p_)

  !> Index of the Lorentz factor
  integer, public     :: lfac_
  !$acc declare create(lfac_)

  !> Index of the inertia
  integer, public     :: xi_
  !$acc declare create(xi_)

  !> Number of tracer species
  integer, public                         :: srhd_n_tracer = 0
  !$acc declare copyin(srhd_n_tracer)

  !> The adiabatic index
  double precision, public                :: srhd_gamma = 5.d0/3.0d0
  !$acc declare copyin(srhd_gamma)

  !> derived values from adiabatic index 
  double precision, public                :: gamma_1,inv_gamma_1,gamma_to_gamma_1
  !$acc declare copyin(gamma_1,inv_gamma_1,gamma_to_gamma_1)

  !> Helium abundance over Hydrogen
  double precision, public  :: He_abundance=0.1d0
  !$acc declare copyin(He_abundance)

  !> Whether particles module is added
  logical, public                         :: srhd_particles = .false.
  !$acc declare copyin(srhd_particles)

  !> switch for source user
  logical, public                         :: srhd_source_usr = .false.
  !$acc declare copyin(srhd_source_usr)

#:enddef

#:def read_params()
    !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srhd_list/ srhd_eos,srhd_gamma,srhd_n_tracer, &
      He_abundance, srhd_source_usr

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, srhd_list, end=111)
111    close(unitpar)
    end do

#ifdef _OPENACC
    !$acc update device(srhd_eos, &
    !$acc&     srhd_gamma, srhd_n_tracer, &
    !$acc&     He_abundance, srhd_source_usr)
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
    phys_gamma = srhd_gamma

    gamma_1=srhd_gamma-1.0d0
    inv_gamma_1=1.0d0/gamma_1
    gamma_to_gamma_1=srhd_gamma/gamma_1
    !$acc update device(gamma_1,inv_gamma_1,gamma_to_gamma_1)

    phys_internal_e=.false.
    phys_partial_ionization=.false.
    need_global_cmax=.false.

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

 !$acc update device(physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization,need_global_cmax,phys_req_diagonal)

    use_particles = srhd_particles

    ! Determine flux variables
    rho_ = var_set_rho()
    d_=rho_
    !$acc update device(rho_,d_)

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)
    !$acc update device(mom)

    ! Set index of energy variable
    e_ = var_set_energy()
    p_ = e_
    !$acc update device(e_,p_)

    ! Register tracer fields
#:if defined('N_TRACER')
    #:for i in range(1, N_TRACER_+1)
        tracer(${i}$) = var_set_fluxvar("trc", "trp", ${i}$, need_bc=.false.)
    #:endfor
    !$acc update device(tracer)
#:endif

    ! Set index for auxiliary variables
    ! MUST be after the possible tracers (which have fluxes)
    xi_  = var_set_auxvar('xi','xi')
    lfac_= var_set_auxvar('lfac','lfac')
    !$acc update device(xi_,lfac_)

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

! use cycle, needs to be dealt with:    
!    ! Initialize particles module
!    if (srhd_particles) then
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

#:def addsource_local()
subroutine addsource_local(qdt, dtfactor, qtC, wCT, wCTprim, qt, wnew, x, dr, &
    qsourcesplit)
  !$acc routine seq
#:if defined('SOURCE_USR')
  use mod_usr, only: addsource_usr
#:endif

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCT(nw_phys), wCTprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim), dr(ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit

  if (.not. qsourcesplit) then 
     !---------------------------------
     ! unsplit sources
     !---------------------------------

#:if defined('SOURCE_USR')
     call addsource_usr(qdt, qt, wCT, wCTprim, wnew, x, .false.)
#:endif

  !!!else
     !---------------------------------
     ! split sources     
     !---------------------------------
        
  end if

end subroutine addsource_local
#:enddef

#:def to_primitive()
  pure subroutine to_primitive(u)
    !$acc routine seq
    use mod_con2prim
    real(dp), intent(inout) :: u(nw_phys)

    real(dp) :: rho,rhoh,pth,E
    real(dp) :: ssqr

    !! begin: call srhd_get_auxiliary(u)
    ssqr=(u(iw_mom(1))**2+u(iw_mom(2))**2+u(iw_mom(3))**2)
    if(srhd_eos)then
       call con2prim_eos(u(lfac_),u(xi_),u(iw_rho),ssqr,u(iw_e))
    else
       call con2prim(u(lfac_),u(xi_),u(iw_rho),ssqr,u(iw_e))
    endif
    !! end: call srhd_get_auxiliary(u)

    rho=u(iw_rho)/u(lfac_)
    rhoh=u(xi_)/u(lfac_)**2
    !! begin: call srhd_get_pressure_eos(rho,rhoh,pth,E)
    if(srhd_eos)then
      E = (rhoh+dsqrt(rhoh**2+(srhd_gamma**2-1.0_dp)*rho**2)) &
         /(srhd_gamma+1.0_dp)
      pth = 0.5d0*gamma_1* (E-rho*(rho/E))
    else
      pth =(rhoh-rho)/gamma_to_gamma_1
    endif 
    !! end: call srhd_get_pressure_eos(rho,rhoh,pth,E)

    u(iw_rho)=rho
    u(iw_mom(1))=u(lfac_)*u(iw_mom(1))/u(xi_)
    u(iw_mom(2))=u(lfac_)*u(iw_mom(2))/u(xi_)
    u(iw_mom(3))=u(lfac_)*u(iw_mom(3))/u(xi_)
    u(iw_e)=pth

  end subroutine to_primitive
#:enddef

#:def to_conservative()  
  pure subroutine to_conservative(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(nw_phys)

    real(dp) :: rho,rhoh,pth
    real(dp) :: E,E_th

    rhoh=(u(iw_mom(1))**2+u(iw_mom(2))**2+u(iw_mom(3))**2)
    u(lfac_)=dsqrt(1.0d0+rhoh)

    rho=u(iw_rho)
    pth=u(iw_e)
    !! begin: call srhd_get_enthalpy_eos(rho,pth,rhoh)
    if(srhd_eos) then
     E_th = pth*inv_gamma_1
     E    = E_th+dsqrt(E_th**2+rho**2)
     ! writing rho/E on purpose, for numerics 
     rhoh = 0.5_dp*((srhd_gamma+1.0_dp)*E &
                   - gamma_1*rho*(rho/E))
    else
     rhoh = rho+gamma_to_gamma_1*pth
    end if
    !! end: call srhd_get_enthalpy_eos(rho,pth,rhoh)

    rhoh=rhoh*u(lfac_)
    u(xi_)=u(lfac_)*rhoh

    u(iw_rho)=u(lfac_)*rho
    u(iw_mom(1))=rhoh*u(iw_mom(1))
    u(iw_mom(2))=rhoh*u(iw_mom(2))
    u(iw_mom(3))=rhoh*u(iw_mom(3))
    u(iw_e)=u(xi_)-pth-u(iw_rho)

  end subroutine to_conservative
#:enddef

#:def get_flux()
  subroutine get_flux(u, xC, flux_dim, flux)
    use mod_global_parameters, only:cmax_global
    !$acc routine seq
    real(dp), intent(in)  :: u(nw_phys)
    real(dp), intent(in)  :: xC(1:ndim)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(nw_flux)
    real(dp)              :: ptotal

    real(dp) :: pth,vel(1:ndim)

    pth=u(iw_e)
    vel(1)=u(iw_mom(1))/u(lfac_)
    vel(2)=u(iw_mom(2))/u(lfac_)
    vel(3)=u(iw_mom(3))/u(lfac_)

    ! Density flux
    flux(iw_rho)=u(iw_rho)*u(lfac_)*vel(flux_dim)

    ! Momentum flux with pressure term
    flux(iw_mom(1))=u(xi_)*u(iw_mom(1))*vel(flux_dim)
    flux(iw_mom(2))=u(xi_)*u(iw_mom(2))*vel(flux_dim)
    flux(iw_mom(3))=u(xi_)*u(iw_mom(3))*vel(flux_dim)
    flux(iw_mom(flux_dim))=flux(iw_mom(flux_dim))+pth

    ! Energy flux
    flux(iw_e)=vel(flux_dim)*(u(xi_)-u(iw_rho))

    ! Tracer flux. Note that tracers stay conservative.
#:if defined('N_TRACER')
  #:for i in range(1, N_TRACER_+1)
      flux(tracer(${i}$)) = u(tracer(${i}$)) * vel(flux_dim)
  #:endfor
#:endif

  end subroutine get_flux
#:enddef

#:def get_cmax()
!> Returns maximum local signal speed from primitive state u in direction flux_dim;
!> used in LLF/TVDLF flux estimation.
pure real(dp) function get_cmax(u, x, flux_dim) result(wC)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim

  real(dp) :: tmp1,tmp2,v2,csound2,vidim,cmax,cmin
  real(dp) :: rho,rhoh,pth,E

    ! input u is in primitive form
    ! auxiliaries are filled here
    rho=u(iw_rho)
    rhoh=u(xi_)/u(lfac_)**2.0d0
    pth=u(iw_e)
    !!tmp1=u(iw_rho)
    !!tmp2=u(xi_)/u(lfac_)**2.0d0
    !!v2=u(iw_e)
    !! begin: call srhd_get_csound2_prim_eos(tmp1,tmp2,v2,csound2)
    if(srhd_eos) then
       E = (rhoh+dsqrt(rhoh**2+(srhd_gamma**2-1.0_dp)&
            *rho**2))/(srhd_gamma+1.0_dp)
       csound2=(pth*((srhd_gamma+1.0_dp)&
                          +gamma_1*(rho/E)**2))&
                      /(2.0d0*rhoh)
    else
       csound2=srhd_gamma*pth/rhoh
    end if
    !! end call srhd_get_csound2_prim_eos(tmp1,tmp2,v2,csound2)

    v2=1.0d0-1.0d0/u(lfac_)**2
    vidim = u(iw_mom(flux_dim))/u(lfac_)
    tmp2=vidim**2.0d0
    tmp1=1.0d0-v2*csound2-tmp2*(1.0d0-csound2)
    tmp2=dsqrt(csound2*(1.0_dp-v2)*tmp1)
    tmp1=vidim*(1.0_dp-csound2)
    cmax=(tmp1+tmp2)/(1.0_dp-v2*csound2)
    cmin=(tmp1-tmp2)/(1.0_dp-v2*csound2)
    ! Limit by speed of light
    cmin = max(cmin, - 1.0d0)
    cmin = min(cmin,   1.0d0)
    cmax = max(cmax, - 1.0d0)
    cmax = min(cmax,   1.0d0)
    ! now take extremal value only for dt limit
    wC = max(dabs(cmax),dabs(cmin))

end function get_cmax
#:enddef  

#:def get_rho()
  pure real(dp) function get_rho(w, x) result(rho)
    !$acc routine seq
    real(dp), intent(in)  :: w(nw_phys)
    real(dp), intent(in)  :: x(1:ndim)

    ! TODO: only for local rad loss
  end function get_rho
#:enddef

#:def get_pthermal()
pure real(dp) function get_pthermal(w, x) result(pth)
  !$acc routine seq
  real(dp), intent(in)  :: w(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)

    ! TODO: only for local rad loss
end function get_pthermal
#:enddef

#:def get_Rfactor()
pure real(dp) function get_Rfactor() result(Rfactor)
  !$acc routine seq

    ! TODO: only for local rad loss
end function get_Rfactor
#:enddef

#:def estimate_speeds_minmax()
!> used in HLL flux estimation.
subroutine estimate_speeds_minmax(uL, uR, xC, flux_dim, wL, wR)
  !$acc routine seq
  real(dp), intent(in)  :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: wL, wR

  ! TODO

end subroutine estimate_speeds_minmax
#:enddef


#:endif
