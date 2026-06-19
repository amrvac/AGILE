#:if PHYS == 'ffhd'
  
#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  !> Base number of physical variables (density, momentum, energy, (q), b1, b2, b3)
  integer, parameter, public              :: nw_phys_base=6
  !> Total number of physical variables (base + conditional additions)
  !> This will be set to nw_phys_base + 1 if HYPERTC is defined
  integer, parameter, public              :: nw_phys=nw_phys_base &
#:if defined('HYPERTC')
    + 1
#:else
    + 0
#:endif

  integer, parameter, public              :: nw_flux_base=3
  integer, parameter, public              :: nw_flux=nw_flux_base &
#:if defined('HYPERTC')
    + 1
#:else
    + 0
#:endif

  !> Whether an energy equation is used
  logical, public                         :: ffhd_energy = .true.
  !$acc declare copyin(ffhd_energy)

  !> Index of the density (in the w array)
  integer, public                         :: rho_
  !$acc declare create(rho_)

  !> Indices of the momentum density
  integer, allocatable, public            :: mom(:)
  !$acc declare create(mom)

  !> Index of the energy density (-1 if not present)
  integer, public                         :: e_
  !$acc declare create(e_)

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public                         :: p_
  !$acc declare create(p_)

  !> Index of the hyperbolic flux variable
  integer, public                         :: q_
  !$acc declare create(q_)

  ! !> Index of the frozen magnetic field
  integer, public                         :: iw_b1
  integer, public                         :: iw_b2
  integer, public                         :: iw_b3
  !$acc declare create(iw_b1, iw_b2, iw_b3)

  !> The adiabatic index
  double precision, public                :: ffhd_gamma = 5.d0/3.0d0
  !$acc declare copyin(ffhd_gamma)

  !> The adiabatic constant
  double precision, public                :: ffhd_adiab = 1.0d0
  !$acc declare copyin(ffhd_adiab)

  !> The helium abundance
  double precision, public                :: He_abundance=0.1d0
  !$acc declare copyin(He_abundance)

  !> Whether p*divb source term is not zero
  logical, public                         :: ffhd_pdivb = .false.
  !$acc declare copyin(ffhd_pdivb)

  !> switch for hyperbolic thermal conduction
  logical, public                         :: ffhd_hyperbolic_thermal_conduction = .false.
  !$acc declare copyin(ffhd_hyperbolic_thermal_conduction)

  !> Whether plasma is partially ionized
  logical, public                         :: ffhd_partial_ionization = .false.
  !$acc declare copyin(ffhd_partial_ionization)

  !> switch for gravity
  logical, public                         :: ffhd_gravity = .false.
  !$acc declare copyin(ffhd_gravity)

  !> switch for radiative cooling
  logical, public                         :: ffhd_radiative_cooling = .false.
  !$acc declare copyin(ffhd_radiative_cooling)

  !> switch for source user
  logical, public                         :: ffhd_source_usr = .false.
  !$acc declare copyin(ffhd_source_usr)

#:enddef

#:def read_params()
    !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /ffhd_list/ ffhd_energy, ffhd_gamma, ffhd_partial_ionization, ffhd_gravity, &
          ffhd_radiative_cooling, ffhd_hyperbolic_thermal_conduction, ffhd_source_usr, ffhd_pdivb, He_abundance

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, ffhd_list, end=111)
111    close(unitpar)
    end do

#ifdef _OPENACC
    !$acc update device( &
    !$acc&         ffhd_energy, ffhd_gamma, ffhd_partial_ionization, &
    !$acc&         ffhd_gravity, ffhd_radiative_cooling, ffhd_hyperbolic_thermal_conduction, &
    !$acc&         ffhd_source_usr, ffhd_pdivb, He_abundance )
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
  end subroutine phys_units
#:enddef
  
#:def phys_init()
    !> Initialize the module
  subroutine phys_init()

    use mod_global_parameters
    #:if defined('COOLING')
    use mod_radiative_cooling, only: rc_fl, radiative_cooling_init_params, radiative_cooling_init
    #:endif

    call read_params(par_files)
    call phys_units()

    phys_energy  = ffhd_energy
    phys_total_energy  = ffhd_energy
    phys_internal_e = .false.
    phys_gamma = ffhd_gamma
    phys_partial_ionization=ffhd_partial_ionization
 !$acc update device(physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization)

    ! Determine flux variables
    rho_ = var_set_rho()
    !$acc update device(rho_)

    allocate(mom(1))
    mom(:) = var_set_momentum(1)
    !$acc update device(mom)

    ! Set index of energy variable
    if (ffhd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if
    !$acc update device(e_,p_)

    ! Set index for heat flux
#:if defined('HYPERTC')
    q_ = var_set_q()
    need_global_cmax = .true.
    !$acc update device(q_)
    !$acc update device(need_global_cmax)
#:endif

    ! Set index of frozen magnetic field
    iw_b1 = var_set_extravar('b1', 'b1')
    iw_b2 = var_set_extravar('b2', 'b2')
    iw_b3 = var_set_extravar('b3', 'b3')
    !$acc update device(iw_b1, iw_b2, iw_b3)

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
    !$acc update device(flux_type)
    
#:if defined('COOLING')
    call radiative_cooling_init_params(phys_gamma,He_abundance)
    call radiative_cooling_init(rc_fl)
    !$acc update device(rc_fl)
    !$acc enter data copyin(rc_fl%tcool,rc_fl%Lcool, rc_fl%Yc)
#:endif

    ! TODO: compile-time nw_phys/nw_flux and runtime nwflux/nw should
    ! be unified so this can't go out of sync. Until then, assert here
    if (nwflux /= nw_flux) then
       write(*,'(A,I0,A,I0)') "ASSERT FAIL: nwflux=", nwflux, " nw_flux=", nw_flux
       call mpistop("nwflux /= nw_flux: var_set_* calls do not match nw_flux parameter")
    end if
    if (nw /= nw_phys) then
       write(*,'(A,I0,A,I0)') "ASSERT FAIL: nw=", nw, " nw_phys=", nw_phys
       call mpistop("nw /= nw_phys: var_set_* calls do not match nw_phys parameter")
    end if

  end subroutine phys_init
#:enddef

#:def phys_get_dt()
  subroutine phys_get_dt(w, x, dx, dtnew)
  !$acc routine seq
#:if defined('GRAVITY')
  use mod_usr, only: gravity_field
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
  use mod_radiative_cooling, only: radiative_cooling_add_source
#:endif

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCT(nw_phys), wCTprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim), dr(ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit
  ! .. local ..
  integer                  :: idim
  real(dp)                 :: field, mag, divb

  if (.not. qsourcesplit) then 
     !---------------------------------
     ! unsplit sources
     !---------------------------------
     
#:if defined('GRAVITY')
     do idim = 1, ndim
        field = gravity_field(wCT, x, idim)
        mag = wCTprim(iw_b1-1+idim)
        wnew(iw_mom(1)) = wnew(iw_mom(1)) + qdt * field * wCT(iw_rho) * mag
        wnew(iw_e)      = wnew(iw_e) + qdt * field * wCT(iw_mom(1)) * mag
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

     ! Not yet implemented

  end if
  
end subroutine addsource_local
#:enddef

#:def addsource_nonlocal()
subroutine addsource_nonlocal(qdt, dtfactor, qtC, wCTprim, qt, wnew, x, dx, idir, &
     qsourcesplit)
  !$acc routine seq
  use mod_global_parameters, only: dt, cmax_global, courantpar, third, &
       unit_temperature, unit_length, unit_density, unit_velocity

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCTprim(nw_phys,5)
  real(dp), intent(in)     :: x(1:ndim), dx(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  integer, intent(in)      :: idir
  logical, intent(in)      :: qsourcesplit
  ! .. local ..
  real(dp)                 :: tau, htc_qrsc, sig_par
  real(dp)                 :: Te(1:5), gradT
  real(dp)                 :: mag5(1:5), divb, mag

  if (.not. qsourcesplit) then 
     !---------------------------------
     ! unsplit sources
     !---------------------------------

#:if defined('PDIVB')
     ! > p*divb 
     mag5(1:5) = wCTprim(iw_b1-1+idir,1:5)
     divb = (8*mag5(4)-8*mag5(2)-mag5(5)+mag5(1))/12.0_dp/dx(idir)
     wnew(iw_mom(1)) = wnew(iw_mom(1)) + qdt*wCTprim(iw_e,3)*divb
#:endif

#:if defined('HYPERTC')
     !> gradient of temperature:
     Te(1:5) = wCTprim(iw_e,1:5) / wCTprim(iw_rho,1:5)
     mag = wCTprim(iw_b1-1+idir,3)

     gradT = (8.d0*(Te(4)-Te(2))-Te(5)+Te(1))/(12.d0*dx(idir))

     sig_par = 8.d-7_dp * unit_temperature**3.5_dp &
             / (unit_length * unit_density * unit_velocity**3.0_dp) &
             * sqrt(Te(3)**5)
     tau = max(4.d0*dt, sig_par*Te(3)*courantpar**2*(phys_gamma-1.0d0)/&
        (wCTprim(iw_e,3)*cmax_global**2))

     htc_qrsc = sig_par * mag * gradT

     wnew(iw_q) = wnew(iw_q) - qdt * (htc_qrsc + wCTprim(iw_q,3)*third) / tau
#:endif

  else
     !---------------------------------
     ! split sources     
     !---------------------------------

     ! Not yet implemented

  end if
  
end subroutine addsource_nonlocal
#:enddef

#:def to_primitive()
pure subroutine to_primitive(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(nw_phys)

  u(iw_mom(1)) = u(iw_mom(1))/u(iw_rho)

  u(iw_e) = (phys_gamma-1.0_dp) * (u(iw_e) - 0.5_dp * u(iw_rho) * &
     u(iw_mom(1))**2 )

end subroutine to_primitive
#:enddef

#:def to_conservative()  
pure subroutine to_conservative(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(nw_phys)
  real(dp)                :: inv_gamma_m1

  inv_gamma_m1 = 1.0_dp/(phys_gamma - 1.0_dp)

  ! Compute energy from pressure and kinetic energy
  u(iw_e) = u(iw_e) * inv_gamma_m1 + 0.5_dp * u(iw_rho) * &
     u(iw_mom(1))**2

  ! Compute momentum from density and velocity components
  u(iw_mom(1)) = u(iw_rho) * u(iw_mom(1))
  
end subroutine to_conservative
#:enddef

#:def get_flux()
subroutine get_flux(u, xC, flux_dim, flux)
  !$acc routine seq

  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: xC(1:ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: flux(nw_flux)
  real(dp)              :: inv_gamma_m1
  real(dp)              :: mag

  inv_gamma_m1 = 1.0_dp/(phys_gamma - 1.0_dp)

  mag = u(iw_b1-1+flux_dim)

  ! Density flux
  flux(iw_rho) = u(iw_rho) * u(iw_mom(1)) * mag

  ! Momentum flux with pressure term
  flux(iw_mom(1)) = (u(iw_rho)*u(iw_mom(1))**2 + u(iw_e)) * mag
  
  ! Energy flux with hyperbolic conduction included
  flux(iw_e) = u(iw_mom(1))*(u(iw_e)*inv_gamma_m1 + &
               0.5_dp*u(iw_rho)*u(iw_mom(1))**2 + u(iw_e)) * mag + &
#:if defined('HYPERTC')
  u(iw_q)*mag

  flux(iw_q) = 0.0_dp
#:else
  0.0_dp
#:endif

end subroutine get_flux
#:enddef

#:def get_cmax()  
pure real(dp) function get_cmax(u, x, flux_dim) result(wC)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim
  real(dp)              :: mag

  mag = u(iw_b1-1+flux_dim)
  
  wC = dsqrt(phys_gamma*u(iw_e)/u(iw_rho)) + abs(u(iw_mom(1))*mag)

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

  pth = (phys_gamma-1.0_dp)*(w(iw_e)-0.5_dp*w(iw_mom(1))**2/w(iw_rho))
end function get_pthermal
#:enddef

#:def get_Rfactor()
pure real(dp) function get_Rfactor() result(Rfactor)
  !$acc routine seq
  Rfactor = 1.0d0
end function get_Rfactor
#:enddef

#:endif
