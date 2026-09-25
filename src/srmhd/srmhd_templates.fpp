#:if PHYS == 'srmhd'


#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter, public              :: nw_phys=2+2*ndim+2
  integer, parameter, public              :: nw_flux=2+2*ndim
  
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
  integer, public                         :: e_
  !$acc declare create(e_)

  !> Index of the gas pressure should equal e_
  integer, public                         :: p_
  !$acc declare create(p_)

!   !> Index of the Lorentz factor
!   integer, public     :: lfac_
!   !$acc declare create(lfac_)

!   !> Index of the inertia
!   integer, public     :: xi_
!   !$acc declare create(xi_)

!   !> Number of tracer species
!   integer, public                         :: srmhd_n_tracer = 0
!   !$acc declare copyin(srmhd_n_tracer)

  !> The adiabatic index
  double precision, public                :: srmhd_gamma = 5.d0/3.0d0
  !$acc declare copyin(srmhd_gamma)

  !> derived values from adiabatic index 
  double precision, public                :: gamma_1,inv_gamma_1,gamma_to_gamma_1
  !$acc declare copyin(gamma_1,inv_gamma_1,gamma_to_gamma_1)

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


#:def to_primitive()
  pure subroutine to_primitive(u)
    !$acc routine seq
    use srmhd_con2prim
    real(dp), intent(inout) :: u(nw_phys)

    real(dp) :: rho, p, d
    real(dp) :: mu
    real(dp) :: vel(3)
    real(dp) :: v_sqr, r_bar_sqr, q_bar

    s_sqr = u(iw_mom(1))**2 + u(iw_mom(2))**2 + u(iw_mom(3))**2

    b_sqr = u(iw_mag(1))**2 + u(iw_mag(2))**2 + u(iw_mag(3))**2

    s_dot_b = u(iw_mom(1))*u(iw_mag(1)) + u(iw_mom(2))*u(iw_mag(2))& 
             +u(iw_mom(3))*u(iw_mag(3))

    d=u(iw_rho)

    tau=u(iw_e)

    call con2prim(mu, u(iw_rho), u(iw_e), s_sqr, b_sqr, s_dot_b)



    x=1/(1 + mu*b_sqr/d) !(26)

    vel(1)= mu*x *(u(iw_mom(1))/d + mu*s_dot_b*u(iw_mag(1))/d**2)

    vel(2)= mu*x *(u(iw_mom(2))/d + mu*s_dot_b*u(iw_mag(2))/d**2)

    vel(3)= mu*x *(u(iw_mom(3))/d + mu*s_dot_b*u(iw_mag(3))/d**2)






  end subroutine to_primitive
#:enddef







