!> Module with finite volume methods for fluxes

#:mute
#:include "physics/mod_physics_templates.fpp"
#:endmute

module mod_finite_volume

  use mod_variables
  use mod_physics_vars
  use mod_global_parameters, only: ndim
  use mod_physicaldata
  implicit none

  private

  public :: finite_volume_local

contains

! instantiate the templated functions here for inlining:
@:addsource_local()
@:addsource_nonlocal()
@:to_primitive()
@:to_conservative()
@:get_cmax()
@:get_flux()
@:estimate_speeds_minmax()
@:estimate_speeds_toro_pvrs()

  ! flux scheme list : (scheme_tag, method_enum, faceflux_proc)
#:set schemes = [ &
  ('muscl_llf',  'fs_tvdlf', 'reconflux_muscl_llf_prim'), &
  ('muscl_hll',  'fs_hll',   'reconflux_muscl_hll_prim'), &
  ('muscl_hllc', 'fs_hllc',  'reconflux_muscl_hllc_prim')]

  subroutine finite_volume_local(qdt, dtfactor, ixImin1,ixImin2,&
    ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
    ixOmax3, idimsmin,idimsmax, qtC, bga, qt, bgb, fC, fE)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    double precision, intent(in)                                       :: qdt,&
      dtfactor, qtC, qt
    integer, intent(in)                                                :: &
      ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
      ixOmin3,ixOmax1,ixOmax2,ixOmax3, idimsmin,idimsmax
    ! remember, old names map as: wCT => bga, wnew => bgb
    type(block_grid_t)                                    :: bga
    type(block_grid_t)                                    :: bgb
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
      ixImin3:ixImax3, 1:nwflux, 1:ndim)  :: fC !not yet provided
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
      ixImin3:ixImax3, sdim:3)            :: fE !not yet provided
    ! .. local ..
    integer                :: n, iigrid, ix1,ix2,ix3
    double precision       :: uprim(nw_phys, ixImin1:ixImax1,ixImin2:ixImax2,&
      ixImin3:ixImax3)
    real(dp)               :: tmp(nw_phys,5)
    real(dp)               :: f(nw_flux, 2)
    real(dp)               :: inv_dr(ndim)
    real(dp)               :: dr(ndim)
    integer                :: typelim
    real(dp)               :: xloc(ndim)
    real(dp)               :: xlocC(ndim,2)
    real(dp)               :: wprim(nw_phys), wCT(nw_phys), wnew(nw_phys)
    integer :: method
    !-----------------------------------------------------------------------------

    method = flux_method(1)  ! TODO: implement per grid level
    select case (method)
#:for scheme_tag, method_enum, faceflux_proc in schemes
    case (${method_enum}$)
      call finite_volume_local_${scheme_tag}$(qdt, dtfactor, ixImin1,ixImin2,&
            ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3, idimsmin,idimsmax, qtC, bga, qt, bgb, fC, fE)
#:endfor
    case default
      call mpistop("finite_volume_local: unknown flux scheme")
    end select
end subroutine finite_volume_local

#:def FV_KERNEL(scheme_tag, faceflux_proc)
  subroutine finite_volume_local_${scheme_tag}$(qdt, dtfactor, ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, idimsmin,idimsmax, qtC, bga, qt, bgb, fC, fE)
    use mod_global_parameters

    double precision, intent(in)                                       :: qdt,&
        dtfactor, qtC, qt
    integer, intent(in)                                                :: &
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3, idimsmin,idimsmax
    ! remember, old names map as: wCT => bga, wnew => bgb
    type(block_grid_t)                                    :: bga
    type(block_grid_t)                                    :: bgb
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nwflux, 1:ndim)  :: fC !not yet provided
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, sdim:3)            :: fE !not yet provided
    ! .. local ..
    integer                :: n, iigrid, ix1,ix2,ix3
    double precision       :: uprim(nw_phys, ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    real(dp)               :: tmp(nw_phys,5)
    real(dp)               :: f(nw_flux, 2)
    real(dp)               :: inv_dr(ndim)
    real(dp)               :: dr(ndim)
    integer                :: typelim
    real(dp)               :: xloc(ndim)
    real(dp)               :: xlocC(ndim,2)
    real(dp)               :: wprim(nw_phys), wCT(nw_phys), wnew(nw_phys)
    integer, parameter     :: max_batch=4096
    integer                :: nbatches, ibatch, igrid_beg, igrid_end
    !-----------------------------------------------------------------------------

    ! batch-launch the kernels (oom issue when many ~30000 blocks):
    nbatches = (igridstail_active + max_batch - 1) / max_batch ! ceiling division
    
    do ibatch = 1, nbatches
       igrid_beg = (ibatch-1) * max_batch + 1
       igrid_end = min(ibatch * max_batch, igridstail_active)

       !$acc parallel loop gang private(uprim, inv_dr, dr, n) default(present)
       do iigrid = igrid_beg, igrid_end
          n = igrids_active(iigrid)

          dr  = rnode(rpdx1_:rnodehi, n)
          inv_dr  = 1/dr
          typelim = type_limiter(node(plevel_, n))

          !$acc loop collapse(ndim) vector
          do ix3=ixImin3,ixImax3 
             do ix2=ixImin2,ixImax2 
                do ix1=ixImin1,ixImax1 
                   ! Convert to primitive
                   uprim(1:nw_phys, ix1,ix2,ix3) = bga%w(ix1,ix2,ix3, 1:nw_phys, n)
                   call to_primitive(uprim(1:nw_phys, ix1,ix2,ix3))
                end do
             end do
          end do
       !!end do
       !!!$acc end parallel loop

       !!!$acc parallel loop gang private(uprim, inv_dr, dr, n) default(present)
       !!do iigrid = igrid_beg, igrid_end
       !!   n = igrids_active(iigrid)

       !!   dr  = rnode(rpdx1_:rnodehi, n)
       !!   inv_dr  = 1/dr
       !!   typelim = type_limiter(node(plevel_, n))

       !!!$acc loop vector collapse(ndim) private(f, wnew, tmp, xlocC, xloc#{if defined('SOURCE_LOCAL')}#, wCT, wprim #{endif}#)
       do ix3=ixOmin3,ixOmax3 
          do ix2=ixOmin2,ixOmax2 
             do ix1=ixOmin1,ixOmax1 
                ! Compute fluxes in all dimensions

                tmp = uprim(1:nw_phys, ix1-2:ix1+2, ix2, ix3)
                xlocC(1:ndim,1) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(1:ndim,2) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(1,1) = xlocC(1,1)-0.5_dp*dr(1)
                xlocC(1,2) = xlocC(1,2)+0.5_dp*dr(1)
                call ${faceflux_proc}$(tmp, xlocC, 1, f, typelim)
                bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = bgb%w(ix1, ix2, ix3, 1:nw_flux,&
                     n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(1)

                tmp = uprim(1:nw_phys, ix1, ix2-2:ix2+2, ix3)
                xlocC(1:ndim,1) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(1:ndim,2) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(2,1) = xlocC(2,1)-0.5_dp*dr(2)
                xlocC(2,2) = xlocC(2,2)+0.5_dp*dr(2)
                call ${faceflux_proc}$(tmp, xlocC, 2, f, typelim)
                bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = bgb%w(ix1, ix2, ix3, 1:nw_flux,&
                     n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(2)

                tmp = uprim(1:nw_phys, ix1, ix2, ix3-2:ix3+2)
                xlocC(1:ndim,1) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(1:ndim,2) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                xlocC(3,1) = xlocC(3,1)-0.5_dp*dr(3)
                xlocC(3,2) = xlocC(3,2)+0.5_dp*dr(3)
                call ${faceflux_proc}$(tmp, xlocC, 3, f, typelim)
                bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = bgb%w(ix1, ix2, ix3, 1:nw_flux,&
                     n) + qdt * (f(:, 1) - f(:, 2)) * inv_dr(3)

#:if defined('SOURCE_LOCAL')
                   ! Add local source terms:
                   xloc(1:ndim) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                   wprim        = uprim(1:nw_phys, ix1, ix2, ix3)
                   wCT          = bga%w(ix1, ix2, ix3, 1:nw_phys, n)
                   wnew         = bgb%w(ix1, ix2, ix3, 1:nw_phys, n)
                   call addsource_local(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                        dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, wCT,&
                        wprim, qt, wnew, xloc, dr, .false. )
                   bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = wnew(1:nw_flux)
#:endif             

#:if defined('SOURCE_NONLOCAL')
                   ! Add non-local (gradient) source terms:
                   xloc(1:ndim) = ps(n)%x(ix1, ix2, ix3, 1:ndim)
                   wnew         = bgb%w(ix1, ix2, ix3, 1:nw_phys, n)

                   tmp = uprim(1:nw_phys, ix1-2:ix1+2, ix2, ix3)
                   call addsource_nonlocal(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                        dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, tmp,&
                        qt, wnew, xloc, dr, 1, .false. )

                   tmp = uprim(1:nw_phys, ix1, ix2-2:ix2+2, ix3)
                   call addsource_nonlocal(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                        dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, tmp,&
                        qt, wnew, xloc, dr, 2, .false. )

                   tmp = uprim(1:nw_phys, ix1, ix2, ix3-2:ix3+2)
                   call addsource_nonlocal(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
                        dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), qtC, tmp,&
                        qt, wnew, xloc, dr, 3, .false. )

                   bgb%w(ix1, ix2, ix3, 1:nw_flux, n) = wnew(1:nw_flux)           
#:endif                
                end do
             end do
          end do
       end do
       !$acc end parallel loop

    end do

  end subroutine finite_volume_local_${scheme_tag}$
  #:enddef


#:for scheme_tag, method_enum, faceflux_proc in schemes
#:call FV_KERNEL(scheme_tag, faceflux_proc)
#:endcall
#:endfor


  !> MUSCL reconstruction in primitive variables for two faces using a 5-point stencil.
  !> Returns uL(:,iface), uR(:,iface) for iface=1 (between cells 2-3) and iface=2 (between 3-4).
  pure subroutine muscl_reconstruct_prim(u, typelim, uL, uR)
    !$acc routine seq
    use mod_limiter, only: limiter_minmod, limiter_vanleer
    real(dp), intent(in)  :: u(nw_phys,5)
    integer,  intent(in)  :: typelim
    real(dp), intent(out) :: uL(nw_phys,2), uR(nw_phys,2)

    real(dp) :: sig(nw_phys,3)   ! slopes at cells 2,3,4
    integer  :: iw

    select case (typelim)
    case (limiter_minmod)
      do iw=1,nw_phys
        sig(iw,1) = minmod (u(iw,2)-u(iw,1), u(iw,3)-u(iw,2))  ! cell 2
        sig(iw,2) = minmod (u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))  ! cell 3
        sig(iw,3) = minmod (u(iw,4)-u(iw,3), u(iw,5)-u(iw,4))  ! cell 4
      end do
    case (limiter_vanleer)
      do iw=1,nw_phys
        sig(iw,1) = vanleer(u(iw,2)-u(iw,1), u(iw,3)-u(iw,2))
        sig(iw,2) = vanleer(u(iw,3)-u(iw,2), u(iw,4)-u(iw,3))
        sig(iw,3) = vanleer(u(iw,4)-u(iw,3), u(iw,5)-u(iw,4))
      end do
    case default  ! Fallback: Godunov (piecewise constant)
      do iw=1,nw_phys
        sig(iw,1)=0._dp
        sig(iw,2)=0._dp
        sig(iw,3)=0._dp
      end do
    end select

    ! Face 1: between cells 2 and 3
    uL(:,1) = u(:,2) + 0.5_dp*sig(:,1)
    uR(:,1) = u(:,3) - 0.5_dp*sig(:,2)

    ! Face 2: between cells 3 and 4
    uL(:,2) = u(:,3) + 0.5_dp*sig(:,2)
    uR(:,2) = u(:,4) - 0.5_dp*sig(:,3)
  end subroutine muscl_reconstruct_prim


  !> One-face LLF/Rusanov numerical flux from primitive L/R states.
  subroutine riemann_llf_prim(uL, uR, xC, flux_dim, F)
    !$acc routine seq
    real(dp), intent(inout) :: uL(nw_phys), uR(nw_phys)
    real(dp), intent(in)    :: xC(ndim)
    integer,  intent(in)    :: flux_dim
    real(dp), intent(out)   :: F(nw_flux)

    real(dp) :: flux_l(nw_flux), flux_r(nw_flux)
    real(dp) :: wL, wR, wmax

    call get_flux(uL, xC, flux_dim, flux_l)
    call get_flux(uR, xC, flux_dim, flux_r)

    wL   = get_cmax(uL, xC, flux_dim)
    wR   = get_cmax(uR, xC, flux_dim)
    wmax = max(abs(wL), abs(wR))

    call to_conservative(uL)
    call to_conservative(uR)

    F = 0.5_dp * ((flux_l + flux_r) - wmax * (uR(1:nw_flux) - uL(1:nw_flux)))
  end subroutine riemann_llf_prim


  !> One-face HLL numerical flux from primitive L/R states.
  subroutine riemann_hll_prim(uL, uR, xC, flux_dim, F)
    !$acc routine seq
    real(dp), intent(inout) :: uL(nw_phys), uR(nw_phys)
    real(dp), intent(in)    :: xC(ndim)
    integer,  intent(in)    :: flux_dim
    real(dp), intent(out)   :: F(nw_flux)

    real(dp) :: flux_l(nw_flux), flux_r(nw_flux)
    real(dp) :: sL, sR, ds, wmax
    real(dp), parameter :: eps = 1e-14_dp

    call get_flux(uL, xC, flux_dim, flux_l)
    call get_flux(uR, xC, flux_dim, flux_r)

    call estimate_speeds_minmax(uL, uR, xC, flux_dim, sL, sR)
    wmax = max(abs(sL), abs(sR))  ! for fallback

    call to_conservative(uL)
    call to_conservative(uR)

    if (sL .ge. 0._dp) then
      F = flux_l
    else if (sR .le. 0._dp) then
      F = flux_r
    else
      ds = sR - sL
      if (ds > eps * (abs(sL) + abs(sR))) then
        F = (sR*flux_l - sL*flux_r + sL*sR*(uR(1:nw_flux) - uL(1:nw_flux))) / ds
      else  ! Fall back to LLF
        F = 0.5_dp * ((flux_l + flux_r) - wmax * (uR(1:nw_flux) - uL(1:nw_flux)))
      end if
    end if
  end subroutine riemann_hll_prim


  !> One-face HLLC numerical flux from primitive L/R states (HD only).
  !> Uses scalar pressure assumption; not valid for MHD (use HLLD instead).
  !> Reference: Toro (2010), chapter 10 (Variant 2)
  subroutine riemann_hllc_prim(uL, uR, xC, flux_dim, F)
    !$acc routine seq
    real(dp), intent(inout) :: uL(nw_phys), uR(nw_phys)
    real(dp), intent(in)    :: xC(ndim)
    integer,  intent(in)    :: flux_dim
    real(dp), intent(out)   :: F(nw_flux)

    real(dp) :: flux_l(nw_flux), flux_r(nw_flux)
    real(dp) :: sL, sR, s_star
    real(dp) :: dL, dR, denom
    real(dp) :: rhoL, rhoR, unL, unR, pL, pR
    real(dp) :: PLR
    real(dp) :: Dstar
    real(dp) :: wmax
    real(dp), parameter :: eps = 1e-14_dp
    integer :: i

    call get_flux(uL, xC, flux_dim, flux_l)
    call get_flux(uR, xC, flux_dim, flux_r)

    call estimate_speeds_toro_pvrs(uL, uR, xC, flux_dim, sL, sR)

    ! cache primitive scalars before uL/uR are overwritten to conservative
    rhoL = uL(iw_rho)
    rhoR = uR(iw_rho)
    unL  = uL(iw_mom(flux_dim))
    unR  = uR(iw_mom(flux_dim))
    pL   = uL(iw_e)
    pR   = uR(iw_e)

    ! we re-use these later (Eq. 10.70, Eq. 10.76)
    dL    = rhoL*(sL - unL)
    dR    = rhoR*(sR - unR)
    denom = dL - dR

    ! needed because (Eq. 10.75) uses U_K
    call to_conservative(uL)
    call to_conservative(uR)

    ! fall back to LLF/Rusanov if denom is ~0
    if (abs(denom) <= eps*(abs(dL) + abs(dR) + 1._dp)) then
      wmax = max(abs(sL), abs(sR))
      F = 0.5_dp * ((flux_l + flux_r) - wmax * (uR(1:nw_flux) - uL(1:nw_flux)))
      return
    end if

    ! contact speed S* (Eq. 10.70)
    s_star = ((pR - pL) + unL*dL - unR*dR) / denom

    ! pressure-like term (Eq. 10.76)
    PLR = 0.5_dp * (pL + pR + dL*(s_star - unL) + dR*(s_star - unR))

    ! HLLC flux (Eq. 10.71) + variant 2 (Eq. 10.75)
    if (0._dp .le. sL) then
      F = flux_l
      return
    elseif (0._dp .ge. sR) then
      F = flux_r
      return
    end if

    if (0._dp .le. s_star) then
      ! F_*L
      do i = 1, nw_flux
        Dstar = 0._dp
        if (i .eq. iw_mom(flux_dim)) Dstar = 1._dp
        if (i .eq. iw_e)             Dstar = s_star
        F(i) = (s_star*(sL*uL(i) - flux_l(i)) + sL*PLR*Dstar) / (sL - s_star)
      end do
    else
      ! F_*R
      do i = 1, nw_flux
        Dstar = 0._dp
        if (i .eq. iw_mom(flux_dim)) Dstar = 1._dp
        if (i .eq. iw_e)             Dstar = s_star
        F(i) = (s_star*(sR*uR(i) - flux_r(i)) + sR*PLR*Dstar) / (sR - s_star)
      end do
    end if

  end subroutine riemann_hllc_prim


  !> MUSCL (primitive-variable) reconstruction with slope limiter; HLL two-wave approximate Riemann flux at faces.
  !> Uses estimated left/right signal speeds (Davis (1988)) for less diffusion than LLF, no contact resolution.
  subroutine reconflux_muscl_hll_prim(u, xlocC, flux_dim, flux, typelim)
    !$acc routine seq
    real(dp), intent(in)  :: u(nw_phys, 5)
    real(dp), intent(in)  :: xlocC(1:ndim, 2)
    integer, intent(in)   :: flux_dim, typelim
    real(dp), intent(out) :: flux(nw_flux, 2)

    real(dp) :: uL(nw_phys,2), uR(nw_phys,2)
    integer  :: iface

    call muscl_reconstruct_prim(u, typelim, uL, uR)

    do iface=1,2
      call riemann_hll_prim(uL(:,iface), uR(:,iface), xlocC(:,iface), flux_dim, flux(:,iface))
    end do
  end subroutine reconflux_muscl_hll_prim


  !> MUSCL (primitive-variable) reconstruction with slope limiter; LLF/Rusanov numerical flux at faces.
  !> Robust and diffusive; uses local max wavespeed for upwinding.
  subroutine reconflux_muscl_llf_prim(u, xlocC, flux_dim, flux, typelim)
    !$acc routine seq
    real(dp), intent(in)  :: u(nw_phys, 5)
    real(dp), intent(in)  :: xlocC(1:ndim, 2)
    integer, intent(in)   :: flux_dim, typelim
    real(dp), intent(out) :: flux(nw_flux, 2)

    real(dp) :: uL(nw_phys,2), uR(nw_phys,2)
    integer  :: iface

    call muscl_reconstruct_prim(u, typelim, uL, uR)

    do iface=1,2
      call riemann_llf_prim(uL(:,iface), uR(:,iface), xlocC(:,iface), flux_dim, flux(:,iface))
    end do

  end subroutine reconflux_muscl_llf_prim  


  !> MUSCL (primitive-variable) reconstruction with slope limiter; HLLC approximate Riemann flux at faces.
  !> Restores the contact wave (and shear in Euler/HD), typically sharper than HLL for similar cost.
  subroutine reconflux_muscl_hllc_prim(u, xlocC, flux_dim, flux, typelim)
    !$acc routine seq
    real(dp), intent(in)  :: u(nw_phys, 5)
    real(dp), intent(in)  :: xlocC(1:ndim, 2)
    integer, intent(in)   :: flux_dim, typelim
    real(dp), intent(out) :: flux(nw_flux, 2)

    real(dp) :: uL(nw_phys,2), uR(nw_phys,2)
    integer  :: iface

    call muscl_reconstruct_prim(u, typelim, uL, uR)
    do iface=1,2
      call riemann_hllc_prim(uL(:,iface), uR(:,iface), xlocC(:,iface), flux_dim, flux(:,iface))
    end do
  end subroutine reconflux_muscl_hllc_prim


  pure real(dp) function vanleer(a, b) result(phi)
    !$acc routine seq
    real(dp), intent(in) :: a, b
    real(dp)             :: ab

    ab = a * b
    if (ab > 0) then
       phi = 2 * ab / (a + b)
    else
       phi = 0
    end if
  end function vanleer

  pure real(dp) function minmod(a, b)
    !$acc routine seq
    real(dp), intent(in) :: a, b

    if (a * b <= 0) then
       minmod = 0.0_dp
    else if (abs(a) < abs(b)) then
       minmod = a
    else
       minmod = b
    end if
  end function minmod

end module mod_finite_volume
