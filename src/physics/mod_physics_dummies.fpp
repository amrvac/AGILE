! Dummy routines which can be overwritten by a physics-dependent implementation
! On non-Cray compilers, these call mpistop with a descriptive message.
! On Cray, STOP cannot be inlined into OpenACC kernels, so dummies just return -1. 

#:def estimate_speeds_minmax()
subroutine estimate_speeds_minmax(uL, uR, xC, flux_dim, wL, wR)
  !$acc routine seq
#ifndef _CRAYFTN
  use mod_comm_lib, only: mpistop
#endif
  real(dp), intent(in)  :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: wL, wR

#ifndef _CRAYFTN
  call mpistop("estimate_speeds_minmax not implemented for this physics module")
#endif
  wL = -1._dp
  wR = -1._dp

end subroutine estimate_speeds_minmax
#:enddef


#:def estimate_speeds_toro_pvrs()
subroutine estimate_speeds_toro_pvrs(uL, uR, xC, flux_dim, sL, sR)
  !$acc routine seq
#ifndef _CRAYFTN
  use mod_comm_lib, only: mpistop
#endif
  real(dp), intent(in)  :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer,  intent(in)  :: flux_dim
  real(dp), intent(out) :: sL, sR

#ifndef _CRAYFTN
  call mpistop("estimate_speeds_toro_pvrs not implemented for this physics module")
#endif
  sL = -1._dp
  sR = -1._dp

end subroutine estimate_speeds_toro_pvrs
#:enddef


#:def addsource_nonlocal()
subroutine addsource_nonlocal(qdt, dtfactor, qtC, wCTprim, qt, wnew, x, dx, idir, &
     qsourcesplit)
  !$acc routine seq

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCTprim(nw_phys,5)
  real(dp), intent(in)     :: x(1:ndim), dx(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  integer, intent(in)      :: idir
  logical, intent(in)      :: qsourcesplit


end subroutine addsource_nonlocal
#:enddef

#:def addsource_compact()
subroutine addsource_compact(qdt, dtfactor, qtC, wCTprim1, wCTprim2, wCTprim3, qt, wnew, x, dx, &
     qsourcesplit)
  !$acc routine seq

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCTprim1(nw_phys,3),wCTprim2(nw_phys,3),wCTprim3(nw_phys,3)
  real(dp), intent(in)     :: x(1:ndim), dx(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit


end subroutine addsource_compact
#:enddef
