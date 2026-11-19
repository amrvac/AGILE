!> Uniform-cell finite-volume Jiang--Shu (1996) reconstruction.
#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute
module mod_weno_js5
  implicit none
  private
  public :: weno_js5_left
contains
  pure subroutine weno_js5_left(v, epsilon, value, candidates, beta, weights)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in) :: v(5), epsilon
    double precision, intent(out) :: value
    double precision, intent(out), optional :: candidates(3), beta(3), weights(3)
    double precision :: q(3), b(3), a(3), d(3), den(3)
    ! v(3) is the donor cell; reconstruct its right face from cell averages.
    q(1) = (2.d0*v(1)-7.d0*v(2)+11.d0*v(3))/6.d0
    q(2) = (-v(2)+5.d0*v(3)+2.d0*v(4))/6.d0
    q(3) = (2.d0*v(3)+5.d0*v(4)-v(5))/6.d0
    b(1) = (13.d0/12.d0)*(v(1)-2.d0*v(2)+v(3))**2 &
         + 0.25d0*(v(1)-4.d0*v(2)+3.d0*v(3))**2
    b(2) = (13.d0/12.d0)*(v(2)-2.d0*v(3)+v(4))**2 &
         + 0.25d0*(v(2)-v(4))**2
    b(3) = (13.d0/12.d0)*(v(3)-2.d0*v(4)+v(5))**2 &
         + 0.25d0*(3.d0*v(3)-4.d0*v(4)+v(5))**2
    d = [0.1d0, 0.6d0, 0.3d0]
    den = epsilon+b
    ! A common factor cancels on normalization: exactly JS p=2 weights,
    ! without overflowing alpha when a user selects a very small epsilon.
    a = d*(minval(den)/den)**2
    a = a/sum(a)
    value = sum(a*q)
    if (present(candidates)) candidates = q
    if (present(beta)) beta = b
    if (present(weights)) weights = a
  end subroutine weno_js5_left
end module mod_weno_js5
