!> Suresh & Huynh (1997), uniform-cell finite-volume MP5 reconstruction.
#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute
module mod_mp5
  implicit none
  private
  public :: mp5_left
contains
  pure double precision function mm(a,b) result(v)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in) :: a,b
    v=0.5d0*(sign(1.d0,a)+sign(1.d0,b))*min(abs(a),abs(b))
  end function mm

  pure subroutine mp5_left(v,tolerance,value)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in) :: v(5),tolerance
    double precision, intent(out) :: value
    double precision :: candidate,mp,dm,d0,dp,dl,dr,ul,md,lc,lower,upper
    candidate=(2.d0*v(1)-13.d0*v(2)+47.d0*v(3)+27.d0*v(4)-3.d0*v(5))/60.d0
    mp=v(3)+mm(v(4)-v(3),4.d0*(v(3)-v(2)))
    if ((candidate-v(3))*(candidate-mp)<=tolerance) then
      value=candidate
      return
    end if
    dm=v(1)-2.d0*v(2)+v(3)
    d0=v(2)-2.d0*v(3)+v(4)
    dp=v(3)-2.d0*v(4)+v(5)
    dl=mm(mm(4.d0*dm-d0,4.d0*d0-dm),mm(dm,d0))
    dr=mm(mm(4.d0*d0-dp,4.d0*dp-d0),mm(d0,dp))
    ul=v(3)+4.d0*(v(3)-v(2))
    md=0.5d0*(v(3)+v(4)-dr)
    lc=v(3)+0.5d0*(v(3)-v(2))+4.d0/3.d0*dl
    lower=max(min(v(3),v(4),md),min(v(3),ul,lc))
    upper=min(max(v(3),v(4),md),max(v(3),ul,lc))
    value=candidate+mm(lower-candidate,upper-candidate)
  end subroutine mp5_left
end module mod_mp5
