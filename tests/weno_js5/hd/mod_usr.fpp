#:mute
#:include "mod_gpu_directives.fpp"
#:endmute
module mod_usr
  use mod_amrvac
  use mod_physics
  implicit none
contains
  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr
    usr_print_log => diagnostics
    call phys_activate()
  end subroutine

  subroutine initonegrid_usr(g1,g2,g3,h1,h2,h3,l1,l2,l3,r1,r2,r3,w,x)
    integer,intent(in) :: g1,g2,g3,h1,h2,h3,l1,l2,l3,r1,r2,r3
    double precision,intent(in) :: x(g1:h1,g2:h2,g3:h3,ndim)
    double precision,intent(inout) :: w(g1:h1,g2:h2,g3:h3,nw)
    integer :: i,j,k,d,axis
    double precision :: rho,p,v(3),sinc,phase,dx(3)
    axis=1
    if(abs(iprob)>=1.and.abs(iprob)<=3) axis=abs(iprob)
    dx=[x(g1+1,g2,g3,1)-x(g1,g2,g3,1), &
        x(g1,g2+1,g3,2)-x(g1,g2,g3,2),x(g1,g2,g3+1,3)-x(g1,g2,g3,3)]
    w=0.d0
    do k=l3,r3
    do j=l2,r2
    do i=l1,r1
      v=0.d0
      v(axis)=1.d0
      if(iprob<0) v(axis)=-1.d0
      phase=x(i,j,k,axis)
      sinc=sin(dpi*dx(axis))/(dpi*dx(axis))
      if(iprob==4) then
        v=1.d0
        phase=sum(x(i,j,k,:))
        sinc=product(sin(dpi*dx)/(dpi*dx))
      end if
      rho=1.d0+0.2d0*sinc*sin(2*dpi*phase)
      p=1.d0
      if(iprob==7) then
        rho=1.d0
        v=0.d0
      else if(iprob==8) then
        rho=1.d0+0.2d0*cos(2*dpi*x(i,j,k,1))
        v=0.d0
        v(1)=0.1d0*sin(2*dpi*x(i,j,k,1))
      else if(iprob==10) then
        v=0.d0
        if(x(i,j,k,1)<0.5d0) then
          rho=1.d0
          p=1.d0
        else
          rho=0.125d0
          p=0.1d0
        end if
      else if(iprob==11) then
        v=0.d0
        if(x(i,j,k,1)<-4.d0) then
          rho=3.857143d0
          v(1)=2.629369d0
          p=10.33333d0
        else
          rho=1.d0+0.2d0*sin(5.d0*x(i,j,k,1))*sin(2.5d0*dx(1))/(2.5d0*dx(1))
        end if
      end if
      w(i,j,k,rho_)=rho
      w(i,j,k,mom)=rho*v
      w(i,j,k,e_)=p/(hd_gamma-1.d0)+0.5d0*rho*sum(v*v)
    end do
    end do
    end do
  end subroutine

  subroutine usr_refine_grid(level,qt,w,x,rf,cf,nrf,ncf)
    ${GPU_ROUTINE_SEQ()}$
    integer,intent(in) :: level
    double precision,intent(in) :: qt,w(nw_phys),x(3)
    logical,intent(inout) :: rf,cf,nrf,ncf
    double precision :: center
    center=0.4d0
    if(iprob==5) center=0.25d0+4.d0*qt
    rf=abs(x(1)-center)<0.16d0
    cf=.not.rf
    nrf=.not.rf
    ncf=rf
  end subroutine

  subroutine diagnostics()
    integer :: ib,n,i,j,k,axis,u
    character(len=80) :: name
    double precision :: local(7),total(7),emax,globmax,pmin,gpmin,vol,rho,exact,dx(3),phase,sinc,p
    axis=1
    if(abs(iprob)>=1.and.abs(iprob)<=3) axis=abs(iprob)
    local=0.d0
    emax=0.d0
    pmin=huge(1.d0)
    write(name,'(a,i4.4,a)') 'cells-',mype,'.txt'
    open(newunit=u,file=trim(name),status='replace')
    do ib=1,igridstail
      n=igrids(ib)
      dx=rnode(rpdx1_:rpdx1_+2,n)
      vol=product(dx)
      do k=ixMlo3,ixMhi3
      do j=ixMlo2,ixMhi2
      do i=ixMlo1,ixMhi1
        rho=ps(n)%w(i,j,k,rho_)
        phase=ps(n)%x(i,j,k,axis)-global_time
        if(iprob<0) phase=ps(n)%x(i,j,k,axis)+global_time
        sinc=sin(dpi*dx(axis))/(dpi*dx(axis))
        if(iprob==4) then
          phase=sum(ps(n)%x(i,j,k,:))-3.d0*global_time
          sinc=product(sin(dpi*dx)/(dpi*dx))
        end if
        exact=1.d0+0.2d0*sinc*sin(2*dpi*phase)
        if(iprob==7) exact=1.d0
        p=(hd_gamma-1.d0)*(ps(n)%w(i,j,k,e_)-0.5d0*sum(ps(n)%w(i,j,k,mom)**2)/rho)
        pmin=min(pmin,p)
        local(1)=local(1)+abs(rho-exact)*vol
        emax=max(emax,abs(rho-exact))
        local(2:6)=local(2:6)+ps(n)%w(i,j,k,1:5)*vol
        local(7)=local(7)+vol
        write(u,'(10(es25.16,1x))') ps(n)%x(i,j,k,:),vol,ps(n)%w(i,j,k,1:5),p
      end do
      end do
      end do
    end do
    close(u)
    call MPI_ALLREDUCE(local,total,7,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(emax,globmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)
    call MPI_ALLREDUCE(pmin,gpmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)
    if(mype==0) then
      open(newunit=u,file='diagnostics.txt',position='append')
      write(u,'(11(es25.16,1x))') global_time,total(1)/total(7),globmax,total(2:7),gpmin,dble(nghostcells)
      close(u)
    end if
  end subroutine
end module mod_usr
