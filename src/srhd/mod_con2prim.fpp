module mod_con2prim

use mod_global_parameters
use mod_physics_vars
implicit none

! Note : we need phys_gamma and use small_pressure

  !> parameters for NR in con2prim
  integer, public                  :: maxitnr   = 100
  double precision, public         :: absaccnr  = 1.0d-8
  double precision, public         :: tolernr   = 1.0d-9
  double precision, public         :: dmaxvel   = 1.0d-7
  
  !$acc declare copyin(maxitnr,absaccnr,tolernr,dmaxvel)


contains 

  !> con2prim: (D,S**2,tau) --> compute auxiliaries lfac and xi
  pure subroutine con2prim_eos(lfac,xi,myd,myssqr,mytau)
    !$acc routine seq
    double precision, intent(in)    :: myd, myssqr, mytau
    double precision, intent(inout) :: lfac, xi

    ! .. local ..
    double precision     :: d,ssqr,tau,vsqr
    double precision     :: f,df,lfacl,dlfac,p,dpdxi
    !------------------------------------------------------------------

    d = myd; ssqr = myssqr; tau = mytau;

    vsqr = ssqr/xi**2
    lfacl = one/dsqrt(one-vsqr)
    dlfac = -lfacl**3*ssqr/(xi**3)
    call FuncPressure_eos(xi,lfacl,d,dlfac,p,dpdxi)
    f  = xi-tau-d-p
    df = one-dpdxi
    if (dabs(f/df)<absaccnr) then
       xi   = xi - f/df
       lfac = lfacl
       return
    end if

    call con2primHydro_eos(lfac,xi,d,ssqr,tau)

  end subroutine con2prim_eos

  !> SRHD iteration solves for p via NR, and then gives xi as output
  pure subroutine con2primHydro_eos(lfac,xi,d,sqrs,tau)
    !$acc routine seq
    double precision, intent(out) :: xi,lfac
    double precision, intent(in)  :: d,sqrs,tau

    ! .. local ..
    integer          :: ni,niiter,nit,n2it,ni3
    double precision :: pcurrent,pnew
    double precision :: er,er1,ff,df,dp,v2
    double precision :: pmin,lfac2inv,pLabs,pRabs,pprev
    double precision :: dv2d2p
    double precision :: xicurrent,h,dhdp
    double precision :: oldff1,oldff2,Nff
    double precision :: pleft,pright
    double precision :: rho,E_th,E
    double precision :: dE_thdp,dEdp
    !---------------------------------------------------------------------

    ! left and right brackets for p-range
    pmin=dsqrt(sqrs)/(one-dmaxvel)-tau-d
    pLabs=max(small_pressure,pmin)
    pRabs=1.0d99
    ! start value from input
    pcurrent=pLabs

    er1=one
    pprev=pcurrent

    ! Fudge Parameters
    oldff1=1.0d7  ! High number
    oldff2=1.0d9  ! High number bigger then oldff1
    n2it = 0
    nit  = 0

    LoopNR:  do ni=1,maxitnr
       nit = nit + 1
       !=== Relax NR iteration accuracy=======!
       if(nit>maxitnr/4)then
          ! mix pressure value for convergence
          pcurrent=half*(pcurrent+pprev)
          ! relax accuracy requirement
          er1=10.0d0*er1
          nit = nit - maxitnr/10
       endif
       !=======================================!

       niiter=ni
       xicurrent=tau+d+pcurrent

       v2=sqrs/xicurrent**2
       lfac2inv=one - v2
          lfac=one/dsqrt(lfac2inv)

       dv2d2p=sqrs/(xicurrent**3)
       !== calculation done using the EOS ==!
       rho=d*dsqrt(lfac2inv)
       E_th = pcurrent/(phys_gamma-1.0d0)
       E = (E_th + dsqrt(E_th**2+rho**2))
       !== Enthalpy ==!
       h = half*((phys_gamma+one)*E-(phys_gamma-1.0d0)*rho*(rho/E))
       !=== Derivative of thermal energy ===!
       dE_thdp = one/(phys_gamma-1.0d0)
       !=== Derivative of internal energy ===!
       dEdp = dE_thdp * (one+E_th/dsqrt(E_th**2+rho**2))&
              +  d**2*dv2d2p/dsqrt(E_th**2+rho**2)
       !====== Derivative of Enthalpy ======!
       dhdp = half*((phys_gamma+one)*dEdp + &
              (phys_gamma-1.0d0)*(rho*(rho/E))*(-2.0d0*dv2d2p/lfac2inv+dEdp/E))
       !=======================================!
       ff=-xicurrent*lfac2inv + h
       df=- two*sqrs/xicurrent**2  + dhdp - lfac2inv

       if (ff*df==zero) then
          if (ff==zero) then
             exit ! zero found
          endif
       else
          pnew=pcurrent-ff/df
          if (ff*df>zero) then
             ! pressure iterate has decreased
             ! restrict to left
             pnew=max(pnew,pLabs)
          else  ! ff*df<0
             ! pressure iterate has increased
             ! restrict to right
             pnew=min(pnew,pRabs)
          endif
       endif

       !===============================================!
       dp=pcurrent-pnew
       er=two*dabs(dp)/(pnew+pcurrent)
       if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
       !===============================================!

       ! For very small values of pressure, NR algorithm is not efficient to
       ! find root, use Euler algorithm to find precise value of pressure
       if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.&
            ff * oldff1 < zero    .and.  dabs(ff)>absaccnr)then

          n2it=n2it+1
          if(n2it<=3) pcurrent=half*(pnew+pcurrent)
          if(n2it>3)then
             pright =pcurrent
             pleft=pprev
             pcurrent=half*(pleft+pright)
             Dicho:  do ni3=1,maxitnr
                !===================!
                xicurrent=tau+d+pcurrent
                v2=sqrs/xicurrent**2
                lfac2inv=one - v2

                   lfac=one/dsqrt(lfac2inv)

                !== calculation done using the EOS ==!
                rho=d*dsqrt(lfac2inv)
                E_th = pnew/(phys_gamma-1.0d0)
                E = (E_th + dsqrt(E_th**2+rho**2))
                !== Enthalpy ==!
                h = half*((phys_gamma+one)*E-(phys_gamma-1.0d0)*rho*(rho/E))
                Nff=-xicurrent*lfac2inv + h
                !=======================================!
                !==== Iterate ====!
                if(ff * Nff < zero)then
                   pleft=pcurrent
                else
                   pright=pcurrent
                endif

                pcurrent=half*(pleft+pright)
                !==================!

                !=== The iteration converged ===!
                if(2.0d0*dabs(pleft-pright)/(pleft+pright)< absaccnr &
                     .or. dabs(ff)<absaccnr)then
                   pnew=pcurrent
                   exit LoopNR
                endif
                !==============================!

                !==============================!

                !=== conserve the last value of Nff ===!
                ff=Nff
                !======================================!
             enddo    Dicho
          endif

       else
          !====== There is no problems, continue the NR iteration ======!
          pprev=pcurrent
          pcurrent=pnew
          !=============================================================!
       endif


       !=== keep the values of the 2 last ff ===!
       oldff2=oldff1
       oldff1=ff
       !========================================!
    enddo LoopNR

    !------------------------------!
    xi=tau+d+pcurrent
    v2=sqrs/xicurrent**2
    lfac2inv=one - v2
       lfac=one/dsqrt(lfac2inv)
  end subroutine con2primHydro_eos

  !> pointwise evaluations used in con2prim
  !> compute pointwise value for pressure p and dpdxi
  pure subroutine FuncPressure_eos(xicurrent,lfac,d,dlfacdxi,p,dpdxi)
    !$acc routine seq
    double precision, intent(in)         :: xicurrent,lfac,d,dlfacdxi
    double precision, intent(out)        :: p,dpdxi
    ! .. local ..
    double precision                     :: rho,h,E,dhdxi,rhotoE
    double precision                     :: dpdchi,dEdxi

    ! rhoh here called h
    h=xicurrent/(lfac**2)
    rho=d/lfac
    E = (h+dsqrt(h**2+(phys_gamma**2-one)*rho**2)) &
              /(phys_gamma+one)
    ! output pressure
    rhotoE = rho/E
    p = half*(phys_gamma-1.0d0)*(E-rho*rhotoE)

    dhdxi = one/(lfac**2)-2.0d0*xicurrent/(lfac**2)*dlfacdxi/lfac

    dEdxi=(dhdxi+(h*dhdxi-(phys_gamma**2-one)*rho**2*dlfacdxi/lfac)&
        /dsqrt(h**2+(phys_gamma**2-one)*rho**2))&
        /(phys_gamma+one)

    ! output pressure derivative to xi
    dpdxi=half*(phys_gamma-1.0d0)*(2.0d0*rho*rhotoE*dlfacdxi/lfac+&
          (one+rhotoE**2)*dEdxi)

  end subroutine FuncPressure_eos

  !> con2prim: (D,S**2,tau) --> compute auxiliaries lfac and xi
  pure subroutine con2prim(lfac,xi,myd,myssqr,mytau)
    !$acc routine seq
    double precision, intent(in)    :: myd, myssqr, mytau
    double precision, intent(inout) :: lfac, xi

    ! .. local ..
    double precision:: f,df,lfacl
    double precision       :: d,ssqr,tau
    !------------------------------------------------------------------

    d = myd; ssqr = myssqr; tau = mytau;

    ! Check if guess is close enough: gives f,df,lfacl
    call funcd(xi,f,df,lfacl,d,ssqr,tau)
    if (dabs(f/df)<absaccnr) then
       xi   = xi - f/df
       lfac = lfacl
       return
    end if

    call con2primHydro(lfac,xi,d,ssqr,tau)

  end subroutine con2prim


  pure subroutine funcd(xi,f,df,mylfac,d,ssqr,tau)
    !$acc routine seq
    double precision, intent(in)  :: xi,d,ssqr,tau
    double precision, intent(out) :: f,df,mylfac

    ! .. local ..
    double precision  :: dlfac
    double precision  :: vsqr,p,dpdxi
    !-----------------------------------------------------------------

    vsqr = ssqr/xi**2

       mylfac = one/dsqrt(one-vsqr)
       dlfac = -mylfac**3*ssqr/(xi**3)
       !===== Pressure, calculate using EOS =====!
       call FuncPressure(xi,mylfac,d,dlfac,p,dpdxi)
       !=========================================!
       f  = xi-tau-d-p
       df = one-dpdxi

  end subroutine funcd

  !> SRHD iteration solves for p via NR, and then gives xi as output
  pure subroutine con2primHydro(lfac,xi,d,sqrs,tau)
    !$acc routine seq
    double precision, intent(out) :: xi,lfac
    double precision, intent(in)  :: d,sqrs,tau

    ! .. local ..
    integer          :: ni,niiter,nit,n2it,ni3
    double precision :: pcurrent,pnew
    double precision :: er,er1,ff,df,dp,v2
    double precision :: pmin,lfac2inv,pLabs,pRabs,pprev
    double precision :: xicurrent,h,dhdp
    double precision :: oldff1,oldff2,Nff
    double precision :: pleft,pright
    double precision :: rho
    !---------------------------------------------------------------------

    ! left and right brackets for p-range
    pmin=dsqrt(sqrs)/(one-dmaxvel)-tau-d
    pLabs=max(small_pressure,pmin)
    pRabs=1.0d99
    ! start value from input
    pcurrent=pLabs

    er1=one
    pprev=pcurrent

    ! Fudge Parameters
    oldff1=1.0d7  ! High number
    oldff2=1.0d9  ! High number bigger then oldff1
    n2it = 0
    nit  = 0

    LoopNR:  do ni=1,maxitnr
       nit = nit + 1
       !=== Relax NR iteration accuracy=======!
       if(nit>maxitnr/4)then
          ! mix pressure value for convergence
          pcurrent=half*(pcurrent+pprev)
          ! relax accuracy requirement
          er1=10.0d0*er1
          nit = nit - maxitnr/10
       endif
       !=======================================!

       niiter=ni
       xicurrent=tau+d+pcurrent

       v2=sqrs/xicurrent**2
       lfac2inv=one - v2
          lfac=one/dsqrt(lfac2inv)

       !== calculation done using the EOS ==!
       rho=d*dsqrt(lfac2inv)
       h = rho + pcurrent*phys_gamma/(phys_gamma-1.0d0)
       dhdp = phys_gamma/(phys_gamma-1.0d0) + d/dsqrt(lfac2inv)*sqrs/xicurrent**3
       !=======================================!
       ff=-xicurrent*lfac2inv + h
       df=- two*sqrs/xicurrent**2  + dhdp - lfac2inv

       if (ff*df==zero) then
          if (ff==zero) then
             exit ! zero found
          endif
       else
          pnew=pcurrent-ff/df
          if (ff*df>zero) then
             ! pressure iterate has decreased
             ! restrict to left
             pnew=max(pnew,pLabs)
          else  ! ff*df<0
             ! pressure iterate has increased
             ! restrict to right
             pnew=min(pnew,pRabs)
          endif
       endif

       !===============================================!
       dp=pcurrent-pnew
       er=two*dabs(dp)/(pnew+pcurrent)
       if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
       !===============================================!

       ! For very small values of pressure, NR algorithm is not efficient to
       ! find root, use Euler algorithm to find precise value of pressure
       if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.&
            ff * oldff1 < zero    .and.  dabs(ff)>absaccnr)then

          n2it=n2it+1
          if(n2it<=3) pcurrent=half*(pnew+pcurrent)
          if(n2it>3)then
             pright =pcurrent
             pleft=pprev
             pcurrent=half*(pleft+pright)
             Dicho:  do ni3=1,maxitnr
                !===================!
                xicurrent=tau+d+pcurrent
                v2=sqrs/xicurrent**2
                lfac2inv=one - v2

                   lfac=one/dsqrt(lfac2inv)

                !== calculation done using the EOS ==!
                rho=d*dsqrt(lfac2inv)
                h = rho + pnew*phys_gamma/(phys_gamma-1.0d0)
                Nff=-xicurrent*lfac2inv + h
                !=======================================!
                !==== Iterate ====!
                if(ff * Nff < zero)then
                   pleft=pcurrent
                else
                   pright=pcurrent
                endif

                pcurrent=half*(pleft+pright)
                !==================!

                !=== The iteration converged ===!
                if(2.0d0*dabs(pleft-pright)/(pleft+pright)< absaccnr &
                     .or. dabs(ff)<absaccnr)then
                   pnew=pcurrent
                   exit LoopNR
                endif
                !==============================!

                !==============================!

                !=== conserve the last value of Nff ===!
                ff=Nff
                !======================================!
             enddo    Dicho
          endif

       else
          !====== There is no problems, continue the NR iteration ======!
          pprev=pcurrent
          pcurrent=pnew
          !=============================================================!
       endif


       !=== keep the values of the 2 last ff ===!
       oldff2=oldff1
       oldff1=ff
       !========================================!
    enddo LoopNR

    !------------------------------!
    xi=tau+d+pcurrent
    v2=sqrs/xicurrent**2
    lfac2inv=one - v2
       lfac=one/dsqrt(lfac2inv)

  end subroutine con2primHydro

  !> pointwise evaluations used in con2prim
  !> compute pointwise value for pressure p and dpdxi
  pure subroutine FuncPressure(xicurrent,lfac,d,dlfacdxi,p,dpdxi)
    !$acc routine seq

    double precision, intent(in)         :: xicurrent,lfac,d,dlfacdxi
    double precision, intent(out)        :: p,dpdxi
    ! .. local ..
    double precision                     :: rho,h,E,dhdxi,rhotoE
    double precision                     :: dpdchi,dEdxi

    ! rhoh here called h
    h=xicurrent/(lfac**2)
    rho=d/lfac
    ! output pressure
    p = (h - rho)*(phys_gamma-1.0d0)/phys_gamma
    dpdchi = one*(phys_gamma-1.0d0)/phys_gamma
    dpdxi = dpdchi * one/lfac**2
    ! zero case dlfacdxi implies zero velocity (ssqr=0)
    if (dlfacdxi /= 0.0d0) &
          dpdxi = dpdxi  + dpdchi * ((d*lfac-2.0d0*xicurrent)/lfac**3) * dlfacdxi

  end subroutine FuncPressure

end module mod_con2prim
