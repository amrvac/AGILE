#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute

module srmhd_con2prim
use mod_physics_vars

implicit none
real(dp), parameter, public :: h0=1+1e-14_dp

${GPU_DECLARE_COPYIN('h0')}$


!this module implements a conserved to primitive variable
!scheme from Kastaun et al. https://arxiv.org/pdf/2005.01821
!the code will be labled with equation and section numbers see the paper for reference

contains 

!=============================================================================
!
! This section implements brents root solver for use in the con2prim algorithm
!
!============================Brent Root Solver=================================
    pure function quad_inter(x1, f1, f2, f3) result(quad)
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: x1, f1, f2, f3
        real(dp) :: quad

        quad=x1*f2*f3/(f1-f2)/(f1-f3)
    end function quad_inter

    pure subroutine swap(a,b)
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(inout) :: a, b
        real(dp) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap 

    #: def brent_template(name, func)
    function ${name}$(mu_plus, params, calls) result(s)
        ${GPU_ROUTINE_SEQ()}$
        implicit none
        real(dp)  :: s
        !real(dp), intent(in) :: d, tau
        !real(dp), intent(in) :: s_sqr, b_sqr, s_dot_b
        real(dp), intent(in) :: params(:)
        real(dp), intent(in) :: mu_plus
        real(dp) :: a, b, c, d
        real(dp) :: fa, fb, fc, fs
        real(dp) :: tol
        logical :: bound_error
        logical :: mflag
        integer :: i, maxiter
        real(dp):: delta
        integer, intent(out), optional :: calls


        calls=0
        delta = 2.0_dp * epsilon(1.0_dp) * abs(b) + 0.5_dp * tol
        bound_error=.false.
        maxiter=100

        a=0
        b=mu_plus

        tol=1e-15_dp

        fa=${func}$(a, params)
        calls=calls+1

        fb=${func}$(b, params)
        calls=calls+1

        if (fa*fb>0) then
            bound_error=.true.
        end if

        if (abs(fa) < abs(fb)) then
            call swap(a, b)
            call swap(fa, fb)
        end if
        
        c=a
        fc=fa
        mflag=.true.

        do i=1, maxiter
            if (bound_error) then
                !print*, bound_error
                exit
            end if

            if (fa/=fc .and. fb/=fc) then
                s= quad_inter(a, fa, fb, fc) + quad_inter(b, fb, fa, fc) + quad_inter(c, fc, fb, fa)
            else
                s = b-fb*(b-a)/(fb-fa)
            end if

            if (.not. ( min((3.0_dp*a+b)/4.0_dp, b) < s .and. s < max((3.0_dp*a+b)/4.0_dp, b) ) ) then
                s=(a+b)/2
                mflag=.true.
            elseif (mflag) then

                if (abs(s-b)>=abs(c-b)/2 .or. abs(c-b)<delta) then
                    s=(b+a)/2
                    mflag=.true.
                else
                    mflag=.false.
                end if
            else

                if (abs(s-b) >= abs(c-d)/2 .or. abs(c-d)<delta) then
                    s=(b+a)/2
                    mflag=.true.
                else
                    mflag=.false.
                end if
            end if

            fs=${func}$(s, params)
            calls=calls+1

            d=c
            c=b
            fc=fb

            if (fa*fs < 0) then
                b=s
                fb=fs
            else
                a=s
                fa=fs
            end if

            if (abs(b-a)/s<tol) then
            !if (fs<tol) then
                exit
            end if

            if (abs(fa) < abs(fb)) then
            call swap(a, b)
            call swap(fa, fb)
            end if
            delta = 2.0_dp * epsilon(1.0_dp) * abs(b) + 0.5_dp * tol

        end do


    end function ${name}$

    #:enddef

    !$:brent_template("brent","master_function")
    !$:brent_template("aux_brent", "aux_f")

!=============================================================================
!
! This section implements a bisection root solver for use in the con2prim algorithm
!
!============================Bisection Root Solver=================================

    #:def bisection_template(name, func)
    function ${name}$(mu_plus, params, iter) result(mu)
        ${GPU_ROUTINE_SEQ()}$
        implicit none
        !real(dp), intent(in) :: d, tau
        !real(dp), intent(in) :: s_sqr, b_sqr, s_dot_b
        real(dp), intent(in) :: mu_plus
        !logical, intent(out), optional :: success
        integer, intent(out), optional :: iter
        real(dp), intent(in) :: params(:)
        real(dp) :: mu
        real(dp) :: a, b, c
        real(dp) :: fa, fb, fc
        real(dp) :: tol
        integer :: maxiter, i
        logical :: bound_error




        bound_error=.false.
        maxiter=100
        tol=1e-15

        a=0

        b=mu_plus

        fa= ${func}$(a, params)

        fb=${func}$(b, params)

        if (fa>0 .eqv. fb>0) then 
            bound_error=.true.
        end if

        do i=1, maxiter
            if (bound_error) then
                exit
            end if

            c=(b+a)/2.0_dp

            fc=${func}$(c, params)

            if (abs((b-a)/c)<tol) then
                mu=c

                iter=i
                exit
            end if

            if (fa>0 .and. fc>0) then
                a=c
                fa=fc

            else if(.not. fa>0 .and. .not. fc>0) then
                a=c
                fa=fc
            else 
                b=c
                fb=fc
            end if


        end do

        if (.not. bound_error) then
            mu=c
        end if
    
    end function ${name}$
    #:enddef

    $:bisection_template("bisection", "master_function")
    $:bisection_template("aux_bis","aux_f")




!=============================================================================
!
! This is the con2prim procedure
!
!==========================Kastaun con2prim===================================



    function ideal_eos(rho,eps) result(p)
        ${GPU_ROUTINE_SEQ()}$
        real(dp), intent(in) :: rho, eps
        real(dp) ::p

        !srmhd_gamma=5.0_dp/3.0_dp
        p=(srmhd_gamma-1)*rho*eps
    end function ideal_eos



    function master_function(mu, params) result(f)
        ${GPU_ROUTINE_SEQ()}$
        implicit none
        real(dp), intent(in) :: mu
        real(dp) :: d, tau, s_sqr, b_sqr, s_dot_b
        real(dp), intent(in) :: params(5)
        real(dp)  ::  q, x
        real(dp) :: r_sqr, myb_sqr, r_dot_b
        real(dp) :: r_bar_sqr, q_bar
        real(dp) :: v0_sqr
        real(dp) :: v_hat_sqr, w_hat
        real(dp) :: rho_hat, eps_hat
        real(dp) :: p_hat
        real(dp) :: a_hat
        real(dp) :: nu_a, nu_b
        real(dp) :: nu_hat
        real(dp) :: f

        
        d = params(1)
        tau = params(2)
        s_sqr = params(3)
        b_sqr = params(4) 
        s_dot_b = params(5)

        r_sqr=s_sqr/d**2 !(23)

        q=tau/d !(23)

        myb_sqr=b_sqr/d !(24)

        r_dot_b=s_dot_b !necessary for (25)

        v0_sqr=r_sqr/h0**2/(1 + r_sqr/h0**2) ! (32)(33)

        x=1.0_dp/(1.0_dp + mu*myb_sqr) !(26)

        r_bar_sqr= r_sqr*x**2 + mu*x*(1+x)*r_dot_b**2/d**3 !(38)

        q_bar= q - 0.5_dp*myb_sqr - 0.5_dp*mu**2*x**2*(r_sqr*myb_sqr-r_dot_b**2/d**3) !(39)

        v_hat_sqr=min(v0_sqr, mu**2*r_bar_sqr) !(40)

        w_hat=1/sqrt(1 - v_hat_sqr) !(40)

        rho_hat=d/w_hat !(41)

        
        if (rho_hat<0) then

            rho_hat=0
        end if

        eps_hat=w_hat*(q_bar - mu*r_bar_sqr) + v_hat_sqr*w_hat**2/(1 + w_hat) !(42)

        !add a limiter for valid eps and rho
        if (eps_hat<0) then
            eps_hat=0
        end if

        if (rho_hat<0) then
            rho_hat=0
        end if

        p_hat = ideal_eos(rho_hat,eps_hat) !(43)

        a_hat=p_hat/rho_hat/(1 + eps_hat) !(43)

        nu_a=(1 + a_hat)*(1 + eps_hat)/w_hat !(46)

        nu_b=(1 + a_hat)*(1 + q_bar - mu*r_bar_sqr) !(47)
         
        nu_hat=max(nu_a, nu_b) !(48)

        f=mu-1.0_dp/(nu_hat + mu*r_bar_sqr)!(44)(45)


    end function master_function



        


    function aux_f(mu, params) result(f)
        ${GPU_ROUTINE_SEQ()}$
        implicit none
        real(dp), intent(in) :: mu
        real(dp), intent(in) :: params(4)
        real(dp)             :: d, s_sqr, s_dot_b, b_sqr
        real(dp)             :: mymu
        real(dp)             :: r_dot_b, myb_sqr, myd
        real(dp)             :: f
        real(dp)             :: r_bar_sqr, r_sqr, x


        d = params(1)
        s_sqr = params(2)
        b_sqr = params(3)
        s_dot_b = params(4)

        mymu=mu
        r_dot_b=s_dot_b
        myd=d
        myb_sqr=b_sqr/myd
        r_sqr=s_sqr/myd**2

        x=1/(1+mu*myb_sqr)

        r_bar_sqr= r_sqr*x**2 + mu*x* (1+x) *r_dot_b**2/myd**3 !(38)

        f=mu*sqrt(h0**2+r_bar_sqr)-1
    end function aux_f


    #:def con2prim_template(name, solver, aux_solver)
    subroutine ${name}$(mu, d, tau, s_sqr, b_sqr, s_dot_b)
        ${GPU_ROUTINE_SEQ()}$
        implicit none
        !real(dp), intent(out) :: lfac
        !real(dp), intent(out), optional:: rho, eps, p, h
        !logical, intent(out), optional :: success
        !integer, intent(out), optional :: iter, aux_iter
        !real(dp), intent(out) :: vi(3)
        real(dp), intent(in)  :: d, tau
        real(dp), intent(in)  :: s_sqr, b_sqr, s_dot_b
        real(dp)              :: myd, mytau
        !real(dp)              :: mysi(3), mybi(3)
        !real(dp)              :: s_sqr, b_sqr, s_dot_b
        real(dp), intent(out):: mu
        real(dp)              :: mu_plus
        real(dp)              :: q_bar, r_bar_sqr
        real(dp)              :: x
        real(dp)              :: params(5)
        real(dp)              :: aux_params(4)
        

        myd=d
        mytau=tau


        aux_params=[myd, s_sqr, b_sqr, s_dot_b]

        if (s_sqr/myd**2 < h0**2) then
            mu_plus=h0
        else
            mu_plus=${aux_solver}$(1/h0, aux_params)
        end if

        
        params=[myd, mytau, s_sqr, b_sqr, s_dot_b]

        mu= ${solver}$(mu_plus, params)


    end subroutine ${name}$
    #:enddef 

    $:con2prim_template("con2prim", "bisection", "aux_bis")

    

end module srmhd_con2prim



