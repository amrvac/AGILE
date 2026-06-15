module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision:: v0,rho0,p0,T0,pbeta0,b0,mach0

!$acc declare create(v0,rho0,p0,T0,pbeta0,b0,mach0)

contains

  subroutine usr_init()
    use mod_global_parameters
    implicit none

    call usr_params_read(par_files)

    unit_length        = 1.d9 ! in cm (1 Mm in cm)
    unit_temperature   = 1.d6 ! in K (1Mk)
    unit_numberdensity = 1.d9 ! in cm^-3

    call set_coordinate_system("Cartesian_3D")
    usr_set_parameters => set_parameters_usr
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    call phys_activate()
  end subroutine usr_init

  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ pbeta0,T0,rho0,mach0

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

    if(mach0<smalldouble) call mpistop("Wrong input values: mach0 must be positive")
    if(rho0<smalldouble) call mpistop("Wrong input values: rho0 must be positive")
    if(T0<smalldouble) call mpistop("Wrong input values: T0 must be positive")
    p0=rho0*T0
    if(pbeta0>smalldouble)then
        b0=dsqrt(2.0d0*p0/pbeta0)
    else
        b0=zero
    endif
    v0=mach0*dsqrt(mhd_gamma*p0/rho0)

!$acc update device(v0,rho0,p0,T0,pbeta0,b0,mach0)

  end subroutine usr_params_read

  subroutine print_params()
    implicit none
    character(len=*), parameter :: fmt = '(A20,1X,SP,ES12.2)'

    print *, '---Orszag-Tang--- PARAMETERS ------------------'
    print fmt,  'mhd_gamma',    mhd_gamma
    print fmt,  'v0',           v0
    print fmt,  'v0',           v0
    print fmt,  'rho0',         rho0
    print fmt,  'p0',           p0
    print fmt,  'T0',           T0
    print fmt,  'pbeta0',       pbeta0
    print fmt,  'b0',           b0
    print fmt,  'mach0',        mach0
    print *, '------------------------------------------------'
  end subroutine print_params

    subroutine set_parameters_usr()
    use mod_global_parameters
    use mod_physics
    implicit none

    if (mype == 0) call print_params()
  end subroutine set_parameters_usr

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    use mod_global_parameters
    use mod_physics
    implicit none
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_) = rho0

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(1)) = &
       -sin(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)) &
       *cos(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,3))

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(2)) = &
        sin(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)) &
       *cos(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,3))

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(3)) = zero

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_) = p0

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(1)) = &
       -b0 * sin(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)) &
       *cos(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,3))
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(2)) = &
        b0 * sin(4.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)) &
       *cos(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,3))
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(3)) = zero

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

! NOTE: This must be named exactly `addsource_usr` in AGILE since during
!   compilation it deals with this function.
! NOTE: This is ran on the GPU.
  subroutine addsource_usr(qdt, qt, wCT, wCTprim, wnew, x, split)
      !$acc routine seq
      use mod_radiative_cooling, only: rc_fl,getvar_cooling_exact
      double precision, intent(in) :: qdt, qt, wCT(nw_phys), wCTprim(nw_phys)
      double precision, intent(inout) :: wnew(nw_phys)
      double precision, intent(in) :: x(1:ndim)
      logical, intent(in) :: split

      double precision :: wCT_nopert(nw_phys),w_nopert(nw_phys)
      double precision :: coolrate

      ! This adds a uniform background heating corresponding to the t=0 uniform unmagnetized background
      wCT_nopert(iw_rho)=rho0
      wCT_nopert(iw_mom(1))=0.0d0
      wCT_nopert(iw_mom(2))=0.0d0
      wCT_nopert(iw_mom(3))=0.0d0
      wCT_nopert(iw_e)=p0/(mhd_gamma-1.0d0)
      wCT_nopert(iw_mag(1))=0.0d0
      wCT_nopert(iw_mag(2))=0.0d0
      wCT_nopert(iw_mag(3))=0.0d0
      w_nopert(1:nw_phys)=wCT_nopert(1:nw_phys)
      call getvar_cooling_exact(qdt,wCT_nopert,w_nopert,x,coolrate,rc_fl)
      wnew(iw_e)=wnew(iw_e)+qdt*coolrate

  end subroutine addsource_usr



  subroutine specialvar_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
    use mod_global_parameters
    use mod_physics
    implicit none
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw+nwauxio)
    double precision             :: normconv(0:nw+nwauxio)

    double precision :: wlocal(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    wlocal(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)= &
       w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    call phys_to_primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wlocal,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1) = &
       wlocal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_e) &
      /wlocal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_rho)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames = 'T'
  end subroutine specialvarnames_output

end module mod_usr

