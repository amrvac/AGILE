module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision :: rho0, p0, b0
  ! Host-side copies of B-field w-indices, saved before device update
  integer :: ib1, ib2, ib3

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output
    usr_init_vector_potential => initvecpot_usr

    call set_coordinate_system("Cartesian_3D")

    call phys_activate()

    rho0 = 25.0_dp / (36.0_dp * dpi)
    p0   =  5.0_dp / (12.0_dp * dpi)
    b0   =  1.0_dp / sqrt(4.0_dp * dpi)

    ! Save B-field indices on host (mag becomes device-only after phys_activate)
    ib1 = mag(1); ib2 = mag(2); ib3 = mag(3)

    if (mype == 0) then
      write(*,'(A)')        ' ========================================'
      write(*,'(A)')        '  Orszag-Tang MHD vortex (3D)'
      write(*,'(A,ES12.5)') '    rho0 = ', rho0
      write(*,'(A,ES12.5)') '    p0   = ', p0
      write(*,'(A,ES12.5)') '    b0   = ', b0
      write(*,'(A)')        ' ========================================'
    end if

  end subroutine usr_init


  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_) = rho0

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(1)) = &
       -sin(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(2)) = &
        sin(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1))

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(3)) = zero

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_) = p0

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(1)) = &
       -b0 * sin(2.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(2)) = &
        b0 * sin(4.0_dp * dpi * x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1))
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(3)) = zero

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  ! Initialise the vectorpotential on the corners - used with staggered grid
  subroutine initvecpot_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, xC, A, idir)
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idir
    double precision, intent(in)       :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)      :: A(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if (idir == 3) then
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = &
         b0 * half / dpi * (half * cos(4.0_dp * dpi * &
         xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)) + &
         cos(2.0_dp * dpi * xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)))
    else
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = zero
    end if

  end subroutine initvecpot_usr


  subroutine specialvar_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
    use mod_functions_bfield
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw+nwauxio)
    double precision             :: normconv(0:nw+nwauxio)

    double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1) = &
       divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine specialvar_output


  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames = 'divb'
  end subroutine specialvarnames_output

end module mod_usr
