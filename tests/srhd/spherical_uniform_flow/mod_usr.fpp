!> Uniform Cartesian relativistic flow on a spherical mesh.
!>
!> A constant density, constant pressure gas moving with a constant
!> (sub-luminal) Cartesian velocity is an exact steady solution of the SRHD
!> equations, exactly as it is for the classical Euler equations. Written in
!> spherical components the velocity is a non-trivial function of theta and
!> phi, so keeping this state exactly uniform exercises every curvilinear
!> flux and every geometric source term. Any error in them shows up
!> immediately as drift in the volume averages the log file reports.
!>
!> The Lorentz factor only depends on |v|, which a spherical rotation of a
!> Cartesian vector leaves invariant, so it is a single global constant here
!> rather than a per-cell quantity.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> ambient state. These are deliberately plain variables rather than
  !> parameters: specialbound_usr runs on the device, and a parameter has no
  !> storage to copy there.
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
  !> the uniform velocity, in Cartesian components (must satisfy |v0| < 1)
  double precision :: v0(3) = [0.3d0, 0.15d0, -0.1d0]
  !$acc declare copyin(rho0, p0, v0)

contains

  subroutine usr_init()

    ! set_coordinate_system is no longer required: the coordinate system comes
    ! from `geometry` in &meshlist. Calling it explicitly is still valid.

    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr

    call phys_activate()

  end subroutine usr_init

  !> spherical components of the uniform Cartesian *ordinary* velocity at x.
  !> A subroutine rather than an array-valued function, because a function
  !> result needs a temporary that not every OpenACC compiler handles inside
  !> a device routine.
  pure subroutine uniform_velocity(x, v)
    !$acc routine seq
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: v(1:3)
    double precision              :: sint, cost, sinp, cosp

    sint = sin(x(2)); cost = cos(x(2))
    sinp = sin(x(3)); cosp = cos(x(3))

    v(1) =  v0(1)*sint*cosp + v0(2)*sint*sinp + v0(3)*cost
    v(2) =  v0(1)*cost*cosp + v0(2)*cost*sinp - v0(3)*sint
    v(3) = -v0(1)*sinp      + v0(2)*cosp

  end subroutine uniform_velocity

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    ! .. local ..
    double precision                :: v(1:3), lfac0
    integer                         :: ix1, ix2, ix3

    lfac0 = 1.0d0 / sqrt(1.0d0 - sum(v0(1:3)**2))

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             call uniform_velocity(x(ix1,ix2,ix3,1:ndim), v)
             w(ix1,ix2,ix3,rho_)   = rho0
             w(ix1,ix2,ix3,p_)     = p0
             ! primitive mom(:) is the spatial four-velocity lfac*v, not v
             w(ix1,ix2,ix3,mom(1)) = lfac0 * v(1)
             w(ix1,ix2,ix3,mom(2)) = lfac0 * v(2)
             w(ix1,ix2,ix3,mom(3)) = lfac0 * v(3)
          end do
       end do
    end do

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  !> keep the exact solution in the ghost cells, so that only the interior
  !> discretisation is under test
  subroutine specialbound_usr(qt, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB, w, x)
    !$acc routine vector
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3
    integer, intent(in)             :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    integer, intent(in)             :: iB
    double precision, intent(in)    :: qt
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    ! .. local ..
    double precision                :: v(1:3), lfac0, xi0
    integer                         :: ix1, ix2, ix3

    lfac0 = 1.0d0 / sqrt(1.0d0 - sum(v0(1:3)**2))
    ! xi = rho*h*lfac^2, gamma-law enthalpy h = 1 + gamma/(gamma-1) * p0/rho0
    xi0 = (rho0 + gamma_to_gamma_1 * p0) * lfac0**2

    !$acc loop collapse(3) vector private(v)
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

             call uniform_velocity(x(ix1,ix2,ix3,1:ndim), v)

             ! conserved D = rho*lfac, S_i = xi*v_i, tau = xi - p - D
             w(ix1,ix2,ix3,iw_rho)    = rho0 * lfac0
             w(ix1,ix2,ix3,iw_mom(1)) = xi0 * v(1)
             w(ix1,ix2,ix3,iw_mom(2)) = xi0 * v(2)
             w(ix1,ix2,ix3,iw_mom(3)) = xi0 * v(3)
             w(ix1,ix2,ix3,iw_e)      = xi0 - p0 - rho0 * lfac0

          end do
       end do
    end do

  end subroutine specialbound_usr

end module mod_usr
