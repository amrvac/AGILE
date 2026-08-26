!> Uniform Cartesian MHD flow on a spherical mesh.
!>
!> A constant density, constant pressure gas moving with a constant Cartesian
!> velocity and threaded by a constant Cartesian magnetic field is an exact
!> steady solution of the MHD equations. Written in spherical components,
!> both the velocity and field are non-trivial functions of theta and phi, so
!> keeping this state exactly uniform exercises every curvilinear flux and
!> every geometric source term in addsource_geometry, including the induction
!> equation and GLM psi coupling terms. Any error in them shows up
!> immediately as drift in the volume averages the log file reports.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> ambient state. These are deliberately plain variables rather than
  !> parameters: specialbound_usr runs on the device, and a parameter has no
  !> storage to copy there.
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
  !> the uniform velocity, in Cartesian components
  double precision :: v0(3) = [1.0d0, 0.5d0, -0.3d0]
  !> the uniform magnetic field, in Cartesian components
  double precision :: b0(3) = [0.2d0, -0.4d0, 0.1d0]
  !$acc declare copyin(rho0, p0, v0, b0)

contains

  subroutine usr_init()

    ! set_coordinate_system is no longer required: the coordinate system comes
    ! from `geometry` in &meshlist. Calling it explicitly is still valid.

    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr

    call phys_activate()

  end subroutine usr_init

  !> spherical components at x of a Cartesian vector v0. Used for both the
  !> velocity and the magnetic field, since this is just a change of basis.
  !> A subroutine rather than an array-valued function, because a function
  !> result needs a temporary that not every OpenACC compiler handles inside
  !> a device routine.
  pure subroutine to_spherical_vector(x, vec0, vec)
    !$acc routine seq
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(in)  :: vec0(1:3)
    double precision, intent(out) :: vec(1:3)
    double precision              :: sint, cost, sinp, cosp

    sint = sin(x(2)); cost = cos(x(2))
    sinp = sin(x(3)); cosp = cos(x(3))

    vec(1) =  vec0(1)*sint*cosp + vec0(2)*sint*sinp + vec0(3)*cost
    vec(2) =  vec0(1)*cost*cosp + vec0(2)*cost*sinp - vec0(3)*sint
    vec(3) = -vec0(1)*sinp      + vec0(2)*cosp

  end subroutine to_spherical_vector

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    ! .. local ..
    double precision                :: v(1:3), b(1:3)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call to_spherical_vector(x_loc, v0, v)
             call to_spherical_vector(x_loc, b0, b)
             w(ix1,ix2,ix3,rho_)   = rho0
             w(ix1,ix2,ix3,p_)     = p0
             w(ix1,ix2,ix3,mom(1)) = v(1)
             w(ix1,ix2,ix3,mom(2)) = v(2)
             w(ix1,ix2,ix3,mom(3)) = v(3)
             w(ix1,ix2,ix3,mag(1)) = b(1)
             w(ix1,ix2,ix3,mag(2)) = b(2)
             w(ix1,ix2,ix3,mag(3)) = b(3)
             w(ix1,ix2,ix3,psi_)   = 0.0d0
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
    double precision                :: v(1:3), b(1:3)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    !$acc loop collapse(3) vector private(v, b, x_loc)
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call to_spherical_vector(x_loc, v0, v)
             call to_spherical_vector(x_loc, b0, b)

             w(ix1,ix2,ix3,iw_rho)    = rho0
             w(ix1,ix2,ix3,iw_mom(1)) = rho0 * v(1)
             w(ix1,ix2,ix3,iw_mom(2)) = rho0 * v(2)
             w(ix1,ix2,ix3,iw_mom(3)) = rho0 * v(3)
             w(ix1,ix2,ix3,iw_mag(1)) = b(1)
             w(ix1,ix2,ix3,iw_mag(2)) = b(2)
             w(ix1,ix2,ix3,iw_mag(3)) = b(3)
             w(ix1,ix2,ix3,psi_)      = 0.0d0
             w(ix1,ix2,ix3,iw_e)      = p0 * inv_gamma_1 + &
                0.5d0 * rho0 * sum(v(1:3)**2) + 0.5d0 * sum(b(1:3)**2)

          end do
       end do
    end do

  end subroutine specialbound_usr

end module mod_usr
