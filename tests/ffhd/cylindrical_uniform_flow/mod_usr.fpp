!> Uniform flow along a uniform frozen field, on a cylindrical mesh.
!>
!> ffhd only evolves scalars (rho, the field-aligned momentum m_par, the
!> field-aligned energy e_hd_par) advected along a fixed field direction
!> b-hat that the user supplies per cell. Take b-hat to be a constant
!> Cartesian direction, and rho, the parallel speed v_par and p spatially
!> uniform: since div(b-hat) = 0 for a uniform Cartesian vector field and
!> every advected quantity is itself spatially constant, every flux
!> divergence in the module's PDEs (see doc/equations.md's "Frozen-field
!> hydrodynamics" section) vanishes identically, in any coordinate system.
!> This is therefore an exact steady solution, exactly as the analogous
!> uniform-flow states are for tests/hd,mhd,srhd/cylindrical_uniform_flow --
!> except here there is no addsource_geometry to exercise (ffhd's conserved
!> quantities are genuine scalars, whose flux divergence by Gauss's theorem
!> needs no curvature correction; see mod_ffhd_templates.fpp). What this does
!> exercise is that b-hat's cylindrical components, which vary with phi even
!> though the physical direction is constant, are handled correctly by the
!> ordinary curvilinear flux machinery.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> ambient state. These are deliberately plain variables rather than
  !> parameters: specialbound_usr runs on the device, and a parameter has no
  !> storage to copy there.
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
  double precision :: vpar0 = 0.4d0
  !> the frozen field direction, in Cartesian components (need not be unit
  !> length: only its direction matters, normalised in to_cylindrical_unit)
  double precision :: b0(3) = [1.0d0, 0.5d0, -0.3d0]
  !$acc declare copyin(rho0, p0, vpar0, b0)

contains

  subroutine usr_init()

    ! set_coordinate_system is no longer required: the coordinate system comes
    ! from `geometry` in &meshlist. Calling it explicitly is still valid.

    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr

    call phys_activate()

  end subroutine usr_init

  !> cylindrical (r, z, phi) components at x of the unit vector along
  !> Cartesian vec0. A subroutine rather than an array-valued function,
  !> because a function result needs a temporary that not every OpenACC
  !> compiler handles inside a device routine.
  pure subroutine to_cylindrical_unit(x, vec0, vec)
    !$acc routine seq
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(in)  :: vec0(1:3)
    double precision, intent(out) :: vec(1:3)
    double precision              :: sinp, cosp, inv_norm
    double precision              :: cart(1:3)

    sinp = sin(x(3)); cosp = cos(x(3))

    inv_norm = 1.0d0 / sqrt(sum(vec0(1:3)**2))
    cart = vec0 * inv_norm

    vec(1) =  cart(1)*cosp + cart(2)*sinp
    vec(2) =  cart(3)
    vec(3) = -cart(1)*sinp + cart(2)*cosp

  end subroutine to_cylindrical_unit

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    ! .. local ..
    double precision                :: bhat(1:3)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call to_cylindrical_unit(x_loc, b0, bhat)
             w(ix1,ix2,ix3,rho_)   = rho0
             w(ix1,ix2,ix3,p_)     = p0
             w(ix1,ix2,ix3,mom(1)) = vpar0
             w(ix1,ix2,ix3,iw_b1)  = bhat(1)
             w(ix1,ix2,ix3,iw_b2)  = bhat(2)
             w(ix1,ix2,ix3,iw_b3)  = bhat(3)
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
    double precision                :: bhat(1:3)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    !$acc loop collapse(3) vector private(bhat, x_loc)
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call to_cylindrical_unit(x_loc, b0, bhat)

             w(ix1,ix2,ix3,iw_rho)    = rho0
             w(ix1,ix2,ix3,iw_mom(1)) = rho0 * vpar0
             w(ix1,ix2,ix3,iw_e)      = p0 * inv_gamma_1 + &
                0.5d0 * rho0 * vpar0**2
             w(ix1,ix2,ix3,iw_b1)     = bhat(1)
             w(ix1,ix2,ix3,iw_b2)     = bhat(2)
             w(ix1,ix2,ix3,iw_b3)     = bhat(3)

          end do
       end do
    end do

  end subroutine specialbound_usr

end module mod_usr
