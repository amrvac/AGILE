!> Uniform field-aligned FFHD flow on a cylindrical mesh onto the axis at r=0.
!>
!> The cylindrical counterpart of tests/ffhd/spherical_pole. A constant
!> density, constant pressure gas moving at a constant speed vpar0 along a
!> uniform Cartesian frozen field b0 is an exact steady solution: rho, the
!> field-aligned momentum m_par and the energy are genuine scalars and stay
!> literally constant, and a uniform Cartesian b-hat has zero divergence.
!> Written in cylindrical (r, z, phi) components b-hat is a non-trivial
!> function of phi, so keeping this state uniform exercises the curvilinear
!> fluxes.
!>
!> This is a pole case: the domain runs onto the singular axis at r=0. FFHD
!> has no boundary condition for the frozen field and does not exchange it
!> through getbc; fill_frozen_field_device re-derives it from
!> usr_set_frozen_field in every cell of every block, the axis ghost cells
!> included. bgeo%x in the ghost layer beyond r=0 carries a negative r, and
!> the analytic fill there reproduces on its own the sign flips a vector picks
!> up across the axis - b_r and b_phi antisymmetric, b_z symmetric - with no
!> entry in the pole sign table.
!>
!> Note the radial entry: upstream MPI-AMRVAC marks only b_phi antisymmetric at
!> the cylindrical axis, but the map (r, phi) -> (r, phi+pi) flips e_r exactly
!> as it flips e_phi, so b_r is odd too. This case pins that down (see "The
!> polar axis" in CLAUDE.md).
!>
!> The fluid variables, being scalars, are simply symmetric across the axis.
!> The volume averages in the log are almost blind to the pole copy, so the
!> real test is check_pole_ghosts below.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision :: rho0  = 1.0d0
  double precision :: p0    = 1.0d0
  double precision :: vpar0 = 0.4d0
  double precision :: b0(3) = [1.0d0, 0.5d0, -0.3d0]
  !$acc declare copyin(rho0, p0, vpar0, b0)

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_process_grid  => check_pole_ghosts

    call phys_activate()

  end subroutine usr_init

  !> cylindrical (r, z, phi) components at x of the unit vector along the
  !> Cartesian vector vec0.
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

  !> The frozen field, in cylindrical components at x. Called by name from
  !> fill_frozen_field_device (device code); see mod_usr_methods.
  pure subroutine usr_set_frozen_field(x, bhat)
    !$acc routine seq
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: bhat(1:3)

    call to_cylindrical_unit(x, b0, bhat)

  end subroutine usr_set_frozen_field

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    integer                         :: ix1, ix2, ix3

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             w(ix1,ix2,ix3,rho_)   = rho0
             w(ix1,ix2,ix3,p_)     = p0
             w(ix1,ix2,ix3,mom(1)) = vpar0
          end do
       end do
    end do

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  !> Keep the exact solution in the ghost cells at the boundaries that are not
  !> the axis, so that only the interior discretisation is under test.
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
    integer                         :: ix1, ix2, ix3

    !$acc loop collapse(3) vector
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
             w(ix1,ix2,ix3,iw_rho)    = rho0
             w(ix1,ix2,ix3,iw_mom(1)) = rho0 * vpar0
             w(ix1,ix2,ix3,iw_e)      = p0 * inv_gamma_1 + &
                0.5d0 * rho0 * vpar0**2
          end do
       end do
    end do

  end subroutine specialbound_usr

  !> Assert that the ghost cells across the r=0 axis carry the exact solution.
  !>
  !> Only meaningful at it = 0. The fluid variables (constants) come from
  !> getbc's pole copy; the frozen field comes from fill_frozen_field_device at
  !> the ghost cell's own negative-r coordinate, and rebuilt in Cartesian it
  !> has to equal the uniform b0. Silence means the check passed.
  subroutine check_pole_ghosts(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    ! .. local ..
    double precision :: err, efluid, sinp, cosp, invb, bcart(3), e0
    integer          :: ix1, ix2, ix3

    if (it /= 0) return
    ! the cylindrical axis is the r = 0 face only
    if (neighbor_pole(-1,0,0,igrid) == 0) return
    if (neighbor_type(-1,0,0,igrid) /= neighbor_sibling) return
    !$acc update host(ps(igrid)%w)

    e0   = p0 * inv_gamma_1 + 0.5d0 * rho0 * vpar0**2
    invb = 1.0d0 / sqrt(sum(b0(1:3)**2))
    err  = 0.0d0

    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixImin1, ixOmin1-1
             efluid = max(abs(w(ix1,ix2,ix3,rho_) - rho0), &
                  abs(w(ix1,ix2,ix3,iw_mom(1)) - rho0*vpar0), &
                  abs(w(ix1,ix2,ix3,iw_e) - e0))
             sinp = sin(x(ix1,ix2,ix3,3)); cosp = cos(x(ix1,ix2,ix3,3))
             bcart(1) = w(ix1,ix2,ix3,iw_b1)*cosp - w(ix1,ix2,ix3,iw_b3)*sinp
             bcart(2) = w(ix1,ix2,ix3,iw_b1)*sinp + w(ix1,ix2,ix3,iw_b3)*cosp
             bcart(3) = w(ix1,ix2,ix3,iw_b2)
             err = max(err, efluid, maxval(abs(bcart - b0*invb)))
          end do
       end do
    end do

    if (err > 1.0d-10) then
       write(*,*) 'pole ghost cells deviate from the exact solution by', err
       call mpistop('pole ghost-cell check failed')
    end if

  end subroutine check_pole_ghosts

end module mod_usr
