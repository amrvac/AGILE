!> Uniform field-aligned FFHD flow on a spherical mesh through both poles.
!>
!> A constant density, constant pressure gas moving at a constant speed vpar0
!> along a uniform Cartesian frozen field b0 is an exact steady solution: the
!> conserved quantities rho, the field-aligned momentum m_par and the energy
!> are all genuine scalars, so they stay literally constant, and a uniform
!> Cartesian b-hat has zero divergence so the scalar fluxes balance. Written in
!> spherical (r, theta, phi) components b-hat is a non-trivial function of
!> theta and phi, so keeping this state uniform exercises every curvilinear
!> flux.
!>
!> This is a pole case: the domain deliberately runs onto the singular axis at
!> both ends. FFHD has no boundary condition for the frozen field and does not
!> exchange it through getbc; fill_nwextra_device re-derives it from
!> usr_set_nwextra in every cell of every block after each grid change,
!> the polar-axis ghost cells included. Because bgeo%x in the ghost layer
!> beyond theta=0 carries the mirrored coordinate (a negative theta), that
!> analytic fill reproduces on its own the sign flips a vector picks up across
!> the axis - b_r symmetric, b_theta and b_phi antisymmetric - with no entry
!> in the pole sign table.
!>
!> The fluid variables, being scalars, are simply symmetric across the axis,
!> which is what src/io/mod_input_output.fpp's read_par_files relies on when it
!> skips the momentum lines of the pole sign table for phys='ffhd'.
!>
!> The volume averages in the log are almost blind to whether the pole copy is
!> right - a cell at the axis has vanishing volume, the face on the axis zero
!> area, and the TVD limiter clips what is left - so the real test is
!> check_pole_ghosts below: it compares the fluid ghost cells against the
!> analytic constant (the getbc pole copy) and the frozen-field ghost cells
!> against the uniform Cartesian b0 they must reduce to (fill_nwextra_device
!> at the axis).
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> ambient state. Plain variables, not parameters: specialbound_usr and
  !> usr_set_nwextra run on the device, and a parameter has no storage to
  !> copy there.
  double precision :: rho0  = 1.0d0
  double precision :: p0    = 1.0d0
  !> field-aligned speed, and the uniform Cartesian frozen-field direction
  !> (any length: fill_nwextra_device normalises it)
  double precision :: vpar0 = 0.4d0
  double precision :: b0(3) = [1.0d0, 0.5d0, -0.3d0]
  !$acc declare copyin(rho0, p0, vpar0, b0)

contains

  subroutine usr_init()

    ! the coordinate system comes from `geometry` in &meshlist

    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_process_grid  => check_pole_ghosts

    call phys_activate()

  end subroutine usr_init

  !> spherical (r, theta, phi) components at x of the unit vector along the
  !> Cartesian vector vec0.
  pure subroutine to_spherical_unit(x, vec0, vec)
    !$acc routine seq
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(in)  :: vec0(1:3)
    double precision, intent(out) :: vec(1:3)
    double precision              :: sint, cost, sinp, cosp, inv_norm
    double precision              :: cart(1:3)

    sint = sin(x(2)); cost = cos(x(2))
    sinp = sin(x(3)); cosp = cos(x(3))

    inv_norm = 1.0d0 / sqrt(sum(vec0(1:3)**2))
    cart = vec0 * inv_norm

    vec(1) =  cart(1)*sint*cosp + cart(2)*sint*sinp + cart(3)*cost
    vec(2) =  cart(1)*cost*cosp + cart(2)*cost*sinp - cart(3)*sint
    vec(3) = -cart(1)*sinp      + cart(2)*cosp

  end subroutine to_spherical_unit

  !> The frozen field, in spherical components at x. Called by name from
  !> fill_nwextra_device (device code); see mod_usr_methods.
  pure subroutine usr_set_nwextra(x, bhat)
    !$acc routine seq
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: bhat(1:3)

    call to_spherical_unit(x, b0, bhat)

  end subroutine usr_set_nwextra

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    integer                         :: ix1, ix2, ix3

    ! rho, the field-aligned momentum and the energy are constants; b1,b2,b3
    ! are filled by fill_nwextra_device from usr_set_nwextra
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

  !> Keep the exact solution in the ghost cells at the radial boundaries, so
  !> that only the interior discretisation is under test.
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

    ! b1,b2,b3 in these ghost cells come from fill_nwextra_device
    ! Vector on the innermost loop, deliberately not collapse(3): nvfortran's
    ! OpenACC miscompiles a collapsed vector loop that calls another !$acc
    ! routine in its body, reached through this !$acc routine vector from the
    ! gang loops in fill_boundary_before_gc / fill_boundary_after_gc. It bit
    ! the hd curvilinear pole cases via analytic_state; the other physics
    ! pole cases keep the same form. See CLAUDE.md - very likely the same
    ! defect as issue #154.
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          !$acc loop vector
          do ix1 = ixOmin1, ixOmax1
             w(ix1,ix2,ix3,iw_rho)    = rho0
             w(ix1,ix2,ix3,iw_mom(1)) = rho0 * vpar0
             w(ix1,ix2,ix3,iw_e)      = p0 * inv_gamma_1 + &
                0.5d0 * rho0 * vpar0**2
          end do
       end do
    end do

  end subroutine specialbound_usr

  !> Assert that the ghost cells across the axis carry the exact solution.
  !>
  !> Only meaningful at it = 0, where the interior is still exactly analytic.
  !> Two things are checked at a same-level pole neighbour:
  !>
  !>  * the fluid variables (rho, m_par, e), which getbc's pole copy fills from
  !>    the block half a revolution away - they are constants, so a correct
  !>    copy reproduces them to round-off;
  !>
  !>  * the frozen field, which fill_nwextra_device fills from
  !>    usr_set_nwextra at the ghost cell's own (mirrored) coordinate -
  !>    converted back to Cartesian it has to equal the uniform b0, which is
  !>    only true if the sign flips across the axis have come out right.
  !>
  !> Silence means the check passed.
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
    double precision :: err, efluid, sint, cost, sinp, cosp, invb, bcart(3)
    double precision :: e0
    integer          :: ix1, ix2, ix3, iside, i2, jxmin2, jxmax2

    if (it /= 0) return
    !$acc update host(ps(igrid)%w)

    e0   = p0 * inv_gamma_1 + 0.5d0 * rho0 * vpar0**2
    invb = 1.0d0 / sqrt(sum(b0(1:3)**2))
    err  = 0.0d0

    do iside = 1, 2
       i2 = 2*iside - 3
       if (neighbor_pole(0,i2,0,igrid) == 0) cycle
       if (neighbor_type(0,i2,0,igrid) /= neighbor_sibling) cycle
       if (iside == 1) then
          jxmin2 = ixImin2;   jxmax2 = ixOmin2-1
       else
          jxmin2 = ixOmax2+1; jxmax2 = ixImax2
       end if
       ! transverse loops run over the whole block (ixI), not just its
       ! interior (ixO), so the pole layer's edges and corners are covered -
       ! in particular the cells where the axis meets the radial boundary
       do ix3 = ixImin3, ixImax3
          do ix2 = jxmin2, jxmax2
             do ix1 = ixImin1, ixImax1
                ! fluid variables: constant everywhere
                efluid = max(abs(w(ix1,ix2,ix3,rho_) - rho0), &
                     abs(w(ix1,ix2,ix3,iw_mom(1)) - rho0*vpar0), &
                     abs(w(ix1,ix2,ix3,iw_e) - e0))
                ! frozen field: rebuild the Cartesian vector from the ghost
                ! cell's spherical components and its own coordinates
                sint = sin(x(ix1,ix2,ix3,2)); cost = cos(x(ix1,ix2,ix3,2))
                sinp = sin(x(ix1,ix2,ix3,3)); cosp = cos(x(ix1,ix2,ix3,3))
                bcart(1) = w(ix1,ix2,ix3,iw_b1)*sint*cosp &
                         + w(ix1,ix2,ix3,iw_b2)*cost*cosp &
                         - w(ix1,ix2,ix3,iw_b3)*sinp
                bcart(2) = w(ix1,ix2,ix3,iw_b1)*sint*sinp &
                         + w(ix1,ix2,ix3,iw_b2)*cost*sinp &
                         + w(ix1,ix2,ix3,iw_b3)*cosp
                bcart(3) = w(ix1,ix2,ix3,iw_b1)*cost &
                         - w(ix1,ix2,ix3,iw_b2)*sint
                err = max(err, efluid, maxval(abs(bcart - b0*invb)))
             end do
          end do
       end do
    end do

    if (err > 1.0d-10) then
       write(*,*) 'pole ghost cells deviate from the exact solution by', err
       call mpistop('pole ghost-cell check failed')
    end if

  end subroutine check_pole_ghosts

end module mod_usr
