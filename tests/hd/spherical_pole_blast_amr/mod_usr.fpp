!> A blast wave that starts off the polar axis and expands across it, on an
!> AMR mesh.
!>
!> This is the demanding pole case. The uniform-flow pole tests keep a state
!> that is smooth across the axis, which is exactly the situation in which a
!> wrong ghost value is harmless: the axis face has zero area, so the flux it
!> would carry is multiplied by zero, and the jump left behind is clipped by
!> the TVD limiter. Here a strong shock sweeps over the pole instead. The
!> state on either side of the axis is genuinely different, nothing is clipped,
!> and the refinement follows the shell straight across - so the same-level,
!> restricting and prolonging pole paths all run, repeatedly, on regridded
!> blocks rather than a static mesh.
!>
!> The hot spot starts a little off the axis and is carried onto it by a
!> uniform Cartesian background velocity pointing mostly along +z. Both
!> matter: the offset means the initial condition itself never touches the
!> axis, so t = 0 is unambiguous, and the background velocity gives every
!> momentum component a non-trivial value there, which is what lets
!> check_pole_ghosts below test the sign flips and not just the mirror.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> ambient state
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
  !> overpressure and radius of the initial hot spot, in length units
  double precision :: pblast = 1.0d1
  double precision :: rblast = 0.25d0
  !> uniform Cartesian background velocity. Mostly along +z, so the blast is
  !> advected onto the north pole rather than only expanding into it, with
  !> transverse components so that nothing about the setup is symmetric.
  double precision :: v0(3) = [0.1d0, -0.05d0, 0.35d0]
  !> Centre of the hot spot in spherical coordinates. theta_c is small but
  !> non-zero: the spot sits near the north pole without straddling it, and
  !> check_blast_fits verifies that.
  double precision :: rc     = 1.0d0
  double precision :: thetac = 0.10d0 * dpi
  double precision :: phic   = 0.50d0 * dpi
  !$acc declare copyin(rho0, p0, pblast, rblast, v0, rc, thetac, phic)

contains

  subroutine usr_init()

    ! the coordinate system comes from `geometry` in &meshlist

    usr_init_one_grid => initonegrid_usr
    usr_process_grid  => check_pole_ghosts

    call phys_activate()

  end subroutine usr_init

  !> Analytic initial state at x, in conserved variables.
  !>
  !> Valid at any coordinates, the ghost cells beyond the axis included, whose
  !> stored theta is negative. That is not a coincidence: the Cartesian point
  !> built from (r, -theta, phi) is the very point at (r, theta, phi + pi) that
  !> the pole copy fetches from, and writing a uniform Cartesian vector in the
  !> ghost cell's own basis reproduces the sign flips the copy applies. So this
  !> one routine is both the initial condition and the reference the pole
  !> ghost cells are checked against.
  pure subroutine analytic_state(x_loc, wpt)
    double precision, intent(in)  :: x_loc(1:ndim)
    double precision, intent(out) :: wpt(1:nw)
    ! .. local ..
    double precision :: sint, cost, sinp, cosp
    double precision :: xcart(3), xc(3), v(3), pres

    sint = sin(x_loc(2)); cost = cos(x_loc(2))
    sinp = sin(x_loc(3)); cosp = cos(x_loc(3))

    ! the blast is spherical in physical space, so place it by Cartesian
    ! distance
    xcart(1) = x_loc(1) * sint * cosp
    xcart(2) = x_loc(1) * sint * sinp
    xcart(3) = x_loc(1) * cost
    xc(1)    = rc * sin(thetac) * cos(phic)
    xc(2)    = rc * sin(thetac) * sin(phic)
    xc(3)    = rc * cos(thetac)

    if (sum((xcart - xc)**2) > rblast**2) then
       pres = p0
    else
       pres = pblast
    end if

    ! spherical components of the uniform Cartesian background velocity
    v(1) =  v0(1)*sint*cosp + v0(2)*sint*sinp + v0(3)*cost
    v(2) =  v0(1)*cost*cosp + v0(2)*cost*sinp - v0(3)*sint
    v(3) = -v0(1)*sinp      + v0(2)*cosp

    wpt(1:nw)      = 0.0d0
    wpt(iw_rho)    = rho0
    wpt(iw_mom(1)) = rho0 * v(1)
    wpt(iw_mom(2)) = rho0 * v(2)
    wpt(iw_mom(3)) = rho0 * v(3)
    wpt(iw_e)      = pres * inv_gamma_1 + 0.5d0 * rho0 * sum(v(1:3)**2)

  end subroutine analytic_state

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    ! .. local ..
    double precision                :: wpt(1:nw), x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    call check_blast_fits

    ! analytic_state returns conserved variables directly, so unlike the other
    ! blast cases there is no phys_to_conserved call here
    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call analytic_state(x_loc, wpt)
             w(ix1,ix2,ix3,1:nw) = wpt(1:nw)
          end do
       end do
    end do

  end subroutine initonegrid_usr

  !> Abort unless the hot spot starts wholly inside the radial domain and
  !> wholly off the axis.
  !>
  !> Reaching the axis is no longer the error it is in spherical_blast_amr -
  !> the whole point here is that the shell crosses it - but the *initial*
  !> spot has to stay clear of it, so that the state at t = 0 is smooth across
  !> the axis and check_pole_ghosts can compare against an unambiguous
  !> analytic value. The shell is free to cross it a few steps later.
  subroutine check_blast_fits
    use mod_comm_lib, only: mpistop
    ! .. local ..
    double precision :: axis_distance

    if (rc - rblast < xprobmin1 .or. rc + rblast > xprobmax1) &
       call mpistop("spherical_pole_blast: hot spot sticks out of the radial &
          &domain")

    ! perpendicular distance from the centre of the spot to the polar axis
    axis_distance = rc * sin(thetac)
    if (axis_distance <= rblast) &
       call mpistop("spherical_pole_blast: hot spot already straddles the &
          &polar axis at t = 0; move it out in theta")

  end subroutine check_blast_fits

  !> Assert that the ghost cells across the axis carry the exact solution.
  !>
  !> Only meaningful at it = 0, where the interior is still exactly the
  !> analytic state, so the pole copy of it has to reproduce that state in the
  !> ghost cells to round-off. Later on there is no closed form to compare
  !> against, and the log takes over as the regression: unlike the smooth
  !> uniform-flow cases, a shock crossing the axis leaves the volume averages
  !> genuinely sensitive to the pole copy. Silence means the check passed.
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
    double precision :: wpt(1:nw), x_loc(1:ndim), err, tol
    integer          :: ix1, ix2, ix3, iside, i2, jxmin2, jxmax2

    if (it /= 0) return
    ! process() runs before the solution is pulled back for output, so the
    ! host copy of w is stale unless we fetch this block ourselves
    !$acc update host(ps(igrid)%w)

    do iside = 1, 2
       i2 = 2*iside - 3
       if (neighbor_pole(0,i2,0,igrid) == 0) cycle
       if (iside == 1) then
          jxmin2 = ixImin2;   jxmax2 = ixOmin2-1
       else
          jxmin2 = ixOmax2+1; jxmax2 = ixImax2
       end if
       ! A pole neighbour at the same level is copied verbatim, so its ghost
       ! cells have to reproduce the analytic state to round-off. Across a
       ! level jump - and the mesh is already refined at it = 0 - the values
       ! are restricted or prolonged on the way and carry the scheme's own
       ! second-order error. The loose tolerance there still catches what
       ! actually goes wrong in that path, a mistaken destination offset or a
       ! missing sign flip, both of which are O(1) on this initial condition.
       if (neighbor_type(0,i2,0,igrid) == neighbor_sibling) then
          tol = 1.0d-10
       else
          tol = 1.0d-1
       end if
       err = 0.0d0
       do ix3 = ixOmin3, ixOmax3
          do ix2 = jxmin2, jxmax2
             do ix1 = ixOmin1, ixOmax1
                x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
                call analytic_state(x_loc, wpt)
                err = max(err, maxval(abs(w(ix1,ix2,ix3,&
                   1:nwflux) - wpt(1:nwflux))))
             end do
          end do
       end do
       if (err > tol) then
          write(*,*) 'pole ghost cells deviate from the exact solution by',&
             err, ' tolerance', tol, ' neighbour type',&
             neighbor_type(0,i2,0,igrid)
          call mpistop('pole ghost-cell check failed')
       end if
    end do

  end subroutine check_pole_ghosts

end module mod_usr
