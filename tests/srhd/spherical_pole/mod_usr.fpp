!> Uniform, sub-luminal Cartesian SRHD flow on a spherical mesh through both
!> poles.
!>
!> A constant density, constant pressure gas moving with a constant Cartesian
!> velocity is an exact steady solution: the Lorentz factor depends only on
!> |v|, which a spherical rotation leaves invariant, so it stays a single
!> global constant however the velocity's components are written. Written in
!> spherical components the velocity itself is a non-trivial function of theta
!> and phi, so keeping this state uniform exercises every curvilinear flux and
!> every geometric source term.
!>
!> This is a pole case: the domain deliberately runs onto the singular axis,
!> where the ghost cells come from the block half a revolution away in phi
!> rather than from a boundary condition. Note that the volume averages in the
!> log are almost blind to whether that copy is right - a cell at the axis has
!> vanishing volume, the face on the axis has zero area, and the TVD limiter
!> clips what is left - so the real test here is check_pole_ghosts below, which
!> compares the ghost cells themselves against the exact solution.
!>
!> SRHD's primitive mom(:) slot holds the spatial four-velocity u^i = lfac*v^i
!> rather than v^i itself (see mod_srhd_templates.fpp's to_primitive /
!> to_conservative). lfac is a scalar, so it is invariant under the pole's
!> pi-rotation, and u^i therefore transforms exactly like an ordinary vector -
!> the same asymm/symm pattern hd's momentum takes at the pole applies here
!> unchanged, which is what src/io/mod_input_output.fpp's read_par_files relies
!> on when it builds the sign table from iw_mom alone.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> ambient state. These are deliberately plain variables rather than
  !> parameters: specialbound_usr runs on the device, and a parameter has no
  !> storage to copy there.
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
  !> the uniform *ordinary* Cartesian velocity, sub-luminal
  double precision :: v0(3) = [0.3d0, 0.15d0, -0.1d0]
  !$acc declare copyin(rho0, p0, v0)

contains

  subroutine usr_init()

    ! the coordinate system comes from `geometry` in &meshlist

    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_process_grid  => check_pole_ghosts

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
    double precision                :: v(1:3), lfac0, x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    lfac0 = 1.0d0 / sqrt(1.0d0 - sum(v0(1:3)**2))

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call uniform_velocity(x_loc, v)
             w(ix1,ix2,ix3,rho_)   = rho0
             w(ix1,ix2,ix3,p_)     = p0
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
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    lfac0 = 1.0d0 / sqrt(1.0d0 - sum(v0(1:3)**2))
    ! xi = rho*h*lfac^2, gamma-law enthalpy h = 1 + gamma/(gamma-1) * p0/rho0
    xi0 = (rho0 + gamma_to_gamma_1 * p0) * lfac0**2

    ! collapse(3) is fine here even though the body calls uniform_velocity:
    ! the nvfortran OpenACC miscompile that hit the hd pole cases needs BOTH a
    ! call AND a runtime-sized private automatic array in the collapsed body
    ! (hd's wpt(1:nw)). srhd's privates (v, x_loc) are all compile-time sized,
    ! so there is nothing to miscompile. See CLAUDE.md ("Bug-hunting notes")
    ! and issue #154.
    !$acc loop collapse(3) vector private(v, x_loc)
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call uniform_velocity(x_loc, v)

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

  !> Analytic conserved state at x, valid at any coordinates - the ghost cells
  !> beyond the axis included, whose stored theta is negative. That those
  !> still come out right is not a coincidence: the mirror and the
  !> per-component sign the pole copy applies are exactly the change of basis
  !> between a point and its partner half a revolution away, so writing the
  !> uniform Cartesian velocity in the ghost cell's own coordinates reproduces
  !> what the copy has to have put there - lfac0 itself does not change, since
  !> it depends only on |v0|.
  pure subroutine analytic_state(x_loc, wpt)
    double precision, intent(in)  :: x_loc(1:ndim)
    double precision, intent(out) :: wpt(1:nw)
    double precision              :: v(1:3), lfac0, xi0

    call uniform_velocity(x_loc, v)
    lfac0 = 1.0d0 / sqrt(1.0d0 - sum(v0(1:3)**2))
    xi0 = (rho0 + gamma_to_gamma_1 * p0) * lfac0**2

    wpt(1:nw)      = 0.0d0
    wpt(iw_rho)    = rho0 * lfac0
    wpt(iw_mom(1)) = xi0 * v(1)
    wpt(iw_mom(2)) = xi0 * v(2)
    wpt(iw_mom(3)) = xi0 * v(3)
    wpt(iw_e)      = xi0 - p0 - rho0 * lfac0

  end subroutine analytic_state

  !> Assert that the ghost cells across the axis carry the exact solution.
  !>
  !> This is the only sharp check of the pole copy the case can make.  The
  !> volume averages the log reports cannot see it: a cell at the axis has
  !> vanishing volume, and the face lying on the axis has zero area, so the
  !> flux that would carry a wrong ghost value into the interior is multiplied
  !> by zero, and the jump a wrong sign leaves behind is clipped by the TVD
  !> limiter.  Reversing every sign at the pole moves mean(w) only in the
  !> fifth digit.  The ghost cells themselves are wrong by O(1) at once.
  !>
  !> Only meaningful at it = 0, where the interior is still exactly the
  !> analytic state, so the pole copy of it has to reproduce the analytic
  !> state in the ghost cells to round-off.  Silence means the check passed.
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
    double precision :: wpt(1:nw), x_loc(1:ndim), err
    integer          :: ix1, ix2, ix3, iside, i2, jxmin2, jxmax2

    if (it /= 0) return
    ! process() runs before the solution is pulled back for output, so the
    ! host copy of w is stale unless we fetch this block ourselves
    !$acc update host(ps(igrid)%w)

    err = 0.0d0
    do iside = 1, 2
       i2 = 2*iside - 3
       if (neighbor_pole(0,i2,0,igrid) == 0) cycle
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
                x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
                call analytic_state(x_loc, wpt)
                err = max(err, maxval(abs(w(ix1,ix2,ix3,&
                   1:nwflux) - wpt(1:nwflux))))
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
