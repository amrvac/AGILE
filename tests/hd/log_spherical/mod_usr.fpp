!> Curvilinear tests on a spherical mesh, away from the polar axis: three runs
!> from one build.
!>
!> Every test directory rebuilds AGILE from scratch, so the suite's cost is
!> dominated by compilation rather than by the runs. Cases that agree on the
!> *compile-time* parameters - `phys`, `geometry` and `specialboundary`, the
!> ones make/config_schema.toml turns into fypp defines - can share one build
!> and differ only in their par file. These do, so they pick their initial
!> condition through the `setup` entry of the &usr_list namelist:
!>
!>   uflow.par       setup = 'uniform'  a uniform Cartesian flow
!>   blast.par       setup = 'blast'    a spherical blast wave, single level
!>   blast_amr.par   setup = 'blast'    the same blast, Lohner-refined
!>
!> The two setups probe the curvilinear machinery from opposite ends.
!>
!> 'uniform' is a constant density, constant pressure gas moving with a
!> constant Cartesian velocity, which is an exact steady solution of the Euler
!> equations. Written in spherical components the velocity is a non-trivial
!> function of theta and phi, so keeping this state exactly uniform exercises
!> every curvilinear flux and every geometric source term at once; any error
!> shows up immediately as drift in the volume averages the log reports.
!>
!> 'blast' is an over-pressured sphere in gas at rest, deliberately off-centre
!> in all three coordinates so that every momentum component is exercised. Its
!> AMR variant re-makes the grid throughout the run as the shell moves outward
!> in r, which is what puts the curvilinear prolongation, coarsening and
!> ghost-cell prolongation paths under test - none of which the single-level
!> run reaches.
!>
!> The domain deliberately stays away from the polar axis and from r = 0; the
!> axis itself is covered by tests/hd/spherical_pole.
!>
!> This directory is the logarithmically stretched variant: `geometry` is
!> 'logSpherical', so the radial mesh is uniform in ln(r) rather than in r and
!> the cell width grows in proportion to r. The domain spans r = 1 to r = 4, a
!> factor of about 3.8 in cell width from the inner edge to the outer one, so
!> the stretching is real rather than nominal.
!>
!> One consequence shows up directly in the code below: read_par_files stores
!> the radial domain bounds as ln(r), so xprobmin1 and xprobmax1 are logarithms
!> here and every use of them for a *physical* radius has to exponentiate.
!>
!> check_log_grid, registered as usr_process_grid, asserts the mesh itself at
!> it = 0. That check is the point of this directory: a stretched mesh that is
!> subtly wrong barely moves a volume average, so the metrics are compared
!> against their closed forms cell by cell instead.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> which initial condition to run: 'uniform' or 'blast'. Read from the
  !> &usr_list namelist of the par file; see params_read_user.
  character(len=20) :: setup = 'uniform'

  !> ambient state, shared by both setups
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
  !> the uniform Cartesian velocity. Zeroed for 'blast' in usr_init, which is
  !> what makes that setup a blast in gas at rest.
  double precision :: v0(3) = [1.0d0, 0.5d0, -0.3d0]

  !> 'blast' only: over-pressure and radius of the initial hot spot, in length
  !> units
  double precision :: pblast = 1.0d2
  double precision :: rblast = 0.2d0
  !> Centre of the hot spot, as a fraction of the domain extent in each
  !> coordinate. Deliberately off-centre in all three directions, so that every
  !> momentum component and every geometric source term is exercised.
  !>
  !> These are fractions rather than absolute (r, theta, phi) on purpose: the
  !> par file spells the angular bounds in units of 2*pi while this module works
  !> in radians, so hard-coding angles here invites placing the blast outside
  !> the domain. Deriving the centre from xprob* cannot get that wrong, and
  !> check_blast_fits below verifies the whole sphere is inside.
  double precision :: fr = 0.45d0, ftheta = 0.4d0, fphi = 0.4d0

  !> set to .true. to keep the uniform state, which must then stay at rest:
  !> the pressure part of the radial and polar fluxes has to cancel the
  !> geometric source terms exactly
  logical :: quiet_start = .false.

  !> the device-side form of `setup`: specialbound_usr runs on the GPU, where
  !> comparing a character string is awkward
  logical :: is_blast = .false.
  !$acc declare copyin(rho0, p0, v0, is_blast)

contains

  subroutine usr_init()
    use mod_comm_lib, only: mpistop

    ! the coordinate system comes from `geometry` in &meshlist

    call params_read_user(par_files)

    select case (trim(setup))
    case ('uniform')
       is_blast = .false.
    case ('blast')
       is_blast = .true.
       v0 = 0.0d0        ! the blast goes off in gas at rest
    case default
       call mpistop("&usr_list setup must be 'uniform' or 'blast'")
    end select
    !$acc update device(is_blast, v0)

    usr_init_one_grid => initonegrid_usr
    usr_process_grid  => check_log_grid
    usr_special_bc    => specialbound_usr

    call phys_activate()

  end subroutine usr_init

  !> Read the &usr_list namelist from every par file given on the command
  !> line. par_files is filled by read_arguments, which runs before usr_init.
  subroutine params_read_user(files)
    use mod_global_parameters, only: unitpar
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ setup, quiet_start

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read_user

  !> spherical velocity components of the uniform Cartesian flow at x.
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
    double precision                :: r, theta, phi, xcart1, xcart2, xcart3
    double precision                :: xc1, xc2, xc3, d2
    double precision                :: rc, thetac, phic
    double precision                :: v(1:3), x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    ! xprobmin2/xprobmax2 and xprobmin3/xprobmax3 are in radians by now:
    ! read_par_files converts them from the par file's units of 2*pi, and it
    ! runs before the grid (and hence this routine) is initialised.
    rc     = exp(xprobmin1) + fr     * (exp(xprobmax1) - exp(xprobmin1))
    thetac = xprobmin2 + ftheta * (xprobmax2 - xprobmin2)
    phic   = xprobmin3 + fphi   * (xprobmax3 - xprobmin3)

    if (is_blast) call check_blast_fits(rc, thetac, phic)

    xc1 = rc * sin(thetac) * cos(phic)
    xc2 = rc * sin(thetac) * sin(phic)
    xc3 = rc * cos(thetac)

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1

             r     = x(ix1,ix2,ix3,1)
             theta = x(ix1,ix2,ix3,2)
             phi   = x(ix1,ix2,ix3,3)

             w(ix1,ix2,ix3,rho_) = rho0
             w(ix1,ix2,ix3,p_)   = p0

             if (is_blast .and. .not.quiet_start) then
                ! the blast is spherical in physical space, so place it using
                ! Cartesian distances
                xcart1 = r * sin(theta) * cos(phi)
                xcart2 = r * sin(theta) * sin(phi)
                xcart3 = r * cos(theta)
                d2 = (xcart1-xc1)**2 + (xcart2-xc2)**2 + (xcart3-xc3)**2
                if (d2 <= rblast**2) w(ix1,ix2,ix3,p_) = pblast
             end if

             ! v0 is zero for the blast, so this leaves it at rest
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call uniform_velocity(x_loc, v)
             w(ix1,ix2,ix3,mom(1)) = v(1)
             w(ix1,ix2,ix3,mom(2)) = v(2)
             w(ix1,ix2,ix3,mom(3)) = v(3)

          end do
       end do
    end do

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  !> Abort if any part of the hot spot falls outside the domain. The angular
  !> half-widths are evaluated at the smallest radius and the smallest
  !> cylindrical radius the sphere reaches, which is the conservative choice.
  subroutine check_blast_fits(rc, thetac, phic)
    use mod_comm_lib, only: mpistop
    double precision, intent(in) :: rc, thetac, phic
    ! .. local ..
    double precision             :: rin, dtheta, dphi, sinmin

    rin = rc - rblast
    if (rin <= 0.0d0 .or. rc + rblast > exp(xprobmax1) .or. rin < exp(xprobmin1)) &
       call mpistop("spherical blast: hot spot sticks out of the radial domain")

    dtheta = asin(min(1.0d0, rblast/rin))
    if (thetac - dtheta < xprobmin2 .or. thetac + dtheta > xprobmax2) &
       call mpistop("spherical blast: hot spot sticks out in theta")

    sinmin = min(sin(thetac-dtheta), sin(thetac+dtheta))
    if (sinmin <= 0.0d0) call mpistop("spherical blast: hot spot reaches the &
       &polar axis; use tests/hd/spherical_pole for that")

    dphi = asin(min(1.0d0, rblast/(rin*sinmin)))
    if (phic - dphi < xprobmin3 .or. phic + dphi > xprobmax3) &
       call mpistop("spherical blast: hot spot sticks out in phi")

  end subroutine check_blast_fits

  !> keep the exact solution in the ghost cells, so that only the interior
  !> discretisation is under test. Used by the uniform setup; the blast leaves
  !> its boundaries open with 'cont' instead.
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
    double precision                :: v(1:3)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    !$acc loop collapse(3) vector private(v, x_loc)
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call uniform_velocity(x_loc, v)

             w(ix1,ix2,ix3,iw_rho)    = rho0
             w(ix1,ix2,ix3,iw_mom(1)) = rho0 * v(1)
             w(ix1,ix2,ix3,iw_mom(2)) = rho0 * v(2)
             w(ix1,ix2,ix3,iw_mom(3)) = rho0 * v(3)
             w(ix1,ix2,ix3,iw_e)      = p0 * inv_gamma_1 + &
                0.5d0 * rho0 * sum(v(1:3)**2)

          end do
       end do
    end do

  end subroutine specialbound_usr


  !> Assert the logarithmically stretched mesh against its closed form, at
  !> it = 0.
  !>
  !> A wrong metric on a stretched grid is close to invisible in the volume
  !> averages the log reports - the same reason the polar-axis cases check
  !> ghost cells directly - so the mesh is checked cell by cell here instead.
  !> Four things, all to round-off:
  !>
  !>   1. consecutive radial positions are in exact geometric progression,
  !>      x(i+1)/x(i) = q = exp(dxi). That is the defining property of a
  !>      logarithmic grid, and it holds for the volume *barycentres* stored in
  !>      x as well as for the faces, because both faces of every cell scale by
  !>      the same q. It is also the property that fails loudly if the
  !>      exponential is ever dropped or applied twice.
  !>   2. ds(1) is the physical width of the cell, and ds(1)/x(1) is therefore
  !>      the same constant everywhere - constant relative resolution.
  !>   3. each stored position is the volume barycentre of its own two faces.
  !>      The faces are rebuilt here from rnode, and the barycentre from the
  !>      raw quotient of quartic (spherical) or cubic (cylindrical)
  !>      differences - deliberately a different expression from the
  !>      cancellation-free one fill_geometry_device evaluates, so that the two
  !>      agreeing means something.
  !>   4. dvolume is the exact integral over those same faces, with the polar
  !>      factor written as cos(theta_L) - cos(theta_R) rather than the
  !>      2 sin(theta_c) sin(dtheta/2) the kernel uses.
  subroutine check_log_grid(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    ! .. local ..
    double precision :: d1, d2, d3, xlo1, xlo2, q, fL, fR, rbar, dvol
    double precision :: tL, tR, err
    double precision, parameter :: tol = 1.0d-11
    integer          :: ix1, ix2, ix3

    if (it /= 0) return
    ! process() syncs the positions for us, but not the metrics
    !$acc update host(ps(igrid)%ds, ps(igrid)%dvolume)

    d1 = rnode(rpdx1_,igrid); d2 = rnode(rpdx2_,igrid); d3 = rnode(rpdx3_,igrid)
    xlo1 = rnode(rpxmin1_,igrid); xlo2 = rnode(rpxmin2_,igrid)
    q = exp(d1)

    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

             ! the two radial faces, from the block corner in ln(r)
             fL = exp(xlo1 + dble(ix1-nghostcells-1)*d1)
             fR = exp(xlo1 + dble(ix1-nghostcells)*d1)
             tL = xlo2 + dble(ix2-nghostcells-1)*d2
             tR = tL + d2

             ! (1) geometric progression of the stored positions
             if (ix1 < ixOmax1) then
                err = abs(x(ix1+1,ix2,ix3,1)/x(ix1,ix2,ix3,1) - q)/q
                call assert_close(err, tol, 'radial positions are not in &
                   &geometric progression', x(ix1,ix2,ix3,1))
             end if

             ! (2) constant relative resolution
             err = abs(ps(igrid)%ds(ix1,ix2,ix3,1) - (fR-fL))/(fR-fL)
             call assert_close(err, tol, 'ds(1) is not the physical cell &
                &width', x(ix1,ix2,ix3,1))

             ! (3) the stored position is the volume barycentre
             rbar = 0.75d0*(fR**4 - fL**4)/(fR**3 - fL**3)
             err = abs(x(ix1,ix2,ix3,1) - rbar)/rbar
             call assert_close(err, tol, 'x(1) is not the volume barycentre &
                &of its faces', x(ix1,ix2,ix3,1))

             ! (4) the cell volume is the exact integral over those faces
             dvol = (fR**3 - fL**3)/3.0d0 * (cos(tL) - cos(tR)) * d3
             err = abs(ps(igrid)%dvolume(ix1,ix2,ix3) - dvol)/dvol
             call assert_close(err, tol, 'dvolume is not the exact cell &
                &volume', x(ix1,ix2,ix3,1))

          end do
       end do
    end do

  end subroutine check_log_grid

  !> Abort with a located message if err exceeds tol.
  subroutine assert_close(err, tol, what, r)
    use mod_global_parameters, only: mype, unitterm
    use mod_comm_lib, only: mpistop
    double precision, intent(in) :: err, tol, r
    character(len=*), intent(in) :: what

    if (err <= tol) return
    if (mype == 0) write(unitterm,*) 'check_log_grid: ', trim(what), &
       ' at r =', r, ' relative error', err
    call mpistop('check_log_grid: the logarithmic mesh is wrong')

  end subroutine assert_close

end module mod_usr
