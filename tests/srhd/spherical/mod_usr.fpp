!> Curvilinear SRHD tests on a spherical mesh, away from the axis: two runs
!> from one build.
!>
!> Cases that agree on the compile-time parameters - `phys`, `geometry` and
!> `specialboundary` - share one build and differ only in their par file,
!> which keeps the suite's cost down since every directory otherwise rebuilds
!> AGILE from scratch. They pick their initial condition through the `setup`
!> entry of the &usr_list namelist:
!>
!>   uflow.par   setup = 'uniform'  a uniform, sub-luminal Cartesian flow
!>   blast.par   setup = 'blast'    a relativistic blast wave in gas at rest
!>
!> 'uniform' is an exact steady solution: the Lorentz factor depends only on
!> |v|, which a rotation of the coordinate frame leaves invariant, so it stays
!> a single global constant however the components are written. That makes it
!> the well-balancing test for SRHD's geometric source terms.
!>
!> The domain deliberately stays away from the axis.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> which initial condition to run: 'uniform' or 'blast'. Read from the
  !> &usr_list namelist of the par file; see params_read_user.
  character(len=20) :: setup = 'uniform'

  !> suppress the hot spot, leaving the ambient state, which must then stay at
  !> rest: the pressure part of the fluxes has to cancel the geometric source
  !> terms exactly
  logical :: quiet_start = .false.

  !> the device-side form of `setup`: specialbound_usr runs on the GPU, where
  !> comparing a character string is awkward
  logical :: is_blast = .false.


  !> ambient state
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
  !> the uniform Cartesian velocity, sub-luminal. Zeroed for 'blast' in
  !> usr_init, which is what makes that setup a blast in gas at rest.
  double precision :: v0(3) = [0.3d0, 0.15d0, -0.1d0]
  !> overpressure and radius of the initial hot spot, in length units
  double precision :: pblast = 1.0d1
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
  !$acc declare copyin(is_blast, v0, rho0, p0, pblast, rblast, fr, ftheta, fphi)

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
    double precision                :: r, theta, phi, xcart1, xcart2, xcart3
    double precision                :: xc1, xc2, xc3, d2
    double precision                :: rc, thetac, phic
    double precision                :: v(1:3), lfac0, x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    ! xprobmin2/xprobmax2 and xprobmin3/xprobmax3 are in radians by now:
    ! read_par_files converts them from the par file's units of 2*pi, and it
    ! runs before the grid (and hence this routine) is initialised.
    rc     = xprobmin1 + fr     * (xprobmax1 - xprobmin1)
    thetac = xprobmin2 + ftheta * (xprobmax2 - xprobmin2)
    phic   = xprobmin3 + fphi   * (xprobmax3 - xprobmin3)

    lfac0 = 1.0d0 / sqrt(1.0d0 - sum(v0(1:3)**2))

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

             ! the blast is spherical in physical space, so place it using
             ! Cartesian distances
             xcart1 = r * sin(theta) * cos(phi)
             xcart2 = r * sin(theta) * sin(phi)
             xcart3 = r * cos(theta)
             d2 = (xcart1-xc1)**2 + (xcart2-xc2)**2 + (xcart3-xc3)**2

             w(ix1,ix2,ix3,rho_) = rho0
             w(ix1,ix2,ix3,p_) = p0

             if (is_blast .and. .not.quiet_start .and. d2 <= rblast**2) &

                w(ix1,ix2,ix3,p_) = pblast

             ! ambient medium at rest: four-velocity zero everywhere
             ! v0 is zero for the blast, so this leaves it at rest
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call uniform_velocity(x_loc, v)
             w(ix1,ix2,ix3,mom(1)) = lfac0 * v(1)
             w(ix1,ix2,ix3,mom(2)) = lfac0 * v(2)
             w(ix1,ix2,ix3,mom(3)) = lfac0 * v(3)

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
    if (rin <= 0.0d0 .or. rc + rblast > xprobmax1 .or. rin < xprobmin1) &
       call mpistop("spherical_blast: hot spot sticks out of the radial domain")

    dtheta = asin(min(1.0d0, rblast/rin))
    if (thetac - dtheta < xprobmin2 .or. thetac + dtheta > xprobmax2) &
       call mpistop("spherical_blast: hot spot sticks out in theta")

    sinmin = min(sin(thetac-dtheta), sin(thetac+dtheta))
    if (sinmin <= 0.0d0) call mpistop("spherical_blast: hot spot reaches the &
       &polar axis, which AGILE cannot handle")

    dphi = asin(min(1.0d0, rblast/(rin*sinmin)))
    if (phic - dphi < xprobmin3 .or. phic + dphi > xprobmax3) &
       call mpistop("spherical_blast: hot spot sticks out in phi")

  end subroutine check_blast_fits

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

end module mod_usr
