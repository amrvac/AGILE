!> Curvilinear FFHD tests on a spherical mesh, away from the axis: two runs
!> from one build.
!>
!> Cases that agree on the compile-time parameters - `phys`, `geometry` and
!> `specialboundary` - share one build and differ only in their par file,
!> which keeps the suite's cost down since every directory otherwise rebuilds
!> AGILE from scratch. They pick their initial condition through the `setup`
!> entry of the &usr_list namelist:
!>
!>   uflow.par   setup = 'uniform'  a uniform state along a uniform frozen field
!>   blast.par   setup = 'blast'    a blast wave along the same field
!>
!> FFHD advects genuine scalars along a user-supplied field direction b-hat, so
!> 'uniform' is a constant state moving with a constant field-aligned speed -
!> an exact steady solution, and the well-balancing test for the geometric
!> source terms. 'blast' releases an over-pressured sphere into the same
!> frozen field from rest.
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
  !> field-aligned velocity. Zeroed for 'blast' in usr_init, which is what
  !> makes that setup a blast in gas at rest.
  double precision :: vpar0 = 0.4d0
  !> overpressure and radius of the initial hot spot, in length units
  double precision :: pblast = 1.0d2
  double precision :: rblast = 0.2d0
  !> the frozen field direction, in Cartesian components (need not be unit
  !> length: only its direction matters, normalised in to_spherical_unit)
  double precision :: b0(3) = [1.0d0, 0.5d0, -0.3d0]
  !> Centre of the hot spot, as a fraction of the domain extent in each
  !> coordinate. Deliberately off-centre in all three directions, so that
  !> the field-aligned momentum equation and the frozen field's curvature
  !> are both exercised.
  !>
  !> These are fractions rather than absolute (r, theta, phi) on purpose: the
  !> par file spells the angular bounds in units of 2*pi while this module works
  !> in radians, so hard-coding angles here invites placing the blast outside
  !> the domain. Deriving the centre from xprob* cannot get that wrong, and
  !> check_blast_fits below verifies the whole sphere is inside.
  double precision :: fr = 0.45d0, ftheta = 0.4d0, fphi = 0.4d0
  !$acc declare copyin(is_blast, rho0, p0, pblast, rblast, b0, fr, ftheta, fphi)

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
       vpar0 = 0.0d0     ! the blast goes off in gas at rest
    case default
       call mpistop("&usr_list setup must be 'uniform' or 'blast'")
    end select
    !$acc update device(is_blast, vpar0)

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

  !> spherical components at x of the unit vector along Cartesian vec0.
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
    double precision                :: bhat(1:3)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    ! xprobmin2/xprobmax2 and xprobmin3/xprobmax3 are in radians by now:
    ! read_par_files converts them from the par file's units of 2*pi, and it
    ! runs before the grid (and hence this routine) is initialised.
    rc     = xprobmin1 + fr     * (xprobmax1 - xprobmin1)
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

             ! the hot spot is spherical in physical space, so place it
             ! using Cartesian distances
             xcart1 = r * sin(theta) * cos(phi)
             xcart2 = r * sin(theta) * sin(phi)
             xcart3 = r * cos(theta)
             d2 = (xcart1-xc1)**2 + (xcart2-xc2)**2 + (xcart3-xc3)**2

             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call to_spherical_unit(x_loc, b0, bhat)

             w(ix1,ix2,ix3,rho_) = rho0
             w(ix1,ix2,ix3,p_) = p0

             if (is_blast .and. .not.quiet_start .and. d2 <= rblast**2) &

                w(ix1,ix2,ix3,p_) = pblast

             ! ambient medium at rest: zero field-aligned momentum everywhere
             ! vpar0 is zero for the blast, so this leaves it at rest
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
    double precision                :: bhat(1:3)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    !$acc loop collapse(3) vector private(bhat, x_loc)
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1

             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call to_spherical_unit(x_loc, b0, bhat)

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
