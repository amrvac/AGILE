!> Relativistic blast wave on a spherical mesh.
!>
!> Same off-centre hot-spot placement as tests/hd/spherical_blast: a dense,
!> hot region is deliberately offset in r, theta and phi so that every
!> momentum component and every geometric source term is exercised. The
!> ambient medium starts at rest (four-velocity zero everywhere), so unlike
!> tests/srhd/spherical_uniform_flow this is a shock/wave-propagation
!> regression test rather than a well-balancing one.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> ambient state
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
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
  !$acc declare copyin(rho0, p0, pblast, rblast, fr, ftheta, fphi)

contains

  subroutine usr_init()

    ! set_coordinate_system is no longer required: the coordinate system comes
    ! from `geometry` in &meshlist. Calling it explicitly is still valid.

    usr_init_one_grid => initonegrid_usr

    call phys_activate()

  end subroutine usr_init

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
    integer                         :: ix1, ix2, ix3

    ! xprobmin2/xprobmax2 and xprobmin3/xprobmax3 are in radians by now:
    ! read_par_files converts them from the par file's units of 2*pi, and it
    ! runs before the grid (and hence this routine) is initialised.
    rc     = xprobmin1 + fr     * (xprobmax1 - xprobmin1)
    thetac = xprobmin2 + ftheta * (xprobmax2 - xprobmin2)
    phic   = xprobmin3 + fphi   * (xprobmax3 - xprobmin3)

    call check_blast_fits(rc, thetac, phic)

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
             if (d2 > rblast**2) then
                w(ix1,ix2,ix3,p_) = p0
             else
                w(ix1,ix2,ix3,p_) = pblast
             end if

             ! ambient medium at rest: four-velocity zero everywhere
             w(ix1,ix2,ix3,mom(1)) = 0.0d0
             w(ix1,ix2,ix3,mom(2)) = 0.0d0
             w(ix1,ix2,ix3,mom(3)) = 0.0d0

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

end module mod_usr
