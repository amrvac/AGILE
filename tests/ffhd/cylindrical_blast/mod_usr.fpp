!> Field-aligned pressure pulse on a cylindrical mesh, threaded by a uniform
!> Cartesian frozen field.
!>
!> Same off-centre hot-spot placement as tests/hd/cylindrical_blast, extended
!> with the uniform frozen field direction from
!> tests/ffhd/cylindrical_uniform_flow (Cartesian, so it is a non-trivial
!> function of phi once written in cylindrical components). Unlike a full
!> HD/MHD blast, ffhd only has a single field-aligned momentum component, so
!> the perturbation can only launch waves along the (curved, spatially
!> varying) field direction rather than in every direction; this is a
!> shock/AMR regression test rather than a well-balancing one
!> (tests/ffhd/cylindrical_uniform_flow covers that, including the frozen
!> field's inter-block ghost-cell exchange, which this relies on too since
!> multiple blocks are used here and the field is not spatially uniform).
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> ambient state
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
  !> overpressure and radius of the initial hot spot, in length units
  double precision :: pblast = 1.0d2
  double precision :: rblast = 0.2d0
  !> the frozen field direction, in Cartesian components (need not be unit
  !> length: only its direction matters, normalised in to_cylindrical_unit)
  double precision :: b0(3) = [1.0d0, 0.5d0, -0.3d0]
  !> Centre of the hot spot, as a fraction of the domain extent in each
  !> coordinate. Deliberately off-centre in all three directions, so that
  !> the field-aligned momentum equation and the frozen field's curvature
  !> are both exercised.
  !>
  !> These are fractions rather than absolute (r, z, phi) on purpose: the par
  !> file spells the phi bound in units of 2*pi while this module works in
  !> radians, so hard-coding an angle here invites placing the blast outside
  !> the domain. Deriving the centre from xprob* cannot get that wrong, and
  !> check_blast_fits below verifies the whole sphere is inside.
  double precision :: fr = 0.45d0, fz = 0.4d0, fphi = 0.4d0
  !$acc declare copyin(rho0, p0, pblast, rblast, b0, fr, fz, fphi)

contains

  subroutine usr_init()

    ! set_coordinate_system is no longer required: the coordinate system comes
    ! from `geometry` in &meshlist. Calling it explicitly is still valid.

    usr_init_one_grid => initonegrid_usr

    call phys_activate()

  end subroutine usr_init

  !> cylindrical (r, z, phi) components at x of the unit vector along
  !> Cartesian vec0.
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
    double precision                :: r, z, phi, xcart1, xcart2, xcart3
    double precision                :: xc1, xc2, xc3, d2
    double precision                :: rc, zc, phic
    double precision                :: bhat(1:3)
    double precision                :: x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

    ! xprobmin3/xprobmax3 is in radians by now: read_par_files converts it
    ! from the par file's units of 2*pi, and it runs before the grid (and
    ! hence this routine) is initialised.
    rc   = xprobmin1 + fr   * (xprobmax1 - xprobmin1)
    zc   = xprobmin2 + fz   * (xprobmax2 - xprobmin2)
    phic = xprobmin3 + fphi * (xprobmax3 - xprobmin3)

    call check_blast_fits(rc, zc, phic)

    xc1 = rc * cos(phic)
    xc2 = rc * sin(phic)
    xc3 = zc

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1

             r   = x(ix1,ix2,ix3,1)
             z   = x(ix1,ix2,ix3,2)
             phi = x(ix1,ix2,ix3,3)

             ! the hot spot is spherical in physical space, so place it
             ! using Cartesian distances
             xcart1 = r * cos(phi)
             xcart2 = r * sin(phi)
             xcart3 = z
             d2 = (xcart1-xc1)**2 + (xcart2-xc2)**2 + (xcart3-xc3)**2

             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call to_cylindrical_unit(x_loc, b0, bhat)

             w(ix1,ix2,ix3,rho_) = rho0
             if (d2 > rblast**2) then
                w(ix1,ix2,ix3,p_) = p0
             else
                w(ix1,ix2,ix3,p_) = pblast
             end if

             ! ambient medium at rest: zero field-aligned momentum everywhere
             w(ix1,ix2,ix3,mom(1)) = 0.0d0
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
  !> half-width is evaluated at the smallest radius the sphere reaches, which
  !> is the conservative choice; z is a plain length so its check is direct.
  subroutine check_blast_fits(rc, zc, phic)
    use mod_comm_lib, only: mpistop
    double precision, intent(in) :: rc, zc, phic
    ! .. local ..
    double precision             :: rin, dphi

    rin = rc - rblast
    if (rin <= 0.0d0 .or. rc + rblast > xprobmax1 .or. rin < xprobmin1) &
       call mpistop("cylindrical_blast: hot spot sticks out of the radial domain")

    if (zc - rblast < xprobmin2 .or. zc + rblast > xprobmax2) &
       call mpistop("cylindrical_blast: hot spot sticks out in z")

    dphi = asin(min(1.0d0, rblast/rin))
    if (phic - dphi < xprobmin3 .or. phic + dphi > xprobmax3) &
       call mpistop("cylindrical_blast: hot spot sticks out in phi")

  end subroutine check_blast_fits

end module mod_usr
