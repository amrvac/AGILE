#:mute
#:include "mod_gpu_directives.fpp"
#:endmute
module mod_usr
  use mod_amrvac, only: usr_init_one_grid, usr_process_grid
  use mod_global_parameters, only: bc_data, levmax, ndim, nw, ps, typeboundary
  use mod_physics, only: phys_activate, rho_

  implicit none
  private

  public :: usr_init, specialbound_usr

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr
    usr_process_grid => check_inlet
    call phys_activate()
  end subroutine usr_init

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    real(dp), intent(in)            :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    real(dp), intent(inout)         :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    w(:,:,:,rho_) = 1.0_dp
  end subroutine initonegrid_usr

  !> Add a UFO to the VTK-filled inlet ghosts for the duration of the sun pulse.
  subroutine specialbound_usr(qt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, &
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB,w,x)
    ${GPU_ROUTINE_VECTOR()}$
    integer, intent(in) :: ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3
    integer, intent(in) :: ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3, iB
    real(dp), intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ndim)
    real(dp), intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,nw)
    real(dp) :: amplitude
    integer :: i, j, k

    amplitude = max(0.0_dp, min(1.0_dp, 16.0_dp * qt, 16.0_dp * (0.25_dp - qt)))
    if (amplitude == 0.0_dp) return

    ${GPU_LOOP_VECTOR('collapse(3)')}$
    do k = ixOmin3, ixOmax3
       do j = ixOmin2, ixOmax2
          do i = ixOmin1, ixOmax1
             w(i,j,k,rho_) = w(i,j,k,rho_) + amplitude * ufo_density(x(i,j,k,2), x(i,j,k,3))
          end do
       end do
    end do
  end subroutine specialbound_usr

  !> Saucer, upper dome, three portholes, and a striped beam in the y-z plane.
  pure function ufo_density(y, z) result(density)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in) :: y, z
    real(dp) :: density, hull, dome, windows, depth, width, beam

    hull = ellipse(y - 0.72_dp, z - 0.82_dp, 0.185_dp, 0.040_dp)
    dome = ellipse(y - 0.72_dp, z - 0.845_dp, 0.091_dp, 0.078_dp) * &
         max(0.0_dp, min(1.0_dp, (z - 0.82_dp) / 0.01_dp))
    windows = max(ellipse(y - 0.63_dp, z - 0.818_dp, 0.017_dp, 0.012_dp), &
                  ellipse(y - 0.72_dp, z - 0.818_dp, 0.017_dp, 0.012_dp), &
                  ellipse(y - 0.81_dp, z - 0.818_dp, 0.017_dp, 0.012_dp))

    depth = (0.79_dp - z) / 0.30_dp
    width = 0.025_dp + 0.12_dp * max(0.0_dp, depth)
    beam = max(0.0_dp, min(1.0_dp, 20.0_dp * depth, 20.0_dp * (1.0_dp - depth), &
               (width - abs(y - 0.72_dp)) / 0.012_dp))
    density = 0.48_dp * max(hull, 0.78_dp * dome) * (1.0_dp - 0.82_dp * windows) + &
         beam * (0.12_dp + 0.06_dp * cos(8.0_dp * pi * depth))
  end function ufo_density

  !> Smooth the ellipse's rim over a narrow band.
  pure function ellipse(y, z, radius_y, radius_z) result(mask)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in) :: y, z, radius_y, radius_z
    real(dp) :: mask

    mask = max(0.0_dp, min(1.0_dp, 8.0_dp * (1.0_dp - (y / radius_y)**2 - (z / radius_z)**2)))
  end function ellipse

  !> Check the ghost cells of the data face against the table, for either
  !> row option (interp_phy_first_row, bc_phy_first_row), with the UFO on top
  !> for bc_icarus.
  subroutine check_inlet(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, &
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters, only: bc_icarus, xprobmin1, xprobmin2, xprobmin3, &
         xprobmax1, xprobmax2, xprobmax3
    use mod_bc_data, only: bc_data_ix, bc_data_get_3d, interp_phy_first_row, &
         bc_phy_first_row
    integer, intent(in) :: igrid, level, ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3
    integer, intent(in) :: ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3
    real(dp), intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ndim)
    real(dp), intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,nw)
    integer :: face, idim, t1, t2, table, a, b, n, nfirst, nlast, pos, inner
    integer :: lo(3), hi(3), ilo(3), ihi(3), c(3), cp(3), ci(3)
    real(dp) :: amplitude, expected, data, xmin(3), xmax(3)
    logical :: icarus

    face = 0
    do n = 1, 2*ndim
       if (bc_data_ix(rho_,n) > 0) face = n
    end do
    if (face == 0) call mpistop('No VTK data face installed')
    table = bc_data_ix(rho_,face)
    icarus = typeboundary(rho_,face) == bc_icarus
    if (qt >= 0.125_dp .and. qt <= 0.25_dp .and. levmax < 2) &
         call mpistop('The VTK density pulse did not trigger AMR')

    idim = (face + 1) / 2
    t1 = merge(2, 1, idim == 1)
    t2 = merge(2, 3, idim == 3)
    lo = [ixOmin1, ixOmin2, ixOmin3]; hi = [ixOmax1, ixOmax2, ixOmax3]
    ilo = [ixImin1, ixImin2, ixImin3]; ihi = [ixImax1, ixImax2, ixImax3]
    xmin = [xprobmin1, xprobmin2, xprobmin3]; xmax = [xprobmax1, xprobmax2, xprobmax3]
    c = lo
    if (mod(face, 2) == 1) then
       c(idim) = ilo(idim)
       if (x(c(1),c(2),c(3),idim) > xmin(idim)) return
       nfirst = ilo(idim); nlast = lo(idim) - 1; pos = ilo(idim); inner = lo(idim)
       if (bc_phy_first_row) nlast = lo(idim)
    else
       c(idim) = ihi(idim)
       if (x(c(1),c(2),c(3),idim) < xmax(idim)) return
       nfirst = hi(idim) + 1; nlast = ihi(idim); pos = hi(idim) + 1; inner = hi(idim)
       if (bc_phy_first_row) nfirst = hi(idim)
    end if

    amplitude = max(0.0_dp, min(1.0_dp, 16.0_dp * qt, 16.0_dp * (0.25_dp - qt)))
    ${GPU_UPDATE_HOST('ps(igrid)%w')}$
    do b = lo(t2), hi(t2)
       do a = lo(t1), hi(t1)
          cp(t1) = a; cp(t2) = b; cp(idim) = pos
          ci = cp; ci(idim) = inner
          data = bc_data_get_3d(table, x(cp(1),cp(2),cp(3),t1), x(cp(1),cp(2),cp(3),t2), qt)
          do n = nfirst, nlast
             c = cp; c(idim) = n
             if (n == inner) then
                expected = data   ! bc_phy_first_row
             else if (interp_phy_first_row) then
                expected = 2 * data - w(ci(1),ci(2),ci(3),rho_)
             else
                expected = data
             end if
             if (icarus .and. n /= inner) expected = expected + amplitude * &
                  ufo_density(x(c(1),c(2),c(3),2), x(c(1),c(2),c(3),3))
             if (abs(w(c(1),c(2),c(3),rho_) - expected) > 1.0e-12_dp) then
                write(*,*) 'face', face, 'cell', c, 'w', w(c(1),c(2),c(3),rho_), 'expected', expected
                call mpistop('VTK inlet density differs from the expected sun or sun-plus-UFO pulse')
             end if
          end do
       end do
    end do
  end subroutine check_inlet

end module mod_usr
