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

  !> Check the sun lookup and the UFO addition on the driven face.
  subroutine check_inlet(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, &
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    use mod_comm_lib, only: mpistop
    use mod_bc_data, only: bc_data_ix, bc_data_get_3d
    integer, intent(in) :: igrid, level, ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3
    integer, intent(in) :: ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3
    real(dp), intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ndim)
    real(dp), intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,nw)
    integer :: face, table, first, last, j, k
    real(dp) :: amplitude, expected, y, z

    ! The bc_icarus run enters from the high side to exercise its face index.
    face = 1
    if (typeboundary(rho_,1) /= bc_data) face = 2
    table = bc_data_ix(rho_,face)
    if (table < 1) call mpistop('VTK inlet installed on the wrong face')
    if (qt >= 0.125_dp .and. qt <= 0.25_dp .and. levmax < 2) &
         call mpistop('The VTK density pulse did not trigger AMR')

    if (face == 2) then
       first = ixOmax1 + 1
       last = ixImax1
       if (x(first,ixOmin2,ixOmin3,1) < 1.0_dp) return
    else
       first = ixImin1
       last = ixOmin1 - 1
       if (x(last,ixOmin2,ixOmin3,1) > 0.0_dp) return
    end if

    amplitude = max(0.0_dp, min(1.0_dp, 16.0_dp * qt, 16.0_dp * (0.25_dp - qt)))
    ${GPU_UPDATE_HOST('ps(igrid)%w')}$
    do k = ixOmin3, ixOmax3
       do j = ixOmin2, ixOmax2
          y = x(first,j,k,2)
          z = x(first,j,k,3)
          expected = bc_data_get_3d(table, y, z, qt)
          if (face == 2) expected = expected + amplitude * ufo_density(y, z)
          if (.not. all(abs(w(first:last,j,k,rho_) - expected) <= 1.0e-12_dp)) &
               call mpistop('VTK inlet density differs from the expected sun or sun-plus-UFO pulse')
       end do
    end do
  end subroutine check_inlet

end module mod_usr
