#:mute
#:include "../../../src/mod_gpu_directives.fpp"
#:endmute

!> Pressure-balanced slab: rho0(x) = rho_c + rho_a cos(kx x), p0(x) = p_c + p_a cos(kx x),
!> B0 = (0, 0, B_z(x)) with B_z^2/2 = b_c^2/2 - p_a cos(kx x), so p0 + B_z^2/2 is uniform
!> and J0 x B0 = grad p0. Periodic, at rest. Holds to round-off only with mhd_equi_split.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision :: rho_c = 1.0d0
  double precision :: rho_a = 0.5d0
  double precision :: p_c = 0.5d0
  double precision :: p_a = 0.25d0
  double precision :: b_c = 1.0d0
  double precision :: kx = 6.283185307179586d0
  ${GPU_DECLARE_COPYIN('rho_c, rho_a, p_c, p_a, b_c, kx')}$

contains

  subroutine usr_init()
    call params_read_user(par_files)
    ${GPU_UPDATE_DEVICE('rho_c, rho_a, p_c, p_a, b_c, kx')}$

    usr_init_one_grid => initonegrid_usr
    usr_print_log     => slab_print_log

    call phys_activate()
  end subroutine usr_init

  subroutine params_read_user(files)
    use mod_global_parameters, only: unitpar
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ rho_c, rho_a, p_c, p_a, b_c, kx

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do
  end subroutine params_read_user

  pure subroutine usr_equilibrium(x, rho0, p0)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: rho0, p0
    rho0 = rho_c + rho_a * cos(kx * x(1))
    p0   = p_c + p_a * cos(kx * x(1))
  end subroutine usr_equilibrium

  pure subroutine usr_b0_field(x, b0)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: b0(3)
    b0(1:2) = 0.0d0
    b0(3) = sqrt(b_c**2 - 2.0d0 * p_a * cos(kx * x(1)))
  end subroutine usr_b0_field

  !> J0 = curl B0 = (0, -dB_z/dx, 0)
  pure subroutine usr_j0_field(x, j0)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: j0(3)
    double precision :: b0(3)
    call usr_b0_field(x, b0)
    j0(1) = 0.0d0
    j0(2) = -p_a * kx * sin(kx * x(1)) / b0(3)
    j0(3) = 0.0d0
  end subroutine usr_j0_field

  !> B0 where the code stores B1, zero otherwise: what diagnostics add to w(mag)
  pure subroutine background_field(x, b)
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: b(1:3)
#:if defined('B0SPLIT')
    call usr_b0_field(x, b)
#:else
    b(1:3) = 0.0d0
#:endif
  end subroutine background_field

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision                :: x_loc(1:ndim), rho0, p0, b(1:3), b0(1:3)
    integer                         :: ix1, ix2, ix3

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw) = 0.0d0
    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call usr_equilibrium(x_loc, rho0, p0)
             call usr_b0_field(x_loc, b0)
             call background_field(x_loc, b)
             w(ix1,ix2,ix3,rho_)       = rho0
             w(ix1,ix2,ix3,p_)         = p0
             w(ix1,ix2,ix3,mag(1:3))   = b0(1:3) - b(1:3)
          end do
       end do
    end do

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
  end subroutine initonegrid_usr

  !> it t dt vmax drho_max dp_max Etot: the deviations from the equilibrium and the total energy
  subroutine slab_print_log()
    use mod_global_parameters
    double precision   :: maxs_ipe(3), maxs(3), etot_ipe, etot
    double precision   :: x_loc(1:ndim), rho0, p0, b(1:3), rho, v2, p, dv
    integer            :: iigrid, igrid, ix1, ix2, ix3
    integer, parameter :: my_unit = 21
    logical, save      :: opened = .false.

    maxs_ipe = 0.0d0
    etot_ipe = 0.0d0
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       dv = rnode(rpdx1_,igrid) * rnode(rpdx2_,igrid) * rnode(rpdx3_,igrid)
       do ix3 = ixMlo3, ixMhi3
          do ix2 = ixMlo2, ixMhi2
             do ix1 = ixMlo1, ixMhi1
                x_loc(1:ndim) = ps(igrid)%x(ix1,ix2,ix3,1:ndim)
                call usr_equilibrium(x_loc, rho0, p0)
                call background_field(x_loc, b)
                rho = ps(igrid)%w(ix1,ix2,ix3,iw_rho)
                v2  = sum(ps(igrid)%w(ix1,ix2,ix3,iw_mom(1:3))**2) / rho**2
                p   = (phys_gamma - 1.0d0) * (ps(igrid)%w(ix1,ix2,ix3,iw_e) - 0.5d0 * rho * v2 &
                      - 0.5d0 * sum(ps(igrid)%w(ix1,ix2,ix3,iw_mag(1:3))**2))
                maxs_ipe(1) = max(maxs_ipe(1), sqrt(v2))
                maxs_ipe(2) = max(maxs_ipe(2), abs(rho - rho0))
                maxs_ipe(3) = max(maxs_ipe(3), abs(p - p0))
                etot_ipe = etot_ipe + (ps(igrid)%w(ix1,ix2,ix3,iw_e) &
                   + sum(b(1:3) * (0.5d0 * b(1:3) + ps(igrid)%w(ix1,ix2,ix3,iw_mag(1:3))))) * dv
             end do
          end do
       end do
    end do

    call MPI_ALLREDUCE(maxs_ipe, maxs, 3, MPI_DOUBLE_PRECISION, MPI_MAX, icomm, ierrmpi)
    call MPI_ALLREDUCE(etot_ipe, etot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierrmpi)

    if (mype == 0) then
       if (.not. opened) then
          if (restart_from_file == undefined) then
             open(unit=my_unit, file=trim(base_filename)//'.log', status='replace')
             write(my_unit,'(A)') '# it t dt vmax drho_max dp_max Etot'
          else
             open(unit=my_unit, file=trim(base_filename)//'.log', status='old', position='append')
          end if
          opened = .true.
       end if
       write(my_unit,'(I8,6ES16.8)') it, global_time, dt, maxs(1), maxs(2), maxs(3), etot
       flush(my_unit)
    end if
  end subroutine slab_print_log

end module mod_usr
