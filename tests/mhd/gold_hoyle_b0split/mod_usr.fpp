#:mute
#:include "../../../src/mod_gpu_directives.fpp"
#:endmute

!> Gold-Hoyle force-free flux rope B_z = B0/(1+q^2 r^2), B_phi = B0 q r/(1+q^2 r^2),
!> uniform rho and p: an exact static equilibrium for any beta. setup = 'rope'
!> holds it. setup = 'bump' and 'helix' are a Gaussian pressure bump in a uniform
!> field and in the periodic force-free field (cos kz, sin kz, 0), run here and in
!> the B0-split twin tests/mhd/gold_hoyle_b0split for comparison.
!> With B0SPLIT the analytic field is B0 and the evolved B1 starts at zero.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  character(len=20) :: setup = 'rope'
  integer :: background = 0
  logical :: is_bump = .false.
  double precision :: rho0 = 1.0d0
  double precision :: b0 = 1.0d0
  double precision :: beta = 0.1d0
  double precision :: qtwist = 1.0d0
  double precision :: p0 = 0.05d0
  double precision :: pbump = 0.05d0
  double precision :: rbump = 0.5d0
  double precision :: bcart(3) = [0.6d0, 0.3d0, 0.8d0]
  double precision :: khelix = 6.283185307179586d0
  ${GPU_DECLARE_COPYIN('background, is_bump, rho0, b0, qtwist, p0, pbump, rbump, bcart, khelix')}$

contains

  subroutine usr_init()
    use mod_global_parameters, only: nwauxio
    use mod_comm_lib, only: mpistop

    call params_read_user(par_files)
    select case (trim(setup))
    case ('rope');  background = 0
    case ('bump');  background = 1
    case ('helix'); background = 2
    case default;   call mpistop("&usr_list setup must be 'rope', 'bump' or 'helix'")
    end select
    is_bump = background > 0
    p0 = 0.5d0 * beta * b0**2
    ${GPU_UPDATE_DEVICE('background, is_bump, rho0, b0, qtwist, p0, pbump, rbump, bcart, khelix')}$

    nwauxio = 3
    usr_init_one_grid    => initonegrid_usr
    usr_special_bc       => specialbound_usr
    usr_before_main_loop => print_setup
    usr_print_log        => gh_print_log
    usr_aux_output       => gh_aux_output
    usr_add_aux_names    => gh_aux_names

    call phys_activate()

  end subroutine usr_init

  subroutine params_read_user(files)
    use mod_global_parameters, only: unitpar
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ setup, rho0, b0, beta, qtwist, pbump, rbump, bcart, khelix

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read_user

  !> the analytic field: (B_x, B_y, B_z), or (B_r, B_z, B_phi) on a cylindrical mesh
  pure subroutine field(x, b)
#ifdef _OPENMP
    ${GPU_ROUTINE_VECTOR()}$
#else
    ${GPU_ROUTINE_SEQ()}$
#endif
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: b(1:3)
    double precision              :: f

    if (background == 1) then
       b(1:3) = bcart(1:3)
       return
    else if (background == 2) then
       b(1:3) = b0 * [cos(khelix * x(3)), sin(khelix * x(3)), 0.0d0]
       return
    end if
#:if GEOM == 'cylindrical'
    f = b0 / (1.0d0 + (qtwist * x(1))**2)
    b(1) = 0.0d0
    b(2) = f
    b(3) = f * qtwist * x(1)
#:else
    f = b0 / (1.0d0 + qtwist**2 * (x(1)**2 + x(2)**2))
    b(1) = -f * qtwist * x(2)
    b(2) =  f * qtwist * x(1)
    b(3) =  f
#:endif

  end subroutine field

  !> curl of the analytic field: alpha B with alpha = 2q/(1+q^2 r^2) for the rope, -k for the helix
  pure subroutine current(x, j)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: j(1:3)
    double precision              :: b(1:3), r2

    call field(x, b)
    if (background == 1) then
       j(1:3) = 0.0d0
       return
    else if (background == 2) then
       j(1:3) = -khelix * b(1:3)
       return
    end if
#:if GEOM == 'cylindrical'
    r2 = x(1)**2
#:else
    r2 = x(1)**2 + x(2)**2
#:endif
    j(1:3) = 2.0d0 * qtwist / (1.0d0 + qtwist**2 * r2) * b(1:3)

  end subroutine current

  pure subroutine usr_b0_field(x, b0f)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: b0f(3)
    call field(x, b0f)
  end subroutine usr_b0_field

  pure subroutine usr_j0_field(x, j0f)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: j0f(3)
    call current(x, j0f)
  end subroutine usr_j0_field

  !> the evolved field at t = 0: the analytic field, or nothing on top of B0
  pure subroutine evolved_field(x, b)
#ifdef _OPENMP
    ${GPU_ROUTINE_VECTOR()}$
#else
    ${GPU_ROUTINE_SEQ()}$
#endif
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: b(1:3)
#:if defined('B0SPLIT')
    b(1:3) = 0.0d0
#:else
    call field(x, b)
#:endif
  end subroutine evolved_field

  !> B0 where the code stores B1, zero otherwise: what diagnostics add to w(mag)
  pure subroutine background_field(x, b)
    double precision, intent(in)  :: x(1:ndim)
    double precision, intent(out) :: b(1:3)
#:if defined('B0SPLIT')
    call field(x, b)
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
    double precision                :: b(1:3), x_loc(1:ndim), xc(1:ndim), d2
    integer                         :: ix1, ix2, ix3

    xc(1) = 0.5d0 * (xprobmin1 + xprobmax1)
    xc(2) = 0.5d0 * (xprobmin2 + xprobmax2)
    xc(3) = 0.5d0 * (xprobmin3 + xprobmax3)
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw) = 0.0d0
    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call evolved_field(x_loc, b)
             d2 = sum((x_loc(1:ndim) - xc(1:ndim))**2)
             w(ix1,ix2,ix3,rho_)   = rho0
             w(ix1,ix2,ix3,p_)     = p0
             if (is_bump) w(ix1,ix2,ix3,p_) = p0 + pbump * exp(-d2 / rbump**2)
             w(ix1,ix2,ix3,mag(1)) = b(1)
             w(ix1,ix2,ix3,mag(2)) = b(2)
             w(ix1,ix2,ix3,mag(3)) = b(3)
          end do
       end do
    end do

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  !> the exact equilibrium held in the ghost cells
  subroutine specialbound_usr(qt, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB, w, x)
    use mod_global_parameters
    ${GPU_ROUTINE_VECTOR()}$
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
    double precision                :: b(1:3), x_loc(1:ndim)
    integer                         :: ix1, ix2, ix3

#ifndef _OPENMP
    ${GPU_LOOP_VECTOR("collapse(3) private(b, x_loc)")}$
#endif
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call evolved_field(x_loc, b)
             w(ix1,ix2,ix3,iw_rho)    = rho0
             w(ix1,ix2,ix3,iw_mom(1)) = 0.0d0
             w(ix1,ix2,ix3,iw_mom(2)) = 0.0d0
             w(ix1,ix2,ix3,iw_mom(3)) = 0.0d0
             w(ix1,ix2,ix3,iw_mag(1)) = b(1)
             w(ix1,ix2,ix3,iw_mag(2)) = b(2)
             w(ix1,ix2,ix3,iw_mag(3)) = b(3)
             w(ix1,ix2,ix3,psi_)      = 0.0d0
             w(ix1,ix2,ix3,iw_e)      = p0 * inv_gamma_1 + 0.5d0 * sum(b(1:3)**2)
          end do
       end do
    end do

  end subroutine specialbound_usr

  subroutine print_setup()
    use mod_global_parameters
    double precision :: lz

#:if GEOM == 'cylindrical'
    lz = xprobmax2 - xprobmin2
#:else
    lz = xprobmax3 - xprobmin3
#:endif
    if (mype /= 0) return
    write(*,'(A,A)')       ' Gold-Hoyle setup = ', trim(setup)
    write(*,'(A,4ES12.5)') '   rho0, B0, p0, q      = ', rho0, b0, p0, qtwist
    write(*,'(A,2ES12.5)') '   tau_A, time_max      = ', sqrt(rho0) / (qtwist * b0), time_max
    write(*,'(A,ES12.5)')  '   total twist / 2pi    = ', qtwist * lz / (2.0d0 * dpi)
#:if defined('B0SPLIT')
    write(*,'(A)')         '   B0 split: B1 = 0 at t = 0'
#:endif

  end subroutine print_setup

  !> |j x B|, |j| |B| and div B on the interior, with B0 and J0 added in a split build
  subroutine gh_lorentz(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,jxb,jb,divb)
    use mod_global_parameters
    use mod_functions_bfield, only: get_current, get_divb
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out) :: jxb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out) :: jb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out) :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: current_w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)
    double precision :: j(1:3), b(1:3), f(1:3), b0f(1:3), j0f(1:3), x_loc(1:ndim)
    integer          :: idirmin, ix1, ix2, ix3

    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current_w)
    call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)

    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
             x_loc(1:ndim) = x(ix1,ix2,ix3,1:ndim)
             call background_field(x_loc, b0f)
             j0f(1:3) = 0.0d0
#:if defined('B0SPLIT')
             call current(x_loc, j0f)
#:endif
             j(1:3) = current_w(ix1,ix2,ix3,1:3) + j0f(1:3)
             b(1:3) = w(ix1,ix2,ix3,iw_mag(1:3)) + b0f(1:3)
             f(1) = j(2)*b(3) - j(3)*b(2)
             f(2) = j(3)*b(1) - j(1)*b(3)
             f(3) = j(1)*b(2) - j(2)*b(1)
             jxb(ix1,ix2,ix3) = sqrt(sum(f(1:3)**2))
             jb(ix1,ix2,ix3)  = sqrt(sum(j(1:3)**2)) * sqrt(sum(b(1:3)**2))
          end do
       end do
    end do

  end subroutine gh_lorentz

  !> it t dt Ekin Emag vmax dB_L2 drho_max sigmaJ jxb_max divb_mean divb_max Etot
  subroutine gh_print_log()
    use mod_global_parameters
    integer, parameter :: nsum = 9, nmax = 4
    double precision   :: jxb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision   :: jb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision   :: divb(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision   :: sums_ipe(nsum), sums(nsum), maxs_ipe(nmax), maxs(nmax)
    double precision   :: dv, rho, v2, b(1:3), b0f(1:3), bex(1:3), x_loc(1:ndim)
    integer            :: iigrid, igrid, ix1, ix2, ix3
    integer, parameter :: my_unit = 21
    logical, save      :: opened = .false.
    character(len=80)  :: filename

    sums_ipe = 0.0d0
    maxs_ipe = 0.0d0

    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       block => ps(igrid)
       dxlevel(1) = rnode(rpdx1_,igrid)
       dxlevel(2) = rnode(rpdx2_,igrid)
       dxlevel(3) = rnode(rpdx3_,igrid)

       call gh_lorentz(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid)%w,ps(igrid)%x,jxb,jb,divb)

       do ix3 = ixMlo3, ixMhi3
          do ix2 = ixMlo2, ixMhi2
             do ix1 = ixMlo1, ixMhi1
                if (slab_uniform) then
                   dv = dxlevel(1) * dxlevel(2) * dxlevel(3)
                else
                   dv = block%dvolume(ix1,ix2,ix3)
                end if
                x_loc(1:ndim) = ps(igrid)%x(ix1,ix2,ix3,1:ndim)
                call background_field(x_loc, b0f)
                call field(x_loc, bex)
                rho    = ps(igrid)%w(ix1,ix2,ix3,iw_rho)
                v2     = sum(ps(igrid)%w(ix1,ix2,ix3,iw_mom(1:3))**2) / rho**2
                b(1:3) = ps(igrid)%w(ix1,ix2,ix3,iw_mag(1:3)) + b0f(1:3)

                sums_ipe(1) = sums_ipe(1) + 0.5d0 * rho * v2 * dv
                sums_ipe(2) = sums_ipe(2) + 0.5d0 * sum(b(1:3)**2) * dv
                sums_ipe(3) = sums_ipe(3) + sum((b(1:3) - bex(1:3))**2) * dv
                sums_ipe(4) = sums_ipe(4) + sum(bex(1:3)**2) * dv
                sums_ipe(5) = sums_ipe(5) + jxb(ix1,ix2,ix3) * dv
                sums_ipe(6) = sums_ipe(6) + jb(ix1,ix2,ix3) * dv
                sums_ipe(7) = sums_ipe(7) + abs(divb(ix1,ix2,ix3)) * dv
                sums_ipe(8) = sums_ipe(8) + dv
                sums_ipe(9) = sums_ipe(9) + (ps(igrid)%w(ix1,ix2,ix3,iw_e) &
                   + sum(b0f(1:3) * (0.5d0 * b0f(1:3) + ps(igrid)%w(ix1,ix2,ix3,iw_mag(1:3))))) * dv
                maxs_ipe(1) = max(maxs_ipe(1), sqrt(v2))
                maxs_ipe(2) = max(maxs_ipe(2), abs(rho - rho0) / rho0)
                maxs_ipe(3) = max(maxs_ipe(3), jxb(ix1,ix2,ix3))
                maxs_ipe(4) = max(maxs_ipe(4), abs(divb(ix1,ix2,ix3)))
             end do
          end do
       end do
    end do

    call MPI_ALLREDUCE(sums_ipe, sums, nsum, MPI_DOUBLE_PRECISION, MPI_SUM,&
       icomm, ierrmpi)
    call MPI_ALLREDUCE(maxs_ipe, maxs, nmax, MPI_DOUBLE_PRECISION, MPI_MAX,&
       icomm, ierrmpi)

    if (mype == 0) then
       if (.not. opened) then
          filename = trim(base_filename) // ".log"
          if (restart_from_file == undefined) then
             open(unit=my_unit, file=filename, status='replace')
             write(my_unit,'(A)') '# it t dt Ekin Emag vmax dB_L2 drho_max ' // &
                'sigmaJ jxb_max divb_mean divb_max Etot'
          else
             open(unit=my_unit, file=filename, status='old', position='append')
          end if
          opened = .true.
       end if
       write(my_unit,'(I8,12ES16.8)') it, global_time, dt, sums(1), sums(2), maxs(1),&
          sqrt(sums(3) / sums(4)), maxs(2), sums(5) / max(sums(6), smalldouble),&
          maxs(3), sums(7) / sums(8), maxs(4), sums(9)
       flush(my_unit)
    end if

  end subroutine gh_print_log

  !> div B, |j x B| and |j x B| / (|j| |B|)
  subroutine gh_aux_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw+nwauxio)
    double precision             :: normconv(0:nw+nwauxio)
    double precision :: jxb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: jb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call gh_lorentz(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw),x,jxb,jb,divb)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1) = &
       divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+2) = &
       jxb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+3) = &
       jxb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) / &
       max(jb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3), smalldouble)

  end subroutine gh_aux_output

  subroutine gh_aux_names(varnames)
    character(len=*) :: varnames
    varnames = 'divb jxb sinjb'
  end subroutine gh_aux_names

end module mod_usr
