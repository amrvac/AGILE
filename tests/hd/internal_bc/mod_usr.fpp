#:mute
#:include "mod_gpu_directives.fpp"
#:endmute

!> Internal boundary conditions (usr_internal_bc), two runs from one build.
!>
!>   fixed.par   setup = 'fixed'   a solid body at rest in a uniform wind
!>   orbit.par   setup = 'orbit'   a body on a circular orbit through gas at
!>                                 rest, so the region moves with time
!>
!> Inside the body of radius rbody the state is held at (rho_b, v_body, p_b)
!> by usr_internal_bc, which getbc applies on the device to every interior
!> cell before each ghost-cell exchange. Both runs use Lohner AMR, so the
!> region is also prolonged and coarsened.
!>
!> The regression log barely notices whether the overwrite reaches the
!> solution: the body is a small fraction of the domain. The real test is
!> check_body, registered as usr_process_grid, which fetches the solution from
!> the device at the start of every step and asserts that every cell well
!> inside the body holds the prescribed state exactly. Before the hook ran on
!> the device, a GPU build wrote the host copy only and this check failed at
!> the first step.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> 'fixed' or 'orbit', from &usr_list
  character(len=20) :: setup = 'fixed'

  !> ambient gas
  double precision :: rho0 = 1.0d0
  double precision :: p0   = 1.0d0
  double precision :: v0(3) = [1.0d0, 0.0d0, 0.0d0]
  !> the body
  double precision :: rho_b = 10.0d0
  double precision :: p_b   = 1.0d0
  double precision :: rbody = 0.2d0
  !> 'orbit' only: orbital radius and angular frequency, about the origin in
  !> the x-y plane
  double precision :: rorbit = 0.4d0
  double precision :: omega  = 2.5d0
  !> the device-side form of `setup`
  logical :: is_orbit = .false.
  ${GPU_DECLARE_COPYIN('rho0, p0, v0, rho_b, p_b, rbody, rorbit, omega, is_orbit')}$

contains

  subroutine usr_init()
    use mod_comm_lib, only: mpistop

    call params_read_user(par_files)
    select case (trim(setup))
    case ('fixed')
       is_orbit = .false.
       v0 = [1.0d0, 0.0d0, 0.0d0]
    case ('orbit')
       is_orbit = .true.
       v0 = 0.0d0
    case default
       call mpistop("&usr_list setup must be 'fixed' or 'orbit'")
    end select
    ${GPU_UPDATE_DEVICE('is_orbit, v0')}$

    usr_init_one_grid => initonegrid_usr
    usr_process_grid  => check_body

    call phys_activate()

  end subroutine usr_init

  subroutine params_read_user(files)
    use mod_global_parameters, only: unitpar
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/ setup

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read_user

  !> Centre and velocity of the body at time qt.
  pure subroutine body_motion(qt, xc, vb)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: qt
    double precision, intent(out) :: xc(3), vb(3)

    if (is_orbit) then
       xc = [rorbit*cos(omega*qt), rorbit*sin(omega*qt), 0.0d0]
       vb = [-omega*xc(2), omega*xc(1), 0.0d0]
    else
       xc = 0.0d0
       vb = 0.0d0
    end if

  end subroutine body_motion

  !> Distance from x to the body's centre at qt.
  pure double precision function body_distance(qt, x)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in) :: qt, x(1:ndim)
    double precision :: xc(3), vb(3)

    call body_motion(qt, xc, vb)
    body_distance = sqrt(sum((x(1:3) - xc)**2))

  end function body_distance

  !> The prescribed conserved state inside the body at qt.
  pure subroutine body_state(qt, w)
    ${GPU_ROUTINE_SEQ()}$
    double precision, intent(in)  :: qt
    double precision, intent(out) :: w(1:nw_phys)
    double precision :: xc(3), vb(3)

    call body_motion(qt, xc, vb)
    w = 0.0d0
    w(iw_rho)       = rho_b
    w(iw_mom(1:3))  = rho_b * vb
    w(iw_e)         = p_b / (hd_gamma - 1.0d0) + 0.5d0 * rho_b * sum(vb**2)

  end subroutine body_state

  !> Internal boundary: hold the body's state inside it. Called by name from
  !> getintbc, on the device, for every interior cell.
  subroutine usr_internal_bc(level, qt, x, w)
    ${GPU_ROUTINE_SEQ()}$
    integer, intent(in)             :: level
    double precision, intent(in)    :: qt
    double precision, intent(in)    :: x(1:ndim)
    double precision, intent(inout) :: w(1:nw_phys)

    if (body_distance(qt, x) <= rbody) call body_state(qt, w)

  end subroutine usr_internal_bc

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision :: wb(1:nw_phys), x_loc(1:ndim)
    integer :: ix1, ix2, ix3

    call body_state(0.0d0, wb)
    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             x_loc = x(ix1,ix2,ix3,1:ndim)
             if (body_distance(0.0d0, x_loc) <= rbody) then
                w(ix1,ix2,ix3,1:nw_phys) = wb
             else
                w(ix1,ix2,ix3,iw_rho)      = rho0
                w(ix1,ix2,ix3,iw_mom(1:3)) = rho0 * v0
                w(ix1,ix2,ix3,iw_e)        = p0 / (hd_gamma - 1.0d0) + &
                     0.5d0 * rho0 * sum(v0**2)
             end if
          end do
       end do
    end do

  end subroutine initonegrid_usr

  !> Assert that every cell well inside the body holds the body's state.
  !> Cells within a small margin of the surface are skipped, because the
  !> time the hook was last applied at and the qt passed here may differ by
  !> round-off, which for the orbiting body moves the surface by as much.
  subroutine check_body(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, parameter :: margin = 1.0d-6, tol = 1.0d-12
    double precision :: wb(1:nw_phys), x_loc(1:ndim), err
    integer :: ix1, ix2, ix3

    ! the solution lives on the device; the host copy is only refreshed for
    ! output, so fetch this block before reading it
    ${GPU_UPDATE_HOST('ps(igrid)%w')}$

    call body_state(qt, wb)
    do ix3 = ixOmin3, ixOmax3
       do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
             x_loc = x(ix1,ix2,ix3,1:ndim)
             if (body_distance(qt, x_loc) > rbody - margin) cycle
             err = maxval(abs(w(ix1,ix2,ix3,1:nw_phys) - wb))
             if (err > tol * maxval(abs(wb))) then
                write(*,*) 'check_body: it =', it, ' igrid =', igrid, &
                     ' x =', x_loc, ' err =', err
                call mpistop('internal boundary check failed')
             end if
          end do
       end do
    end do

  end subroutine check_body

end module mod_usr
