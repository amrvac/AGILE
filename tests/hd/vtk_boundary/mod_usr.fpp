#:mute
#:include "mod_gpu_directives.fpp"
#:endmute

!> Data-driven boundaries (typeboundary 'bc_data') in hd, two runs from one
!> build:
!>
!>   agile.par      setup = 'mixed'      face min1: only the density comes from
!>                                       time-varying data (mixed.vtk), with
!>                                       the default interp_phy_first_row = T;
!>                                       momentum and energy are 'cont'
!>   primitive.par  setup = 'primitive'  face min3: all five variables come
!>                                       from static primitive data
!>                                       (primitive.vtk), written directly
!>
!> The data files are made by make_vtk.py. The real test is check_face,
!> registered as usr_process_grid: it fetches each block touching the data
!> face from the device and compares its ghost cells with what the table
!> prescribes, converting primitive data to conserved by hand rather than
!> through the physics module. The rho/vtk_boundary test covers the other row
!> option (bc_phy_first_row), bc_icarus and AMR.
module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  !> 'mixed' or 'primitive', from &usr_list
  character(len=20) :: setup = 'mixed'
  !> initial velocity, towards the interior from the data face
  double precision :: v0(3)

contains

  subroutine usr_init()
    use mod_comm_lib, only: mpistop

    call params_read_user(par_files)
    select case (trim(setup))
    case ('mixed')
       v0 = [0.5d0, 0.0d0, 0.0d0]
    case ('primitive')
       v0 = [0.0d0, 0.0d0, 0.5d0]
    case default
       call mpistop("&usr_list setup must be 'mixed' or 'primitive'")
    end select

    usr_init_one_grid => initonegrid_usr
    usr_process_grid  => check_face

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

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    integer :: idir

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,iw_rho) = 1.0d0
    do idir = 1, 3
       w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,iw_mom(idir)) = v0(idir)
    end do
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,iw_e) = &
         1.0d0 / (hd_gamma - 1.0d0) + 0.5d0 * sum(v0**2)
  end subroutine initonegrid_usr

  !> Compare the ghost cells of the data face with the table: data-driven
  !> variables against the data (2*data - first interior cell with
  !> interp_phy_first_row), 'cont' ones against the first interior cell.
  subroutine check_face(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    use mod_bc_data, only: bc_data_ix, bc_data_get_3d, interp_phy_first_row, &
         bc_phy_first_row, boundary_data_primitive
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, parameter :: tol = 1.0d-12
    double precision :: dat(nw_phys), cons(nw_phys), expected, xmin(3), xmax(3)
    integer :: face, idim, t1, t2, a, b, n, iw, nfirst, nlast, pos, inner
    integer :: lo(3), hi(3), ilo(3), ihi(3), c(3), cp(3), ci(3)

    face = 0
    do n = 1, 2*ndim
       if (any(bc_data_ix(:,n) > 0)) face = n
    end do
    if (face == 0) call mpistop('No data face installed')
    if (bc_phy_first_row) call mpistop('check_face does not cover bc_phy_first_row')

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
    else
       c(idim) = ihi(idim)
       if (x(c(1),c(2),c(3),idim) < xmax(idim)) return
       nfirst = hi(idim) + 1; nlast = ihi(idim); pos = hi(idim) + 1; inner = hi(idim)
    end if

    ! the solution lives on the device
    ${GPU_UPDATE_HOST('ps(igrid)%w')}$

    do b = lo(t2), hi(t2)
       do a = lo(t1), hi(t1)
          cp(t1) = a; cp(t2) = b; cp(idim) = pos
          ci = cp; ci(idim) = inner
          dat = 0.0d0
          do iw = 1, nw_phys
             if (bc_data_ix(iw,face) > 0) dat(iw) = bc_data_get_3d(bc_data_ix(iw,face), &
                  x(cp(1),cp(2),cp(3),t1), x(cp(1),cp(2),cp(3),t2), qt)
          end do
          if (boundary_data_primitive) then
             cons(iw_rho) = dat(iw_rho)
             cons(iw_mom(1:3)) = dat(iw_rho) * dat(iw_mom(1:3))
             cons(iw_e) = dat(iw_e) / (hd_gamma - 1.0d0) + &
                  0.5d0 * dat(iw_rho) * sum(dat(iw_mom(1:3))**2)
          end if
          do n = nfirst, nlast
             c = cp; c(idim) = n
             do iw = 1, nw_phys
                if (boundary_data_primitive) then
                   expected = cons(iw)
                else if (bc_data_ix(iw,face) < 1) then
                   expected = w(ci(1),ci(2),ci(3),iw)     ! 'cont'
                else if (interp_phy_first_row) then
                   expected = 2 * dat(iw) - w(ci(1),ci(2),ci(3),iw)
                else
                   expected = dat(iw)
                end if
                if (abs(w(c(1),c(2),c(3),iw) - expected) > tol * max(1.0d0, abs(expected))) then
                   write(*,*) 'face', face, 'cell', c, 'iw', iw, 'w', &
                        w(c(1),c(2),c(3),iw), 'expected', expected
                   call mpistop('ghost cell differs from the boundary data')
                end if
             end do
          end do
       end do
    end do
  end subroutine check_face

end module mod_usr
