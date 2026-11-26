#:mute
#:include "mod_gpu_directives.fpp"
#:endmute

!> Boundary conditions from user data ('bc_data' and 'bc_icarus').
!>
!> The VTK file is read once into a single table, bc_values, whose axes are
!> the face's two tangential coordinates and time, with one column per
!> data-driven (variable, face) pair, numbered by bc_data_ix. It is copied to
!> the device once, after which it never changes, and host and device read it
!> through the same routines (bc_locate, bc_interp). It is kept in plain
!> module arrays rather than a derived type with an allocatable component:
!> nvhpc 25.11 faulted reading the latter on the device, where 26.9 did not.
module mod_bc_data
  use mod_global_parameters, only: std_len

  implicit none
  private

  integer, parameter :: dp = kind(1.0d0)

  type bc_data_t
     integer                        :: n_variables
     double precision               :: origin(3)
     double precision               :: dx(3)
     integer                        :: n_points(3)
     character(len=40), allocatable :: names(:)
     double precision, allocatable  :: values(:, :, :, :)
  end type bc_data_t

  !> All boundary data: axes (tangent1, tangent2, time), one column per
  !> data-driven variable and face
  real(dp), allocatable :: bc_values(:, :, :, :)
  !> Table origin, inverse spacing and size along (tangent1, tangent2, time)
  real(dp) :: bc_xmin(3), bc_inv_dx(3)
  integer  :: bc_npts(3)
  ${GPU_DECLARE_CREATE('bc_values, bc_xmin, bc_inv_dx, bc_npts')}$

  !> Column of bc_values for each variable and face, -1 where not data driven
  integer, public, protected, allocatable :: bc_data_ix(:, :)
  ${GPU_DECLARE_CREATE('bc_data_ix')}$

  !> data file name
  character(len=std_len), public, protected :: boundary_data_file_name = ''
  !> Set the data at the physical boundary (the face) instead of the ghost
  !> cell centres: the ghost cells get 2*data - (first interior cell)
  logical, public, protected :: interp_phy_first_row=.true.
  !> Also overwrite the first interior row with the data
  logical, public, protected :: bc_phy_first_row=.false.
  !> The data are primitive variables (every physical variable of a data
  !> face must then be data driven)
  logical, public, protected :: boundary_data_primitive=.false.
  ${GPU_DECLARE_CREATE('interp_phy_first_row, bc_phy_first_row, boundary_data_primitive')}$

  public :: bc_data_init
  public :: bc_data_set
  public :: bc_data_get_3d

contains

  subroutine bc_data_init()
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    use mod_physics, only: nw_phys

    integer                :: iw, ib, n_bc
    type(bc_data_t)        :: bc

    call bc_read_params(par_files)

    allocate(bc_data_ix(nwfluxbc, 2*ndim))
    bc_data_ix(:, :) = -1

    if (.not. any(typeboundary(1:nwfluxbc,:)==bc_data .or. &
                  typeboundary(1:nwfluxbc,:)==bc_icarus)) return
#:if not defined('SPECIALBOUNDARY')
    if (any(typeboundary(1:nwfluxbc,:)==bc_icarus)) &
         call mpistop("bc_icarus requires a build with specialboundary=.true. and specialbound_usr")
#:endif
    if (len_trim(boundary_data_file_name) == 0) call mpistop("bc_data requires boundary_data_file_name")
    if (interp_phy_first_row .and. bc_phy_first_row) &
         call mpistop("bc_data: interp_phy_first_row and bc_phy_first_row exclude each other")
    call read_vtk_structured_points(trim(boundary_data_file_name), bc)
    if (any(bc%n_points(1:2) < 2) .or. bc%n_points(3) < 1 .or. any(bc%dx <= 0.d0)) &
         call mpistop("bc_data: invalid VTK dimensions or spacing")
    if (count(typeboundary(1:nwfluxbc,:)==bc_data .or. &
              typeboundary(1:nwfluxbc,:)==bc_icarus) /= bc%n_variables) &
         call mpistop("bc_data: scalar arrays must follow boundary then variable order")

    n_bc = 0
    do ib = 1, 2 * ndim
       do iw = 1, nwfluxbc
          if (typeboundary(iw, ib)==bc_data .or. typeboundary(iw, ib)==bc_icarus) then
             n_bc               = n_bc + 1
             bc_data_ix(iw, ib) = n_bc
          end if
       end do
       ! Primitive data are converted per cell, which needs the whole state
       if (boundary_data_primitive .and. any(bc_data_ix(:, ib) > 0) .and. &
           any(bc_data_ix(1:nw_phys, ib) < 0)) &
            call mpistop("bc_data: primitive data need every physical variable of a data face")
    end do

    bc_xmin   = bc%origin
    bc_inv_dx = 1 / bc%dx
    bc_npts   = bc%n_points
    if (bc%n_points(3) > 1) then
       allocate(bc_values, source=bc%values)
    else
       ! Static data: two identical time levels, so time is clamped away
       bc_npts(3)   = 2
       bc_inv_dx(3) = 1
       allocate(bc_values(bc_npts(1), bc_npts(2), 2, bc%n_variables))
       bc_values(:, :, 1, :) = bc%values(:, :, 1, :)
       bc_values(:, :, 2, :) = bc%values(:, :, 1, :)
    end if

    ${GPU_UPDATE_DEVICE('bc_values, bc_xmin, bc_inv_dx, bc_npts')}$
    ${GPU_UPDATE_DEVICE('bc_data_ix, interp_phy_first_row, bc_phy_first_row, boundary_data_primitive')}$
  end subroutine bc_data_init

  !> Read this module"s parameters from a file
  subroutine bc_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /bd_list/ boundary_data_file_name, interp_phy_first_row,&
        bc_phy_first_row,boundary_data_primitive

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, bd_list, end=111)
111    close(unitpar)
    end do

  end subroutine bc_read_params

  !> Boundary value of table column n_bc at tangential coordinates (x1, x2)
  !> and time qt
  pure function bc_data_get_3d(n_bc, x1, x2, qt) result(val)
    ${GPU_ROUTINE_SEQ()}$
    integer, intent(in)  :: n_bc
    real(dp), intent(in) :: x1, x2, qt
    real(dp)             :: val

    integer  :: ix(3)
    real(dp) :: frac(3)

    call bc_locate(x1, x2, qt, ix, frac)
    val = bc_interp(n_bc, ix, frac)
  end function bc_data_get_3d

  !> Cell of the table containing (x1, x2, qt): the lower corner ix and the
  !> weight frac of that corner along each axis. Coordinates outside the
  !> table are clamped to its edge.
  pure subroutine bc_locate(x1, x2, qt, ix, frac)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: x1, x2, qt
    integer, intent(out)  :: ix(3)
    real(dp), intent(out) :: frac(3)
    real(dp) :: s
    integer  :: i

    do i = 1, 3
       select case (i)
       case (1)
          s = (x1 - bc_xmin(1)) * bc_inv_dx(1)
       case (2)
          s = (x2 - bc_xmin(2)) * bc_inv_dx(2)
       case default
          s = (qt - bc_xmin(3)) * bc_inv_dx(3)
       end select
       ! table point j sits at s = j - 1
       if (s <= 0) then
          ix(i) = 1; frac(i) = 1
       else if (s >= bc_npts(i) - 1) then
          ix(i) = bc_npts(i) - 1; frac(i) = 0
       else
          ix(i) = int(s) + 1; frac(i) = ix(i) - s
       end if
    end do
  end subroutine bc_locate

  !> Trilinear interpolation of column col at a location from bc_locate
  pure function bc_interp(col, ix, frac) result(val)
    ${GPU_ROUTINE_SEQ()}$
    integer, intent(in)  :: col, ix(3)
    real(dp), intent(in) :: frac(3)
    real(dp) :: val, w1, w2, w3
    integer  :: i, j, k

    val = 0.0_dp
    do k = 0, 1
       w3 = merge(frac(3), 1 - frac(3), k == 0)
       do j = 0, 1
          w2 = merge(frac(2), 1 - frac(2), j == 0)
          do i = 0, 1
             w1 = merge(frac(1), 1 - frac(1), i == 0)
             val = val + w1 * w2 * w3 * bc_values(ix(1)+i, ix(2)+j, ix(3)+k, col)
          end do
       end do
    end do
  end function bc_interp

  !> Fill the ghost cells ixO of face iB from the boundary data, after the
  !> ordinary boundary types. Runs per cell, on the host or inside bc_phys on
  !> the device. Variables without data on this face are left alone.
  !>
  !> The data are looked up at the tangential position of the ghost row
  !> ixOmin along the face normal. With interp_phy_first_row the data are the
  !> face value, and every ghost cell gets 2*data - (first interior cell);
  !> with bc_phy_first_row the first interior row is overwritten as well.
  subroutine bc_data_set(qt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB,w,x)
    use mod_global_parameters
    use mod_physics, only: nw_phys, to_primitive, to_conservative
    ${GPU_ROUTINE_VECTOR()}$

    integer, intent(in)     :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB
    real(dp), intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    real(dp), intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    real(dp)        :: u(nw_phys), u_in(nw_phys), val
    real(dp)        :: frac(3)
    integer         :: loc(3)
    integer         :: idim, tangent1, tangent2, pos_row, inner_row, iw, n_bc
    integer         :: lo1, lo2, lo3, hi1, hi2, hi3, ix1, ix2, ix3
    integer         :: jx1, jx2, jx3, kx1, kx2, kx3

    idim = (iB + 1) / 2
    select case (idim)
    case (1)
       tangent1 = 2; tangent2 = 3
    case (2)
       tangent1 = 1; tangent2 = 3
    case default
       tangent1 = 1; tangent2 = 2
    end select

    ! Cells to fill: the ghost layer, plus the first interior row on request
    lo1 = ixOmin1; lo2 = ixOmin2; lo3 = ixOmin3
    hi1 = ixOmax1; hi2 = ixOmax2; hi3 = ixOmax3
    select case (iB)
    case (1)
       pos_row = ixOmin1; inner_row = ixOmax1 + 1
       if (bc_phy_first_row) hi1 = hi1 + 1
    case (2)
       pos_row = ixOmin1; inner_row = ixOmin1 - 1
       if (bc_phy_first_row) lo1 = lo1 - 1
    case (3)
       pos_row = ixOmin2; inner_row = ixOmax2 + 1
       if (bc_phy_first_row) hi2 = hi2 + 1
    case (4)
       pos_row = ixOmin2; inner_row = ixOmin2 - 1
       if (bc_phy_first_row) lo2 = lo2 - 1
    case (5)
       pos_row = ixOmin3; inner_row = ixOmax3 + 1
       if (bc_phy_first_row) hi3 = hi3 + 1
    case default
       pos_row = ixOmin3; inner_row = ixOmin3 - 1
       if (bc_phy_first_row) lo3 = lo3 - 1
    end select

    ${GPU_LOOP_VECTOR('collapse(3) private(u, u_in, val, loc, frac, iw, n_bc, jx1, jx2, jx3, kx1, kx2, kx3)')}$
    do ix3 = lo3, hi3
       do ix2 = lo2, hi2
          do ix1 = lo1, hi1
             ! jx: where the data are looked up; kx: the first interior cell
             jx1 = ix1; jx2 = ix2; jx3 = ix3
             kx1 = ix1; kx2 = ix2; kx3 = ix3
             select case (idim)
             case (1)
                jx1 = pos_row; kx1 = inner_row
             case (2)
                jx2 = pos_row; kx2 = inner_row
             case default
                jx3 = pos_row; kx3 = inner_row
             end select
             call bc_locate(x(jx1,jx2,jx3,tangent1), x(jx1,jx2,jx3,tangent2), &
                  qt, loc, frac)

             if (boundary_data_primitive) then
                ! bc_data_init made sure every physical variable has data here
                if (interp_phy_first_row) then
                   u_in = w(kx1,kx2,kx3,1:nw_phys)
                   call to_primitive(u_in)
                end if
                do iw = 1, nw_phys
                   u(iw) = bc_interp(bc_data_ix(iw,iB), loc, frac)
                   if (interp_phy_first_row) u(iw) = 2 * u(iw) - u_in(iw)
                end do
                call to_conservative(u)
                w(ix1,ix2,ix3,1:nw_phys) = u
             end if
             ! Conserved data, and any variables beyond the physical ones
             do iw = merge(nw_phys + 1, 1, boundary_data_primitive), nwfluxbc
                n_bc = bc_data_ix(iw,iB)
                if (n_bc < 1) cycle
                val = bc_interp(n_bc, loc, frac)
                if (interp_phy_first_row) val = 2 * val - w(kx1,kx2,kx3,iw)
                w(ix1,ix2,ix3,iw) = val
             end do
          end do
       end do
    end do
  end subroutine bc_data_set

  subroutine read_vtk_structured_points(fname, bc)
    character(len=*), intent(in)  :: fname
    type(bc_data_t), intent(out)  :: bc
    double precision, allocatable :: tmp_data(:, :, :, :)
    integer, parameter            :: max_variables = 10
    character(len=40)             :: tmp_names(max_variables)
    character(len=256)            :: line
    character(len=40)             :: word, typename
    integer, parameter            :: my_unit = 123
    integer                       :: n, n_points_total
    integer                       :: n_components

    open(my_unit, file=trim(fname), status="old", action="read")

    ! Header, e.g. # vtk DataFile Version 2.0
    read(my_unit, "(A)") line

    ! Dataset name
    read(my_unit, "(A)") line

    ! ASCII / BINARY
    read(my_unit, "(A)") line

    if (line /= "ASCII") then
       print *, "line: ", trim(line)
       error stop "read_vtk: not an ASCII file"
    end if

    ! DATASET STRUCTURED_POINTS
    read(my_unit, "(A)") line

    if (line /= "DATASET STRUCTURED_POINTS") then
       print *, "line: ", trim(line)
       error stop "read_vtk must have: DATASET STRUCTURED_POINTS"
    end if

    ! DIMENSIONS NX NY NZ
    read(my_unit, "(A)") line
    read(line, *) word, bc%n_points

    if (word /= "DIMENSIONS") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: DIMENSIONS"
    end if

    ! SPACING DX DY DZ
    read(my_unit, *) word, bc%dx

    if (word /= "SPACING") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: SPACING"
    end if

    ! ORIGIN 0 0 0
    read(my_unit, *) word, bc%origin
    if (word /= "ORIGIN") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: ORIGIN"
    end if

    ! POINT_DATA NPOINTS
    read(my_unit, *) word, n_points_total

    if (word /= "POINT_DATA") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: POINT_DATA n_points"
    end if

    if (n_points_total /= product(bc%n_points)) error stop "read_vtk: n_points not equal to NX*NY*NZ"

    allocate(tmp_data(bc%n_points(1), bc%n_points(2), bc%n_points(3),&
        max_variables))

    ! Read all scalar variables
    do n = 1, max_variables

       ! SCALARS name type ncomponents
       read(my_unit, *, end=900) word, tmp_names(n), typename, n_components

       if (word /= "SCALARS") then
          print *, "line: ", trim(line)
          error stop "read_vtk expects: SCALARS name type ncomponents"
       end if

       if (n_components /= 1) error stop "read_vtk: ncomponents should be 1"

       ! LOOKUP_TABLE default
       read(my_unit, *) word, typename

       if (word /= "LOOKUP_TABLE") then
          print *, "line: ", trim(line)
          error stop "read_vtk expects: LOOKUP_TABLE name"
       end if

       ! Read list of data values
       read(my_unit, *) tmp_data(:, :, :, n)
    end do

    ! Done reading variables
900 continue

    close(my_unit)

    if (n == max_variables + 1) error stop "read_vtk: increase max_variables"

    ! Loop index is one higher than number of variables
    bc%n_variables = n-1
    bc%values      = tmp_data(:, :, :, 1:n-1)
    bc%names       = tmp_names(1:n-1)
  end subroutine read_vtk_structured_points

end module mod_bc_data
