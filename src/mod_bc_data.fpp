#:mute
#:include "mod_gpu_directives.fpp"
#:endmute

!> Boundary conditions from user data ('bc_data' and 'bc_icarus').
!>
!> The VTK file is read once into a single 3D lookup table whose axes are the
!> face's two tangential coordinates and time, with one column per
!> data-driven (variable, face) pair, numbered by bc_data_ix. The table is
!> copied to the device once, after which it never changes, and bc_data_set
!> reads it through the same lookup routines on host and device.
module mod_bc_data
  use mod_lookup_table
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
  type(LT3_t) :: bc_table
  ${GPU_DECLARE_CREATE('bc_table')}$

  !> Column of bc_table for each variable and face, -1 where not data driven
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
    double precision       :: xmin(3), xmax(3)
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

    xmin = bc%origin
    xmax = bc%origin + (bc%n_points-1) * bc%dx
    if (bc%n_points(3) > 1) then
       bc_table = LT3_create_from_data(xmin, xmax, bc%values)
    else
       ! Static data: two identical time levels, so time is clamped away
       xmax(3) = xmin(3) + 1.d0
       bc_table = LT3_create_from_data(xmin, xmax, &
            reshape([bc%values, bc%values], [bc%n_points(1), bc%n_points(2), &
            2, bc%n_variables], order=[1, 2, 4, 3]))
    end if

    ! Shallow copy first, then attach the table itself (as for bgeo)
    ${GPU_UPDATE_DEVICE('bc_table')}$
    ${GPU_ENTER_DATA_COPYIN('bc_table%rows_cols')}$
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

    val = LT3_get_col(bc_table, n_bc, x1, x2, qt)
  end function bc_data_get_3d

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
    type(LT3_loc_t) :: loc
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

    ${GPU_LOOP_VECTOR('collapse(3) private(u, u_in, val, loc, iw, n_bc, jx1, jx2, jx3, kx1, kx2, kx3)')}$
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
             loc = LT3_get_loc(bc_table, x(jx1,jx2,jx3,tangent1), &
                  x(jx1,jx2,jx3,tangent2), qt)

             if (boundary_data_primitive) then
                ! bc_data_init made sure every physical variable has data here
                if (interp_phy_first_row) then
                   u_in = w(kx1,kx2,kx3,1:nw_phys)
                   call to_primitive(u_in)
                end if
                do iw = 1, nw_phys
                   u(iw) = LT3_get_col_at_loc(bc_table, bc_data_ix(iw,iB), loc)
                   if (interp_phy_first_row) u(iw) = 2 * u(iw) - u_in(iw)
                end do
                call to_conservative(u)
                w(ix1,ix2,ix3,1:nw_phys) = u
             end if
             ! Conserved data, and any variables beyond the physical ones
             do iw = merge(nw_phys + 1, 1, boundary_data_primitive), nwfluxbc
                n_bc = bc_data_ix(iw,iB)
                if (n_bc < 1) cycle
                val = LT3_get_col_at_loc(bc_table, n_bc, loc)
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
