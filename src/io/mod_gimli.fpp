!> Legolas eigenmodes exported by Pylbo's GIMLI (.ldat files), added to the
!> primitive variables of an initial condition. Host only. Same file layout,
!> routine names and conventions as MPI-AMRVAC's mod_gimli; the one addition is
!> that a cylindrical Legolas mode can be put on a Cartesian mesh or on the
!> (r, z, phi) cylindrical mesh, with the vector components mapped accordingly.
module mod_gimli
#:if PHYS == 'hd' or PHYS == 'mhd'
  use mod_global_parameters
  use mod_geometry, only: coordinate, Cartesian, cylindrical
  use mod_physics_vars
  use mod_comm_lib, only: mpistop

  implicit none
  private

  !> rho, v1, v2, v3, p, B1, B2, B3 in Legolas' (r, theta, z)
  integer, parameter :: n_ef = 8
  integer :: mhd_bool, ef_gridpts
  double precision :: k2, k3
  double precision :: file_units(8)
  double precision, allocatable :: ef_grid(:)
  complex(kind(1.0d0)), allocatable :: ef(:,:)
  logical :: units_checked = .false.

  public :: read_legolas_data, add_perturbation_to_w_array, legolas_perturbation

contains

  subroutine read_legolas_data(legolas_file, file_id)
    character(len=*), intent(in) :: legolas_file
    integer, intent(in) :: file_id
    integer :: i_ef

    if (coordinate /= Cartesian .and. coordinate /= cylindrical) &
       call mpistop("mod_gimli: needs a Cartesian or cylindrical mesh")

    open(unit=file_id+mype, file=trim(legolas_file), form='unformatted', &
       status='old', action='read')
    read(file_id+mype) mhd_bool
    read(file_id+mype) ef_gridpts
    read(file_id+mype) k2, k3
    allocate(ef_grid(ef_gridpts), ef(ef_gridpts, n_ef))
    ef = (0.0d0, 0.0d0)
    read(file_id+mype) ef_grid
    do i_ef = 1, merge(8, 5, mhd_bool == 1)
       read(file_id+mype) ef(:, i_ef)
    end do
    file_units = 0.0d0
    read(file_id+mype) file_units(1:merge(8, 7, mhd_bool == 1))
    close(file_id+mype)

  end subroutine read_legolas_data

  !> The file's units must be the case's: the perturbation is in code units.
  !> Checked at the first use, after the case and phys_units have set them, to 1e-3,
  !> since a case may round a unit the Legolas run derived (n0 from rho0, say).
  subroutine check_units()
    character(len=18), parameter :: names(8) = [character(len=18) :: 'unit_length', &
       'unit_numberdensity', 'unit_temperature', 'unit_density', 'unit_pressure', &
       'unit_velocity', 'unit_magneticfield', 'unit_time']
    double precision :: units(8)
    integer :: i

    units = [unit_length, unit_numberdensity, unit_temperature, unit_density, &
       unit_pressure, unit_velocity, unit_magneticfield, unit_time]
    if (mhd_bool /= 1) then
       units(7) = units(8)
       file_units(7) = file_units(8)
    end if
    do i = 1, merge(8, 7, mhd_bool == 1)
       if (abs(units(i) - file_units(i)) > 1.0d-3 * abs(units(i))) then
          if (mype == 0) write(*,'(A,ES16.8,A,ES16.8)') 'mod_gimli: ' // trim(names(i)) // &
             ' is ', units(i), ' here but ', file_units(i), ' in the Legolas file'
          call mpistop("mod_gimli: the Legolas run and this case use different units")
       end if
    end do
    units_checked = .true.

  end subroutine check_units

  !> w(w_index) = w(w_index) + Re[f(r) exp(i (k2 theta + k3 z))], w primitive;
  !> nothing for a w_index that is not rho, mom, p or mag
  subroutine add_perturbation_to_w_array(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
     ixGmax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,w_index,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    integer, intent(in)             :: w_index
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision                :: pert(n_ef)
    integer                         :: w_of_ef(n_ef), i_ef, ix1, ix2, ix3

    if (.not. units_checked) call check_units()
    w_of_ef = -1
    w_of_ef(1) = rho_
    w_of_ef(2:4) = mom(1:3)
    w_of_ef(5) = p_
#:if PHYS == 'mhd'
    w_of_ef(6:8) = mag(1:3)
#:endif
    i_ef = findloc(w_of_ef, w_index, dim=1)
    if (i_ef == 0) return

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             pert = legolas_perturbation(x(ix1,ix2,ix3,1:ndim))
             w(ix1,ix2,ix3,w_index) = w(ix1,ix2,ix3,w_index) + pert(i_ef)
          end do
       end do
    end do

  end subroutine add_perturbation_to_w_array

  !> rho, v, p, B of the eigenmode at x, vectors in the components of the mesh
  function legolas_perturbation(x) result(pert)
    double precision, intent(in) :: x(1:ndim)
    double precision             :: pert(n_ef)
    double precision             :: r, theta, z

    if (coordinate == Cartesian) then
       r     = hypot(x(1), x(2))
       theta = atan2(x(2), x(1))
       z     = x(3)
    else
       r     = x(r_)
       theta = x(phi_)
       z     = x(z_)
    end if
    pert = dble(ef_amplitude(r) * exp(cmplx(0.0d0, k2 * theta + k3 * z, kind(1.0d0))))
    pert(2:4) = mesh_components(pert(2:4), theta)
    pert(6:8) = mesh_components(pert(6:8), theta)

  end function legolas_perturbation

  !> Eigenfunctions at r, as MPI-AMRVAC: linear between Legolas points, constant
  !> inside the first and beyond the last one within the domain, zero outside it;
  !> on a cylindrical mode the inner extension ramps to zero towards the axis
  function ef_amplitude(r) result(a)
    double precision, intent(in) :: r
    complex(kind(1.0d0))         :: a(n_ef)
    integer                      :: idl, idu

    if (r < xprobmin1 .or. r > xprobmax1) then
       a = (0.0d0, 0.0d0)
    else if (r <= ef_grid(1)) then
       a = ef(1,:)
       if (ef_grid(1) > 0.0d0) a = ef(1,:) * max(2.0d0 * r / ef_grid(1) - 1.0d0, 0.0d0)
    else if (r >= ef_grid(ef_gridpts)) then
       a = ef(ef_gridpts,:)
    else
       idl = maxloc(ef_grid, mask=(ef_grid < r), dim=1)
       idu = minloc(ef_grid, mask=(ef_grid > r), dim=1)
       a = ef(idl,:) + (r - ef_grid(idl)) * (ef(idu,:) - ef(idl,:)) / (ef_grid(idu) - ef_grid(idl))
    end if

  end function ef_amplitude

  !> (r, theta, z) components of a vector as components of the mesh
  pure function mesh_components(v, theta) result(u)
    double precision, intent(in) :: v(3), theta
    double precision             :: u(3)

    if (coordinate == Cartesian) then
       u(1) = cos(theta) * v(1) - sin(theta) * v(2)
       u(2) = sin(theta) * v(1) + cos(theta) * v(2)
       u(3) = v(3)
    else
       u(r_)   = v(1)
       u(phi_) = v(2)
       u(z_)   = v(3)
    end if

  end function mesh_components
#:endif
end module mod_gimli
