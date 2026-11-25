module mod_input_output_helper

  implicit none

  ! public methods
  public :: count_ix
  public :: create_output_file
  public :: datfile_write_header1
  public :: block_shape_io
  public :: get_names_from_string

  !> whether a manually inserted datfile is saved
  logical :: save_now
  !> Version number of the .dat file output
  integer, parameter :: version_number = 5
  !> The code that wrote a file, from version 6 on: 1 for AGILE
  integer, parameter :: origin_agile = 1
  contains

   function get_names_from_string(aux_variable_names,nwc) result(names) 
    use mod_global_parameters
    character(len=*),intent(in):: aux_variable_names
    integer, intent(in) :: nwc
    character(len=name_len)   :: names(1:nwc)

    integer::  space_position,iw
    character(len=name_len)::  wname
    character(len=std_len)::  scanstring

    ! copied from subroutine getheadernames in calculate_xw
    !TODO check for errors 
    scanstring=TRIM(adjustl(aux_variable_names))
    space_position=index(scanstring,' ')
    do iw=1,nwc
       do while (space_position==1)
         scanstring=scanstring(2:)
         space_position=index(scanstring,' ')
       enddo
       wname=scanstring(:space_position-1)
       scanstring=scanstring(space_position+1:)
       space_position=index(scanstring,' ')

       names(iw)=TRIM(adjustl(wname))
    enddo
  end function get_names_from_string  

  !> Compute number of elements in index range
  pure integer function count_ix(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3)
    integer, intent(in) :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3

    count_ix = product([ ixOmax1,ixOmax2,ixOmax3 ] - [ ixOmin1,ixOmin2,&
       ixOmin3 ] + 1)
  end function count_ix

  !> Determine the shape of a block for output (whether to include ghost cells,
  !> and on which sides)
  subroutine block_shape_io(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
     ixOmax2,ixOmax3, n_values, nghost_uniform)
    use mod_global_parameters

    integer, intent(in) :: igrid
    !> nghost(1:ndim) contains the number of ghost cells on the block's minimum
    !> boundaries, and nghost(ndim+1:2*ndim) on the block's maximum boundaries
    integer, intent(out) :: n_ghost(2*ndim)
    !> Index range on output block
    integer, intent(out) :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    !> Number of cells/values in output
    integer, intent(out) :: n_values
    !> If present (snapshot mode), include this many ghost-cell layers on every
    !> block face (capped at nghostcells)
    integer, intent(in), optional :: nghost_uniform

    integer            :: idim, ng

    n_ghost(:) = 0

    if(present(nghost_uniform)) then
      ng = max(0, min(nghost_uniform, nghostcells))
      n_ghost(:) = ng
    else if(save_physical_boundary) then
      do idim=1,ndim
        ! Include ghost cells on lower boundary
        if(ps(igrid)%is_physical_boundary(2*idim-1)) n_ghost(idim)=nghostcells
        ! Include ghost cells on upper boundary
        if(ps(igrid)%is_physical_boundary(2*idim)) &
           n_ghost(ndim+idim)=nghostcells
      end do
    end if

    ixOmin1 = ixMlo1 - n_ghost(1)
    ixOmin2 = ixMlo2 - n_ghost(2)
    ixOmin3 = ixMlo3 - n_ghost(3)
    ixOmax1 = ixMhi1 + n_ghost(ndim+1)
    ixOmax2 = ixMhi2 + n_ghost(ndim+2)
    ixOmax3 = ixMhi3 + n_ghost(ndim+3)

    n_values = count_ix(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) * nw

  end subroutine block_shape_io

  !> Standard method for creating a new output file
  subroutine create_output_file(fh, ix, extension, suffix)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    integer, intent(out)         :: fh !< File handle
    integer, intent(in)          :: ix !< Index of file
    character(len=*), intent(in) :: extension !< Extension of file
    character(len=*), intent(in), optional :: suffix !< Optional suffix
    character(len=std_len)       :: filename
    integer :: amode

    if (ix >= 10000) then
      call mpistop("Number of output files is limited to 10000 (0...9999)")
    end if

    if (present(suffix)) then
       write(filename,"(a,a,i4.4,a)") trim(base_filename), trim(suffix), ix,&
           extension
    else
       write(filename,"(a,i4.4,a)") trim(base_filename), ix, extension
    end if

    ! MPI cannot easily replace existing files
    open(unit=unitdatfile,file=filename,status='replace')
    close(unitdatfile, status='delete')

    amode = ior(MPI_MODE_CREATE, MPI_MODE_WRONLY)
    call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode, MPI_INFO_NULL, fh,&
        ierrmpi)

    if (ierrmpi /= 0) then
      print *, "Error, cannot create file ", trim(filename)
      call mpistop("Fatal error")
    end if

  end subroutine create_output_file

  !> Write header for a datfile, generalize cons_wnames and nw
  !>
  !> If you edit the header, don't forget to update: datfile_write_header(),
  !> datfile_read_header(), doc/fileformat.md, and
  !> tools/python/amrvac_pytools/datfiles/reading/datfile_utilities.py
  subroutine datfile_write_header1(fh, offset_tree, offset_block,&
      dataset_names, nw_vars, version, is_snap)
    use mod_forest
    use mod_physics
    use mod_global_parameters
    integer, intent(in)                       :: fh           !< File handle
    integer(kind=MPI_OFFSET_KIND), intent(in) :: offset_tree !< Offset of tree info
    integer(kind=MPI_OFFSET_KIND), intent(in) :: offset_block !< Offset of block data
    character(len=*), intent(in) :: dataset_names(:)
    integer, intent(in) :: nw_vars
    !> File format version to write (defaults to the v5 format)
    integer, intent(in), optional :: version
    !> An analysis snapshot: its first counter is the next snapshot index
    logical, intent(in), optional :: is_snap
    integer, dimension(MPI_STATUS_SIZE)       :: st
    integer                                   :: iw, er

    character(len=name_len) :: dname
    integer                 :: file_version
    logical                 :: snap_mode

    file_version = version_number
    if (present(version)) file_version = version
    snap_mode = .false.
    if (present(is_snap)) snap_mode = is_snap

    call MPI_FILE_WRITE(fh, file_version, 1, MPI_INTEGER, st, er)
    if (file_version > 5)&
      call MPI_FILE_WRITE(fh, origin_agile, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, int(offset_tree), 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, int(offset_block), 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, nw_vars, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, ndir, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, ndim, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, levmax, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, nleafs, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, nparents, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, it, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, global_time, 1, MPI_DOUBLE_PRECISION, st, er)

    call MPI_FILE_WRITE(fh, [ xprobmin1,xprobmin2,xprobmin3 ], ndim,&
        MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, [ xprobmax1,xprobmax2,xprobmax3 ], ndim,&
        MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, [ domain_nx1,domain_nx2,domain_nx3 ], ndim,&
        MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, [ block_nx1,block_nx2,block_nx3 ], ndim,&
        MPI_INTEGER, st, er)

    ! Periodicity (assume all variables are periodic if one is)
    call MPI_FILE_WRITE(fh, periodB, ndim, MPI_LOGICAL, st, er)

    ! Geometry
    call MPI_FILE_WRITE(fh, geometry_name(1:name_len), name_len, MPI_CHARACTER,&
        st, er)

    ! Write stagger grid mark
    call MPI_FILE_WRITE(fh, stagger_grid, 1, MPI_LOGICAL, st, er)

    ! From version 6: bytes per real in the block data (self-describing
    ! precision, so readers need not map version number to precision)
    if (file_version > 5) then
      call MPI_FILE_WRITE(fh, snap_size_real, 1, MPI_INTEGER, st, er)

      ! Compression scheme name, 'none' if uncompressed.
      ! Scheme-specific metadata follows the name.
      if (snap_compress_zfp) then
        dname = 'zfp'
      else
        dname = 'none'
      end if
      call MPI_FILE_WRITE(fh, dname, name_len, MPI_CHARACTER, st, er)
      if (snap_compress_zfp) then
        call MPI_FILE_WRITE(fh, snap_zfp_rel_tol, 1, MPI_DOUBLE_PRECISION,&
            st, er)
        call MPI_FILE_WRITE(fh, snap_zfp_rel_eps, 1, MPI_DOUBLE_PRECISION,&
            st, er)
      end if
    end if

    do iw = 1, nw_vars
      ! using directly trim(adjustl((dataset_names(iw)))) in MPI_FILE_WRITE call 
      ! does not work, there will be trailing characters
      dname = trim(adjustl((dataset_names(iw))))
      call MPI_FILE_WRITE(fh, dname, name_len, MPI_CHARACTER, st, er)
    end do

    ! Physics related information
    call MPI_FILE_WRITE(fh, physics_type, name_len, MPI_CHARACTER, st, er)

    ! Format:
    ! integer :: n_par
    ! double precision :: values(n_par)
    ! character(n_par * name_len) :: names
    call phys_write_info(fh)

    ! The output counters, useful for restarting.
    ! datfilenext is the next file index in the series, which can either be the
    ! series of resume snapshots or analysis snapshots.
    ! Counters for snapshots are advanced *after* this procedure, hence the +1.
    if (snap_mode) then
      call MPI_FILE_WRITE(fh, snapnext+1, 1, MPI_INTEGER, st, er)
    else if (pass_wall_time.or.save_now) then
      call MPI_FILE_WRITE(fh, datfilenext, 1, MPI_INTEGER, st, er)
    else
      call MPI_FILE_WRITE(fh, datfilenext+1, 1, MPI_INTEGER, st, er)
    end if
    call MPI_FILE_WRITE(fh, slicenext, 1, MPI_INTEGER, st, er)
    call MPI_FILE_WRITE(fh, collapsenext, 1, MPI_INTEGER, st, er)
    ! The analysis snapshot counter, to point a resume snapshot to the next
    ! analysis snapshot. Has no meaning for analysis snapshots themselves.
    if (file_version > 5) then
      call MPI_FILE_WRITE(fh, merge(-1, snapnext, snap_mode), 1, MPI_INTEGER,&
          st, er)
    end if

  end subroutine datfile_write_header1

end module mod_input_output_helper
