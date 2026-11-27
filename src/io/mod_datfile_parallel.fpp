#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute
!> Parallel MPI-IO for restart datfiles (the v5 format).
!>
!> The default path in mod_input_output funnels every block through rank 0,
!> which writes (or reads) it with separate small calls. Here every rank
!> writes and reads its own blocks with MPI-IO instead, producing the same
!> file byte for byte. Selected by datfile_parallel_io in &filelist.
!>
!> A rank's leaves are consecutive in Morton order and the file stores the
!> block records in Morton order, so each rank's records form one contiguous
!> byte range. Every record is a 2*ndim integer ghost-cell header followed by
!> the block's values in double precision. The header is 8*ndim/2 = 24 bytes,
!> exactly 3 doubles, so a rank packs its records into one double precision
!> buffer and moves it in one call. The buffer is capped at
!> datfile_io_buffer_mb, so a rank with more data than that needs several
!> rounds.
!>
!> By default (datfile_io_collective = F) each rank moves its own rounds with
!> independent MPI_FILE_WRITE_AT / MPI_FILE_READ_AT. With
!> datfile_io_collective = T it uses MPI_FILE_WRITE_AT_ALL /
!> MPI_FILE_READ_AT_ALL, and every rank joins every round, with a count of 0
!> once it has nothing left. The data of a rank is one contiguous range, so
!> collective buffering has nothing to merge: with Open MPI's OMPIO it only
!> adds an exchange through aggregator ranks and a synchronisation per round,
!> which made reads several times slower than the serial path (KH3D 512^3 on
!> 4 A100s: 43 s collective vs 7 s independent at 16 MiB, serial 10 s).
!>
!> The metadata (header, tree, block levels, indices and offsets) is small
!> and stays with rank 0, as in the default path.
!>
!> Analysis snapshots (v6, write_snapshot_parallel) work the same way, with
!> the v6 differences: no per-record header (the ghost-cell counts go in the
!> tree), values in snap_size_real bytes, optionally snap_nghost ghost layers,
!> and optionally ZFP compression. ZFP streams have sizes known only after
!> compressing, so a rank first compresses all its blocks (the same adaptive
!> per-field search, in the same order and with the same precision chain as
!> the default path, hence the same bytes), keeps the streams in memory, and
!> only then learns its offset and writes them.
module mod_datfile_parallel
#ifdef USE_MPIWRAPPERS
  use mod_mpi_wrapper, only: mpi_file_write_at_wrapper,&
     mpi_file_write_at_all_wrapper
#else
#define mpi_file_write_at_wrapper MPI_FILE_WRITE_AT
#define mpi_file_write_at_all_wrapper MPI_FILE_WRITE_AT_ALL
#endif
  use mod_comm_lib, only: mpistop

  implicit none
  private

  !> Doubles in a record's ghost-cell header (2*ndim integers)
  integer, parameter :: nhdr_dbl = 3

  public :: write_datfile_parallel
  public :: write_snapshot_parallel
  public :: read_datfile_blocks_parallel
  public :: store_datfile_block

contains

  !> Write a v5 restart datfile with every rank writing its own blocks.
  subroutine write_datfile_parallel()
    use mod_forest
    use mod_global_parameters
    use mod_variables, only: cons_wnames
    use mod_input_output_helper, only: count_ix, block_shape_io,&
       create_output_file, datfile_write_header1
    use mod_functions_forest, only: write_forest

    integer                       :: fh, igrid, ipe, Morton_no, ib, ib0
    integer                       :: n_ghost(2*ndim), n_values, n_rec
    integer                       :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    integer                       :: nrounds, iround, nfill, m_next
    integer                       :: offset_size, amode
    integer                       :: istatus(MPI_STATUS_SIZE)
    integer                       :: counts(0:npe-1), displs(0:npe-1)
    integer(kind=8)               :: cap_dbl
    integer(kind=MPI_OFFSET_KIND) :: offset_tree_info, offset_offsets,&
       offset_block_data, offset, offset_end, pos
    integer(kind=8), allocatable  :: rec_dbl(:), my_rec_dbl(:)
    integer, allocatable          :: block_ig(:, :), block_lvl(:)
    integer(kind=MPI_OFFSET_KIND), allocatable :: block_offset(:)
    double precision, allocatable :: buf(:)
    type(tree_node), pointer      :: pnode
    character(len=std_len)        :: filename

    call check_header_size()

    ! Size of each local record in doubles
    allocate(my_rec_dbl(max(1, Morton_stop(mype)-Morton_start(mype)+1)))
    do Morton_no = Morton_start(mype), Morton_stop(mype)
      igrid = sfc_to_igrid(Morton_no)
      call record_shape(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
         ixOmax2,ixOmax3, n_values)
      my_rec_dbl(Morton_no-Morton_start(mype)+1) = nhdr_dbl + n_values
    end do

    ! Every rank learns every record size, and with it every file offset
    do ipe = 0, npe-1
      counts(ipe) = Morton_stop(ipe) - Morton_start(ipe) + 1
      displs(ipe) = Morton_start(ipe) - 1
    end do
    allocate(rec_dbl(nleafs))
    call MPI_ALLGATHERV(my_rec_dbl, counts(mype), MPI_INTEGER8, rec_dbl,&
       counts, displs, MPI_INTEGER8, icomm, ierrmpi)

    ! Rank 0 writes the metadata, with the final offsets, and closes the file
    if (mype == 0) then
      call create_output_file(fh, datfilenext, ".dat")

      ! The header holds offsets not known yet; it is written again below
      offset_tree_info = -1
      offset_block_data = -1
      call datfile_write_header1(fh, offset_tree_info, offset_block_data,&
         cons_wnames, nw, 5, .false.)
      call MPI_FILE_GET_POSITION(fh, offset_tree_info, ierrmpi)

      call write_forest(fh)

      allocate(block_ig(ndim, nleafs), block_lvl(nleafs))
      do Morton_no = 1, nleafs
        pnode => igrid_to_node(sfc(1, Morton_no), sfc(2, Morton_no))%node
        block_ig(:, Morton_no) = [ pnode%ig1,pnode%ig2,pnode%ig3 ]
        block_lvl(Morton_no) = pnode%level
      end do
      call MPI_FILE_WRITE(fh, block_lvl, nleafs, MPI_INTEGER, istatus,&
         ierrmpi)
      call MPI_FILE_WRITE(fh, block_ig, size(block_ig), MPI_INTEGER,&
         istatus, ierrmpi)

      call MPI_FILE_GET_POSITION(fh, offset_offsets, ierrmpi)
      call MPI_TYPE_SIZE(MPI_OFFSET, offset_size, ierrmpi)
      offset_block_data = offset_offsets + int(nleafs, MPI_OFFSET_KIND)*&
         offset_size

      allocate(block_offset(nleafs))
      block_offset(1) = offset_block_data
      do Morton_no = 2, nleafs
        block_offset(Morton_no) = block_offset(Morton_no-1) +&
           rec_dbl(Morton_no-1)*size_double
      end do
      call MPI_FILE_WRITE(fh, block_offset, nleafs, MPI_OFFSET, istatus,&
         ierrmpi)

      call MPI_FILE_GET_POSITION(fh, pos, ierrmpi)
      if (pos /= offset_block_data) call mpistop(&
         "write_datfile_parallel: metadata size differs from expected")

      call MPI_FILE_SEEK(fh, 0_MPI_OFFSET_KIND, MPI_SEEK_SET, ierrmpi)
      call datfile_write_header1(fh, offset_tree_info, offset_block_data,&
         cons_wnames, nw, 5, .false.)
      call MPI_FILE_CLOSE(fh, ierrmpi)
    end if

    ! After rank 0 closed the file, so it exists for everyone
    call MPI_BCAST(offset_block_data, 1, MPI_OFFSET, 0, icomm, ierrmpi)

    filename = output_filename(datfilenext, "")
    amode = MPI_MODE_WRONLY
    call MPI_FILE_OPEN(icomm, filename, amode, MPI_INFO_NULL, fh, ierrmpi)
    if (ierrmpi /= MPI_SUCCESS) call mpistop(&
       "write_datfile_parallel: cannot open "//trim(filename))

    offset = offset_block_data
    do Morton_no = 1, Morton_start(mype)-1
      offset = offset + rec_dbl(Morton_no)*size_double
    end do
    offset_end = offset
    do Morton_no = Morton_start(mype), Morton_stop(mype)
      offset_end = offset_end + rec_dbl(Morton_no)*size_double
    end do

    cap_dbl = io_buffer_doubles()
    nrounds = io_rounds(rec_dbl(Morton_start(mype):Morton_stop(mype)),&
       cap_dbl)
    allocate(buf(buffer_length(rec_dbl(Morton_start(mype):&
       Morton_stop(mype)), cap_dbl)))

    m_next = Morton_start(mype)
    do iround = 1, nrounds
      ! Pack as many whole records as fit (always at least one)
      ib = 0
      do while (m_next <= Morton_stop(mype))
        n_rec = int(rec_dbl(m_next))
        if (ib > 0 .and. int(ib, 8) + n_rec > cap_dbl) exit
        igrid = sfc_to_igrid(m_next)
        ib0 = ib
        call pack_record(igrid, buf, ib)
        if (ib - ib0 /= n_rec) call mpistop(&
           "write_datfile_parallel: record size differs from expected")
        m_next = m_next + 1
      end do
      nfill = ib

      if (datfile_io_collective) then
        call MPI_FILE_WRITE_AT_ALL(fh, offset, buf, nfill,&
           MPI_DOUBLE_PRECISION, istatus, ierrmpi)
      else
        call MPI_FILE_WRITE_AT(fh, offset, buf, nfill, MPI_DOUBLE_PRECISION,&
           istatus, ierrmpi)
      end if
      if (ierrmpi /= MPI_SUCCESS) call mpistop(&
         "write_datfile_parallel: writing the block data failed")
      offset = offset + int(nfill, MPI_OFFSET_KIND)*size_double
    end do

    if (m_next /= Morton_stop(mype)+1 .or. offset /= offset_end) call&
       mpistop("write_datfile_parallel: not every block was written")

    call MPI_FILE_CLOSE(fh, ierrmpi)
    call MPI_BARRIER(icomm, ierrmpi)

  end subroutine write_datfile_parallel

  !> Write a v6 analysis snapshot, <base_filename>_snapNNNN.dat, with every
  !> rank writing its own blocks. The file is the same, byte for byte, as the
  !> one write_datfile(is_snap=.true.) writes.
  !>
  !> 1. Each rank sizes its records: from the block shape when uncompressed;
  !>    with ZFP by compressing all its blocks (unchanged adaptive search and
  !>    precision chain) into a host buffer that it keeps.
  !> 2. Rank 0 gathers the record sizes, ghost-cell counts and (ZFP) per-field
  !>    byte counts, writes the metadata and closes the file; each rank's
  !>    start offset is an exclusive scan of the ranks' byte totals.
  !> 3. Every rank writes its contiguous range in rounds of at most
  !>    datfile_io_buffer_mb.
  subroutine write_snapshot_parallel()
    use mod_forest
    use mod_global_parameters
    use mod_variables, only: cons_wnames
    use mod_input_output_helper, only: create_output_file,&
       datfile_write_header1
    use mod_functions_forest, only: write_forest
#:if defined('COMPRESS_ZFP')
    use mod_zfp, only: zfp_max_field_nbytes,&
       zfp_compress_field_adaptive_dp, zfp_compress_field_adaptive_sp
#:endif

    integer                       :: fh, igrid, ipe, Morton_no, iloc, nloc
    integer                       :: n_ghost(2*ndim), n_values
    integer                       :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    integer                       :: nrounds, iround, nfill, ib, ib0, m_next
    integer                       :: offset_size, size_real_out
    integer                       :: istatus(MPI_STATUS_SIZE)
    integer                       :: counts(0:npe-1), displs(0:npe-1)
    integer                       :: counts_k(0:npe-1), displs_k(0:npe-1)
    integer(kind=8)               :: cap, my_bytes, my_start, n_rec, zpos
    integer(kind=MPI_OFFSET_KIND) :: offset_tree_info, offset_offsets,&
       offset_block_data, offset, offset_end, pos
    integer(kind=8), allocatable  :: my_rec(:), my_val(:), rec_bytes(:)
    integer, allocatable          :: my_nghost(:, :), block_nghost(:, :)
    integer, allocatable          :: my_field_nbytes(:, :),&
       zfp_block_nbytes(:, :)
    integer, allocatable          :: block_ig(:, :), block_lvl(:)
    integer(kind=MPI_OFFSET_KIND), allocatable :: block_offset(:)
    double precision, allocatable :: buf(:)
    real(kind=4), allocatable     :: buf_sp(:)
    integer(kind=1), allocatable  :: zbuf(:)
    logical                       :: write_sp, compress
    type(tree_node), pointer      :: pnode
    character(len=std_len)        :: filename
#:if defined('COMPRESS_ZFP')
    integer                       :: iw, nx1, nx2, nx3
    ! ZFP precision of the previous block, per field: the start guess of the
    ! adaptive search, as in write_datfile
    integer                       :: zfp_prec(nw)
    integer(kind=8)               :: zfp_nbytes, field_cap
    integer(kind=1), allocatable  :: zgrow(:)
#:endif

    size_real_out = snap_size_real
    write_sp = (size_real_out == 4)
    compress = .false.
#:if defined('COMPRESS_ZFP')
    compress = snap_compress_zfp
#:endif

    ! 1. Shapes and sizes of this rank's records
    nloc = Morton_stop(mype) - Morton_start(mype) + 1
    allocate(my_rec(max(1, nloc)), my_nghost(2*ndim, max(1, nloc)))
    my_rec = 0
    my_nghost = 0
    do iloc = 1, nloc
      igrid = sfc_to_igrid(Morton_start(mype) + iloc - 1)
      call snap_record_shape(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
         ixOmax2,ixOmax3, n_values)
      my_nghost(:, iloc) = n_ghost
      my_rec(iloc) = int(n_values, 8) * size_real_out
    end do

    zpos = 0
#:if defined('COMPRESS_ZFP')
    if (compress) then
      ! Every block's nw streams fit in field_cap bytes (compression ratio 1
      ! at the maximum precision), as in write_datfile
      field_cap = nw * zfp_max_field_nbytes(block_nx1 + 2 * nghostcells,&
         block_nx2 + 2 * nghostcells, block_nx3 + 2 * nghostcells,&
         .not. write_sp, 64)
      allocate(my_field_nbytes(nw, max(1, nloc)))
      my_field_nbytes = 0
      allocate(zbuf(field_cap * min(max(nloc, 1), 4)))
      zfp_prec = 0
      do iloc = 1, nloc
        igrid = sfc_to_igrid(Morton_start(mype) + iloc - 1)
        call snap_record_shape(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,&
           ixOmax1,ixOmax2,ixOmax3, n_values)
        if (zpos + field_cap > size(zbuf, kind=8)) then
          allocate(zgrow(max(2 * size(zbuf, kind=8), zpos + field_cap)))
          zgrow(1:zpos) = zbuf(1:zpos)
          call move_alloc(zgrow, zbuf)
        end if
        nx1 = ixOmax1 - ixOmin1 + 1
        nx2 = ixOmax2 - ixOmin2 + 1
        nx3 = ixOmax3 - ixOmin3 + 1
        do iw = 1, nw
          if (write_sp) then
            zfp_nbytes = zfp_compress_field_adaptive_sp(real(ps(igrid)%w(&
               ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, iw), 4),&
               nx1, nx2, nx3, snap_zfp_rel_tol, snap_zfp_rel_eps,&
               zfp_prec(iw), zbuf, zpos + 1, size(zbuf, kind=8))
          else
            zfp_nbytes = zfp_compress_field_adaptive_dp(ps(igrid)%w(&
               ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, iw),&
               nx1, nx2, nx3, snap_zfp_rel_tol, snap_zfp_rel_eps,&
               zfp_prec(iw), zbuf, zpos + 1, size(zbuf, kind=8))
          end if
          my_field_nbytes(iw, iloc) = int(zfp_nbytes)
          zpos = zpos + zfp_nbytes
        end do
        my_rec(iloc) = sum(int(my_field_nbytes(:, iloc), 8))
      end do
    end if
#:endif

    ! 2. Offsets and metadata
    my_bytes = sum(my_rec(1:nloc))
    my_start = 0
    call MPI_EXSCAN(my_bytes, my_start, 1, MPI_INTEGER8, MPI_SUM, icomm,&
       ierrmpi)
    if (mype == 0) my_start = 0

    do ipe = 0, npe-1
      counts(ipe) = Morton_stop(ipe) - Morton_start(ipe) + 1
      displs(ipe) = Morton_start(ipe) - 1
    end do
    if (mype == 0) then
      allocate(rec_bytes(nleafs), block_nghost(2*ndim, nleafs))
      if (compress) allocate(zfp_block_nbytes(nw, nleafs))
    else
      allocate(rec_bytes(1), block_nghost(2*ndim, 1))
      if (compress) allocate(zfp_block_nbytes(nw, 1))
    end if
    call MPI_GATHERV(my_rec, nloc, MPI_INTEGER8, rec_bytes, counts, displs,&
       MPI_INTEGER8, 0, icomm, ierrmpi)
    counts_k = 2*ndim * counts
    displs_k = 2*ndim * displs
    call MPI_GATHERV(my_nghost, 2*ndim*nloc, MPI_INTEGER, block_nghost,&
       counts_k, displs_k, MPI_INTEGER, 0, icomm, ierrmpi)
    if (compress) then
      counts_k = nw * counts
      displs_k = nw * displs
      call MPI_GATHERV(my_field_nbytes, nw*nloc, MPI_INTEGER,&
         zfp_block_nbytes, counts_k, displs_k, MPI_INTEGER, 0, icomm, ierrmpi)
    end if

    ! Rank 0 writes the metadata, with the final offsets, and closes the file.
    ! block_nghost(2*ndim, nleafs) has the memory layout of the file's
    ! (ndim, 2, nleafs) table: minimum faces, then maximum faces.
    if (mype == 0) then
      call create_output_file(fh, snapnext, ".dat", "_snap")

      offset_tree_info = -1
      offset_block_data = -1
      call datfile_write_header1(fh, offset_tree_info, offset_block_data,&
         cons_wnames, nw, 6, .true.)
      call MPI_FILE_GET_POSITION(fh, offset_tree_info, ierrmpi)

      call write_forest(fh)

      allocate(block_ig(ndim, nleafs), block_lvl(nleafs))
      do Morton_no = 1, nleafs
        pnode => igrid_to_node(sfc(1, Morton_no), sfc(2, Morton_no))%node
        block_ig(:, Morton_no) = [ pnode%ig1,pnode%ig2,pnode%ig3 ]
        block_lvl(Morton_no) = pnode%level
      end do
      call MPI_FILE_WRITE(fh, block_lvl, nleafs, MPI_INTEGER, istatus,&
         ierrmpi)
      call MPI_FILE_WRITE(fh, block_ig, size(block_ig), MPI_INTEGER,&
         istatus, ierrmpi)
      call MPI_FILE_WRITE(fh, block_nghost, size(block_nghost), MPI_INTEGER,&
         istatus, ierrmpi)
      if (compress) call MPI_FILE_WRITE(fh, zfp_block_nbytes,&
         size(zfp_block_nbytes), MPI_INTEGER, istatus, ierrmpi)

      call MPI_FILE_GET_POSITION(fh, offset_offsets, ierrmpi)
      call MPI_TYPE_SIZE(MPI_OFFSET, offset_size, ierrmpi)
      offset_block_data = offset_offsets + int(nleafs, MPI_OFFSET_KIND)*&
         offset_size

      allocate(block_offset(nleafs))
      block_offset(1) = offset_block_data
      do Morton_no = 2, nleafs
        block_offset(Morton_no) = block_offset(Morton_no-1) +&
           rec_bytes(Morton_no-1)
      end do
      call MPI_FILE_WRITE(fh, block_offset, nleafs, MPI_OFFSET, istatus,&
         ierrmpi)

      call MPI_FILE_GET_POSITION(fh, pos, ierrmpi)
      if (pos /= offset_block_data) call mpistop(&
         "write_snapshot_parallel: metadata size differs from expected")

      call MPI_FILE_SEEK(fh, 0_MPI_OFFSET_KIND, MPI_SEEK_SET, ierrmpi)
      call datfile_write_header1(fh, offset_tree_info, offset_block_data,&
         cons_wnames, nw, 6, .true.)
      call MPI_FILE_CLOSE(fh, ierrmpi)
    end if

    ! After rank 0 closed the file, so it exists for everyone
    call MPI_BCAST(offset_block_data, 1, MPI_OFFSET, 0, icomm, ierrmpi)

    ! 3. Block data
    filename = output_filename(snapnext, "_snap")
    call MPI_FILE_OPEN(icomm, filename, MPI_MODE_WRONLY, MPI_INFO_NULL, fh,&
       ierrmpi)
    if (ierrmpi /= MPI_SUCCESS) call mpistop(&
       "write_snapshot_parallel: cannot open "//trim(filename))

    offset = offset_block_data + my_start
    offset_end = offset + my_bytes

    if (compress) then
      ! The streams are already contiguous: write them in slices
      cap = io_buffer_items(1)
      nrounds = sync_rounds(int((zpos + cap - 1) / cap))
      n_rec = 0
      do iround = 1, nrounds
        nfill = int(min(cap, zpos - n_rec))
        ! A slice of zbuf: through the MPI wrappers (USE_MPIWRAPPERS), which
        ! NVHPC 25.3 needs for buffers with an offset
        if (datfile_io_collective) then
          call mpi_file_write_at_all_wrapper(fh, offset,&
             zbuf(n_rec+1:n_rec+nfill), nfill, MPI_BYTE, istatus, ierrmpi)
        else
          call mpi_file_write_at_wrapper(fh, offset,&
             zbuf(n_rec+1:n_rec+nfill), nfill, MPI_BYTE, istatus, ierrmpi)
        end if
        if (ierrmpi /= MPI_SUCCESS) call mpistop(&
           "write_snapshot_parallel: writing the block data failed")
        offset = offset + nfill
        n_rec = n_rec + nfill
      end do
      if (n_rec /= zpos .or. offset /= offset_end) call mpistop(&
         "write_snapshot_parallel: not every stream was written")
    else
      ! Records in values of snap_size_real bytes
      allocate(my_val(max(1, nloc)))
      my_val = my_rec / size_real_out
      cap = io_buffer_items(size_real_out)
      nrounds = io_rounds(my_val(1:nloc), cap)
      allocate(buf(buffer_length(my_val(1:nloc), cap)))
      if (write_sp) allocate(buf_sp(size(buf)))

      m_next = Morton_start(mype)
      do iround = 1, nrounds
        ! Pack as many whole records as fit (always at least one)
        ib = 0
        do while (m_next <= Morton_stop(mype))
          n_rec = my_val(m_next - Morton_start(mype) + 1)
          if (ib > 0 .and. int(ib, 8) + n_rec > cap) exit
          ib0 = ib
          call pack_snap_values(sfc_to_igrid(m_next), buf, ib)
          if (ib - ib0 /= n_rec) call mpistop(&
             "write_snapshot_parallel: record size differs from expected")
          m_next = m_next + 1
        end do
        nfill = ib

        if (write_sp) then
          buf_sp(1:nfill) = real(buf(1:nfill), 4)
          if (datfile_io_collective) then
            call MPI_FILE_WRITE_AT_ALL(fh, offset, buf_sp, nfill, MPI_REAL4,&
               istatus, ierrmpi)
          else
            call MPI_FILE_WRITE_AT(fh, offset, buf_sp, nfill, MPI_REAL4,&
               istatus, ierrmpi)
          end if
        else
          if (datfile_io_collective) then
            call MPI_FILE_WRITE_AT_ALL(fh, offset, buf, nfill,&
               MPI_DOUBLE_PRECISION, istatus, ierrmpi)
          else
            call MPI_FILE_WRITE_AT(fh, offset, buf, nfill,&
               MPI_DOUBLE_PRECISION, istatus, ierrmpi)
          end if
        end if
        if (ierrmpi /= MPI_SUCCESS) call mpistop(&
           "write_snapshot_parallel: writing the block data failed")
        offset = offset + int(nfill, MPI_OFFSET_KIND) * size_real_out
      end do

      if (m_next /= Morton_stop(mype)+1 .or. offset /= offset_end) call&
         mpistop("write_snapshot_parallel: not every block was written")
    end if

    call MPI_FILE_CLOSE(fh, ierrmpi)
    call MPI_BARRIER(icomm, ierrmpi)

  end subroutine write_snapshot_parallel

  !> Read the block data of a v5 restart datfile with every rank reading its
  !> own blocks. Called by read_datfile after the header and the tree have been
  !> read and the local blocks allocated.
  !>
  !> fh0 is the restart file opened on MPI_COMM_SELF by rank 0 (unused on
  !> other ranks); it is closed here.
  subroutine read_datfile_blocks_parallel(fh0, offset_block_data, nw_found,&
     stagger_mark_dat)
    use mod_forest
    use mod_global_parameters
    use mod_input_output_helper, only: count_ix

    integer, intent(inout)                   :: fh0
    integer(kind=MPI_OFFSET_KIND), intent(in) :: offset_block_data
    integer, intent(in)                      :: nw_found
    logical, intent(in)                      :: stagger_mark_dat

    integer                       :: fh, igrid, Morton_no, ib, n_rec
    integer                       :: n_ghost(2*ndim), n_values
    integer                       :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    integer                       :: nrounds, iround, nfill, m_first,&
       m_next, offset_size
    integer                       :: istatus(MPI_STATUS_SIZE)
    integer(kind=8)               :: cap_dbl
    integer(kind=MPI_OFFSET_KIND) :: file_size
    integer(kind=MPI_OFFSET_KIND), allocatable :: block_offset(:)
    integer(kind=8), allocatable  :: rec_dbl(:)
    double precision, allocatable :: buf(:)

    call check_header_size()

    ! The offset table ends where the block data starts; the last block ends
    ! where the file does
    allocate(block_offset(nleafs+1))
    if (mype == 0) then
      call MPI_TYPE_SIZE(MPI_OFFSET, offset_size, ierrmpi)
      call MPI_FILE_READ_AT(fh0, offset_block_data - int(nleafs,&
         MPI_OFFSET_KIND)*offset_size, block_offset, nleafs, MPI_OFFSET,&
         istatus, ierrmpi)
      call MPI_FILE_GET_SIZE(fh0, file_size, ierrmpi)
      block_offset(nleafs+1) = file_size
      call MPI_FILE_CLOSE(fh0, ierrmpi)
    end if
    call MPI_BCAST(block_offset, nleafs+1, MPI_OFFSET, 0, icomm, ierrmpi)

    allocate(rec_dbl(nleafs))
    do Morton_no = 1, nleafs
      if (mod(block_offset(Morton_no+1) - block_offset(Morton_no),&
         int(size_double, MPI_OFFSET_KIND)) /= 0) call mpistop(&
         "read_datfile_blocks_parallel: block offsets are not whole doubles")
      rec_dbl(Morton_no) = (block_offset(Morton_no+1) -&
         block_offset(Morton_no)) / size_double
    end do

    call MPI_FILE_OPEN(icomm, restart_from_file, MPI_MODE_RDONLY,&
       MPI_INFO_NULL, fh, ierrmpi)
    if (ierrmpi /= MPI_SUCCESS) call mpistop(&
       "read_datfile_blocks_parallel: cannot open "//trim(restart_from_file))

    cap_dbl = io_buffer_doubles()
    nrounds = io_rounds(rec_dbl(Morton_start(mype):Morton_stop(mype)),&
       cap_dbl)
    allocate(buf(buffer_length(rec_dbl(Morton_start(mype):&
       Morton_stop(mype)), cap_dbl)))

    m_next = Morton_start(mype)
    do iround = 1, nrounds
      ! The same grouping of whole records as the writer
      m_first = m_next
      nfill = 0
      do while (m_next <= Morton_stop(mype))
        n_rec = int(rec_dbl(m_next))
        if (nfill > 0 .and. int(nfill, 8) + n_rec > cap_dbl) exit
        nfill = nfill + n_rec
        m_next = m_next + 1
      end do

      if (.not. datfile_io_collective) then
        ! Independent IO: a rank only has rounds with data
        call MPI_FILE_READ_AT(fh, block_offset(m_first), buf, nfill,&
           MPI_DOUBLE_PRECISION, istatus, ierrmpi)
      else if (m_first <= Morton_stop(mype)) then
        call MPI_FILE_READ_AT_ALL(fh, block_offset(m_first), buf, nfill,&
           MPI_DOUBLE_PRECISION, istatus, ierrmpi)
      else
        call MPI_FILE_READ_AT_ALL(fh, 0_MPI_OFFSET_KIND, buf, 0,&
           MPI_DOUBLE_PRECISION, istatus, ierrmpi)
      end if
      if (ierrmpi /= MPI_SUCCESS) call mpistop(&
         "read_datfile_blocks_parallel: reading the block data failed")

      ib = 0
      do Morton_no = m_first, m_next-1
        igrid = sfc_to_igrid(Morton_no)
        n_ghost = transfer(buf(ib+1:ib+nhdr_dbl), n_ghost)
        ixOmin1 = ixMlo1 - n_ghost(1)
        ixOmin2 = ixMlo2 - n_ghost(2)
        ixOmin3 = ixMlo3 - n_ghost(3)
        ixOmax1 = ixMhi1 + n_ghost(ndim+1)
        ixOmax2 = ixMhi2 + n_ghost(ndim+2)
        ixOmax3 = ixMhi3 + n_ghost(ndim+3)
        n_values = count_ix(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)&
           * nw_found
        if (stagger_mark_dat) n_values = n_values + count_ix(ixOmin1-1,&
           ixOmin2-1,ixOmin3-1,ixOmax1,ixOmax2,ixOmax3) * nws
        if (nhdr_dbl + n_values /= rec_dbl(Morton_no)) call mpistop(&
           "read_datfile_blocks_parallel: record size differs from its header")

        call store_datfile_block(igrid, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
           ixOmax2,ixOmax3, buf(ib+nhdr_dbl+1:ib+nhdr_dbl+n_values),&
           nw_found, stagger_mark_dat)
        ib = ib + nhdr_dbl + n_values
      end do
    end do

    if (m_next /= Morton_stop(mype)+1) call mpistop(&
       "read_datfile_blocks_parallel: not every block was read")

    call MPI_FILE_CLOSE(fh, ierrmpi)

  end subroutine read_datfile_blocks_parallel

  !> Put one block read from a v5 datfile into ps(igrid) and on the device.
  !>
  !> wbuf holds the cell-centred values on ixO (nw_found variables), then, if
  !> the file has them, the staggered ones. Variables the run has but the file
  !> lacks (or the other way round) are left to usr_transform_w, if set.
  subroutine store_datfile_block(igrid, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
     ixOmax2,ixOmax3, wbuf, nw_found, stagger_mark_dat)
    use mod_global_parameters
    use mod_usr_methods, only: usr_transform_w
    use mod_input_output_helper, only: count_ix

    integer, intent(in)          :: igrid
    integer, intent(in)          :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    double precision, intent(in) :: wbuf(:)
    integer, intent(in)          :: nw_found
    logical, intent(in)          :: stagger_mark_dat

    double precision, allocatable :: w(:,:,:,:)
    integer :: n_values_w

    allocate(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw_found))
    block=>ps(igrid)

    n_values_w = count_ix(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) *&
       nw_found
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw_found) =&
       reshape(wbuf(1:n_values_w), [ ixOmax1-ixOmin1+1, ixOmax2-ixOmin2+1,&
       ixOmax3-ixOmin3+1, nw_found ])
    if (stagger_mark_dat .and. stagger_grid) then
      ps(igrid)%ws(ixOmin1-1:ixOmax1,ixOmin2-1:ixOmax2,ixOmin3-1:ixOmax3,&
         1:nws) = reshape(wbuf(n_values_w+1:), [ ixOmax1-ixOmin1+2,&
         ixOmax2-ixOmin2+2, ixOmax3-ixOmin3+2, nws ])
    end if

    if (nw_found /= nw .and. associated(usr_transform_w)) then
      call usr_transform_w(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
         ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,nw_found,w,ps(igrid)%x,&
         ps(igrid)%w)
    else
      ps(igrid)%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:min(nw, nw_found)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:min(nw, nw_found))
    end if

#:if defined('FILL_NWEXTRA_ANALYTIC')
    ${GPU_UPDATE_DEVICE('bg(1)%w(:,:,:,1:nwgc,igrid)')}$
#:else
    ${GPU_UPDATE_DEVICE('bg(1)%w(:,:,:,:,igrid)')}$
#:endif

  end subroutine store_datfile_block

  !> Index range and value count (cell-centred plus staggered) of a block's
  !> record, as the default writer stores it. nghost_uniform is passed on to
  !> block_shape_io (snapshots with snap_nghost > 0).
  subroutine record_shape(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
     ixOmax2,ixOmax3, n_values, nghost_uniform)
    use mod_global_parameters
    use mod_input_output_helper, only: count_ix, block_shape_io

    integer, intent(in)  :: igrid
    integer, intent(out) :: n_ghost(2*ndim)
    integer, intent(out) :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    integer, intent(out) :: n_values
    integer, intent(in), optional :: nghost_uniform

    call block_shape_io(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3, n_values, nghost_uniform)
    if (stagger_grid) n_values = n_values + count_ix(ixOmin1-1,ixOmin2-1,&
       ixOmin3-1,ixOmax1,ixOmax2,ixOmax3) * nws

  end subroutine record_shape

  !> The record shape of a snapshot block, as write_datfile(is_snap=.true.)
  !> chooses it: a uniform halo of snap_nghost layers if snap_nghost > 0,
  !> otherwise the restart-file rule (physical-boundary ghost cells if
  !> save_physical_boundary)
  subroutine snap_record_shape(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,&
     ixOmax1,ixOmax2,ixOmax3, n_values)
    use mod_global_parameters

    integer, intent(in)  :: igrid
    integer, intent(out) :: n_ghost(2*ndim)
    integer, intent(out) :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    integer, intent(out) :: n_values

    if (snap_nghost > 0) then
      call record_shape(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
         ixOmax2,ixOmax3, n_values, snap_nghost)
    else
      call record_shape(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
         ixOmax2,ixOmax3, n_values)
    end if

  end subroutine snap_record_shape

  !> Append one snapshot block's values (w, then ws), with no header, to buf
  !> in double precision; the caller casts them to snap_size_real
  subroutine pack_snap_values(igrid, buf, ib)
    use mod_global_parameters

    integer, intent(in)             :: igrid
    double precision, intent(inout) :: buf(:)
    integer, intent(inout)          :: ib

    integer :: n_ghost(2*ndim), n_values, n
    integer :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3

    call snap_record_shape(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3, n_values)

    n = (ixOmax1-ixOmin1+1)*(ixOmax2-ixOmin2+1)*(ixOmax3-ixOmin3+1)*nw
    buf(ib+1:ib+n) = reshape(ps(igrid)%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1:nw), [ n ])
    ib = ib + n

    if (stagger_grid) then
      n = (ixOmax1-ixOmin1+2)*(ixOmax2-ixOmin2+2)*(ixOmax3-ixOmin3+2)*nws
      buf(ib+1:ib+n) = reshape(ps(igrid)%ws(ixOmin1-1:ixOmax1,&
         ixOmin2-1:ixOmax2,ixOmin3-1:ixOmax3,1:nws), [ n ])
      ib = ib + n
    end if

  end subroutine pack_snap_values

  !> Name of output file number ix: <base_filename><suffix>NNNN.dat, as
  !> create_output_file names it
  function output_filename(ix, suffix) result(filename)
    use mod_global_parameters

    integer, intent(in)          :: ix
    character(len=*), intent(in) :: suffix
    character(len=std_len)       :: filename

    write(filename, "(a,a,i4.4,a)") trim(base_filename), trim(suffix), ix,&
       ".dat"
  end function output_filename

  !> Append one block's record (ghost-cell header, w, then ws) to buf
  subroutine pack_record(igrid, buf, ib)
    use mod_global_parameters

    integer, intent(in)             :: igrid
    double precision, intent(inout) :: buf(:)
    integer, intent(inout)          :: ib

    integer :: n_ghost(2*ndim), n_values, n
    integer :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3

    call record_shape(igrid, n_ghost, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3, n_values)

    buf(ib+1:ib+nhdr_dbl) = transfer(n_ghost, buf(1:1), nhdr_dbl)
    ib = ib + nhdr_dbl

    n = (ixOmax1-ixOmin1+1)*(ixOmax2-ixOmin2+1)*(ixOmax3-ixOmin3+1)*nw
    buf(ib+1:ib+n) = reshape(ps(igrid)%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1:nw), [ n ])
    ib = ib + n

    if (stagger_grid) then
      n = (ixOmax1-ixOmin1+2)*(ixOmax2-ixOmin2+2)*(ixOmax3-ixOmin3+2)*nws
      buf(ib+1:ib+n) = reshape(ps(igrid)%ws(ixOmin1-1:ixOmax1,&
         ixOmin2-1:ixOmax2,ixOmin3-1:ixOmax3,1:nws), [ n ])
      ib = ib + n
    end if

  end subroutine pack_record

  !> The 2*ndim integer header of a record must be a whole number of doubles
  subroutine check_header_size()
    use mod_global_parameters

    if (2*ndim*size_int /= nhdr_dbl*size_double) call mpistop(&
       "mod_datfile_parallel: record header is not a whole number of doubles")
  end subroutine check_header_size

  !> Pack buffer capacity in doubles, from datfile_io_buffer_mb
  integer(kind=8) function io_buffer_doubles()
    use mod_global_parameters

    io_buffer_doubles = int(datfile_io_buffer_mb, 8) * 1024_8 * 1024_8 /&
       size_double
  end function io_buffer_doubles

  !> Rounds this rank runs: its own with independent IO; with collective IO
  !> every rank runs as many as the rank that needs the most
  integer function io_rounds(rec_dbl, cap)
    use mod_global_parameters

    integer(kind=8), intent(in) :: rec_dbl(:)
    integer(kind=8), intent(in) :: cap

    io_rounds = sync_rounds(count_rounds(rec_dbl, cap))
  end function io_rounds

  !> Rounds to run, given this rank's own nloc: nloc with independent IO, the
  !> largest nloc of any rank with collective IO
  integer function sync_rounds(nloc)
    use mod_global_parameters

    integer, intent(in) :: nloc

    if (datfile_io_collective) then
      call MPI_ALLREDUCE(nloc, sync_rounds, 1, MPI_INTEGER, MPI_MAX, icomm,&
         ierrmpi)
    else
      sync_rounds = nloc
    end if
  end function sync_rounds

  !> Pack buffer capacity in items of item_bytes bytes, from
  !> datfile_io_buffer_mb, capped so that one MPI call's count fits a default
  !> integer
  integer(kind=8) function io_buffer_items(item_bytes)
    use mod_global_parameters

    integer, intent(in) :: item_bytes

    io_buffer_items = min(int(datfile_io_buffer_mb, 8) * 1024_8 * 1024_8 /&
       item_bytes, int(huge(1), 8))
  end function io_buffer_items

  !> Number of rounds needed for these records: consecutive whole records
  !> are grouped up to cap (in the records' unit: doubles, reals or bytes),
  !> and a record bigger than cap gets a round of its own
  integer function count_rounds(rec_dbl, cap)
    integer(kind=8), intent(in) :: rec_dbl(:)
    integer(kind=8), intent(in) :: cap

    integer(kind=8) :: acc
    integer         :: i

    count_rounds = 0
    acc = 0
    do i = 1, size(rec_dbl)
      if (acc > 0 .and. acc + rec_dbl(i) > cap) then
        count_rounds = count_rounds + 1
        acc = 0
      end if
      acc = acc + rec_dbl(i)
    end do
    if (acc > 0) count_rounds = count_rounds + 1
  end function count_rounds

  !> Buffer length that holds any one round: the cap, but no more than all the
  !> records together and no less than the largest one (at least 1)
  integer function buffer_length(rec_dbl, cap)
    integer(kind=8), intent(in) :: rec_dbl(:)
    integer(kind=8), intent(in) :: cap

    integer(kind=8) :: total, biggest

    total = sum(rec_dbl)
    biggest = 0
    if (size(rec_dbl) > 0) biggest = maxval(rec_dbl)
    buffer_length = int(max(1_8, min(total, max(cap, biggest))))
  end function buffer_length

end module mod_datfile_parallel
