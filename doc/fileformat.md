title: AGILE data file format

> this doc describes the latest version of the data file format, which is numbered "6"

AGILE writes two kinds of `.dat` file, whose layouts are described below:

* **datfiles** (`<base_filename>NNNN.dat`, file type 2, `dtsave_dat`) are
  always written in double precision, version 5. They are the full-fidelity
  output and the files to restart from.
* **snapshots** (`<base_filename>_snapNNNN.dat`, file type 6, `dtsave_snap`)
  are written in version 6, by default in double precision (`snap_size_real = 4`
  for single precision, roughly half the file size), and support lossy compression.
  They are intended for analysis; a run cannot restart from them.

Each file holds a single time level, and can be converted to other data formats
directly usable for visualization. Note that restart is possible on a differing
number of CPUs, and may suddenly allow more refinement levels. Also, note that
individual files will typically have different lengths, as the number of grid
blocks will vary dynamically. You can find the exact implementation in
`src/io/mod_input_output.fpp`, more specifically in `datfile_write_header()`
and `datfile_read_header()`.

A `.dat` file contains a header, mesh tree information, and block
data, in the following order:

# Header (v5)

```fortran
integer           :: Version number
integer           :: Byte offset where tree information starts
integer           :: Byte offset where block data starts
integer           :: nw
integer           :: ndir
integer           :: ndim
integer           :: levmax
integer           :: nleafs
integer           :: nparents
integer           :: it
double precision  :: global_time
double precision  :: xprobmin(ndim)
double precision  :: xprobmax(ndim)
integer           :: domain_nx(ndim)
integer           :: block_nx(ndim)
logical           :: periodic(ndim)
character(len=16) :: geometry
logical           :: staggered
character(len=16) :: w_names(nw)
character(len=16) :: physics_type
! The physics parameters, such as the adiabatic gamma
integer           :: n_params
double precision  :: parameters(n_params)
character(len=16) :: parameter_names(n_params)
! Indexes for file output (for restarting)
integer           :: datfilenext ! the next datfile index
integer           :: slicenext
integer           :: collapsenext
```

# Header (v6)


```fortran
integer           :: Version number
integer           :: Data origin (1 for AGILE)
integer           :: Byte offset where tree information starts
integer           :: Byte offset where block data starts
integer           :: nw
integer           :: ndir
integer           :: ndim
integer           :: levmax
integer           :: nleafs
integer           :: nparents
integer           :: it
double precision  :: global_time
double precision  :: xprobmin(ndim)
double precision  :: xprobmax(ndim)
integer           :: domain_nx(ndim)
integer           :: block_nx(ndim)
logical           :: periodic(ndim)
character(len=16) :: geometry
logical           :: staggered
integer           :: size_real
character(len=16) :: compression
double precision  :: zfp_rel_tol ! present if compression == 'zfp'
double precision  :: zfp_rel_eps ! present if compression == 'zfp'
character(len=16) :: w_names(nw)
character(len=16) :: physics_type
! The physics parameters, such as the adiabatic gamma
integer           :: n_params
double precision  :: parameters(n_params)
character(len=16) :: parameter_names(n_params)
! Indexes for file output (for restarting)
integer           :: datfilenext ! next index in the datfile or snapshot series
integer           :: slicenext
integer           :: collapsenext
integer           :: snapnext ! in a datfile, the next snapshot
```

# Tree information (v5)

```fortran
logical :: leaf(nleafs+nparents)
integer :: refinement_level(nleafs)
integer :: spatial_index(ndim, nleafs)
integer(kind=MPI_OFFSET_KIND) :: offset_block(nleafs)
```

# Tree information (v6)

```fortran
logical :: leaf(nleafs+nparents)
integer :: refinement_level(nleafs)
integer :: spatial_index(ndim, nleafs)
integer :: n_ghost(ndim, 2, nleafs)
integer :: field_nbytes(nw, nleafs) ! present if compression == 'zfp'
integer(kind=MPI_OFFSET_KIND) :: offset_block(nleafs)
```

# Block 1 to nleafs (v5)

```fortran
integer :: n_ghost_lo(ndim)
integer :: n_ghost_hi(ndim)
! block_shape = 1-n_ghost_lo:block_nx+n_ghost_hi
double precision :: w(block_shape, nw)
```

# Block 1 to nleafs (v6, uncompressed)

```fortran
! block_shape = 1-n_ghost(:,1):block_nx+n_ghost(:,2), n_ghost from the tree
double precision :: w(block_shape, nw) ! if size_real = 8
real             :: w(block_shape, nw) ! if size_real = 4
```

Ghost cells are nonzero for physical-boundary blocks when `save_physical_boundary=T`,
and (v6 snapshots only) on every block face when `snap_nghost>0`.

# Block 1 to nleafs (v6, ZFP compressed)

A block record is the concatenation of `nw` independent canonical
fixed-precision ZFP streams (with headers), one per field, their lengths given
by `field_nbytes` in the tree. Before compression each field is padded in each
direction to a multiple of 4 (ZFP's 4^3 chunks); after decompression that
padding is cropped, and the fields give the same layout as an uncompressed v6
block:
```fortran
double precision :: w(block_shape, nw) ! if size_real = 8
real             :: w(block_shape, nw) ! if size_real = 4
```

# Version history

## Version 1

Version 1 contained the following information

    1. block data (nw variables)
    2. leaf/parent logical array
    3. header:
        nx^D
        domain_nx^D
        xprobmin^D
        xprobmax^D
        nleafs
        levmax
        ndim
        ndir
        nw
        it
        global_time

The idea is that you can reconstruct the full grid when you know the Morton
order used for the leaf/parent logical array.

## Version 2

Version 2 had the same information as version 1, but changes were made to the
Morton order on the coarse grid, causing incompatibility.

## Version 5 (current for datfiles)

Version 5 introduced the `geometry` parameter (e.g. "polar_2D",
"cartesian_1.75D"... ), as well as a `periodic` and a `staggered` flags.
The periodic flag is an `ndim`-long 1D boolean array. For each direction,
it defaults to `.false.` but is set to `.true.` if at least one
quantity is using a periodic boundary condition in the corresponding
direction. This change was motivated by improving the (upcoming)
compatibility with `yt`.  The `staggered` flag is added to support
staggered grid for constrained transport MHD.

## Version 6 (current for analysis snapshots)

Version 6 is written by the snapshot output stream only currently
(file type 6, `<base_filename>_snapNNNN.dat`) and is intended for analysis;
use datfiles (always v5, and always double precision) for restarts.

Version 6's header adds to version 5's (see both above): `origin`,
`size_real`, the compression scheme and its parameters, and `snapnext`.
`origin` describes the origin of the file, 1 for AGILE.
`size_real` is 4 or 8, describing single or double precision block data.

Version 6 stores the per-block ghost-cell counts in the tree section as
`integer :: n_ghost(ndim, 2, nleafs)`, as opposed to in the blocks.

Version 6 snapshots support LLNL ZFP lossy compression in fixed-precision mode
for 3D data. Each field in each block is separately compressed as a standard
ZFP stream (containing a header and data) satisfying some modified relative
error criterion, see [Analysis snapshots](par.html#par_snapshots).
The block constituting contiguous fields is replaced by a contiguous set of
ZFP streams decompressing into said fields, plus some padding due to ZFP's 4^3 chunk size.
