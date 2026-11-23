title: AGILE data file format

> this doc describes the latest version of the data file format, which is numbered "6"

AGILE writes two kinds of `.dat` file, which share the layout described below:

* **datfiles** (`<base_filename>NNNN.dat`, file type 2, `dtsave_dat`) are
  always written in double precision, version 5. They are the full-fidelity
  output and the files to restart from.
* **snapshots** (`<base_filename>_snapNNNN.dat`, file type 6, `dtsave_snap`)
  are written in version 6, by default in single precision (roughly half the
  file size, and see `snap_size_real`).
  They are intended for analysis; you cannot restart from them.

Each file holds a single time level, and can be converted to other data formats
directly usable for visualization. Note that restart is possible on a differing
number of CPUs, and may suddenly allow more refinement levels. Also, note that
individual files will typically have different lengths, as the number of grid
blocks will vary dynamically. You can find the exact implementation in
`src/io/mod_input_output.fpp`, more specifically in `datfile_write_header()`
and `datfile_read_header()`.

A `.dat` file contains a header, mesh tree information, and block
data, in the following order:

# Header

Comments here mean present under said condition.

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
integer           :: size_real ! version >= 6
character(len=16) :: compression ! version >= 6
integer           :: zfp_precision ! compression == 'zfp'
character(len=16) :: w_names(nw)
character(len=16) :: physics_type
! The physics parameters, such as gamma
integer           :: n_params
double precision  :: parameters(n_params)
character(len=16) :: parameter_names(n_params)
! Indexes for file output (for restarting)
integer           :: datfilenext
integer           :: slicenext
integer           :: collapsenext
```

# Tree information

```fortran
logical :: leaf(nleafs+nparents)
integer :: refinement_level(nleafs)
integer :: spatial_index(ndim, nleafs)
integer :: n_ghost(ndim, 2, nleafs) ! version >= 6
integer :: field_nbytes(nw, nleafs) ! version >= 6, present when compression /= 'none'
integer(kind=MPI_OFFSET_KIND) :: offset_block(nleafs)
```

# Block 1 to nleafs (uncompressed)

```fortran
integer :: n_ghost_lo(ndim) ! version <= 5, number of ghost cells on lower boundaries
integer :: n_ghost_hi(ndim) ! version <= 5, number of ghost cells on upper boundaries
! block_shape = 1-n_ghost_lo:block_nx+n_ghost_hi
! double precision up to version 5, single precision from version 6
double precision :: w(block_shape, nw)
```

Ghost cells are nonzero for physical-boundary blocks when `save_physical_boundary=T`,
and (v6 snapshots only) on every block face when `snap_nghost>0`. This may help in
post processing, avoids having to reconstruct ghost cells.

# Block 1 to nleafs (ZFP compressed, version >=6)

A block record is the concatenation of `nw` independent canonical
fixed-precision ZFP streams (with headers). The fields are padded up to
multiples of 4^3, as required by ZFP. After decompression and concatenation,
the layout is the same as for any version >=6 file.
```
double precision :: w(block_shape, nw)
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

## Version 5 (current)

Version 5 introduced the `geometry` parameter (e.g. "polar_2D",
"cartesian_1.75D"... ), as well as a `periodic` and a `staggered` flags.
The periodic flag a `ndir`-long 1D boolean array. For each direction,
it defaults to `.false.` but is set to `.true.` if at least one
quantity is using a periodic boundary condition in the corresponding
direction. This change was motivated by improving the (upcoming)
compatibility with `yt`.  The `staggered` flag is added to support
staggered grid for constrained transport MHD.

## Version 6 (current)

Version 6 is written by the snapshot output stream (file type 6,
`<base_filename>_snapNNNN.dat`) and is intended for analysis; use datfiles
(always v5, and always double precision) for restarts.

Version 6 has the same header as version 5, but stores the block data in
single precision by default, roughly halving the file size.
(`snap_size_real=8` in `filelist` keeps double precision.)
A `size_real` integer (bytes per real: 4 or 8) is added to
the header, so this defines the actual block real kind for version 6 files.

Version 6 stores the per-block ghost-cell counts in the tree section as
`integer :: n_ghost(ndim, 2, nleafs)`, as opposed to in the blocks.

Version 6 supports storing ghost cells per block more generally, so that they
don't have to be reconstructed and reshaped during further analysis.

Version 6 snapshots support LLNL ZFP compression in lossy fixed-precision mode.
Currently unsupported in tandem with `stagger_grid`.
See https://doi.org/10.1137/18M1168832 in regard to the parameter
`zfp_precision`, to get a sense of how this translates into relative error
statistics. It is recommended to not use compression in combination with
`snap_size_real=4` since this is likely to hurt (especially block-average)
accuracy over `snap_size_real=8`, but it will work.
