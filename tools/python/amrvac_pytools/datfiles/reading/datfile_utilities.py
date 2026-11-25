"""
Module containing reading and processing methods for an MPI-AMRVAC .datfiles file.

@author: Jannis Theunissen (original)
         Clément Robert (extensions, modifications)
         Adrian Kelly (v6)
         Olaf Willocx (v6)
Last edit: 1 August 2026
"""

import struct
import numpy as np

# Size of basic types (in bytes)
SIZE_LOGICAL = 4
SIZE_INT = 4
SIZE_DOUBLE = 8
NAME_LEN = 16

# For un-aligned data, use '=' (for aligned data set to '')
ALIGN = '='


def get_header(istream):
    """Read header from an MPI-AMRVAC 2.1 snapshot. This is compatible with versions down to 2.0.
    :param: istream     open datfile buffer with 'rb' mode
    :return: h          header information contained in a dictionary
    """
    istream.seek(0)
    h = {}

    fmt = ALIGN + 'i'
    [h['datfile_version']] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    if h['datfile_version'] < 3:
        raise IOError("Unsupported AMRVAC .datfiles file version: %d", h['datfile_version'])

    if h['datfile_version'] >= 6:
        fmt = ALIGN + 'i'
        # File origin.
        [h['origin']] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    h['compression'] = 'none'

    # Read scalar data at beginning of file
    fmt = ALIGN + 9 * 'i' + 'd'
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    [h['offset_tree'], h['offset_blocks'], h['nw'],
     h['ndir'], h['ndim'], h['levmax'], h['nleafs'], h['nparents'],
     h['it'], h['time']] = hdr

    # Read min/max coordinates
    fmt = ALIGN + h['ndim'] * 'd'
    h['xmin'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))
    h['xmax'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))

    # Read domain and block size (in number of cells)
    fmt = ALIGN + h['ndim'] * 'i'
    h['domain_nx'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))
    h['block_nx'] = np.array(
        struct.unpack(fmt, istream.read(struct.calcsize(fmt))))

    if h['datfile_version'] <= 5:
        h['size_real'] = SIZE_DOUBLE

    if h['datfile_version'] >= 5:
        # Read periodicity
        fmt = ALIGN + h['ndim'] * 'i' # Fortran logical is 4 byte int
        h['periodic'] = np.array(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))), dtype=bool)

        # Read geometry name
        fmt = ALIGN + NAME_LEN * 'c'
        hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        h['geometry'] = b''.join(hdr).strip().decode()

        # Read staggered flag
        fmt = ALIGN + 'i' # Fortran logical is 4 byte int
        h['staggered'] = bool(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt)))[0])

        if h['datfile_version'] >= 6:
            # Bytes per real in the block data (4 for single precision)
            fmt = ALIGN + 'i'
            [h['size_real']] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

            # Compression scheme name ('none', 'zfp', ...).
            # Scheme-specific metadata follows the name.
            fmt = ALIGN + NAME_LEN * 'c'
            hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
            h['compression'] = b''.join(hdr).strip().decode()
            if h['compression'] == 'zfp':
                # ZFP compression error criterion params.
                fmt = ALIGN + 'dd'
                h['zfp_rel_tol'], h['zfp_rel_eps'] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    # Read w_names
    w_names = []
    for i in range(h['nw']):
        fmt = ALIGN + NAME_LEN * 'c'
        hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        w_names.append(b''.join(hdr).strip().decode())
    h['w_names'] = w_names

    # Read physics type
    fmt = ALIGN + NAME_LEN * 'c'
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    h['physics_type'] = b''.join(hdr).strip().decode()

    # Read number of physics-defined parameters
    fmt = ALIGN + 'i'
    [n_pars] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    # First physics-parameter values are given, then their names
    fmt = ALIGN + n_pars * 'd'
    vals = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    fmt = ALIGN + n_pars * NAME_LEN * 'c'
    names = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    # Split and join the name strings (from one character array)
    names = [b''.join(names[i:i+NAME_LEN]).strip().decode()
             for i in range(0, len(names), NAME_LEN)]

    # Store the values corresponding to the names
    for val, name in zip(vals, names):
        h[name] = val
    return h


def get_tree_info(istream):
    """
    Read levels, morton-curve indices, and byte offsets for each block as stored in the datfile
    :param:  istream            Open datfile buffer with 'rb' mode.
    :return: block_lvls         NumPy array with block levels.
             block_ixs          NumPy array with block morton-curve indices
             block_offsets      NumPy array with block (data) offset in the datfile.
             block_nghost       NumPy array (nleafs, 2*ndim) of per-block ghost-cell
                                counts: [n_ghost_lo(ndim), n_ghost_hi(ndim)].
             block_field_bytes  NumPy array (nleafs, nw) with the byte count of
                                each block's per-field compressed stream, or None
                                for uncompressed files.
    """
    istream.seek(0)
    header = get_header(istream)
    nleafs = header['nleafs']
    nparents = header['nparents']
    ndim = header['ndim']
    block_field_bytes = None

    # Read tree info. Skip 'leaf' array
    istream.seek(header['offset_tree'] + (nleafs+nparents) * SIZE_LOGICAL)

    # Read block levels
    fmt = ALIGN + nleafs * 'i'
    block_lvls = np.array(struct.unpack(fmt, istream.read(struct.calcsize(fmt))))

    # Read block indices
    fmt = ALIGN + nleafs * header['ndim'] * 'i'
    block_ixs = np.reshape(struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
                           [nleafs, header['ndim']])

    gfmt = ALIGN + (2 * ndim) * 'i'
    bcsize = struct.calcsize(gfmt)

    if header['datfile_version'] >= 6:
        # From version 6, the per-block ghost-cell counts are stored in the tree
        # as Fortran n_ghost(ndim, 2, nleafs).
        fmt = ALIGN + nleafs * 2 * ndim * 'i'
        block_nghost = np.reshape(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            [nleafs, 2 * ndim])

        if header['compression'] != 'none':
            # Per-block per-field compressed byte counts, stored in the tree
            # as Fortran field_nbytes(nw, nleafs)
            fmt = ALIGN + nleafs * header['nw'] * 'i'
            block_field_bytes = np.reshape(
                struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
                [nleafs, header['nw']])

        # Offsets point directly at the block data.
        fmt = ALIGN + nleafs * 'q'
        block_offsets = np.array(struct.unpack(fmt, istream.read(struct.calcsize(fmt))))
    else:
        # Up to version 5, each stored offset points at that block's ghost-count,
        # `2*ndim` ints.
        fmt = ALIGN + nleafs * 'q'
        raw_offsets = np.array(struct.unpack(fmt, istream.read(struct.calcsize(fmt))))
        block_offsets = raw_offsets + bcsize          # start of block data

        # Per-block ghost-cell counts, read from each block's header.
        block_nghost = np.zeros((nleafs, 2 * ndim), dtype=int)
        for i, raw in enumerate(raw_offsets):
            istream.seek(int(raw))
            block_nghost[i] = struct.unpack(gfmt, istream.read(bcsize))
    return block_lvls, block_ixs, block_offsets, block_nghost, block_field_bytes


def block_shape_with_ghost(header, nghost):
    """Per-block array shape including its ghost halo: (block_nx + lo + hi, ..., nw)."""
    ndim = header['ndim']
    lo = np.asarray(nghost[:ndim])
    hi = np.asarray(nghost[ndim:])
    return np.append(np.asarray(header['block_nx']) + lo + hi, header['nw'])


def strip_ghost(block, header, nghost):
    """Return the interior (block_nx) part of a ghosted block, dropping the halo."""
    ndim = header['ndim']
    lo = np.asarray(nghost[:ndim])
    hi = np.asarray(nghost[ndim:])
    sl = tuple(slice(int(lo[k]), block.shape[k] - int(hi[k])) for k in range(ndim))
    return block[sl + (slice(None),)]


def get_single_block_data(istream, byte_offset, block_shape, size_real=SIZE_DOUBLE):
    """"
    Retrieve a specific block from the datfile
    :param: istream       open datfile buffer in 'rb' mode
    :param: byte_offset   offset of the given block in the datfile
    :param: block_shape   the shape of the block (list containing dimensions + number of variables)
    :param: size_real     bytes per real in the block data (header['size_real'], 4 or 8)
    :return: block_data   numpy array containing the block data, with dimensions equal to block_shape
    """
    istream.seek(byte_offset)
    # Read actual data
    fmt = ALIGN + int(np.prod(block_shape)) * ('f' if size_real == 4 else 'd')
    d = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    # Fortran ordering
    block_data = np.reshape(d, block_shape, order='F')
    return block_data


def get_single_block_data_zfp(istream, byte_offset, field_bytes, block_shape):
    """
    Retrieve a specific ZFP-compressed block from the datfile
    :param: istream       open datfile buffer in 'rb' mode
    :param: byte_offset   offset of the given block in the datfile
    :param: field_bytes   per-field compressed byte counts (block_field_bytes row)
    :param: block_shape   the shape of the block (list containing dimensions + number of variables)
    :return: block_data   numpy array with dimensions equal to block_shape
    """
    import zfpy
    ndim = len(block_shape) - 1
    # There is padding due to e.g. blocks 4^3 in ndim = 3, so slice with this.
    unpad = tuple(slice(0, int(n)) for n in block_shape[:ndim])
    fields = []
    pos = int(byte_offset)
    for nb in field_bytes:
        istream.seek(pos)
        # Float kind is already deduced from the per-field ZFP header.
        arr = zfpy.decompress_numpy(istream.read(int(nb)))
        # Fortran ordering
        fields += [arr.T[unpad]]
        pos += int(nb)
    return np.stack(fields, axis=-1)


def read_block(dataset, ileaf, keep_ghost=False):
    """
    Read leaf block `ileaf` at its true (ghosted) shape, stripping the halo
    unless keep_ghost=True. Handles both ghosted and ghost-free files.
    """
    nghost = dataset.block_nghost[ileaf]
    shape = block_shape_with_ghost(dataset.header, nghost)
    match dataset.header['compression']:
      case 'zfp':
          block = get_single_block_data_zfp(
              dataset.file, dataset.block_offsets[ileaf], dataset.block_field_bytes[ileaf], shape)
      case 'none' | _:
          size_real = dataset.header.get('size_real', SIZE_DOUBLE)
          block = get_single_block_data(dataset.file, dataset.block_offsets[ileaf], shape, size_real)
    if not keep_ghost:
        block = strip_ghost(block, dataset.header, nghost)
    return block


def get_blocks(dataset, keep_ghost=False):
    """
    Reads all block data from an MPI-AMRVAC 2.0 snapshot.
    :param dataset      instance of 'amrvac_reader.load_file' class
    :param keep_ghost   if True, return each block with its ghost halo (for
                        per-block gradients / periodic-edge stencils); if False
                        (default) strip the halo so blocks tile the domain.
    :return list containing block data as dictionaries with level, morton index,
            data ('w') and ghost counts ('nghost').
    """
    size_real = dataset.header.get('size_real', SIZE_DOUBLE)
    hdr = dataset.header
    blocks = []
    for ileaf, offset in enumerate(dataset.block_offsets):
        nghost = dataset.block_nghost[ileaf]
        shape = block_shape_with_ghost(hdr, nghost)
        if hdr['compression'] == 'zfp':
            block = get_single_block_data_zfp(dataset.file, offset,
                                              dataset.block_field_bytes[ileaf], shape)
        else:
            block = get_single_block_data(dataset.file, offset, shape, size_real)
        if not keep_ghost:
            block = strip_ghost(block, hdr, nghost)
        lvl = dataset.block_lvls[ileaf]
        ix = dataset.block_ixs[ileaf]
        b = {'lvl': lvl, 'ix': ix, 'w': block, 'nghost': nghost}
        blocks.append(b)

    return blocks



def get_uniform_data(dataset):
    """
    Retrieves the data for a uniform data set.
    :param dataset: instance of 'amrvac_reader.load_file' class
    :return The raw data as a NumPy array.
    """
    blocks = get_blocks(dataset)
    hdr = dataset.header

    refined_nx = 2 ** (hdr['levmax'] - 1) * hdr['domain_nx']
    domain_shape = np.append(refined_nx, hdr['nw'])
    d = np.zeros(domain_shape, order='F')

    for b in blocks:
        i0 = (b['ix'] - 1) * hdr['block_nx']
        i1 = i0 + hdr['block_nx']
        if hdr['ndim'] == 1:
            d[i0[0]:i1[0], :] = b['w']
        elif hdr['ndim'] == 2:
            d[i0[0]:i1[0], i0[1]:i1[1], :] = b['w']
        elif hdr['ndim'] == 3:
            d[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2], :] = b['w']
    return d


def get_upsampled_data(istream):
    """
    Retrieves the data of a uniform or AMR data set on the finest level's grid.
    A block on a coarser level is upsampled piecewise-constant, each of its cells
    repeated 2**(levmax-level) times per direction, so the result shows the cells
    the solver actually had rather than an interpolation of them. Unlike
    get_uniform_data, which assumes every block is on levmax, this places
    multi-level data correctly; on single-level data the two agree exactly.
    Works for 1D, 2D and 3D files, including the (ndim-1)-dimensional slices
    written with slice_type = 'dat'.
    :param: istream     open datfile buffer with 'rb' mode
    :return The raw data as a NumPy array of shape
            (domain_nx * 2**(levmax-1), nw), domain_nx having ndim entries.
    """
    hdr = get_header(istream)
    block_lvls, block_ixs, block_offsets = get_tree_info(istream)
    block_nx = hdr['block_nx']
    block_shape = np.append(block_nx, hdr['nw'])

    refined_nx = 2 ** (hdr['levmax'] - 1) * hdr['domain_nx']
    d = np.full(np.append(refined_nx, hdr['nw']), np.nan)

    for lvl, ix, offset in zip(block_lvls, block_ixs, block_offsets):
        block = get_single_block_data(istream, offset, block_shape)
        factor = 2 ** (hdr['levmax'] - lvl)
        for axis in range(hdr['ndim']):
            block = np.repeat(block, factor, axis=axis)
        # ix is the 1-based block index on the block's own level
        i0 = (ix - 1) * block_nx * factor
        d[tuple(slice(a, a + n) for a, n in zip(i0, block.shape[:-1]))] = block

    if np.isnan(d).any():
        raise ValueError("the leaf blocks do not cover the domain")
    return d
