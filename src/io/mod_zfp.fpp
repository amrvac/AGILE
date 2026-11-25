module mod_zfp

#:if defined('COMPRESS_ZFP')
  use, intrinsic :: iso_c_binding
  use zfp
  use mod_comm_lib, only: mpistop
  implicit none
  private

  public :: zfp_padded_dim
  public :: zfp_max_field_nbytes
  public :: zfp_compress_field_adaptive_dp
  public :: zfp_compress_field_adaptive_sp

  !> Staging arrays, dynamically reallocated.
  double precision, allocatable, target, save :: pad_dp(:,:,:), dec_dp(:,:,:)
  real(kind=4),     allocatable, target, save :: pad_sp(:,:,:), dec_sp(:,:,:)

  !> Scratch streams of the precision search: the current trial, and the
  !> smallest-precision stream that satisfied the tolerance so far (so that it
  !> does not have to be recomputed).
  integer(kind=1), allocatable, target, save :: trial_buf(:), best_buf(:)
  !> Their size: the worst-case stream of a field of this shape and type.
  integer(kind=8), save :: stream_nbytes = 0
  integer, save         :: stream_shape(3) = 0
  logical, save         :: stream_dp = .false.

  !> Corrected zFORp interfaces. They lack the value attribute, not usable.
  interface
    function zfp_write_header_c(stream, field, mask) result(num_bits)&
        bind(c, name="zfp_write_header")
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value    :: stream, field
      integer(c_int), value :: mask
      integer(c_size_t)     :: num_bits
    end function zfp_write_header_c
    function zfp_read_header_c(stream, field, mask) result(num_bits)&
        bind(c, name="zfp_read_header")
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value    :: stream, field
      integer(c_int), value :: mask
      integer(c_size_t)     :: num_bits
    end function zfp_read_header_c
  end interface

contains

  !> Need to pad, ZFP uses 4^3 blocks.
  elemental function zfp_padded_dim(n) result(np)
    integer, intent(in) :: n
    integer             :: np
    np = ((n+3)/4)*4
  end function zfp_padded_dim

  !> Worst case (compression ratio 1) number of bytes of a compressed field,
  !> with header. (nx, ny, nz) are unpadded.
  function zfp_max_field_nbytes(nx, ny, nz, write_dp, prec) result(nbytes)
    integer, intent(in) :: nx, ny, nz, prec
    logical, intent(in) :: write_dp
    integer(kind=8)     :: nbytes

    integer(kind=1), target :: scratch(64)
    type(zFORp_bitstream)   :: bitstream
    type(zFORp_stream)      :: stream
    type(zFORp_field)       :: field
    integer                 :: iprec

    field = zFORp_field_3d(c_null_ptr,&
        merge(zFORp_type_double, zFORp_type_float, write_dp),&
        zfp_padded_dim(nx), zfp_padded_dim(ny), zfp_padded_dim(nz))
    bitstream = zFORp_bitstream_stream_open(c_loc(scratch), int(size(scratch), 8))
    stream = zFORp_stream_open(bitstream)
    iprec = zFORp_stream_set_precision(stream, prec)

    ! zFORp_stream_maximum_size describes the field data. Add the header bits,
    ! but round to the word size. (x+63)/64 is ceiling.
    nbytes = zFORp_stream_maximum_size(stream, field)+&
        ((zfORp_header_max_bits+63)/64)*8

    call zFORp_stream_close(stream)
    call zFORp_bitstream_stream_close(bitstream)
    call zFORp_field_free(field)
  end function zfp_max_field_nbytes

  !> Size trial_buf and best_buf for the worst-case stream of an nx*ny*nz
  !> field (compression ratio 1 and precision 64) when the shape changes.
  subroutine reserve_streams(nx, ny, nz, write_dp)
    integer, intent(in) :: nx, ny, nz
    logical, intent(in) :: write_dp

    if (all(stream_shape == [nx,ny,nz]) .and. (stream_dp .eqv. write_dp))&
        return
    stream_shape = [nx,ny,nz]
    stream_dp = write_dp
    stream_nbytes = zfp_max_field_nbytes(nx, ny, nz, write_dp, 64)
    if (allocated(trial_buf)) deallocate(trial_buf, best_buf)
    allocate(trial_buf(stream_nbytes), best_buf(stream_nbytes))
  end subroutine reserve_streams

  !> Compress an already padded field. The header goes first.
  function compress_padded(field_ptr, scalar_type, px, py, pz, prec, buf,&
      buf_start, buf_size) result(nbytes)
    type(c_ptr), intent(in)     :: field_ptr
    integer, intent(in)         :: scalar_type, px, py, pz, prec
    integer(kind=8), intent(in) :: buf_start, buf_size
    integer(kind=1), target, intent(inout) :: buf(buf_size)
    integer(kind=8)             :: nbytes

    type(zFORp_bitstream) :: bitstream
    type(zFORp_stream)    :: stream
    type(zFORp_field)     :: field
    integer               :: iprec
    integer(kind=8)       :: header_bits

    field = zFORp_field_3d(field_ptr, scalar_type, px, py, pz)
    bitstream = zFORp_bitstream_stream_open(&
      c_loc(buf(buf_start)), buf_size-buf_start+1)
    stream = zFORp_stream_open(bitstream)
    iprec = zFORp_stream_set_precision(stream, prec)
    header_bits = zfp_write_header_c(&
        transfer(stream, c_null_ptr),&
        transfer(field, c_null_ptr),&
        int(zFORp_header_full, c_int))
    if (header_bits == 0) call mpistop("mod_zfp: writing zfp header failed")
    nbytes = zFORp_compress(stream, field)
    if (nbytes == 0) call mpistop("mod_zfp: zfp compression failed")

    call zFORp_stream_close(stream)
    call zFORp_bitstream_stream_close(bitstream)
    call zFORp_field_free(field)
  end function compress_padded

#:for K, FTYPE, ZTYPE, IS_DP in [&
  ('dp', 'double precision', 'zFORp_type_double', '.true.'), &
  ('sp', 'real(kind=4)', 'zFORp_type_float', '.false.')]

  !> Compress a field q as Q into buf(buf_start:) at the smallest* ZFP precision
  !> p in [1, 64] satisfying
  !>   |Q - q| <= rel_tol*(|q| + rel_eps*max|q|)
  !> in every cell of the block, max|q| over the block, with rel_tol and
  !> rel_eps user-side parameters.
  !> Returns the stream's byte count.
  !> The given prec is the first guess for p on entry if in [1, 64], otherwise
  !> we make a guess based on rel_tol. On return prec is p.
  !> Each step compresses, decompresses and measures the error, taken as
  !> nonincreasing in p. From the first guess the search steps away by 1, 2,
  !> 4, etc., until the outcome flips, then bisects. It takes 2 steps when the
  !> guess is right, a few more when it is close, at most ~13. Converges faster
  !> than simple bisection, typically takes 3 steps, although bisection has a
  !> better worst case of log2(64) = 6 steps.
  !> When rel_tol is below about epsilon/rel_eps it may not be
  !> attainable in every block in a worst case scenario. Then p = 64 is used.
  !> *: The smallest p is found when compressed size increases monotonically
  !>    with p, which is almost always true. If false, a solution to the
  !>    inequality is still found, but may not be optimal.
  function zfp_compress_field_adaptive_${K}$(w, nx, ny, nz, rel_tol, rel_eps,&
      prec, buf, buf_start, buf_size) result(nbytes)
    integer, intent(in)          :: nx, ny, nz
    ${FTYPE}$, intent(in)        :: w(nx,ny,nz)
    double precision, intent(in) :: rel_tol, rel_eps
    integer, intent(inout)       :: prec
    integer(kind=8), intent(in)  :: buf_start, buf_size
    integer(kind=1), target, intent(inout) :: buf(buf_size)
    integer(kind=8)              :: nbytes

    integer(kind=1), allocatable :: swap(:)
    integer          :: i, j, k, px, py, pz, p, p0, lo, hi, step
    integer(kind=8)  :: trial_nbytes, best_nbytes
    double precision :: l_max, rel_tol_measured
    logical          :: pass

    px = zfp_padded_dim(nx)
    py = zfp_padded_dim(ny)
    pz = zfp_padded_dim(nz)

    if (allocated(pad_${K}$)) then
      if (any(shape(pad_${K}$) /= [px,py,pz]))&
        deallocate(pad_${K}$, dec_${K}$)
    end if
    if (.not. allocated(pad_${K}$))&
      allocate(pad_${K}$(px,py,pz), dec_${K}$(px, py, pz))
    call reserve_streams(nx,ny,nz, ${IS_DP}$)

    ! Pad by edge replication on the right. It's important to keep the
    ! exponent bits uniform over the data, it improves compressability.
    do k = 1, pz
    do j = 1, py
    do i = 1, px
      pad_${K}$(i,j,k) = w(min(i, nx), min(j, ny), min(k, nz))
    end do
    end do
    end do

    l_max = maxval(abs(dble(w)))

    ! lo: the largest  precision known to fail
    ! hi: the smallest precision known to pass
    lo = 0
    hi = 65
    p0 = prec
    if (p0 < 1 .or. p0 > 64) then
      ! We estimate p0. We expect rel_tol >= c*2^-p for p bitplanes and c some
      ! constant. Pattern appears to hold quite well for reasonable values of
      ! rel_tol, and we find c around 6 for the doubleGEM3D MHD setup, so using
      ! that for now. This only really matters if no estimate is given by prec.
      p0 = min(64, max(1, ceiling(log(6.0d0/rel_tol)/log(2.0d0))))
    end if
    p = p0
    step = 1
    best_nbytes = 0

    do
      trial_nbytes = compress_padded(c_loc(pad_${K}$), ${ZTYPE}$,&
          px, py, pz, p, trial_buf, 1_8, stream_nbytes)
      rel_tol_measured = rel_tol_measured_${K}$(trial_nbytes, w, nx, ny, nz,&
          l_max, rel_eps)
      pass = rel_tol_measured <= rel_tol

      if (pass .or. p == 64) then
        ! Keep the stream of the last one that passes, or 64 if none pass, so
        ! that we don't have to recompress at the end. Instead of copying, just
        ! swap the pointers.
        best_nbytes = trial_nbytes
        call move_alloc(trial_buf, swap)
        call move_alloc(best_buf, trial_buf)
        call move_alloc(swap, best_buf)
      end if
      ! 3 modes while converging:
      !   1: hi = 65, moving -> geometrically
      !   2: lo = 0,  moving <- geometrically
      !   3: [lo,hi] in [1,64], bisecting
      if (pass) then
        ! Passed, so at most this p.
        hi = p
      else
        ! Failed, so at least this p.
        lo = p
      end if
      ! converged
      if (hi == lo+1) exit
      if (hi == 65) then
        ! mode 1
        p = min(64, p0+step)
        step = 2*step
      else if (lo == 0) then
        ! mode 2
        p = max(1, p0-step)
        step = 2*step
      else
        ! mode 3
        p = (lo+hi)/2
      end if
    end do

    if (buf_start+best_nbytes-1 > buf_size)&
        call mpistop("mod_zfp: compressed stream exceeds the field buffer")
    buf(buf_start:buf_start+best_nbytes-1) = best_buf(1:best_nbytes)
    nbytes = best_nbytes
    prec = min(hi, 64)
  end function zfp_compress_field_adaptive_${K}$

  !> Decompress the trial_buf stream and measure the smallest effective
  !> rel_tol it satisfies, i.e.
  !>   |Q-q|/(|q|+rel_eps*max|q|),
  !> the max taken over the block. This is then to be compared to the desired
  !> rel_tol.
  function rel_tol_measured_${K}$(stream_nbytes, w, nx, ny, nz, l_max, rel_eps)&
      result(rel_tol_measured)
    integer(kind=8), intent(in)  :: stream_nbytes
    integer, intent(in)          :: nx, ny, nz
    ${FTYPE}$, intent(in)        :: w(nx,ny,nz)
    double precision, intent(in) :: l_max, rel_eps
    double precision             :: rel_tol_measured

    type(zFORp_bitstream) :: bitstream
    type(zFORp_stream)    :: stream
    type(zFORp_field)     :: field
    integer(kind=8)       :: header_bits, nb
    integer               :: i, j, k
    double precision      :: q, q_measured

    bitstream = zFORp_bitstream_stream_open(c_loc(trial_buf), stream_nbytes)
    stream = zFORp_stream_open(bitstream)
    field = zFORp_field_alloc()
    header_bits = zfp_read_header_c(&
        transfer(stream, c_null_ptr),&
        transfer(field, c_null_ptr),&
        int(zFORp_header_full, c_int))
    if (header_bits == 0) call mpistop("mod_zfp: reading zfp header failed")
    call zFORp_field_set_pointer(field, c_loc(dec_${K}$))
    nb = zFORp_decompress(stream, field)
    if (nb == 0) call mpistop("mod_zfp: zfp decompression failed")
    call zFORp_stream_close(stream)
    call zFORp_bitstream_stream_close(bitstream)
    call zFORp_field_free(field)

    rel_tol_measured = 0.0d0
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
      q = dble(w(i,j,k))
      q_measured = dble(dec_${K}$(i,j,k))
      ! Careful with l_max = 0.
      rel_tol_measured = max(rel_tol_measured,&
          abs(q_measured-q)/max(abs(q)+rel_eps*l_max, tiny(1.0d0)))
    end do
    end do
    end do
  end function rel_tol_measured_${K}$

#:endfor

#:endif

end module mod_zfp
