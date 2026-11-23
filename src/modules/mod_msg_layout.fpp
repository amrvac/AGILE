!> Laying an exchange out as one message per destination rank.
!>
!> An exchange that moves many small chunks between the same pair of ranks does
!> not need one message per chunk. If both ranks can put a peer's chunks in the
!> same order without talking to each other, every chunk bound for that peer can
!> occupy one contiguous run of a single buffer and travel in a single Isend.
!> The message count is then the number of neighbouring ranks rather than the
!> number of chunks, and the tag no longer has to encode a block index, so it
!> stays trivially inside MPI_TAG_UB.
!>
!> The ordering is supplied by a key that both ranks compute for the same chunk
!> and that totally orders a peer's run - typically built from the *receiving*
!> block's index, which the sender reads out of its neighbour tables and the
!> receiver reads off its own igrid. layout_runs sorts by it; rank_keys does the
!> sorting.
!>
!> Two things to keep in mind when using this:
!>
!>  - both ranks must enumerate from tree state that is already consistent,
!>    before any message is posted, and the chunk sizes must agree chunk for
!>    chunk;
!>  - aggregation gives up MPI's own mismatch detection, because a disagreement
!>    about which chunks a run holds is no longer a length mismatch on any
!>    single message. Check each arrival's MPI_GET_COUNT against the run length
!>    that was laid out. That catches a difference in the *set* of chunks the
!>    two ranks enumerated; a difference in order cannot arise, since both sort
!>    the same run by the same key.
!>
!> Used by mod_fix_conserve (refluxing) and mod_coarsen_refine (coarsening).
module mod_msg_layout

  implicit none
  private

  public :: layout_runs
  public :: rank_keys

contains

  !> Give every peer one contiguous run of the exchange buffer, with the chunks
  !> inside a run ordered by the key. Both ranks compute the same key for the
  !> same chunk, so sorting by it makes the sender's run and the receiver's run
  !> agree element for element without any handshake.
  !>
  !> size_of(k) is the extent of chunk k in buffer elements; peers come out
  !> ascending, peer_off/peer_len delimit their runs, and off_of(k) is where
  !> chunk k starts.
  subroutine layout_runs(n, pe_of, key_of, size_of, off_of, npeer, peer,&
     peer_off, peer_len)
    use mod_global_parameters, only: npe

    integer, intent(in)  :: n, pe_of(:), key_of(:), size_of(:)
    integer, intent(out) :: off_of(:), npeer, peer(:), peer_off(:), peer_len(:)

    integer :: k, p, j
    integer :: nper(0:npe-1), cursor(0:npe-1)
    integer, allocatable :: perm(:)

    npeer = 0
    if (n == 0) return

    ! how long each peer's run is, in buffer elements
    nper = 0
    do k = 1, n
      nper(pe_of(k)) = nper(pe_of(k)) + size_of(k)
    end do

    ! peers ascending, each run following the previous one
    j = 1
    do p = 0, npe-1
      if (nper(p) == 0) cycle
      npeer           = npeer + 1
      peer(npeer)     = p
      peer_off(npeer) = j
      peer_len(npeer) = nper(p)
      cursor(p)       = j
      j               = j + nper(p)
    end do

    ! walk the chunks in ascending key order and hand each one the next slot
    ! in its peer's run, so the two ranks fill a run identically
    allocate(perm(n))
    call rank_keys(n, key_of, perm)
    do k = 1, n
      j         = perm(k)
      p         = pe_of(j)
      off_of(j) = cursor(p)
      cursor(p) = cursor(p) + size_of(j)
    end do
    deallocate(perm)

  end subroutine layout_runs

  !> Rank key(1:n) ascending: perm(k) is the index of the k-th smallest key.
  !> Heapsort, so O(n log n) worst case and no scratch beyond perm itself.
  !> Local rather than the equivalent mrgrnk vendored with octree-mg, to keep
  !> the callers independent of the multigrid solver.
  subroutine rank_keys(n, key, perm)
    integer, intent(in)  :: n, key(:)
    integer, intent(out) :: perm(:)

    integer :: i, j, l, ir, tmp

    do i = 1, n
      perm(i) = i
    end do
    if (n < 2) return

    l  = n/2 + 1
    ir = n
    do
      if (l > 1) then
        l   = l - 1
        tmp = perm(l)
      else
        tmp      = perm(ir)
        perm(ir) = perm(1)
        ir       = ir - 1
        if (ir == 1) then
          perm(1) = tmp
          exit
        end if
      end if
      i = l
      j = l + l
      do while (j <= ir)
        if (j < ir) then
          if (key(perm(j)) < key(perm(j+1))) j = j + 1
        end if
        if (key(tmp) < key(perm(j))) then
          perm(i) = perm(j)
          i = j
          j = j + j
        else
          j = ir + 1
        end if
      end do
      perm(i) = tmp
    end do

  end subroutine rank_keys

end module mod_msg_layout
