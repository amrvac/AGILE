module mod_load_balance
#ifdef USE_MPIWRAPPERS
  use mod_mpi_wrapper
#else
#define mpi_irecv_wrapper MPI_IRECV
#define mpi_isend_wrapper MPI_ISEND
#endif
  implicit none
  private
  !> MPI recv send variables for AMR.  One request per neighbouring rank, so
  !> these are sized by npe.
  integer :: itag, irecv, isend
  integer :: recv_igrid, recv_ipe, send_igrid, send_ipe, igrid
  integer, dimension(:), allocatable :: recvrequest, sendrequest
  integer, dimension(:,:), allocatable :: recvstatus, sendstatus
  !> MPI recv send variables for staggered-variable AMR
  integer, dimension(:), allocatable :: recvrequest_stg, sendrequest_stg

  !> MPI buffers to send blocks.
  !>
  !> The migration is aggregated per destination rank: every block bound for
  !> one peer occupies a contiguous run of snd_buff_lb and travels in a single
  !> Isend.  Both ranks order a run by the key
  !>
  !>     ikey = recv_igrid
  !>
  !> the block index the block will occupy on its new owner.  That is unique
  !> among the blocks arriving at a given receiver, because each comes from its
  !> own getnode(), and the sender knows it because load_balance walks the
  !> globally replicated Morton list in the same order on every rank and every
  !> rank calls getnode() for every migration - so igrid_inuse, and with it the
  !> destination index, stays replicated.
  !>
  !> Sizing matters more here than the message count.  These buffers used to be
  !> block_nx^3 * nw * max_blocks doubles each, allocated on the first call
  !> whether or not a single block ever moved: 1.05 GB apiece at block_nx=16,
  !> nw=8, max_blocks=4000, and device-resident.  They are now sized by what
  !> actually migrates, and only ever grow, so the steady state neither wastes
  !> the memory nor reallocates every regrid.
  !>
  !> Deliberately NOT declare create: these are resized, and an OpenACC
  !> declare-create allocatable does not reliably get a fresh device copy when
  !> the host array is reallocated - an update device afterwards then writes
  !> through a stale mapping.  They are managed explicitly with enter data
  !> create / exit data delete instead, the same way mod_fix_conserve handles
  !> its own resizable exchange buffers.
  double precision, allocatable, dimension(:) :: snd_buff_lb, rcv_buff_lb
  !> high-water mark of the two buffers, in chunks
  integer :: max_buff = 0

  !> Per migrating block, appended during the walk and turned into per-peer
  !> runs by exchange_migrated_blocks.  lb_snd_igrid / lb_snd_ibuf and the rcv
  !> counterparts are read by the pack and unpack kernels and are therefore
  !> device-resident; the rest is the host-side working set.
  integer :: n_lb_snd, n_lb_rcv
  integer, allocatable, dimension(:) :: lb_snd_igrid, lb_snd_ibuf
  integer, allocatable, dimension(:) :: lb_snd_dest, lb_snd_key
  integer, allocatable, dimension(:) :: lb_rcv_igrid, lb_rcv_ibuf
  integer, allocatable, dimension(:) :: lb_rcv_src, lb_rcv_key
  !$acc declare create(lb_snd_igrid,lb_snd_ibuf,lb_rcv_igrid,lb_rcv_ibuf)

  !> The peers themselves, ascending, with the extent of each one's run.  Host
  !> only: these drive the MPI calls and nothing else.
  integer :: n_lb_send_pe, n_lb_recv_pe
  integer, allocatable, dimension(:) :: lb_send_pe, lb_send_pe_off,&
     lb_send_pe_len
  integer, allocatable, dimension(:) :: lb_recv_pe, lb_recv_pe_off,&
     lb_recv_pe_len

  public :: load_balance

contains
  !> reallocate blocks into processors for load balance
  subroutine load_balance
    use mod_forest
    use mod_global_parameters
    use mod_space_filling_curve
    use mod_amr_solution_node, only: getnode,putnode
    use mod_functions_forest, only: change_ipe_tree_leaf

    integer :: Morton_no, ipe
    !> MPI recv send variables for staggered-variable AMR
    integer :: itag_stg
    integer, dimension(:,:), allocatable :: recvstatus_stg, sendstatus_stg

    ! Jannis: for now, not using version for passive/active blocks
    call get_Morton_range()

    if (npe==1) then
       sfc_to_igrid(:)=sfc(1,Morton_start(mype):Morton_stop(mype))
       return
    end if

    irecv=0
    isend=0
    n_lb_snd=0
    n_lb_rcv=0
    if (.not.allocated(recvrequest)) then
       allocate(recvstatus(MPI_STATUS_SIZE,npe),recvrequest(npe),&
           sendstatus(MPI_STATUS_SIZE,npe),sendrequest(npe))
       allocate(lb_send_pe(npe),lb_send_pe_off(npe),lb_send_pe_len(npe),&
           lb_recv_pe(npe),lb_recv_pe_off(npe),lb_recv_pe_len(npe))
    end if
    recvrequest=MPI_REQUEST_NULL
    sendrequest=MPI_REQUEST_NULL

    if(stagger_grid) then
      allocate(recvstatus_stg(MPI_STATUS_SIZE,max_blocks*3),&
         recvrequest_stg(max_blocks*3), sendstatus_stg(MPI_STATUS_SIZE,&
         max_blocks*3),sendrequest_stg(max_blocks*3))
      recvrequest_stg=MPI_REQUEST_NULL
      sendrequest_stg=MPI_REQUEST_NULL
    end if

    ! Only the descriptors are allocated up front, and they are integers.  The
    ! data buffers are sized by what actually migrates, in
    ! exchange_migrated_blocks, once the walk below has counted it.
    if ( .not. allocated(lb_snd_igrid) ) then
       allocate( lb_snd_igrid(max_blocks), lb_snd_ibuf(max_blocks), &
            lb_snd_dest(max_blocks), lb_snd_key(max_blocks), &
            lb_rcv_igrid(max_blocks), lb_rcv_ibuf(max_blocks), &
            lb_rcv_src(max_blocks), lb_rcv_key(max_blocks) )
    end if

    do ipe=0,npe-1; do Morton_no=Morton_start(ipe),Morton_stop(ipe)
       recv_ipe=ipe

       send_igrid=sfc(1,Morton_no)
       send_ipe=sfc(2,Morton_no)

       if (recv_ipe/=send_ipe) then
          ! get an igrid number for the new node in recv_ipe processor
          recv_igrid=getnode(recv_ipe)
          ! update node igrid and ipe on the tree
          call change_ipe_tree_leaf(recv_igrid,recv_ipe,send_igrid,send_ipe)
          ! receive physical data of the new node in recv_ipe processor
          if (recv_ipe==mype) call lb_recv
          ! send physical data of the old node in send_ipe processor
          if (send_ipe==mype) call lb_send
       end if
       if (recv_ipe==mype) then
          if (recv_ipe==send_ipe) then
             sfc_to_igrid(Morton_no)=send_igrid
          else
             sfc_to_igrid(Morton_no)=recv_igrid
          end if
       end if
    end do; end do

    ! Lay the blocks out per peer, exchange them in one message each, and apply
    ! them.  The walk above only recorded descriptors.
    call exchange_migrated_blocks

    ! The staggered path still posts per block; see the note in lb_send.  It is
    ! unreachable here (stagger_grid is false and fix_edges is not called), so
    ! it is left as it was found.
    if(stagger_grid) then
      if (irecv>0) call MPI_WAITALL(irecv,recvrequest_stg,recvstatus_stg,&
         ierrmpi)
      if (isend>0) call MPI_WAITALL(isend,sendrequest_stg,sendstatus_stg,&
         ierrmpi)
   end if

    if(stagger_grid) deallocate(recvstatus_stg,recvrequest_stg,sendstatus_stg,&
       sendrequest_stg)

    ! post processing
    do ipe=0,npe-1; do Morton_no=Morton_start(ipe),Morton_stop(ipe)
       recv_ipe=ipe

       send_igrid=sfc(1,Morton_no)
       send_ipe=sfc(2,Morton_no)

       if (recv_ipe/=send_ipe) then
          !if (send_ipe==mype) call dealloc_node(send_igrid)
          call putnode(send_igrid,send_ipe)
       end if
    end do; end do


    ! Update sfc array: igrid and ipe info in space filling curve
    call amr_Morton_order()

  end subroutine load_balance

  !> Exchange and apply the migrating blocks, one message per neighbouring
  !> rank, with the buffers sized by what actually moves.
  !>
  !> The walk in load_balance has already allocated every arriving block with
  !> alloc_node and recorded a descriptor for each departure and arrival. Here
  !> those are laid out as one contiguous run per peer, ordered by the shared
  !> key, packed by a single kernel and sent as a single message each.
  subroutine exchange_migrated_blocks
    use mod_global_parameters
    use mod_msg_layout, only: layout_runs
    use mod_comm_lib, only: mpistop

    integer :: k, nchunk, nbuf, nneed
    integer :: igrid, ibuf, iw, ix1, ix2, ix3
    integer, allocatable :: chunksize(:)

    ! Only 1:nwgc travels. The analytic extras past it - for ffhd the frozen
    ! field - were already rebuilt by fill_nwextra_device inside the
    ! alloc_node that lb_recv calls for every arriving block, so putting them
    ! on the wire would only overwrite identical values.
    nchunk = block_nx1*block_nx2*block_nx3*nwgc

    ! Every chunk is one block interior, so they are all the same size.
    allocate(chunksize(max(n_lb_snd,n_lb_rcv,1)))
    chunksize = nchunk
    call layout_runs(n_lb_snd, lb_snd_dest, lb_snd_key, chunksize, lb_snd_ibuf,&
       n_lb_send_pe, lb_send_pe, lb_send_pe_off, lb_send_pe_len)
    call layout_runs(n_lb_rcv, lb_rcv_src, lb_rcv_key, chunksize, lb_rcv_ibuf,&
       n_lb_recv_pe, lb_recv_pe, lb_recv_pe_off, lb_recv_pe_len)
    deallocate(chunksize)

    ! Grow the buffers to fit, never shrink: load_balance runs on every regrid,
    ! so resizing down would thrash the device heap for nothing.  The device
    ! copy is torn down and recreated around the reallocation rather than
    ! updated afterwards - see the note on the declarations.
    nneed = max(n_lb_snd, n_lb_rcv)
    if (nneed > max_buff) then
      if (allocated(snd_buff_lb)) then
        !$acc exit data delete(snd_buff_lb, rcv_buff_lb)
        deallocate(snd_buff_lb, rcv_buff_lb)
      end if
      max_buff = nneed
      allocate(snd_buff_lb(max_buff*nchunk), rcv_buff_lb(max_buff*nchunk))
      !$acc enter data create(snd_buff_lb, rcv_buff_lb)
    end if
    if (n_lb_snd == 0 .and. n_lb_rcv == 0) return

    if (n_lb_snd > 0) then
      !$acc update device(lb_snd_igrid(1:n_lb_snd), lb_snd_ibuf(1:n_lb_snd))
    end if
    if (n_lb_rcv > 0) then
      !$acc update device(lb_rcv_igrid(1:n_lb_rcv), lb_rcv_ibuf(1:n_lb_rcv))
    end if

    ! One Irecv per peer, straight into that peer's run. A peer sends exactly
    ! one message, so (communicator, source) already disambiguates and the tag
    ! carries nothing.
    itag = 0
    if (n_lb_recv_pe > 0) then
#ifndef NOGPUDIRECT
      !$acc host_data use_device(rcv_buff_lb)
#endif
      do k = 1, n_lb_recv_pe
        call mpi_irecv_wrapper(rcv_buff_lb(lb_recv_pe_off(k)),&
           lb_recv_pe_len(k),MPI_DOUBLE_PRECISION,lb_recv_pe(k),itag,icomm,&
           recvrequest(k),ierrmpi)
      end do
#ifndef NOGPUDIRECT
      !$acc end host_data
#endif
    end if

    ! Pack every departing block in one kernel. bg(1)%w carries the grid index
    ! last, so the block can be selected by a device-side index.
    if (n_lb_snd > 0) then
      !$acc parallel loop gang default(present) private(igrid,ibuf)
      do k = 1, n_lb_snd
         igrid = lb_snd_igrid(k)
         ibuf  = lb_snd_ibuf(k)
         !$acc loop collapse(4) vector
         do iw = 1, nwgc
            do ix3 = 1, block_nx3
               do ix2 = 1, block_nx2
                  do ix1 = 1, block_nx1
                     snd_buff_lb(ibuf + (ix1-1) + block_nx1*(ix2-1) &
                          + block_nx1*block_nx2*(ix3-1) &
                          + block_nx1*block_nx2*block_nx3*(iw-1)) = &
                          bg(1)%w(ixMlo1-1+ix1, ixMlo2-1+ix2, ixMlo3-1+ix3,&
                                  iw, igrid)
                  end do
               end do
            end do
         end do
      end do
    end if

    if (n_lb_send_pe > 0) then
#ifdef NOGPUDIRECT
      !$acc update host(snd_buff_lb(1:n_lb_snd*nchunk))
#else
      !$acc host_data use_device(snd_buff_lb)
#endif
      do k = 1, n_lb_send_pe
        call mpi_isend_wrapper(snd_buff_lb(lb_send_pe_off(k)),&
           lb_send_pe_len(k),MPI_DOUBLE_PRECISION,lb_send_pe(k),itag,icomm,&
           sendrequest(k),ierrmpi)
      end do
#ifndef NOGPUDIRECT
      !$acc end host_data
#endif
    end if

    if (n_lb_recv_pe > 0) then
      call MPI_WAITALL(n_lb_recv_pe,recvrequest,recvstatus,ierrmpi)
      ! With one aggregated message per peer a disagreement about which blocks
      ! the run holds is no longer an MPI mismatch, so check the arrival
      ! lengths. A difference in order cannot arise: both ranks sort the same
      ! run by the same key.
      do k = 1, n_lb_recv_pe
        call MPI_GET_COUNT(recvstatus(:,k),MPI_DOUBLE_PRECISION,nbuf,ierrmpi)
        if (nbuf /= lb_recv_pe_len(k)) call mpistop( &
           "exchange_migrated_blocks: message length disagrees with the layout")
      end do
#ifdef NOGPUDIRECT
      !$acc update device(rcv_buff_lb(1:n_lb_rcv*nchunk))
#endif
    end if

    ! Apply every arriving block in one kernel.
    if (n_lb_rcv > 0) then
      !$acc parallel loop gang default(present) private(igrid,ibuf)
      do k = 1, n_lb_rcv
         igrid = lb_rcv_igrid(k)
         ibuf  = lb_rcv_ibuf(k)
         !$acc loop collapse(4) vector
         do iw = 1, nwgc
            do ix3 = 1, block_nx3
               do ix2 = 1, block_nx2
                  do ix1 = 1, block_nx1
                     bg(1)%w(ixMlo1-1+ix1, ixMlo2-1+ix2, ixMlo3-1+ix3,&
                             iw, igrid) &
                          = rcv_buff_lb(ibuf + (ix1-1) + block_nx1*(ix2-1) &
                            + block_nx1*block_nx2*(ix3-1) &
                            + block_nx1*block_nx2*block_nx3*(iw-1))
                  end do
               end do
            end do
         end do
      end do
    end if

    if (n_lb_send_pe > 0) then
      call MPI_WAITALL(n_lb_send_pe,sendrequest,sendstatus,ierrmpi)
    end if

  end subroutine exchange_migrated_blocks

  subroutine lb_recv
    use mod_global_parameters
    use mod_amr_solution_node, only: alloc_node
    use mod_comm_lib, only: mpistop

    ! alloc_node also re-derives the analytic extras (for ffhd, the frozen
    ! field) through fill_nwextra_device, which is what lets the transfer below
    ! carry only 1:nwgc rather than 1:nw.
    call alloc_node(recv_igrid)

    ! Record the block; exchange_migrated_blocks posts them all once the walk
    ! is over, one message per peer.  The key is the index this block will
    ! occupy here, which the sender knows too.
    irecv=irecv+1
    n_lb_rcv=n_lb_rcv+1
    if (n_lb_rcv > size(lb_rcv_igrid)) then
       call mpistop('load_balance: more arrivals than max_blocks')
    end if
    lb_rcv_igrid(n_lb_rcv) = recv_igrid
    lb_rcv_src(n_lb_rcv)   = send_ipe
    lb_rcv_key(n_lb_rcv)   = recv_igrid
    if(stagger_grid) then
       itag=recv_igrid+max_blocks
       call mpi_irecv_wrapper(ps(recv_igrid)%ws,1,type_block_io_stg,send_ipe,itag,&
             icomm,recvrequest_stg(irecv),ierrmpi)
    end if
  end subroutine lb_recv

  subroutine lb_send
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    ! Record the block; the pack and the post both happen in
    ! exchange_migrated_blocks.  The key is the receiving block index, so it
    ! matches the one the receiver derives from its own getnode().
    isend=isend+1
    n_lb_snd=n_lb_snd+1
    if (n_lb_snd > size(lb_snd_igrid)) then
       call mpistop('load_balance: more departures than max_blocks')
    end if
    lb_snd_igrid(n_lb_snd) = send_igrid
    lb_snd_dest(n_lb_snd)  = recv_ipe
    lb_snd_key(n_lb_snd)   = recv_igrid

    ! The staggered branch below is NOT aggregated and is left as found: it is
    ! unreachable (stagger_grid is false throughout this fork and fix_edges has
    ! no caller).  Its tag, recv_igrid+max_blocks, is sound but doubles the tag
    ! range, so it would be the first thing to hit MPI_TAG_UB.
    if(stagger_grid) then
       itag=recv_igrid+max_blocks
       call mpi_isend_wrapper(ps(send_igrid)%ws,1,type_block_io_stg,recv_ipe,itag,&
             icomm,sendrequest_stg(isend),ierrmpi)
    end if
  end subroutine lb_send

end module mod_load_balance
