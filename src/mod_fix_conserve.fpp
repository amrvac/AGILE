!> Module for flux conservation near refinement boundaries
module mod_fix_conserve
#ifdef USE_MPIWRAPPERS
  use mod_mpi_wrapper
#else
#define mpi_irecv_wrapper MPI_IRECV
#define mpi_isend_wrapper MPI_ISEND
#endif
  implicit none
  private

  type fluxalloc
     double precision, dimension(:,:,:,:,:), allocatable :: flux
     !!double precision, dimension(:,:,:,:), pointer:: flux => null()
     !!double precision, dimension(:,:,:,:), pointer:: edge => null()
  end type fluxalloc
  !> store flux to fix conservation
  type(fluxalloc), dimension(:,:), allocatable, public :: pflux

  integer, save                        :: nrecv, nsend
  !> Extent of recvbuffer / sendbuffer actually used by the current exchange.
  !> The buffers are grown and never shrunk, so size() can exceed this, and the
  !> host/device transfers have to be sliced by it rather than moving the whole
  !> array.
  integer, save                        :: fc_recvsize, fc_sendsize
  double precision, allocatable, save  :: recvbuffer(:), sendbuffer(:)
  integer, dimension(:), allocatable   :: fc_recvreq, fc_sendreq
  integer, dimension(:,:), allocatable :: fc_recvstat, fc_sendstat
  integer, dimension(3), save        :: isize
  
  !The flux exchange is aggregated per destination rank: every chunk bound for
  !one peer occupies a single contiguous run of sendbuffer and travels in one
  !Isend, so the message count is the number of neighbouring ranks rather than
  !the number of (block, face, child) triples. That is what keeps the tag space
  !trivially inside MPI_TAG_UB - the tag no longer encodes a block index at all
  !- and it replaces the earlier scheme of striping tags across several
  !duplicated communicators to widen the range.
  !
  !It works because both ranks can order a peer's run identically without
  !talking to each other. The ordering key is the one ibuf_offset is already
  !indexed by,
  !
  !    ikey = 4**3*(ineighbor-1) + inc1 + 4*inc2 + 16*inc3
  !
  !which the sender builds from the receiving block's index on *its* rank and
  !the receiver from its own igrid, both arriving at the same number. Sorting
  !each run by it makes the two layouts agree element for element. The key also
  !fixes the chunk's size: inc_d lies in {1,2} exactly when i_d is zero, so the
  !one direction with inc_d in {0,3} is idims, and the size is isize(idims).
  !
  !ibuf_offset therefore keeps its meaning - key -> offset in recvbuffer, read
  !by fix_conserve - and only its values change: a chunk now sits at its peer's
  !run start plus the sizes of all lower-keyed chunks from that peer.
  integer, allocatable, save         :: ibuf_offset(:)
  !
  !The peers themselves, ascending, with the extent of each one's run. Host
  !only: these drive the MPI calls and nothing else.
  integer, save                      :: n_send_pe, n_recv_pe
  integer, allocatable, save         :: send_pe(:), send_pe_off(:), send_pe_len(:)
  integer, allocatable, save         :: recv_pe(:), recv_pe_off(:), recv_pe_len(:)
  !
  !Per outgoing chunk, all built by build_message_layout. snd_igrid and
  !snd_ibuf are read by the pack kernels and are therefore device-resident;
  !snd_key, snd_dest and snd_dims are the host-side working set the layout is
  !derived from. snd_group lists the chunks that read one (iside, idims) slot
  !of pflux, which is what lets a single kernel serve them all: the slot is
  !then loop-invariant, and a derived-type component holding an allocatable
  !array cannot be selected by a device-side index.
  integer, allocatable, save         :: snd_igrid(:), snd_ibuf(:)
  integer, allocatable, save         :: snd_key(:), snd_dest(:), snd_dims(:)
  integer, allocatable, save         :: snd_group(:,:,:)
  integer, save                      :: snd_ngroup(2,3)
  !
  !The same working set for the incoming chunks, from which ibuf_offset and the
  !recv_pe_* runs are built.
  integer, allocatable, save         :: rcv_key(:), rcv_src(:), rcv_dims(:)
  integer, allocatable, save         :: rcv_off(:)
  integer, save                      :: nwflux_fc
  integer, dimension(3,3), save      :: nxCo_fc
  !$acc declare create(isize, nxCo_fc, nwflux_fc)

  integer                              :: ibuf, ibuf_send
  ! ct for corner total
  integer, save                        :: nrecv_ct, nsend_ct
  ! buffer for corner coarse
  double precision, allocatable, save  :: recvbuffer_cc(:), sendbuffer_cc(:)
  integer, dimension(:), allocatable   :: cc_recvreq, cc_sendreq
  integer, dimension(:,:), allocatable :: cc_recvstat, cc_sendstat
  integer, dimension(3), save        :: isize_stg
  integer                              :: ibuf_cc, ibuf_cc_send
  integer                              :: itag, itag_cc, isend, isend_cc,&
      irecv, irecv_cc

  public :: init_comm_fix_conserve
  public :: sendflux
  public :: recvflux
  public :: store_flux
  public :: store_edge
  public :: fix_conserve
  public :: fix_edges

 contains

   subroutine init_comm_fix_conserve(idimmin,idimmax,nwfluxin)
     use mod_global_parameters
     use mod_comm_lib, only: mpistop

     integer, intent(in) :: idimmin,idimmax,nwfluxin

     integer :: iigrid, igrid, idims, iside, i1,i2,i3, nxCo1,nxCo2,nxCo3
     integer :: ic1,ic2,ic3, inc1,inc2,inc3, ipe_neighbor
     integer :: recvsize, sendsize
     integer :: recvsize_cc, sendsize_cc
     ! MPI tag out of bounds safeguard

     ! JESSENEW
     nwflux_fc = nwfluxin

     if (.not.allocated(pflux(1,1)%flux)) call mpistop(&
        "init_comm_fix_conserve: pflux%flux is not allocated yet")

     nxCo_fc(1,1)=1
     nxCo_fc(2,1)=size(pflux(1,1)%flux,2)/2
     nxCo_fc(3,1)=size(pflux(1,1)%flux,3)/2

     nxCo_fc(1,2)=size(pflux(1,2)%flux,1)/2
     nxCo_fc(2,2)=1
     nxCo_fc(3,2)=size(pflux(1,2)%flux,3)/2

     nxCo_fc(1,3)=size(pflux(1,3)%flux,1)/2
     nxCo_fc(2,3)=size(pflux(1,3)%flux,2)/2
     nxCo_fc(3,3)=1

     if (nxCo_fc(2,1)/=(ixMhi2-ixMlo2+1)/2 .or. nxCo_fc(3,1)/=(ixMhi3-ixMlo3+&
        1)/2 .or. nxCo_fc(1,2)/=(ixMhi1-ixMlo1+1)/2 .or. nxCo_fc(3,&
        2)/=(ixMhi3-ixMlo3+1)/2 .or. nxCo_fc(1,3)/=(ixMhi1-ixMlo1+1)/2 .or. &
        nxCo_fc(2,3)/=(ixMhi2-ixMlo2+1)/2) call mpistop(&
        "init_comm_fix_conserve: pflux%flux shape disagrees with the mesh")

     if (nwfluxin > size(pflux(1,1)%flux,4)) call mpistop(&
        "init_comm_fix_conserve: nwfluxin exceeds the pflux%flux w extent")

     nsend    = 0
     nrecv    = 0
     recvsize = 0
     sendsize = 0
     if(stagger_grid) then
       ! Special communication for diagonal 'coarse corners'
       nsend_ct=0
       nrecv_ct=0
       recvsize_cc=0
       sendsize_cc=0
     end if

     do idims= idimmin,idimmax
       select case (idims)
         case (1)
         nrecv=nrecv+nrecv_fc(1)
         nsend=nsend+nsend_fc(1)
         isize(1)=nxCo_fc(1,1)*nxCo_fc(2,1)*nxCo_fc(3,1)*(nwfluxin)
         recvsize=recvsize+nrecv_fc(1)*isize(1)
         sendsize=sendsize+nsend_fc(1)*isize(1)
         if(stagger_grid) then
           ! This does not consider the 'coarse corner' case
           nxCo1=1;nxCo2=ixGhi2/2-nghostcells+1;nxCo3=ixGhi3/2-nghostcells+1;
           isize_stg(1)=nxCo1*nxCo2*nxCo3*(3-1)
           ! the whole size is used (cell centered and staggered)
           isize(1)=isize(1)+isize_stg(1)
           recvsize=recvsize+nrecv_fc(1)*isize_stg(1)
           sendsize=sendsize+nsend_fc(1)*isize_stg(1)
           ! Coarse corner case
           nrecv_ct=nrecv_ct+nrecv_cc(1)
           nsend_ct=nsend_ct+nsend_cc(1)
           recvsize_cc=recvsize_cc+nrecv_cc(1)*isize_stg(1)
           sendsize_cc=sendsize_cc+nsend_cc(1)*isize_stg(1)
         end if

         case (2)
         nrecv=nrecv+nrecv_fc(2)
         nsend=nsend+nsend_fc(2)
         isize(2)=nxCo_fc(1,2)*nxCo_fc(2,2)*nxCo_fc(3,2)*(nwfluxin)
         recvsize=recvsize+nrecv_fc(2)*isize(2)
         sendsize=sendsize+nsend_fc(2)*isize(2)
         if(stagger_grid) then
           ! This does not consider the 'coarse corner' case
           nxCo1=ixGhi1/2-nghostcells+1;nxCo2=1;nxCo3=ixGhi3/2-nghostcells+1;
           isize_stg(2)=nxCo1*nxCo2*nxCo3*(3-1)
           ! the whole size is used (cell centered and staggered)
           isize(2)=isize(2)+isize_stg(2)
           recvsize=recvsize+nrecv_fc(2)*isize_stg(2)
           sendsize=sendsize+nsend_fc(2)*isize_stg(2)
           ! Coarse corner case
           nrecv_ct=nrecv_ct+nrecv_cc(2)
           nsend_ct=nsend_ct+nsend_cc(2)
           recvsize_cc=recvsize_cc+nrecv_cc(2)*isize_stg(2)
           sendsize_cc=sendsize_cc+nsend_cc(2)*isize_stg(2)
         end if

         case (3)
         nrecv=nrecv+nrecv_fc(3)
         nsend=nsend+nsend_fc(3)
         isize(3)=nxCo_fc(1,3)*nxCo_fc(2,3)*nxCo_fc(3,3)*(nwfluxin)
         recvsize=recvsize+nrecv_fc(3)*isize(3)
         sendsize=sendsize+nsend_fc(3)*isize(3)
         if(stagger_grid) then
           ! This does not consider the 'coarse corner' case
           nxCo1=ixGhi1/2-nghostcells+1;nxCo2=ixGhi2/2-nghostcells+1;nxCo3=1;
           isize_stg(3)=nxCo1*nxCo2*nxCo3*(3-1)
           ! the whole size is used (cell centered and staggered)
           isize(3)=isize(3)+isize_stg(3)
           recvsize=recvsize+nrecv_fc(3)*isize_stg(3)
           sendsize=sendsize+nsend_fc(3)*isize_stg(3)
           ! Coarse corner case
           nrecv_ct=nrecv_ct+nrecv_cc(3)
           nsend_ct=nsend_ct+nsend_cc(3)
           recvsize_cc=recvsize_cc+nrecv_cc(3)*isize_stg(3)
           sendsize_cc=sendsize_cc+nsend_cc(3)*isize_stg(3)
         end if

       end select
     end do

     fc_recvsize = recvsize
     fc_sendsize = sendsize

     ! Grow the buffers to fit, never shrink.  init_comm_fix_conserve runs
     ! every timestep and the number of coarse-fine interfaces oscillates as
     ! the mesh moves, so reallocating whenever the size merely *differs* meant
     ! a device free and malloc on most steps: measured on four ranks in
     ! tests/hd/spherical/blast_amr.par, 164 of 216 calls on the busiest rank.
     ! Those are synchronising on CUDA.  Holding the high-water mark costs
     ! almost nothing here - it was 4640 doubles, 37 kB, in that same run.
     !
     ! Everything downstream is sized by recvsize / sendsize rather than by
     ! size(buffer), so a buffer larger than the exchange is harmless; the two
     ! host/device transfers are sliced for exactly this reason.
     if (allocated(recvbuffer)) then
       if (recvsize > size(recvbuffer)) then
         !$acc exit data delete(recvbuffer)
         deallocate(recvbuffer)
         allocate(recvbuffer(recvsize))
         !$acc enter data create(recvbuffer)
       end if
     else
       allocate(recvbuffer(max(recvsize,1)))
       !$acc enter data create(recvbuffer)
     end if

     ! Key -> offset table for the incoming chunks, read by fix_conserve.
     ! Allocated once: max_blocks is fixed, and so therefore is the key range.
     if (.not.allocated(ibuf_offset)) then
       allocate(ibuf_offset(4**3*max_blocks))
       ibuf_offset = -1
       !$acc enter data copyin(ibuf_offset)
     end if

     if (allocated(sendbuffer)) then
       if (sendsize > size(sendbuffer)) then
         !$acc exit data delete(sendbuffer)
         deallocate(sendbuffer)
         allocate(sendbuffer(sendsize))
         !$acc enter data create(sendbuffer)
       end if
     else
       allocate(sendbuffer(max(sendsize,1)))
       !$acc enter data create(sendbuffer)
     end if

     ! Per-chunk working set, sized by nsend / nrecv, and grown on the same
     ! never-shrink rule as the buffers above - it churned in lockstep with
     ! them, and three of these carry a device copy. Only 1:nsend / 1:nrecv is
     ! ever read, so a longer array is harmless. npe bounds the peer lists
     ! trivially.
     if (allocated(snd_ibuf)) then
       if (max(nsend,1) > size(snd_ibuf)) then
         !$acc exit data delete(snd_igrid, snd_ibuf, snd_group)
         deallocate(snd_igrid, snd_ibuf, snd_key, snd_dest, snd_dims, snd_group)
       end if
     end if
     if (.not.allocated(snd_ibuf)) then
       allocate(snd_igrid(max(nsend,1)), snd_ibuf(max(nsend,1)),&
          snd_key(max(nsend,1)), snd_dest(max(nsend,1)),&
          snd_dims(max(nsend,1)), snd_group(max(nsend,1),2,3))
       !$acc enter data create(snd_igrid, snd_ibuf, snd_group)
     end if

     if (allocated(rcv_key)) then
       if (max(nrecv,1) > size(rcv_key)) deallocate(rcv_key, rcv_src,&
          rcv_dims, rcv_off)
     end if
     if (.not.allocated(rcv_key)) then
       allocate(rcv_key(max(nrecv,1)), rcv_src(max(nrecv,1)),&
          rcv_dims(max(nrecv,1)), rcv_off(max(nrecv,1)))
     end if

     if (.not.allocated(send_pe)) then
       allocate(send_pe(npe), send_pe_off(npe), send_pe_len(npe))
       allocate(recv_pe(npe), recv_pe_off(npe), recv_pe_len(npe))
     end if

     ! Only 1:n_recv_pe / 1:n_send_pe of these is ever posted or waited on, and
     ! that is bounded by the chunk count, so growing without shrinking is safe.
     if (allocated(fc_recvreq)) then
       if (nrecv > size(fc_recvreq)) then
         deallocate(fc_recvreq, fc_recvstat)
         allocate(fc_recvstat(MPI_STATUS_SIZE,max(nrecv,1)),&
            fc_recvreq(max(nrecv,1)))
       end if
     else
       allocate(fc_recvstat(MPI_STATUS_SIZE,max(nrecv,1)),&
          fc_recvreq(max(nrecv,1)))
     end if

     if (allocated(fc_sendreq)) then
       if (nsend > size(fc_sendreq)) then
         deallocate(fc_sendreq, fc_sendstat)
         allocate(fc_sendstat(MPI_STATUS_SIZE,max(nsend,1)),&
            fc_sendreq(max(nsend,1)))
       end if
     else
       allocate(fc_sendstat(MPI_STATUS_SIZE,max(nsend,1)),&
          fc_sendreq(max(nsend,1)))
     end if

     if(stagger_grid) then

       if (allocated(recvbuffer_cc)) then
         if (recvsize_cc /= size(recvbuffer_cc)) then
           deallocate(recvbuffer_cc)
           allocate(recvbuffer_cc(recvsize_cc))
         end if
       else
         allocate(recvbuffer_cc(recvsize_cc))
       end if

       if (allocated(cc_recvreq)) then
         if (nrecv_ct /= size(cc_recvreq)) then
           deallocate(cc_recvreq, cc_recvstat)
           allocate(cc_recvstat(MPI_STATUS_SIZE,nrecv_ct),&
               cc_recvreq(nrecv_ct))
         end if
       else
         allocate(cc_recvstat(MPI_STATUS_SIZE,nrecv_ct), cc_recvreq(nrecv_ct))
       end if

       if (allocated(sendbuffer_cc)) then
         if (sendsize_cc /= size(sendbuffer_cc)) then
           deallocate(sendbuffer_cc)
           allocate(sendbuffer_cc(sendsize_cc))
         end if
       else
         allocate(sendbuffer_cc(sendsize_cc))
       end if

       if (allocated(cc_sendreq)) then
         if (nsend_ct /= size(cc_sendreq)) then
           deallocate(cc_sendreq, cc_sendstat)
           allocate(cc_sendstat(MPI_STATUS_SIZE,nsend_ct),&
               cc_sendreq(nsend_ct))
         end if
       else
         allocate(cc_sendstat(MPI_STATUS_SIZE,nsend_ct), cc_sendreq(nsend_ct))
       end if
     end if

     !$acc update device(isize, nxCo_fc, nwflux_fc)


     call build_message_layout(idimmin,idimmax)

   end subroutine init_comm_fix_conserve

   !> Enumerate every chunk this rank sends and receives, and lay each peer's
   !> chunks out as one contiguous run of the exchange buffer.  Fills snd_igrid
   !> / snd_ibuf / snd_group for the pack kernels, ibuf_offset for
   !> fix_conserve, and the send_pe_* / recv_pe_* runs for the MPI calls.
   subroutine build_message_layout(idimmin,idimmax)
     use mod_global_parameters
     use mod_comm_lib, only: mpistop
     use mod_msg_layout, only: layout_runs

     integer, intent(in) :: idimmin,idimmax

     integer :: iigrid, igrid, idims, iside, i1,i2,i3, ic1,ic2,ic3
     integer :: inc1,inc2,inc3, ineighbor, ipe_neighbor, n, m, k

     ! Outgoing: this rank's fine blocks that face a coarse neighbour
     ! elsewhere.  The key is built from the *receiver's* block index, which is
     ! what neighbor() holds, so it matches the key the receiver derives from
     ! its own igrid below.
     n = 0
     snd_ngroup = 0
     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims=idimmin,idimmax
         do iside=1,2
           i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
           i3=kr(3,idims)*(2*iside-3);
           if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
           if (neighbor_type(i1,i2,i3,igrid)/=neighbor_coarse) cycle
           ineighbor   =neighbor(1,i1,i2,i3,igrid)
           ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
           if (ipe_neighbor==mype) cycle
           ic1=1+modulo(node(pig1_,igrid)-1,2)
           ic2=1+modulo(node(pig2_,igrid)-1,2)
           ic3=1+modulo(node(pig3_,igrid)-1,2);
           inc1=-2*i1+ic1;inc2=-2*i2+ic2;inc3=-2*i3+ic3;
           n = n + 1
           if (n > nsend) call mpistop(&
              "build_message_layout: more sends than nsend_fc counted")
           snd_igrid(n) = igrid
           snd_dest(n)  = ipe_neighbor
           snd_dims(n)  = idims
           snd_key(n)   = 4**3*(ineighbor-1)+inc1+4*inc2+16*inc3
           snd_ngroup(iside,idims) = snd_ngroup(iside,idims) + 1
           snd_group(snd_ngroup(iside,idims),iside,idims) = n
         end do
       end do
     end do
     if (n /= nsend) call mpistop(&
        "build_message_layout: fewer sends than nsend_fc counted")

     ! the key fixes the chunk's size: inc_d lies in {1,2} exactly when i_d is
     ! zero, so the one direction with inc_d in {0,3} is idims - but the walk
     ! above already carried idims along, so hand layout_runs the sizes directly
     call layout_runs(n, snd_dest, snd_key, isize(snd_dims(1:n)), snd_ibuf,&
        n_send_pe, send_pe, send_pe_off, send_pe_len)
     !$acc update device(snd_igrid, snd_ibuf, snd_group)

     ! Incoming: this rank's coarse blocks that face fine children elsewhere.
     m = 0
     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims=idimmin,idimmax
         do iside=1,2
           i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
           i3=kr(3,idims)*(2*iside-3);
           if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
           if (neighbor_type(i1,i2,i3,igrid)/=neighbor_fine) cycle
           do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
             inc3=2*i3+ic3
           do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
             inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
             inc1=2*i1+ic1
             ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
             if (ipe_neighbor==mype) cycle
             m = m + 1
             if (m > nrecv) call mpistop(&
                "build_message_layout: more receives than nrecv_fc counted")
             rcv_src(m)  = ipe_neighbor
             rcv_dims(m) = idims
             rcv_key(m)  = 4**3*(igrid-1)+inc1+4*inc2+16*inc3
           end do
           end do
           end do
         end do
       end do
     end do
     if (m /= nrecv) call mpistop(&
        "build_message_layout: fewer receives than nrecv_fc counted")

     call layout_runs(m, rcv_src, rcv_key, isize(rcv_dims(1:m)), rcv_off,&
        n_recv_pe, recv_pe, recv_pe_off, recv_pe_len)

     ! Scatter into the key -> offset table fix_conserve reads.  Entries for
     ! keys not in this exchange keep whatever they held; nothing reads them.
     do k = 1, m
       ibuf_offset(rcv_key(k)+1) = rcv_off(k)
     end do
     !$acc update device(ibuf_offset)

   end subroutine build_message_layout

   subroutine recvflux(idimmin,idimmax)
     use mod_global_parameters
     use mod_comm_lib, only: mpistop

     integer, intent(in) :: idimmin,idimmax

     integer :: iigrid, igrid, idims, iside, i1,i2,i3, nxCo1,nxCo2,nxCo3
     integer :: ic1,ic2,ic3, inc1,inc2,inc3, ipe_neighbor
     integer :: pi1,pi2,pi3,mi1,mi2,mi3,ph1,ph2,ph3,mh1,mh2,mh3,idir

     ! One receive per peer.  build_message_layout has already decided where
     ! each peer's run starts and how long it is, and filled ibuf_offset so
     ! fix_conserve can find an individual chunk inside it.  The tag no longer
     ! identifies a chunk - (communicator, source) plus one tag per exchange is
     ! enough now that a peer sends exactly one message - so it only has to
     ! separate two exchanges with different dimension ranges.
     if (n_recv_pe>0) then
       fc_recvreq=MPI_REQUEST_NULL
       itag=idimmin+4*idimmax
#ifndef NOGPUDIRECT
       !$acc host_data use_device(recvbuffer)
#endif
       do irecv=1,n_recv_pe
         call mpi_irecv_wrapper(recvbuffer(recv_pe_off(irecv)),&
            recv_pe_len(irecv),MPI_DOUBLE_PRECISION,recv_pe(irecv),itag,icomm,&
            fc_recvreq(irecv),ierrmpi)
       end do
#ifndef NOGPUDIRECT
       !$acc end host_data
#endif
     end if

     if(stagger_grid) then
     ! receive corners
       if (nrecv_ct>0) then
         cc_recvreq=MPI_REQUEST_NULL
         ibuf_cc=1
         irecv_cc=0

         do iigrid=1,igridstail; igrid=igrids(iigrid);
           do idims= idimmin,idimmax
             do iside=1,2
               i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
               i3=kr(3,idims)*(2*iside-3);
               ! Check if there are special corners
               ! (Coarse block diagonal to a fine block)
               ! If there are, receive.
               ! Tags are calculated in the same way as for
               ! normal fluxes, but should not overlap because
               ! inc^D are different
               if (neighbor_type(i1,i2,i3,igrid)==3) then
                 do idir=idims+1,ndim
                   pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
                   mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
                   ph1=pi1-kr(idims,1)*(2*iside-3)
                   ph2=pi2-kr(idims,2)*(2*iside-3)
                   ph3=pi3-kr(idims,3)*(2*iside-3);
                   mh1=mi1-kr(idims,1)*(2*iside-3)
                   mh2=mi2-kr(idims,2)*(2*iside-3)
                   mh3=mi3-kr(idims,3)*(2*iside-3);

                   if (neighbor_type(pi1,pi2,pi3,&
                      igrid)==4.and.neighbor_type(ph1,ph2,ph3,&
                      igrid)==3.and.neighbor_pole(pi1,pi2,pi3,igrid)==0) then
                      ! Loop on children (several in 3D)
                    do ic3=1+int((1-pi3)/2),2-int((1+pi3)/2)
                       inc3=2*pi3+ic3
                    do ic2=1+int((1-pi2)/2),2-int((1+pi2)/2)
                       inc2=2*pi2+ic2
                    do ic1=1+int((1-pi1)/2),2-int((1+pi1)/2)
                       inc1=2*pi1+ic1
                       ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                       if (mype/=ipe_neighbor) then
                         irecv_cc=irecv_cc+1
                         itag_cc=4**3*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                            inc3*4**(3-1)
                         call mpi_irecv_wrapper(recvbuffer_cc(ibuf_cc),&
                            isize_stg(idims),MPI_DOUBLE_PRECISION,ipe_neighbor,&
                            itag_cc,icomm,cc_recvreq(irecv_cc),ierrmpi)
                         ibuf_cc=ibuf_cc+isize_stg(idims)
                       end if
                    end do
                    end do
                    end do
                   end if

                   if (neighbor_type(mi1,mi2,mi3,&
                      igrid)==4.and.neighbor_type(mh1,mh2,mh3,&
                      igrid)==3.and.neighbor_pole(mi1,mi2,mi3,igrid)==0) then
                      ! Loop on children (several in 3D)
                    do ic3=1+int((1-mi3)/2),2-int((1+mi3)/2)
                        inc3=2*mi3+ic3
                    do ic2=1+int((1-mi2)/2),2-int((1+mi2)/2)
                        inc2=2*mi2+ic2
                    do ic1=1+int((1-mi1)/2),2-int((1+mi1)/2)
                        inc1=2*mi1+ic1
                       ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                       if (mype/=ipe_neighbor) then
                         irecv_cc=irecv_cc+1
                         itag_cc=4**3*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
                            inc3*4**(3-1)
                         call mpi_irecv_wrapper(recvbuffer_cc(ibuf_cc),&
                            isize_stg(idims),MPI_DOUBLE_PRECISION,ipe_neighbor,&
                            itag_cc,icomm,cc_recvreq(irecv_cc),ierrmpi)
                         ibuf_cc=ibuf_cc+isize_stg(idims)
                       end if
                    end do
                    end do
                    end do
                   end if
                 end do
               end if
             end do
           end do
         end do
       end if
     end if ! end if stagger grid

   end subroutine recvflux

   !> Pack the outgoing flux chunks and post one message per peer.  Where each
   !> chunk goes was decided by build_message_layout; nothing here walks the
   !> tree.
   subroutine sendflux(idimmin,idimmax)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax

     integer :: idims, iside, ix1,ix2,ix3, iw
     integer :: k, ng, imsg

     fc_sendreq = MPI_REQUEST_NULL

     ! One kernel per (iside, idims) group rather than one per chunk: pflux's
     ! slot is then loop-invariant, and a derived-type component holding an
     ! allocatable array cannot be selected by a device-side index.
     do idims = idimmin,idimmax
       do iside = 1,2
         ng = snd_ngroup(iside,idims)
         if (ng == 0) cycle
         select case (idims)
         case (1)
           !$acc parallel loop gang private(imsg) default(present)
           do k = 1,ng
             imsg = snd_group(k,iside,1)
             !$acc loop vector collapse(3)
             do iw=1,nwflux_fc
               do ix3=1,nxCo_fc(3,1)
                 do ix2=1,nxCo_fc(2,1)
                   sendbuffer(snd_ibuf(imsg)+(ix2-1)+(ix3-1)*nxCo_fc(2,1) &
                      +(iw-1)*nxCo_fc(2,1)*nxCo_fc(3,1)) = &
                      pflux(iside,1)%flux(1,ix2,ix3,iw,snd_igrid(imsg))
                 end do
               end do
             end do
           end do
         case (2)
           !$acc parallel loop gang private(imsg) default(present)
           do k = 1,ng
             imsg = snd_group(k,iside,2)
             !$acc loop vector collapse(3)
             do iw=1,nwflux_fc
               do ix3=1,nxCo_fc(3,2)
                 do ix1=1,nxCo_fc(1,2)
                   sendbuffer(snd_ibuf(imsg)+(ix1-1)+(ix3-1)*nxCo_fc(1,2) &
                      +(iw-1)*nxCo_fc(1,2)*nxCo_fc(3,2)) = &
                      pflux(iside,2)%flux(ix1,1,ix3,iw,snd_igrid(imsg))
                 end do
               end do
             end do
           end do
         case (3)
           !$acc parallel loop gang private(imsg) default(present)
           do k = 1,ng
             imsg = snd_group(k,iside,3)
             !$acc loop vector collapse(3)
             do iw=1,nwflux_fc
               do ix2=1,nxCo_fc(2,3)
                 do ix1=1,nxCo_fc(1,3)
                   sendbuffer(snd_ibuf(imsg)+(ix1-1)+(ix2-1)*nxCo_fc(1,3) &
                      +(iw-1)*nxCo_fc(1,3)*nxCo_fc(2,3)) = &
                      pflux(iside,3)%flux(ix1,ix2,1,iw,snd_igrid(imsg))
                 end do
               end do
             end do
           end do
         end select
       end do
     end do

     ! One transfer for the whole buffer, then one Isend per peer.  Each peer's
     ! chunks are contiguous and in the order that peer expects them, so a
     ! single message carries all of them and the tag carries no block index.
#ifdef NOGPUDIRECT
     if (n_send_pe > 0) then
       ! sliced, not whole-array: the buffer is grown and not shrunk, so
       ! size(sendbuffer) can exceed what this exchange actually uses
       !$acc update host(sendbuffer(1:fc_sendsize))
     end if
#else
     !$acc host_data use_device(sendbuffer)
#endif
     itag=idimmin+4*idimmax
     do k = 1,n_send_pe
       call mpi_isend_wrapper(sendbuffer(send_pe_off(k)),send_pe_len(k),&
          MPI_DOUBLE_PRECISION,send_pe(k),itag,icomm,fc_sendreq(k),ierrmpi)
     end do
#ifndef NOGPUDIRECT
     !$acc end host_data
#endif

   end subroutine sendflux

   !  ------------------------------------------------------------------------
   !  Staggered-grid (constrained transport) fragments from the pre-aggregation
   !  sendflux, kept for reference.
   !
   !  AGILE has no stagger_grid support: pflux's `edge` component is commented
   !  out of the fluxalloc type, store_edge and fix_edges are inert, and the
   !  corner buffers sendbuffer_cc / cc_sendreq are sized but never written.
   !  What follows is upstream MPI-AMRVAC's send side for that path as it stood
   !  in this file, and is the most useful starting point if constrained
   !  transport is implemented here. The matching receive side is still present,
   !  and live, in recvflux's if(stagger_grid) block.
   !
   !  Two things need adapting before any of it can be used again:
   !    - it packs and posts one message at a time, inline in a tree walk that
   !      sendflux no longer has. The equivalent now is to enumerate the corner
   !      chunks in build_message_layout and give them their own per-peer runs,
   !      exactly as the face fluxes get theirs;
   !    - isize(idims) used to absorb isize_stg(idims), so reviving the edge
   !      payload means sizing and laying that part out too.
   !
   !  Per direction: first the staggered variant of the face-flux send, then
   !  the corner ('coarse corner') exchange for a fine block surrounded by
   !  coarse ones.
   !  ------------------------------------------------------------------------

   !  ---- direction 1 ----

                 !!  ibuf_send_next=ibuf_send+isize(1)
                 !!  sendbuffer(ibuf_send:ibuf_send_next-isize_stg(1)-&
                 !!     1)=reshape(pflux(iside,1,igrid)%flux,&
                 !!     (/isize(1)-isize_stg(1)/))

                 !!  sendbuffer(ibuf_send_next-isize_stg(1):ibuf_send_next-&
                 !!     1)=reshape(pflux(iside,1,igrid)%edge,(/isize_stg(1)/))
                 !!  call mpi_isend_wrapper(sendbuffer(ibuf_send),isize(1),&
                 !!      MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                 !!     fc_sendreq(isend),ierrmpi)
                 !!  ibuf_send=ibuf_send_next

               !!if(stagger_grid) then
               !!  ! If we are in a fine block surrounded by coarse blocks
               !!  do idir=idims+1,ndim
               !!    pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
               !!    mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
               !!    ph1=pi1-kr(idims,1)*(2*iside-3)
               !!    ph2=pi2-kr(idims,2)*(2*iside-3)
               !!    ph3=pi3-kr(idims,3)*(2*iside-3);
               !!    mh1=mi1-kr(idims,1)*(2*iside-3)
               !!    mh2=mi2-kr(idims,2)*(2*iside-3)
               !!    mh3=mi3-kr(idims,3)*(2*iside-3);

               !!    if (neighbor_type(pi1,pi2,pi3,&
               !!       igrid)==2.and.neighbor_type(ph1,ph2,ph3,&
               !!       igrid)==2.and.mype/=neighbor(2,pi1,pi2,pi3,&
               !!       igrid).and.neighbor_pole(pi1,pi2,pi3,igrid)==0) then
               !!      ! Get relative position in the grid for tags
               !!      ineighbor=neighbor(1,pi1,pi2,pi3,igrid)
               !!      ipe_neighbor=neighbor(2,pi1,pi2,pi3,igrid)
               !!      ic1=1+modulo(node(pig1_,igrid)-1,2)
               !!      ic2=1+modulo(node(pig2_,igrid)-1,2)
               !!      ic3=1+modulo(node(pig3_,igrid)-1,2);
               !!      inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
               !!      itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
               !!         inc3*4**(3-1)
               !!      ! Reshape to buffer and send
               !!      isend_cc=isend_cc+1
               !!      ibuf_cc_send_next=ibuf_cc_send+isize_stg(1)
               !!      sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
               !!         1)=reshape(pflux(iside,1,igrid)%edge,&
               !!         shape=(/isize_stg(1)/))
               !!      call mpi_isend_wrapper(sendbuffer_cc(ibuf_cc_send),isize_stg(1),&
               !!         MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
               !!         cc_sendreq(isend_cc),ierrmpi)
               !!      ibuf_cc_send=ibuf_cc_send_next
               !!    end if

               !!    if (neighbor_type(mi1,mi2,mi3,&
               !!       igrid)==2.and.neighbor_type(mh1,mh2,mh3,&
               !!       igrid)==2.and.mype/=neighbor(2,mi1,mi2,mi3,&
               !!       igrid).and.neighbor_pole(mi1,mi2,mi3,igrid)==0) then
               !!      ! Get relative position in the grid for tags
               !!      ineighbor=neighbor(1,mi1,mi2,mi3,igrid)
               !!      ipe_neighbor=neighbor(2,mi1,mi2,mi3,igrid)
               !!      ic1=1+modulo(node(pig1_,igrid)-1,2)
               !!      ic2=1+modulo(node(pig2_,igrid)-1,2)
               !!      ic3=1+modulo(node(pig3_,igrid)-1,2);
               !!      inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
               !!      inc1=-2*mi1+ic1;inc2=-2*mi2+ic2;inc3=-2*mi3+ic3;
               !!      itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
               !!         inc3*4**(3-1)
               !!      ! Reshape to buffer and send
               !!      isend_cc=isend_cc+1
               !!      ibuf_cc_send_next=ibuf_cc_send+isize_stg(1)
               !!      sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
               !!         1)=reshape(pflux(iside,1,igrid)%edge,&
               !!         shape=(/isize_stg(1)/))
               !!      call mpi_isend_wrapper(sendbuffer_cc(ibuf_cc_send),isize_stg(1),&
               !!         MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
               !!         cc_sendreq(isend_cc),ierrmpi)
               !!      ibuf_cc_send=ibuf_cc_send_next
               !!    end if
               !!  end do
               !!end if ! end if stagger grid

   !  ---- direction 2 ----

                 !!  ibuf_send_next=ibuf_send+isize(2)
                 !!  sendbuffer(ibuf_send:ibuf_send_next-isize_stg(2)-&
                 !!     1)=reshape(pflux(iside,2,igrid)%flux,&
                 !!     (/isize(2)-isize_stg(2)/))

                 !!  sendbuffer(ibuf_send_next-isize_stg(2):ibuf_send_next-&
                 !!     1)=reshape(pflux(iside,2,igrid)%edge,(/isize_stg(2)/))
                 !!  call mpi_isend_wrapper(sendbuffer(ibuf_send),isize(2),&
                 !!      MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                 !!     fc_sendreq(isend),ierrmpi)
                 !!  ibuf_send=ibuf_send_next

               !!if(stagger_grid) then
               !!  ! If we are in a fine block surrounded by coarse blocks
               !!  do idir=idims+1,ndim
               !!    pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
               !!    mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
               !!    ph1=pi1-kr(idims,1)*(2*iside-3)
               !!    ph2=pi2-kr(idims,2)*(2*iside-3)
               !!    ph3=pi3-kr(idims,3)*(2*iside-3);
               !!    mh1=mi1-kr(idims,1)*(2*iside-3)
               !!    mh2=mi2-kr(idims,2)*(2*iside-3)
               !!    mh3=mi3-kr(idims,3)*(2*iside-3);

               !!    if (neighbor_type(pi1,pi2,pi3,&
               !!       igrid)==2.and.neighbor_type(ph1,ph2,ph3,&
               !!       igrid)==2.and.mype/=neighbor(2,pi1,pi2,pi3,&
               !!       igrid).and.neighbor_pole(pi1,pi2,pi3,igrid)==0) then
               !!      ! Get relative position in the grid for tags
               !!      ineighbor=neighbor(1,pi1,pi2,pi3,igrid)
               !!      ipe_neighbor=neighbor(2,pi1,pi2,pi3,igrid)
               !!      ic1=1+modulo(node(pig1_,igrid)-1,2)
               !!      ic2=1+modulo(node(pig2_,igrid)-1,2)
               !!      ic3=1+modulo(node(pig3_,igrid)-1,2);
               !!      inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
               !!      itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
               !!         inc3*4**(3-1)
               !!      ! Reshape to buffer and send
               !!      isend_cc=isend_cc+1
               !!      ibuf_cc_send_next=ibuf_cc_send+isize_stg(2)
               !!      sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
               !!         1)=reshape(pflux(iside,2,igrid)%edge,&
               !!         shape=(/isize_stg(2)/))
               !!      call mpi_isend_wrapper(sendbuffer_cc(ibuf_cc_send),isize_stg(2),&
               !!         MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
               !!         cc_sendreq(isend_cc),ierrmpi)
               !!      ibuf_cc_send=ibuf_cc_send_next
               !!    end if

               !!    if (neighbor_type(mi1,mi2,mi3,&
               !!       igrid)==2.and.neighbor_type(mh1,mh2,mh3,&
               !!       igrid)==2.and.mype/=neighbor(2,mi1,mi2,mi3,&
               !!       igrid).and.neighbor_pole(mi1,mi2,mi3,igrid)==0) then
               !!      ! Get relative position in the grid for tags
               !!      ineighbor=neighbor(1,mi1,mi2,mi3,igrid)
               !!      ipe_neighbor=neighbor(2,mi1,mi2,mi3,igrid)
               !!      ic1=1+modulo(node(pig1_,igrid)-1,2)
               !!      ic2=1+modulo(node(pig2_,igrid)-1,2)
               !!      ic3=1+modulo(node(pig3_,igrid)-1,2);
               !!      inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
               !!      inc1=-2*mi1+ic1;inc2=-2*mi2+ic2;inc3=-2*mi3+ic3;
               !!      itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
               !!         inc3*4**(3-1)
               !!      ! Reshape to buffer and send
               !!      isend_cc=isend_cc+1
               !!      ibuf_cc_send_next=ibuf_cc_send+isize_stg(2)
               !!      sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
               !!         1)=reshape(pflux(iside,2,igrid)%edge,&
               !!         shape=(/isize_stg(2)/))
               !!      call mpi_isend_wrapper(sendbuffer_cc(ibuf_cc_send),isize_stg(2),&
               !!         MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
               !!         cc_sendreq(isend_cc),ierrmpi)
               !!      ibuf_cc_send=ibuf_cc_send_next
               !!    end if
               !!  end do
               !!end if ! end if stagger grid

   !  ---- direction 3 ----

                 !!  ibuf_send_next=ibuf_send+isize(3)
                 !!  sendbuffer(ibuf_send:ibuf_send_next-isize_stg(3)-&
                 !!     1)=reshape(pflux(iside,3,igrid)%flux,&
                 !!     (/isize(3)-isize_stg(3)/))

                 !!  sendbuffer(ibuf_send_next-isize_stg(3):ibuf_send_next-&
                 !!     1)=reshape(pflux(iside,3,igrid)%edge,(/isize_stg(3)/))
                 !!  call mpi_isend_wrapper(sendbuffer(ibuf_send),isize(3),&
                 !!      MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                 !!     fc_sendreq(isend),ierrmpi)
                 !!  ibuf_send=ibuf_send_next

               !!if(stagger_grid) then
               !!  ! If we are in a fine block surrounded by coarse blocks
               !!  do idir=idims+1,ndim
               !!    pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
               !!    mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
               !!    ph1=pi1-kr(idims,1)*(2*iside-3)
               !!    ph2=pi2-kr(idims,2)*(2*iside-3)
               !!    ph3=pi3-kr(idims,3)*(2*iside-3);
               !!    mh1=mi1-kr(idims,1)*(2*iside-3)
               !!    mh2=mi2-kr(idims,2)*(2*iside-3)
               !!    mh3=mi3-kr(idims,3)*(2*iside-3);

               !!    if (neighbor_type(pi1,pi2,pi3,&
               !!       igrid)==2.and.neighbor_type(ph1,ph2,ph3,&
               !!       igrid)==2.and.mype/=neighbor(2,pi1,pi2,pi3,&
               !!       igrid).and.neighbor_pole(pi1,pi2,pi3,igrid)==0) then
               !!      ! Get relative position in the grid for tags
               !!      ineighbor=neighbor(1,pi1,pi2,pi3,igrid)
               !!      ipe_neighbor=neighbor(2,pi1,pi2,pi3,igrid)
               !!      ic1=1+modulo(node(pig1_,igrid)-1,2)
               !!      ic2=1+modulo(node(pig2_,igrid)-1,2)
               !!      ic3=1+modulo(node(pig3_,igrid)-1,2);
               !!      inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
               !!      itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
               !!         inc3*4**(3-1)
               !!      ! Reshape to buffer and send
               !!      isend_cc=isend_cc+1
               !!      ibuf_cc_send_next=ibuf_cc_send+isize_stg(3)
               !!      sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
               !!         1)=reshape(pflux(iside,3,igrid)%edge,&
               !!         shape=(/isize_stg(3)/))
               !!      call mpi_isend_wrapper(sendbuffer_cc(ibuf_cc_send),isize_stg(3),&
               !!         MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
               !!         cc_sendreq(isend_cc),ierrmpi)
               !!      ibuf_cc_send=ibuf_cc_send_next
               !!    end if

               !!    if (neighbor_type(mi1,mi2,mi3,&
               !!       igrid)==2.and.neighbor_type(mh1,mh2,mh3,&
               !!       igrid)==2.and.mype/=neighbor(2,mi1,mi2,mi3,&
               !!       igrid).and.neighbor_pole(mi1,mi2,mi3,igrid)==0) then
               !!      ! Get relative position in the grid for tags
               !!      ineighbor=neighbor(1,mi1,mi2,mi3,igrid)
               !!      ipe_neighbor=neighbor(2,mi1,mi2,mi3,igrid)
               !!      ic1=1+modulo(node(pig1_,igrid)-1,2)
               !!      ic2=1+modulo(node(pig2_,igrid)-1,2)
               !!      ic3=1+modulo(node(pig3_,igrid)-1,2);
               !!      inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;inc3=-2*pi3+ic3;
               !!      inc1=-2*mi1+ic1;inc2=-2*mi2+ic2;inc3=-2*mi3+ic3;
               !!      itag_cc=4**3*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)+&
               !!         inc3*4**(3-1)
               !!      ! Reshape to buffer and send
               !!      isend_cc=isend_cc+1
               !!      ibuf_cc_send_next=ibuf_cc_send+isize_stg(3)
               !!      sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
               !!         1)=reshape(pflux(iside,3,igrid)%edge,&
               !!         shape=(/isize_stg(3)/))
               !!      call mpi_isend_wrapper(sendbuffer_cc(ibuf_cc_send),isize_stg(3),&
               !!         MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
               !!         cc_sendreq(isend_cc),ierrmpi)
               !!      ibuf_cc_send=ibuf_cc_send_next
               !!    end if
               !!  end do
               !!end if ! end if stagger grid

   !> Correct the coarse cells abutting a refinement boundary: take out the
   !> coarse face flux the update used and put back the sum of the fine ones.
   !>
   !> The stored flux is Cartesian qdt*f/dx or curvilinear qdt*f*A, see the
   !> flux-fixing stores in mod_finite_volume.  The Cartesian form is already a
   !> per-volume quantity but of the *fine* cell, so it is rescaled by
   !> CoFiratio; the curvilinear form is extensive and is divided here by the
   !> coarse cell's own volume, with no ratio, because the four fine face areas
   !> already sum to the coarse face's.
   subroutine fix_conserve(psb,idimmin,idimmax,nw0,nwfluxin)
     use mod_global_parameters
     use mod_comm_lib, only: mpistop

     integer, intent(in) :: idimmin,idimmax, nw0, nwfluxin
     type(state) :: psb(max_blocks)

     integer :: iigrid, igrid, idims, iside, iotherside, i1,i2,i3, ic1,ic2,ic3,&
         inc1,inc2,inc3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
     integer :: ix1, ix2, ix3 !JESSE ADDED
     integer :: nxCo1,nxCo2,nxCo3, iw, ix, ipe_neighbor, ineighbor, nbuf,&
         ibufnext, nw1
#:if GEOM == 'Cartesian'
     double precision :: CoFiratio
#:endif

     nw1=nw0-1+nwfluxin
#:if GEOM == 'Cartesian'
     ! The flux is divided by volume of fine cell. We need, however,
     ! to divide by volume of coarse cell => muliply by volume ratio
     CoFiratio=one/dble(2**ndim)
#:endif

     if (n_recv_pe>0) then
       call MPI_WAITALL(n_recv_pe,fc_recvreq,fc_recvstat,ierrmpi)
       ! With one aggregated message per peer, a disagreement about which
       ! chunks the run holds no longer shows up as an MPI mismatch, so check
       ! the arrival lengths. This catches any difference in the *set* of
       ! chunks the two ranks enumerated; a difference in their *order* cannot
       ! arise, since both sort the same run by the same key.
       do irecv=1,n_recv_pe
         call MPI_GET_COUNT(fc_recvstat(:,irecv),MPI_DOUBLE_PRECISION,nbuf,&
            ierrmpi)
         if (nbuf /= recv_pe_len(irecv)) call mpistop(&
            "fix_conserve: flux message length disagrees with the layout")
       end do
#ifdef NOGPUDIRECT
       ! Without GPU-direct the IRECVs landed in host memory; the unpack below
       ! runs on the device, so push the payload across.  Sliced, not
       ! whole-array: the buffer is grown and not shrunk.
       !$acc update device(recvbuffer(1:fc_recvsize))
#endif
     end if

     nxCo1=(ixMhi1-ixMlo1+1)/2
     nxCo2=(ixMhi2-ixMlo2+1)/2
     nxCo3=(ixMhi3-ixMlo3+1)/2

     ! for all grids: perform flux update at Coarse-Fine interfaces
     ! Everything assigned per gang is named here rather than left to the
     ! implicit firstprivate OpenACC gives an unlisted scalar in a parallel
     ! region: the rule does cover them, but a half-populated list is
     ! indistinguishable from an oversight.
     !$acc parallel loop gang default(present) &
     !$acc& private(igrid, idims, iside, i1,i2,i3, ix, ic1,ic2,ic3, &
     !$acc&         inc1,inc2,inc3, ineighbor, ipe_neighbor, iotherside, &
     !$acc&         ixmin1,ixmin2,ixmin3, ixmax1,ixmax2,ixmax3, &
     !$acc&         ix1,ix2,ix3, iw)
     do iigrid=1,igridstail
       igrid=igrids(iigrid)

       do idims=idimmin,idimmax
         select case (idims)
           case (1)
           do iside=1,2
             i1=kr(1,1)*(2*iside-3)
             i2=kr(2,1)*(2*iside-3)
             i3=kr(3,1)*(2*iside-3)

             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

             if (neighbor_type(i1,i2,i3,igrid)/=4) cycle

 !opedit: skip over active/passive interface since flux for passive ones is
             ! not computed, keep the buffer counter up to date:
          !   if (.not.neighbor_active(i1,i2,i3,&
          !      igrid).or..not.neighbor_active(0,0,0,igrid) ) then
          !     do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
          !     inc3=2*i3+ic3
          ! do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
          !     inc2=2*i2+ic2
          ! do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
          !     inc1=2*i1+ic1
          !     ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
          !     if (ipe_neighbor/=mype) then
          !       ibufnext=ibuf+isize(1)
          !       ibuf=ibufnext
          !     end if
          !     end do
          ! end do
          ! end do
          !     cycle
          !   end if
             !

             select case (iside)
             case (1)
               ix=ixMlo1
             case (2)
               ix=ixMhi1
             end select

             ! remove coarse flux
#:if GEOM == 'Cartesian'
                !$acc loop vector collapse(3)
                do ix3=ixMlo3,ixMhi3
                  do ix2=ixMlo2,ixMhi2 
                    do iw=1,nwfluxin
                      psb(igrid)%w(ix,ix2,ix3,nw0+iw-1) = &
                        psb(igrid)%w(ix,ix2,ix3,nw0+iw-1) - &
                        pflux(iside,1)%flux(1,ix2-nghostcells,ix3-nghostcells,&
                                            iw,igrid)
                    end do
                  end do
                end do
#:else
                !$acc loop vector collapse(3)
                do ix3=ixMlo3,ixMhi3
                  do ix2=ixMlo2,ixMhi2
                    do iw=1,nwfluxin
                      psb(igrid)%w(ix,ix2,ix3,nw0+iw-1) = &
                        psb(igrid)%w(ix,ix2,ix3,nw0+iw-1) - &
                        pflux(iside,1)%flux(1,ix2-nghostcells,ix3-nghostcells,&
                                            iw,igrid) &
                        / bgeo%dvolume(ix,ix2,ix3,igrid)
                    end do
                  end do
                end do
#:endif


             ! add fine flux

             do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
             do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
             do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               ixmin1=ix
               ixmin2=ixMlo2+(ic2-1)*nxCo2
               ixmin3=ixMlo3+(ic3-1)*nxCo3
               ixmax1=ix
               ixmax2=ixmin2-1+nxCo2
               ixmax3=ixmin3-1+nxCo3
               if (ipe_neighbor==mype) then
                 iotherside=3-iside
#:if GEOM == 'Cartesian'
                     !$acc loop vector collapse(3)
                     do ix3=1,nxCo3 
                        do ix2=1,nxCo2 
                          do iw=1,nwfluxin
                             psb(igrid)%w(ix,ixmin2+ix2-1,ixmin3+ix3-1,nw0+iw-1) = &
                              psb(igrid)%w(ix,ixmin2+ix2-1,ixmin3+ix3-1,nw0+iw-1) + &
                              pflux(iotherside,1)%flux(1,ix2,ix3,iw,&
                                ineighbor) * CoFiratio
                          end do
                        end do
                     end do
#:else
                     ! Direction 1, so loop runs over directions 2 and 3
                     !$acc loop vector collapse(3)
                     do ix3=1,nxCo3
                        do ix2=1,nxCo2
                          do iw=1,nwfluxin
                             psb(igrid)%w(ix,ixmin2+ix2-1,ixmin3+ix3-1,nw0+iw-1) = &
                              psb(igrid)%w(ix,ixmin2+ix2-1,ixmin3+ix3-1,nw0+iw-1) + &
                              pflux(iotherside,1)%flux(1,ix2,ix3,iw,&
                                ineighbor) &
                              / bgeo%dvolume(ix,ixmin2+ix2-1,ixmin3+ix3-1,igrid)
                          end do
                        end do
                     end do
#:endif
               !else
               !  if (slab_uniform) then
               !    ibufnext=ibuf+isize(1)
               !    if(stagger_grid) ibufnext=ibufnext-isize_stg(1)
               !    psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
               !       nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !       ixmin3:ixmax3,nw0:nw1)+CoFiratio &
               !       *reshape(source=recvbuffer(ibuf:ibufnext-1),&
               !        shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !       ixmin3:ixmax3,nw0:nw1)))
               !    ibuf=ibuf+isize(1)
               !  else
               !    ibufnext=ibuf+isize(1)
               !    if(stagger_grid) then
               !      nbuf=(isize(1)-isize_stg(1))/nwfluxin
               !    else
               !      nbuf=isize(1)/nwfluxin
               !    end if
               !    do iw=nw0,nw1
               !      psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
               !         iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !         ixmin3:ixmax3,iw) &
               !         +reshape(source=recvbuffer(ibuf:ibufnext-1),&
               !          shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !         ixmin3:ixmax3,iw))) /ps(igrid)%dvolume(ixmin1:ixmax1,&
               !         ixmin2:ixmax2,ixmin3:ixmax3)
               !      ibuf=ibuf+nbuf
               !    end do
               !    ibuf=ibufnext
               !  end if
               else
#:if GEOM == 'Cartesian'
                   ! Two transverse indices plus the variable index: every iteration
                   ! lands in a distinct cell of a distinct variable, and the buffer
                   ! offset is a pure function of the three, so all three collapse.
                   !$acc loop vector collapse(3)
                   do ix3=1,nxCo_fc(3,1)
                     do ix2=1,nxCo_fc(2,1)
                       do iw=1,nwfluxin
                         psb(igrid)%w(ix,ixmin2+ix2-1,ixmin3+ix3-1,nw0+iw-1) = &
                           psb(igrid)%w(ix,ixmin2+ix2-1,ixmin3+ix3-1,nw0+iw-1) + &
                           recvbuffer(ibuf_offset(4**3*(igrid-1)+inc1+4*inc2+16*inc3+1) &
                              +(ix2-1)+(ix3-1)*nxCo_fc(2,1)+(iw-1)*nxCo_fc(2,1)*nxCo_fc(3,1)) * CoFiratio
                       end do
                     end do
                   end do
#:else
                   ! Two transverse indices plus the variable index: every iteration
                   ! lands in a distinct cell of a distinct variable, and the buffer
                   ! offset is a pure function of the three, so all three collapse.
                   !$acc loop vector collapse(3)
                   do ix3=1,nxCo_fc(3,1)
                     do ix2=1,nxCo_fc(2,1)
                       do iw=1,nwfluxin
                         psb(igrid)%w(ix,ixmin2+ix2-1,ixmin3+ix3-1,nw0+iw-1) = &
                           psb(igrid)%w(ix,ixmin2+ix2-1,ixmin3+ix3-1,nw0+iw-1) + &
                           recvbuffer(ibuf_offset(4**3*(igrid-1)+inc1+4*inc2+16*inc3+1) &
                              +(ix2-1)+(ix3-1)*nxCo_fc(2,1)+(iw-1)*nxCo_fc(2,1)*nxCo_fc(3,1)) &
                           / bgeo%dvolume(ix,ixmin2+ix2-1,ixmin3+ix3-1,igrid)
                       end do
                     end do
                   end do
#:endif
               end if
            end do
           end do
           end do
           end do
           case (2)
           do iside=1,2
             i1=kr(1,2)*(2*iside-3)
             i2=kr(2,2)*(2*iside-3)
             i3=kr(3,2)*(2*iside-3)

             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

             if (neighbor_type(i1,i2,i3,igrid)/=4) cycle

 !opedit: skip over active/passive interface since flux for passive ones is
             ! not computed, keep the buffer counter up to date:
            ! if (.not.neighbor_active(i1,i2,i3,&
            !    igrid).or..not.neighbor_active(0,0,0,igrid) ) then
            !   do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
            !   inc3=2*i3+ic3
            !   do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
            !   inc2=2*i2+ic2
            !   do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
            !   inc1=2*i1+ic1
            !   ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
            !   if (ipe_neighbor/=mype) then
            !     ibufnext=ibuf+isize(2)
            !     ibuf=ibufnext
            !   end if
            !   end do
            !  end do
            !  end do
            !   cycle
            ! end if
             !

             select case (iside)
             case (1)
               ix=ixMlo2
             case (2)
               ix=ixMhi2
             end select

             ! remove coarse flux
#:if GEOM == 'Cartesian'
                !$acc loop vector collapse(3)
                do ix3=ixMlo3,ixMhi3
                  do ix1=ixMlo1,ixMhi1 
                    do iw=1,nwfluxin
                      psb(igrid)%w(ix1,ix,ix3,nw0+iw-1) = &
                       psb(igrid)%w(ix1,ix,ix3,nw0+iw-1) - &
                       pflux(iside,2)%flux(ix1-nghostcells,1,ix3-nghostcells,&
                                    iw,igrid)
                    end do
                  end do
                end do
#:else
                !$acc loop vector collapse(3)
                do ix3=ixMlo3,ixMhi3
                  do ix1=ixMlo1,ixMhi1
                    do iw=1,nwfluxin
                      psb(igrid)%w(ix1,ix,ix3,nw0+iw-1) = &
                       psb(igrid)%w(ix1,ix,ix3,nw0+iw-1) - &
                       pflux(iside,2)%flux(ix1-nghostcells,1,ix3-nghostcells,&
                                    iw,igrid) &
                       / bgeo%dvolume(ix1,ix,ix3,igrid)
                    end do
                  end do
                end do
#:endif


             ! add fine flux
             do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
             do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
             do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               ixmin1=ixMlo1+(ic1-1)*nxCo1
               ixmin2=ix
               ixmin3=ixMlo3+(ic3-1)*nxCo3
               ixmax1=ixmin1-1+nxCo1
               ixmax2=ix
               ixmax3=ixmin3-1+nxCo3

               if (ipe_neighbor==mype) then
                 iotherside=3-iside

#:if GEOM == 'Cartesian'
                   !$acc loop vector collapse(3)
                   do ix3=1,nxCo3 
                     do ix1=1,nxCo1 
                       do iw=1,nwfluxin
                         psb(igrid)%w(ixmin1+ix1-1,ix,ixmin3+ix3-1,nw0+iw-1) = &
                           psb(igrid)%w(ixmin1+ix1-1,ix,ixmin3+ix3-1,nw0+iw-1) + &
                           pflux(iotherside,2)%flux(ix1,1,ix3,&
                              iw,ineighbor) * CoFiratio
                       end do
                     end do
                   end do
#:else
                   !$acc loop vector collapse(3)
                   do ix3=1,nxCo3
                     do ix1=1,nxCo1
                       do iw=1,nwfluxin
                         psb(igrid)%w(ixmin1+ix1-1,ix,ixmin3+ix3-1,nw0+iw-1) = &
                           psb(igrid)%w(ixmin1+ix1-1,ix,ixmin3+ix3-1,nw0+iw-1) + &
                           pflux(iotherside,2)%flux(ix1,1,ix3,&
                              iw,ineighbor) &
                           / bgeo%dvolume(ixmin1+ix1-1,ix,ixmin3+ix3-1,igrid)
                       end do
                     end do
                   end do
#:endif
               !else
               !  if (slab_uniform) then
               !    ibufnext=ibuf+isize(2)
               !    if(stagger_grid) ibufnext=ibufnext-isize_stg(2)
               !    psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
               !       nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !       ixmin3:ixmax3,nw0:nw1)+CoFiratio &
               !       *reshape(source=recvbuffer(ibuf:ibufnext-1),&
               !        shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !       ixmin3:ixmax3,nw0:nw1)))
               !    ibuf=ibuf+isize(2)
               !  else
               !    ibufnext=ibuf+isize(2)
               !    if(stagger_grid) then
               !      nbuf=(isize(2)-isize_stg(2))/nwfluxin
               !    else
               !      nbuf=isize(2)/nwfluxin
               !    end if
               !    do iw=nw0,nw1
               !      psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
               !         iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !         ixmin3:ixmax3,iw) &
               !         +reshape(source=recvbuffer(ibuf:ibufnext-1),&
               !          shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !         ixmin3:ixmax3,iw))) /ps(igrid)%dvolume(ixmin1:ixmax1,&
               !         ixmin2:ixmax2,ixmin3:ixmax3)
               !      ibuf=ibuf+nbuf
               !    end do
               !    ibuf=ibufnext
               !  end if
               else
#:if GEOM == 'Cartesian'
                   ! Two transverse indices plus the variable index: every iteration
                   ! lands in a distinct cell of a distinct variable, and the buffer
                   ! offset is a pure function of the three, so all three collapse.
                   !$acc loop vector collapse(3)
                   do ix3=1,nxCo_fc(3,2)
                     do ix1=1,nxCo_fc(1,2)
                       do iw=1,nwfluxin
                         psb(igrid)%w(ixmin1+ix1-1,ix,ixmin3+ix3-1,nw0+iw-1) = &
                           psb(igrid)%w(ixmin1+ix1-1,ix,ixmin3+ix3-1,nw0+iw-1) + &
                           recvbuffer(ibuf_offset(4**3*(igrid-1)+inc1+4*inc2+16*inc3+1) &
                              +(ix1-1)+(ix3-1)*nxCo_fc(1,2)+(iw-1)*nxCo_fc(1,2)*nxCo_fc(3,2)) * CoFiratio
                       end do
                     end do
                   end do
#:else
                   ! Two transverse indices plus the variable index: every iteration
                   ! lands in a distinct cell of a distinct variable, and the buffer
                   ! offset is a pure function of the three, so all three collapse.
                   !$acc loop vector collapse(3)
                   do ix3=1,nxCo_fc(3,2)
                     do ix1=1,nxCo_fc(1,2)
                       do iw=1,nwfluxin
                         psb(igrid)%w(ixmin1+ix1-1,ix,ixmin3+ix3-1,nw0+iw-1) = &
                           psb(igrid)%w(ixmin1+ix1-1,ix,ixmin3+ix3-1,nw0+iw-1) + &
                           recvbuffer(ibuf_offset(4**3*(igrid-1)+inc1+4*inc2+16*inc3+1) &
                              +(ix1-1)+(ix3-1)*nxCo_fc(1,2)+(iw-1)*nxCo_fc(1,2)*nxCo_fc(3,2)) &
                           / bgeo%dvolume(ixmin1+ix1-1,ix,ixmin3+ix3-1,igrid)
                       end do
                     end do
                   end do
#:endif
               end if
             end do
             end do
             end do
           end do

           case (3)
           do iside=1,2
             i1=kr(1,3)*(2*iside-3)
             i2=kr(2,3)*(2*iside-3)
             i3=kr(3,3)*(2*iside-3)

             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle

             if (neighbor_type(i1,i2,i3,igrid)/=4) cycle

 !opedit: skip over active/passive interface since flux for passive ones is
             ! not computed, keep the buffer counter up to date:
            !   if (.not.neighbor_active(i1,i2,i3,&
            !      igrid).or..not.neighbor_active(0,0,0,igrid) ) then
            !     do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
            !     inc3=2*i3+ic3
            ! do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
            !     inc2=2*i2+ic2
            ! do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
            !     inc1=2*i1+ic1
            !     ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
            !     if (ipe_neighbor/=mype) then
            !       ibufnext=ibuf+isize(3)
            !       ibuf=ibufnext
            !     end if
            !     end do
            ! end do
            ! end do
            !     cycle
            !   end if

             select case (iside)
             case (1)
               ix=ixMlo3
             case (2)
               ix=ixMhi3
             end select

             ! remove coarse flux
#:if GEOM == 'Cartesian'
               !$acc loop vector collapse(3)
               do ix2=ixMlo2,ixMhi2
                 do ix1=ixMlo1,ixMhi1 
                   do iw=1,nwfluxin
                     psb(igrid)%w(ix1,ix2,ix,nw0+iw-1) = &
                       psb(igrid)%w(ix1,ix2,ix,nw0+iw-1) - &
                       pflux(iside,3)%flux(ix1-nghostcells,ix2-nghostcells,&
                          1,iw,igrid)
                   end do
                 end do
               end do
#:else
               !$acc loop vector collapse(3)
               do ix2=ixMlo2,ixMhi2
                 do ix1=ixMlo1,ixMhi1
                   do iw=1,nwfluxin
                     psb(igrid)%w(ix1,ix2,ix,nw0+iw-1) = &
                       psb(igrid)%w(ix1,ix2,ix,nw0+iw-1) - &
                       pflux(iside,3)%flux(ix1-nghostcells,ix2-nghostcells,&
                          1,iw,igrid) &
                       / bgeo%dvolume(ix1,ix2,ix,igrid)
                   end do
                 end do
               end do
#:endif


             ! add fine flux
             do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
             do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
             do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
               ixmin1=ixMlo1+(ic1-1)*nxCo1
               ixmin2=ixMlo2+(ic2-1)*nxCo2
               ixmin3=ix
               ixmax1=ixmin1-1+nxCo1
               ixmax2=ixmin2-1+nxCo2
               ixmax3=ix
               if (ipe_neighbor==mype) then
                 iotherside=3-iside
#:if GEOM == 'Cartesian'
                   !$acc loop vector collapse(3)
                   do ix2=1,nxCo2 
                     do ix1=1,nxCo1 
                       do iw=1,nwfluxin
                         psb(igrid)%w(ixmin1+ix1-1,ixmin2+ix2-1,ix,nw0+iw-1) = &
                           psb(igrid)%w(ixmin1+ix1-1,ixmin2+ix2-1,ix,nw0+iw-1) + &
                           pflux(iotherside,3)%flux(ix1,ix2,1,iw,&
                              ineighbor)* CoFiratio
                       end do
                     end do
                   end do
#:else
                   !$acc loop vector collapse(3)
                   do ix2=1,nxCo2
                     do ix1=1,nxCo1
                       do iw=1,nwfluxin
                         psb(igrid)%w(ixmin1+ix1-1,ixmin2+ix2-1,ix,nw0+iw-1) = &
                           psb(igrid)%w(ixmin1+ix1-1,ixmin2+ix2-1,ix,nw0+iw-1) + &
                           pflux(iotherside,3)%flux(ix1,ix2,1,iw,&
                              ineighbor) &
                           / bgeo%dvolume(ixmin1+ix1-1,ixmin2+ix2-1,ix,igrid)
                       end do
                     end do
                   end do
#:endif
               !else
               !  if (slab_uniform) then
               !    ibufnext=ibuf+isize(3)
               !    if(stagger_grid) ibufnext=ibufnext-isize_stg(3)
               !    psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
               !       nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !       ixmin3:ixmax3,nw0:nw1)+CoFiratio &
               !       *reshape(source=recvbuffer(ibuf:ibufnext-1),&
               !        shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !       ixmin3:ixmax3,nw0:nw1)))
               !    ibuf=ibuf+isize(3)
               !  else
               !    ibufnext=ibuf+isize(3)
               !    if(stagger_grid) then
               !      nbuf=(isize(3)-isize_stg(3))/nwfluxin
               !    else
               !      nbuf=isize(3)/nwfluxin
               !    end if
               !    do iw=nw0,nw1
               !      psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
               !         iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !         ixmin3:ixmax3,iw) &
               !         +reshape(source=recvbuffer(ibuf:ibufnext-1),&
               !          shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               !         ixmin3:ixmax3,iw))) /ps(igrid)%dvolume(ixmin1:ixmax1,&
               !         ixmin2:ixmax2,ixmin3:ixmax3)
               !      ibuf=ibuf+nbuf
               !    end do
               !    ibuf=ibufnext
               !  end if
               else
#:if GEOM == 'Cartesian'
                   ! Two transverse indices plus the variable index: every iteration
                   ! lands in a distinct cell of a distinct variable, and the buffer
                   ! offset is a pure function of the three, so all three collapse.
                   !$acc loop vector collapse(3)
                   do ix2=1,nxCo_fc(2,3)
                     do ix1=1,nxCo_fc(1,3)
                       do iw=1,nwfluxin
                         psb(igrid)%w(ixmin1+ix1-1,ixmin2+ix2-1,ix,nw0+iw-1) = &
                           psb(igrid)%w(ixmin1+ix1-1,ixmin2+ix2-1,ix,nw0+iw-1) + &
                           recvbuffer(ibuf_offset(4**3*(igrid-1)+inc1+4*inc2+16*inc3+1) &
                              +(ix1-1)+(ix2-1)*nxCo_fc(1,3)+(iw-1)*nxCo_fc(1,3)*nxCo_fc(2,3)) * CoFiratio
                       end do
                     end do
                   end do
#:else
                   ! Two transverse indices plus the variable index: every iteration
                   ! lands in a distinct cell of a distinct variable, and the buffer
                   ! offset is a pure function of the three, so all three collapse.
                   !$acc loop vector collapse(3)
                   do ix2=1,nxCo_fc(2,3)
                     do ix1=1,nxCo_fc(1,3)
                       do iw=1,nwfluxin
                         psb(igrid)%w(ixmin1+ix1-1,ixmin2+ix2-1,ix,nw0+iw-1) = &
                           psb(igrid)%w(ixmin1+ix1-1,ixmin2+ix2-1,ix,nw0+iw-1) + &
                           recvbuffer(ibuf_offset(4**3*(igrid-1)+inc1+4*inc2+16*inc3+1) &
                              +(ix1-1)+(ix2-1)*nxCo_fc(1,3)+(iw-1)*nxCo_fc(1,3)*nxCo_fc(2,3)) &
                           / bgeo%dvolume(ixmin1+ix1-1,ixmin2+ix2-1,ix,igrid)
                       end do
                     end do
                   end do
#:endif
               end if
             end do
             end do
             end do
           end do
         end select
       end do
     end do

     if (n_send_pe>0) then
       call MPI_WAITALL(n_send_pe,fc_sendreq,fc_sendstat,ierrmpi)
     end if

   end subroutine fix_conserve

   subroutine store_flux(igrid,fC,idimmin,idimmax,nwfluxin)
     use mod_global_parameters

     integer, intent(in)          :: igrid, idimmin,idimmax, nwfluxin
     double precision, intent(in) :: fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
        ixGlo3:ixGhi3,1:nwfluxin,1:ndim)

     integer :: idims, iside, i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3, ix1,ix2,&
        ix3, ixCo1,ixCo2,ixCo3, nxCo1,nxCo2,nxCo3, iw
!!
!!     do idims = idimmin,idimmax
!!       select case (idims)
!!         case (1)
!!         do iside=1,2
!!           i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3)
!!           i3=kr(3,1)*(2*iside-3);
!!
!!           if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
!!
!!           select case (neighbor_type(i1,i2,i3,igrid))
!!           case (neighbor_fine)
!!             select case (iside)
!!             case (1)
!!               pflux(iside,1,igrid)%flux(1,:,:,1:nwfluxin) = -fC(nghostcells,&
!!                  ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:nwfluxin,1)
!!             case (2)
!!               pflux(iside,1,igrid)%flux(1,:,:,1:nwfluxin) = fC(ixMhi1,&
!!                  ixMlo2:ixMhi2,ixMlo3:ixMhi3,1:nwfluxin,1)
!!             end select
!!           case (neighbor_coarse)
!!             nxCo1=1;nxCo2=ixGhi2/2-nghostcells;nxCo3=ixGhi3/2-nghostcells;
!!             select case (iside)
!!             case (1)
!!               do iw=1,nwfluxin
!!                do ixCo3=1,nxCo3
!!         do ixCo2=1,nxCo2
!!         do ixCo1=1,nxCo1
!!                   ix1=nghostcells;ix2=ixMlo2+2*(ixCo2-1)
!!                   ix3=ixMlo3+2*(ixCo3-1);
!!                   pflux(iside,1,igrid)%flux(ixCo1,ixCo2,ixCo3,&
!!                      iw) = sum(fC(ix1,ix2:ix2+1,ix3:ix3+1,iw,1))
!!                end do
!!         end do
!!         end do
!!               end do
!!             case (2)
!!               do iw=1,nwfluxin
!!                do ixCo3=1,nxCo3
!!         do ixCo2=1,nxCo2
!!         do ixCo1=1,nxCo1
!!                   ix1=ixMhi1;ix2=ixMlo2+2*(ixCo2-1);ix3=ixMlo3+2*(ixCo3-1);
!!                   pflux(iside,1,igrid)%flux(ixCo1,ixCo2,ixCo3,&
!!                      iw) =-sum(fC(ix1,ix2:ix2+1,ix3:ix3+1,iw,1))
!!                end do
!!         end do
!!         end do
!!               end do
!!             end select
!!           end select
!!         end do
!!         case (2)
!!         do iside=1,2
!!           i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3)
!!           i3=kr(3,2)*(2*iside-3);
!!
!!           if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
!!
!!           select case (neighbor_type(i1,i2,i3,igrid))
!!           case (neighbor_fine)
!!             select case (iside)
!!             case (1)
!!               pflux(iside,2,igrid)%flux(:,1,:,1:nwfluxin) = -fC(ixMlo1:ixMhi1,&
!!                  nghostcells,ixMlo3:ixMhi3,1:nwfluxin,2)
!!             case (2)
!!               pflux(iside,2,igrid)%flux(:,1,:,1:nwfluxin) = fC(ixMlo1:ixMhi1,&
!!                  ixMhi2,ixMlo3:ixMhi3,1:nwfluxin,2)
!!             end select
!!           case (neighbor_coarse)
!!             nxCo1=ixGhi1/2-nghostcells;nxCo2=1;nxCo3=ixGhi3/2-nghostcells;
!!             select case (iside)
!!             case (1)
!!               do iw=1,nwfluxin
!!                do ixCo3=1,nxCo3
!!         do ixCo2=1,nxCo2
!!         do ixCo1=1,nxCo1
!!                   ix1=ixMlo1+2*(ixCo1-1);ix2=nghostcells
!!                   ix3=ixMlo3+2*(ixCo3-1);
!!                   pflux(iside,2,igrid)%flux(ixCo1,ixCo2,ixCo3,&
!!                      iw) = sum(fC(ix1:ix1+1,ix2,ix3:ix3+1,iw,2))
!!                end do
!!         end do
!!         end do
!!               end do
!!             case (2)
!!               do iw=1,nwfluxin
!!                do ixCo3=1,nxCo3
!!         do ixCo2=1,nxCo2
!!         do ixCo1=1,nxCo1
!!                   ix1=ixMlo1+2*(ixCo1-1);ix2=ixMhi2;ix3=ixMlo3+2*(ixCo3-1);
!!                   pflux(iside,2,igrid)%flux(ixCo1,ixCo2,ixCo3,&
!!                      iw) =-sum(fC(ix1:ix1+1,ix2,ix3:ix3+1,iw,2))
!!                end do
!!         end do
!!         end do
!!               end do
!!             end select
!!           end select
!!         end do
!!         case (3)
!!         do iside=1,2
!!           i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3)
!!           i3=kr(3,3)*(2*iside-3);
!!
!!           if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
!!
!!           select case (neighbor_type(i1,i2,i3,igrid))
!!           case (neighbor_fine)
!!             select case (iside)
!!             case (1)
!!               pflux(iside,3,igrid)%flux(:,:,1,1:nwfluxin) = -fC(ixMlo1:ixMhi1,&
!!                  ixMlo2:ixMhi2,nghostcells,1:nwfluxin,3)
!!             case (2)
!!               pflux(iside,3,igrid)%flux(:,:,1,1:nwfluxin) = fC(ixMlo1:ixMhi1,&
!!                  ixMlo2:ixMhi2,ixMhi3,1:nwfluxin,3)
!!             end select
!!           case (neighbor_coarse)
!!             nxCo1=ixGhi1/2-nghostcells;nxCo2=ixGhi2/2-nghostcells;nxCo3=1;
!!             select case (iside)
!!             case (1)
!!               do iw=1,nwfluxin
!!                do ixCo3=1,nxCo3
!!         do ixCo2=1,nxCo2
!!         do ixCo1=1,nxCo1
!!                   ix1=ixMlo1+2*(ixCo1-1);ix2=ixMlo2+2*(ixCo2-1)
!!                   ix3=nghostcells;
!!                   pflux(iside,3,igrid)%flux(ixCo1,ixCo2,ixCo3,&
!!                      iw) = sum(fC(ix1:ix1+1,ix2:ix2+1,ix3,iw,3))
!!                end do
!!         end do
!!         end do
!!               end do
!!             case (2)
!!               do iw=1,nwfluxin
!!                do ixCo3=1,nxCo3
!!         do ixCo2=1,nxCo2
!!         do ixCo1=1,nxCo1
!!                   ix1=ixMlo1+2*(ixCo1-1);ix2=ixMlo2+2*(ixCo2-1);ix3=ixMhi3;
!!                   pflux(iside,3,igrid)%flux(ixCo1,ixCo2,ixCo3,&
!!                      iw) =-sum(fC(ix1:ix1+1,ix2:ix2+1,ix3,iw,3))
!!                end do
!!         end do
!!         end do
!!               end do
!!             end select
!!           end select
!!         end do
!!       end select
!!     end do

   end subroutine store_flux

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ALL OF THE FOLLOWING IS FOR STAGGERED GRIDS !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine store_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      fE,idimmin,idimmax)
     use mod_global_parameters

     integer, intent(in)          :: igrid, ixImin1,ixImin2,ixImin3,ixImax1,&
        ixImax2,ixImax3, idimmin,idimmax
     double precision, intent(in) :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
        ixImin3:ixImax3,sdim:3)

     integer :: idims, idir, iside, i1,i2,i3
     integer :: pi1,pi2,pi3, mi1,mi2,mi3, ph1,ph2,ph3, mh1,mh2,mh3 !To detect corners
     integer :: ixMcmin1,ixMcmin2,ixMcmin3,ixMcmax1,ixMcmax2,ixMcmax3

     !!do idims = idimmin,idimmax  !loop over face directions
     !!  !! Loop over block faces
     !!  do iside=1,2
     !!    i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
     !!    i3=kr(3,idims)*(2*iside-3);
     !!    if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
     !!    select case (neighbor_type(i1,i2,i3,igrid))
     !!    case (neighbor_fine)
     !!      ! The neighbour is finer
     !!      ! Face direction, side (left or right), restrict ==ired?, fE
     !!      call flux_to_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     !!         ixImax3,idims,iside,.false.,fE)
     !!    case(neighbor_coarse)
     !!      ! The neighbour is coarser
     !!      call flux_to_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     !!         ixImax3,idims,iside,.true.,fE)
     !!    case(neighbor_sibling)
     !!      ! If the neighbour is at the same level,
     !!      ! check if there are corners
     !!      ! If there is any corner, store the fluxes from that side
     !!      do idir=idims+1,ndim
     !!        pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
     !!        mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
     !!        ph1=pi1-kr(idims,1)*(2*iside-3);ph2=pi2-kr(idims,2)*(2*iside-3)
     !!        ph3=pi3-kr(idims,3)*(2*iside-3);
     !!        mh1=mi1-kr(idims,1)*(2*iside-3);mh2=mi2-kr(idims,2)*(2*iside-3)
     !!        mh3=mi3-kr(idims,3)*(2*iside-3);
     !!        if (neighbor_type(pi1,pi2,pi3,igrid)==4.and.neighbor_type(ph1,ph2,&
     !!           ph3,igrid)==3) then
     !!          call flux_to_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     !!             ixImax3,idims,iside,.false.,fE)
     !!        end if
     !!        if (neighbor_type(mi1,mi2,mi3,igrid)==4.and.neighbor_type(mh1,mh2,&
     !!           mh3,igrid)==3) then
     !!          call flux_to_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     !!             ixImax3,idims,iside,.false.,fE)
     !!        end if
     !!      end do
     !!    end select
     !!  end do
     !!end do

   end subroutine store_edge

   subroutine flux_to_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
      ixImax3,idims,iside,restrict,fE)
     use mod_global_parameters

     integer                      :: igrid,ixImin1,ixImin2,ixImin3,ixImax1,&
        ixImax2,ixImax3,idims,iside
     logical                      :: restrict
     double precision, intent(in) :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
        ixImin3:ixImax3,sdim:3)

     integer                      :: idir1,idir2
     integer                      :: ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,&
        ixEmax3,ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3, jxFmin1,&
        jxFmin2,jxFmin3,jxFmax1,jxFmax2,jxFmax3, nx1,nx2,nx3,nxCo1,nxCo2,nxCo3

     nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
     nxCo1=nx1/2;nxCo2=nx2/2;nxCo3=nx3/2;
     ! ixE are the indices on the 'edge' array.
     ! ixF are the indices on the 'fE' array
     ! jxF are indices advanced to perform the flux restriction (sum) in 3D
     ! A line integral of the electric field on the coarse side
     ! lies over two edges on the fine side. So, in 3D we restrict by summing
     ! over two cells on the fine side.

     !!do idir1=1,ndim-1
     !!  ! 3D: rotate indices among 1 and 2 to save space
     !!  idir2=mod(idir1+idims-1,3)+1


     !!  if (restrict) then
     !!    ! Set up indices for restriction
     !!    ixFmin1=ixMlo1-1+kr(1,idir2);ixFmin2=ixMlo2-1+kr(2,idir2)
     !!    ixFmin3=ixMlo3-1+kr(3,idir2);
     !!    ixFmax1=ixMhi1-kr(1,idir2);ixFmax2=ixMhi2-kr(2,idir2)
     !!    ixFmax3=ixMhi3-kr(3,idir2);

     !!    jxFmin1=ixFmin1+kr(1,idir2);jxFmin2=ixFmin2+kr(2,idir2)
     !!    jxFmin3=ixFmin3+kr(3,idir2);jxFmax1=ixFmax1+kr(1,idir2)
     !!    jxFmax2=ixFmax2+kr(2,idir2);jxFmax3=ixFmax3+kr(3,idir2);

     !!    ixEmin1=0+kr(1,idir2);ixEmin2=0+kr(2,idir2);ixEmin3=0+kr(3,idir2);
     !!    ixEmax1=nxCo1;ixEmax2=nxCo2;ixEmax3=nxCo3;
     !!    select case(idims)
     !!   case(1)
     !!      ixEmin1=1;ixEmax1=1;
     !!      select case(iside)
     !!      case(1)
     !!        ixFmax1=ixFmin1
     !!        jxFmax1=ixFmin1
     !!      case(2)
     !!        ixFmin1=ixFmax1
     !!        jxFmin1=ixFmax1
     !!      end select

     !!   case(2)
     !!      ixEmin2=1;ixEmax2=1;
     !!      select case(iside)
     !!      case(1)
     !!        ixFmax2=ixFmin2
     !!        jxFmax2=ixFmin2
     !!      case(2)
     !!        ixFmin2=ixFmax2
     !!        jxFmin2=ixFmax2
     !!      end select

     !!   case(3)
     !!      ixEmin3=1;ixEmax3=1;
     !!      select case(iside)
     !!      case(1)
     !!        ixFmax3=ixFmin3
     !!        jxFmax3=ixFmin3
     !!      case(2)
     !!        ixFmin3=ixFmax3
     !!        jxFmin3=ixFmax3
     !!      end select

     !!    end select

     !!  pflux(iside,idims,igrid)%edge(ixEmin1:ixEmax1,ixEmin2:ixEmax2,&
     !!     ixEmin3:ixEmax3,idir1)=fE(ixFmin1:ixFmax1:2,ixFmin2:ixFmax2:2,&
     !!     ixFmin3:ixFmax3:2,idir2) +fE(jxFmin1:jxFmax1:2,jxFmin2:jxFmax2:2,&
     !!     jxFmin3:jxFmax3:2,idir2);

     !!  else
     !!    ! Set up indices for copying
     !!    ixFmin1=ixMlo1-1+kr(1,idir2);ixFmin2=ixMlo2-1+kr(2,idir2)
     !!    ixFmin3=ixMlo3-1+kr(3,idir2);
     !!    ixFmax1=ixMhi1;ixFmax2=ixMhi2;ixFmax3=ixMhi3;
     !!    ixEmin1=0+kr(1,idir2);ixEmin2=0+kr(2,idir2);ixEmin3=0+kr(3,idir2);
     !!    ixEmax1=nx1;ixEmax2=nx2;ixEmax3=nx3;

     !!    select case(idims)
     !!   case(1)
     !!      ixEmin1=1;ixEmax1=1;
     !!      select case(iside)
     !!      case(1)
     !!        ixFmax1=ixFmin1
     !!      case(2)
     !!        ixFmin1=ixFmax1
     !!      end select

     !!   case(2)
     !!      ixEmin2=1;ixEmax2=1;
     !!      select case(iside)
     !!      case(1)
     !!        ixFmax2=ixFmin2
     !!      case(2)
     !!        ixFmin2=ixFmax2
     !!      end select

     !!   case(3)
     !!      ixEmin3=1;ixEmax3=1;
     !!      select case(iside)
     !!      case(1)
     !!        ixFmax3=ixFmin3
     !!      case(2)
     !!        ixFmin3=ixFmax3
     !!      end select

     !!    end select

     !!    pflux(iside,idims,igrid)%edge(ixEmin1:ixEmax1,ixEmin2:ixEmax2,&
     !!       ixEmin3:ixEmax3,idir1)=fE(ixFmin1:ixFmax1,ixFmin2:ixFmax2,&
     !!       ixFmin3:ixFmax3,idir2)

     !!  end if

     !!end do

   end subroutine flux_to_edge

   subroutine fix_edges(psuse,idimmin,idimmax)
     use mod_global_parameters

     type(state) :: psuse(max_blocks)
     integer, intent(in) :: idimmin,idimmax

     integer :: iigrid, igrid, idims, iside, iotherside, i1,i2,i3, ic1,ic2,ic3,&
         inc1,inc2,inc3, ixMcmin1,ixMcmin2,ixMcmin3,ixMcmax1,ixMcmax2,ixMcmax3
     integer :: nbuf, ibufnext
     integer :: ibufnext_cc
     integer :: pi1,pi2,pi3, mi1,mi2,mi3, ph1,ph2,ph3, mh1,mh2,mh3 !To detect corners
     integer :: ixEmin1(1:3),ixEmin2(1:3),ixEmin3(1:3),ixEmax1(1:3),&
        ixEmax2(1:3),ixEmax3(1:3), ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,&
        ixtEmax2,ixtEmax3, ixFmin1(1:ndim),ixFmin2(1:ndim),ixFmin3(1:ndim),&
        ixFmax1(1:ndim),ixFmax2(1:ndim),ixFmax3(1:ndim), ixfEmin1(1:3),&
        ixfEmin2(1:3),ixfEmin3(1:3),ixfEmax1(1:3),ixfEmax2(1:3),ixfEmax3(1:3)
     integer :: nx1,nx2,nx3, idir, ix, ipe_neighbor, ineighbor
     logical :: pcorner(1:ndim),mcorner(1:ndim)

     !!if (nrecv_ct>0) then
     !!   call MPI_WAITALL(nrecv_ct,cc_recvreq,cc_recvstat,ierrmpi)
     !!end if

     !!! Initialise buffer counter again
     !!ibuf=1
     !!ibuf_cc=1
     !!do iigrid=1,igridstail; igrid=igrids(iigrid);
     !!  do idims= idimmin,idimmax
     !!    do iside=1,2
     !!      i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
     !!      i3=kr(3,idims)*(2*iside-3);
     !!      if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
     !!      select case(neighbor_type(i1,i2,i3,igrid))
     !!      case(neighbor_fine)
     !!        ! The first neighbour is finer
     !!        if (.not.neighbor_active(i1,i2,i3,&
     !!           igrid).or..not.neighbor_active(0,0,0,igrid) ) then
     !!          do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
     !!             inc3=2*i3+ic3
     !!          do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
     !!             inc2=2*i2+ic2
     !!          do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
     !!             inc1=2*i1+ic1
     !!             ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
     !!             !! When the neighbour is in a different process
     !!             if (ipe_neighbor/=mype) then
     !!                ibufnext=ibuf+isize(idims)
     !!                ibuf=ibufnext
     !!                end if
     !!          end do
     !!          end do
     !!          end do
     !!           cycle
     !!        end if

     !!        ! Check if there are corners
     !!        pcorner=.false.
     !!        mcorner=.false.
     !!        do idir=1,ndim
     !!          pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
     !!          mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
     !!          ph1=pi1-kr(idims,1)*(2*iside-3)
     !!          ph2=pi2-kr(idims,2)*(2*iside-3)
     !!          ph3=pi3-kr(idims,3)*(2*iside-3);
     !!          mh1=mi1-kr(idims,1)*(2*iside-3)
     !!          mh2=mi2-kr(idims,2)*(2*iside-3)
     !!          mh3=mi3-kr(idims,3)*(2*iside-3);
     !!          if (neighbor_type(ph1,ph2,ph3,&
     !!             igrid)==neighbor_fine) pcorner(idir)=.true.
     !!          if (neighbor_type(mh1,mh2,mh3,&
     !!             igrid)==neighbor_fine) mcorner(idir)=.true.
     !!        end do
     !!        ! Calculate indices range
     !!        call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
     !!           ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
     !!           ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
     !!           ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,igrid,idims,iside,.false.,&
     !!           .false.,0,0,0,pcorner,mcorner)
     !!        ! Remove coarse part of circulation
     !!        call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
     !!           ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,&
     !!           ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,&
     !!           ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,pflux(iside,idims,&
     !!           igrid)%edge,idims,iside,.false.,psuse(igrid))
     !!        ! Add fine part of the circulation
     !!       do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
     !!          inc3=2*i3+ic3
     !!       do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
     !!          inc2=2*i2+ic2
     !!       do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
     !!          inc1=2*i1+ic1
     !!          ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
     !!          ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
     !!          iotherside=3-iside
     !!          nx1=(ixMhi1-ixMlo1+1)/2;nx2=(ixMhi2-ixMlo2+1)/2
     !!          nx3=(ixMhi3-ixMlo3+1)/2;
     !!          call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!             ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!             ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
     !!             ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,igrid,&
     !!             idims,iside,.true.,.false.,inc1,inc2,inc3,pcorner,mcorner)
     !!          if (ipe_neighbor==mype) then
     !!            call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!               ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!               ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
     !!               ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
     !!               pflux(iotherside,idims,ineighbor)%edge,idims,iside,.true.,&
     !!               psuse(igrid))
     !!          else
     !!            ibufnext=ibuf+isize(idims)
     !!            call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!               ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!               ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
     !!               ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
     !!               reshape(source=recvbuffer(ibufnext-&
     !!               isize_stg(idims):ibufnext-1),shape=(/ ixtEmax1-ixtEmin1+1,&
     !!               ixtEmax2-ixtEmin2+1,ixtEmax3-ixtEmin3+1 ,3-1 /)),idims,&
     !!               iside,.true.,psuse(igrid))
     !!            ibuf=ibufnext
     !!          end if
     !!       end do
     !!       end do
     !!       end do

     !!      case(neighbor_sibling)
     !!        ! The first neighbour is at the same level
     !!        ! Check if there are corners
     !!        do idir=idims+1,ndim
     !!          pcorner=.false.
     !!          mcorner=.false.
     !!          pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
     !!          mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
     !!          ph1=pi1-kr(idims,1)*(2*iside-3)
     !!          ph2=pi2-kr(idims,2)*(2*iside-3)
     !!          ph3=pi3-kr(idims,3)*(2*iside-3);
     !!          mh1=mi1-kr(idims,1)*(2*iside-3)
     !!          mh2=mi2-kr(idims,2)*(2*iside-3)
     !!          mh3=mi3-kr(idims,3)*(2*iside-3);
     !!          if (neighbor_type(pi1,pi2,pi3,&
     !!             igrid)==neighbor_fine.and.neighbor_type(ph1,ph2,ph3,&
     !!             igrid)==neighbor_sibling.and.neighbor_pole(pi1,pi2,pi3,&
     !!             igrid)==0) then
     !!            pcorner(idir)=.true.
     !!            ! Remove coarse part
     !!            ! Set indices
     !!            call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!               ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!               ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
     !!               ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
     !!               igrid,idims,iside,.false.,.true.,0,0,0,pcorner,mcorner)
     !!            ! Remove
     !!            call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!               ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!               ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
     !!               ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
     !!               pflux(iside,idims,igrid)%edge,idims,iside,.false.,&
     !!               psuse(igrid))
     !!            ! Add fine part
     !!            ! Find relative position of finer grid
     !!  do ix=1,2
     !!            inc1=kr(idims,1)*3*(iside-1)+3*kr(idir,1)+kr(6-idir-idims,&
     !!               1)*ix
     !!            inc2=kr(idims,2)*3*(iside-1)+3*kr(idir,2)+kr(6-idir-idims,&
     !!               2)*ix
     !!            inc3=kr(idims,3)*3*(iside-1)+3*kr(idir,3)+kr(6-idir-idims,&
     !!               3)*ix;
     !!            ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
     !!            ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
     !!            iotherside=3-iside
     !!            ! Set indices
     !!            call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!               ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!               ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
     !!               ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
     !!               igrid,idims,iside,.true.,.true.,inc1,inc2,inc3,pcorner,&
     !!               mcorner)
     !!            ! add
     !!            if (ipe_neighbor==mype) then
     !!              call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!                 ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!                 ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
     !!                 ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
     !!                 pflux(iotherside,idims,ineighbor)%edge,idims,iside,&
     !!                 .true.,psuse(igrid))
     !!            else
     !!              ibufnext_cc=ibuf_cc+isize_stg(idims)
     !!              call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!                 ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!                 ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
     !!                 ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
     !!                 reshape(source=recvbuffer_cc(ibuf_cc:ibufnext_cc-1),&
     !!                 shape=(/ ixtEmax1-ixtEmin1+1,ixtEmax2-ixtEmin2+1,&
     !!                 ixtEmax3-ixtEmin3+1 ,3-1 /)),idims,iside,.true.,&
     !!                 psuse(igrid))
     !!              ibuf_cc=ibufnext_cc
     !!            end if
     !!  end do
     !!          ! Set CoCorner to false again for next step
     !!            pcorner(idir)=.false.
     !!          end if

     !!          if (neighbor_type(mi1,mi2,mi3,&
     !!             igrid)==neighbor_fine.and.neighbor_type(mh1,mh2,mh3,&
     !!             igrid)==neighbor_sibling.and.neighbor_pole(mi1,mi2,mi3,&
     !!             igrid)==0) then
     !!              mcorner(idir)=.true.
     !!              ! Remove coarse part
     !!              ! Set indices
     !!              call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!                 ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!                 ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
     !!                 ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
     !!                 igrid,idims,iside,.false.,.true.,0,0,0,pcorner,mcorner)
     !!              ! Remove
     !!              call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!                 ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!                 ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
     !!                 ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
     !!                 pflux(iside,idims,igrid)%edge,idims,iside,.false.,&
     !!                 psuse(igrid))
     !!              ! Add fine part
     !!              ! Find relative position of finer grid
     !!    do ix=1,2
     !!              inc1=kr(idims,1)*3*(iside-1)+kr(6-idir-idims,1)*ix
     !!              inc2=kr(idims,2)*3*(iside-1)+kr(6-idir-idims,2)*ix
     !!              inc3=kr(idims,3)*3*(iside-1)+kr(6-idir-idims,3)*ix;
     !!              ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
     !!              ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
     !!              iotherside=3-iside
     !!              ! Set indices
     !!              call set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!                 ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!                 ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,ixEmax3,&
     !!                 ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,ixfEmax3,&
     !!                 igrid,idims,iside,.true.,.true.,inc1,inc2,inc3,pcorner,&
     !!                 mcorner)
     !!              ! add
     !!              if (ipe_neighbor==mype) then
     !!                call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!                   ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!                   ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,&
     !!                   ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,&
     !!                   ixfEmax3,pflux(iotherside,idims,ineighbor)%edge,idims,&
     !!                   iside,.true.,psuse(igrid))
     !!              else
     !!                ibufnext_cc=ibuf_cc+isize_stg(idims)
     !!                call add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,&
     !!                   ixFmax3,ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,&
     !!                   ixtEmax3,ixEmin1,ixEmin2,ixEmin3,ixEmax1,ixEmax2,&
     !!                   ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,ixfEmax2,&
     !!                   ixfEmax3,reshape(source=recvbuffer_cc(&
     !!                   ibuf_cc:ibufnext_cc-1),shape=(/ ixtEmax1-ixtEmin1+1,&
     !!                   ixtEmax2-ixtEmin2+1,ixtEmax3-ixtEmin3+1 ,3-1 /)),idims,&
     !!                   iside,.true.,psuse(igrid))
     !!                ibuf_cc=ibufnext_cc
     !!              end if
     !!    end do
     !!            ! Set CoCorner to false again for next step
     !!             mcorner(idir)=.false.
     !!          end if
     !!        end do
     !!      end select
     !!    end do
     !!  end do
     !!end do

     !!if (nsend_ct>0) call MPI_WAITALL(nsend_ct,cc_sendreq,cc_sendstat,ierrmpi)

   end subroutine fix_edges

   !> This routine sets the indexes for the correction
   !> of the circulation according to several different
   !> cases, as grids located in different cpus,
   !> presence of corners, and different relative locations
   !> of the fine grid respect to the coarse one
   subroutine set_ix_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
      ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,ixEmin2,&
      ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,&
      ixfEmax2,ixfEmax3,igrid,idims,iside,add,CoCorner,inc1,inc2,inc3,pcorner,&
      mcorner)
     use mod_global_parameters

     integer,intent(in)    :: igrid,idims,iside,inc1,inc2,inc3
     logical,intent(in)    :: add,CoCorner
     logical,intent(inout) :: pcorner(1:ndim),mcorner(1:ndim)
     integer,intent(out)   :: ixFmin1(1:ndim),ixFmin2(1:ndim),ixFmin3(1:ndim),&
        ixFmax1(1:ndim),ixFmax2(1:ndim),ixFmax3(1:ndim),ixtEmin1,ixtEmin2,&
        ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1(1:3),ixEmin2(1:3),&
        ixEmin3(1:3),ixEmax1(1:3),ixEmax2(1:3),ixEmax3(1:3),ixfEmin1(1:3),&
        ixfEmin2(1:3),ixfEmin3(1:3),ixfEmax1(1:3),ixfEmax2(1:3),ixfEmax3(1:3) !Indices for faces and edges
     integer               :: icor1,icor2,icor3,idim1,idir,nx1,nx2,nx3,middle1,&
        middle2,middle3
     integer               :: ixtfEmin1,ixtfEmin2,ixtfEmin3,ixtfEmax1,&
        ixtfEmax2,ixtfEmax3

     ! ixF -> Indices for the _F_aces, and
     ! depends on the field component
     ! ixtE -> are the _t_otal range of the 'edge' array
     ! ixE -> are the ranges of the edge array,
     ! depending on the component
     ! ixfE -> are the ranges of the fE array (3D),
     ! and also depend on the component

     ! ... General ...
     ! Assign indices for the size of the E field array

     ixtfEmin1=ixMlo1-1;ixtfEmin2=ixMlo2-1;ixtfEmin3=ixMlo3-1;
     ixtfEmax1=ixMhi1;ixtfEmax2=ixMhi2;ixtfEmax3=ixMhi3;

     if(add) then
       nx1=(ixMhi1-ixMlo1+1)/2;nx2=(ixMhi2-ixMlo2+1)/2
       nx3=(ixMhi3-ixMlo3+1)/2;
     else
       nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
     end if

     do idim1=1,ndim
       ixtEmin1=0;ixtEmin2=0;ixtEmin3=0;
       ixtEmax1=nx1;ixtEmax2=nx2;ixtEmax3=nx3;
       select case(idims)
       case(1)
         ixtEmin1=1;ixtEmax1=1;
         if (iside==1) ixtfEmax1=ixtfEmin1;
         if (iside==2) ixtfEmin1=ixtfEmax1;

       case(2)
         ixtEmin2=1;ixtEmax2=1;
         if (iside==1) ixtfEmax2=ixtfEmin2;
         if (iside==2) ixtfEmin2=ixtfEmax2;

       case(3)
         ixtEmin3=1;ixtEmax3=1;
         if (iside==1) ixtfEmax3=ixtfEmin3;
         if (iside==2) ixtfEmin3=ixtfEmax3;

       end select
     end do

     ! Assign indices, considering only the face
     ! (idims and iside)
     do idim1=1,ndim
       ixFmin1(idim1)=ixMlo1-kr(idim1,1);ixFmin2(idim1)=ixMlo2-kr(idim1,2)
       ixFmin3(idim1)=ixMlo3-kr(idim1,3);
       ixFmax1(idim1)=ixMhi1;ixFmax2(idim1)=ixMhi2;ixFmax3(idim1)=ixMhi3;
       select case(idims)
       case(1)
          select case(iside)
          case(1)
          ixFmax1(idim1)=ixFmin1(idim1)
          case(2)
          ixFmin1(idim1)=ixFmax1(idim1)
          end select

       case(2)
          select case(iside)
          case(1)
          ixFmax2(idim1)=ixFmin2(idim1)
          case(2)
          ixFmin2(idim1)=ixFmax2(idim1)
          end select

       case(3)
          select case(iside)
          case(1)
          ixFmax3(idim1)=ixFmin3(idim1)
          case(2)
          ixFmin3(idim1)=ixFmax3(idim1)
          end select

       end select
     end do
     ! ... Relative position ...
     ! Restrict range using relative position
     if(add) then
       middle1=(ixMhi1+ixMlo1)/2;middle2=(ixMhi2+ixMlo2)/2
       middle3=(ixMhi3+ixMlo3)/2;

       if(inc1==1) then
         ixFmax1(:)=middle1
         ixtfEmax1=middle1
       end if
       if(inc1==2) then
         ixFmin1(:)=middle1+1
         ixtfEmin1=middle1
       end if


       if(inc2==1) then
         ixFmax2(:)=middle2
         ixtfEmax2=middle2
       end if
       if(inc2==2) then
         ixFmin2(:)=middle2+1
         ixtfEmin2=middle2
       end if


       if(inc3==1) then
         ixFmax3(:)=middle3
         ixtfEmax3=middle3
       end if
       if(inc3==2) then
         ixFmin3(:)=middle3+1
         ixtfEmin3=middle3
       end if

     end if
     ! ... Adjust ranges of edges according to direction ...
     do idim1=1,3
       ixfEmax1(idim1)=ixtfEmax1;ixfEmax2(idim1)=ixtfEmax2
       ixfEmax3(idim1)=ixtfEmax3;
       ixEmax1(idim1)=ixtEmax1;ixEmax2(idim1)=ixtEmax2
       ixEmax3(idim1)=ixtEmax3;
       ixfEmin1(idim1)=ixtfEmin1+kr(idim1,1)
       ixfEmin2(idim1)=ixtfEmin2+kr(idim1,2)
       ixfEmin3(idim1)=ixtfEmin3+kr(idim1,3);
       ixEmin1(idim1)=ixtEmin1+kr(idim1,1)
       ixEmin2(idim1)=ixtEmin2+kr(idim1,2)
       ixEmin3(idim1)=ixtEmin3+kr(idim1,3);
     end do
     ! ... Corners ...
     ! 'Coarse' corners
     if (CoCorner) then
       do idim1=idims+1,ndim
         if (pcorner(idim1)) then
           do idir=1,3!Index arrays have size ndim
             if (idir==6-idim1-idims) then
              !!! Something here has to change
              !!! Array ixfE must have size 3, while
              !!! ixE must have size ndim
              if (1==idim1) then
                 ixfEmin1(idir)=ixfEmax1(idir)
                 if (add) then
                   ixEmax1(idir) =ixEmin1(idir)
                 else
                   ixEmin1(idir) =ixEmax1(idir)
                 end if
               end if
              if (2==idim1) then
                 ixfEmin2(idir)=ixfEmax2(idir)
                 if (add) then
                   ixEmax2(idir) =ixEmin2(idir)
                 else
                   ixEmin2(idir) =ixEmax2(idir)
                 end if
               end if
              if (3==idim1) then
                 ixfEmin3(idir)=ixfEmax3(idir)
                 if (add) then
                   ixEmax3(idir) =ixEmin3(idir)
                 else
                   ixEmin3(idir) =ixEmax3(idir)
                 end if
               end if
             else
               ixEmin1(idir)=1;ixEmin2(idir)=1;ixEmin3(idir)=1;
               ixEmax1(idir)=0;ixEmax2(idir)=0;ixEmax3(idir)=0;
               ixfEmin1(idir)=1;ixfEmin2(idir)=1;ixfEmin3(idir)=1;
               ixfEmax1(idir)=0;ixfEmax2(idir)=0;ixfEmax3(idir)=0;
             end if
           end do
         end if
         if (mcorner(idim1)) then
           do idir=1,3
             if (idir==6-idim1-idims) then
              if (1==idim1) then
                 ixfEmax1(idir)=ixfEmin1(idir)
                 if (add) then
                   ixEmin1(idir) =ixEmax1(idir)
                 else
                   ixEmax1(idir) =ixEmin1(idir)
                 end if
               end if
              if (2==idim1) then
                 ixfEmax2(idir)=ixfEmin2(idir)
                 if (add) then
                   ixEmin2(idir) =ixEmax2(idir)
                 else
                   ixEmax2(idir) =ixEmin2(idir)
                 end if
               end if
              if (3==idim1) then
                 ixfEmax3(idir)=ixfEmin3(idir)
                 if (add) then
                   ixEmin3(idir) =ixEmax3(idir)
                 else
                   ixEmax3(idir) =ixEmin3(idir)
                 end if
               end if
             else
               ixEmin1(idir)=1;ixEmin2(idir)=1;ixEmin3(idir)=1;
               ixEmax1(idir)=0;ixEmax2(idir)=0;ixEmax3(idir)=0;
               ixfEmin1(idir)=1;ixfEmin2(idir)=1;ixfEmin3(idir)=1;
               ixfEmax1(idir)=0;ixfEmax2(idir)=0;ixfEmax3(idir)=0;
             end if
           end do
         end if
       end do
     else
     ! Other kinds of corners
     ! Crop ranges to account for corners
     ! When the fine fluxes are added, we consider
     ! whether they come from the same cpu or from
     ! a different one, in order to minimise the
     ! amount of communication
     ! Case for different processors still not implemented!!!
      if((idims.gt.1).and.pcorner(1)) then
         if((.not.add).or.(inc1==2)) then
 !ixFmax1(:)=ixFmax1(:)-kr(1,1);!ixFmax2(:)=ixFmax2(:)-kr(1,2);!ixFmax3(:)=ixFmax3(:)-kr(1,3);
           do idir=1,3
             if ((idir==idims).or.(idir==1)) cycle
               ixfEmax1(idir)=ixfEmax1(idir)-1
               ixEmax1(idir)=ixEmax1(idir)-1
           end do
         end if
       end if
      if((idims.gt.2).and.pcorner(2)) then
         if((.not.add).or.(inc2==2)) then
 !ixFmax1(:)=ixFmax1(:)-kr(2,1);!ixFmax2(:)=ixFmax2(:)-kr(2,2);!ixFmax3(:)=ixFmax3(:)-kr(2,3);
           do idir=1,3
             if ((idir==idims).or.(idir==2)) cycle
               ixfEmax2(idir)=ixfEmax2(idir)-1
               ixEmax2(idir)=ixEmax2(idir)-1
           end do
         end if
       end if
      if((idims.gt.3).and.pcorner(3)) then
         if((.not.add).or.(inc3==2)) then
 !ixFmax1(:)=ixFmax1(:)-kr(3,1);!ixFmax2(:)=ixFmax2(:)-kr(3,2);!ixFmax3(:)=ixFmax3(:)-kr(3,3);
           do idir=1,3
             if ((idir==idims).or.(idir==3)) cycle
               ixfEmax3(idir)=ixfEmax3(idir)-1
               ixEmax3(idir)=ixEmax3(idir)-1
           end do
         end if
       end if
      if((idims>1).and.mcorner(1)) then
         if((.not.add).or.(inc1==1)) then
 !ixFmin1(:)=ixFmin1(:)+kr(1,1);!ixFmin2(:)=ixFmin2(:)+kr(1,2);!ixFmin3(:)=ixFmin3(:)+kr(1,3);
           do idir=1,3
             if ((idir==idims).or.(idir==1)) cycle
               ixfEmin1(idir)=ixfEmin1(idir)+1
               ixEmin1(idir)=ixEmin1(idir)+1
           end do
         end if
       end if
      if((idims>2).and.mcorner(2)) then
         if((.not.add).or.(inc2==1)) then
 !ixFmin1(:)=ixFmin1(:)+kr(2,1);!ixFmin2(:)=ixFmin2(:)+kr(2,2);!ixFmin3(:)=ixFmin3(:)+kr(2,3);
           do idir=1,3
             if ((idir==idims).or.(idir==2)) cycle
               ixfEmin2(idir)=ixfEmin2(idir)+1
               ixEmin2(idir)=ixEmin2(idir)+1
           end do
         end if
       end if
      if((idims>3).and.mcorner(3)) then
         if((.not.add).or.(inc3==1)) then
 !ixFmin1(:)=ixFmin1(:)+kr(3,1);!ixFmin2(:)=ixFmin2(:)+kr(3,2);!ixFmin3(:)=ixFmin3(:)+kr(3,3);
           do idir=1,3
             if ((idir==idims).or.(idir==3)) cycle
               ixfEmin3(idir)=ixfEmin3(idir)+1
               ixEmin3(idir)=ixEmin3(idir)+1
           end do
         end if
       end if
     end if

   end subroutine set_ix_circ

   subroutine add_sub_circ(ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3,&
      ixtEmin1,ixtEmin2,ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1,ixEmin2,&
      ixEmin3,ixEmax1,ixEmax2,ixEmax3,ixfEmin1,ixfEmin2,ixfEmin3,ixfEmax1,&
      ixfEmax2,ixfEmax3,edge,idims,iside,add,s)
     use mod_global_parameters

     type(state)        :: s
     integer,intent(in) :: idims,iside
     integer            :: ixFmin1(1:ndim),ixFmin2(1:ndim),ixFmin3(1:ndim),&
        ixFmax1(1:ndim),ixFmax2(1:ndim),ixFmax3(1:ndim),ixtEmin1,ixtEmin2,&
        ixtEmin3,ixtEmax1,ixtEmax2,ixtEmax3,ixEmin1(1:3),ixEmin2(1:3),&
        ixEmin3(1:3),ixEmax1(1:3),ixEmax2(1:3),ixEmax3(1:3),ixfEmin1(1:3),&
        ixfEmin2(1:3),ixfEmin3(1:3),ixfEmax1(1:3),ixfEmax2(1:3),ixfEmax3(1:3)
     double precision   :: edge(ixtEmin1:ixtEmax1,ixtEmin2:ixtEmax2,&
        ixtEmin3:ixtEmax3,1:ndim-1)
     logical,intent(in) :: add

     integer            :: idim1,idim2,idir,middle1,middle2,middle3
     integer            :: ixfECmin1,ixfECmin2,ixfECmin3,ixfECmax1,ixfECmax2,&
        ixfECmax3,ixECmin1,ixECmin2,ixECmin3,ixECmax1,ixECmax2,ixECmax3
     double precision   :: fE(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
        sdim:3)
     double precision   :: circ(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
        1:ndim)
     integer            :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,hxmin1,&
        hxmin2,hxmin3,hxmax1,hxmax2,hxmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
        ixCmax2,ixCmax3,hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,hxCmax3 !Indices for edges

     ! ixF -> Indices for the faces, depends on the field component
     ! ixE -> Total range for the edges
     ! ixfE -> Edges in fE (3D) array
     ! ix,hx,ixC,hxC -> Auxiliary indices
     ! Assign quantities stored ad edges to make it as similar as
     ! possible to the routine updatefaces.
     fE(:,:,:,:)=zero
     do idim1=1,ndim-1
        ! 3D: rotate indices (see routine flux_to_edge)
       idir=mod(idim1+idims-1,3)+1

       ixfECmin1=ixfEmin1(idir);ixfECmin2=ixfEmin2(idir)
       ixfECmin3=ixfEmin3(idir);ixfECmax1=ixfEmax1(idir)
       ixfECmax2=ixfEmax2(idir);ixfECmax3=ixfEmax3(idir);
       ixECmin1=ixEmin1(idir);ixECmin2=ixEmin2(idir);ixECmin3=ixEmin3(idir)
       ixECmax1=ixEmax1(idir);ixECmax2=ixEmax2(idir);ixECmax3=ixEmax3(idir);
       fE(ixfECmin1:ixfECmax1,ixfECmin2:ixfECmax2,ixfECmin3:ixfECmax3,&
          idir)=edge(ixECmin1:ixECmax1,ixECmin2:ixECmax2,ixECmin3:ixECmax3,&
          idim1)
     end do

     ! Calculate part of circulation needed
     circ=zero
     do idim1=1,ndim
        do idim2=1,ndim
           do idir=sdim,3
             if (lvc(idim1,idim2,idir)==0) cycle
             ! Assemble indices
             ixCmin1=ixFmin1(idim1);ixCmin2=ixFmin2(idim1)
             ixCmin3=ixFmin3(idim1);ixCmax1=ixFmax1(idim1)
             ixCmax2=ixFmax2(idim1);ixCmax3=ixFmax3(idim1);
             hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
             hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
             hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
             if(idim1==idims) then
               circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)-fE(hxCmin1:hxCmax1,&
                  hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir))
             else
               select case(iside)
               case(2)
                 circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                    idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3,idim1)+lvc(idim1,idim2,&
                    idir)*fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                    idir)
               case(1)
                 circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                    idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    ixCmin3:ixCmax3,idim1)-lvc(idim1,idim2,&
                    idir)*fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
                    idir)
               end select
             end if
           end do
        end do
     end do

     ! Divide circulation by surface and add
     do idim1=1,ndim
        ixCmin1=ixFmin1(idim1);ixCmin2=ixFmin2(idim1);ixCmin3=ixFmin3(idim1)
        ixCmax1=ixFmax1(idim1);ixCmax2=ixFmax2(idim1);ixCmax3=ixFmax3(idim1);
        where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)>1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3))
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)
        elsewhere
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
        end where
        ! Add/subtract to field at face
        if (add) then
          s%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)=s%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)
        else
          s%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)=s%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)+circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)
        end if
     end do

   end subroutine add_sub_circ

end module mod_fix_conserve
