module mod_load_balance
#ifdef USE_MPIWRAPPERS
  use mod_mpi_wrapper
#else
#define mpi_irecv_wrapper MPI_IRECV
#define mpi_isend_wrapper MPI_ISEND
#endif
  implicit none
  private
  !> MPI recv send variables for AMR
  integer :: itag, irecv, isend
  integer :: recv_igrid, recv_ipe, send_igrid, send_ipe, igrid
  integer, dimension(:), allocatable :: recvrequest, sendrequest
  !> MPI recv send variables for staggered-variable AMR
  integer, dimension(:), allocatable :: recvrequest_stg, sendrequest_stg
  !> MPI buffers to send blocks
  double precision, allocatable, dimension(:,:,:,:,:)  :: snd_buff_lb, rcv_buff_lb
  integer, allocatable, dimension(:)  :: rcv_info_lb
  !$acc declare create(snd_buff_lb,rcv_buff_lb,rcv_info_lb)
  !$omp declare target(snd_buff_lb,rcv_buff_lb,rcv_info_lb)
  !> maximum number of blocks to send
  integer, parameter :: max_buff=1024

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
    !> MPI recv send variables for AMR
    integer, dimension(:,:), allocatable :: recvstatus, sendstatus
    !> MPI recv send variables for staggered-variable AMR
    integer :: itag_stg
    integer, dimension(:,:), allocatable :: recvstatus_stg, sendstatus_stg
    integer :: ix1, ix2, ix3, iw, ibuff

    ! Jannis: for now, not using version for passive/active blocks
    call get_Morton_range()

    if (npe==1) then
       sfc_to_igrid(:)=sfc(1,Morton_start(mype):Morton_stop(mype))
       return
    end if

    irecv=0
    isend=0
    allocate(recvstatus(MPI_STATUS_SIZE,max_blocks),recvrequest(max_blocks),&
        sendstatus(MPI_STATUS_SIZE,max_blocks),sendrequest(max_blocks))
    recvrequest=MPI_REQUEST_NULL
    sendrequest=MPI_REQUEST_NULL

    if(stagger_grid) then
      allocate(recvstatus_stg(MPI_STATUS_SIZE,max_blocks*3),&
         recvrequest_stg(max_blocks*3), sendstatus_stg(MPI_STATUS_SIZE,&
         max_blocks*3),sendrequest_stg(max_blocks*3))
      recvrequest_stg=MPI_REQUEST_NULL
      sendrequest_stg=MPI_REQUEST_NULL
    end if

    ! Allocate the send and receive buffers
    if ( .not. allocated(snd_buff_lb) ) then
       allocate( snd_buff_lb(block_nx1, block_nx2, block_nx3, nw, max_buff), &
            rcv_buff_lb(block_nx1, block_nx2, block_nx3, nw, max_buff), &
            rcv_info_lb(max_buff) )
       !$acc update device(snd_buff_lb, rcv_buff_lb, rcv_info_lb)
       !$omp target update to(snd_buff_lb, rcv_buff_lb, rcv_info_lb)
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

    if (irecv>0) then
      call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
      if(stagger_grid) call MPI_WAITALL(irecv,recvrequest_stg,recvstatus_stg,&
         ierrmpi)
    end if
    if (isend>0) then
      call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
      if(stagger_grid) call MPI_WAITALL(isend,sendrequest_stg,sendstatus_stg,&
         ierrmpi)
   end if

    ! unpack the receive buffers on GPU
#ifdef NOGPUDIRECT
   !$acc update device(rcv_buff_lb(:,:,:,:,1:irecv))
   !$omp target update to(rcv_buff_lb(:,:,:,:,1:irecv))
#endif
    !$acc update device(rcv_info_lb(1:irecv))
    !$acc parallel loop gang
    !$omp target update to(rcv_info_lb(1:irecv))
    !$omp target teams distribute
    do ibuff = 1, irecv
       recv_igrid = rcv_info_lb(ibuff)
       !$acc loop collapse(4) vector
       !$omp parallel do collapse(4)
       do iw = 1, nw
          do ix3 = 1, block_nx3
             do ix2 = 1, block_nx2
                do ix1 = 1, block_nx1
                   bg(1)%w(ixMlo1-1 + ix1, ixMlo2-1 + ix2, ixMlo3-1 + ix3, iw, recv_igrid) &
                        = rcv_buff_lb(ix1, ix2, ix3, iw, ibuff)
                end do
             end do
          end do
       end do
    end do

    deallocate(recvstatus,recvrequest,sendstatus,sendrequest)
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

  subroutine lb_recv
    use mod_global_parameters
    use mod_amr_solution_node, only: alloc_node
    use mod_comm_lib, only: mpistop

    call alloc_node(recv_igrid)

    itag=recv_igrid
    irecv=irecv+1
    if (irecv > max_buff) then
       call mpistop('load_balance: max_buff too small in receive')
    end if
#ifndef NOGPUDIRECT
    !$acc host_data use_device(rcv_buff_lb)
    !$omp target data use_device_addr(rcv_buff_lb)
#endif
    call mpi_irecv_wrapper(rcv_buff_lb(:,:,:,:,irecv), &
                   block_nx1*block_nx2*block_nx3*nw, MPI_DOUBLE_PRECISION, &
          send_ipe,itag, icomm, &
          recvrequest(irecv),ierrmpi)
#ifndef NOGPUDIRECT
    !$acc end host_data
    !$omp end target data
#endif
    rcv_info_lb(irecv) = recv_igrid
    if(stagger_grid) then
       itag=recv_igrid+max_blocks
       call mpi_irecv_wrapper(ps(recv_igrid)%ws,1,type_block_io_stg,send_ipe,itag,&
             icomm,recvrequest_stg(irecv),ierrmpi)
    end if
  end subroutine lb_recv

  subroutine lb_send
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer :: ix1, ix2, ix3, iw, ibuff

    itag=recv_igrid
    isend=isend+1
    if (isend > max_buff) then
       call mpistop('load_balance: max_buff too small in send')
    end if
    !$acc parallel loop gang default(present)
    !$omp target teams distribute
    do iw = 1, nw
       !$acc loop collapse(3) vector
       !$omp parallel do collapse(3)
       do ix3 = 1, block_nx3
          do ix2 = 1, block_nx2
             do ix1 = 1, block_nx1
                snd_buff_lb(ix1, ix2, ix3, iw, isend) = &
                      bg(1)%w(ixMlo1-1+ix1, ixMlo2-1+ix2, ixMlo3-1+ix3, iw, send_igrid)
             end do
          end do
       end do
    end do

#ifndef NOGPUDIRECT
    !$acc host_data use_device(snd_buff_lb)
    !$omp target data use_device_addr(snd_buff_lb)
#else
    !$acc update host(snd_buff_lb(:,:,:,:,isend))
    !$omp target update from(snd_buff_lb(:,:,:,:,isend))
#endif
    call mpi_isend_wrapper(snd_buff_lb(:,:,:,:,isend), &
                   block_nx1*block_nx2*block_nx3*nw, MPI_DOUBLE_PRECISION, &
          recv_ipe,itag, icomm, sendrequest(isend),ierrmpi)
#ifndef NOGPUDIRECT
    !$acc end host_data
    !$omp end target data
#endif
    if(stagger_grid) then
       itag=recv_igrid+max_blocks
       call mpi_isend_wrapper(ps(send_igrid)%ws,1,type_block_io_stg,recv_ipe,itag,&
             icomm,sendrequest_stg(isend),ierrmpi)
    end if
  end subroutine lb_send

end module mod_load_balance
