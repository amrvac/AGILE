!> Module to coarsen and refine grids for AMR
module mod_coarsen_refine
#ifdef USE_MPIWRAPPERS
  use mod_mpi_wrapper
#else
#define mpi_irecv_wrapper MPI_IRECV
#define mpi_isend_wrapper MPI_ISEND
#endif
  implicit none
  private
  !> MPI recv send variables for AMR.  One request per neighbouring rank now,
  !> not one per chunk, so these are sized by npe.
  integer :: itag, irecv, isend
  integer, dimension(:), allocatable :: recvrequest, sendrequest
  integer, dimension(:,:), allocatable :: recvstatus, sendstatus

  !> MPI buffers to send non-local coarsened grids.
  !>
  !> The exchange is aggregated per destination rank: every coarsened child
  !> bound for one peer occupies a contiguous run of snd_buff_cf and travels in
  !> a single Isend, so the message count is the number of neighbouring ranks
  !> rather than the number of (parent, child) pairs.  That is also what lets
  !> the tag be a constant; it used to be itag=ipeFi+igridFi, which is unique
  !> only because the matching Irecv names ipeFi as its source, and was itself a
  !> workaround for MPI_TAG_UB.
  !>
  !> Both ranks order a peer's run by the key
  !>
  !>     ikey = 8*(igridCo-1) + (ic1-1) + 2*(ic2-1) + 4*(ic3-1)
  !>
  !> built from the *receiving* parent's block index: the sender has it as the
  !> igrid argument of coarsen_grid_siblings, the receiver as its own igrid.
  !> Both know it because amr_coarsen_refine walks the globally replicated
  !> coarsen() table in the same order on every rank and calls getnode() for
  !> every event, so the free-node bookkeeping is replicated too.  Sorting by
  !> the key makes the two layouts agree element for element with no handshake.
  !>
  !> Aggregation gives up MPI's own mismatch detection, so the arrival lengths
  !> are checked against the layout in exchange_coarsened_blocks.
  double precision, allocatable, dimension(:) :: snd_buff_cf, rcv_buff_cf
  !$acc declare create(snd_buff_cf,rcv_buff_cf)

  !> Per outgoing / incoming chunk, appended by coarsen_grid_siblings during
  !> the walk and turned into per-peer runs by exchange_coarsened_blocks.
  !> cf_snd_igrid / cf_snd_ibuf and the rcv counterparts are read by the pack
  !> and unpack kernels and are therefore device-resident; cf_snd_dest,
  !> cf_snd_key, cf_rcv_src and cf_rcv_key are the host-side working set the
  !> layout is derived from.
  integer :: n_cf_snd, n_cf_rcv
  integer, allocatable, dimension(:) :: cf_snd_igrid, cf_snd_ibuf
  integer, allocatable, dimension(:) :: cf_snd_dest, cf_snd_key
  integer, allocatable, dimension(:) :: cf_rcv_igrid, cf_rcv_ibuf
  integer, allocatable, dimension(:) :: cf_rcv_src, cf_rcv_key
  integer, allocatable, dimension(:,:) :: cf_rcv_ic
  !$acc declare create(cf_snd_igrid,cf_snd_ibuf,cf_rcv_igrid,cf_rcv_ibuf,&
  !$acc&               cf_rcv_ic)

  !> The peers themselves, ascending, with the extent of each one's run.  Host
  !> only: these drive the MPI calls and nothing else.
  integer :: n_cf_send_pe, n_cf_recv_pe
  integer, allocatable, dimension(:) :: cf_send_pe, cf_send_pe_off,&
     cf_send_pe_len
  integer, allocatable, dimension(:) :: cf_recv_pe, cf_recv_pe_off,&
     cf_recv_pe_len

  !> maximum number of coarse blocks that can be sent after coarsening
  integer, parameter :: max_buff=1024
  !$acc declare copyin(max_buff)
  !> MPI recv send variables for staggered-variable AMR
  integer :: itag_stg
  integer, dimension(:), allocatable :: recvrequest_stg, sendrequest_stg
  integer, dimension(:,:), allocatable :: recvstatus_stg, sendstatus_stg

  ! Public subroutines
  public :: amr_coarsen_refine

contains

  !> coarsen and refine blocks to update AMR grid
  subroutine amr_coarsen_refine
    use mod_forest
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_usr_methods, only: usr_after_refine
    use mod_amr_fct
    use mod_space_filling_curve
    use mod_load_balance
    use mod_functions_connectivity, only: get_level_range,getigrids,&
         build_connectivity
    use mod_amr_solution_node, only: getnode, putnode
    use mod_functions_forest, only: coarsen_tree_leaf,refine_tree_leaf
    use mod_selectgrids, only: selectgrids
    use mod_refine, only: refine_grids
    use mod_multigrid_coupling


    integer :: iigrid, igrid, ipe, igridCo, ipeCo, level, ic1,ic2,ic3
    integer, dimension(2,2,2) :: igridFi, ipeFi
    integer :: n_coarsen, n_refine
    type(tree_node_ptr) :: tree, sibling
    logical             :: active
    integer             :: ibuff, iw, ix1, ix2, ix3

    call proper_nesting

    if(stagger_grid) then
       call store_faces
       call comm_faces
    end if

    n_coarsen = count(coarsen(:, :))
    n_refine = count(refine(:, :))

    ! to save memory: first coarsen then refine
    irecv=0
    isend=0
    n_cf_snd=0
    n_cf_rcv=0
    if (.not.allocated(recvrequest)) then
       allocate(recvstatus(MPI_STATUS_SIZE,npe),recvrequest(npe),&
            sendstatus(MPI_STATUS_SIZE,npe),sendrequest(npe))
       allocate(cf_send_pe(npe),cf_send_pe_off(npe),cf_send_pe_len(npe),&
            cf_recv_pe(npe),cf_recv_pe_off(npe),cf_recv_pe_len(npe))
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

    ! Allocate the chunk descriptors and the exchange buffers.  Both are
    ! capped by max_buff chunks, which coarsen_grid_siblings checks against.
    if ( .not. allocated(snd_buff_cf) ) then
       allocate( snd_buff_cf((block_nx1/2)*(block_nx2/2)*(block_nx3/2)*nw &
            *max_buff), &
            rcv_buff_cf((block_nx1/2)*(block_nx2/2)*(block_nx3/2)*nw &
            *max_buff) )
       allocate( cf_snd_igrid(max_buff), cf_snd_ibuf(max_buff), &
            cf_snd_dest(max_buff), cf_snd_key(max_buff), &
            cf_rcv_igrid(max_buff), cf_rcv_ibuf(max_buff), &
            cf_rcv_src(max_buff), cf_rcv_key(max_buff), cf_rcv_ic(3,max_buff) )
       !$acc update device(snd_buff_cf, rcv_buff_cf)
    end if

    do ipe=0,npe-1
       do igrid=1,max_blocks
          if (coarsen(igrid,ipe)) then
             if (.not.associated(igrid_to_node(igrid,ipe)%node)) cycle

             tree%node => igrid_to_node(igrid,ipe)%node%parent%node
             do ic3=1,2
                do ic2=1,2
                   do ic1=1,2
                      sibling%node => tree%node%child(ic1,ic2,ic3)%node
                      ipeFi(ic1,ic2,ic3)=sibling%node%ipe
                      igridFi(ic1,ic2,ic3)=sibling%node%igrid
                   end do
                end do
             end do

             ipeCo=ipeFi(1,1,1)
             igridCo=getnode(ipeCo)

             call coarsen_tree_leaf(igridCo,ipeCo,igridFi,ipeFi,active)

             call coarsen_grid_siblings(igridCo,ipeCo,igridFi,ipeFi,active)

             ! local coarsening done
             do ic3=1,2
                do ic2=1,2
                   do ic1=1,2
                      if (ipeFi(ic1,ic2,ic3)==ipeCo) then
                         call putnode(igridFi(ic1,ic2,ic3),ipeFi(ic1,ic2,ic3))
                         coarsen(igridFi(ic1,ic2,ic3),ipeFi(ic1,ic2,ic3))=.false.
                      end if
                   end do
                end do
             end do
          end if
       end do
    end do

    ! Lay the chunks out per peer, exchange them in one message each, and
    ! apply them.  The walk above only recorded descriptors.
    call exchange_coarsened_blocks

    ! The staggered path still posts per chunk; see the note in
    ! coarsen_grid_siblings.  It is unreachable here (stagger_grid is false and
    ! fix_edges is not called), so it is left as it was found.
    if(stagger_grid) then
       if (irecv>0) call MPI_WAITALL(irecv,recvrequest_stg,recvstatus_stg,&
            ierrmpi)
       if (isend>0) call MPI_WAITALL(isend,sendrequest_stg,sendstatus_stg,&
            ierrmpi)
       deallocate(recvstatus_stg,recvrequest_stg,sendstatus_stg,&
            sendrequest_stg)
    end if

    ! non-local coarsening done
    do ipe=0,npe-1
       do igrid=1,max_blocks
          if (coarsen(igrid,ipe)) then
             !if (ipe==mype) call dealloc_node(igrid) ! do not deallocate node
             ! memory preventing fragmentization of system memory as a result
             ! of frequent allocating and deallocating memory

             ! put the node (igrid number) into unused.
             call putnode(igrid,ipe)
             coarsen(igrid,ipe)=.false.
          end if
       end do
    end do

    do ipe=0,npe-1
       do igrid=1,max_blocks
          if (refine(igrid,ipe)) then

             do ic3=1,2
                do ic2=1,2
                   do ic1=1,2
                      igridFi(ic1,ic2,ic3)=getnode(ipe)
                      ipeFi(ic1,ic2,ic3)=ipe
                   end do
                end do
             end do

             call refine_tree_leaf(igridFi,ipeFi,igrid,ipe,active)

             if (ipe==mype) call refine_grids(igridFi,ipeFi,igrid,ipe,active)

             ! refinement done
             call putnode(igrid,ipe)
             refine(igrid,ipe)=.false.
          end if
       end do
    end do

    ! A crash occurs in later MPI_WAITALL when initial condition comsumes too
    ! much time to filling new blocks with both gfortran and intel fortran compiler.
    ! This barrier cure this problem
    !TODO to find the reason
    if(.not.time_advance) call MPI_BARRIER(icomm,ierrmpi)

    if(stagger_grid) call end_comm_faces

    call get_level_range

    ! Update sfc array: igrid and ipe info in space filling curve
    call amr_Morton_order

    call load_balance

    ! Rebuild tree connectivity
    call getigrids
    call build_connectivity

    ! Update the list of active grids
    call selectgrids
    !  grid structure now complete again.

    ! since we only filled mesh values, and advance assumes filled
    ! ghost cells, do boundary filling for the new levels
    if (time_advance) then
       call getbc(global_time+dt,0.d0,ps,iwstart,nwgc)
    else
       call getbc(global_time,0.d0,ps,iwstart,nwgc)
    end if

    if (use_multigrid) call mg_update_refinement(n_coarsen, n_refine)

    if (associated(usr_after_refine)) then
       call usr_after_refine(n_coarsen, n_refine)
    end if

    !$acc update device(coarsen, refine)

  end subroutine amr_coarsen_refine

  !> For all grids on all processors, do a check on refinement flags. Make
  !> sure that neighbors will not differ more than one level of refinement.
  subroutine proper_nesting
    use mod_forest
    use mod_global_parameters
    use mod_amr_neighbors, only: find_neighbor

    logical, dimension(:,:), allocatable :: refine2
    integer :: iigrid, igrid, level, ic1,ic2,ic3, inp1,inp2,inp3, i1,i2,i3,&
         my_neighbor_type,ipe
    logical :: coarsening, pole(ndim), sendbuf(max_blocks)
    type(tree_node_ptr) :: tree, p_neighbor, my_parent, sibling, my_neighbor,&
         neighborchild

    if (nbufferx1/=0.or.nbufferx2/=0.or.nbufferx3/=0) then
       allocate(refine2(max_blocks,npe))
       call MPI_ALLREDUCE(refine,refine2,max_blocks*npe,MPI_LOGICAL,MPI_LOR,&
            icomm,ierrmpi)
       refine=refine2
    else
       sendbuf(:)=refine(:,mype)
       call MPI_ALLGATHER(sendbuf,max_blocks,MPI_LOGICAL,refine,max_blocks,&
            MPI_LOGICAL,icomm,ierrmpi)
    end if

    do level=min(levmax,refine_max_level-1),levmin+1,-1
       tree%node => level_head(level)%node
       do
          if (.not.associated(tree%node)) exit

          if (refine(tree%node%igrid,tree%node%ipe)) then
             ic1=1+modulo(tree%node%ig1-1,2);ic2=1+modulo(tree%node%ig2-1,2)
             ic3=1+modulo(tree%node%ig3-1,2);
             do inp3=ic3-2,ic3-1
                do inp2=ic2-2,ic2-1
                   do inp1=ic1-2,ic1-1
                      if (inp1==0.and.inp2==0.and.inp3==0) cycle
                      p_neighbor%node => tree%node%parent%node
                      if (inp1/=0) then
                         p_neighbor%node => p_neighbor%node%neighbor(ic1,1)%node
                         if (.not.associated(p_neighbor%node)) cycle
                      end if
                      if (inp2/=0) then
                         p_neighbor%node => p_neighbor%node%neighbor(ic2,2)%node
                         if (.not.associated(p_neighbor%node)) cycle
                      end if
                      if (inp3/=0) then
                         p_neighbor%node => p_neighbor%node%neighbor(ic3,3)%node
                         if (.not.associated(p_neighbor%node)) cycle
                      end if
                      if (p_neighbor%node%leaf) then
                         refine(p_neighbor%node%igrid,p_neighbor%node%ipe)=.true.
                      end if
                   end do
                end do
             end do
          end if

          tree%node => tree%node%next%node
       end do
    end do

    ! On each processor locally, check if grids set for coarsening are already
    ! set for refinement.

    do iigrid=1,igridstail; igrid=igrids(iigrid);
       if (refine(igrid,mype).and.coarsen(igrid,mype)) coarsen(igrid,&
            mype)=.false.
    end do

    ! For all grids on all processors, do a check on coarse refinement flags
    sendbuf(:)=coarsen(:,mype)
    call MPI_ALLGATHER(sendbuf,max_blocks,MPI_LOGICAL,coarsen,max_blocks,&
         MPI_LOGICAL,icomm,ierrmpi)

    do level=levmax,max(2,levmin),-1
       tree%node => level_head(level)%node
       do
          if (.not.associated(tree%node)) exit

          if (coarsen(tree%node%igrid,tree%node%ipe)) then
             coarsening=.true.
             my_parent%node => tree%node%parent%node

             ! are all siblings flagged for coarsen ?
             check1:  do ic3=1,2
                do ic2=1,2
                   do ic1=1,2
                      sibling%node => my_parent%node%child(ic1,ic2,ic3)%node
                      if (sibling%node%leaf) then
                         if (coarsen(sibling%node%igrid,sibling%node%ipe)) cycle
                      end if
                      call unflag_coarsen_siblings
                      exit check1
                   end do
                end do
             end do check1

             ! Make sure that neighbors will not differ more than one level of
             ! refinement, otherwise unflag all siblings
             if (coarsening) then
                check2:     do ic3=1,2
                   do ic2=1,2
                      do ic1=1,2
                         sibling%node => my_parent%node%child(ic1,ic2,ic3)%node
                         do i3=ic3-2,ic3-1
                            do i2=ic2-2,ic2-1
                               do i1=ic1-2,ic1-1
                                  if (i1==0.and.i2==0.and.i3==0) cycle
                                  call find_neighbor(my_neighbor,my_neighbor_type, sibling,&
                                       i1,i2,i3,pole)
                                  select case (my_neighbor_type)
                                  case (neighbor_sibling)
                                     if (refine(my_neighbor%node%igrid,&
                                          my_neighbor%node%ipe)) then
                                        call unflag_coarsen_siblings
                                        exit check2
                                     else
                                        cycle
                                     end if
                                  case (neighbor_fine)
                                     neighborchild%node=>my_neighbor%node%child(1,1,&
                                          1)%node
                                     if (neighborchild%node%leaf) then
                                        if (coarsen(neighborchild%node%igrid,&
                                             neighborchild%node%ipe)) then
                                           cycle
                                        end if
                                     end if
                                     call unflag_coarsen_siblings
                                     exit check2
                                  end select
                               end do
                            end do
                         end do
                      end do
                   end do
                end do check2
             end if

          end if

          tree%node => tree%node%next%node
       end do
    end do

  contains

    subroutine unflag_coarsen_siblings

      integer :: ic1,ic2,ic3
      type(tree_node_ptr) :: sibling

      do ic3=1,2
         do ic2=1,2
            do ic1=1,2
               sibling%node => my_parent%node%child(ic1,ic2,ic3)%node
               if (sibling%node%leaf) then
                  coarsen(sibling%node%igrid,sibling%node%ipe)=.false.
               end if
            end do
         end do
      end do
      coarsening=.false.

    end subroutine unflag_coarsen_siblings

  end subroutine proper_nesting

  !> coarsen sibling blocks into one block
  subroutine coarsen_grid_siblings(igrid,ipe,child_igrid,child_ipe,active)
    use mod_global_parameters
    use mod_coarsen, only: coarsen_grid
    use mod_initialize_amr, only: initial_condition
    use mod_amr_solution_node, only: alloc_node
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: igrid, ipe
    integer, dimension(2,2,2), intent(in) :: child_igrid, child_ipe
    logical, intent(in) :: active

    integer :: igridFi, ipeFi, ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
         ixComax3, ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
         ixCoMmin1,ixCoMmin2,ixCoMmin3,ixCoMmax1,ixCoMmax2,ixCoMmax3, ic1,ic2,&
         ic3, idir
    integer :: ix1, ix2, ix3, iw

    if (ipe==mype) call alloc_node(igrid)

    ! New passive cell, coarsen from initial condition:
    if (.not. active) then
       if (ipe == mype) then
          ! initial_condition fetches this block's positions back itself
          call initial_condition(igrid)
          do ic3=1,2
             do ic2=1,2
                do ic1=1,2
                   igridFi=child_igrid(ic1,ic2,ic3)
                   ipeFi=child_ipe(ic1,ic2,ic3)
                   !if (ipeFi==mype) then
                   !   ! remove solution space of child
                   !   call dealloc_node(igridFi)
                   !end if
                end do
             end do
          end do
       end if
       return
    end if

    do ic3=1,2
       do ic2=1,2
          do ic1=1,2
             igridFi=child_igrid(ic1,ic2,ic3)
             ipeFi=child_ipe(ic1,ic2,ic3)

             if (ipeFi==mype) then
                dxlevel(1)=rnode(rpdx1_,igridFi);dxlevel(2)=rnode(rpdx2_,igridFi)
                dxlevel(3)=rnode(rpdx3_,igridFi);
                if (ipe==mype) then
                   ixComin1=ixMlo1+(ic1-1)*(ixMhi1-ixMlo1+1)/2
                   ixComin2=ixMlo2+(ic2-1)*(ixMhi2-ixMlo2+1)/2
                   ixComin3=ixMlo3+(ic3-1)*(ixMhi3-ixMlo3+1)/2;
                   ixComax1=ixMhi1+(ic1-2)*(ixMhi1-ixMlo1+1)/2
                   ixComax2=ixMhi2+(ic2-2)*(ixMhi2-ixMlo2+1)/2
                   ixComax3=ixMhi3+(ic3-2)*(ixMhi3-ixMlo3+1)/2;

                   call coarsen_grid(ps(igridFi),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
                        ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid),&
                        ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixComin1,ixComin2,&
                        ixComin3,ixComax1,ixComax2,ixComax3,bgeo,igridFi,bgeo,igrid)
                   ! remove solution space of child
                   !call dealloc_node(igridFi)
                else
                   ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1;
                   ixCoGmax1=ixGhi1/2+nghostcells;ixCoGmax2=ixGhi2/2+nghostcells
                   ixCoGmax3=ixGhi3/2+nghostcells;
                   ixCoMmin1=ixCoGmin1+nghostcells;ixCoMmin2=ixCoGmin2+nghostcells
                   ixCoMmin3=ixCoGmin3+nghostcells;ixCoMmax1=ixCoGmax1-nghostcells
                   ixCoMmax2=ixCoGmax2-nghostcells;ixCoMmax3=ixCoGmax3-nghostcells;
                   call coarsen_grid(ps(igridFi),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
                        ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,psc(igridFi),&
                        ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
                        ixCoMmin1,ixCoMmin2,ixCoMmin3,ixCoMmax1,ixCoMmax2,ixCoMmax3,&
                        bgeo,igridFi,bgeoc,igridFi)

                   ! Record the chunk; exchange_coarsened_blocks packs and
                   ! posts them all once the walk is over, one message per
                   ! peer.  The key is built from the receiving parent's block
                   ! index, igrid, which the receiver derives from its own.
                   isend=isend+1
                   n_cf_snd=n_cf_snd+1
                   if (n_cf_snd > max_buff) then
                      call mpistop('coarsen_grid_siblings: max_buff too small in send')
                   end if
                   cf_snd_igrid(n_cf_snd) = igridFi
                   cf_snd_dest(n_cf_snd)  = ipe
                   cf_snd_key(n_cf_snd)   = 8*(igrid-1) + (ic1-1) + 2*(ic2-1)&
                       + 4*(ic3-1)

                   ! The staggered branch below is NOT aggregated, and carries
                   ! two defects that are recorded rather than repaired because
                   ! nothing reaches it (stagger_grid is false throughout this
                   ! fork and fix_edges has no caller): the do idir loop stores
                   ! every request into the single slot sendrequest_stg(isend),
                   ! leaking ndim-1 handles that are never waited on; and
                   ! itag_stg multiplies igridFi by 3, 4 or 5 for ndir=3, so
                   ! (igridFi=4,idir=1) collides with (igridFi=3,idir=2).
                   if(stagger_grid) then
                      do idir=1,ndim
                         itag_stg=(npe+ipeFi+1)+igridFi*(ndir-1+idir)
                         call mpi_isend_wrapper(psc(igridFi)%ws,1,type_coarse_block_stg(idir,&
                              ic1,ic2,ic3),ipe,itag_stg, icomm,sendrequest_stg(isend),&
                              ierrmpi)
                      end do
                   end if
                end if
             else
                if (ipe==mype) then
                   irecv=irecv+1
                   n_cf_rcv=n_cf_rcv+1
                   if (n_cf_rcv > max_buff) then
                      call mpistop('coarsen_grid_siblings: max_buff too small in receive')
                   end if
                   cf_rcv_igrid(n_cf_rcv) = igrid
                   cf_rcv_src(n_cf_rcv)   = ipeFi
                   cf_rcv_ic(1,n_cf_rcv)  = ic1
                   cf_rcv_ic(2,n_cf_rcv)  = ic2
                   cf_rcv_ic(3,n_cf_rcv)  = ic3
                   ! the sender builds this same key from the block index it
                   ! was handed for this parent
                   cf_rcv_key(n_cf_rcv)   = 8*(igrid-1) + (ic1-1) + 2*(ic2-1)&
                       + 4*(ic3-1)
                   if(stagger_grid) then
                      do idir=1,ndim
                         itag_stg=(npe+ipeFi+1)+igridFi*(ndir-1+idir)
                         call mpi_irecv_wrapper(ps(igrid)%ws,1,type_sub_block_stg(idir,ic1,ic2,&
                              ic3),ipeFi,itag_stg, icomm,recvrequest_stg(irecv),ierrmpi)
                      end do
                   end if
                end if
             end if
          end do
       end do
    end do

  end subroutine coarsen_grid_siblings

  !> Exchange and apply the non-local half of a coarsening step, one message
  !> per neighbouring rank.
  !>
  !> coarsen_grid_siblings has already written each off-rank child into its own
  !> coarse representative psc(igridFi) - which is bgc(1)%w(...,igridFi), so it
  !> is on the device - and recorded a descriptor for it.  Here those
  !> descriptors are laid out as one contiguous run per peer, ordered by the
  !> shared key, packed by a single kernel and sent as a single message each.
  subroutine exchange_coarsened_blocks
    use mod_global_parameters
    use mod_msg_layout, only: layout_runs
    use mod_comm_lib, only: mpistop

    integer :: k, nxCo1, nxCo2, nxCo3, nchunk, nbuf
    integer :: ixCoMmin1, ixCoMmin2, ixCoMmin3
    integer :: igrid, ibuf, ic1, ic2, ic3, iw, ix1, ix2, ix3
    integer, allocatable :: chunksize(:)

    nxCo1 = block_nx1/2; nxCo2 = block_nx2/2; nxCo3 = block_nx3/2
    ! only 1:nwgc is coarsened and only 1:nwgc is applied - the analytic extras
    ! past it are re-derived by alloc_node - so nwgc, not nw, is what travels
    nchunk = nxCo1*nxCo2*nxCo3*nwgc
    ! psc is indexed from 1 in every direction, like bgc(1)%w
    ixCoMmin1 = 1+nghostcells
    ixCoMmin2 = 1+nghostcells
    ixCoMmin3 = 1+nghostcells

    ! Every chunk is one coarsened child, so they are all the same size.
    allocate(chunksize(max(n_cf_snd,n_cf_rcv,1)))
    chunksize = nchunk
    call layout_runs(n_cf_snd, cf_snd_dest, cf_snd_key, chunksize, cf_snd_ibuf,&
       n_cf_send_pe, cf_send_pe, cf_send_pe_off, cf_send_pe_len)
    call layout_runs(n_cf_rcv, cf_rcv_src, cf_rcv_key, chunksize, cf_rcv_ibuf,&
       n_cf_recv_pe, cf_recv_pe, cf_recv_pe_off, cf_recv_pe_len)
    deallocate(chunksize)

    if (n_cf_snd > 0) then
      !$acc update device(cf_snd_igrid(1:n_cf_snd), cf_snd_ibuf(1:n_cf_snd))
    end if
    if (n_cf_rcv > 0) then
      !$acc update device(cf_rcv_igrid(1:n_cf_rcv), cf_rcv_ibuf(1:n_cf_rcv),&
      !$acc&              cf_rcv_ic(:,1:n_cf_rcv))
    end if

    ! One Irecv per peer, straight into that peer's run.  A peer sends exactly
    ! one message, so (communicator, source) already disambiguates and the tag
    ! carries nothing.
    itag = 0
    if (n_cf_recv_pe > 0) then
#ifndef NOGPUDIRECT
      !$acc host_data use_device(rcv_buff_cf)
#endif
      do k = 1, n_cf_recv_pe
        call mpi_irecv_wrapper(rcv_buff_cf(cf_recv_pe_off(k)),&
           cf_recv_pe_len(k),MPI_DOUBLE_PRECISION,cf_recv_pe(k),itag,icomm,&
           recvrequest(k),ierrmpi)
      end do
#ifndef NOGPUDIRECT
      !$acc end host_data
#endif
    end if

    ! Pack every outgoing chunk in one kernel.  bgc(1)%w carries the grid index
    ! last, so the block can be selected by a device-side index.
    if (n_cf_snd > 0) then
      !$acc parallel loop gang default(present) private(igrid,ibuf)
      do k = 1, n_cf_snd
         igrid = cf_snd_igrid(k)
         ibuf  = cf_snd_ibuf(k)
         !$acc loop collapse(4) vector
         do iw = 1, nwgc
            do ix3 = 1, nxCo3
               do ix2 = 1, nxCo2
                  do ix1 = 1, nxCo1
                     snd_buff_cf(ibuf + (ix1-1) + nxCo1*(ix2-1) &
                          + nxCo1*nxCo2*(ix3-1) &
                          + nxCo1*nxCo2*nxCo3*(iw-1)) = &
                          bgc(1)%w(ixCoMmin1-1+ix1, ixCoMmin2-1+ix2,&
                                   ixCoMmin3-1+ix3, iw, igrid)
                  end do
               end do
            end do
         end do
      end do
    end if

    if (n_cf_send_pe > 0) then
#ifdef NOGPUDIRECT
      !$acc update host(snd_buff_cf(1:n_cf_snd*nchunk))
#else
      !$acc host_data use_device(snd_buff_cf)
#endif
      do k = 1, n_cf_send_pe
        call mpi_isend_wrapper(snd_buff_cf(cf_send_pe_off(k)),&
           cf_send_pe_len(k),MPI_DOUBLE_PRECISION,cf_send_pe(k),itag,icomm,&
           sendrequest(k),ierrmpi)
      end do
#ifndef NOGPUDIRECT
      !$acc end host_data
#endif
    end if

    if (n_cf_recv_pe > 0) then
      call MPI_WAITALL(n_cf_recv_pe,recvrequest,recvstatus,ierrmpi)
      ! With one aggregated message per peer a disagreement about which chunks
      ! the run holds is no longer an MPI mismatch, so check the arrival
      ! lengths.  A difference in order cannot arise: both ranks sort the same
      ! run by the same key.
      do k = 1, n_cf_recv_pe
        call MPI_GET_COUNT(recvstatus(:,k),MPI_DOUBLE_PRECISION,nbuf,ierrmpi)
        if (nbuf /= cf_recv_pe_len(k)) call mpistop( &
           "exchange_coarsened_blocks: message length disagrees with the layout")
      end do
#ifdef NOGPUDIRECT
      !$acc update device(rcv_buff_cf(1:n_cf_rcv*nchunk))
#endif
    end if

    ! Apply every incoming chunk in one kernel.
    if (n_cf_rcv > 0) then
      !$acc parallel loop gang default(present) private(igrid,ibuf,ic1,ic2,ic3)
      do k = 1, n_cf_rcv
         igrid = cf_rcv_igrid(k)
         ibuf  = cf_rcv_ibuf(k)
         ic1   = cf_rcv_ic(1,k)
         ic2   = cf_rcv_ic(2,k)
         ic3   = cf_rcv_ic(3,k)
         !$acc loop collapse(4) vector
         do iw = 1, nwgc  ! analytic extras past nwgc are set in alloc_node
            do ix3 = 1, nxCo3
               do ix2 = 1, nxCo2
                  do ix1 = 1, nxCo1
                     bg(1)%w( &
                          ixMlo1-1+(ic1-1)*nxCo1 + ix1, &
                          ixMlo2-1+(ic2-1)*nxCo2 + ix2, &
                          ixMlo3-1+(ic3-1)*nxCo3 + ix3, &
                          iw, igrid) &
                          = rcv_buff_cf(ibuf + (ix1-1) + nxCo1*(ix2-1) &
                            + nxCo1*nxCo2*(ix3-1) &
                            + nxCo1*nxCo2*nxCo3*(iw-1))
                  end do
               end do
            end do
         end do
      end do
    end if

    if (n_cf_send_pe > 0) then
      call MPI_WAITALL(n_cf_send_pe,sendrequest,sendstatus,ierrmpi)
    end if

  end subroutine exchange_coarsened_blocks

end module mod_coarsen_refine
