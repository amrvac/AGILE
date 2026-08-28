module mod_amr_solution_node
  use mod_comm_lib, only: mpistop

  implicit none
  private

  public :: getnode, putnode
  public :: alloc_node, alloc_state
  public :: dealloc_node

contains


  !> Get first available igrid on processor ipe
  integer function getnode(ipe)
    use mod_forest, only: igrid_inuse
    use mod_global_parameters

    integer, intent(in) :: ipe
    integer :: igrid, igrid_available
  
    igrid_available=0
  
    do igrid=1,max_blocks
       if (igrid_inuse(igrid,ipe)) cycle
  
       igrid_available=igrid
       exit
    end do
  
    if (igrid_available == 0) then
       getnode = -1
       print *, "Current maximum number of grid blocks:", max_blocks
       call mpistop&
          ("Insufficient grid blocks; increase max_blocks in meshlist")
    else
       getnode=igrid_available
       igrid_inuse(igrid,ipe)=.true.
    end if
  
    if (ipe==mype) then
       ! initialize node on host and device
       node(1:nodehi,getnode) = 0
       !$acc update device(node(1:nodehi,getnode))
       rnode(1:rnodehi,getnode) = zero
       !$acc update device(rnode(1:rnodehi,getnode))
    end if
  
  end function getnode
  
  ! put igrid on processor ipe to be not in use
  subroutine putnode(igrid,ipe)
    use mod_forest
    implicit none
  
    integer, intent(in) :: igrid, ipe
  
    igrid_inuse(igrid,ipe)=.false.
  
  end subroutine putnode
  
  !> allocate arrays on igrid node
  subroutine alloc_node(igrid)
#ifdef _OPENACC
    use acc_utils
#endif    
    use mod_forest
    use mod_global_parameters
    use mod_geometry
    use mod_physics, only: phys_set_equi_vars
    use mod_b0, only: set_B0_grid 
    
    integer, intent(in) :: igrid
  
    integer :: level, ig1,ig2,ig3, ign1,ign2,ign3, ixCoGmin1,ixCoGmin2,&
       ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3, i1,i2,i3

    ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1;
    ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells
    ixCoGmax2=(ixGhi2-2*nghostcells)/2+2*nghostcells
    ixCoGmax3=(ixGhi3-2*nghostcells)/2+2*nghostcells;
  
    ! set level information
    level=igrid_to_node(igrid,mype)%node%level
  
    if(.not. associated(ps(igrid)%w)) then
       
       ! allocate arrays for solution and space
       call alloc_state(igrid, ps(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
          ixGhi3, .true.)
       ! allocate arrays for one level coarser solution
       call alloc_state_coarse(igrid, psc(igrid), ixCoGmin1,ixCoGmin2,&
          ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3)
       if(.not.convert) then
        ! allocate arrays for temp solution 1
        call alloc_state(igrid, ps1(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
           ixGhi3, .false.)

        ! allocate temporary solution space
        select case (t_integrator)
        case(ssprk3,ssprk4,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
          call alloc_state(igrid, ps2(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, .false.)
        case(RK3_BT,rk4,ssprk5,IMEX_CB3a)
          call alloc_state(igrid, ps2(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, .false.)
          call alloc_state(igrid, ps3(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, .false.)
        case(IMEX_ARS3,IMEX_232)
          call alloc_state(igrid, ps2(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, .false.)
          call alloc_state(igrid, ps3(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, .false.)
          call alloc_state(igrid, ps4(igrid), ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3, .false.)
        end select
      end if
  
    end if

    ! avoid dividing by zero rho in skipped corner ghostcells when phys_req_diagonal=F
    ps(igrid)%w(:,:,:,1)=1.d0
    ps(igrid)%level=level
    psc(igrid)%level=level-1
    ! avoid dividing by zero rho in skipped corner ghostcells when phys_req_diagonal=F
    psc(igrid)%w(:,:,:,1)=1.d0
    if(phys_trac) ps(igrid)%special_values=0.d0
    if(.not.convert) then
      ps1(igrid)%level=level
      select case (t_integrator)
      case(ssprk3,ssprk4,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
        ps2(igrid)%level=level
      case(RK3_BT,rk4,ssprk5,IMEX_CB3a)
        ps2(igrid)%level=level
        ps3(igrid)%level=level
      case(IMEX_ARS3,IMEX_232)
        ps2(igrid)%level=level
        ps3(igrid)%level=level
        ps4(igrid)%level=level
      end select
    end if

    ! block pointer to current block
    block=>ps(igrid)
    ig1=igrid_to_node(igrid,mype)%node%ig1
    ig2=igrid_to_node(igrid,mype)%node%ig2
    ig3=igrid_to_node(igrid,mype)%node%ig3;
    node(plevel_,igrid)=level
    node(pig1_,igrid)=ig1
    node(pig2_,igrid)=ig2
    node(pig3_,igrid)=ig3
 !$acc update device(node(plevel_,igrid),node(pig1_,igrid),node(pig2_,igrid),node(pig3_,igrid))
    
    ! set dx information
    rnode(rpdx1_,igrid)=dx(1,level)
    rnode(rpdx2_,igrid)=dx(2,level)
    rnode(rpdx3_,igrid)=dx(3,level)
    dxlevel(:)=dx(:,level)
 !$acc update device(rnode(rpdx1_,igrid),rnode(rpdx2_,igrid),rnode(rpdx3_,igrid), dxlevel)

    ! uniform cartesian case as well as all unstretched coordinates
    ! determine the minimal and maximal corners
    rnode(rpxmin1_,igrid)=xprobmin1+dble(ig1-1)*dg1(level)
    rnode(rpxmin2_,igrid)=xprobmin2+dble(ig2-1)*dg2(level)
    rnode(rpxmin3_,igrid)=xprobmin3+dble(ig3-1)*dg3(level)
    rnode(rpxmax1_,igrid)=xprobmin1+dble(ig1)*dg1(level)
    rnode(rpxmax2_,igrid)=xprobmin2+dble(ig2)*dg2(level)
    rnode(rpxmax3_,igrid)=xprobmin3+dble(ig3)*dg3(level)
   if(rnode(rpxmax1_,igrid)>xprobmax1) rnode(rpxmax1_,igrid)=xprobmax1
   if(rnode(rpxmax2_,igrid)>xprobmax2) rnode(rpxmax2_,igrid)=xprobmax2
   if(rnode(rpxmax3_,igrid)>xprobmax3) rnode(rpxmax3_,igrid)=xprobmax3

 !$acc update device( rnode(rpxmax1_,igrid),rnode(rpxmax2_,igrid),rnode(rpxmax3_,igrid), rnode(rpxmin1_,igrid),rnode(rpxmin2_,igrid),rnode(rpxmin3_,igrid) )
   
    ! Fill this block's cell metrics. A uniform block is an analytic function
    ! of rnode alone, so this runs on the device and only the positions come
    ! back; see sync_geometry_host for the rest.
    call fill_geometry_device(igrid)

#:if defined('FILL_NWEXTRA_ANALYTIC')
    ! Fill this block's analytic extra variables (past nwgc) the same way, from
    ! the positions fill_geometry_device just wrote. This is the only place
    ! they are set: getbc, prolongation and coarsening all stop at nwgc, so
    ! every block that reaches alloc_node - a fresh root, a refined child, a
    ! coarsened parent, a load-balanced arrival - gets them here and nothing
    ! overwrites them afterwards.
    call fill_nwextra_device(igrid)
#:endif

    ! initialize background non-evolving solution; both of these read
    ! ps(igrid)%x on the host, so this block is the one place in alloc_node
    ! that still needs the positions back
    if (B0field .or. number_equi_vars>0) call sync_positions_host(igrid)
    if (B0field) call set_B0_grid(igrid)
    if (number_equi_vars>0) call phys_set_equi_vars(igrid)
  
    ! find the blocks on the boundaries
    ps(igrid)%is_physical_boundary=.false.
    
    do i1=-1,1
      if(i1==0) cycle
      ign1=ig1+i1
      ! blocks at periodic boundary have neighbors in the physical domain
      ! thus threated at internal blocks with no physical boundary
      if (periodB(1)) ign1=1+modulo(ign1-1,ng1(level))
      if (ign1 > ng1(level)) then
         if(phi_ > 0 .and. poleB(2,1)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*1)=.false.
         else
           ps(igrid)%is_physical_boundary(2*1)=.true.
         end if
      else if (ign1 < 1) then
         if(phi_ > 0 .and. poleB(1,1)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*1-1)=.false.
         else
           ps(igrid)%is_physical_boundary(2*1-1)=.true.
         end if
      end if
    end do
    
    do i2=-1,1
      if(i2==0) cycle
      ign2=ig2+i2
      ! blocks at periodic boundary have neighbors in the physical domain
      ! thus threated at internal blocks with no physical boundary
      if (periodB(2)) ign2=1+modulo(ign2-1,ng2(level))
      if (ign2 > ng2(level)) then
         if(phi_ > 0 .and. poleB(2,2)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*2)=.false.
         else
           ps(igrid)%is_physical_boundary(2*2)=.true.
         end if
      else if (ign2 < 1) then
         if(phi_ > 0 .and. poleB(1,2)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*2-1)=.false.
         else
           ps(igrid)%is_physical_boundary(2*2-1)=.true.
         end if
      end if
    end do
    
    
    do i3=-1,1
      if(i3==0) cycle
      ign3=ig3+i3
      ! blocks at periodic boundary have neighbors in the physical domain
      ! thus threated at internal blocks with no physical boundary
      if (periodB(3)) ign3=1+modulo(ign3-1,ng3(level))
      if (ign3 > ng3(level)) then
         if(phi_ > 0 .and. poleB(2,3)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*3)=.false.
         else
           ps(igrid)%is_physical_boundary(2*3)=.true.
         end if
      else if (ign3 < 1) then
         if(phi_ > 0 .and. poleB(1,3)) then
           ! if at a pole, the boundary is not physical boundary
           ps(igrid)%is_physical_boundary(2*3-1)=.false.
         else
           ps(igrid)%is_physical_boundary(2*3-1)=.true.
         end if
      end if
    end do
    
    if(any(ps(igrid)%is_physical_boundary)) then
      phyboundblock(igrid)=.true.
    else
      phyboundblock(igrid)=.false.
   end if

   !$acc update device( phyboundblock(igrid) )
#ifdef _OPENACC
   call copy_or_update(ps(igrid)%igrid) 
   call copy_or_update(ps1(igrid)%igrid) 
   call copy_or_update(ps2(igrid)%igrid)

   call copy_or_update_pointer(ps(igrid)%w, no_update=.true.)
   call copy_or_update_pointer(ps1(igrid)%w, no_update=.true.)
   call copy_or_update_pointer(ps2(igrid)%w, no_update=.true.)

   call copy_or_update_pointer(ps(igrid)%is_physical_boundary)
   call copy_or_update_pointer(ps1(igrid)%is_physical_boundary, no_update=.true.)
   call copy_or_update_pointer(ps2(igrid)%is_physical_boundary, no_update=.true.)
   call copy_or_update_pointer(psc(igrid)%w, no_update=.true.)
   ! The metrics are already resident on the device for every block (bgeo and
   ! bgeoc are allocated once in initialize_vars) and were filled above, on
   ! whichever side fill_geometry_* chose; ps1/ps2/psc share them, so there is
   ! nothing to do per state here.
#endif

 end subroutine alloc_node
  
  !> Build block igrid's cell metrics directly on the device.
  !>
  !> A uniform block is an analytic function of nothing but its corner and its
  !> spacing, both of which alloc_node has just put in rnode, so there is no
  !> reason to compute a dozen block-sized arrays on the CPU and ship them
  !> across the bus every time a node is (re)allocated - which, with AMR, is
  !> most regrids. Everything is built here instead, for whichever coordinate
  !> system the build selected, and nothing is handed back to the host - see
  !> sync_positions_host and sync_geometry_host in mod_geometry for that.
  !>
  !> The expressions reproduce what get_surface_area and alloc_node's
  !> select(coordinate) used to do term for term, including the two distinct
  !> radial coordinates: x uses the two-sided form (measured from rnode's lower
  !> corner in the mesh and its upper corner in the outer ghost layer, so the
  !> overlapping ghost cells of neighbouring blocks agree to the last bit),
  !> while the volumes and cell sizes use the one-sided form over the extended
  !> range that the old code called xext.
  subroutine fill_geometry_device(igrid)
    use mod_global_parameters

    integer, intent(in) :: igrid

    integer          :: ix1,ix2,ix3
    integer          :: ixCoGmax1,ixCoGmax2,ixCoGmax3
    ! last cell in each direction whose centre is measured from the lower corner
    integer          :: iup1,iup2,iup3
    ! cell spacing of the block and of its coarse representative
    double precision :: d1,d2,d3, c1,c2,c3
    double precision :: xmin1,xmin2,xmin3, xmax1,xmax2,xmax3
    ! cell centre (p), first-cell centre (q) and extended-range centre (e)
    double precision :: p1,p2,p3, q1,q2, e1,e2

    d1=rnode(rpdx1_,igrid); d2=rnode(rpdx2_,igrid); d3=rnode(rpdx3_,igrid)
    c1=2.0d0*d1; c2=2.0d0*d2; c3=2.0d0*d3
    xmin1=rnode(rpxmin1_,igrid); xmin2=rnode(rpxmin2_,igrid)
    xmin3=rnode(rpxmin3_,igrid)
    xmax1=rnode(rpxmax1_,igrid); xmax2=rnode(rpxmax2_,igrid)
    xmax3=rnode(rpxmax3_,igrid)

    iup1=ixMhi1-nghostcells; iup2=ixMhi2-nghostcells; iup3=ixMhi3-nghostcells

    ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells
    ixCoGmax2=(ixGhi2-2*nghostcells)/2+2*nghostcells
    ixCoGmax3=(ixGhi3-2*nghostcells)/2+2*nghostcells

    ! q1 and q2 are the centres of the first cell in the radial and polar
    ! directions.  They do not depend on the loop indices, so they are scalars
    ! computed once here and read by the kernel, rather than recomputed inside
    ! it.  (ixGlo is always inside the mesh, so only the lower-corner form of
    ! the two-sided rule below can apply.)
    q1 = xmin1 + (dble(ixGlo1-nghostcells)-half)*d1
    q2 = xmin2 + (dble(ixGlo2-nghostcells)-half)*d2

    ! ---- positions, and (curvilinear only) the face areas built from them ----
    ! The three "if (ix==ixGlo)" branches fill the lower face of the first
    ! cell in each direction, which the flux update reaches as surfaceC(ix-1,.)
    ! A Cartesian build has no face areas to fill: its flux update divides by
    ! the uniform rnode spacing instead, so surfaceC is never allocated.
    ! The cell centre is measured from the block's lower corner inside the
    ! mesh, and from its upper corner in the outer ghost layer, so that the
    ! overlapping ghost cells of two neighbouring blocks come out identical to
    ! the last bit.  Written out here rather than called as a helper: as an
    ! !$acc routine seq inside this loop, nvfortran's -Minline hoisted the call
    ! and gave every cell of the block the position of the first one.
    !$acc parallel loop collapse(3) default(present) private(p1,p2,p3)
    do ix3=ixGlo3,ixGhi3
       do ix2=ixGlo2,ixGhi2
          do ix1=ixGlo1,ixGhi1
             if (ix1 <= iup1) then
                p1 = xmin1 + (dble(ix1-nghostcells)-half)*d1
             else
                p1 = xmax1 + (dble(ix1-ixMhi1)-half)*d1
             end if
             if (ix2 <= iup2) then
                p2 = xmin2 + (dble(ix2-nghostcells)-half)*d2
             else
                p2 = xmax2 + (dble(ix2-ixMhi2)-half)*d2
             end if
             if (ix3 <= iup3) then
                p3 = xmin3 + (dble(ix3-nghostcells)-half)*d3
             else
                p3 = xmax3 + (dble(ix3-ixMhi3)-half)*d3
             end if
             bgeo%x(ix1,ix2,ix3,1,igrid) = p1
             bgeo%x(ix1,ix2,ix3,2,igrid) = p2
             bgeo%x(ix1,ix2,ix3,3,igrid) = p3

#:if GEOM != 'Cartesian'
#:if GEOM == 'cylindrical'
             ! 3D cylindrical is (r, z, phi): set_coordinate_system fixes
             ! r_=1, z_=2, phi_=3, so the z and phi faces are named directly
             bgeo%surfaceC(ix1,ix2,ix3,1,igrid) = dabs(p1+half*d1)*d2*d3
             bgeo%surfaceC(ix1,ix2,ix3,2,igrid) = p1*d1*d3
             bgeo%surfaceC(ix1,ix2,ix3,3,igrid) = d1*d2
             if (ix1==ixGlo1) bgeo%surfaceC(ixGlo1-1,ix2,ix3,1,igrid) = &
                  dabs(q1-half*d1)*d2*d3
             ! the z and phi face areas do not depend on the index being
             ! stepped back, so the extra plane repeats the same expression
             if (ix2==ixGlo2) bgeo%surfaceC(ix1,ixGlo2-1,ix3,2,igrid) = p1*d1*d3
             if (ix3==ixGlo3) bgeo%surfaceC(ix1,ix2,ixGlo3-1,3,igrid) = d1*d2
#:else
             bgeo%surfaceC(ix1,ix2,ix3,1,igrid) = (p1+half*d1)**2 &
                  *two*dsin(p2)*dsin(half*d2)*d3
             bgeo%surfaceC(ix1,ix2,ix3,2,igrid) = p1*d1*dsin(p2+half*d2)*d3
             bgeo%surfaceC(ix1,ix2,ix3,3,igrid) = p1*d1*d2
             if (ix1==ixGlo1) bgeo%surfaceC(ixGlo1-1,ix2,ix3,1,igrid) = &
                  (q1-half*d1)**2*two*dsin(p2)*dsin(half*d2)*d3
             if (ix2==ixGlo2) bgeo%surfaceC(ix1,ixGlo2-1,ix3,2,igrid) = &
                  p1*d1*dsin(q2-half*d2)*d3
             if (ix3==ixGlo3) bgeo%surfaceC(ix1,ix2,ixGlo3-1,3,igrid) = p1*d1*d2
#:endif
#:endif
          end do
       end do
    end do

#:if GEOM != 'Cartesian'
    ! ---- volumes and cell sizes, over the extended range ----
    ! Curvilinear only: in a Cartesian build every one of these is the uniform
    ! rnode spacing, which its kernels read from rnode directly, so none of
    ! dx, dvolume or ds is allocated there.
    !$acc parallel loop collapse(3) default(present) private(e1,e2)
    do ix3=ixGextmin3,ixGextmax3
       do ix2=ixGextmin2,ixGextmax2
          do ix1=ixGextmin1,ixGextmax1
             e1 = xmin1 + (dble(ix1-nghostcells)-half)*d1
             e2 = xmin2 + (dble(ix2-nghostcells)-half)*d2

             bgeo%dx(ix1,ix2,ix3,1,igrid) = d1
             bgeo%dx(ix1,ix2,ix3,2,igrid) = d2
             bgeo%dx(ix1,ix2,ix3,3,igrid) = d3

#:if GEOM == 'cylindrical'
             bgeo%dvolume(ix1,ix2,ix3,igrid) = dabs(e1)*d1*d2*d3
             bgeo%ds(ix1,ix2,ix3,1,igrid)  = d1
             bgeo%ds(ix1,ix2,ix3,2,igrid)  = d2
             bgeo%ds(ix1,ix2,ix3,3,igrid)  = e1*d3
#:else
             bgeo%dvolume(ix1,ix2,ix3,igrid) = (e1**2+d1**2/12.0d0)*d1 &
                  *two*dabs(dsin(e2))*dsin(half*d2)*d3
             bgeo%ds(ix1,ix2,ix3,1,igrid)  = d1
             bgeo%ds(ix1,ix2,ix3,2,igrid)  = e1*d2
             bgeo%ds(ix1,ix2,ix3,3,igrid)  = e1*dsin(e2)*d3
#:endif
          end do
       end do
    end do
#:endif

    ! ---- the one-level-coarser representative ----
    ! Its positions are measured from the lower corner throughout - there is no
    ! neighbouring coarse block whose ghost cells have to match - and its
    ! spacing is doubled, so one loop covers positions, volumes and areas.
    ! As above, only the positions exist in a Cartesian build.
    !$acc parallel loop collapse(3) default(present) private(p1,p2,q1,q2)
    do ix3=1,ixCoGmax3
       do ix2=1,ixCoGmax2
          do ix1=1,ixCoGmax1
             p1 = xmin1 + (dble(ix1-nghostcells)-half)*c1
             p2 = xmin2 + (dble(ix2-nghostcells)-half)*c2
             bgeoc%x(ix1,ix2,ix3,1,igrid) = p1
             bgeoc%x(ix1,ix2,ix3,2,igrid) = p2
             bgeoc%x(ix1,ix2,ix3,3,igrid) = &
                  xmin3 + (dble(ix3-nghostcells)-half)*c3

#:if GEOM != 'Cartesian'
             bgeoc%dx(ix1,ix2,ix3,1,igrid) = c1
             bgeoc%dx(ix1,ix2,ix3,2,igrid) = c2
             bgeoc%dx(ix1,ix2,ix3,3,igrid) = c3
             ! nothing reads the coarse ds today, but leaving it
             ! uninitialised would be a trap for whoever first does
             bgeoc%ds(ix1,ix2,ix3,1,igrid)  = c1
             bgeoc%ds(ix1,ix2,ix3,2,igrid)  = c2
             bgeoc%ds(ix1,ix2,ix3,3,igrid)  = c3
#:endif

#:if GEOM != 'Cartesian'
#:if GEOM == 'cylindrical'
             bgeoc%dvolume(ix1,ix2,ix3,igrid)    = dabs(p1)*c1*c2*c3
             bgeoc%surfaceC(ix1,ix2,ix3,1,igrid) = dabs(p1+half*c1)*c2*c3
             bgeoc%surfaceC(ix1,ix2,ix3,2,igrid) = p1*c1*c3
             bgeoc%surfaceC(ix1,ix2,ix3,3,igrid) = c1*c2
             if (ix1==1) then
                q1 = xmin1 + (dble(1-nghostcells)-half)*c1
                bgeoc%surfaceC(0,ix2,ix3,1,igrid) = dabs(q1-half*c1)*c2*c3
             end if
             if (ix2==1) bgeoc%surfaceC(ix1,0,ix3,2,igrid) = p1*c1*c3
             if (ix3==1) bgeoc%surfaceC(ix1,ix2,0,3,igrid) = c1*c2
#:else
             bgeoc%dvolume(ix1,ix2,ix3,igrid) = (p1**2+c1**2/12.0d0)*c1 &
                  *two*dabs(dsin(p2))*dsin(half*c2)*c3
             bgeoc%surfaceC(ix1,ix2,ix3,1,igrid) = (p1+half*c1)**2 &
                  *two*dsin(p2)*dsin(half*c2)*c3
             bgeoc%surfaceC(ix1,ix2,ix3,2,igrid) = p1*c1*dsin(p2+half*c2)*c3
             bgeoc%surfaceC(ix1,ix2,ix3,3,igrid) = p1*c1*c2
             if (ix1==1) then
                q1 = xmin1 + (dble(1-nghostcells)-half)*c1
                bgeoc%surfaceC(0,ix2,ix3,1,igrid) = (q1-half*c1)**2 &
                     *two*dsin(p2)*dsin(half*c2)*c3
             end if
             if (ix2==1) then
                q2 = xmin2 + (dble(1-nghostcells)-half)*c2
                bgeoc%surfaceC(ix1,0,ix3,2,igrid) = p1*c1*dsin(q2-half*c2)*c3
             end if
             if (ix3==1) bgeoc%surfaceC(ix1,ix2,0,3,igrid) = p1*c1*c2
#:endif
#:endif
          end do
       end do
    end do

    ! Nothing is copied back here. alloc_node runs on every regrid, mostly for
    ! blocks no host routine will look at, so the fetch belongs at the sites
    ! that actually read the geometry: sync_positions_host before host code
    ! that needs ps(igrid)%x, sync_geometry_host before output.

  end subroutine fill_geometry_device

#:if defined('FILL_NWEXTRA_ANALYTIC')
  !> Fill this block's extra w-variables on the device from usr_set_nwextra.
  !>
  !> The `nwextra` variables (registered by var_set_extravar) are not advected,
  !> carry no boundary condition, and here are taken to be analytic functions
  !> of position alone. Rather than exchange them through getbc or interpolate
  !> them in prolongation and coarsening - all of which stop at nwgc - they are
  !> re-derived from the user's usr_set_nwextra in every cell of every block,
  !> the full ixG range with ghost cells included. `alloc_node` is the sole
  !> caller: it runs for every block that comes into existence (a fresh root, a
  !> refined child, a coarsened parent, a load-balanced arrival), right after
  !> fill_geometry_device has written this block's positions, and nothing
  !> overwrites the extra slots afterwards.
  !>
  !> Filling the ghost cells this way is also what lets a build with such a
  !> variable reach the polar axis: bgeo%x in the ghost layer beyond theta=0
  !> (or r=0) carries the mirrored coordinate, so evaluating the user's
  !> analytic field there reproduces on its own the sign flips a vector picks
  !> up across the axis, with no entry in typeboundary and no pole branch in
  !> getbc.
  !>
  !> This flag is set by config_schema.toml's `implies` for phys='ffhd', whose
  !> frozen field b1,b2,b3 is exactly such a variable.
  !>
  !> Runs on the device: bgeo%x and bg(1)%w are both resident for all
  !> max_blocks from initialize_vars, so the kernel indexes them directly.
  !> mod_variables tracks nwextra and, via var_set_extravar, the w indices
  !> iw_extra(1:nwextra) the extra variables actually occupy - this routine
  !> makes no assumption about where in w they live.
  subroutine fill_nwextra_device(igrid)
    use mod_global_parameters
    use mod_usr, only: usr_set_nwextra

    integer, intent(in) :: igrid

    integer          :: ix1, ix2, ix3, iwx
    ! wx is sized by the compile-time max_nw, not the runtime nwextra: a
    ! variable-length private array in an !$acc parallel loop is a fragile
    ! corner across OpenACC compilers. Only wx(1:nwextra) is passed and used.
    double precision :: xloc(1:ndim), wx(max_nw)

    !$acc parallel loop collapse(3) default(present) private(xloc, wx)
    do ix3 = ixGlo3, ixGhi3
       do ix2 = ixGlo2, ixGhi2
          do ix1 = ixGlo1, ixGhi1
             xloc(1:ndim) = bgeo%x(ix1,ix2,ix3,1:ndim,igrid)
             call usr_set_nwextra(xloc, wx(1:nwextra))
             do iwx = 1, nwextra
                bg(1)%w(ix1,ix2,ix3, iw_extra(iwx), igrid) = wx(iwx)
             end do
          end do
       end do
    end do

  end subroutine fill_nwextra_device
#:endif

  !> allocate memory to physical state of igrid node
  subroutine alloc_state(igrid, s, ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
     ixGmax3, alloc_once_for_ps)
    use mod_global_parameters
    type(state) :: s
    integer, intent(in) :: igrid, ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3
    logical, intent(in) :: alloc_once_for_ps
    integer             :: ixGsmin1,ixGsmin2,ixGsmin3,ixGsmax1,ixGsmax2,&
       ixGsmax3
  
    !allocate(s%w(ixG^S,1:nw))
    s%w => bg(s%istep)%w(:,:,:,:,igrid)
    s%igrid=igrid
    s%w=0.d0
    s%ixGmin1=ixGmin1;s%ixGmin2=ixGmin2;s%ixGmin3=ixGmin3;s%ixGmax1=ixGmax1
    s%ixGmax2=ixGmax2;s%ixGmax3=ixGmax3;
     ixGsmin1 = ixGmin1-1; ixGsmax1 = ixGmax1; ixGsmin2 = ixGmin2-1
     ixGsmax2 = ixGmax2; ixGsmin3 = ixGmin3-1; ixGsmax3 = ixGmax3
    if(stagger_grid) then
      allocate(s%ws(ixGsmin1:ixGsmax1,ixGsmin2:ixGsmax2,ixGsmin3:ixGsmax3,&
         1:nws))
      s%ws=0.d0
      if(record_electric_field) then
        allocate(s%we(ixGsmin1:ixGsmax1,ixGsmin2:ixGsmax2,ixGsmin3:ixGsmax3,&
           1:nws))
        s%we=0.d0
      end if
      s%ixGsmin1=ixGsmin1;s%ixGsmin2=ixGsmin2;s%ixGsmin3=ixGsmin3
      s%ixGsmax1=ixGsmax1;s%ixGsmax2=ixGsmax2;s%ixGsmax3=ixGsmax3;
    end if
    if(alloc_once_for_ps) then
      ! allocate extra variables for ps state
      if(nw_extra>0) allocate(s%wextra(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1:nw_extra))
      ! All of this block's geometry lives in bgeo, one array per quantity for
      ! all blocks; point this block's views at its slice of them.
      call point_at_geometry(s, bgeo, igrid, ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
         ixGmax2,ixGmax3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,&
         ixGextmax2,ixGextmax3)
      ! allocate physical boundary flag
      allocate(s%is_physical_boundary(2*ndim))
      if(local_timestep) then
        allocate(s%dt(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3))
      endif
  
      if(B0field) then
        allocate(s%B0(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndir,&
           0:ndim))
        allocate(s%J0(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
           7-2*ndir:3))
      end if
      if(number_equi_vars > 0) then
        allocate(s%equi_vars(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
           1:number_equi_vars,0:ndim))
      endif
  
      ! allocate space for special values for each block state
      if(phys_trac) then
        ! special_values(1) Tcoff local
        ! special_values(2) Tmax local
        ! special_values(3:2+ndim) Bdir local
        allocate(s%special_values(ndim+2))
      end if
    else
      ! share common info from ps states to save memory
      if(nw_extra>0) s%wextra=>ps(igrid)%wextra
      s%x=>ps(igrid)%x
      s%dx=>ps(igrid)%dx
      s%ds=>ps(igrid)%ds
      s%dvolume=>ps(igrid)%dvolume
      s%surfaceC=>ps(igrid)%surfaceC
      s%is_physical_boundary=>ps(igrid)%is_physical_boundary
      if(B0field) then
        s%B0=>ps(igrid)%B0
        s%J0=>ps(igrid)%J0
      end if
      if(number_equi_vars > 0) then
        s%equi_vars=>ps(igrid)%equi_vars
      endif
      if(phys_trac) s%special_values=>ps(igrid)%special_values
   end if

  end subroutine alloc_state

  !> Point the metric components of state s at block igrid's slice of geo.
  !> Bounds remapping is what keeps the index ranges (which start below 1 for
  !> the face areas, and for the volumes when nghostcells is odd) intact: a
  !> plain pointer assignment would silently rebase them to 1.
  subroutine point_at_geometry(s, geo, igrid, ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
     ixGmax2,ixGmax3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
     ixGextmax3)
    use mod_physicaldata, only: state, geo_t
    use mod_global_parameters, only: ndim

    type(state), intent(inout) :: s
    type(geo_t), intent(inout), target :: geo
    integer, intent(in) :: igrid, ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixGextmin1,ixGextmin2,ixGextmin3,ixGextmax1,ixGextmax2,&
       ixGextmax3

    s%x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim) => &
       geo%x(:,:,:,:,igrid)

    ! A Cartesian build allocates nothing but the positions, so the metric
    ! views are left unassociated rather than remapped onto an unallocated
    ! array.  That is deliberate: any host reader that turns out to need one
    ! after all faults here instead of quietly reading nonsense.
#:if GEOM != 'Cartesian'
    s%ds(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
       1:ndim) => geo%ds(:,:,:,:,igrid)
    s%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
       ixGextmin3:ixGextmax3) => geo%dvolume(:,:,:,igrid)
    s%surfaceC(ixGmin1-1:ixGmax1,ixGmin2-1:ixGmax2,ixGmin3-1:ixGmax3,1:ndim) &
       => geo%surfaceC(:,:,:,:,igrid)
    s%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ixGextmin3:ixGextmax3,&
       1:ndim) => geo%dx(:,:,:,:,igrid)
#:endif

  end subroutine point_at_geometry

  !> allocate memory to one-level coarser physical state of igrid node
  subroutine alloc_state_coarse(igrid, s, ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
     ixGmax2,ixGmax3)
    use mod_global_parameters
    type(state) :: s
    integer, intent(in) :: igrid, ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3
    integer             :: ixGsmin1,ixGsmin2,ixGsmin3,ixGsmax1,ixGsmax2,&
       ixGsmax3
  
!    allocate(s%w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw))
    s%w => bgc(1)%w(:,:,:,:,igrid)


    s%igrid=igrid
    s%w=0.d0
    s%ixGmin1=ixGmin1;s%ixGmin2=ixGmin2;s%ixGmin3=ixGmin3;s%ixGmax1=ixGmax1
    s%ixGmax2=ixGmax2;s%ixGmax3=ixGmax3;
     ixGsmin1 = ixGmin1-1; ixGsmax1 = ixGmax1; ixGsmin2 = ixGmin2-1
     ixGsmax2 = ixGmax2; ixGsmin3 = ixGmin3-1; ixGsmax3 = ixGmax3
    if(stagger_grid) then
      allocate(s%ws(ixGsmin1:ixGsmax1,ixGsmin2:ixGsmax2,ixGsmin3:ixGsmax3,&
         1:nws))
      s%ws=0.d0
      s%ixGsmin1=ixGsmin1;s%ixGsmin2=ixGsmin2;s%ixGsmin3=ixGsmin3
      s%ixGsmax1=ixGsmax1;s%ixGsmax2=ixGsmax2;s%ixGsmax3=ixGsmax3;
    end if
    if(B0fieldAllocCoarse) then
      allocate(s%B0(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndir,&
         0:ndim))
    end if
    ! coarse representatives take their geometry from bgeoc, see alloc_state;
    ! they need no extra layer, so their "extended" range is just ixG
    call point_at_geometry(s, bgeoc, igrid, ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
       ixGmax2,ixGmax3, ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3)
    ! allocate physical boundary flag
    allocate(s%is_physical_boundary(2*ndim))

  end subroutine alloc_state_coarse

  subroutine dealloc_state(igrid, s,dealloc_x)
    use mod_global_parameters
    integer, intent(in) :: igrid
    type(state) :: s
    logical, intent(in) :: dealloc_x
  
    deallocate(s%w)
    if(stagger_grid) then
      deallocate(s%ws)
      if(record_electric_field) deallocate(s%we)
    end if
    if(dealloc_x) then
      if(nw_extra>0) deallocate(s%wextra)
      ! the geometry in bgeo outlives the block, only drop the views on it
      nullify(s%x,s%ds,s%dvolume,s%surfaceC,s%dx)
      deallocate(s%is_physical_boundary)
      if(B0field) then
        deallocate(s%B0)
        deallocate(s%J0)
      end if
      if(number_equi_vars > 0) then
        deallocate(s%equi_vars)
      end if
    else
      nullify(s%x,s%dx,s%ds,s%dvolume,s%surfaceC)
      nullify(s%is_physical_boundary)
      if(B0field) nullify(s%B0,s%J0)
      if(number_equi_vars > 0) then
        nullify(s%equi_vars)
      end if
      if(nw_extra>0) nullify(s%wextra)
    end if
  end subroutine dealloc_state
  
  subroutine dealloc_state_coarse(igrid, s)
    use mod_global_parameters
    integer, intent(in) :: igrid
    type(state) :: s
  
    deallocate(s%w)
    if(stagger_grid) then
      deallocate(s%ws)
    end if
    if(B0fieldAllocCoarse) then
      deallocate(s%B0)
    end if
    ! the geometry in bgeoc outlives the block, only drop the views on it
    nullify(s%x,s%ds,s%dvolume,s%surfaceC,s%dx)
    deallocate(s%is_physical_boundary)
  end subroutine dealloc_state_coarse
  
  subroutine dealloc_node(igrid)
    use mod_global_parameters
  
    integer, intent(in) :: igrid
  
    if (igrid==0) then
       call mpistop("trying to delete a non-existing grid in dealloc_node")
    end if
  
    call dealloc_state(igrid, ps(igrid),.true.)
    call dealloc_state_coarse(igrid, psc(igrid))
    if(.not.convert) then
      call dealloc_state(igrid, ps1(igrid),.false.)
      ! deallocate temporary solution space
      select case (t_integrator)
      case(ssprk3,ssprk4,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
        call dealloc_state(igrid, ps2(igrid),.false.)
      case(RK3_BT,rk4,ssprk5,IMEX_CB3a)
        call dealloc_state(igrid, ps2(igrid),.false.)
        call dealloc_state(igrid, ps3(igrid),.false.)
      case(IMEX_ARS3,IMEX_232)
        call dealloc_state(igrid, ps2(igrid),.false.)
        call dealloc_state(igrid, ps3(igrid),.false.)
        call dealloc_state(igrid, ps4(igrid),.false.)
      end select
    end if
  
  end subroutine dealloc_node

end module mod_amr_solution_node
