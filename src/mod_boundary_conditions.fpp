#:mute
#:include "mod_gpu_directives.fpp"
#:endmute
module mod_boundary_conditions

  implicit none
  private

  public :: bc_phys
  public :: getintbc

contains

  !> fill ghost cells at a physical boundary
  !> x is this block's slice of bgeo/bgeoc rather than s%x, so that device
  !> code never has to dereference a per-block pointer to reach the positions.
  subroutine bc_phys(iside,idims,time,qdt,s,x,ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
       ixGmax2,ixGmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3)
#:if defined('SPECIALBOUNDARY')    
    use mod_usr, only: specialbound_usr
#:endif
    use mod_bc_data, only: bc_data_set
    use mod_global_parameters
    ${GPU_ROUTINE_VECTOR()}$

    integer, intent(in) :: iside, idims, ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
       ixGmax2,ixGmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3
    double precision, intent(in) :: time,qdt
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    type(state), intent(inout) :: s

    integer :: idir, is
    integer :: ixOsmin1,ixOsmin2,ixOsmin3,ixOsmax1,ixOsmax2,ixOsmax3,hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,jxOmin1,jxOmin2,jxOmin3,jxOmax1,&
       jxOmax2,jxOmax3
    integer :: iw, iB, ix1,ix2,ix3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3, ixMmin1,ixMmin2,ixMmin3,ixMmax1,ixMmax2,ixMmax3, nghostcellsi,&
       iib1,iib2,iib3
    logical  :: isphysbound

    associate(w=>s%w,ws=>s%ws)
    select case (idims)
    case (1)
       if (iside==2) then
          ! maximal boundary
          iB=2*1
          ixOmin1=ixBmax1+1-nghostcells;ixOmin2=ixBmin2;ixOmin3=ixBmin3;
          ixOmax1=ixBmax1;ixOmax2=ixBmax2;ixOmax3=ixBmax3;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix1=ixOmin1,ixOmax1
                   w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) = w(ixOmin1-1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+1)then
                  do ix1=ixOmin1,ixOmax1
                      w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                         iw) = max(w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                         iw),zero)
                  end do
                else
                  do ix1=ixOmin1,ixOmax1
                      w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) = w(ixOmin1-1,&
                         ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)
                  end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       else
          ! minimal boundary
          iB=2*1-1
          ixOmin1=ixBmin1;ixOmin2=ixBmin2;ixOmin3=ixBmin3;
          ixOmax1=ixBmin1-1+nghostcells;ixOmax2=ixBmax2;ixOmax3=ixBmax3;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix1=ixOmin1,ixOmax1
                   w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) = w(ixOmax1+1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+1)then
                   do ix1=ixOmin1,ixOmax1
                     w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                        iw) = min(w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                        iw),zero)
                   end do
                else
                   do ix1=ixOmin1,ixOmax1
                     w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) = w(ixOmax1+1,&
                        ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)
                   end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       end if
    case (2)
       if (iside==2) then
          ! maximal boundary
          iB=2*2
          ixOmin1=ixBmin1;ixOmin2=ixBmax2+1-nghostcells;ixOmin3=ixBmin3;
          ixOmax1=ixBmax1;ixOmax2=ixBmax2;ixOmax3=ixBmax3;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix2=ixOmin2,ixOmax2
                   w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                      iw) = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,&
                   ixOmin3:ixOmax3,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,&
                   ixOmin3:ixOmax3,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+2)then
                  do ix2=ixOmin2,ixOmax2
                      w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                         iw) = max(w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,&
                         iw),zero)
                  end do
                else
                  do ix2=ixOmin2,ixOmax2
                      w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                         iw) = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,iw)
                  end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       else
          ! minimal boundary
          iB=2*2-1
          ixOmin1=ixBmin1;ixOmin2=ixBmin2;ixOmin3=ixBmin3;
          ixOmax1=ixBmax1;ixOmax2=ixBmin2-1+nghostcells;ixOmax3=ixBmax3;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix2=ixOmin2,ixOmax2
                   w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                      iw) = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,&
                   ixOmin3:ixOmax3,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,&
                   ixOmin3:ixOmax3,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+2)then
                   do ix2=ixOmin2,ixOmax2
                     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                        iw) = min(w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,&
                        iw),zero)
                   end do
                else
                   do ix2=ixOmin2,ixOmax2
                     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                        iw) = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,iw)
                   end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       end if 
    case (3)
       if (iside==2) then
          ! maximal boundary
          iB=2*3
          ixOmin1=ixBmin1;ixOmin2=ixBmin2;ixOmin3=ixBmax3+1-nghostcells;
          ixOmax1=ixBmax1;ixOmax2=ixBmax2;ixOmax3=ixBmax3;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix3=ixOmin3,ixOmax3
                   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                      iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3-1:ixOmin3-nghostcells:-1,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3-1:ixOmin3-nghostcells:-1,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+3)then
                  do ix3=ixOmin3,ixOmax3
                      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                         iw) = max(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1,&
                         iw),zero)
                  end do
                else
                  do ix3=ixOmin3,ixOmax3
                      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                         iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1,iw)
                  end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       else
          ! minimal boundary
          iB=2*3-1
          ixOmin1=ixBmin1;ixOmin2=ixBmin2;ixOmin3=ixBmin3;
          ixOmax1=ixBmax1;ixOmax2=ixBmax2;ixOmax3=ixBmin3-1+nghostcells;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix3=ixOmin3,ixOmax3
                   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                      iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+1,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmax3+nghostcells:ixOmax3+1:-1,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmax3+nghostcells:ixOmax3+1:-1,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+3)then
                   do ix3=ixOmin3,ixOmax3
                     w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                        iw) = min(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+1,&
                        iw),zero)
                   end do
                else
                   do ix3=ixOmin3,ixOmax3
                     w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                        iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+1,iw)
                   end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       end if 
    end select

    ! Fill external data after the ordinary boundary conditions.
    if (any(typeboundary(1:nwflux+nwaux,iB)==bc_data .or. &
            typeboundary(1:nwflux+nwaux,iB)==bc_icarus)) then
       call bc_data_set(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2, &
            ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB,w,x)
    end if

#:if defined('SPECIALBOUNDARY')    
    ! do user defined special boundary conditions
    if (any(typeboundary(1:nwflux+nwaux,iB)==bc_special .or. &
            typeboundary(1:nwflux+nwaux,iB)==bc_icarus)) then
       call specialbound_usr(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
            ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB,w,x)

    end if
#:endif
    
  end associate
  end subroutine bc_phys

  !> Fill internal boundary values: on the device, apply the user's by-name
  !> per-cell hook usr_internal_bc (see mod_usr_methods) to every interior
  !> cell of every block of time level bgstep. Only variables
  !> nwhead:nwtail are written back, honouring getbc's variable range. Ghost
  !> cells are not touched: getbc's exchange, which follows, fills them from
  !> the overwritten interior.
  subroutine getintbc(time,bgstep,nwhead,nwtail)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
#:if defined('INTERNAL_BC')
    use mod_usr, only: usr_internal_bc
    use mod_physics_vars, only: nw_phys
#:endif

    double precision, intent(in)   :: time
    integer, intent(in)            :: bgstep, nwhead, nwtail

#:if defined('INTERNAL_BC')
    integer :: iigrid, igrid, level, ix1, ix2, ix3
    double precision :: wloc(1:nw_phys), xloc(1:ndim)

    ${GPU_PARALLEL_LOOP_GANG("private(igrid, level)")}$
    do iigrid = 1, igridstail; igrid = igrids(iigrid)
       level = node(plevel_,igrid)
       ${GPU_LOOP_VECTOR("collapse(3) private(wloc, xloc)")}$
       do ix3 = ixMlo3, ixMhi3
          do ix2 = ixMlo2, ixMhi2
             do ix1 = ixMlo1, ixMhi1
                ! fixed-size private copies, as in forcedrefine_grid: a
                ! strided section would be a device-heap temporary per cell
                wloc = bg(bgstep)%w(ix1,ix2,ix3,1:nw_phys,igrid)
                xloc = bgeo%x(ix1,ix2,ix3,1:ndim,igrid)
                call usr_internal_bc(level, time, xloc, wloc)
                bg(bgstep)%w(ix1,ix2,ix3,nwhead:min(nwtail,nw_phys),igrid) = &
                     wloc(nwhead:min(nwtail,nw_phys))
             end do
          end do
       end do
    end do
#:else
    call mpistop('internalboundary=T needs a build with internalboundary=T '//&
         'in &boundlist of agile.par (the INTERNAL_BC compile-time flag)')
#:endif

  end subroutine getintbc

end module mod_boundary_conditions
