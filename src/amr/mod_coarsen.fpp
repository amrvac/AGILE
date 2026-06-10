module mod_coarsen
  use mod_global_parameters
  use mod_physics

  implicit none
  private

  public :: coarsen_grid

contains


  !> coarsen one grid to its coarser representative
  subroutine coarsen_grid(sFi,ixFiGmin1,ixFiGmin2,ixFiGmin3,ixFiGmax1,&
     ixFiGmax2,ixFiGmax3,ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,&
     sCo,ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,ixComin1,&
     ixComin2,ixComin3,ixComax1,ixComax2,ixComax3)

    type(state), intent(inout)      :: sFi, sCo
    integer, intent(in) :: ixFiGmin1,ixFiGmin2,ixFiGmin3,ixFiGmax1,ixFiGmax2,&
       ixFiGmax3, ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,&
        ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3, ixComin1,&
       ixComin2,ixComin3,ixComax1,ixComax2,ixComax3

    integer :: ixCo1,ixCo2,ixCo3, ixFi1,ixFi2,ixFi3, iw
    double precision :: CoFiratio
    double precision :: B_energy_change(ixCoGmin1:ixCoGmax1,&
       ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3)

    ! coarsen by 2 in every direction - conservatively

! AGILE: tbd
!    if(coarsenprimitive) call phys_to_primitive(ixFiGmin1,ixFiGmin2,ixFiGmin3,&
!       ixFiGmax1,ixFiGmax2,ixFiGmax3,ixFimin1,ixFimin2,ixFimin3,ixFimax1,&
!       ixFimax2,ixFimax3,wFi,sFi%x)

    if(slab_uniform) then
       CoFiratio=one/dble(2**ndim)
       !$acc parallel loop gang
       !$omp target teams loop
       do iw=1,nw
          !$acc loop collapse(3) vector
          !$omp loop collapse(3)
          do ixCo3 = ixComin3,ixComax3
             do ixCo2 = ixComin2,ixComax2
                do ixCo1 = ixComin1,ixComax1
                   ixFi3=2*(ixCo3-ixComin3)+ixFimin3+ixFiGmin3-1
                   ixFi2=2*(ixCo2-ixComin2)+ixFimin2+ixFiGmin2-1
                   ixFi1=2*(ixCo1-ixComin1)+ixFimin1+ixFiGmin1-1
                   sCo%w(ixCo1+ixCoGmin1-1,ixCo2+ixCoGmin2-1,ixCo3+ixCoGmin3-1,iw)=sum(sFi%w(ixFi1:ixFi1+1,ixFi2:ixFi2+1,&
                        ixFi3:ixFi3+1,iw))*CoFiratio
                end do
             end do
          end do
      end do
    else
       !$acc parallel loop gang
       !$omp target teams loop
       do iw=1,nw
          !$acc loop collapse(3) vector
          !$omp loop collapse(3)
          do ixCo3 = ixComin3,ixComax3
             do ixCo2 = ixComin2,ixComax2
                do ixCo1 = ixComin1,ixComax1
                   ixFi3=2*(ixCo3-ixComin3)+ixFimin3+ixFiGmin3-1
                   ixFi2=2*(ixCo2-ixComin2)+ixFimin2+ixFiGmin2-1
                   ixFi1=2*(ixCo1-ixComin1)+ixFimin1+ixFiGmin1-1
                   sCo%w(ixCo1+ixCoGmin1-1,ixCo2+ixCoGmin2-1,ixCo3+ixCoGmin3-1,iw)= sum(sFi%dvolume(ixFi1:ixFi1+1,&
                        ixFi2:ixFi2+1,ixFi3:ixFi3+1)*sFi%w(ixFi1:ixFi1+1,ixFi2:ixFi2+1,&
                        ixFi3:ixFi3+1,iw)) /sCo%dvolume(ixCo1,ixCo2,ixCo3)
                end do
             end do
          end do
       end do
    end if

! AGILE: tbd
    ! if(stagger_grid) then
    !   do iw=1,nws
    !     ! Start one layer before
    !     do ixCo3 = ixComin3-kr(3,iw),ixComax3
    !        ixFi3=2*(ixCo3-ixComin3+kr(3,iw))+ixFimin3-kr(3,iw)
    !     do ixCo2 = ixComin2-kr(2,iw),ixComax2
    !        ixFi2=2*(ixCo2-ixComin2+kr(2,iw))+ixFimin2-kr(2,iw)
    !     do ixCo1 = ixComin1-kr(1,iw),ixComax1
    !        ixFi1=2*(ixCo1-ixComin1+kr(1,iw))+ixFimin1-kr(1,iw)
    !        ! This if statement catches the axis where surface is zero:
    !        if (sCo%surfaceC(ixCo1,ixCo2,ixCo3,iw)>1.0d-9*sCo%dvolume(ixCo1,&
    !           ixCo2,ixCo3)) then !Normal case
    !          sCo%w(ixCo1,ixCo2,ixCo3,iw)=sum(sFi%surfaceC(ixFi1:ixFi1+1-kr(iw,&
    !             1),ixFi2:ixFi2+1-kr(iw,2),ixFi3:ixFi3+1-kr(iw,3),&
    !             iw)*sFi%ws(ixFi1:ixFi1+1-kr(iw,1),ixFi2:ixFi2+1-kr(iw,2),&
    !             ixFi3:ixFi3+1-kr(iw,3),iw)) /sCo%surfaceC(ixCo1,ixCo2,ixCo3,&
    !             iw)
    !        else ! On axis
    !          sCo%w(ixCo1,ixCo2,ixCo3,iw)=zero
    !        end if
    !     end do
    !     end do
    !     end do
    !   end do
    !   if(phys_total_energy.and. .not.coarsenprimitive) then
    !     B_energy_change(ixComin1:ixComax1,ixComin2:ixComax2,&
    !        ixComin3:ixComax3)=0.5d0*sum(wCo(ixComin1:ixComax1,&
    !        ixComin2:ixComax2,ixComin3:ixComax3,iw_mag(:))**2,dim=ndim+1)
    !   end if
    !   ! average to fill cell-centred values
    !   call phys_face_to_center(ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
    !      ixComax3,sCo)
    !   if(phys_total_energy.and. .not.coarsenprimitive) then
    !     wCo(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,&
    !        iw_e)=wCo(ixComin1:ixComax1,ixComin2:ixComax2,ixComin3:ixComax3,&
    !        iw_e)-B_energy_change(ixComin1:ixComax1,ixComin2:ixComax2,&
    !        ixComin3:ixComax3)+0.5d0*sum(wCo(ixComin1:ixComax1,&
    !        ixComin2:ixComax2,ixComin3:ixComax3,iw_mag(:))**2,dim=ndim+1)
    !   end if
    ! end if

! AGILE: tbd
!    if(coarsenprimitive) then
!      call phys_to_conserved(ixFiGmin1,ixFiGmin2,ixFiGmin3,ixFiGmax1,ixFiGmax2,&
!         ixFiGmax3,ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,wFi,&
!         sFi%x)
!      call phys_to_conserved(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,&
!         ixCoGmax3,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,wCo,&
!         sCo%x)
!    end if
  end subroutine coarsen_grid

end module mod_coarsen
