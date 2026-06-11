module mod_errest
  use mod_comm_lib, only: mpistop
  implicit none
  private

  public :: errest, forcedrefine_grid_io

contains

  !> Do all local error estimation which determines (de)refinement
  subroutine errest
    use mod_forest, only: refine, buffer, coarsen
    use mod_global_parameters

    integer :: igrid, iigrid

    if (igridstail==0) return

    select case (refine_criterion)
    case (1)
       ! all refinement solely based on user routine usr_refine_grid
    case (3)
       ! Error estimation is based on Lohner's scheme
       !$acc parallel loop gang private(igrid)
       !$omp target teams distribute private(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          call lohner_grid(igrid)
       end do
    case default
       call mpistop("Unknown error estimator")
    end select

    if ( refine_usr ) then
       !$acc parallel loop gang private(igrid)
       !$omp target teams distribute parallel do private(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          call forcedrefine_grid(igrid)
       end do
    end if

    !$acc update host(refine, coarsen)
    !$omp target update from(refine, coarsen)

  end subroutine errest

  subroutine lohner_grid(igrid)
    use mod_forest, only: coarsen, refine
    use mod_global_parameters
    !$acc routine vector
    !$omp declare target

    integer, intent(in) :: igrid

    integer                            :: iflag, idims1, idims2, level
    integer                            :: ix1, ix2, ix3
    double precision                   :: threshold, error, numerator, denominator
    logical                            :: refineflag, coarsenflag
    double precision, parameter        :: epsilon=1.0d-6

      level       = node(plevel_,igrid)
      threshold   = refine_threshold(level)

      refineflag  = .false.
      coarsenflag = .true.
      !$acc loop vector collapse(3) reduction(.or.:refineflag) reduction(.and.:coarsenflag)
      !$omp loop bind(parallel) collapse(3) reduction(.or.:refineflag) reduction(.and.:coarsenflag)
      do ix3 = ixMlo3, ixMhi3
         do ix2 = ixMlo2, ixMhi2
            do ix1 = ixMlo1, ixMhi1

               error = zero
               !$acc loop seq reduction(+:error)
               !ToDo: check if omp is needed here?
               do iflag = 1, nw
                  if(w_refine_weight(iflag)==0.d0) cycle

                  numerator   = zero
                  denominator = zero
                  !$acc loop seq reduction(+:numerator, denominator)
                  !$omp loop bind(thread) reduction(+:numerator, denominator)
                  do idims1 = 1, ndim
                     do idims2 = 1, ndim

                        numerator = numerator + &
                             ( &
                             ( bg(1)%w(ix1+kr(1,idims2)+kr(1,idims1), &
                             ix2+kr(2,idims2)+kr(2,idims1), &
                             ix3+kr(3,idims2)+kr(3,idims1), iflag, igrid)    &
                             - bg(1)%w(ix1-kr(1,idims2)+kr(1,idims1), &
                             ix2-kr(2,idims2)+kr(2,idims1), &
                             ix3-kr(3,idims2)+kr(3,idims1), iflag, igrid) )  &
                             - &
                             ( bg(1)%w(ix1+kr(1,idims2)-kr(1,idims1), &
                             ix2+kr(2,idims2)-kr(2,idims1), &
                             ix3+kr(3,idims2)-kr(3,idims1), iflag, igrid)    &
                             - bg(1)%w(ix1-kr(1,idims2)-kr(1,idims1), &
                             ix2-kr(2,idims2)-kr(2,idims1), &
                             ix3-kr(3,idims2)-kr(3,idims1), iflag, igrid) )  &
                             )**2

                        denominator = denominator + &
                             ( &
                             abs( &
                             bg(1)%w(ix1+2*kr(1,idims1), ix2+2*kr(2,idims1), ix3+2*kr(3,idims1), iflag, igrid) &
                             - bg(1)%w(ix1, ix2, ix3, iflag, igrid) &
                             ) &
                             + abs( &
                             bg(1)%w(ix1, ix2, ix3, iflag, igrid) &
                             - bg(1)%w(ix1-2*kr(1,idims1), ix2-2*kr(2,idims1), ix3-2*kr(3,idims1), iflag, igrid) &
                             ) &
                             + amr_wavefilter(level) * ( &
                             ( abs( bg(1)%w(ix1+kr(1,idims1)+kr(1,idims2), &
                             ix2+kr(2,idims1)+kr(2,idims2), &
                             ix3+kr(3,idims1)+kr(3,idims2), iflag, igrid) )   &
                             + abs( bg(1)%w(ix1-kr(1,idims1)+kr(1,idims2), &
                             ix2-kr(2,idims1)+kr(2,idims2), ix3-kr(3,idims1)+kr(3,idims2), iflag, igrid) ) ) &
                             + &
                             ( abs( bg(1)%w(ix1+kr(1,idims1)-kr(1,idims2), &
                             ix2+kr(2,idims1)-kr(2,idims2), &
                             ix3+kr(3,idims1)-kr(3,idims2), iflag, igrid) )   &
                             + abs( bg(1)%w(ix1-kr(1,idims1)-kr(1,idims2), &
                             ix2-kr(2,idims1)-kr(2,idims2), &
                             ix3-kr(3,idims1)-kr(3,idims2), iflag, igrid) ) ) &
                             ) &
                             )**2
                     end do
                  end do

                  error = error + w_refine_weight(iflag) * sqrt( numerator / max( denominator, epsilon ) )

               end do

               if (error > threshold) then
                  refineflag = .true.
               end if
               if (error > derefine_ratio(level) * threshold) then
                  coarsenflag = .false.
               end if

            end do
         end do
      end do

      if (refineflag .and. level < refine_max_level) refine(igrid,mype)=.true.
      if (coarsenflag .and. level > 1) coarsen(igrid,mype)=.true.

  end subroutine lohner_grid

  subroutine forcedrefine_grid(igrid)
    #:if defined('REFINE_USR')
    use mod_usr, only: usr_refine_grid
    #:endif
    use mod_forest, only: coarsen, refine, buffer
    use mod_global_parameters
    !$acc routine vector
    !$omp declare target

    integer, intent(in) :: igrid

    integer :: level
    integer :: my_refine, my_coarsen
    double precision :: qt

    level=node(plevel_,igrid)

    ! initialize to 0
    my_refine   = 0
    my_coarsen  = 0

    if (time_advance) then
       qt=global_time+dt
    else
       qt=global_time
    end if

#:if defined('REFINE_USR')
    call usr_refine_grid(igrid,level,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2, &
         ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,qt, &
         bg(1)%w(:,:,:,:, igrid), ps(igrid)%x, &
         my_refine,my_coarsen)
#:endif

    if (my_coarsen==1) then
       if (level>1) then
          refine(igrid,mype)=.false.
          coarsen(igrid,mype)=.true.
       else
          refine(igrid,mype)=.false.
          coarsen(igrid,mype)=.false.
       end if
    end if

    if (my_coarsen==-1)then
       coarsen(igrid,mype)=.false.
    end if

    if (my_refine==1) then
       if (level<refine_max_level) then
          refine(igrid,mype)=.true.
          coarsen(igrid,mype)=.false.
       else
          refine(igrid,mype)=.false.
          coarsen(igrid,mype)=.false.
       end if
    end if

    if (my_refine==-1) then
      refine(igrid,mype)=.false.
    end if

  end subroutine forcedrefine_grid

  subroutine forcedrefine_grid_io(igrid,w)
    use mod_forest, only: coarsen, refine
    use mod_global_parameters

    integer, intent(in)          :: igrid
    double precision, intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nw)

    integer                   :: level, my_levmin, my_levmax
    logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3) :: refineflag

    level=node(plevel_,igrid)

    if (level_io > 0) then
       my_levmin = level_io
       my_levmax = level_io
    else
       my_levmin = max(1,level_io_min)
       my_levmax = min(refine_max_level,level_io_max)
    end if

    if (level>my_levmax) then
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.true.
    elseif (level<my_levmin) then
      refine(igrid,mype)=.true.
      coarsen(igrid,mype)=.false.
    end if

    if (level==my_levmin .or. level==my_levmax) then
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.false.
    end if

    if(refine(igrid,mype).and.level>=refine_max_level)refine(igrid,&
       mype)=.false.
    if(coarsen(igrid,mype).and.level<=1)coarsen(igrid,mype)=.false.

  end subroutine forcedrefine_grid_io

end module mod_errest
