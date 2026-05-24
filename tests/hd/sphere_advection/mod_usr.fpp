module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision :: rhodens  = 10.0d0
  double precision :: rholight = 1.0d0
  !$acc declare copyin(rhodens, rholight)
  double precision :: p        = 10.0d0
  double precision :: rsphere  = 0.25d0
  double precision :: v1 = 1.0d0, v2 = 1.0d0, v3 = 1.0d0

contains

  subroutine usr_init()

    call set_coordinate_system("Cartesian_3D")

    usr_init_one_grid => initonegrid_usr
    usr_aux_output => special_output
    usr_add_aux_names => special_aux_names

    call phys_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
         ixGmin3:ixGmax3,1:nw)
    ! .. local ..
    double precision                :: r2
    integer                         :: ix1, ix2, ix3

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1

             r2 = x(ix1,ix2,ix3,1)**2 + x(ix1,ix2,ix3,2)**2 + x(ix1,ix2,ix3,3)**2

             if ( r2 < rsphere**2 ) then
                w(ix1,ix2,ix3,rho_) = rhodens
             else
                w(ix1,ix2,ix3,rho_) = rholight

             end if
             
             w(ix1,ix2,ix3,p_)     = p
             w(ix1,ix2,ix3,mom(1)) = v1
             w(ix1,ix2,ix3,mom(2)) = v2
             w(ix1,ix2,ix3,mom(3)) = v3
                
          end do
       end do
    end do
    
    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  subroutine usr_refine_grid(igrid,level,ixGmin1,ixGmin2,ixGmin3,&
    ixGmax1,ixGmax2,ixGmax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
    qt,w,x,refine,coarsen)

    use mod_global_parameters
    !$acc routine vector
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

    ! you must set consistent values for integers refine/coarsen:

    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement

    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen

    integer, intent(in)             :: igrid, level, ixGmin1,ixGmin2,&
        ixGmin3,ixGmax1,ixGmax2,ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,&
        ixmax2,ixmax3
    double precision, intent(in)    :: qt, x(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
    double precision, intent(in)    :: w(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw)
    integer, intent(inout) :: refine, coarsen
    ! .. local ..
    integer                         :: ix1, ix2, ix3
    logical                         :: has_sphere

    has_sphere = .false.
    !$acc loop collapse(3) reduction(.or.:has_sphere)
    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             if ( w(ix1, ix2, ix3, rho_) > (rhodens-rholight)/2.0d0 ) then
                has_sphere = .true.
             end if
          end do
       end do
    end do

    if (has_sphere) then
       coarsen = -1
       refine  = 1
    else 
       coarsen = 1
       refine  = -1
    end if

  end subroutine usr_refine_grid

    subroutine special_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
      double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,nw+nwauxio)
      double precision             :: normconv(0:nw+nwauxio)

      double precision                   :: threshold, error, numerator, denominator
      double precision, parameter        :: epsilon=1.0d-6
      integer                            :: iflag, idims1, idims2, level
      integer                            :: ix1, ix2, ix3, igrid

      igrid       = block%igrid
      level       = node(plevel_,igrid)
      threshold   = refine_threshold(level)

      do ix3 = ixMlo3, ixMhi3
         do ix2 = ixMlo2, ixMhi2
            do ix1 = ixMlo1, ixMhi1

               error = zero
               do iflag = 1, nw
                  if(w_refine_weight(iflag)==0.d0) cycle

                  numerator   = zero
                  denominator = zero
                  do idims1 = 1, ndim
                     do idims2 = 1, ndim

                        numerator = numerator + &
                             ( &
                             ( w(ix1+kr(1,idims2)+kr(1,idims1), &
                             ix2+kr(2,idims2)+kr(2,idims1), &
                             ix3+kr(3,idims2)+kr(3,idims1), iflag)    &
                             - w(ix1-kr(1,idims2)+kr(1,idims1), &
                             ix2-kr(2,idims2)+kr(2,idims1), &
                             ix3-kr(3,idims2)+kr(3,idims1), iflag) )  &
                             - &
                             ( w(ix1+kr(1,idims2)-kr(1,idims1), &
                             ix2+kr(2,idims2)-kr(2,idims1), &
                             ix3+kr(3,idims2)-kr(3,idims1), iflag)    &
                             - w(ix1-kr(1,idims2)-kr(1,idims1), &
                             ix2-kr(2,idims2)-kr(2,idims1), &
                             ix3-kr(3,idims2)-kr(3,idims1), iflag) )  &
                             )**2

                        denominator =  denominator + &
                             ( &
                             abs( &
                             w(ix1+2*kr(1,idims1), ix2+2*kr(2,idims1), ix3+2*kr(3,idims1), iflag) &
                             - w(ix1, ix2, ix3, iflag) &
                             ) &
                             + abs( &
                             w(ix1, ix2, ix3, iflag) &
                             - w(ix1-2*kr(1,idims1), ix2-2*kr(2,idims1), ix3-2*kr(3,idims1), iflag) &
                             ) &
                             + amr_wavefilter(level) * ( &
                             ( abs( w(ix1+kr(1,idims1)+kr(1,idims2), &
                             ix2+kr(2,idims1)+kr(2,idims2), &
                             ix3+kr(3,idims1)+kr(3,idims2), iflag) )   &
                             + abs( w(ix1-kr(1,idims1)+kr(1,idims2), &
                             ix2-kr(2,idims1)+kr(2,idims2), ix3-kr(3,idims1)+kr(3,idims2), iflag) ) ) &
                             + &
                             ( abs( w(ix1+kr(1,idims1)-kr(1,idims2), &
                             ix2+kr(2,idims1)-kr(2,idims2), &
                             ix3+kr(3,idims1)-kr(3,idims2), iflag) )   &
                             + abs( w(ix1-kr(1,idims1)-kr(1,idims2), &
                             ix2-kr(2,idims1)-kr(2,idims2), &
                             ix3-kr(3,idims1)-kr(3,idims2), iflag) ) ) &
                             ) &
                             )**2

                     end do
                  end do

                  error = error + w_refine_weight(iflag) * sqrt( numerator / max( denominator, epsilon ) )

               end do

               w(ix1, ix2, ix3 ,nw+1) = error

            end do
         end do
      end do

      w(:,:,:,nw+2) = mype

    end subroutine special_output

    !> Add names for the auxiliary variables
    subroutine special_aux_names(varnames)
      use mod_global_parameters
      character(len=*) :: varnames
      
      varnames = 'error cpu'
      
    end subroutine special_aux_names

  
end module mod_usr
