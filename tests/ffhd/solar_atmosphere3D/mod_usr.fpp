#:mute
#:include "../../../src/mod_gpu_directives.fpp"
#:endmute

module mod_usr
    use mod_amrvac
    use mod_physics

    implicit none

    double precision, allocatable :: pbc(:),rbc(:)
    ${GPU_DECLARE_CREATE('pbc,rbc')}$
    double precision :: usr_grav,trelax
    double precision :: heatunit,gzone,SRadius,dya
    double precision, allocatable :: pa(:),ra(:)
    integer, parameter :: jmax=5000
    double precision :: dh, Bh, dv, Bv, xv, dp1
    double precision :: ch, cv, dh1, dv1
    ${GPU_DECLARE_CREATE('dh,Bh,dv,Bv,xv,dp1,ch,cv,dh1,dv1')}$
    ${GPU_DECLARE_CREATE('usr_grav,SRadius,pbc,rbc')}$
  
  contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ dh,Bh,dv,Bv,xv,dp1 

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read
  
    subroutine usr_init()

      call set_coordinate_system("Cartesian_3D")
      call usr_params_read(par_files)
      ${GPU_UPDATE_DEVICE('dh,Bh,dv,Bv,xv,dp1')}$
  
      unit_length        = 1.d9 !< cm
      unit_temperature   = 1.d6 !< K
      unit_numberdensity = 1.d9 !< cm^-3
  
      usr_init_one_grid   => initonegrid_usr
      usr_set_parameters  => initglobaldata_usr
  
      call phys_activate()

    end subroutine usr_init
  
    subroutine initglobaldata_usr()

      heatunit=unit_pressure/unit_time !< 3.697693390805347E-003 erg*cm^-3/s
      gzone=0.2d0 !< thickness of a ghostzone below the bottom boundary
      dya=(2.d0*gzone+xprobmax3-xprobmin3)/dble(jmax) !< cells size of high-resolution 1D solar atmosphere

      usr_grav=-2.74d4*unit_length/unit_velocity**2 !< solar gravity
      ${GPU_UPDATE_DEVICE('usr_grav')}$

      SRadius=69.61d0 !< Solar radius
      ${GPU_UPDATE_DEVICE('SRadius')}$

      ch = (Bh*dh**3.d0)*0.5d0
      cv = (Bv*dv**3.d0)*0.5d0
      dh1 = dh + dp1
      dv1 = dv + dp1
      ${GPU_UPDATE_DEVICE('ch,cv,dh1,dv1')}$

      call inithdstatic

    end subroutine initglobaldata_usr
  
    subroutine inithdstatic
      !> hydrostatic vertical stratification of density, temperature, pressure
      integer :: i,j,na,ibc
      integer, parameter :: n_val=49
      double precision :: res
      double precision :: rpho,htra,Ttr,Fc,invT,kappa
      double precision :: h_val(n_val),t_val(n_val)
      double precision, allocatable :: ya(:),Ta(:),gg(:)
  
      rpho=0.71d15/unit_numberdensity !< bottom density at y=-gzone
      Fc=2.d5/heatunit/unit_length !< constant
      !> constant kappa but actually should be related to Ta, so not necessary to be the same with kappa in Thermal Conduction
      kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      !> VAL-C table, only temperature
      data h_val / 0, 50, 100, 150, 250, &
                   350, 450, 515, 555, 605, &
                   655, 705, 755, 855, 905, &
                   980, 1065, 1180, 1280, 1380, &
                   1515, 1605, 1785, 1925, 1990, &
                   2016, 2050, 2070, 2080, 2090, &
                   2104, 2107, 2109, 2113, 2115, &
                   2120, 2129, 2160, 2200, 2230, &
                   2255, 2263, 2267, 2271, 2274, &
                   2280, 2290, 2298, 2543 /
      data t_val / 6420, 5840, 5455, 5180, 4780, &
                   4465, 4220, 4170, 4230, 4420, &
                   4730, 5030, 5280, 5650, 5755, &
                   5925, 6040, 6150, 6220, 6280, &
                   6370, 6440, 6630, 6940, 7160, &
                   7360, 7660, 7940, 8180, 8440, &
                   9500, 10700, 12300, 18500, 21000, &
                   22500, 23000, 23500, 24000, 24200, &
                   24500, 25500, 28000, 32000, 37000, &
                   50000, 89100, 141000, 447000 /
      h_val(1:n_val)=h_val(1:n_val)*1.d5/unit_length
      t_val(1:n_val)=t_val(1:n_val)/unit_temperature
      htra=maxval(h_val)
      Ttr=maxval(t_val)
      allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))
      do j=1,jmax
        ya(j)=(dble(j)-0.5d0)*dya-gzone
        if(ya(j)>=htra) then
          Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
        else
          do i=1,n_val
            if(ya(j)<h_val(i+1)) then
              Ta(j)=t_val(i)+(ya(j)-h_val(i))*(t_val(i+1)-t_val(i))/(h_val(i+1)-h_val(i))
              exit
            endif
          enddo
        endif
        !> keep gg the same with the settings in gravity
        gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
      enddo
      ra(1)=rpho
      pa(1)=rpho*Ta(1)
      invT=0.d0
      do j=2,jmax
        invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
        pa(j)=pa(1)*dexp(invT*dya)
        ra(j)=pa(j)/Ta(j)
      end do
      allocate(rbc(nghostcells))
      allocate(pbc(nghostcells))
      do ibc=nghostcells,1,-1
        na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dya+0.5d0)
        res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
        rbc(ibc)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        pbc(ibc)=pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
      end do
      deallocate(ya,gg,Ta)
      ${GPU_UPDATE_DEVICE('rbc,pbc')}$
    end subroutine inithdstatic

    subroutine initonegrid_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
        integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
        double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3,1:ndim)
        double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3,1:nw)
    
        double precision :: res
        integer :: ix1, ix2, ix3 ,na

        double precision :: xlocal(3)
        double precision :: bx, by, bz
        double precision :: xh, yh, zh, yv, zv

        xh = 0.d0
        yh = 0.d0
        zh = -dh1
        yv = 0.d0
        zv = -dv1

  
        do ix3=ixOmin3,ixOmax3
        do ix2=ixOmin2,ixOmax2
        do ix1=ixOmin1,ixOmax1

            na=floor((x(ix1,ix2,ix3,3)-xprobmin3+gzone)/dya+0.5d0)
            res=x(ix1,ix2,ix3,3)-xprobmin3+gzone-(dble(na)-0.5d0)*dya
            w(ix1,ix2,ix3,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
            w(ix1,ix2,ix3,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))

        end do
        end do
        end do

        do ix3=ixImin3,ixImax3
        do ix2=ixImin2,ixImax2
        do ix1=ixImin1,ixImax1

          xlocal(1) = x(ix1,ix2,ix3,1)
          xlocal(2) = x(ix1,ix2,ix3,2)
          xlocal(3) = x(ix1,ix2,ix3,3)
          bx = -3.d0*ch*((xlocal(3)-zh)**2+(xlocal(2)-yh)**2) / ((xlocal(1)-xh)**2+(xlocal(2)-yh)**2+(xlocal(3)-zh)**2)**(5./2.) + &
            3.d0*cv*(xlocal(1)-xv)*(xlocal(3)-zv) / ((xlocal(1)-xv)**2+(xlocal(2)-yv)**2+(xlocal(3)-zv)**2)**(5./2.) + &
            2*ch / ((xlocal(1)-xh)**2+(xlocal(2)-yh)**2+(xlocal(3)-zh)**2)**(3./2.)

          by = 3.d0*ch*(xlocal(1)-xh)*(xlocal(2)-yh) / ((xlocal(1)-xh)**2+(xlocal(2)-yh)**2+(xlocal(3)-zh)**2)**(5./2.) + &
            3.d0*cv*(xlocal(2)-yv)*(xlocal(3)-zv) / ((xlocal(1)-xv)**2+(xlocal(2)-yv)**2+(xlocal(3)-zv)**2)**(5./2.)

          bz = 3.d0*ch*(xlocal(1)-xh)*(xlocal(3)-zh) / ((xlocal(1)-xh)**2+(xlocal(2)-yh)**2+(xlocal(3)-zh)**2)**(5./2.) + & 
            (-3.d0)*cv*((xlocal(2)-yv)**2+(xlocal(1)-xv)**2) / ((xlocal(1)-xv)**2+(xlocal(2)-yv)**2+(xlocal(3)-zv)**2)**(5./2.) + & 
            2.d0*cv/((xlocal(1)-xv)**2+(xlocal(2)-yv)**2+(xlocal(3)-zv)**2)**(3./2.)

          w(ix1,ix2,ix3,iw_b1) = bx/(bx**2+by**2+bz**2)**(1./2.)
          w(ix1,ix2,ix3,iw_b2) = by/(bx**2+by**2+bz**2)**(1./2.)
          w(ix1,ix2,ix3,iw_b3) = bz/(bx**2+by**2+bz**2)**(1./2.)
        end do
        end do
        end do

        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(:))=zero
        #:if defined('HYPERTC')
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,q_)=zero
        #:endif
        call phys_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    end subroutine initonegrid_usr

    subroutine specialbound_usr(qt, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
        ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB, w, x)
        ${GPU_ROUTINE_VECTOR()}$

        integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
        ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB
        double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
        ixImin3:ixImax3, 1:ndim)
        double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        ixImin3:ixImax3, 1:nw)

        integer                         :: ix1,ix2,ix3
        double precision                :: inv_gamma_m1
        double precision                :: pth, invT, gravity, tmp

        inv_gamma_m1 = 1.0d0/(phys_gamma - 1.0d0)

        select case(iB)
        case(5)

        ${GPU_LOOP_VECTOR("collapse(3)")}$
        do ix3=ixOmin3,ixOmax3
            do ix2=ixOmin2,ixOmax2
            do ix1=ixOmin1,ixOmax1
            w(ix1,ix2,ix3,iw_mom(1)) = &
                -w(ix1,ix2,2*ixOmax3+1-ix3,iw_mom(1))&
                /w(ix1,ix2,2*ixOmax3+1-ix3,rho_)
            w(ix1,ix2,ix3,rho_)=rbc(ix3)
            w(ix1,ix2,ix3,e_)=pbc(ix3)
            #:if defined('HYPERTC')
            w(ix1,ix2,ix3,q_)=zero
            #:endif

            ! inline to_conservative() for performance:
            w(ix1,ix2,ix3, e_)=w(ix1,ix2,ix3, e_)*inv_gamma_m1+0.5_dp*w(ix1,ix2,ix3, iw_rho)*&
                w(ix1,ix2,ix3, iw_mom(1))**2
            w(ix1,ix2,ix3, iw_mom(1)) = w(ix1,ix2,ix3, iw_rho)*w(ix1,ix2,ix3, iw_mom(1))

            end do
            end do
        enddo

        ! case(6)

        ! ! to include ix3 into loop, perhaps precompute?
        ! ${GPU_LOOP_VECTOR("collapse(2)")}$
        ! do ix2=ixOmin2,ixOmax2
        ! do ix1=ixOmin1,ixOmax1

        !   ${GPU_ROUTINE_SEQ()}$
        !   tmp = 0.d0
        !   do ix3=ixOmin3,ixOmax3
        !     gravity = 0.5d0*(gravity_field(w(ix1,ix2,ix3,:), x(ix1,ix2,ix3,:), 3) + &
        !                 gravity_field(w(ix1,ix2,ix3-1,:), x(ix1,ix2,ix3-1,:), 3))
        !     tmp = tmp + gravity*invT
        !     tmp_array(ix3) = tmp
        !   end do

        ! end do 
        ! end do
        ! ${GPU_LOOP_VECTOR("collapse(2)")}$
        !   do ix2=ixOmin2,ixOmax2
        !     do ix1=ixOmin1,ixOmax1

        !     ! inline get pthermal to improve performance
        !     pth = (phys_gamma-1.0_dp)*(w(ix1,ix2,ixOmin3-1, e_)- &
        !       0.5_dp*w(ix1,ix2,ixOmin3-1, mom(1))**2/w(ix1,ix2,ixOmin3-1, rho_))
        !     invT     = w(ix1,ix2,ixOmin3-1,rho_)/pth

        !     tmp = 0.d0
        !     do ix3=ixOmin3,ixOmax3
        !       gravity = 0.5d0*(gravity_field(w(ix1,ix2,ix3,:), x(ix1,ix2,ix3,:), 3) + &
        !                 gravity_field(w(ix1,ix2,ix3-1,:), x(ix1,ix2,ix3-1,:), 3))
        !       tmp = tmp + gravity*invT
        !       ! repeated dxlevel(3) calculation to be optimized
        !       w(ix1,ix2,ix3,p_)=pth*dexp(tmp*(x(ix1,ix2,ix3,3)-x(ix1,ix2,ix3-1,3)))
        !       w(ix1,ix2,ix3,rho_) = w(ix1,ix2,ix3,p_)*invT
        !       #:if defined('HYPERTC')
        !       w(ix1,ix2,ix3,q_)=zero
        !       #:endif
        !       w(ix1,ix2,ix3,iw_mom(1)) = -w(ix1,ix2,2*ixOmin3-1-ix3,iw_mom(1)) / &
        !                                   w(ix1,ix2,2*ixOmin3-1-ix3,rho_)  

        !       ! inline to_conservative() for performance:
        !       w(ix1,ix2,ix3, e_)=w(ix1,ix2,ix3, e_)*inv_gamma_m1 +0.5_dp*w(ix1,ix2,ix3, rho_)*&
        !           w(ix1,ix2,ix3, iw_mom(1))**2
        !       w(ix1,ix2,ix3, iw_mom(1)) = w(ix1,ix2,ix3, iw_rho)*w(ix1,ix2,ix3, iw_mom(1))
        !       end do
        !       end do
        !   end do
      end select

   end subroutine specialbound_usr

   pure real(dp) function gravity_field(wCT, x, idim) result(field)
     ${GPU_ROUTINE_SEQ()}$
     real(dp), intent(in)    :: wCT(nw_phys)
     real(dp), intent(in)    :: x(1:ndim)
     integer, value, intent(in)     :: idim

     if (idim == 1) field =  0.0_dp
     if (idim == 2) field =  0.0_dp
     if (idim == 3) then
        field = -usr_grav*(SRadius/(SRadius+x(3)))**2
     end if

   end function gravity_field

  end module mod_usr
