!> This module contains variables that describe the connectivity of the mesh and
!> also data structures for connectivity-related communication.
module mod_connectivity
   implicit none

   integer, parameter :: neighbor_boundary = 1
   integer, parameter :: neighbor_coarse = 2
   integer, parameter :: neighbor_sibling = 3
   integer, parameter :: neighbor_fine = 4

   integer, dimension(:,:,:,:,:), allocatable :: neighbor
   integer, dimension(:,:,:,:,:), allocatable :: neighbor_child
   integer, dimension(:,:,:,:), allocatable :: neighbor_type
   logical, dimension(:,:,:,:), allocatable :: neighbor_active
   integer, dimension(:,:,:,:), allocatable :: neighbor_pole
   !$acc declare create(neighbor, neighbor_type, neighbor_pole, neighbor_child)

   ! grid number array per processor
   integer, dimension(:), allocatable :: igrids
   integer, dimension(:), allocatable :: igrids_active
   integer, dimension(:), allocatable :: igrids_passive
   !$acc declare create(igrids, igrids_active, igrids_passive)

   ! phys boundary indices
   integer, dimension(:,:), allocatable :: idphyb
   !$acc declare create(idphyb)

   ! number of grids on current processor
   integer :: igridstail
   integer :: igridstail_active
   integer :: igridstail_passive
   !$acc declare create(igridstail, igridstail_active, igridstail_passive)

   integer, dimension(3) :: nrecv_fc, nsend_fc
   ! cc for corner coarse
   integer, dimension(3) :: nrecv_cc, nsend_cc

   ! srl neighbor info
   type nbinfo_srl_t
      integer                            :: nigrids=0
      integer, pointer, dimension(:) :: igrid, iencode, isize, ibuf_start
   end type nbinfo_srl_t

   ! c and f neighbor info
   type nbinfo_cf_t
      integer                            :: nigrids=0
      integer, allocatable, dimension(:) :: igrid, inc1, inc2, inc3, isize, ibuf_start, i1, i2, i3
    contains
      procedure, non_overridable         :: init => init_cf
   end type nbinfo_cf_t

   ! buffer
   type nbinfo_buffer_t
      double precision, allocatable, dimension(:) :: buffer
      integer                                     :: size = -1
   end type nbinfo_buffer_t

   type nbinfo_buffer_i_t
      integer, allocatable, dimension(:) :: buffer
      integer                            :: size = -1
   end type nbinfo_buffer_i_t

   type course_nb_t
      type(nbinfo_cf_t)                    :: info            ! info for each nbproc
      type(nbinfo_buffer_t)                :: send            ! double precision send data. One for each nb proc
      type(nbinfo_buffer_t)                :: rcv             ! double precision receive data
      type(nbinfo_buffer_i_t)              :: info_send       ! info package send
      type(nbinfo_buffer_i_t)              :: info_rcv        ! info package receive
   contains
      procedure, non_overridable           :: init => init_course_nb
   end type course_nb_t
   type srl_nb_t
      type(nbinfo_srl_t)                   :: info            ! list of the ipelist for each nbproc
      type(nbinfo_buffer_t)                :: send            ! double precision send data. One for each nb proc
      type(nbinfo_buffer_t)                :: rcv             ! double precision receive data
      type(nbinfo_buffer_i_t)              :: info_send       ! info package send
      type(nbinfo_buffer_i_t)              :: info_rcv        ! info package receive
   contains
      procedure, non_overridable           :: init => init_srl_nb
    end type srl_nb_t

   ! neighbor cpu info structure
   type nbprocs_info_t
      ! SRL (neighbor is at same level)
      integer                              :: nbprocs_srl=0        ! number of neighboring processes at srl
      integer                              :: imaxigrids_srl=0     ! max number of igrids over all neighbor processors (for loop collasing)
      integer, allocatable                 :: nbprocs_srl_list(:)  ! list of neighboring ipe at srl
      integer, allocatable                 :: ipe_to_inbpe_srl(:)  ! inverse to nbprocs_srl_list
      ! F (neighbor is finer)
      integer                              :: nbprocs_f=0          ! number of neighboring processes with finer grids
      integer                              :: imaxigrids_f=0       ! max number of igrids over all neighbor processors (for loop collasing)
      integer, allocatable                 :: nbprocs_f_list(:)    ! list of neighboring ipe at f
      integer, allocatable                 :: ipe_to_inbpe_f(:)    ! inverse to nbprocs_f_list
      type(nbinfo_cf_t), allocatable       :: f(:)                 ! list of the ipelist for each nbproc
      type(nbinfo_buffer_t), allocatable   :: f_send(:)            ! double precision send data. One for each nb proc
      type(nbinfo_buffer_t), allocatable   :: f_rcv(:)             ! double precision receive data
      type(nbinfo_buffer_i_t), allocatable :: f_info_send(:)       ! info package send
      type(nbinfo_buffer_i_t), allocatable :: f_info_rcv(:)        ! info package receive
      ! C (neighbor is coarser)
      integer                              :: nbprocs_c=0          ! number of neighboring processes with coarser grids
      integer                              :: imaxigrids_c=0       ! max number of igrids over all neighbor processors (for loop collasing)
      integer, pointer                     :: nbprocs_c_list(:) => null()    ! list of neighboring ipe at f
      integer, allocatable                 :: ipe_to_inbpe_c(:)    ! inverse to nbprocs_f_list
      type(course_nb_t), allocatable       :: course_nb(:)         ! Info about course neighbours
      type(srl_nb_t), allocatable          :: srl_nb(:)            ! Info about same refinment level neighbours
      integer                              :: max_size = -1        ! maximum buffer size (initialized in init)
      integer                              :: max_nbprocs = 2     ! maximum nr of neighbor procs
      integer                              :: max_igrids = 2    ! maximum nr of igrids per neighbor proc
    contains
      procedure, non_overridable :: init, reset
      procedure, non_overridable :: add_igrid_to_srl, add_igrid_to_c, add_igrid_to_f
      procedure, non_overridable :: add_to_srl, add_to_f, add_to_c
      procedure, non_overridable :: alloc_buffers_srl, alloc_buffers_f, alloc_buffers_c
      procedure, non_overridable :: iencode, idecode
   end type nbprocs_info_t

   type(nbprocs_info_t) :: nbprocs_info
   !$acc declare create(nbprocs_info)

   public  :: nbprocs_info_t, nbprocs_info, nbprocs_update_device
   private :: reset, init, add_to_srl, add_to_f, add_to_c
   private :: init_course_nb
   private :: add_igrid_to_srl, alloc_buffers_srl, iencode, idecode
   private :: add_igrid_to_c, add_igrid_to_f


 contains
   subroutine init_srl_info(self, nigrids)
      class(nbinfo_srl_t) :: self
      integer :: nigrids
      integer, pointer, dimension(:) :: igrid, iencode, isize, ibuf_start
      logical :: realloc
      realloc =  allocated(self%igrid)

      if (realloc) then
        allocate(igrid(nigrids), iencode(nigrids), ibuf_start(nigrids), isize(nigrids))

        igrid(:size(self%igrid)) = self%igrid
        iencode(:size(self%iencode)) = self%iencode
        ibuf_start(:size(self%ibuf_start)) = self%ibuf_start
        isize(:size(self%isize)) = self%isize

        deallocate(self%igrid, self%iencode, self%ibuf_start, self%isize)

        self%igrid => igrid
        self%iencode => iencode
        self%ibuf_start => ibuf_start
        self%isize => isize
      else
        allocate(self%igrid(nigrids), self%iencode(nigrids), self%ibuf_start(nigrids), self%isize(nigrids))
      end if
   end subroutine

   subroutine init_srl_nb(self, nigrids, bufsize)
      class(srl_nb_t) :: self
      integer :: nigrids, bufsize
      call init_srl_info(self%info, nigrids)
      allocate(self%info_send%buffer( 3 * bufsize ))
      allocate(self%info_rcv%buffer( 3 *bufsize ))
      allocate(self%send%buffer(bufsize))
      allocate(self%rcv%buffer(bufsize))
   end subroutine

   subroutine init_course_nb(self, nigrids, bufsize)
      class(course_nb_t) :: self
      integer :: nigrids, bufsize
      call self%info%init(nigrids)
      allocate(self%info_send%buffer( 5 * bufsize ))
      allocate(self%info_rcv%buffer( 5 * bufsize ))
      allocate(self%send%buffer(bufsize))
      allocate(self%rcv%buffer(bufsize))
   end subroutine

   subroutine init_cf(self, nigrids)
     class(nbinfo_cf_t)   :: self
     integer, intent(in)  :: nigrids

     allocate(self%igrid(nigrids), self%inc1(nigrids), self%inc2(nigrids), &
          self%inc3(nigrids), self%ibuf_start(nigrids), self%isize(nigrids), &
          self%i1(nigrids), self%i2(nigrids), self%i3(nigrids))

   end subroutine init_cf

   subroutine reset(self)
     class(nbprocs_info_t) :: self
     ! .. local ..
     integer               :: i

     do i=1, self%nbprocs_srl
        self%srl_nb(i)%info%nigrids=0
     end do

     do i=1, self%nbprocs_c
        self%course_nb(i)%info%nigrids=0
     end do

     do i=1, self%nbprocs_f
        self%f(i)%nigrids=0
     end do

     self%nbprocs_srl         = 0
     self%imaxigrids_srl      = 0
     self%ipe_to_inbpe_srl(:) = -1

     self%nbprocs_c           = 0
     self%imaxigrids_c        = 0
     self%ipe_to_inbpe_c(:)   = -1

     self%nbprocs_f           = 0
     self%imaxigrids_f        = 0
     self%ipe_to_inbpe_f(:)   = -1

   end subroutine reset

   subroutine init(self, npe, max_size)
     class(nbprocs_info_t)  :: self
     integer, intent(in)   :: npe, max_size
     ! .. local ..
     integer               :: i

     self%max_size    = max_size

     allocate(self%nbprocs_srl_list(self%max_nbprocs), &
          self%srl_nb(self%max_nbprocs))

     allocate(self%nbprocs_f_list(self%max_nbprocs), &
          self%f(self%max_nbprocs))

     allocate(self%nbprocs_c_list(self%max_nbprocs), &
          self%course_nb(self%max_nbprocs))

     allocate(self%ipe_to_inbpe_srl(0:npe-1))
     self%ipe_to_inbpe_srl(:) = -1

     allocate(self%ipe_to_inbpe_f(0:npe-1))
     self%ipe_to_inbpe_f(:) = -1

     allocate(self%ipe_to_inbpe_c(0:npe-1))
     self%ipe_to_inbpe_c(:) = -1

     do i = 1, self%max_nbprocs
        call self%srl_nb(i)%init(self%max_igrids, self%max_size)
        call self%f(i)%init(self%max_igrids)
        call self%course_nb(i)%init(self%max_igrids, self%max_size)
     end do

     ! pre-allocate buffers
     allocate( &
          self%f_send(1:self%max_nbprocs), &
          self%f_rcv(1:self%max_nbprocs), &
          self%f_info_send(1:self%max_nbprocs), &
          self%f_info_rcv(1:self%max_nbprocs) &
          )

     do i = 1, self%max_nbprocs
        allocate( &
             self%f_info_send(i)%buffer( 5 * self%max_igrids ), &
             self%f_info_rcv(i)%buffer( 5 * self%max_igrids ), &
             self%f_send(i)%buffer(self%max_size), &
             self%f_rcv(i)%buffer(self%max_size) &
        )
     end do

   end subroutine init

   subroutine add_to_srl(self, ipe, igrid, i1, i2, i3)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: ipe, igrid, i1, i2, i3

     ! add the procesor id if not already present:
     if (self%ipe_to_inbpe_srl(ipe) == -1) then

        ! enlarge counter
        self%nbprocs_srl = self%nbprocs_srl + 1

        if (self%nbprocs_srl > self%max_nbprocs) then
           print *, 'add_ipe_to_srl_list: nbprocs_srl exceeding max_nbprocs', self%nbprocs_srl, self%max_nbprocs
           error stop

        end if

        ! add the process
        self%nbprocs_srl_list(self%nbprocs_srl) = ipe
        ! add to inverse list
        self%ipe_to_inbpe_srl(ipe) = self%nbprocs_srl

     end if

     call self%add_igrid_to_srl(ipe, igrid, i1, i2, i3)

   end subroutine add_to_srl

   subroutine add_to_c(self, ipe, igrid, i1, i2, i3, inc1, inc2, inc3)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: ipe, igrid, i1, i2, i3, inc1, inc2, inc3

     ! add the procesor id if not already present:
     if (self%ipe_to_inbpe_c(ipe) == -1) then

        ! enlarge counter
        self%nbprocs_c = self%nbprocs_c + 1

        if (self%nbprocs_c > self%max_nbprocs) then

           !allocate(self%nbprocs_c_list(self%max_nbprocs), &
           !     self%course_nb%info(self%max_nbprocs))
           print *, 'add_ipe_to_c_list: nbprocs_c exceeding max_nbprocs', self%nbprocs_c, self%max_nbprocs
        end if

        ! add the process
        self%nbprocs_c_list(self%nbprocs_c) = ipe
        ! add to inverse list
        self%ipe_to_inbpe_c(ipe) = self%nbprocs_c

     end if

     call self%add_igrid_to_c(ipe, igrid, i1, i2, i3, inc1, inc2, inc3)

   end subroutine add_to_c

   subroutine add_to_f(self, ipe, igrid, inc1, inc2, inc3)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: ipe, igrid, inc1, inc2, inc3

     ! add the procesor id if not already present:
     if (self%ipe_to_inbpe_f(ipe) == -1) then

        ! enlarge counter
        self%nbprocs_f = self%nbprocs_f + 1

        if (self%nbprocs_f > self%max_nbprocs) then
           print *, 'add_ipe_to_c_list: nbprocs_f exceeding max_nbprocs', self%nbprocs_f, self%max_nbprocs
           stop
        end if

        ! add the process
        self%nbprocs_f_list(self%nbprocs_f) = ipe
        ! add to inverse list
        self%ipe_to_inbpe_f(ipe) = self%nbprocs_f

     end if

     call self%add_igrid_to_f(ipe, igrid, inc1, inc2, inc3)

   end subroutine add_to_f

   subroutine add_igrid_to_srl(self, ipe, igrid, i1, i2, i3)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: ipe, igrid, i1, i2, i3
     ! .. local ..
     integer               :: inbpe, i, nigrids
     integer, allocatable, dimension(:) :: oldigrid, oldiencode

     ! translate to neighbor processor index for srl
     inbpe = self%ipe_to_inbpe_srl(ipe)

     ! enlarge counter
     self%srl_nb(inbpe)%info%nigrids = self%srl_nb(inbpe)%info%nigrids + 1
     nigrids = self%srl_nb(inbpe)%info%nigrids   

     ! need to enlarge storage
     if ( self%srl_nb(inbpe)%info%nigrids > size(self%srl_nb(inbpe)%info%igrid) ) then
        call init_srl_info(self%srl_nb(inbpe)%info, self%srl_nb(inbpe)%info%nigrids)
     end if


     call self%iencode(i1,i2,i3,i)
     ! add the data at the counter
     self%srl_nb(inbpe)%info%igrid( self%srl_nb(inbpe)%info%nigrids )    = igrid
     self%srl_nb(inbpe)%info%iencode( self%srl_nb(inbpe)%info%nigrids )  = i

   end subroutine add_igrid_to_srl

   subroutine add_igrid_to_c(self, ipe, igrid, i1, i2, i3, inc1, inc2, inc3)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: ipe, igrid, i1, i2, i3, inc1, inc2, inc3
     ! .. local ..
     integer               :: inbpe, i

     ! translate to neighbor processor index for c
     inbpe = self%ipe_to_inbpe_c(ipe)

     ! enlarge counter
     self%course_nb(inbpe)%info%nigrids = self%course_nb(inbpe)%info%nigrids + 1

     ! need to enlarge storage
     if ( self%course_nb(inbpe)%info%nigrids > self%max_igrids ) then
        error stop 'add_igrid_to_c: nigrids exceeding max_igrids'
     end if

     ! add the data at the counter
     self%course_nb(inbpe)%info%igrid( self%course_nb(inbpe)%info%nigrids )    = igrid
     self%course_nb(inbpe)%info%i1( self%course_nb(inbpe)%info%nigrids )       = i1
     self%course_nb(inbpe)%info%i2( self%course_nb(inbpe)%info%nigrids )       = i2
     self%course_nb(inbpe)%info%i3( self%course_nb(inbpe)%info%nigrids )       = i3
     self%course_nb(inbpe)%info%inc1( self%course_nb(inbpe)%info%nigrids )     = inc1
     self%course_nb(inbpe)%info%inc2( self%course_nb(inbpe)%info%nigrids )     = inc2
     self%course_nb(inbpe)%info%inc3( self%course_nb(inbpe)%info%nigrids )     = inc3

   end subroutine add_igrid_to_c

   subroutine add_igrid_to_f(self, ipe, igrid, inc1, inc2, inc3)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: ipe, igrid, inc1, inc2, inc3
     ! .. local ..
     integer               :: inbpe, i

     ! translate to neighbor processor index for f
     inbpe = self%ipe_to_inbpe_f(ipe)

     ! enlarge counter
     self%f(inbpe)%nigrids = self%f(inbpe)%nigrids + 1

     ! need to enlarge storage
     if ( self%f(inbpe)%nigrids > self%max_igrids ) then
        print *, 'add_igrid_to_f: nigrids exceeding max_igrids'
        stop
     end if

     ! add the data at the counter
     self%f(inbpe)%igrid( self%f(inbpe)%nigrids )    = igrid
     self%f(inbpe)%inc1( self%f(inbpe)%nigrids )     = inc1
     self%f(inbpe)%inc2( self%f(inbpe)%nigrids )     = inc2
     self%f(inbpe)%inc3( self%f(inbpe)%nigrids )     = inc3

   end subroutine add_igrid_to_f

   subroutine alloc_buffers_srl(self, nwgc, &
         ixS_srl_min1, ixS_srl_max1, &
         ixS_srl_min2, ixS_srl_max2, &
         ixS_srl_min3, ixS_srl_max3, &
         ixR_srl_min1, ixR_srl_max1, &
         ixR_srl_min2, ixR_srl_max2, &
         ixR_srl_min3, ixR_srl_max3)

     class(nbprocs_info_t)                     :: self
     integer, intent(in)                       :: nwgc
     integer, dimension(-1:2,-1:1), intent(in) :: &
         ixS_srl_min1, ixS_srl_max1, &
         ixS_srl_min2, ixS_srl_max2, &
         ixS_srl_min3, ixS_srl_max3, &
         ixR_srl_min1, ixR_srl_max1, &
         ixR_srl_min2, ixR_srl_max2, &
         ixR_srl_min3, ixR_srl_max3
     ! .. local ..
     integer         :: inb, igrid, isize_S, isize_R
     integer         :: n_i1, n_i2, n_i3, i1, i2, i3
     integer         :: ixSmin1, ixSmin2, ixSmin3
     integer         :: ixRmin1, ixRmin2, ixRmin3
     integer         :: ixSmax1, ixSmax2, ixSmax3
     integer         :: ixRmax1, ixRmax2, ixRmax3
     integer         :: iib1, iib2, iib3, ibuf_start

     self%imaxigrids_srl = 0
     do inb = 1, self%nbprocs_srl
        self%imaxigrids_srl = max(self%srl_nb(inb)%info%nigrids, self%imaxigrids_srl)

        isize_S = 0; isize_R = 0; ibuf_start = 1
        do igrid = 1, self%srl_nb(inb)%info%nigrids
           iib1 = idphyb(1,self%srl_nb(inb)%info%igrid(igrid))
           iib2 = idphyb(2,self%srl_nb(inb)%info%igrid(igrid))
           iib3 = idphyb(3,self%srl_nb(inb)%info%igrid(igrid))
           call self%idecode( i1, i2, i3, self%srl_nb(inb)%info%iencode(igrid) )

           n_i1=-i1; n_i2=-i2; n_i3=-i3

           ixSmin1=ixS_srl_min1(iib1,i1);ixSmin2=ixS_srl_min2(iib2,i2)
           ixSmin3=ixS_srl_min3(iib3,i3);ixSmax1=ixS_srl_max1(iib1,i1)
           ixSmax2=ixS_srl_max2(iib2,i2);ixSmax3=ixS_srl_max3(iib3,i3)
           ixRmin1=ixR_srl_min1(iib1,n_i1);ixRmin2=ixR_srl_min2(iib2,n_i2)
           ixRmin3=ixR_srl_min3(iib3,n_i3);ixRmax1=ixR_srl_max1(iib1,n_i1)
           ixRmax2=ixR_srl_max2(iib2,n_i2);ixRmax3=ixR_srl_max3(iib3,n_i3)

           isize_S = isize_S + (ixSmax1-ixSmin1+1) * (ixSmax2-ixSmin2+1) * (ixSmax3-ixSmin3+1)
           isize_R = isize_R + (ixRmax1-ixRmin1+1) * (ixRmax2-ixRmin2+1) * (ixRmax3-ixRmin3+1)

           self%srl_nb(inb)%info%ibuf_start(igrid) = ibuf_start
           self%srl_nb(inb)%info%isize(igrid) = (ixSmax1-ixSmin1+1) * (ixSmax2-ixSmin2+1) * (ixSmax3-ixSmin3+1) * nwgc

           ibuf_start = ibuf_start + self%srl_nb(inb)%info%isize(igrid)

        end do

        if (isize_S*nwgc > self%max_size .or. isize_R*nwgc > self%max_size) then
           print *, 'alloc_buffers_srl: exceeding max_size for ghost-cell buffer', isize_S*nwgc, isize_R*nwgc, self%max_size
           stop
        end if

        self%srl_nb(inb)%rcv%size       = isize_R * nwgc
        self%srl_nb(inb)%send%size      = isize_S * nwgc
        self%srl_nb(inb)%info_send%size = 3 * self%srl_nb(inb)%info%nigrids
        self%srl_nb(inb)%info_rcv%size  = 3 * self%srl_nb(inb)%info%nigrids

     end do

   end subroutine alloc_buffers_srl

   subroutine alloc_buffers_f(self, nwgc, &
         ixS_p_min1, ixS_p_max1, &
         ixS_p_min2, ixS_p_max2, &
         ixS_p_min3, ixS_p_max3, &
         ixR_r_min1, ixR_r_max1, &
         ixR_r_min2, ixR_r_max2, &
         ixR_r_min3, ixR_r_max3)
     ! neighbor is fine
     ! required buffers:
     ! send: psb(ixS_p_) with indices type_send_p
     ! receive: psb(ixR_r_) with indices type_recv_r

     class(nbprocs_info_t)                     :: self
     integer, intent(in)                       :: nwgc
     integer, dimension(-1:1,0:3), intent(in)  :: &
         ixS_p_min1, ixS_p_max1, &
         ixS_p_min2, ixS_p_max2, &
         ixS_p_min3, ixS_p_max3, &
         ixR_r_min1, ixR_r_max1, &
         ixR_r_min2, ixR_r_max2, &
         ixR_r_min3, ixR_r_max3
     ! .. local ..
     integer         :: inb, isize_S, isize_R
     integer         :: igrid, inc1, inc2, inc3
     integer         :: iib1, iib2, iib3, ibuf_start
     integer         :: ixSmin1, ixSmin2, ixSmin3
     integer         :: ixRmin1, ixRmin2, ixRmin3
     integer         :: ixSmax1, ixSmax2, ixSmax3
     integer         :: ixRmax1, ixRmax2, ixRmax3

     self%imaxigrids_f = 0
     do inb = 1, self%nbprocs_f
        self%imaxigrids_f = max(self%f(inb)%nigrids, self%imaxigrids_f)

        isize_S = 0; isize_R = 0; ibuf_start = 1
        do igrid = 1, self%f(inb)%nigrids
           iib1 = idphyb(1,self%f(inb)%igrid(igrid))
           iib2 = idphyb(2,self%f(inb)%igrid(igrid))
           iib3 = idphyb(3,self%f(inb)%igrid(igrid))
           inc1 = self%f(inb)%inc1(igrid)
           inc2 = self%f(inb)%inc2(igrid)
           inc3 = self%f(inb)%inc3(igrid)

           ixSmin1=ixS_p_min1(iib1,inc1);ixSmin2=ixS_p_min2(iib2,inc2)
           ixSmin3=ixS_p_min3(iib3,inc3);ixSmax1=ixS_p_max1(iib1,inc1)
           ixSmax2=ixS_p_max2(iib2,inc2);ixSmax3=ixS_p_max3(iib3,inc3)

           ixRmin1=ixR_r_min1(iib1,inc1);ixRmin2=ixR_r_min2(iib2,inc2)
           ixRmin3=ixR_r_min3(iib3,inc3);ixRmax1=ixR_r_max1(iib1,inc1)
           ixRmax2=ixR_r_max2(iib2,inc2);ixRmax3=ixR_r_max3(iib3,inc3)

           isize_S = isize_S + (ixSmax1-ixSmin1+1) * (ixSmax2-ixSmin2+1) * (ixSmax3-ixSmin3+1)
           isize_R = isize_R + (ixRmax1-ixRmin1+1) * (ixRmax2-ixRmin2+1) * (ixRmax3-ixRmin3+1)

           self%f(inb)%ibuf_start(igrid) = ibuf_start
           self%f(inb)%isize(igrid) = (ixSmax1-ixSmin1+1) * (ixSmax2-ixSmin2+1) * (ixSmax3-ixSmin3+1) * nwgc

           ibuf_start = ibuf_start + self%f(inb)%isize(igrid)

        end do

        if (isize_S*nwgc > self%max_size .or. isize_R*nwgc > self%max_size) then
           print *, 'alloc_buffers_f: exceeding max_size for ghost-cell buffer', isize_S*nwgc, isize_R*nwgc, self%max_size
           stop
        end if

        self%f_rcv(inb)%size       = isize_R*nwgc
        self%f_send(inb)%size      = isize_S*nwgc
        self%f_info_send(inb)%size = 5 * self%f(inb)%nigrids
        self%f_info_rcv(inb)%size  = 5 * self%f(inb)%nigrids

     end do

   end subroutine alloc_buffers_f

   subroutine alloc_buffers_c(self, nwgc, &
         ixS_r_min1, ixS_r_max1, &
         ixS_r_min2, ixS_r_max2, &
         ixS_r_min3, ixS_r_max3, &
         ixR_p_min1, ixR_p_max1, &
         ixR_p_min2, ixR_p_max2, &
         ixR_p_min3, ixR_p_max3 &
         )
     ! neighbor is coarse
     ! required buffers:
     ! send: psc(ixS_r_) with indices type_send_r
     ! receive: psc(ixR_p_) with indices type_recv_p

     class(nbprocs_info_t)                     :: self
     integer, intent(in)                       :: nwgc
     integer, dimension(-1:1,-1:1), intent(in) :: &
         ixS_r_min1, ixS_r_max1, &
         ixS_r_min2, ixS_r_max2, &
         ixS_r_min3, ixS_r_max3
     integer, dimension(-1:1,0:3), intent(in)  :: &
         ixR_p_min1, ixR_p_max1, &
         ixR_p_min2, ixR_p_max2, &
         ixR_p_min3, ixR_p_max3
     ! .. local ..
     integer         :: inb, isize_S, isize_R
     integer         :: igrid, i1, i2, i3
     integer         :: inc1, inc2, inc3
     integer         :: iib1, iib2, iib3, ibuf_start
     integer         :: ixSmin1, ixSmin2, ixSmin3
     integer         :: ixRmin1, ixRmin2, ixRmin3
     integer         :: ixSmax1, ixSmax2, ixSmax3
     integer         :: ixRmax1, ixRmax2, ixRmax3

     self%imaxigrids_c = 0
     do inb = 1, self%nbprocs_c
        self%imaxigrids_c = max(self%course_nb(inb)%info%nigrids, self%imaxigrids_c)

        isize_S = 0; isize_R = 0; ibuf_start = 1
        do igrid = 1, self%course_nb(inb)%info%nigrids
           iib1 = idphyb(1,self%course_nb(inb)%info%igrid(igrid))
           iib2 = idphyb(2,self%course_nb(inb)%info%igrid(igrid))
           iib3 = idphyb(3,self%course_nb(inb)%info%igrid(igrid))
           i1   = self%course_nb(inb)%info%i1(igrid)
           i2   = self%course_nb(inb)%info%i2(igrid)
           i3   = self%course_nb(inb)%info%i3(igrid)
           inc1 = self%course_nb(inb)%info%inc1(igrid)
           inc2 = self%course_nb(inb)%info%inc2(igrid)
           inc3 = self%course_nb(inb)%info%inc3(igrid)

           ixSmin1=ixS_r_min1(iib1,i1);ixSmin2=ixS_r_min2(iib2,i2)
           ixSmin3=ixS_r_min3(iib3,i3);ixSmax1=ixS_r_max1(iib1,i1)
           ixSmax2=ixS_r_max2(iib2,i2);ixSmax3=ixS_r_max3(iib3,i3)

           ixRmin1=ixR_p_min1(iib1,inc1);ixRmin2=ixR_p_min2(iib2,inc2)
           ixRmin3=ixR_p_min3(iib3,inc3);ixRmax1=ixR_p_max1(iib1,inc1)
           ixRmax2=ixR_p_max2(iib2,inc2);ixRmax3=ixR_p_max3(iib3,inc3)

           isize_S = isize_S + (ixSmax1-ixSmin1+1) * (ixSmax2-ixSmin2+1) * (ixSmax3-ixSmin3+1)
           isize_R = isize_R + (ixRmax1-ixRmin1+1) * (ixRmax2-ixRmin2+1) * (ixRmax3-ixRmin3+1)

           self%course_nb(inb)%info%ibuf_start(igrid) = ibuf_start
           self%course_nb(inb)%info%isize(igrid) = (ixSmax1-ixSmin1+1) * (ixSmax2-ixSmin2+1) * (ixSmax3-ixSmin3+1) * nwgc

           ibuf_start = ibuf_start + self%course_nb(inb)%info%isize(igrid)

        end do

        if (isize_S*nwgc > self%max_size .or. isize_R*nwgc > self%max_size) then
           print *, 'alloc_buffers_c: exceeding max_size for ghost-cell buffer', isize_S*nwgc, isize_R*nwgc, self%max_size
           stop
        end if

        self%course_nb(inb)%rcv%size       = isize_R*nwgc
        self%course_nb(inb)%send%size      = isize_S*nwgc
        self%course_nb(inb)%info_send%size = 5 * self%course_nb(inb)%info%nigrids
        self%course_nb(inb)%info_rcv%size  = 5 * self%course_nb(inb)%info%nigrids

     end do

   end subroutine alloc_buffers_c

   subroutine iencode(self, i1, i2, i3, i)
     class(nbprocs_info_t)                     :: self
     integer, intent(in)                       :: i1, i2, i3
     integer, intent(out)                      :: i
     ! .. local ..
     integer, parameter                        :: i1min=-1, i2min=-1, i3min=-1

     i = 1 + (i1-i1min) + 3 * (i2-i2min) + 9 * (i3-i3min)

   end subroutine iencode

   subroutine idecode(self, i1, i2, i3, i)
     class(nbprocs_info_t)                     :: self
     integer, intent(in)                       :: i
     integer, intent(out)                      :: i1, i2, i3
     ! .. local ..
     integer, parameter                        :: i1min=-1, i2min=-1, i3min=-1
     integer                                   :: id, idd

     i3  = ceiling(dble(i)/9.0d0) - 1 + i3min
     id  = i - 9 * (i3-i3min)
     i2  = ceiling(dble(id)/3.0d0) -1 + i2min
     idd = id - 3 * (i2-i2min)
     i1  = idd + i1min - 1

   end subroutine idecode

   subroutine nbprocs_update_device
     ! .. local ..
     logical               :: first = .true.
     integer               :: ipe_neighbor

    if (first) then
       !$acc update device(nbprocs_info)
       !$acc enter data copyin(nbprocs_info%nbprocs_srl_list, nbprocs_info%ipe_to_inbpe_srl)
       !$acc enter data copyin(nbprocs_info%nbprocs_f_list, nbprocs_info%ipe_to_inbpe_f)
       !$acc enter data copyin(nbprocs_info%nbprocs_c_list, nbprocs_info%ipe_to_inbpe_c)
       !$acc enter data copyin(nbprocs_info%course_nb)
       !$acc enter data copyin(nbprocs_info%srl_nb)
       !$acc enter data copyin(nbprocs_info%f_rcv, nbprocs_info%f_send, nbprocs_info%f, nbprocs_info%f_info_rcv, nbprocs_info%f_info_send)
       do ipe_neighbor = 1, nbprocs_info%max_nbprocs
          !$acc enter data copyin(nbprocs_info%srl_nb(ipe_neighbor)%rcv%buffer, nbprocs_info%srl_nb(ipe_neighbor)%send%buffer)
          !$acc enter data copyin( nbprocs_info%srl_nb(ipe_neighbor)%info_rcv%buffer, nbprocs_info%srl_nb(ipe_neighbor)%info_send%buffer)
          !$acc enter data copyin(nbprocs_info%srl_nb(ipe_neighbor)%info%igrid)
          !$acc enter data copyin(nbprocs_info%srl_nb(ipe_neighbor)%info%iencode)
          !$acc enter data copyin(nbprocs_info%srl_nb(ipe_neighbor)%info%isize)
          !$acc enter data copyin(nbprocs_info%srl_nb(ipe_neighbor)%info%ibuf_start)
          !$acc enter data copyin(nbprocs_info%srl_nb(ipe_neighbor)%info%nigrids)
       end do
       do ipe_neighbor = 1, nbprocs_info%max_nbprocs
          !$acc enter data copyin(nbprocs_info%f_rcv(ipe_neighbor)%buffer, nbprocs_info%f_info_rcv(ipe_neighbor)%buffer)
          !$acc enter data copyin(nbprocs_info%f_send(ipe_neighbor)%buffer, nbprocs_info%f_info_send(ipe_neighbor)%buffer)
          !$acc enter data copyin(nbprocs_info%f(ipe_neighbor)%igrid)
          !$acc enter data copyin(nbprocs_info%f(ipe_neighbor)%i1)
          !$acc enter data copyin(nbprocs_info%f(ipe_neighbor)%i2)
          !$acc enter data copyin(nbprocs_info%f(ipe_neighbor)%i3)
          !$acc enter data copyin(nbprocs_info%f(ipe_neighbor)%inc1)
          !$acc enter data copyin(nbprocs_info%f(ipe_neighbor)%inc2)
          !$acc enter data copyin(nbprocs_info%f(ipe_neighbor)%inc3)
          !$acc enter data copyin(nbprocs_info%f(ipe_neighbor)%isize)
          !$acc enter data copyin(nbprocs_info%f(ipe_neighbor)%ibuf_start)
          !$acc enter data copyin(nbprocs_info%f(ipe_neighbor)%nigrids)
       end do
       do ipe_neighbor = 1, nbprocs_info%max_nbprocs
          !$acc enter data copyin(nbprocs_info%course_nb(ipe_neighbor)%rcv%buffer, nbprocs_info%course_nb(ipe_neighbor)%send%buffer)
          !$acc enter data copyin( nbprocs_info%course_nb(ipe_neighbor)%info_rcv%buffer, nbprocs_info%course_nb(ipe_neighbor)%info_send%buffer)
          !$acc enter data copyin(nbprocs_info%course_nb(ipe_neighbor)%info%igrid)
          !$acc enter data copyin(nbprocs_info%course_nb(ipe_neighbor)%info%i1)
          !$acc enter data copyin(nbprocs_info%course_nb(ipe_neighbor)%info%i2)
          !$acc enter data copyin(nbprocs_info%course_nb(ipe_neighbor)%info%i3)
          !$acc enter data copyin(nbprocs_info%course_nb(ipe_neighbor)%info%inc1)
          !$acc enter data copyin(nbprocs_info%course_nb(ipe_neighbor)%info%inc2)
          !$acc enter data copyin(nbprocs_info%course_nb(ipe_neighbor)%info%inc3)
          !$acc enter data copyin(nbprocs_info%course_nb(ipe_neighbor)%info%isize)
          !$acc enter data copyin(nbprocs_info%course_nb(ipe_neighbor)%info%ibuf_start)
          !$acc enter data copyin(nbprocs_info%course_nb(ipe_neighbor)%info%nigrids)
       end do

       first = .false.

    else
       !$acc update device(nbprocs_info%nbprocs_srl, nbprocs_info%imaxigrids_srl)
       !$acc update device(nbprocs_info%nbprocs_f, nbprocs_info%imaxigrids_f)
       !$acc update device(nbprocs_info%nbprocs_c, nbprocs_info%imaxigrids_c)
       !$acc update device(nbprocs_info%nbprocs_srl_list, nbprocs_info%ipe_to_inbpe_srl)
       !$acc update device(nbprocs_info%nbprocs_f_list, nbprocs_info%ipe_to_inbpe_f)
       !$acc update device(nbprocs_info%nbprocs_c_list, nbprocs_info%ipe_to_inbpe_c)
       do ipe_neighbor = 1, nbprocs_info%nbprocs_srl
          !$acc update device(nbprocs_info%srl_nb(ipe_neighbor)%info%nigrids)
          !$acc update device(nbprocs_info%srl_nb(ipe_neighbor)%info%iencode)
          !$acc update device(nbprocs_info%srl_nb(ipe_neighbor)%info%isize)
          !$acc update device(nbprocs_info%srl_nb(ipe_neighbor)%info%ibuf_start)
          !$acc update device(nbprocs_info%srl_nb(ipe_neighbor)%send%size)
          !$acc update device(nbprocs_info%srl_nb(ipe_neighbor)%rcv%size)
          !$acc update device(nbprocs_info%srl_nb(ipe_neighbor)%info_send%size)
          !$acc update device(nbprocs_info%srl_nb(ipe_neighbor)%info_rcv%size)
       end do
       !$acc update device(nbprocs_info%nbprocs_f, nbprocs_info%imaxigrids_f)
       do ipe_neighbor = 1, nbprocs_info%nbprocs_f
          !$acc update device(nbprocs_info%f(ipe_neighbor)%nigrids)
          !$acc update device(nbprocs_info%f(ipe_neighbor)%igrid)
          !$acc update device(nbprocs_info%f(ipe_neighbor)%i1)
          !$acc update device(nbprocs_info%f(ipe_neighbor)%i2)
          !$acc update device(nbprocs_info%f(ipe_neighbor)%i3)
          !$acc update device(nbprocs_info%f(ipe_neighbor)%inc1)
          !$acc update device(nbprocs_info%f(ipe_neighbor)%inc2)
          !$acc update device(nbprocs_info%f(ipe_neighbor)%inc3)
          !$acc update device(nbprocs_info%f(ipe_neighbor)%isize)
          !$acc update device(nbprocs_info%f(ipe_neighbor)%ibuf_start)
          !$acc update device(nbprocs_info%f_send(ipe_neighbor)%size)
          !$acc update device(nbprocs_info%f_rcv(ipe_neighbor)%size)
          !$acc update device(nbprocs_info%f_info_send(ipe_neighbor)%size)
          !$acc update device(nbprocs_info%f_info_rcv(ipe_neighbor)%size)
       end do
       !$acc update device(nbprocs_info%nbprocs_c, nbprocs_info%imaxigrids_c)
       do ipe_neighbor = 1, nbprocs_info%nbprocs_c
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info%nigrids)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info%igrid)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info%i1)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info%i2)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info%i3)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info%inc1)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info%inc2)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info%inc3)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info%isize)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info%ibuf_start)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%send%size)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%rcv%size)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info_send%size)
          !$acc update device(nbprocs_info%course_nb(ipe_neighbor)%info_rcv%size)
       end do

    end if

    end subroutine nbprocs_update_device

 end module mod_connectivity
