#if defined(_OPENACC) || defined(_OPENMP)
!=====================================================================
! Transfer wrappers for copy or update if present
!=====================================================================
module gpu_utils
#:if defined('OPENACC')
  use openacc
#:elif defined('OPENMP')
  use omp_lib
#:endif
  implicit none

  private
  public :: copy_or_update, copy_or_update_pointer, copy_or_update_alloc

  ! Adjust this as needed
  #:set MAXRANK = 6
  #:set TYPES = [("double", "double precision"), ("logical", "logical"), ("integer", "integer")]

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Generic interface for arrays and scalars (non-pointer, non-allocatable)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface copy_or_update
  #:for tname, tdecl in TYPES
  #:for rank in range(1, MAXRANK+1)
     module procedure copy_or_update_${tname}$_nonpointer_r${rank}$
  #:endfor
  #:endfor

  #:for tname, tdecl in TYPES
    module procedure copy_or_update_${tname}$_nonpointer
  #:endfor
  end interface copy_or_update
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Generic interface for data with pointer attribute
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface copy_or_update_pointer
  #:for tname, tdecl in TYPES
  #:for rank in range(1, MAXRANK+1)
     module procedure copy_or_update_${tname}$_pointer_r${rank}$
  #:endfor
  #:endfor

  #:for tname, tdecl in TYPES
    module procedure copy_or_update_${tname}$_pointer
  #:endfor
  end interface copy_or_update_pointer
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Generic interface for arrays with allocatable attribute
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface copy_or_update_alloc
  #:for tname, tdecl in TYPES
  #:for rank in range(1, MAXRANK+1)
     module procedure copy_or_update_${tname}$_alloc_r${rank}$
  #:endfor
  #:endfor
  end interface copy_or_update_alloc

contains

  ! Versions for scalars, pointers
  #:for tname, tdecl in TYPES

  subroutine copy_or_update_${tname}$_pointer(scalar, no_update)
    implicit none
    ${tdecl}$, pointer, intent(inout) :: scalar
    logical, optional, intent(in)     :: no_update

    logical                           :: no_update_

    if (present(no_update)) then
       no_update_ = no_update
    else
       no_update_ = .false.
    end if

    if (.not. associated(scalar)) then
       print *, 'copy_or_update: null-pointer encountered'
       return
    end if

#:if defined('OPENACC')
    if (.not. acc_is_present(scalar, sizeof(scalar))) then
       !$acc enter data copyin(scalar)
       !$acc enter data attach(scalar)
    else if (.not. no_update_) then
       !$acc update device(scalar)
    end if
    call acc_detach(scalar)
    call acc_attach(scalar)
#:elif defined('OPENMP')
    ! No presence check (unlike the OpenACC branch above): OpenMP has no
    ! Fortran-friendly equivalent of acc_is_present for arbitrary
    ! pointer/array/scalar arguments, and pointer attachment happens
    ! implicitly when both the pointer and its target are mapped. This is
    ! always correct for the transferred *value*, at the cost of an
    ! ever-growing device refcount on repeated calls (map(to:) on an
    ! already-present item only bumps the refcount, it does not recopy).
    !$omp target enter data map(to: scalar)
    if (.not. no_update_) then
       !$omp target update to(scalar)
    end if
#:endif

  end subroutine copy_or_update_${tname}$_pointer

  ! Versions for scalars, non-pointers
  subroutine copy_or_update_${tname}$_nonpointer(scalar, no_update)
    implicit none
    ${tdecl}$, intent(inout) :: scalar
    logical, optional, intent(in)     :: no_update

    logical                           :: no_update_

    if (present(no_update)) then
       no_update_ = no_update
    else
       no_update_ = .false.
    end if

#:if defined('OPENACC')
    if (.not. acc_is_present(scalar, sizeof(scalar))) then
       !$acc enter data copyin(scalar)
    else if (.not. no_update_) then
       !$acc update device(scalar)
    end if
#:elif defined('OPENMP')
    !$omp target enter data map(to: scalar)
    if (.not. no_update_) then
       !$omp target update to(scalar)
    end if
#:endif

  end subroutine copy_or_update_${tname}$_nonpointer

  #:endfor

  ! Versions for various array ranks

  #:for tname, tdecl in TYPES
  ! Generate non-pointer and pointer versions for ranks 1..MAXRANK
  #:for rank in range(1, MAXRANK+1)
  #:set shp = '(' + ','.join([':']*rank) + ')'

  !------------------------------
  ! Non-pointer, non-allocatable
  !------------------------------
  subroutine copy_or_update_${tname}$_nonpointer_r${rank}$(array, no_update)
    implicit none
    ${tdecl}$, intent(inout) :: array${shp}$
    logical, optional, intent(in)     :: no_update

    logical                           :: no_update_

    if (present(no_update)) then
       no_update_ = no_update
    else
       no_update_ = .false.
    end if

#:if defined('OPENACC')
    if (.not. acc_is_present(array)) then
       !$acc enter data copyin(array)
    else if (.not. no_update_) then
       !$acc update device(array)
    end if
#:elif defined('OPENMP')
    !$omp target enter data map(to: array)
    if (.not. no_update_) then
       !$omp target update to(array)
    end if
#:endif

  end subroutine copy_or_update_${tname}$_nonpointer_r${rank}$

  !-------------
  ! Pointer case
  !-------------
  subroutine copy_or_update_${tname}$_pointer_r${rank}$(array, no_update)
    implicit none
    ${tdecl}$, pointer, intent(inout) :: array${shp}$
    logical, optional, intent(in)     :: no_update

    logical                           :: no_update_

    if (present(no_update)) then
       no_update_ = no_update
    else
       no_update_ = .false.
    end if

    if (.not. associated(array)) then
       print *, 'copy_or_update: null-pointer encountered (rank ${rank}$)'
       return
    end if

#:if defined('OPENACC')
    if (.not. acc_is_present(array)) then
       !$acc enter data copyin(array)
       !$acc enter data attach(array)
    else if (.not. no_update_) then
       !$acc update device(array)
    end if
    call acc_detach(array)
    call acc_attach(array)
#:elif defined('OPENMP')
    ! See the comment in copy_or_update_${tname}$_pointer above: no presence
    ! check, and attach/detach is dropped since OpenMP attaches a mapped
    ! pointer to its (also mapped) target implicitly.
    !$omp target enter data map(to: array)
    if (.not. no_update_) then
       !$omp target update to(array)
    end if
#:endif

  end subroutine copy_or_update_${tname}$_pointer_r${rank}$

  !----------------
  ! Allocatable case
  !----------------
  subroutine copy_or_update_${tname}$_alloc_r${rank}$(array, no_update)
    implicit none
    ${tdecl}$, allocatable, intent(inout) :: array${shp}$
    logical, optional, intent(in)     :: no_update

    logical                           :: no_update_

    if (present(no_update)) then
       no_update_ = no_update
    else
       no_update_ = .false.
    end if

    if (.not. allocated(array)) then
       print *, 'copy_or_update: unallocated allocatable (rank ${rank}$)'
       return
    end if

#:if defined('OPENACC')
    if (.not. acc_is_present(array)) then
       !$acc enter data copyin(array)
    else if (.not. no_update_) then
       !$acc update device(array)
    end if
#:elif defined('OPENMP')
    !$omp target enter data map(to: array)
    if (.not. no_update_) then
       !$omp target update to(array)
    end if
#:endif

  end subroutine copy_or_update_${tname}$_alloc_r${rank}$

  #:endfor
  #:endfor

end module gpu_utils
#endif
