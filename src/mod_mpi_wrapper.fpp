module mod_mpi_wrapper
  ! This module is a temporary workaround for a bug in NVHPC 25.3:
  ! The runtime can fail when an offset is given to a buffer used in MPI
  ! commands. By using these wrappers, the actual MPI call has no offset.

  use mpi
  implicit none
  private

  public :: mpi_file_read_wrapper
  public :: mpi_file_read_at_wrapper
  public :: mpi_file_write_wrapper
  public :: mpi_file_write_at_wrapper
  public :: mpi_file_write_at_all_wrapper
  public :: mpi_irecv_wrapper
  public :: mpi_isend_wrapper
  public :: mpi_recv_wrapper
  public :: mpi_send_wrapper

contains

  subroutine mpi_file_read_wrapper(fh, buf, count, datatype, status, ierror)
    integer, intent(in)                              :: fh, count, datatype
    type(*), intent(inout)                           :: buf(..)
    integer, dimension(MPI_STATUS_SIZE), intent(out) :: status
    integer, optional, intent(out)                   :: ierror

    call MPI_FILE_READ(fh, buf, count, datatype, status, ierror)
  end subroutine mpi_file_read_wrapper

  subroutine mpi_file_read_at_wrapper(fh, offset, buf, count, datatype, status, ierror)
    integer, intent(in)                              :: fh, count, datatype
    integer(MPI_OFFSET_KIND), intent(in)             :: offset
    type(*), intent(inout)                           :: buf(..)
    integer, dimension(MPI_STATUS_SIZE), intent(out) :: status
    integer, optional, intent(out)                   :: ierror

    call MPI_FILE_READ_AT(fh, offset, buf, count, datatype, status, ierror)
  end subroutine mpi_file_read_at_wrapper

  subroutine mpi_file_write_wrapper(fh, buf, count, datatype, status, ierror)
    integer, intent(in)                              :: fh, count, datatype
    type(*), intent(in)                              :: buf(..)
    integer, dimension(MPI_STATUS_SIZE), intent(out) :: status
    integer, optional, intent(out)                   :: ierror

    call MPI_FILE_WRITE(fh, buf, count, datatype, status, ierror)
  end subroutine mpi_file_write_wrapper

  subroutine mpi_file_write_at_wrapper(fh, offset, buf, count, datatype, status, ierror)
    integer, intent(in)                              :: fh, count, datatype
    integer(MPI_OFFSET_KIND), intent(in)             :: offset
    type(*), intent(in)                              :: buf(..)
    integer, dimension(MPI_STATUS_SIZE), intent(out) :: status
    integer, optional, intent(out)                   :: ierror

    call MPI_FILE_WRITE_AT(fh, offset, buf, count, datatype, status, ierror)
  end subroutine mpi_file_write_at_wrapper

  subroutine mpi_file_write_at_all_wrapper(fh, offset, buf, count, datatype, status, ierror)
    integer, intent(in)                              :: fh, count, datatype
    integer(MPI_OFFSET_KIND), intent(in)             :: offset
    type(*), intent(in)                              :: buf(..)
    integer, dimension(MPI_STATUS_SIZE), intent(out) :: status
    integer, optional, intent(out)                   :: ierror

    call MPI_FILE_WRITE_AT_ALL(fh, offset, buf, count, datatype, status, ierror)
  end subroutine mpi_file_write_at_all_wrapper

  subroutine mpi_irecv_wrapper(buf, count, datatype, source, tag, comm, request, ierror)
    integer, intent(in)            :: count, source, tag, datatype, comm
    type(*), intent(inout)         :: buf(..)
    integer, intent(out)           :: request
    integer, optional, intent(out) :: ierror

    call MPI_IRECV(buf, count, datatype, source, tag, comm, request, ierror)
  end subroutine mpi_irecv_wrapper

  subroutine mpi_isend_wrapper(buf, count, datatype, dest, tag, comm, request, ierror)
    integer, intent(in)            :: count, dest, tag, datatype, comm
    type(*), intent(in)            :: buf(..)
    integer, intent(out)           :: request
    integer, optional, intent(out) :: ierror

    call MPI_ISEND(buf, count, datatype, dest, tag, comm, request, ierror)
  end subroutine mpi_isend_wrapper

  subroutine mpi_recv_wrapper(buf, count, datatype, source, tag, comm, status, ierror)
    integer, intent(in)                              :: count, source, tag, datatype, comm
    type(*), intent(inout)                           :: buf(..)
    integer, dimension(MPI_STATUS_SIZE), intent(out) :: status
    integer, optional, intent(out)                   :: ierror

    call MPI_RECV(buf, count, datatype, source, tag, comm, status, ierror)
  end subroutine mpi_recv_wrapper

  subroutine mpi_send_wrapper(buf, count, datatype, dest, tag, comm, ierror)
    integer, intent(in)            :: count, dest, tag, datatype, comm
    type(*), intent(in)            :: buf(..)
    integer, optional, intent(out) :: ierror

    call MPI_SEND(buf, count, datatype, dest, tag, comm, ierror)
  end subroutine mpi_send_wrapper

end module mod_mpi_wrapper