PROGRAM DATASET

  ! Run this program with 6 cores

  use hdf5_utils_mpi, only : hdf_open_file, hdf_write_dataset, &
    hdf_read_dataset, hdf_close_file, hdf_set_dims, hdf_get_rank, hdf_read_attribute, hdf_set_even_offset

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(kind=8) :: file_id       ! File identifier

  INTEGER, DIMENSION(2) :: dimsf = (/5, 8/) ! Dataset dimensions.
  INTEGER, DIMENSION(2) :: dimsfi

  INTEGER, ALLOCATABLE :: data0(:, :)   ! Data to write
  real(kind=8), ALLOCATABLE :: data1(:)   ! Data to write
  integer, allocatable :: data2(:, :)
  complex(kind=8), allocatable :: data3(:, :)

  CHARACTER(len=20) :: data4
  character(len=:), dimension(:, :), allocatable :: data5
  complex(kind=8) :: data6
  CHARACTER(len=20) :: data7
  character(len=:), dimension(:), allocatable :: data8
  integer, allocatable :: dims(:), offset_global(:), count_global(:)
  integer :: axis_write, mpi_old_size

  integer, allocatable :: offset(:)

  INTEGER :: i, j, ii, jj
  !
  ! MPI definitions and calls.
  !
  INTEGER :: mpierror       ! MPI error flag
  INTEGER :: comm, info
  INTEGER :: mpi_size, mpi_rank
  integer :: mpi_status(MPI_STATUS_SIZE)
  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL
  CALL MPI_INIT(mpierror)
  CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
  CALL MPI_COMM_RANK(comm, mpi_rank, mpierror)


  allocate(offset(mpi_size))
  ! ================ Read ==============================
  ! read test
  call hdf_open_file(file_id, "test_hl.h5", STATUS='OLD', ACTION='READ')

  allocate(data2(10, 5))
  data2 = -1

  call hdf_set_even_offset(file_id, file_id, "data2", offset)
  if (mpi_rank == 0) then
    write(*,*) offset
  end if

  call hdf_read_dataset(file_id, "data2", data2, offset)

  call hdf_close_file(file_id)

  call hdf_open_file(file_id, "test_hl_2.h5", STATUS='NEW')


  do i = 1 , mpi_size
    write(data4, '("data",I1)') i
    call hdf_write_dataset(file_id, trim(data4), data2, processor=i-1)
  end do

  call hdf_close_file(file_id)

  deallocate(data2)

  CALL MPI_FINALIZE(mpierror)

END PROGRAM DATASET
