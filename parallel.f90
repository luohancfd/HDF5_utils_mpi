PROGRAM DATASET

  USE HDF5 ! This module contains all necessary modules
  use hdf5_utils_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'
  CHARACTER(LEN=10), PARAMETER :: filename = "sds.h5"  ! File name
  CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  INTEGER(HID_T) :: plist_id      ! Property list identifier

  INTEGER(HSIZE_T), DIMENSION(2) :: dimsf = (/5, 8/) ! Dataset dimensions.
!     INTEGER, DIMENSION(7) :: dimsfi = (/5,8,0,0,0,0,0/)
!     INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi = (/5,8/)
  INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi

  INTEGER, ALLOCATABLE :: data0(:, :)   ! Data to write
  real(kind=8), ALLOCATABLE :: data1(:)   ! Data to write
  integer, allocatable :: data2(:)
  complex(kind=8), allocatable :: data3(:, :)

  CHARACTER(len=20) :: data4
  character(len=:), dimension(:, :), allocatable :: data5
  complex(kind=8) :: data6

  integer, allocatable :: dims(:), offset(:)

  INTEGER :: rank = 2 ! Dataset rank

  INTEGER :: error, error_n  ! Error flags
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

  !
  ! Initialize data buffer with trivial data.
  !

  write (*, *) "Rank = ", mpi_rank

  ! data0
  ALLOCATE (data0(dimsf(1), dimsf(2)))
  do i = 1, dimsf(2)
    do j = 1, dimsf(1)
      data0(j, i) = j - 1 + (i - 1)*dimsf(1)
    end do
  end do
  if (mpi_rank == 0) then
    data0(1, 1) = -2
  end if

  ! data1
  allocate (data1(mpi_size*5))
  do i = 1, mpi_size*5
    data1(i) = dble(i)
  end do

  ! data2
  allocate (data2((mpi_rank + 2)*2))
  do i = 1, (mpi_rank + 2)*2
    data2(i) = i
  end do

  ! data3
  allocate (data3(0:10, (mpi_rank + 2)*2))
  do i = 1, (mpi_rank + 2)*2
    do j = 0, 10
      data3(j, i) = cmplx(dble(mpi_rank), dble(i), kind=8)
    end do
  end do

  ! data4
  write (data4, '("Rank=", I2)') mpi_rank

  ! data5: extremely complex array of character
  i = 20 + mpi_rank
  allocate (character(len=i)::data5(0:10, (mpi_rank + 2)*2))
  do i = 1, (mpi_rank + 2)*2
    do j = 0, 10
      write (data5(j, i), '("Rank=", I2,"_", I2)') mpi_rank, i
    end do
  end do

  ! data6 : complex number
  data6 = cmplx(mpi_rank, 1.0)
  ! ================ Write ==============================

  call hdf_open_file(file_id, "test_hl.h5", STATUS='NEW')

  ! write data0 with value on processor0
  call hdf_write_dataset(file_id, "data0", data0, processor=0)

  ! data1 has the same value across processors, write it
  call hdf_write_dataset(file_id, "data1", data1)

  ! data2 is a 1D array with different size and values on each processor
  ! we stack the value along axis = 1 and write it to the file
  call hdf_write_dataset(file_id, "data2", data2, axis=1)

  ! data3 is a 2D array of complex with the same size on the first dimension, but different size
  ! on the second dimension
  !
  ! The following call will fail because data3 doesn't have the same dimension along
  !    axises other than axis=1
  ! call hdf_write_dataset(file_id, "data3", data3, axis=1)
  !
  ! The following call stacks data3 along axis=2 and write the result
  !
  ! For example:
  ! On processor 1, data3 is [[1, 2, 3], [4, 5, 6]], dimension is (2, 3)
  ! On processor 2, data3 is [[-1, -2, -3, -4], [-5, -6, -7, -8]], dimension is (2, 4)
  ! After the call the following data will be written to the file
  ! [[1,2,3, -1,-2,-3,-4], [4,5,6, -5,-6,-7,-8]]
  call hdf_write_dataset(file_id, "data3", data3, axis=2)

  ! data4 is a CHARACTER(len=20) with different values on each processor
  !
  ! The following is a bad call since data4 has different values across processors
  ! there is no guarantee what will be saved. It is only acceptable if there is only one core
  ! call hdf_write_dataset(file_id, "data4", data4)
  !
  ! The following will stack along axis 1, thus an array with size = num_of_processors and type
  ! = character(len=20) will be writtent
  call hdf_write_dataset(file_id, "data4", data4, axis=1)

  ! data5 is an extremely complex array of characters
  ! character(len=20+mpi_rank)::data5(0:10, (mpi_rank+2)*2)
  ! The code will first find the max length of string
  ! Then use that as the fundamental type and stack data5 along axis = 2, similar to data3
  call hdf_write_dataset(file_id, "data5", data5, axis=2)


  ! data6 is a single complex number
  call hdf_write_dataset(file_id, "data6", data6, axis=1)


  ! call hdf_write_attribute(file_id, "", "test", 1)
  call hdf_close_file(file_id)

  ! ================ Read ==============================
  ! read test
  call hdf_open_file(file_id, "test_hl.h5", STATUS='OLD', ACTION='READ')

  ! data0
  call hdf_read_dataset(file_id, "data0", data0)
  if (mpi_rank > 0) then
    call MPI_Recv(ii, 1, MPI_INTEGER, mpi_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, mpi_status, mpierror)
  end if
  write (*, *) "Rank=", mpi_rank, "data0=", data0
  call flush (6)
  if (mpi_rank < mpi_size - 1) then
    call MPI_Send(ii, 1, MPI_INTEGER, mpi_rank + 1, 1, MPI_COMM_WORLD, mpierror)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, mpierror)

  ! data1
  call hdf_read_dataset(file_id, "data1", data1)
  if (mpi_rank > 0) then
    call MPI_Recv(ii, 1, MPI_INTEGER, mpi_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, mpi_status, mpierror)
  end if
  write (*, '(A, I2, A, 20(F6.2, 2X))') "Rank=", mpi_rank, "data1=", data1
  call flush (6)
  if (mpi_rank < mpi_size - 1) then
    call MPI_Send(ii, 1, MPI_INTEGER, mpi_rank + 1, 1, MPI_COMM_WORLD, mpierror)
  end if

  ! data2
  call hdf_read_dataset(file_id, "data2", data2)
  if (mpi_rank > 0) then
    call MPI_Recv(ii, 1, MPI_INTEGER, mpi_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, mpi_status, mpierror)
  end if
  write (*, *) "Rank=", mpi_rank, "data2=", data2
  call flush (6)
  if (mpi_rank < mpi_size - 1) then
    call MPI_Send(ii, 1, MPI_INTEGER, mpi_rank + 1, 1, MPI_COMM_WORLD, mpierror)
  end if

  ! data3
  call hdf_read_dataset(file_id, "data3", data3)
  if (mpi_rank > 0) then
    call MPI_Recv(ii, 1, MPI_INTEGER, mpi_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, mpi_status, mpierror)
  end if
  write (*, *) "Rank=", mpi_rank, "data3:"
  do ii = 0,10
    write(*,'("  ", 20(F6.2,"+",F5.2,"i",",", 2X))') data3(ii,:)
  end do
  call flush (6)
  if (mpi_rank < mpi_size - 1) then
    call MPI_Send(ii, 1, MPI_INTEGER, mpi_rank + 1, 1, MPI_COMM_WORLD, mpierror)
  end if

  ! data6
  call hdf_read_dataset(file_id, "data6", data6)
  if (mpi_rank > 0) then
    call MPI_Recv(ii, 1, MPI_INTEGER, mpi_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, mpi_status, mpierror)
  end if
  write (*, *) "Rank=", mpi_rank, "data6:", data6
  call flush (6)
  if (mpi_rank < mpi_size - 1) then
    call MPI_Send(ii, 1, MPI_INTEGER, mpi_rank + 1, 1, MPI_COMM_WORLD, mpierror)
  end if




  ! ! data2 = 0
  ! ! call hdf_read_dataset(file_id, "data2", data2)

  ! deallocate(data2)
  ! ii = 0
  ! do jj = 1, mpi_size
  !      ii = ii + (jj-1+2)*2
  ! end do
  ! ii = ii / mpi_size
  ! allocate(data2(ii), offset(mpi_size))
  ! offset(1) = 0
  ! do jj = 2, mpi_size
  !      offset(jj) = offset(jj-1) + ii
  ! end do
  ! call hdf_read_dataset(file_id, "data2", data2, offset=offset)

  ! write(*,*) "Rank=", mpi_rank, "data2 = ", data2

  ! call hdf_set_dims(file_id, 'data3', dims)

  ! call hdf_read_dataset(file_id, "data3", data3)

  call hdf_close_file(file_id)

  DEALLOCATE (data0, data1, data2, data3, data5)

  CALL MPI_FINALIZE(mpierror)

END PROGRAM DATASET
