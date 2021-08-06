PROGRAM DATASET

  ! Run this program with 6 cores

  use hdf5_utils_mpi, only : hdf_open_file, hdf_write_dataset, &
    hdf_read_dataset, hdf_close_file, hdf_set_dims, hdf_get_rank, hdf_read_attribute, &
    hdf_set_even_offset

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(kind=8) :: file_id       ! File identifier

  INTEGER, DIMENSION(2) :: dimsf = (/5, 8/) ! Dataset dimensions.
  INTEGER, DIMENSION(2) :: dimsfi

  INTEGER, ALLOCATABLE :: data0(:, :)   ! Data to write
  real(kind=8), ALLOCATABLE :: data1(:)   ! Data to write
  integer, allocatable :: data2(:,:), data22(:,:)
  complex(kind=8), allocatable :: data3(:, :)

  CHARACTER(len=20) :: data4
  character(len=:), dimension(:, :), allocatable :: data5
  complex(kind=8) :: data6
  CHARACTER(len=20) :: data7
  character(len=:), dimension(:), allocatable :: data8
  integer, allocatable :: dims(:), offset_global(:), count_global(:)
  integer :: axis_write, mpi_old_size

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

  if (mpi_size .ne. 6) then
    if (mpi_rank == 0) then
      write(*,*) "Run this program with 6 cores"
      call MPI_Abort(comm, mpierror)
    end if
    call MPI_Barrier(comm, mpierror)
  end if

  allocate(offset_global(mpi_size), count_global(mpi_size))
  ! ================ Read ==============================
  ! read test
  call hdf_open_file(file_id, "test_hl.h5", STATUS='OLD', ACTION='READ')

  ! get the number of processors used to write the file
  call hdf_read_attribute(file_id, "", "mpi_nrank", mpi_old_size)
  if (mpi_old_size .ne. mpi_size) then
    if (mpi_rank == 0) then
      write(*,*) "We are reading with different number of cores"
    end if
  end if
  call MPI_Barrier(comm, mpierror)

  ! data2 ====================================
  ! This method dynamically set data2 buffer

  ! get the original dimension
  call hdf_set_dims(file_id, "data2", dims)
  ! get the axis that data is stacked, it should be 1 here
  call hdf_read_attribute(file_id, "data2", "axis_write", axis_write)
  if (axis_write .ne. -1) then
    if (mpi_rank .eq. 0) then
      write(*,*) "we need to do load balancing"
    end if
  end if

  ! load balancing
  if (mpi_rank == 0) then
    ! calculate offset of the stacked axis on processor 0
    jj = floor(dble(dims(axis_write))/dble(mpi_size))
    ! write(*,*) dims(axis_write), jj
    offset_global(1) = 0
    offset_global(2) = dims(axis_write) - jj*(mpi_size-1)
    do ii = 3, mpi_size
      offset_global(ii) = offset_global(ii-1) + jj
    end do

    count_global(1:mpi_size-1) = offset_global(2:mpi_size) - offset_global(1:mpi_size-1)

    count_global(mpi_size) = dims(axis_write) - offset_global(mpi_size)
    write(*,*) offset_global
  end if

  ! Broadcast offset_global
  call MPI_Bcast(offset_global, mpi_size, MPI_INTEGER, 0, comm, mpierror)

  ! Scatter memory size to each processor
  call MPI_Scatter(count_global, 1, MPI_INTEGER, jj, 1, MPI_INTEGER, 0, comm, mpierror)

  ! we know axis_write=1 and data2 has rank=1
  allocate(data2(jj, dims(2)))

  ! read the datatset
  call hdf_read_dataset(file_id, "data2", data2, offset=offset_global)
  if (mpi_rank > 0) then
    call MPI_Recv(ii, 1, MPI_INTEGER, mpi_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, mpi_status, mpierror)
  end if
  write (*, *) "Rank=", mpi_rank, "data2="
  do i = 1, jj
    write(*,*) (data2(i,j), j=1,dims(2))
  end do

  call flush (6)
  if (mpi_rank < mpi_size - 1) then
    call MPI_Send(ii, 1, MPI_INTEGER, mpi_rank + 1, 1, MPI_COMM_WORLD, mpierror)
  end if


  ! data22 ================ this method preset data22
  call hdf_set_dims(file_id, "data2", dims)
  call hdf_read_attribute(file_id, "data2", "axis_write", axis_write)

  if (axis_write .eq. 1) then
    ALLOCATE(data22(10, dims(2)))
  end if
  if (axis_write == 2) then
     allocate(data22(dims(1), 10))
  end if
  data22 = -1
  call hdf_set_even_offset(file_id, file_id, "data2", offset_global, 10)
  call hdf_read_dataset(file_id, "data2", data22, offset=offset_global)
  if (mpi_rank > 0) then
    call MPI_Recv(ii, 1, MPI_INTEGER, mpi_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, mpi_status, mpierror)
  end if
  write (*, *) "Rank=", mpi_rank, "data2="
  do i = 1, 10
    write(*,*) (data22(i,j), j=1,dims(2))
  end do
  call flush (6)
  if (mpi_rank < mpi_size - 1) then
    call MPI_Send(ii, 1, MPI_INTEGER, mpi_rank + 1, 1, MPI_COMM_WORLD, mpierror)
  end if




  call hdf_close_file(file_id)

  DEALLOCATE (data2, data22, offset_global, count_global)

  CALL MPI_FINALIZE(mpierror)

END PROGRAM DATASET
