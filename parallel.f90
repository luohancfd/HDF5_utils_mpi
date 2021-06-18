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

     INTEGER(HSIZE_T), DIMENSION(2) :: dimsf = (/5,8/) ! Dataset dimensions.
!     INTEGER, DIMENSION(7) :: dimsfi = (/5,8,0,0,0,0,0/)
!     INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi = (/5,8/)
     INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi

     INTEGER, ALLOCATABLE :: data0(:,:)   ! Data to write
     INTEGER, ALLOCATABLE :: data1(:)   ! Data to write

     INTEGER :: rank = 2 ! Dataset rank

     INTEGER :: error, error_n  ! Error flags
     INTEGER :: i, j
     !
     ! MPI definitions and calls.
     !
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank
     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL
     CALL MPI_INIT(mpierror)
     CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
     CALL MPI_COMM_RANK(comm, mpi_rank, mpierror)
     !
     ! Initialize data buffer with trivial data.
     !

     write(*,* ) "Rank = ", mpi_rank
     ALLOCATE ( data0(dimsf(1),dimsf(2)))
     allocate(data1(mpi_size*5))
     do i = 1, dimsf(2)
     do j = 1, dimsf(1)
        data0(j,i) = j - 1 + (i-1)*dimsf(1)
     enddo
     enddo
     if (mpi_rank == 0) then
          data0(1,1) = -2
     end if
     if (mpi_rank == 1) then
          data0(1,1) = -1
     end if

     do i = 1, mpi_size*5
          data1(i) = i
     end do
     !
     ! Initialize FORTRAN interface

     call hdf_open_file(file_id, "test_hl.h5", STATUS='NEW')
     call hdf_write_dataset(file_id, "data0", data0, processor=1)
     call hdf_write_dataset(file_id, "data1", data1)

     call hdf_write_dataset(file_id, "data2", mpi_rank, processor=2)
     call hdf_close_file(file_id)

     ! read test
     call hdf_open_file(file_id, "test_hl.h5", STATUS='OLD', ACTION='READ')
     call hdf_read_dataset(file_id, "data0", data0)
     call hdf_read_dataset(file_id, "data1", data1)
     if (mpi_rank == 2) then
          write(*,*) 'data0:'
          write(*,*) data0
          write(*,*) 'data1:'
          write(*,*) data1
     end if
     call hdf_close_file(file_id)


     DEALLOCATE(data0, data1)


     CALL MPI_FINALIZE(mpierror)

     END PROGRAM DATASET
