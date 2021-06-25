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
     real(kind=8), ALLOCATABLE :: data1(:)   ! Data to write
     integer, allocatable :: data2(:)
     complex(kind=8), allocatable :: data3(:,:)

     CHARACTER(len=20) :: data4
     character(len=:), dimension(:, :), allocatable :: data5

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
     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL
     CALL MPI_INIT(mpierror)
     CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
     CALL MPI_COMM_RANK(comm, mpi_rank, mpierror)


     !
     ! Initialize data buffer with trivial data.
     !

     write(*,* ) "Rank = ", mpi_rank

     ! data0
     ALLOCATE ( data0(dimsf(1),dimsf(2)))
     do i = 1, dimsf(2)
          do j = 1, dimsf(1)
          data0(j,i) = j - 1 + (i-1)*dimsf(1)
          enddo
     enddo
     if (mpi_rank == 0) then
          data0(1,1) = -2
     end if

     ! data1
     allocate(data1(mpi_size*5))
     do i = 1, mpi_size*5
        data1(i) = dble(i)
     end do

     ! data2
     allocate(data2((mpi_rank+2)*2))
     do i = 1, (mpi_rank+2)*2
       data2(i) = mpi_rank + i
     end do

     ! data3
     allocate(data3(0:10, (mpi_rank+2)*2))
     do i = 1, (mpi_rank+2)*2
          do j = 0,10
               data3(j,i) = cmplx(dble(mpi_rank), dble(i), kind=8)
          end do
     end do

     ! data4
     write(data4, '("Rank=", I2)') mpi_rank

     ! data5: extremely complex array
     i = 20+mpi_rank
     ! i = 20
     allocate (character(len=i)::data5(0:10, (mpi_rank+2)*2))
     do i = 1, (mpi_rank+2)*2
          do j = 0,10
               write(data5(j, i), '("Rank=", I2,"_", I2)') mpi_rank, i
          end do
     end do


     !
     ! Initialize FORTRAN interface

     call hdf_open_file(file_id, "test_hl.h5", STATUS='NEW')
     call hdf_write_dataset(file_id, "data0", data0, processor=0)
     call hdf_write_dataset(file_id, "data1", data1)
     call hdf_write_dataset(file_id, "data2", data2, axis=1)

     ! The following call will fail because data3 doesn't have the same dimension along
     ! axises other than axis=1
     ! call hdf_write_dataset(file_id, "data3", data3, axis=1)

     ! The following call stacks data3 along axis=2 and write the result
     call hdf_write_dataset(file_id, "data3", data3, axis=2)

     ! The following call writes data3 on processor = 0 to the file
     call hdf_write_dataset(file_id, "data3_2", data3, processor=0)

     ! ! This is not a bad call since data4 has different values across processors
     ! ! there is no guarantee what will be saved
     call hdf_write_dataset(file_id, "data4", data4)
     call hdf_write_dataset(file_id, "data4_0", data4, processor=0)
     call hdf_write_dataset(file_id, "data4_1", data4, axis=1)

     call hdf_write_dataset(file_id, "data5", data5, axis=2)


     ! call hdf_write_attribute(file_id, "", "test", 1)
     call hdf_close_file(file_id)

     ! read test
     call hdf_open_file(file_id, "test_hl.h5", STATUS='OLD', ACTION='READ')
     call hdf_read_dataset(file_id, "data0", data0)
     call hdf_read_dataset(file_id, "data1", data1)

     ! data2 = 0
     ! call hdf_read_dataset(file_id, "data2", data2)

     deallocate(data2)
     ii = 0
     do jj = 1, mpi_size
          ii = ii + (jj-1+2)*2
     end do
     ii = ii / mpi_size
     allocate(data2(ii), offset(mpi_size))
     offset(1) = 0
     do jj = 2, mpi_size
          offset(jj) = offset(jj-1) + ii
     end do
     call hdf_read_dataset(file_id, "data2", data2, offset=offset)



     write(*,*) "Rank=", mpi_rank, "data2 = ", data2

     call hdf_set_dims(file_id, 'data3', dims)
     deallocate(data3)
     allocate(data3(dims(1), dims(2)))
     call hdf_read_dataset(file_id, "data3", data3)
     ! if (mpi_rank == 2) then
     !      write(*,*) 'data0:'
     !      write(*,*) data0
     !      write(*,*) 'data1:'
     !      write(*,*) data1

     ! end if
     call hdf_close_file(file_id)


     DEALLOCATE(data0, data1, data2, data3, data5)


     CALL MPI_FINALIZE(mpierror)

     END PROGRAM DATASET
