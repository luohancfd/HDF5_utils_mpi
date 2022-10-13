program main
	use hdf5_utils_mpi
	implicit none

	include "mpif.h"

	integer, allocatable :: data(:,:)
	integer :: ipart_particle, num_ch_pic
	integer :: num_particles(200), ii, jj
    integer :: mpierror, comm, mpi_size, mpi_rank
    integer :: mpi_status(MPI_STATUS_SIZE)
	integer(HID_T) :: h5file_id, group_id
	character(len=50) :: file_restart, group_name

	comm = MPI_COMM_WORLD

	call mpi_init(mpierror)

	call mpi_comm_size(comm, mpi_size, mpierror)

	call mpi_comm_rank(comm, mpi_rank, mpierror)

	ipart_particle = (mpi_size+1)*2 + num_ch_pic + 10
    num_ch_pic = 3

	allocate(data(ipart_particle, num_ch_pic))
	data = 0

	num_particles = 0
	do ii = 1,num_ch_pic
		! some random value
		num_particles(ii) = (mpi_rank+1)*2 + ii
		do jj = 1, num_particles(ii)
			data(jj, ii) = mpi_rank * 100 + jj
		end do
	end do

	if (mpi_rank > 0) then
		call MPI_Recv(ii, 1, MPI_INTEGER, mpi_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, mpi_status, mpierror)
	end if
	write (*, *) "Rank=", mpi_rank, "data:"
	do ii = 1, ipart_particle
		write(*,*) (data(ii, jj), jj = 1, num_ch_pic)
	end do
	call flush (6)
	if (mpi_rank < mpi_size - 1) then
		call MPI_Send(ii, 1, MPI_INTEGER, mpi_rank + 1, 1, MPI_COMM_WORLD, mpierror)
	end if


	file_restart = "astra_res.h5"
    call hdf_open_file(h5file_id, file_restart, STATUS='NEW')
    call hdf_write_dataset(h5file_id, 'num_ch_pic', num_ch_pic)
    call hdf_write_dataset(h5file_id, 'ipart_particle', ipart_particle)
	do ii = 1, num_ch_pic
		write(group_name,'("prtcl_",I3.3)') ii
		call hdf_create_group(h5file_id, group_name)
		call hdf_open_group(h5file_id, group_name, group_id)

		jj = num_particles(ii)
		! only write valid particles
		call hdf_write_dataset(group_id, "data", data(1:jj,ii), axis=1)
        call hdf_close_group(group_id)
	end do
    call hdf_close_file(h5file_id)
    deallocate(data)
    CALL MPI_FINALIZE(mpierror)
end program


