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
	integer, allocatable :: prtcl_offset(:)

	integer :: mpi_size_old, total_prtcl(1)

	comm = MPI_COMM_WORLD

	call mpi_init(mpierror)

	call mpi_comm_size(comm, mpi_size, mpierror)

	call mpi_comm_rank(comm, mpi_rank, mpierror)

	ipart_particle = (mpi_size+1)*2 + num_ch_pic + 10
    num_ch_pic = 3

	num_particles = 0

	file_restart = "astra_res.h5"
    call hdf_open_file(h5file_id, file_restart, STATUS='OLD', ACTION='READ')
    call hdf_read_dataset(h5file_id, 'num_ch_pic', num_ch_pic)
    call hdf_read_dataset(h5file_id, 'ipart_particle', ipart_particle)
    call hdf_read_attribute(h5file_id, "", "mpi_nrank", mpi_size_old)

	ipart_particle = ceiling(dble(ipart_particle * mpi_size_old) / dble(mpi_size))

	allocate(data(ipart_particle, num_ch_pic))
    data = 0

	allocate(prtcl_offset(mpi_size))
	do ii = 1, num_ch_pic
		write(group_name,'("prtcl_",I3.3)') ii
		call hdf_open_group(h5file_id, group_name, group_id)

		call hdf_get_dims(group_id, "data", total_prtcl)
		if (mpi_size .eq. mpi_size_old) then
			call hdf_read_attribute(group_id, "data", "offset", prtcl_offset)
		else
			call hdf_set_even_offset(h5file_id, group_id, "data", prtcl_offset, ipart_particle)
		end if

		if (mpi_rank+1 == mpi_size) then
			num_particles(ii) = total_prtcl(1) - prtcl_offset(mpi_rank+1)
		else
			num_particles(ii) = prtcl_offset(mpi_rank+2) - prtcl_offset(mpi_rank+1)
		end if
		jj = num_particles(ii)

		call hdf_read_dataset(group_id, "data", data(1:jj,ii), offset=prtcl_offset)
        call hdf_close_group(group_id)
	end do

	if (mpi_rank > 0) then
		call MPI_Recv(ii, 1, MPI_INTEGER, mpi_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, mpi_status, mpierror)
	end if
	write (*, *) "Rank=", mpi_rank, "data:"
	! write(*,*) (num_particles(ii), ii = 1, num_ch_pic)
	do ii = 1, ipart_particle
		write(*,*) (data(ii, jj), jj = 1, num_ch_pic)
	end do
	call flush (6)
	if (mpi_rank < mpi_size - 1) then
		call MPI_Send(ii, 1, MPI_INTEGER, mpi_rank + 1, 1, MPI_COMM_WORLD, mpierror)
	end if

    call hdf_close_file(h5file_id)
    deallocate(data, prtcl_offset)
    CALL MPI_FINALIZE(mpierror)
end program


