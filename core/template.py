
write_single_value_template = '''  subroutine hdf_write_dataset_{ftype_name}_0(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
{declaration}
    integer, optional, intent(in) :: chunks           ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis          ! axis used to stack the data among processors

    integer :: rank
    integer(HSIZE_T) :: dimsf(1), dimsm(1)
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: hdferror
    integer :: processor_write, axis_write
    integer(HSIZE_T), dimension(1) :: offset, count
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: ii
{additonal_delcaration}
    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_{ftype_name}_0: "//trim(dset_name)
    end if

    !
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_{ftype_name}_0: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_{ftype_name}_0: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      processor_write = processor
    end if

    ! set rank and dims
    rank = 1
    dimsm = (/0/)
    dimsf = (/mpi_nrank/)
{dtype_creation}
    ! set axis, might be changed later
    axis_write = -1
    if (processor_write == -1) then
      if (present(axis)) then
        axis_write = axis
      end if
    end if

    ! create dataspace and dataset
    call h5screate_f(H5S_SCALAR_F, mem_space_id, hdferror)
    if (axis_write .eq. -1) then
      file_space_id = mem_space_id
    else
      axis_write = 1
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, {h5type}, file_space_id, dset_id, hdferror)

    ! select hyperslab if needed
    if (axis_write .eq. -1) then
      mem_space_id = H5S_ALL_F
    else
      call h5sclose_f(file_space_id, hdferror)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      offset(1) = mpi_irank
      count(1)  = 1
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset, count, hdferror)

      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
      do ii = 1, mpi_nrank
        offset_glob(ii) = ii - 1
      end do
      count_glob = 1
    end if

    ! write dataset
{write_string}
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if
{additional_attribute}
    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)
{additional_close}
  end subroutine hdf_write_dataset_{ftype_name}_0
'''

write_array_template = '''  subroutine hdf_write_dataset_{ftype_name}_{rank}(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
{declaration}
    integer, optional, intent(in) :: chunks({rank})        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension({rank}) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
{additonal_delcaration}
    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_{ftype_name}_{rank}: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_{ftype_name}_{rank}: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_{ftype_name}_{rank}: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_{ftype_name}_{rank}: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = {rank}
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm
{dtype_creation}
    ! set axis_write, might be changed later
    axis_write = -1
    if (processor_write == -1) then
      if (present(axis)) then
        axis_write = axis
      end if
    end if

    ! create dataspace
    call h5screate_simple_f(rank, dimsm, mem_space_id, hdferror)
    if (axis_write .eq. -1) then
      file_space_id = mem_space_id
    else if (axis_write > {rank} .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", {rank}
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, {h5type}, file_space_id, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! select hyperslab if needed
    if (axis_write .eq. -1) then
      mem_space_id = H5S_ALL_F
    else
      call h5sclose_f(file_space_id, hdferror)
      call h5dget_space_f(dset_id, file_space_id, hdferror)

      ! set offset information
      offset = 0
      if (mpi_irank .ne. 0) then
        call MPI_Recv(offset_end, 1, mpi_hsize_t, mpi_irank-1, MPI_ANY_TAG, mpi_comm, status, mpi_ierr)
        offset(axis_write) = offset_end
      end if
      offset_end = offset(axis_write) + dimsm(axis_write)
      if (mpi_irank < mpi_nrank - 1) then
        call MPI_Send(offset_end, 1, mpi_hsize_t, mpi_irank+1, 0, mpi_comm, mpi_ierr)
      end if

      ! gather offset information
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))

      ! first check if other axises have the same dimension
      if (mpi_nrank > 1) then
        do ii = 1, rank
          if (ii .ne. axis_write) then
            call MPI_Gather(dimsm(ii),      1, mpi_hsize_t,   &
                            offset_glob,  1, mpi_hsize_t,   &
                            0, mpi_comm, mpi_ierr)
            if (mpi_irank == 0) then
              do jj = 2, mpi_nrank
                if (offset_glob(jj) .ne. offset_glob(1)) then
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_{ftype_name}_{rank}: "//trim(dset_name), &
                    " axis=", ii, " doesn't have the same size across processors"
                  call MPI_Abort(mpi_comm, 1, mpi_ierr)
                end if
              end do
            end if
          end if
        end do
      end if

      call MPI_Allgather(offset(axis_write), 1, mpi_hsize_t,   &
                         offset_glob,        1, mpi_hsize_t,   &
                         mpi_comm, mpi_ierr)
      call MPI_Allgather(dimsm(axis_write),  1, mpi_hsize_t,   &
                         count_glob,         1, mpi_hsize_t,   &
                         mpi_comm, mpi_ierr)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset, dimsm, hdferror)
    end if


    ! write dataset
{write_string}
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if
{additional_attribute}
    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)
{additional_close}
  end subroutine hdf_write_dataset_{ftype_name}_{rank}
'''

write_regular_dataset_template = '''    if (processor_write == -1) then
      call h5dwrite_f(dset_id, {h5type}, {data_name}, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, {h5type}, {data_name}, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if'''

write_complex_dataset_template = '''    if (processor_write == -1) then
      call h5dwrite_f(dset_id, complexd_field_id(1), real(array),  &
                      dimsf, hdferror,                             &
                      file_space_id=file_space_id,                 &
                      mem_space_id=mem_space_id,                   &
                      xfer_prp=dplist_complex_collective)
      call h5dwrite_f(dset_id, complexd_field_id(2), aimag(array), &
                      dimsf, hdferror,                             &
                      file_space_id=file_space_id,                 &
                      mem_space_id=mem_space_id,                   &
                      xfer_prp=dplist_complex_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, complexd_field_id(1), real(array),  &
                      dimsf, hdferror,                             &
                      file_space_id=file_space_id,                 &
                      mem_space_id=mem_space_id,                   &
                      xfer_prp=dplist_independent)
      call h5dwrite_f(dset_id, complexd_field_id(2), aimag(array), &
                      dimsf, hdferror,                             &
                      file_space_id=file_space_id,                 &
                      mem_space_id=mem_space_id,                   &
                      xfer_prp=dplist_independent)
    end if'''

read_array_template = '''  subroutine hdf_read_dataset_{ftype_name}_{rank}(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
{declaration}
{additional_declaration}
    integer :: rank
{dims_declaration}
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_{ftype_name}_{rank}: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = {rank}
{set_dims}

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
{read_string}

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_{ftype_name}_{rank}
'''

read_regular_dataset_template = '''    call h5dread_f(dset_id, {h5type}, array, dims, hdferror, xfer_prp=dplist_collective)'''

read_complex_dataset_template = '''    call h5dread_f(dset_id, complexd_field_id(1), buffer({buffer_indexing}1), dims, hdferror, xfer_prp=dplist_collective)
    call h5dread_f(dset_id, complexd_field_id(2), buffer({buffer_indexing}2), dims, hdferror, xfer_prp=dplist_collective)
    array = cmplx(buffer({buffer_indexing}1), buffer({buffer_indexing}2), kind=dp)'''