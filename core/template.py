
write_single_value_template = '''  subroutine hdf_write_dataset_{ftype_name}_0(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
{declaration}
    integer, optional, intent(in) :: chunks           ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: dset_id, dspace_id
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_{ftype_name}_0: "//trim(dset_name)
    end if

    ! set rank and dims
    dims = (/0/)

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

    ! create dataspace
    call h5screate_f(H5S_SCALAR_F, dspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_f: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, {h5type}, dspace_id, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
{write_string}
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_{ftype_name}_0
'''

write_array_template = '''  subroutine hdf_write_dataset_{ftype_name}_{rank}(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
{declaration}
    integer, optional, intent(in) :: chunks({rank})        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims({rank}), cdims({rank})
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_{ftype_name}_{rank}: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = {rank}
    dims = shape(array, KIND=HID_T)

    !
    if (present(filter)) then
      filter_case = filter
    else
      filter_case = hdf_default_filter
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      cdims = int(chunks, SIZE_T)
    else
      cdims = 0
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      processor_write = processor
    end if

    ! create and set property list
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, hdferror)
    call hdf_set_property_list(plist_id, rank, dims, cdims, filter_case)

    ! create dataspace
    call h5screate_simple_f(rank, dims, dspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, {h5type}, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
{write_string}
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_{ftype_name}_{rank}
'''

write_regular_dataset_template = '''    if (processor_write == -1) then
      call h5dwrite_f(dset_id, {h5type}, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, {h5type}, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if'''

write_complex_dataset_template = '''    if (processor_write == -1) then
      call h5dwrite_f(dset_id, complexd_field_id(1), real(array), &
                      dims, hdferror, xfer_prp=dplist_complex_collective)
      call h5dwrite_f(dset_id, complexd_field_id(2), aimag(array), &
                      dims, hdferror, xfer_prp=dplist_complex_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, complexd_field_id(1), real(array), &
                    dims, hdferror, xfer_prp=dplist_independent)
        call h5dwrite_f(dset_id, complexd_field_id(2), aimag(array), &
                    dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if'''


write_char0_template = '''  subroutine hdf_write_dataset_character_0(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    character(len=*), intent(in) :: array             ! data to be written
    integer, optional, intent(in) :: chunks           ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer(SIZE_T) :: dims(1), length
    integer(HID_T) :: dset_id, dspace_id, dtype_id
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_character_0: "//trim(dset_name)
    end if

    ! set rank, dims, dtype
    dims = (/0/)
    length = len(array)
    call h5tcopy_f(H5T_FORTRAN_S1, dtype_id, hdferror)
    call h5tset_size_f(dtype_id, length, hdferror)

    !
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_character_0: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_character_0: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      processor_write = processor
    end if

    ! create dataspace
    call h5screate_f(H5S_SCALAR_F, dspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_f: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, dtype_id, dspace_id, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
{write_string}
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)
    call h5tclose_f(dtype_id, hdferror)

  end subroutine hdf_write_dataset_{ftype_name}_0
'''

write_char_array_template = '''  subroutine hdf_write_dataset_{ftype_name}_{rank}(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
{declaration}
    integer, optional, intent(in) :: chunks({rank})        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims({rank}), cdims({rank}), length
    integer(HID_T) :: dset_id, dspace_id, plist_id, dtype_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_{ftype_name}_{rank}: "//trim(dset_name)
    end if

    ! set rank, dims, dtype
    rank = {rank}
    dims = shape(array, KIND=HID_T)
{length_calculation}
    call h5tcopy_f(H5T_FORTRAN_S1, dtype_id, hdferror)
    call h5tset_size_f(dtype_id, length, hdferror)

    !
    if (present(filter)) then
      filter_case = filter
    else
      filter_case = hdf_default_filter
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      cdims = int(chunks, SIZE_T)
    else
      cdims = 0
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      processor_write = processor
    end if

    ! create and set property list
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, hdferror)
    call hdf_set_property_list(plist_id, rank, dims, cdims, filter_case)

    ! create dataspace
    call h5screate_simple_f(rank, dims, dspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, {h5type}, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
{write_string}
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)
    call h5tclose_f(dtype_id, hdferror)

  end subroutine hdf_write_dataset_{ftype_name}_{rank}
'''