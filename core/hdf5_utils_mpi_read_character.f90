  !!---------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_character--------------------------------
  !!---------------------------------------------------------------------------------------

  !  \brief reads a string from an hdf5 file
  subroutine hdf_read_dataset_character_0(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    character(len=*), intent(out) :: array      ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(1), length
    integer(HID_T) :: dset_id, dspace_id, dtype_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_character_0: "//trim(dset_name)
    end if

    ! set rank and dims
    dims = (/0/)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dopen: ", hdferror

    ! open dataspace and check dimension
    call h5dget_space_f(dset_id, dspace_id, hdferror)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
    if (rank .ne. 0) then
      write (*, '(A)') "hdf_read_dataset_character_0: inconsistent rank"
      stop
    end if
    call h5sclose_f(dspace_id, hdferror)

    ! open datatype
    call h5dget_type_f(dset_id, dtype_id, hdferror)
    call h5tget_size_f(dtype_id, length, hdferror)

    if (length > len(array)) then
      write (*, '(A, I3, A, I3, A)') "hdf_read_dataset_character_0: buffer string length (", &
        len(array), " < ", length, ")"
      stop
    end if

    ! write dataset
    call h5dread_f(dset_id, dtype_id, array, dims, hdferror, xfer_prp=dplist_collective)
    array = array(1:length)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5tclose_f(dtype_id, hdferror)
    !write(*,'(A20,I0)') "h5tclose: ", hdferror
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_character_0

  !  \brief reads an array of strings from an hdf5 file
  subroutine hdf_read_dataset_character_1(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    character(len=*), intent(out) :: array(:)   ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims_in(1), dims(1), max_dims(1), ii, length
    integer(HID_T) :: dset_id, dspace_id, dtype_id
    character(len=:), dimension(:), allocatable :: buffer
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_character_1: "//trim(dset_name)
    end if

    ! set rank and dims
    dims_in = shape(array)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dopen: ", hdferror

    ! open dataspace and check dimension
    call h5dget_space_f(dset_id, dspace_id, hdferror)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
    if (rank .ne. 1) then
      write (*, '(A)') "hdf_read_dataset_character_1: inconsistent rank"
      stop
    end if
    call h5sget_simple_extent_dims_f(dspace_id, dims, max_dims, hdferror)
    do ii = 1, rank
      if (dims_in(ii) .ne. dims(ii)) then
        write (*, '(A, I2, A, I2, A, I2)') &
          "hdf_read_dataset_character_1: inconsistent dimension(", ii, ") ", &
          dims(ii), "vs", dims_in(ii)
        stop
      end if
    end do
    call h5sclose_f(dspace_id, hdferror)

    ! open data type and check string length
    call h5dget_type_f(dset_id, dtype_id, hdferror)
    call h5tget_size_f(dtype_id, length, hdferror)
    if (length > len(array(1))) then
      write (*, '(A, I3, A, I3, A)') "hdf_read_dataset_character_1: buffer string length (", &
        len(array(1)), " < ", length, ")"
      stop
    end if

    ! read dataset
    if (length == len(array(1))) then
      call h5dread_f(dset_id, dtype_id, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      allocate (character(len=length)::buffer(dims(1)))
      call h5dread_f(dset_id, dtype_id, buffer, dims, hdferror, xfer_prp=dplist_collective)
      do ii = 1, dims(1)
        array(ii) = buffer(ii)
      end do
      deallocate (buffer)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5tclose_f(dtype_id, hdferror)
    !write(*,'(A20,I0)') "h5tclose: ", hdferror
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_character_1

  !  \brief reads an array of strings from an hdf5 file
  subroutine hdf_read_dataset_character_2(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    character(len=*), intent(out) :: array(:, :)   ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims_in(2), dims(2), max_dims(2), ii, jj, length
    integer(HID_T) :: dset_id, dspace_id, dtype_id
    character(len=:), dimension(:, :), allocatable :: buffer
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_character_2: "//trim(dset_name)
    end if

    ! set rank and dims
    dims_in = shape(array)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dopen: ", hdferror

    ! open dataspace and check dimension
    call h5dget_space_f(dset_id, dspace_id, hdferror)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
    if (rank .ne. 2) then
      write (*, '(A)') "hdf_read_dataset_character_2: inconsistent rank"
      stop
    end if
    call h5sget_simple_extent_dims_f(dspace_id, dims, max_dims, hdferror)
    do ii = 1, rank
      if (dims_in(ii) .ne. dims(ii)) then
        write (*, '(A, I2, A, I2, A, I2)') &
          "hdf_read_dataset_character_2: inconsistent dimension(", ii, ") ", &
          dims(ii), "vs", dims_in(ii)
        stop
      end if
    end do
    call h5sclose_f(dspace_id, hdferror)

    ! open data type and check string length
    call h5dget_type_f(dset_id, dtype_id, hdferror)
    call h5tget_size_f(dtype_id, length, hdferror)
    if (length > len(array(1, 1))) then
      write (*, '(A, I3, A, I3, A)') "hdf_read_dataset_character_2: buffer string length (", &
        len(array(1, 1)), " < ", length, ")"
      stop
    end if

    ! read dataset
    if (length == len(array(1, 1))) then
      call h5dread_f(dset_id, dtype_id, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      allocate (character(len=length)::buffer(dims(1), dims(2)))
      call h5dread_f(dset_id, dtype_id, buffer, dims, hdferror, xfer_prp=dplist_collective)
      do ii = 1, dims(1)
        do jj = 1, dims(2)
          array(ii, jj) = buffer(ii, jj)
        end do
      end do
      deallocate (buffer)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5tclose_f(dtype_id, hdferror)
    !write(*,'(A20,I0)') "h5tclose: ", hdferror
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_character_2