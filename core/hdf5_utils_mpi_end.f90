
  !!---------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_integer--------------------------------
  !!---------------------------------------------------------------------------------------

  !  \brief reads a scalar from an hdf5 file
  subroutine hdf_read_dataset_integer_0(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(out) :: array                ! data to be written

    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_0: "//trim(dset_name)
    end if

    ! set rank and dims
    dims = (/0/)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_integer_0

  !  \brief reads a 1d array from an hdf5 file
  subroutine hdf_read_dataset_integer_1(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(out) :: array(:)            ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_integer_1

  !  \brief reads a 2d array from an hdf5 file
  subroutine hdf_read_dataset_integer_2(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(out) :: array(:, :)           ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(2)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_integer_2

  !  \brief reads a 3d array from an hdf5 file
  subroutine hdf_read_dataset_integer_3(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(out) :: array(:, :, :)         ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(3)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_integer_3

  !  \brief reads a 4d array from an hdf5 file
  subroutine hdf_read_dataset_integer_4(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(out) :: array(:, :, :, :)       ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(4)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_integer_4

  !  \brief reads a 5d array from an hdf5 file
  subroutine hdf_read_dataset_integer_5(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(out) :: array(:, :, :, :, :)     ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(5)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_integer_5

  !  \brief reads a 6d array from an hdf5 file
  subroutine hdf_read_dataset_integer_6(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(out) :: array(:, :, :, :, :, :)   ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(6)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 6
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_integer_6

  !!---------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_real--------------------------------
  !!---------------------------------------------------------------------------------------

  !  \brief reads a scalar from an hdf5 file
  subroutine hdf_read_dataset_real_0(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(sp), intent(out) :: array              ! data to be written

    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_0: "//trim(dset_name)
    end if

    ! set rank and dims
    dims = (/0/)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_real_0

  !  \brief reads a 1d array from an hdf5 file
  subroutine hdf_read_dataset_real_1(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(sp), intent(out) :: array(:)            ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_real_1

  !  \brief reads a 2d array from an hdf5 file
  subroutine hdf_read_dataset_real_2(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(sp), intent(out) :: array(:, :)          ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(2)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_real_2

  !  \brief reads a 3d array from an hdf5 file
  subroutine hdf_read_dataset_real_3(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(sp), intent(out) :: array(:, :, :)        ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(3)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_real_3

  !  \brief reads a 4d array from an hdf5 file
  subroutine hdf_read_dataset_real_4(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(sp), intent(out) :: array(:, :, :, :)      ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(4)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_real_4

  !  \brief reads a 5d array from an hdf5 file
  subroutine hdf_read_dataset_real_5(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(sp), intent(out) :: array(:, :, :, :, :)    ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(5)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_real_5

  !  \brief reads a 6d array from an hdf5 file
  subroutine hdf_read_dataset_real_6(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(sp), intent(out) :: array(:, :, :, :, :, :)  ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(6)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 6
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_real_6

  !!---------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_double--------------------------------
  !!---------------------------------------------------------------------------------------

  !  \brief reads a scalar from an hdf5 file
  subroutine hdf_read_dataset_double_0(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(dp), intent(out) :: array               ! data to be written

    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_0: "//trim(dset_name)
    end if

    ! set rank and dims
    dims = (/0/)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_double_0

  !  \brief reads a 1d array from an hdf5 file
  subroutine hdf_read_dataset_double_1(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(dp), intent(out) :: array(:)            ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_double_1

  !  \brief reads a 2d array from an hdf5 file
  subroutine hdf_read_dataset_double_2(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(dp), intent(out) :: array(:, :)          ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(2)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_double_2

  !  \brief reads a 3d array from an hdf5 file
  subroutine hdf_read_dataset_double_3(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(dp), intent(out) :: array(:, :, :)        ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(3)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_double_3

  !  \brief reads a 4d array from an hdf5 file
  subroutine hdf_read_dataset_double_4(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(dp), intent(out) :: array(:, :, :, :)      ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(4)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_double_4

  !  \brief reads a 5d array from an hdf5 file
  subroutine hdf_read_dataset_double_5(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(dp), intent(out) :: array(:, :, :, :, :)    ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(5)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_double_5

  !  \brief reads a 6d array from an hdf5 file
  subroutine hdf_read_dataset_double_6(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    real(dp), intent(out) :: array(:, :, :, :, :, :)  ! data to be written

    integer :: rank
    integer(HSIZE_T) :: dims(6)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 6
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_double_6

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

  !!---------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_double_complex------------------------
  !!---------------------------------------------------------------------------------------

  !  \brief reads a scalar from an hdf5 file
  subroutine hdf_read_dataset_complex_double_0(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    complex(dp), intent(out) :: array           ! data to be read
    real(dp) :: buffer(2)                       ! buffer to save real and imag part

    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_0: "//trim(dset_name)
    end if

    ! set rank and dims
    dims = (/0/)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! read dataset
    call h5dread_f(dset_id, complexd_field_id(1), buffer(1), dims, hdferror, xfer_prp=dplist_collective)
    call h5dread_f(dset_id, complexd_field_id(2), buffer(2), dims, hdferror, xfer_prp=dplist_collective)
    array = cmplx(buffer(1), buffer(2), kind=dp)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_complex_double_0

  !  \brief reads a 1d array from an hdf5 file
  subroutine hdf_read_dataset_complex_double_1(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    complex(dp), intent(out) :: array(:)        ! data to be written
    real(dp), allocatable, dimension(:, :) :: buffer ! buffer to save real and imag part

    integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! read dataset
    allocate (buffer(dims(1), 2))
    call h5dread_f(dset_id, complexd_field_id(1), buffer(:, 1), dims, hdferror, xfer_prp=dplist_collective)
    call h5dread_f(dset_id, complexd_field_id(2), buffer(:, 2), dims, hdferror, xfer_prp=dplist_collective)
    array = cmplx(buffer(:, 1), buffer(:, 2), kind=dp)
    deallocate (buffer)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_complex_double_1

  !  \brief reads a 2d array from an hdf5 file
  subroutine hdf_read_dataset_complex_double_2(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    complex(dp), intent(out) :: array(:, :)        ! data to be written
    real(dp), allocatable, dimension(:, :, :) :: buffer ! buffer to save real and imag part

    integer :: rank
    integer(HSIZE_T) :: dims(2)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! read dataset
    allocate (buffer(dims(1), dims(2), 2))
    call h5dread_f(dset_id, complexd_field_id(1), buffer(:, :, 1), dims, hdferror, xfer_prp=dplist_collective)
    call h5dread_f(dset_id, complexd_field_id(2), buffer(:, :, 2), dims, hdferror, xfer_prp=dplist_collective)
    array = cmplx(buffer(:, :, 1), buffer(:, :, 2), kind=dp)
    deallocate (buffer)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_complex_double_2

  !  \brief reads a 3d array from an hdf5 file
  subroutine hdf_read_dataset_complex_double_3(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    complex(dp), intent(out) :: array(:, :, :)        ! data to be written
    real(dp), allocatable, dimension(:, :, :, :) :: buffer ! buffer to save real and imag part

    integer :: rank
    integer(HSIZE_T) :: dims(3)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! read dataset
    allocate (buffer(dims(1), dims(2), dims(3), 2))
    call h5dread_f(dset_id, complexd_field_id(1), buffer(:, :, :, 1), dims, hdferror, xfer_prp=dplist_collective)
    call h5dread_f(dset_id, complexd_field_id(2), buffer(:, :, :, 2), dims, hdferror, xfer_prp=dplist_collective)
    array = cmplx(buffer(:, :, :, 1), buffer(:, :, :, 2), kind=dp)
    deallocate (buffer)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_complex_double_3

  !  \brief reads a 4d array from an hdf5 file
  subroutine hdf_read_dataset_complex_double_4(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    complex(dp), intent(out) :: array(:, :, :, :)        ! data to be written
    real(dp), allocatable, dimension(:, :, :, :, :) :: buffer ! buffer to save real and imag part

    integer :: rank
    integer(HSIZE_T) :: dims(4)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! read dataset
    allocate (buffer(dims(1), dims(2), dims(3), dims(4), 2))
    call h5dread_f(dset_id, complexd_field_id(1), buffer(:, :, :, :, 1), dims, hdferror, xfer_prp=dplist_collective)
    call h5dread_f(dset_id, complexd_field_id(2), buffer(:, :, :, :, 2), dims, hdferror, xfer_prp=dplist_collective)
    array = cmplx(buffer(:, :, :, :, 1), buffer(:, :, :, :, 2), kind=dp)
    deallocate (buffer)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_complex_double_4

  !  \brief reads a 5d array from an hdf5 file
  subroutine hdf_read_dataset_complex_double_5(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    complex(dp), intent(out) :: array(:, :, :, :, :)        ! data to be written
    real(dp), allocatable, dimension(:, :, :, :, :, :) :: buffer ! buffer to save real and imag part

    integer :: rank
    integer(HSIZE_T) :: dims(5)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! read dataset
    allocate (buffer(dims(1), dims(2), dims(3), dims(4), dims(5), 2))
    call h5dread_f(dset_id, complexd_field_id(1), buffer(:, :, :, :, :, 1), dims, hdferror, xfer_prp=dplist_collective)
    call h5dread_f(dset_id, complexd_field_id(2), buffer(:, :, :, :, :, 2), dims, hdferror, xfer_prp=dplist_collective)
    array = cmplx(buffer(:, :, :, :, :, 1), buffer(:, :, :, :, :, 2), kind=dp)
    deallocate (buffer)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_complex_double_5

  !  \brief reads a 2d array from an hdf5 file
  subroutine hdf_read_dataset_complex_double_6(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    complex(dp), intent(out) :: array(:, :, :, :, :, :)        ! data to be written
    real(dp), allocatable, dimension(:, :, :, :, :, :, :) :: buffer ! buffer to save real and imag part

    integer :: rank
    integer(HSIZE_T) :: dims(6)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
    dims = shape(array, KIND=HSIZE_T)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! read dataset
    allocate (buffer(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), 2))
    call h5dread_f(dset_id, complexd_field_id(1), buffer(:, :, :, :, :, :, 1), dims, hdferror, xfer_prp=dplist_collective)
    call h5dread_f(dset_id, complexd_field_id(2), buffer(:, :, :, :, :, :, 2), dims, hdferror, xfer_prp=dplist_collective)
    array = cmplx(buffer(:, :, :, :, :, :, 1), buffer(:, :, :, :, :, :, 2), kind=dp)
    deallocate (buffer)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5dclose_f(dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror

  end subroutine hdf_read_dataset_complex_double_6

  !!---------------------------------------------------------------------------------------
  !!---------------------------------------hdf_write_attr*---------------------------------
  !!---------------------------------------------------------------------------------------

  !  \brief writes a scalar attribute
  subroutine hdf_write_attr_integer_0(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    integer, intent(in) :: array                 ! data to write to attribute

    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror
    logical :: attr_exists

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_integer_0: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create dataspace
    dims = (/0/)
    call h5screate_f(H5S_SCALAR_F, aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! delete attribute if exists
    call h5aexists_f(obj_id, attr_name, attr_exists, hdferror)
    if (attr_exists) then
      call h5adelete_f(obj_id, attr_name, hdferror)
    end if

    ! create attribute
    call h5acreate_f(obj_id, attr_name, H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    call h5sclose_f(aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5sclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_write_attr_integer_0

  !  \brief writes 1d array attribute
  subroutine hdf_write_attr_integer_1(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    integer, intent(in) :: array(:)              ! data to write to attribute

    integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror
    logical :: attr_exists

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_integer_1: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! delete attribute if exists
    call h5aexists_f(obj_id, attr_name, attr_exists, hdferror)
    if (attr_exists) then
      call h5adelete_f(obj_id, attr_name, hdferror)
    end if

    ! create dataspace
    rank = 1
    dims = shape(array, KIND=HSIZE_T)
    call h5screate_simple_f(rank, dims, aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create attribute
    call h5acreate_f(obj_id, attr_name, H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    call h5sclose_f(aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5sclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_write_attr_integer_1

  subroutine hdf_write_attr_integer_1_8(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    integer(HSIZE_T), intent(in) :: array(:)    ! data to write to attribute

    integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror
    logical :: attr_exists

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_integer_1_8: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! delete attribute if exists
    call h5aexists_f(obj_id, attr_name, attr_exists, hdferror)
    if (attr_exists) then
      call h5adelete_f(obj_id, attr_name, hdferror)
    end if

    ! create dataspace
    rank = 1
    dims = shape(array, KIND=HSIZE_T)
    call h5screate_simple_f(rank, dims, aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create attribute
    call h5acreate_f(obj_id, attr_name, H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset, the data is write as 4 bytes integer
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, int(array), dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    call h5sclose_f(aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5sclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if
  end subroutine hdf_write_attr_integer_1_8

  !  \brief writes a scalar attribute
  subroutine hdf_write_attr_real_0(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    real(sp), intent(in) :: array                ! data to write to attribute

    !integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror
    logical :: attr_exists

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_real_0: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! delete attribute if exists
    call h5aexists_f(obj_id, attr_name, attr_exists, hdferror)
    if (attr_exists) then
      call h5adelete_f(obj_id, attr_name, hdferror)
    end if

    ! create dataspace
    dims = (/0/)
    call h5screate_f(H5S_SCALAR_F, aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create attribute
    call h5acreate_f(obj_id, attr_name, H5T_NATIVE_REAL, aspace_id, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5awrite_f(attr_id, H5T_NATIVE_REAL, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    call h5sclose_f(aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5sclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_write_attr_real_0

  !  \brief writes 1d array attribute
  subroutine hdf_write_attr_real_1(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    real(sp), intent(in) :: array(:)             ! data to write to attribute

    integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror
    logical :: attr_exists

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_real_1: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! delete attribute if exists
    call h5aexists_f(obj_id, attr_name, attr_exists, hdferror)
    if (attr_exists) then
      call h5adelete_f(obj_id, attr_name, hdferror)
    end if

    ! create dataspace
    rank = 1
    dims = shape(array, KIND=HSIZE_T)
    call h5screate_simple_f(rank, dims, aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create attribute
    call h5acreate_f(obj_id, attr_name, H5T_NATIVE_REAL, aspace_id, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5awrite_f(attr_id, H5T_NATIVE_REAL, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    call h5sclose_f(aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5sclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_write_attr_real_1

  !  \brief writes a scalar attribute
  subroutine hdf_write_attr_double_0(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    real(dp), intent(in) :: array                ! data to write to attribute

    !integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror
    logical :: attr_exists

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_double_0: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! delete attribute if exists
    call h5aexists_f(obj_id, attr_name, attr_exists, hdferror)
    if (attr_exists) then
      call h5adelete_f(obj_id, attr_name, hdferror)
    end if

    ! create dataspace
    dims = (/0/)
    call h5screate_f(H5S_SCALAR_F, aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create attribute
    call h5acreate_f(obj_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    call h5sclose_f(aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5sclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_write_attr_double_0

  !  \brief writes 1d array attribute
  subroutine hdf_write_attr_double_1(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    real(dp), intent(in) :: array(:)             ! data to write to attribute

    integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror
    logical :: attr_exists

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_double_1: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! delete attribute if exists
    call h5aexists_f(obj_id, attr_name, attr_exists, hdferror)
    if (attr_exists) then
      call h5adelete_f(obj_id, attr_name, hdferror)
    end if

    ! create dataspace
    rank = 1
    dims = shape(array, KIND=HSIZE_T)
    call h5screate_simple_f(rank, dims, aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create attribute
    call h5acreate_f(obj_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    call h5sclose_f(aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5sclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_write_attr_double_1

  !  \brief writes a string attribute
  subroutine hdf_write_attr_string(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    character(len=*), intent(in) :: array        ! data to write to attribute

    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, type_id, aspace_id, attr_id
    integer :: hdferror
    logical :: attr_exists

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_string: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! delete attribute if exists
    call h5aexists_f(obj_id, attr_name, attr_exists, hdferror)
    if (attr_exists) then
      call h5adelete_f(obj_id, attr_name, hdferror)
    end if

    ! create type_id and aspace_id
    dims(1) = len(array, KIND=HSIZE_T)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, hdferror)
    !write(*,*) 'h5tcopy_f returns', type_id
    call h5tset_size_f(type_id, dims(1), hdferror)
    !write(*,*) 'h5tset_size_f returns', hdferror
    call h5screate_f(H5S_SCALAR_F, aspace_id, hdferror)

    ! create attribute
    call h5acreate_f(obj_id, attr_name, type_id, aspace_id, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5awrite_f(attr_id, type_id, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5tclose_f(type_id, hdferror)
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    call h5sclose_f(aspace_id, hdferror)
    !write(*,'(A20,I0)') "h5sclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_write_attr_string

  !!---------------------------------------------------------------------------------------
  !!-------------------------------------hdf_read_attr-------------------------------------
  !!---------------------------------------------------------------------------------------

  !  \brief writes a scalar attribute
  subroutine hdf_read_attr_integer_0(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    integer, intent(out) :: array                ! data to write to attribute

    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_attr_integer_0: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create attribute
    call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_read_attr_integer_0

  !  \brief writes a scalar attribute
  subroutine hdf_read_attr_integer_1(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    integer, intent(out) :: array(:)             ! data to write to attribute

    integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_attr_integer_1: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create dataspace
    rank = 1
    dims = shape(array, KIND=HSIZE_T)

    ! create attribute
    call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_read_attr_integer_1

  !  \brief writes a scalar attribute
  subroutine hdf_read_attr_real_0(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    real(sp), intent(out) :: array               ! data to write to attribute

    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_attr_real_0: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create attribute
    call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5aread_f(attr_id, H5T_NATIVE_REAL, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_read_attr_real_0

  !  \brief reads 1d array attribute
  subroutine hdf_read_attr_real_1(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    real(sp), intent(out) :: array(:)            ! data to write to attribute

    integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_attr_real_1: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create dataspace
    rank = 1
    dims = shape(array, KIND=HSIZE_T)

    ! create attribute
    call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5aread_f(attr_id, H5T_NATIVE_REAL, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_read_attr_real_1

  !  \brief writes a scalar attribute
  subroutine hdf_read_attr_double_0(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    real(dp), intent(out) :: array               ! data to write to attribute

    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_attr_double_0: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create attribute
    call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_read_attr_double_0

  !  \brief reads 1d array attribute
  subroutine hdf_read_attr_double_1(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    real(dp), intent(out) :: array(:)            ! data to write to attribute

    integer :: rank
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_attr_double_1: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create dataspace
    rank = 1
    dims = shape(array, KIND=HSIZE_T)

    ! create attribute
    call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_read_attr_double_1

  !  \brief writes a string attribute
  subroutine hdf_read_attr_string(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    character(len=*), intent(out) :: array       ! data to write to attribute

    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, type_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_attr_string: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create type_id
    dims(1) = len(array, KIND=HSIZE_T)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, hdferror)
    !write(*,*) 'h5tcopy_f returns', type_id
    call h5tset_size_f(type_id, dims(1), hdferror)

    ! create attribute
    call h5aopen_f(obj_id, attr_name, attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    call h5aread_f(attr_id, type_id, array, dims, hdferror)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_read_attr_string

end module hdf5_utils_mpi
