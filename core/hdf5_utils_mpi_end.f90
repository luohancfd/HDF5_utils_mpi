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

    ! write dataset
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, int(array, kind=4), dims, hdferror)
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
  subroutine hdf_read_attr_integer_1_8(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    integer(HSIZE_T), intent(out) :: array(:)   ! data to write to attribute
    integer, allocatable :: array_buffer(:)

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
    allocate(array_buffer(dims(1)))
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, array_buffer, dims, hdferror)
    array = int(array_buffer, kind=HSIZE_T)
    deallocate(array_buffer)
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5aclose_f(attr_id, hdferror)
    !write(*,'(A20,I0)') "h5dclose: ", hdferror
    if (obj_name /= "") then
      call h5oclose_f(obj_id, hdferror)
    end if

  end subroutine hdf_read_attr_integer_1_8

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
