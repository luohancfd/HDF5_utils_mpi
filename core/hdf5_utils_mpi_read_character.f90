  !!---------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_character--------------------------------
  !!---------------------------------------------------------------------------------------

  !  \brief reads a string from an hdf5 file
  subroutine hdf_read_dataset_character_0(loc_id, dset_name, array, offset)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading, useless
    character(len=*), intent(out) :: array      ! data to be written

    integer :: rank
    integer(HSIZE_T),dimension(1) :: dimsf, dimsm, offset_local
    integer(HSIZE_T) :: length
    integer(HID_T) :: dset_id, file_space_id, mem_space_id, dtype_id
    character(len=:), dimension(:), allocatable :: buffer
    integer :: hdferror
    integer :: processor_write, axis_write
    logical :: is_parallel

    rank = 0

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_character_0: "//trim(dset_name)
    end if

    ! get file_space dimension
    ! call hdf_get_dims(loc_id, dset_name, dimsf)

    ! get stacked axis
    axis_write = -1
    is_parallel = .false.
    call hdf_read_attribute(loc_id, dset_name, 'processor', processor_write)
    if (processor_write .eq. -1) then
      call hdf_read_attribute(loc_id, dset_name, 'axis_write', axis_write)
      if (axis_write .ne. -1) then
        is_parallel = .true.
        if (axis_write .ne. 1) then
          write(*, '(A)') "hdf_read_dataset_character_0("//trim(dset_name)// &
                        "): illegal axis_write"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end if
    end if
    dimsm = (/1/)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! open datatype
    call h5dget_type_f(dset_id, dtype_id, hdferror)
    call h5tget_size_f(dtype_id, length, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      mem_space_id = H5S_ALL_F
      file_space_id = H5S_ALL_F
    else
      call h5screate_f(H5S_SCALAR_F, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = mpi_irank
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, dimsm, hdferror)
    end if

    if (length .le. len(array)) then
      call h5dread_f(dset_id, dtype_id, array, dimsm, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    else
      allocate (character(len=length)::buffer(dimsm(1)))
      call h5dread_f(dset_id, dtype_id, buffer, dimsm, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
      array = buffer(1)(1:len(array))
      deallocate(buffer)
    end if

    if (is_parallel) then
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5tclose_f(dtype_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_character_0

  !  \brief reads a 1D array of strings from a hdf5 file
  subroutine hdf_read_dataset_character_1(loc_id, dset_name, array, offset)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    character(len=*), intent(out) :: array(:)   ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(1) :: dimsf, dimsm, offset_local
    integer(HSIZE_T) :: length
    integer(HID_T) :: dset_id, file_space_id, mem_space_id, dtype_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    character(len=:), dimension(:), allocatable :: buffer
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 1

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_character_1: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)

    ! get stacked axis
    axis_write = -1
    is_parallel = .false.
    call hdf_read_attribute(loc_id, dset_name, 'processor', processor_write)
    if (processor_write .eq. -1) then
      call hdf_read_attribute(loc_id, dset_name, 'axis_write', axis_write)
      if (axis_write .ne. -1) then
        is_parallel = .true.
      end if
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_character_1("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_character_1("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_character_1("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_character_1("//trim(dset_name)//&
                    "): illegal offset value"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do

      offset_glob = int(offset, kind=HSIZE_T)
      do ii = 1, mpi_nrank-1
        count_glob(ii) = offset_glob(ii+1) - offset_glob(ii)
      end do
      count_glob(mpi_nrank) = dimsf(axis_write) - offset_glob(mpi_nrank)
    else if (is_parallel) then
      if (mpi_nrank .ne. mpi_nrank_old) then
        write(*, '(A)') "hdf_read_dataset_character_1("//trim(dset_name)// &
                        "): different number of processors, offset needs to be explicitly specified"
        call MPI_Abort(mpi_comm, mpi_ierr)
      else
        call hdf_read_attribute(loc_id, dset_name, 'count',  count_glob)
        call hdf_read_attribute(loc_id, dset_name, 'offset', offset_glob)
      end if
    end if

    dimsm = shape(array, KIND=HSIZE_T)
    if (.not. is_parallel) then
      do ii = 1, rank
        if (dimsm(ii) .ne. dimsf(ii)) then
          write(*, '(A)') "hdf_read_dataset_character_1 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      do ii = 1, rank
        jj = dimsf(ii)
        if (ii == axis_write) jj = count_glob(mpi_irank+1)
        if (dimsm(ii) .ne. jj) then
          write(*, '(A, I2, I8)') "hdf_read_dataset_character_1 ("//trim(dset_name)// &
                          "): array size is wrong", dimsm(ii), jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! open datatype
    call h5dget_type_f(dset_id, dtype_id, hdferror)
    call h5tget_size_f(dtype_id, length, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      mem_space_id = H5S_ALL_F
      file_space_id = H5S_ALL_F
    else
      call h5screate_simple_f(rank, dimsm, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, dimsm, hdferror)
    end if

    if (length .le. len(array)) then
      call h5dread_f(dset_id, dtype_id, array, dimsm, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    else
      allocate (character(len=length)::buffer(dimsm(1)))
      call h5dread_f(dset_id, dtype_id, buffer, dimsm, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
      do ii = 1, dimsm(1)
        array(ii) = buffer(ii)
      end do
      deallocate(buffer)
    end if

    if (is_parallel) then
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5tclose_f(dtype_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_character_1

  !  \brief reads a 2D array of strings from a hdf5 file
  subroutine hdf_read_dataset_character_2(loc_id, dset_name, array, offset)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    character(len=*), intent(out) :: array(:,:) ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(2) :: dimsf, dimsm, offset_local
    integer(HSIZE_T) :: length
    integer(HID_T) :: dset_id, file_space_id, mem_space_id, dtype_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable, dimension(:) :: offset_glob, count_glob
    character(len=:), dimension(:,:), allocatable :: buffer
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 2

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_character_2: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)

    ! get stacked axis
    axis_write = -1
    is_parallel = .false.
    call hdf_read_attribute(loc_id, dset_name, 'processor', processor_write)
    if (processor_write .eq. -1) then
      call hdf_read_attribute(loc_id, dset_name, 'axis_write', axis_write)
      if (axis_write .ne. -1) then
        is_parallel = .true.
      end if
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_character_2("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_character_2("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_character_2("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_character_2("//trim(dset_name)//&
                    "): illegal offset value"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do

      offset_glob = int(offset, kind=HSIZE_T)
      do ii = 1, mpi_nrank-1
        count_glob(ii) = offset_glob(ii+1) - offset_glob(ii)
      end do
      count_glob(mpi_nrank) = dimsf(axis_write) - offset_glob(mpi_nrank)
    else if (is_parallel) then
      if (mpi_nrank .ne. mpi_nrank_old) then
        write(*, '(A)') "hdf_read_dataset_character_2("//trim(dset_name)// &
                        "): different number of processors, offset needs to be explicitly specified"
        call MPI_Abort(mpi_comm, mpi_ierr)
      else
        call hdf_read_attribute(loc_id, dset_name, 'count',  count_glob)
        call hdf_read_attribute(loc_id, dset_name, 'offset', offset_glob)
      end if
    end if

    dimsm = shape(array, KIND=HSIZE_T)
    if (.not. is_parallel) then
      do ii = 1, rank
        if (dimsm(ii) .ne. dimsf(ii)) then
          write(*, '(A)') "hdf_read_dataset_character_2 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      do ii = 1, rank
        jj = dimsf(ii)
        if (ii == axis_write) jj = count_glob(mpi_irank+1)
        if (dimsm(ii) .ne. jj) then
          write(*, '(A, I2, I8)') "hdf_read_dataset_character_2 ("//trim(dset_name)// &
                          "): array size is wrong", dimsm(ii), jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! open datatype
    call h5dget_type_f(dset_id, dtype_id, hdferror)
    call h5tget_size_f(dtype_id, length, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      mem_space_id = H5S_ALL_F
      file_space_id = H5S_ALL_F
    else
      call h5screate_simple_f(rank, dimsm, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, dimsm, hdferror)
    end if

    if (length .le. len(array)) then
      call h5dread_f(dset_id, dtype_id, array, dimsm, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    else
      allocate (character(len=length)::buffer(dimsm(1), dimsm(2)))
      call h5dread_f(dset_id, dtype_id, buffer, dimsm, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
      do ii = 1, dimsm(1)
        do jj = 1, dimsm(2)
          array(ii,jj) = buffer(ii,jj)
        end do
      end do
      deallocate(buffer)
    end if

    if (is_parallel) then
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5tclose_f(dtype_id, hdferror)
    call h5dclose_f(dset_id, hdferror)


  end subroutine hdf_read_dataset_character_2