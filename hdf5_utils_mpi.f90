!>  \brief A set of high level wrapper subroutines for HDF5
!>
!>  \par \b Features:
!>   - opening and closing files
!>   - creating/opening/closing groups
!>   - get rank and dimensions of dataset
!>   - reading and writing dataset (integer, real, double)
!>     - uses a generic interface to switch on rank and kind
!>   - writing/reading attributes (integer, real, double, string)
!>     - uses a generic interface to switch on rank and kind
!>
!>  \todo
!>   - reading and writing ( strings )
!>   - hdf_exists  (h5o_exist_by_name or h5l_exists)
!>   - hdf_get_*
!>     - hdf_get_obj_name  (h5iget_name_f)
!>     - hdf_get_obj_type  (h5iget_type_f)
!>     - hdf_get_dset_type   (H5Dget_type)
!>     - hdf_get_obj_id (not needed)
!>   - error checking,
!>     - check dims when reading
!>     - check dataset/attribute name when reading
!>     - check group name when reading/writing
!>     - stop on error vs return error flag vs global error flag
!>
!>  \note I might use H5T_STD_I32LE, H5T_IEEE_F32LE, H5T_IEEE_F64LE instead of H5T_NATIVE_DOUBLE when
!>    creating a dataset. This would make the hdf5 file more portable?
!>  \note should I support other integer types (ie 16 and 64)? I am not sure if hdf5 fortran supports these
!>
module hdf5_utils_mpi

  use hdf5
  implicit none
  include 'mpif.h'

  private

  public :: hdf_open_file, hdf_close_file
  public :: hdf_create_group, hdf_open_group, hdf_close_group
  public :: hdf_exists, hdf_get_rank, hdf_get_dims
  public :: hdf_write_dataset, hdf_read_dataset
  public :: hdf_write_attribute, hdf_read_attribute
  public :: hdf_create_dataset
  public :: hdf_write_vector_to_dataset, hdf_read_vector_from_dataset
  public :: HID_T, hdf_set_print_messages, hdf_set_default_filter

  private :: hdf_preset_prop, hdf_close_prop

  !>  \brief Generic interface to write a dataset
  !>
  !>  Supported types
  !>   - integers (scalar and 1d-6d arrays)
  !>   - doubles (scalar and 1d-6d arrays)
  !>   - reals (scalar and 1d-6d arrays)
  !>   - string (scalar and 1d-2d arrays)
  !>
  !>  \param[in] loc_id     local id in file
  !>  \param[in] dset_name  name of dataset
  !>  \param[in] array      data array to be written
  !>  \param[in] chunks     (optional) chunk size for dataset
  !>  \param[in] filter     (optional) filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
  interface hdf_write_dataset
    module procedure hdf_write_dataset_integer_0
    module procedure hdf_write_dataset_integer_1
    module procedure hdf_write_dataset_integer_2
    module procedure hdf_write_dataset_integer_3
    module procedure hdf_write_dataset_integer_4
    module procedure hdf_write_dataset_integer_5
    module procedure hdf_write_dataset_integer_6
    module procedure hdf_write_dataset_real_0
    module procedure hdf_write_dataset_real_1
    module procedure hdf_write_dataset_real_2
    module procedure hdf_write_dataset_real_3
    module procedure hdf_write_dataset_real_4
    module procedure hdf_write_dataset_real_5
    module procedure hdf_write_dataset_real_6
    module procedure hdf_write_dataset_double_0
    module procedure hdf_write_dataset_double_1
    module procedure hdf_write_dataset_double_2
    module procedure hdf_write_dataset_double_3
    module procedure hdf_write_dataset_double_4
    module procedure hdf_write_dataset_double_5
    module procedure hdf_write_dataset_double_6
    module procedure hdf_write_dataset_character_0
    module procedure hdf_write_dataset_character_1
    module procedure hdf_write_dataset_character_2
    module procedure hdf_write_dataset_complex_double_0
    module procedure hdf_write_dataset_complex_double_1
    module procedure hdf_write_dataset_complex_double_2
    module procedure hdf_write_dataset_complex_double_3
    module procedure hdf_write_dataset_complex_double_4
    module procedure hdf_write_dataset_complex_double_5
    module procedure hdf_write_dataset_complex_double_6
  end interface hdf_write_dataset

  !>  \brief Generic interface to write a vector to dataset
  !>
  !>  The vector is is written out in along the fast dimension
  !>  (column oriented in FORTRAN, row oriented in the HDF5 file).
  !>  So the vector should have the same length as the first dimension
  !>  of the dataset, and the offset should agree with dims(2:rank) of
  !>  of the dataset.
  !>
  !>  Supported types
  !>   - integers (scalar and 1d-6d arrays)
  !>   - doubles (scalar and 1d-6d arrays)
  !>   - reals (scalar and 1d-6d arrays)
  !   - string (scalar and 1d-6d arrays)
  !>
  !>  \param[in] loc_d      local id in file
  !>  \param[in] dset_name  name of dataset
  !>  \param[in] offset     position within the dataset
  !>  \param[in] vector     data array to be written
  interface hdf_write_vector_to_dataset
    module procedure hdf_write_vector_to_dataset_integer
    module procedure hdf_write_vector_to_dataset_real
    module procedure hdf_write_vector_to_dataset_double
  end interface hdf_write_vector_to_dataset

  !>  \brief Generic interface to read a dataset of doubles
  !>
  !>  Supported types
  !>   - integers (scalar and 1d-6d arrays)
  !>   - doubles (scalar and 1d-6d arrays)
  !>   - reals (scalar and 1d-6d arrays)
  !>   - string (scalar and 1d-2d arrays)
  !>
  !>  \param[in]  loc_d      local id in file
  !>  \param[in]  dset_name  name of dataset
  !>  \param[out] array      data array to be read
  interface hdf_read_dataset
    module procedure hdf_read_dataset_integer_0
    module procedure hdf_read_dataset_integer_1
    module procedure hdf_read_dataset_integer_2
    module procedure hdf_read_dataset_integer_3
    module procedure hdf_read_dataset_integer_4
    module procedure hdf_read_dataset_integer_5
    module procedure hdf_read_dataset_integer_6
    module procedure hdf_read_dataset_real_0
    module procedure hdf_read_dataset_real_1
    module procedure hdf_read_dataset_real_2
    module procedure hdf_read_dataset_real_3
    module procedure hdf_read_dataset_real_4
    module procedure hdf_read_dataset_real_5
    module procedure hdf_read_dataset_real_6
    module procedure hdf_read_dataset_double_0
    module procedure hdf_read_dataset_double_1
    module procedure hdf_read_dataset_double_2
    module procedure hdf_read_dataset_double_3
    module procedure hdf_read_dataset_double_4
    module procedure hdf_read_dataset_double_5
    module procedure hdf_read_dataset_double_6
    module procedure hdf_read_dataset_character_0
    module procedure hdf_read_dataset_character_1
    module procedure hdf_read_dataset_character_2
    module procedure hdf_read_dataset_complex_double_0
    module procedure hdf_read_dataset_complex_double_1
    module procedure hdf_read_dataset_complex_double_2
    module procedure hdf_read_dataset_complex_double_3
    module procedure hdf_read_dataset_complex_double_4
    module procedure hdf_read_dataset_complex_double_5
    module procedure hdf_read_dataset_complex_double_6
  end interface hdf_read_dataset

  !>  \brief Generic interface to read a vector from a dataset
  !>
  !>  The vector is is read in along the fast dimension
  !>  (column oriented in FORTRAN, row oriented in the HDF5 file).
  !>  So the vector should have the same length as the first dimension
  !>  of the dataset, and the offset should agree with dims(2:rank) of
  !>  of the dataset.
  !>
  !>  Supported types
  !>   - integers (scalar and 1d-6d arrays)
  !>   - doubles (scalar and 1d-6d arrays)
  !>   - reals (scalar and 1d-6d arrays)
  !   - string (scalar and 1d-6d arrays)
  !>
  !>  \param[in] loc_d      local id in file
  !>  \param[in] dset_name  name of dataset
  !>  \param[in] offset     position within the dataset
  !>  \param[out] vector    data array to be written
  interface hdf_read_vector_from_dataset
    module procedure hdf_read_vector_from_dataset_integer
    module procedure hdf_read_vector_from_dataset_real
    module procedure hdf_read_vector_from_dataset_double
  end interface hdf_read_vector_from_dataset

  !>  \brief Generic interface to write attribute
  !>
  !>  Supported types
  !>   - integers (scalar and 1d arrays)
  !>   - reals (scalar and 1d arrays)
  !>   - doubles (scalar and 1d arrays)
  !>   - string (1d array of characters)
  !>
  !>  \param[in] loc_id     local id in file
  !>  \param[in] obj_name   name of object to be attached to (if left blank, just use loc_id)
  !>  \param[in] attr_name  name of attribute to be added
  !>  \param[in] array      attribute data to be written
  interface hdf_write_attribute
    module procedure hdf_write_attr_integer_0
    module procedure hdf_write_attr_integer_1
    module procedure hdf_write_attr_real_0
    module procedure hdf_write_attr_real_1
    module procedure hdf_write_attr_double_0
    module procedure hdf_write_attr_double_1
    module procedure hdf_write_attr_string
  end interface hdf_write_attribute

  !>  \brief Generic interface to read attribute
  !>
  !>  Supported types
  !>   - integers (scalar and 1d arrays)
  !>   - reals (scalar and 1d arrays)
  !>   - doubles (scalar and 1d arrays)
  !>   - string (1d array of characters)
  !>
  !>  \param[in] loc_id     local id in file
  !>  \param[in] obj_name   name of object to be attached to (if left blank, just use loc_id)
  !>  \param[in] attr_name  name of attribute to be added
  !>  \param[in] array      attribute data to be written
  interface hdf_read_attribute
    module procedure hdf_read_attr_integer_0
    module procedure hdf_read_attr_integer_1
    module procedure hdf_read_attr_real_0
    module procedure hdf_read_attr_real_1
    module procedure hdf_read_attr_double_0
    module procedure hdf_read_attr_double_1
    module procedure hdf_read_attr_string
  end interface hdf_read_attribute

  ! precision for this file
  integer, parameter :: sp = kind(1.0)     ! single precision
  integer, parameter :: dp = kind(1.0d0)   ! double precision

  !
  logical :: hdf_print_messages = .false.

  !
  character(len=32) :: hdf_default_filter = 'none'
  integer :: hdf_gzip_level = 6  ! 0-9
  integer :: hdf_szip_pixels_per_block = 8    ! should be even number less than 32 (https://portal.hdfgroup.org/display/HDF5/H5P_SET_SZIP)
  character(len=2) :: hdf_szip_options = 'NN'  ! H5_SZIP_NN_OM_F or H5_SZIP_EC_OM_F

  !
  ! compound type for double precision complex number
  integer(HID_T) :: complexd_type_id, complexd_field_id(2)
  character(len=1), parameter :: complexd_field_name(2) = (/'r', 'i'/)
  integer(HID_T) :: complexd_pid
  integer(HID_T) :: dplist_complex_collective, dplist_complex_independent

  !
  ! property list for parallel API
  integer :: mpi_comm, mpi_irank, mpi_nrank, mpi_ierr
  integer(HID_T) :: file_plist_id   !< parallel file access property
  integer(HID_T) :: dplist_collective, dplist_independent !< dataset access property


contains

  !>  \brief Sets the value of hdf_print_messages
  !>
  !>  By default, hdf_print_messages = .false. By setting it
  !>  to .true., some messages are printed detailing what hdf_utils
  !>  is doing.
  subroutine hdf_set_print_messages(val_print_messages)

    logical, intent(in) :: val_print_messages  !<  new value for hdf_print_messages

    hdf_print_messages = val_print_messages

  end subroutine hdf_set_print_messages

  !>  \brief Sets the value of hdf_print_messages
  !>
  !>
  subroutine hdf_set_default_filter(filter, gzip_level, szip_pixels_per_block, szip_options)

    character(len=*), intent(in) :: filter  !<  new value for hdf_print_messages
    integer, optional, intent(in) :: gzip_level
    integer, optional, intent(in) :: szip_pixels_per_block
    character(len=2), optional, intent(in) :: szip_options

    hdf_default_filter = filter

    if (present(gzip_level)) then
      hdf_gzip_level = gzip_level
    end if

    if (present(szip_pixels_per_block)) then
      hdf_szip_pixels_per_block = szip_pixels_per_block
    end if

    if (present(szip_options)) then
      !hdf_szip_options = szip_options
      if (szip_options == 'NN') then
        hdf_szip_options = 'NN'
      elseif (szip_options == 'EC') then
        hdf_szip_options = 'EC'
      else
        write (*, '(A)') "-->hdf_set_default_filter: warning szip_options "//trim(szip_options)//" not recognized"
        stop
      end if
    end if
    !write(*,'(A,I0,A,I0)') ' H5_SZIP_NN_OM_F=', H5_SZIP_NN_OM_F, ', H5_SZIP_EC_OM_F=', H5_SZIP_EC_OM_F

    write (*, *) filter, gzip_level, szip_pixels_per_block, szip_options

  end subroutine hdf_set_default_filter

  !>  \brief Check if location exists.
  !>
  !>  Also checks is intemediate paths exists in a safe way.
  subroutine hdf_exists(loc_id, obj_name, exists)

    integer(HID_T), intent(in) :: loc_id       !< local id
    character(len=*), intent(in) :: obj_name   !< relative path to object
    logical, intent(out) :: exists  !< .TRUE. if everything exists, .FALSE. otherwise

    integer :: hdferror, pos, cpos, str_len

    ! check intermediate paths (subgroups)
    str_len = len_trim(obj_name)
    cpos = 0
    exists = .false.
    do
      !start = cpos + 1
      !write(*,*) start, str_len, obj_name(start:str_len)

      pos = index(obj_name(cpos + 1:str_len), "/")

      ! no subgroup found
      if (pos == 0) exit

      ! check subgroup
      cpos = cpos + pos
      call h5lexists_f(loc_id, obj_name(1:cpos - 1), exists, hdferror)
      ! write(*,*) obj_name(1:cpos-1), exists

      ! return if intermediate path fails
      if (exists .eqv. .false.) then
        if (hdf_print_messages) then
          write (*, '(A,A,A)') "--->hdf_exists: subpath '", trim(obj_name(1:cpos - 1)), "' does not exist, return false"
        end if
        exists = .false.
        return
      end if

    end do

    ! check object (unless obj_name ended with "/"
    if (cpos /= str_len) then
      call h5lexists_f(loc_id, obj_name, exists, hdferror)
      ! write(*,*) obj_name, exists
      if (exists .eqv. .false.) then
        if (hdf_print_messages) then
          write (*, '(A,A,A)') "--->hdf_exists: object '", trim(obj_name), "' does not exist, return false"
        end if
        exists = .false.
        return
      end if
    end if

    exists = .true.
    return

  end subroutine hdf_exists

  !>  \brief Opens file and return identifier
  !>
  !>  \todo
  !>   - case insentive STATUS and ACTION
  !>   - delete file for REPLACE case
  !>
  !>  | STATUS  | ACTION    | Description                          |
  !>  | :-----: | :-------: | :----------------------------------- |
  !>  | NEW     | na        | calls h5fcreate with H5F_ACC_TRUNC_F |
  !>  | REPLACE | na        | calls h5fcreate with H5F_ACC_EXCL_F  |
  !>  | OLD     | READ      | calls h5fopen with H5F_ACC_RDONLY_F  |
  !>  | OLD     | WRITE     | calls h5fopen with H5F_ACC_RDWR_F    |
  !>  | OLD     | READWRITE | calls h5fopen with H5F_ACC_RDWR_F    |
  !>
  subroutine hdf_open_file(file_id, filename, STATUS, ACTION)

    integer(HID_T), intent(out) :: file_id            !< HDF5 id of the file
    character(len=*), intent(in) :: filename          !< the HDF5 filename
    character(len=*), optional, intent(in) :: STATUS  !< file status (OLD, NEW, REPLACE)
    character(len=*), optional, intent(in) :: ACTION  !< file action (READ, WRITE, READWRITE)

    integer :: hdferror
    character(len=16) :: status2, action2

    if (hdf_print_messages) then
      write (*, '(A)') "->hdf_open_file: "//trim(filename)
    end if

    ! open hdf5 interface
    call h5open_f(hdferror)
    !write(*,'(A20,I0)') "h5open: ", hdferror

    ! set defaults
    status2 = 'NEW'
    if (present(STATUS)) status2 = STATUS
    action2 = 'READWRITE'
    if (present(STATUS)) action2 = ACTION

    ! set mpi related
    mpi_comm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_nrank, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_irank, mpi_ierr)

    ! set up file access property list with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, file_plist_id, hdferror)
    call h5pset_fapl_mpio_f(file_plist_id, mpi_comm, MPI_INFO_NULL, hdferror)

    ! open/create hdf5 file
    if (status2 == 'OLD') then
      if (action2 == 'READ') then
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferror, access_prp=file_plist_id)
      elseif ((action2 == 'WRITE') .or. (action2 == 'READWRITE')) then
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdferror, access_prp=file_plist_id)
      else
        write (*, *) "hdf_open: action = ", action2, " not supported."
        stop
      end if
    elseif (status2 == 'NEW') then
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdferror, access_prp=file_plist_id)
    elseif (status2 == 'REPLACE') then
      call MPI_FILE_DELETE(filename, MPI_INFO_NULL, mpi_ierr)
      call h5fcreate_f(filename, H5F_ACC_EXCL_F, file_id, hdferror, access_prp=file_plist_id)
    else
      write (*, *) "hdf_open: status = ", status2, " not supported."
      stop
    end if

    if (hdferror .ne. 0) then
      write (*, *) "Fail to open HDF5 file"
      stop
    end if

    call hdf_preset_prop()

    !write(*,'(A20,I0)') "h5fcreate: ", hdferror

  end subroutine hdf_open_file

  !>  \brief Closes a hdf5 file
  subroutine hdf_close_file(file_id)

    integer(HID_T), intent(in) :: file_id  !< file id to be closed

    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "->hdf_close_file"
    end if

    call hdf_close_prop()

    call h5pclose_f(file_plist_id, hdferror)
    call h5fclose_f(file_id, hdferror)
    !write(*,'(A20,I0)') "h5fclose: ", hdferror

    call h5close_f(hdferror)

  end subroutine hdf_close_file

  !>  \brief Preset some properties
  subroutine hdf_preset_prop()
    implicit none
    integer :: ii
    integer(SIZE_T) ::offset
    INTEGER(SIZE_T)     ::   type_sizei  ! Size of the integer datatype
    INTEGER(SIZE_T)     ::   type_sized  ! Size of the double precision datatype
    INTEGER(SIZE_T)     ::   type_sizer  ! Size of the real datatype
    integer :: hdferror

    integer(HID_T) :: temp_int(3)

    call h5tget_size_f(H5T_NATIVE_INTEGER, type_sizei, hdferror)
    call h5tget_size_f(H5T_NATIVE_DOUBLE, type_sized, hdferror)
    call h5tget_size_f(H5T_NATIVE_REAL, type_sizer, hdferror)

    ! set dataset transfer property to preserve partially initialized fields
    ! during write/read to/from dataset with compound datatype.
    temp_int = (/complexd_pid, dplist_complex_collective, dplist_complex_collective/)
    do ii = 1, 3
      call h5pcreate_f(H5P_DATASET_XFER_F, temp_int(ii), hdferror)
      call h5pset_preserve_f(temp_int(ii), .TRUE., hdferror)
    end do

    !
    ! create compound type for double complex
    !
    call h5tcreate_f(H5T_COMPOUND_F, type_sized*2, complexd_type_id, hdferror)
    offset = 0
    do ii = 1, 2
      call h5tinsert_f(complexd_type_id, complexd_field_name(ii), offset, &
                       H5T_NATIVE_DOUBLE, hdferror)
      offset = offset + type_sized
    end do
    offset = 0
    do ii = 1, 2
      call h5tcreate_f(H5T_COMPOUND_F, type_sized, complexd_field_id(ii), hdferror)
      call h5tinsert_f(complexd_field_id(ii), complexd_field_name(ii), offset, &
                       H5T_NATIVE_DOUBLE, hdferror)
    end do

    !
    ! set dataset access property
    call h5pcreate_f(H5P_DATASET_XFER_F, dplist_collective, hdferror)
    call h5pset_dxpl_mpio_f(dplist_collective, H5FD_MPIO_COLLECTIVE_F, hdferror)
    call h5pcreate_f(H5P_DATASET_XFER_F, dplist_independent, hdferror)
    call h5pset_dxpl_mpio_f(dplist_independent, H5FD_MPIO_INDEPENDENT_F, hdferror)
    call h5pcreate_f(H5P_DATASET_XFER_F, dplist_complex_collective, hdferror)
    call h5pset_dxpl_mpio_f(dplist_complex_collective, H5FD_MPIO_COLLECTIVE_F, hdferror)
    call h5pcreate_f(H5P_DATASET_XFER_F, dplist_complex_independent, hdferror)
    call h5pset_dxpl_mpio_f(dplist_complex_independent, H5FD_MPIO_INDEPENDENT_F, hdferror)

  end subroutine hdf_preset_prop

  !>  \brief Preset some properties
  subroutine hdf_close_prop()
    implicit none
    integer :: ii
    integer :: hdferror

    call h5pclose_f(complexd_pid, hdferror)

    call h5tclose_f(complexd_type_id, hdferror)
    do ii = 1, 2
      call h5tclose_f(complexd_field_id(ii), hdferror)
    end do

    call h5pclose_f(dplist_complex_collective, hdferror)
    call h5pclose_f(dplist_complex_independent, hdferror)
    call h5pclose_f(dplist_collective, hdferror)
    call h5pclose_f(dplist_independent, hdferror)

  end subroutine hdf_close_prop

  !>  \brief Create a new group
  subroutine hdf_create_group(loc_id, group_name)

    integer(HID_T), intent(in) :: loc_id         !< location id where to put the group
    character(len=*), intent(in) :: group_name   !< name of the group

    integer(HID_T) :: grp_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "->hdf_create_group: "//trim(group_name)
    end if

    call h5gcreate_f(loc_id, group_name, grp_id, hdferror)
    !write(*,'(A20,I0)') "h5gcreate: ", hdferror

    call h5gclose_f(grp_id, hdferror)
    !write(*,'(A20,I0)') "h5gclose: ", hdferror

  end subroutine hdf_create_group

  !>  \brief Opens a group and returns the identifier
  subroutine hdf_open_group(loc_id, group_name, group_id)

    integer(HID_T), intent(in) :: loc_id         !< location id where to put the group
    character(len=*), intent(in) :: group_name   !< name of the group
    integer(HID_T), intent(out) :: group_id      !< id for the group
    logical :: exist

    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A,A,A)') "->hdf_open_group: '"//trim(group_name)//"'"
    end if

    call hdf_exists(loc_id, group_name, exist)
    if (exist) then
      if (hdf_print_messages) then
        write (*, '(A,A,A)') "->hdf_open_group: opening group '"//trim(group_name)//"'"
      end if
      call h5gopen_f(loc_id, group_name, group_id, hdferror)
    else
      if (hdf_print_messages) then
        write (*, '(A,A,A)') "->hdf_open_group: group '"//trim(group_name)//"' does not exist, return with error"
      end if
      hdferror = -1
    end if

  end subroutine hdf_open_group

  !>  \brief Close a group by identifier
  subroutine hdf_close_group(group_id)

    integer(HID_T), intent(in) :: group_id   !< id for the group

    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "->hdf_close_group"
    end if

    call h5gclose_f(group_id, hdferror)
    !write(*,'(A20,I0)') "h5gclose: ", hdferror

  end subroutine hdf_close_group

  !>  \brief Get the rank of a dataset
  subroutine hdf_get_rank(loc_id, dset_name, rank)

    integer(HID_T), intent(in) :: loc_id        !< location id
    character(len=*), intent(in) :: dset_name   !< dataset name
    integer, intent(out) :: rank                !< rank of the dataset

    integer(HID_T) :: dset_id, dspace_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "->hdf_get_rank"
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! get dataspace
    call h5dget_space_f(dset_id, dspace_id, hdferror)

    ! get rank (ndims)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)

    ! close id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_get_rank

  !>  \brief get the dimensions of a dataset
  subroutine hdf_get_dims(loc_id, dset_name, dims)

    integer(HID_T), intent(in) :: loc_id        !< location id
    character(len=*), intent(in) :: dset_name   !< name of dataset
    integer, intent(out) :: dims(:)             !< dimensions of the dataset

    integer(HID_T) :: dset_id, dspace_id
    integer :: rank
    integer(HSIZE_T) :: dset_dims(6), max_dims(6)
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "->hdf_get_dims"
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! get dataspace
    call h5dget_space_f(dset_id, dspace_id, hdferror)

    ! get rank (ndims)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)

    ! get dims
    call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)
    dims(1:rank) = int(dset_dims(1:rank))

    ! close id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_get_dims

  !     - hdf_get_kind   (H5Dget_type)

  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_vector_to_dataset-----------------------------
  !!----------------------------------------------------------------------------------------

  !>  \brief create a dataset 'dset_name' with shape 'dset_dims' of type 'dset_type'
  subroutine hdf_create_dataset(loc_id, dset_name, dset_dims, dset_type)

    integer(HID_T), intent(in) :: loc_id        !< local id in file
    character(len=*), intent(in) :: dset_name   !< name of dataset
    integer, intent(in) :: dset_dims(:)         !< dimensions of the dataset
    character(len=*), intent(in) :: dset_type   !< type of dataset (integer or double)

    integer :: rank
    integer(SIZE_T) :: dims(8)
    integer(HID_T) :: dset_id, dspace_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_create_dataset: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = size(dset_dims, 1)
    dims(1:rank) = int(dset_dims, SIZE_T)

    ! create dataspace
    call h5screate_simple_f(rank, dims, dspace_id, hdferror)
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    select case (dset_type)
    case ('integer')
      call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror)
      call h5dclose_f(dset_id, hdferror)
    case ('real')
      call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror)
      call h5dclose_f(dset_id, hdferror)
    case ('double')
      call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror)
      call h5dclose_f(dset_id, hdferror)
    case default
      write (*, '(A,A,A)') "---> ERROR: dset_type ", dset_type, " not supported"
    end select

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)

  end subroutine hdf_create_dataset

  !
  subroutine hdf_write_vector_to_dataset_integer(loc_id, dset_name, offset, vector)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(in) :: offset(:)            ! position within dataset
    integer, intent(out) :: vector(:)           ! data to be written

    integer(HID_T) :: dset_id, dspace_id, mspace_id
    integer :: rank
    integer(HSIZE_T) :: dset_dims(6), max_dims(6), mdims(1)
    integer(HSIZE_T) :: hs_count(6), hs_offset(6) ! Hyperslab offset
    integer :: hdferror

    integer :: i
    character(len=32) :: format_string

    vector = 0
    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_vector_to_dataset_integer: "//trim(dset_name)
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! get dataspace
    call h5dget_space_f(dset_id, dspace_id, hdferror)

    ! get rank (ndims), and dims
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
    call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)

    ! check size and offset
    if (size(vector, 1) == dset_dims(1)) then

      if (all(offset(1:rank - 1) <= dset_dims(2:rank)) .and. all(offset(1:rank - 1) > 0)) then

        ! select hyperslab in dataset (note: convert FORTRAN offset to C offset by subtracting 1)
        hs_count(1) = dset_dims(1)
        hs_count(2:rank) = 1
        hs_offset(1) = 0
        hs_offset(2:rank) = offset(1:rank - 1) - 1
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_offset(1:rank), hs_count(1:rank), hdferror)

        ! set mspace to a vector
        mdims(1) = size(vector, 1)
        call h5screate_simple_f(1, mdims, mspace_id, hdferror)

        ! write out vector
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, vector, mdims, hdferror, mspace_id, dspace_id)

        ! close mspace_id
        call h5sclose_f(mspace_id, hdferror)

      else
        write (format_string, '(A,I0,A,I0,A)') '(A,', rank - 1, '(I0,A),A,', rank - 1, '(I0,A),A)'
        write (*, format_string) "--->ERROR: offset=(", (offset(i), ',', i=1, rank - 1), &
          "), is not constent with dset_dims(2:rank)=(", (dset_dims(i), ',', i=2, rank), ")"
      end if

    else
      write (*, '(A,I0,A,I0)') "--->ERROR: size(vector)=", size(vector), &
        ", is not constent with dset_dims(1)=", dset_dims(1)
    end if

    ! close id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_vector_to_dataset_integer

  !
  subroutine hdf_write_vector_to_dataset_real(loc_id, dset_name, offset, vector)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(in) :: offset(:)            ! position within dataset
    real(sp), intent(in) :: vector(:)           ! data to be written

    integer(HID_T) :: dset_id, dspace_id, mspace_id
    integer :: rank
    integer(HSIZE_T) :: dset_dims(6), max_dims(6), mdims(1)
    integer(HSIZE_T) :: hs_count(6), hs_offset(6)
    integer :: hdferror

    integer :: i
    character(len=32) :: format_string

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_vector_to_dataset_real: "//trim(dset_name)
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! get dataspace
    call h5dget_space_f(dset_id, dspace_id, hdferror)

    ! get rank (ndims), and dims
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
    call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)

    ! check size and offset
    if (size(vector, 1) == dset_dims(1)) then

      if (all(offset(1:rank - 1) <= dset_dims(2:rank)) .and. all(offset(1:rank - 1) > 0)) then

        ! select hyperslab in dataset (note: convert FORTRAN offset to C offset by subtracting 1)
        hs_count(1) = dset_dims(1)
        hs_count(2:rank) = 1
        hs_offset(1) = 0
        hs_offset(2:rank) = offset(1:rank - 1) - 1
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_offset(1:rank), hs_count(1:rank), hdferror)

        ! set mspace to a vector
        mdims(1) = size(vector, 1)
        call h5screate_simple_f(1, mdims, mspace_id, hdferror)

        ! write out vector
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, vector, mdims, hdferror, mspace_id, dspace_id)

        ! close mspace_id
        call h5sclose_f(mspace_id, hdferror)

      else
        write (format_string, '(A,I0,A,I0,A)') '(A,', rank - 1, '(I0,A),A,', rank - 1, '(I0,A),A)'
        write (*, format_string) "--->ERROR: offset=(", (offset(i), ',', i=1, rank - 1), &
          "), is not constent with dset_dims(2:rank)=(", (dset_dims(i), ',', i=2, rank), ")"
      end if

    else
      write (*, '(A,I0,A,I0)') "--->ERROR: size(vector)=", size(vector), &
        ", is not constent with dset_dims(1)=", dset_dims(1)
    end if

    ! close id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_vector_to_dataset_real

  !
  subroutine hdf_write_vector_to_dataset_double(loc_id, dset_name, offset, vector)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(in) :: offset(:)            ! position within dataset
    real(dp), intent(in) :: vector(:)           ! data to be written

    integer(HID_T) :: dset_id, dspace_id, mspace_id
    integer :: rank
    integer(HSIZE_T) :: dset_dims(6), max_dims(6), mdims(1)
    integer(HSIZE_T) :: hs_count(6), hs_offset(6)
    integer :: hdferror

    integer :: i
    character(len=32) :: format_string

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_vector_to_dataset_double: "//trim(dset_name)
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! get dataspace
    call h5dget_space_f(dset_id, dspace_id, hdferror)

    ! get rank (ndims), and dims
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
    call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)

    ! check size and offset
    if (size(vector, 1) == dset_dims(1)) then

      if (all(offset(1:rank - 1) <= dset_dims(2:rank)) .and. all(offset(1:rank - 1) > 0)) then

        ! select hyperslab in dataset (note: convert FORTRAN offset to C offset by subtracting 1)
        hs_count(1) = dset_dims(1)
        hs_count(2:rank) = 1
        hs_offset(1) = 0
        hs_offset(2:rank) = offset(1:rank - 1) - 1
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_offset(1:rank), hs_count(1:rank), hdferror)

        ! set mspace to a vector
        mdims(1) = size(vector, 1)
        call h5screate_simple_f(1, mdims, mspace_id, hdferror)

        ! write out vector
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vector, mdims, hdferror, mspace_id, dspace_id)

        ! close mspace_id
        call h5sclose_f(mspace_id, hdferror)

      else
        write (format_string, '(A,I0,A,I0,A)') '(A,', rank - 1, '(I0,A),A,', rank - 1, '(I0,A),A)'
        write (*, format_string) "--->ERROR: offset=(", (offset(i), ',', i=1, rank - 1), &
          "), is not constent with dset_dims(2:rank)=(", (dset_dims(i), ',', i=2, rank), ")"
      end if

    else
      write (*, '(A,I0,A,I0)') "--->ERROR: size(vector)=", size(vector), &
        ", is not constent with dset_dims(1)=", dset_dims(1)
    end if

    ! close id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_vector_to_dataset_double

  !
  subroutine hdf_read_vector_from_dataset_integer(loc_id, dset_name, offset, vector)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(in) :: offset(:)            ! position within dataset
    integer, intent(out) :: vector(:)           ! data to be written

    integer(HID_T) :: dset_id, dspace_id, mspace_id
    integer :: rank
    integer(HSIZE_T) :: dset_dims(6), max_dims(6), mdims(1)
    integer(HSIZE_T) :: hs_count(6), hs_offset(6) ! Hyperslab offset
    integer :: hdferror

    integer :: i
    character(len=32) :: format_string

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_vector_from_dataset_integer: "//trim(dset_name)
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! get dataspace
    call h5dget_space_f(dset_id, dspace_id, hdferror)

    ! get rank (ndims), and dims
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
    call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)

    ! check size and offset
    if (size(vector, 1) == dset_dims(1)) then

      if (all(offset(1:rank - 1) <= dset_dims(2:rank)) .and. all(offset(1:rank - 1) > 0)) then

        ! select hyperslab in dataset (note: convert FORTRAN offset to C offset by subtracting 1)
        hs_count(1) = dset_dims(1)
        hs_count(2:rank) = 1
        hs_offset(1) = 0
        hs_offset(2:rank) = offset(1:rank - 1) - 1
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_offset(1:rank), hs_count(1:rank), hdferror)

        ! set mspace to a vector
        mdims(1) = size(vector, 1)
        call h5screate_simple_f(1, mdims, mspace_id, hdferror)

        ! write out vector
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, vector, mdims, hdferror, mspace_id, dspace_id)

        ! close mspace_id
        call h5sclose_f(mspace_id, hdferror)

      else
        write (format_string, '(A,I0,A,I0,A)') '(A,', rank - 1, '(I0,A),A,', rank - 1, '(I0,A),A)'
        write (*, format_string) "--->ERROR: offset=(", (offset(i), ',', i=1, rank - 1), &
          "), is not constent with dset_dims(2:rank)=(", (dset_dims(i), ',', i=2, rank), ")"
      end if

    else
      write (*, '(A,I0,A,I0)') "--->ERROR: size(vector)=", size(vector), &
        ", is not constent with dset_dims(1)=", dset_dims(1)
    end if

    ! close id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_vector_from_dataset_integer

  !
  subroutine hdf_read_vector_from_dataset_real(loc_id, dset_name, offset, vector)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(in) :: offset(:)            ! position within dataset
    real(sp), intent(out) :: vector(:)          ! data to be written

    integer(HID_T) :: dset_id, dspace_id, mspace_id
    integer :: rank
    integer(HSIZE_T) :: dset_dims(6), max_dims(6), mdims(1)
    integer(HSIZE_T) :: hs_count(6), hs_offset(6)
    integer :: hdferror

    integer :: i
    character(len=32) :: format_string

    vector = 0.0
    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_vector_from_dataset_real: "//trim(dset_name)
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! get dataspace
    call h5dget_space_f(dset_id, dspace_id, hdferror)

    ! get rank (ndims), and dims
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
    call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)

    ! check size and offset
    if (size(vector, 1) == dset_dims(1)) then

      if (all(offset(1:rank - 1) <= dset_dims(2:rank)) .and. all(offset(1:rank - 1) > 0)) then

        ! select hyperslab in dataset (note: convert FORTRAN offset to C offset by subtracting 1)
        hs_count(1) = dset_dims(1)
        hs_count(2:rank) = 1
        hs_offset(1) = 0
        hs_offset(2:rank) = offset(1:rank - 1) - 1
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_offset(1:rank), hs_count(1:rank), hdferror)

        ! set mspace to a vector
        mdims(1) = size(vector, 1)
        call h5screate_simple_f(1, mdims, mspace_id, hdferror)

        ! write out vector
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, vector, mdims, hdferror, mspace_id, dspace_id)

        ! close mspace_id
        call h5sclose_f(mspace_id, hdferror)

      else
        write (format_string, '(A,I0,A,I0,A)') '(A,', rank - 1, '(I0,A),A,', rank - 1, '(I0,A),A)'
        write (*, format_string) "--->ERROR: offset=(", (offset(i), ',', i=1, rank - 1), &
          "), is not constent with dset_dims(2:rank)=(", (dset_dims(i), ',', i=2, rank), ")"
      end if

    else
      write (*, '(A,I0,A,I0)') "--->ERROR: size(vector)=", size(vector), &
        ", is not constent with dset_dims(1)=", dset_dims(1)
    end if

    ! close id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_vector_from_dataset_real

  !
  subroutine hdf_read_vector_from_dataset_double(loc_id, dset_name, offset, vector)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(in) :: offset(:)            ! position within dataset
    real(dp), intent(out) :: vector(:)          ! data to be written

    integer(HID_T) :: dset_id, dspace_id, mspace_id
    integer :: rank
    integer(HSIZE_T) :: dset_dims(6), max_dims(6), mdims(1)
    integer(HSIZE_T) :: hs_count(6), hs_offset(6)
    integer :: hdferror

    integer :: i
    character(len=32) :: format_string

    vector = 0.0
    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_vector_from_dataset_double: "//trim(dset_name)
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! get dataspace
    call h5dget_space_f(dset_id, dspace_id, hdferror)

    ! get rank (ndims), and dims
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)
    call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)

    ! check size and offset
    if (size(vector, 1) == dset_dims(1)) then

      if (all(offset(1:rank - 1) <= dset_dims(2:rank)) .and. all(offset(1:rank - 1) > 0)) then

        ! select hyperslab in dataset (note: convert FORTRAN offset to C offset by subtracting 1)
        hs_count(1) = dset_dims(1)
        hs_count(2:rank) = 1
        hs_offset(1) = 0
        hs_offset(2:rank) = offset(1:rank - 1) - 1
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_offset(1:rank), hs_count(1:rank), hdferror)

        ! set mspace to a vector
        mdims(1) = size(vector, 1)
        call h5screate_simple_f(1, mdims, mspace_id, hdferror)

        ! write out vector
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vector, mdims, hdferror, mspace_id, dspace_id)

        ! close mspace_id
        call h5sclose_f(mspace_id, hdferror)

      else
        write (format_string, '(A,I0,A,I0,A)') '(A,', rank - 1, '(I0,A),A,', rank - 1, '(I0,A),A)'
        write (*, format_string) "--->ERROR: offset=(", (offset(i), ',', i=1, rank - 1), &
          "), is not constent with dset_dims(2:rank)=(", (dset_dims(i), ',', i=2, rank), ")"
      end if

    else
      write (*, '(A,I0,A,I0)') "--->ERROR: size(vector)=", size(vector), &
        ", is not constent with dset_dims(1)=", dset_dims(1)
    end if

    ! close id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_vector_from_dataset_double

  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_set_property_list--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief set property list for chunking and compression
  subroutine hdf_set_property_list(plist_id, rank, dims, cdims, filter)

    integer(HID_T), intent(inout) :: plist_id
    integer, intent(in) :: rank
    integer(SIZE_T), intent(in) :: dims(:)
    integer(SIZE_T), intent(inout) :: cdims(:)
    character(len=*), intent(in) :: filter

    integer :: hdferror
    integer :: szip_pixels_per_block, szip_options_mask

    ! set chunk (if needed)
    if (filter == 'none') then
      if (all(cdims .ne. 0)) then
        call h5pset_chunk_f(plist_id, rank, cdims, hdferror)
      end if
    else
      if (any(cdims == 0)) then
        cdims = dims
      end if
      call h5pset_chunk_f(plist_id, rank, cdims, hdferror)
    end if

    ! set filter
    select case (filter)
    case ('none')
      continue
    case ('szip')
      if (hdf_szip_options == 'NN') then
        szip_options_mask = H5_SZIP_NN_OM_F
      elseif (hdf_szip_options == 'EC') then
        szip_options_mask = H5_SZIP_EC_OM_F
      else
        write (*, '(A)') "-->hdf_set_property_list: warning hdf_szip_options "//trim(hdf_szip_options)//" not recognized"
        continue
      end if
      szip_pixels_per_block = min(int(product(cdims), 4), hdf_szip_pixels_per_block)
      call H5Pset_szip_f(plist_id, szip_options_mask, szip_pixels_per_block, hdferror)
    case ('gzip')
      call h5pset_deflate_f(plist_id, hdf_gzip_level, hdferror)
    case ('gzip+shuffle')
      call h5pset_shuffle_f(plist_id, hdferror)
      call h5pset_deflate_f(plist_id, hdf_gzip_level, hdferror)
    case default
      write (*, '(A)') "-->hdf_set_property_list: warning filter "//trim(filter)//" not recognized"
    end select

  end subroutine hdf_set_property_list


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_dataset_integer--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief write a scalar to a hdf5 file
  subroutine hdf_write_dataset_integer_0(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array                      ! data to be written
    integer, optional, intent(in) :: chunks           ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: dset_id, dspace_id
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_0: "//trim(dset_name)
    end if

    ! set rank and dims
    dims = (/0/)

    !
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_0: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_0: warning chunks not used"
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_0

  !  \brief write a 1 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_1(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:)                   ! data to be written
    integer, optional, intent(in) :: chunks(1)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(1), cdims(1)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_1

  !  \brief write a 2 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_2(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:,:)                 ! data to be written
    integer, optional, intent(in) :: chunks(2)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(2), cdims(2)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_2

  !  \brief write a 3 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_3(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:,:,:)               ! data to be written
    integer, optional, intent(in) :: chunks(3)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(3), cdims(3)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_3

  !  \brief write a 4 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_4(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:,:,:,:)             ! data to be written
    integer, optional, intent(in) :: chunks(4)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(4), cdims(4)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_4

  !  \brief write a 5 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_5(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:,:,:,:,:)           ! data to be written
    integer, optional, intent(in) :: chunks(5)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(5), cdims(5)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_5

  !  \brief write a 6 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_6(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:,:,:,:,:,:)         ! data to be written
    integer, optional, intent(in) :: chunks(6)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(6), cdims(6)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 6
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_6


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_dataset_real--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief write a scalar to a hdf5 file
  subroutine hdf_write_dataset_real_0(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array                     ! data to be written
    integer, optional, intent(in) :: chunks           ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: dset_id, dspace_id
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_0: "//trim(dset_name)
    end if

    ! set rank and dims
    dims = (/0/)

    !
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_real_0: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_real_0: warning chunks not used"
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, dspace_id, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_0

  !  \brief write a 1 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_1(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:)                  ! data to be written
    integer, optional, intent(in) :: chunks(1)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(1), cdims(1)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_1

  !  \brief write a 2 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_2(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:,:)                ! data to be written
    integer, optional, intent(in) :: chunks(2)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(2), cdims(2)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_2

  !  \brief write a 3 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_3(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:,:,:)              ! data to be written
    integer, optional, intent(in) :: chunks(3)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(3), cdims(3)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_3

  !  \brief write a 4 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_4(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:,:,:,:)            ! data to be written
    integer, optional, intent(in) :: chunks(4)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(4), cdims(4)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_4

  !  \brief write a 5 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_5(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:,:,:,:,:)          ! data to be written
    integer, optional, intent(in) :: chunks(5)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(5), cdims(5)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_5

  !  \brief write a 6 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_6(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:,:,:,:,:,:)        ! data to be written
    integer, optional, intent(in) :: chunks(6)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(6), cdims(6)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 6
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_6


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_dataset_double--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief write a scalar to a hdf5 file
  subroutine hdf_write_dataset_double_0(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array                     ! data to be written
    integer, optional, intent(in) :: chunks           ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: dset_id, dspace_id
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_0: "//trim(dset_name)
    end if

    ! set rank and dims
    dims = (/0/)

    !
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_double_0: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_double_0: warning chunks not used"
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_0

  !  \brief write a 1 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_1(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:)                  ! data to be written
    integer, optional, intent(in) :: chunks(1)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(1), cdims(1)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_1

  !  \brief write a 2 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_2(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:,:)                ! data to be written
    integer, optional, intent(in) :: chunks(2)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(2), cdims(2)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_2

  !  \brief write a 3 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_3(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:,:,:)              ! data to be written
    integer, optional, intent(in) :: chunks(3)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(3), cdims(3)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_3

  !  \brief write a 4 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_4(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:,:,:,:)            ! data to be written
    integer, optional, intent(in) :: chunks(4)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(4), cdims(4)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_4

  !  \brief write a 5 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_5(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:,:,:,:,:)          ! data to be written
    integer, optional, intent(in) :: chunks(5)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(5), cdims(5)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_5

  !  \brief write a 6 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_6(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:,:,:,:,:,:)        ! data to be written
    integer, optional, intent(in) :: chunks(6)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(6), cdims(6)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 6
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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_6


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_dataset_character--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief write a scalar to a hdf5 file
  subroutine hdf_write_dataset_character_0(loc_id, dset_name, array, chunks, filter, processor)

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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, dtype_id, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, dtype_id, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)
    call h5tclose_f(dtype_id, hdferror)

  end subroutine hdf_write_dataset_character_0

  !  \brief write a 1 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_character_1(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    character(len=*), intent(in) :: array(:)          ! data to be written
    integer, optional, intent(in) :: chunks(1)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(1), cdims(1), length
    integer(HID_T) :: dset_id, dspace_id, plist_id, dtype_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_character_1: "//trim(dset_name)
    end if

    ! set rank, dims, dtype
    rank = 1
    dims = shape(array, KIND=HID_T)
    length=len(array(1))
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
    call h5dcreate_f(loc_id, dset_name, dtype_id, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, dtype_id, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, dtype_id, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)
    call h5tclose_f(dtype_id, hdferror)

  end subroutine hdf_write_dataset_character_1

  !  \brief write a 2 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_character_2(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    character(len=*), intent(in) :: array(:,:)        ! data to be written
    integer, optional, intent(in) :: chunks(2)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(2), cdims(2), length
    integer(HID_T) :: dset_id, dspace_id, plist_id, dtype_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_character_2: "//trim(dset_name)
    end if

    ! set rank, dims, dtype
    rank = 2
    dims = shape(array, KIND=HID_T)
    length=len(array(1,1))
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
    call h5dcreate_f(loc_id, dset_name, dtype_id, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, dtype_id, array, dims, hdferror, xfer_prp=dplist_collective)
    else
      if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, dtype_id, array, dims, hdferror, xfer_prp=dplist_independent)
      end if
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)
    call h5tclose_f(dtype_id, hdferror)

  end subroutine hdf_write_dataset_character_2


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_dataset_complex_double--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief write a scalar to a hdf5 file
  subroutine hdf_write_dataset_complex_double_0(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array                  ! data to be written
    integer, optional, intent(in) :: chunks           ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: dset_id, dspace_id
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_0: "//trim(dset_name)
    end if

    ! set rank and dims
    dims = (/0/)

    !
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_0: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_0: warning chunks not used"
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
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, dspace_id, dset_id, hdferror)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_0

  !  \brief write a 1 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_1(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:)               ! data to be written
    integer, optional, intent(in) :: chunks(1)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(1), cdims(1)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
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
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_1

  !  \brief write a 2 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_2(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:,:)             ! data to be written
    integer, optional, intent(in) :: chunks(2)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(2), cdims(2)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
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
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_2

  !  \brief write a 3 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_3(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:,:,:)           ! data to be written
    integer, optional, intent(in) :: chunks(3)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(3), cdims(3)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
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
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_3

  !  \brief write a 4 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_4(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:,:,:,:)         ! data to be written
    integer, optional, intent(in) :: chunks(4)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(4), cdims(4)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
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
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_4

  !  \brief write a 5 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_5(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:,:,:,:,:)       ! data to be written
    integer, optional, intent(in) :: chunks(5)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(5), cdims(5)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
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
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_5

  !  \brief write a 6 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_6(loc_id, dset_name, array, chunks, filter, processor)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:,:,:,:,:,:)     ! data to be written
    integer, optional, intent(in) :: chunks(6)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that write the data. (-1: collective write, i.e. same data on each processor_write)

    integer :: rank
    integer(SIZE_T) :: dims(6), cdims(6)
    integer(HID_T) :: dset_id, dspace_id, plist_id
    character(len=32) :: filter_case
    integer :: hdferror
    integer :: processor_write

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 6
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
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
    !write(*,'(A20,I0)') "h5dcreate: ", hdferror

    ! write dataset
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! close all id's
    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_6


  !!---------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_integer--------------------------------
  !!---------------------------------------------------------------------------------------

  !  \brief reads a scalar from an hdf5 file
  subroutine hdf_read_dataset_integer_0(loc_id, dset_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, intent(out) :: array                ! data to be written

    integer(SIZE_T) :: dims(1)
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
    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(2)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(3)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(4)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(5)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(6)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 6
    dims = shape(array, KIND=HID_T)

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

    integer(SIZE_T) :: dims(1)
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
    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(2)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(3)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(4)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(5)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(6)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 6
    dims = shape(array, KIND=HID_T)

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

    integer(SIZE_T) :: dims(1)
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
    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(2)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(3)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(4)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(5)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(6)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 6
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(1), length
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
    integer(SIZE_T) :: dims_in(1), dims(1), max_dims(1), ii, length
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
    integer(SIZE_T) :: dims_in(2), dims(2), max_dims(2), ii, jj, length
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

    integer(SIZE_T) :: dims(1)
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
    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_1: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 1
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(2)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_2: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(3)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_3: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 3
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(4)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_4: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 4
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(5)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_5: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 5
    dims = shape(array, KIND=HID_T)

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
    integer(SIZE_T) :: dims(6)
    integer(HID_T) :: dset_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_6: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = 2
    dims = shape(array, KIND=HID_T)

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

    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror

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
    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_integer_1: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create dataspace
    rank = 1
    dims = shape(array, KIND=HID_T)
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

  !  \brief writes a scalar attribute
  subroutine hdf_write_attr_real_0(loc_id, obj_name, attr_name, array)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: obj_name    ! object name attribute will be attached to (if "" use loc_id)
    character(len=*), intent(in) :: attr_name   ! name of attribute
    real(sp), intent(in) :: array                ! data to write to attribute

    !integer :: rank
    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_real_0: "//trim(obj_name)//"/"//trim(attr_name)
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
    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_real_1: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create dataspace
    rank = 1
    dims = shape(array, KIND=HID_T)
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
    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_double_0: "//trim(obj_name)//"/"//trim(attr_name)
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
    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, aspace_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_double_1: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create dataspace
    rank = 1
    dims = shape(array, KIND=HID_T)
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

    integer(SIZE_T) :: dims(1)
    integer(HID_T) :: obj_id, type_id, aspace_id, attr_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_attr_string: "//trim(obj_name)//"/"//trim(attr_name)
    end if

    ! open object
    if (obj_name == "") then
      obj_id = loc_id
    else
      call h5oopen_f(loc_id, obj_name, obj_id, hdferror)
    end if

    ! create type_id and aspace_id
    dims(1) = len(array, KIND=HID_T)
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

    integer(SIZE_T) :: dims(1)
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
    integer(SIZE_T) :: dims(1)
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
    dims = shape(array, KIND=HID_T)

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

    integer(SIZE_T) :: dims(1)
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
    integer(SIZE_T) :: dims(1)
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
    dims = shape(array, KIND=HID_T)

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

    integer(SIZE_T) :: dims(1)
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
    integer(SIZE_T) :: dims(1)
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
    dims = shape(array, KIND=HID_T)

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

    integer(SIZE_T) :: dims(1)
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
    dims(1) = len(array, KIND=HID_T)
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
