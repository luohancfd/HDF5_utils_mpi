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
  public :: hdf_exists, hdf_get_rank, hdf_get_dims, hdf_set_dims
  public :: hdf_write_dataset, hdf_read_dataset
  public :: hdf_write_attribute, hdf_read_attribute
  public :: hdf_create_dataset
  public :: hdf_write_vector_to_dataset, hdf_read_vector_from_dataset
  public :: HID_T, hdf_set_print_messages, hdf_set_default_filter
  public :: hdf_set_even_offset

  private :: hdf_preset_prop, hdf_close_prop

  !>  \brief Generic interface to write a dataset
  !>  Supported types
  !>   - integers (scalar and 1d-6d arrays)
  !>   - doubles (scalar and 1d-6d arrays)
  !>   - reals (scalar and 1d-6d arrays)
  !>   - string (scalar and 1d-2d arrays)
  !>   - complex double number (compound data type, "r"/"i" for real and imaginary part, scalar and 1d-6d arrays)
  !>
  !>  \param[in] loc_id     local id in file, e.g. file_id
  !>  \param[in] dset_name  name of dataset, NOTE: HDF5 assumes the dataset doesn't exist before !!
  !>  \param[in] array      data array to be written
  !>  \param[in] chunks     (optional) chunk size for dataset
  !>  \param[in] filter     (optional) filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
  !>  \param[in] processor  (optional, default=-1) processor that provides the data, -1 if the data is the same on all processors.
  !>
  !>  if processor != -1, the following options become useless
  !>  \param[in] axis       (optional, default=-1) dimension on which the data will be stacked, starting from 1
  !>
  !>  \details
  !>    - a variable with the name "num_reacts". All processors have the same data. We only want one copy saved in the final file
  !>          call hdf_write_dataset(file_id, "num_reacts", num_reacts)
  !>    - a variable with the name "num_reacts_0". All processors have this variable but different shape and content.
  !>          we only want to save the copy on processor 0 (first processor)
  !>          call hdf_write_dataset(file_id, "num_reacts_0", num_reacts_0, processor=0)
  !>    - a scalar variable with the name "num_particle". All processors have this variable with different values. We want to save all of them.
  !>          call hdf_write_dataset(file_id, "num_particle", num_particle, axis=1)
  !>       (Note: the result of this call is an array with its size equal to the number of processors will be save to the hdf5 file. If you want to add
  !>        them up and then save the result, you need to do the calculation by yourself)
  !>    - an array "x_loc(XSIZE, YSIZE, ZSIZE)". All processors have this variable with different values. In addition, the meaningful data is not the whole
  !>      array but "x_loc(xstart:xend, ystart:yend, 1:num)" (xstart, xend, ystart, yend are the same across processors but num is has different values
  !>      on each processor). To stack the variable on the third dimension, use the following call. The API will set hyperslab based on the shape of array provided
  !>          call hdf_write_dataset(file_id, "x_loc", x_loc(xstart:xend, ystart:yend, 1:num), axis=3)
  !>    The function will also save the number of rows contributed by each processor as an attribute.


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

  interface hdf_set_even_offset
    module procedure hdf_set_even_offset_from_dataset
    module procedure hdf_set_even_offset_from_number
  end interface hdf_set_even_offset

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
    module procedure hdf_write_attr_integer_1_8  ! kind = 8 integer
    module procedure hdf_write_attr_real_0
    module procedure hdf_write_attr_real_1
    module procedure hdf_write_attr_double_0
    module procedure hdf_write_attr_double_1
    module procedure hdf_write_attr_string
  end interface hdf_write_attribute

  !>  \brief Get the appropriate mpi integer type
  !>  \param[in] int_number  an integer number
  !>  \return    MPI_TYPE for MPI communication
  interface hdf_get_mpi_int
    module procedure hdf_get_mpi_int_4
    module procedure hdf_get_mpi_int_8
    ! module procedure hdf_get_mpi_int_16
  end interface hdf_get_mpi_int

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
    module procedure hdf_read_attr_integer_1_8
    module procedure hdf_read_attr_real_0
    module procedure hdf_read_attr_real_1
    module procedure hdf_read_attr_double_0
    module procedure hdf_read_attr_double_1
    module procedure hdf_read_attr_string
  end interface hdf_read_attribute

  interface hdf_get_dims
    module procedure hdf_get_dims_4
    module procedure hdf_get_dims_8
  end interface hdf_get_dims

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
  integer :: mpi_comm, mpi_irank, mpi_nrank, mpi_ierr, mpi_hsize_t
  integer(HID_T) :: file_plist_id   !< parallel file access property
  integer(HID_T) :: dplist_collective, dplist_independent !< dataset access property
  integer :: mpi_nrank_old

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
    integer(HSIZE_T) :: temp

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
    mpi_hsize_t = hdf_get_mpi_int(temp)

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

    mpi_nrank_old = -1
    if ((status2 .ne. 'OLD') .or. (action2 == 'WRITE') .or. (action2 == 'READWRITE')) then
      call hdf_preset_file_attribute(file_id)
    else
      call hdf_read_attribute(file_id, "", 'mpi_nrank', mpi_nrank_old)
    end if

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

  !>  \brief Preset some file-level attributes
  subroutine hdf_preset_file_attribute(file_id)
    implicit none
    integer(HID_T), intent(in) :: file_id            !< HDF5 id of the file

    call hdf_write_attribute(file_id, "", 'mpi_nrank', mpi_nrank)
  end subroutine hdf_preset_file_attribute

  !>  \brief Preset some properties
  subroutine hdf_preset_prop()
    implicit none
    integer :: ii
    integer(HSIZE_T) ::offset
    INTEGER(HSIZE_T)     ::   type_sizei  ! Size of the integer datatype
    INTEGER(HSIZE_T)     ::   type_sized  ! Size of the double precision datatype
    INTEGER(HSIZE_T)     ::   type_sizer  ! Size of the real datatype
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

  !>  \brief Set even offset for read_dataset
  subroutine hdf_set_even_offset_from_dataset(file_id, loc_id, dset_name, offset, new_size)

    integer(HID_T), intent(in) :: file_id, loc_id  !< location id
    character(len=*), intent(in) :: dset_name      !< dataset name
    !< size for the read buffer each processor can have different value
    integer :: offset(:)
    integer, intent(in), optional :: new_size

    integer :: mpi_nrank_, ii, jj, kk
    integer, allocatable :: count_old(:), total_count

    if (size(offset) .ne. mpi_nrank) then
      write(*,'(A)') "-->hdf_set_even_offset: size of offset is wrong"
      call MPI_Abort(mpi_comm, 1, mpi_ierr)
    end if

    call hdf_read_attribute(file_id, "", "mpi_nrank", mpi_nrank_)

    ! sanity check
    call hdf_read_attribute(loc_id, dset_name, "axis_write", ii)
    if (ii == -1) then
      write(*,'(A)') "-->hdf_set_even_offset: dataset ["//dset_name//"] is not written in parallel"
      call MPI_Abort(mpi_comm, 1, mpi_ierr)
    end if

    ! get old data size
    allocate(count_old(mpi_nrank_))
    call hdf_read_attribute(loc_id, dset_name, "count", count_old)
    total_count = sum(count_old)
    deallocate(count_old)

    ! calculate offset of the stacked axis
    jj = floor(dble(total_count)/dble(mpi_nrank))
    if (jj == 0) then
        write(*,'(A)') "-->hdf_set_even_offset: some processor may get nothing"
        call MPI_Abort(mpi_comm, 1, mpi_ierr)
    end if

    kk = total_count - jj * mpi_nrank
    offset(1) = 0
    do ii = 1,  kk
      offset(ii+1) = offset(ii) + jj + 1
    end do
    do ii = kk + 1, mpi_nrank-1
      offset(ii+1) = offset(ii) + jj
    end do

    ! sanity check whether we have enough space
    if (present(new_size)) then
      if (mpi_irank == mpi_nrank - 1) then
        ii = total_count - offset(mpi_irank+1)
      else
        ii = offset(mpi_irank + 2) - offset(mpi_irank + 1)
      end if
      if (ii > new_size) then
        write(*,'(A)') "-->hdf_set_even_offset: new_size is not enough for dataset ["//dset_name//"]"
        call MPI_Abort(mpi_comm, 1, mpi_ierr)
      end if
    end if
  end subroutine hdf_set_even_offset_from_dataset

  subroutine hdf_set_even_offset_from_number(total_count, nprocs, offset)
    integer, intent(in) :: total_count
    integer, intent(in) :: nprocs
    integer :: offset(nprocs)

    integer :: ii, jj, kk

    kk = total_count - jj * nprocs
    offset(1) = 0
    do ii = 1,  kk
      offset(ii+1) = offset(ii) + jj + 1
    end do
    do ii = kk + 1, nprocs-1
      offset(ii+1) = offset(ii) + jj
    end do
  end subroutine hdf_set_even_offset_from_number

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
  subroutine hdf_get_dims_4(loc_id, dset_name, dims)

    integer(HID_T), intent(in) :: loc_id        !< location id
    character(len=*), intent(in) :: dset_name   !< name of dataset
    integer, intent(out) :: dims(:)             !< dimensions of the dataset

    integer(HID_T) :: dset_id, dspace_id
    integer :: rank
    integer(HSIZE_T) :: dset_dims(6), max_dims(6)
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "->hdf_get_dims_4"
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

  end subroutine hdf_get_dims_4

  !>  \brief get the dimensions of a dataset
  subroutine hdf_get_dims_8(loc_id, dset_name, dims)

    integer(HID_T), intent(in) :: loc_id        !< location id
    character(len=*), intent(in) :: dset_name   !< name of dataset
    integer(HSIZE_T), intent(out) :: dims(:)             !< dimensions of the dataset

    integer(HID_T) :: dset_id, dspace_id
    integer :: rank
    integer(HSIZE_T) :: dset_dims(6), max_dims(6)
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "->hdf_get_dims_8"
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! get dataspace
    call h5dget_space_f(dset_id, dspace_id, hdferror)

    ! get rank (ndims)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)

    ! get dims
    call h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), hdferror)
    dims(1:rank) = dset_dims(1:rank)

    ! close id's
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_get_dims_8

  !>  \brief set the dimensions of a dataset
  subroutine hdf_set_dims(loc_id, dset_name, dims)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, allocatable :: dims(:)

    integer :: rank, hdferror
    integer(HID_T) :: dset_id, dspace_id
    integer(HSIZE_T) :: dims_data(6), max_dims(6)

    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)
    call h5dget_space_f(dset_id, dspace_id, hdferror)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferror)

    if (allocated(dims)) then
      deallocate(dims)
    end if
    allocate(dims(rank))

    call h5sget_simple_extent_dims_f(dspace_id, dims_data(1:rank), max_dims(1:rank), hdferror)
    dims = int(dims_data(1:rank))
    call h5sclose_f(dspace_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_set_dims

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
    integer(HSIZE_T) :: dims(8)
    integer(HID_T) :: dset_id, dspace_id
    integer :: hdferror

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_create_dataset: "//trim(dset_name)
    end if

    ! set rank and dims
    rank = size(dset_dims, 1)
    dims(1:rank) = int(dset_dims, HSIZE_T)

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
    integer(HSIZE_T), intent(in) :: dims(:)
    integer(HSIZE_T), intent(inout) :: cdims(:)
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
  !!--------------------------------hdf_get_mpi_int-----------------------------------------
  !!----------------------------------------------------------------------------------------

  function hdf_get_mpi_int_4(int_number) result(r)
    implicit none
    integer(kind=4), INTENT(IN) :: int_number
    integer :: r
    r = MPI_INTEGER4
  end function hdf_get_mpi_int_4

  function hdf_get_mpi_int_8(int_number) result(r)
    implicit none
    integer(kind=8), INTENT(IN) :: int_number
    integer :: r
    r = MPI_INTEGER8
  end function hdf_get_mpi_int_8

  ! function hdf_get_mpi_int_16(int_number) result(r)
  !   implicit none
  !   integer(kind=16), INTENT(IN) :: int_number
  !   integer :: r
  !   r = MPI_INTEGER16
  ! end function hdf_get_mpi_int_16