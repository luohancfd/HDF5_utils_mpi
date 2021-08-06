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
  subroutine hdf_set_even_offset(file_id, loc_id, dset_name, offset, new_size)

    integer(HID_T), intent(in) :: file_id, loc_id  !< location id
    character(len=*), intent(in) :: dset_name      !< dataset name
    !< size for the read buffer each processor can have different value
    integer, allocatable :: offset(:)
    integer, intent(in), optional :: new_size

    integer :: mpi_nrank_, ii, jj, kk
    integer, allocatable :: count_old(:), total_count

    ! ignore whatever is set originally
    if (allocated(offset)) then
      deallocate(offset)
    end if
    allocate(offset(mpi_nrank))

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
  end subroutine hdf_set_even_offset

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
  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_dataset_integer--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief write a scalar to a hdf5 file
  subroutine hdf_write_dataset_integer_0(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array                      ! data to be written
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

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_0: "//trim(dset_name)
    end if

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

    ! set rank and dims
    rank = 1
    dimsm = (/0/)
    dimsf = (/mpi_nrank/)

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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, file_space_id, dset_id, hdferror)

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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',     axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_0

  !  \brief write a 1 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_1(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:)                   ! data to be written
    integer, optional, intent(in) :: chunks(1)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(1) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_1: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_1: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_1: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_integer_1: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 1
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 1 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 1
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_integer_1: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_1

  !  \brief write a 2 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_2(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:,:)                 ! data to be written
    integer, optional, intent(in) :: chunks(2)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(2) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_2: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_2: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_2: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_integer_2: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 2
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 2 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 2
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_integer_2: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_2

  !  \brief write a 3 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_3(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:,:,:)               ! data to be written
    integer, optional, intent(in) :: chunks(3)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(3) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_3: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_3: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_3: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_integer_3: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 3
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 3 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 3
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_integer_3: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_3

  !  \brief write a 4 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_4(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:,:,:,:)             ! data to be written
    integer, optional, intent(in) :: chunks(4)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(4) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_4: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_4: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_4: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_integer_4: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 4
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 4 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 4
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_integer_4: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_4

  !  \brief write a 5 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_5(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:,:,:,:,:)           ! data to be written
    integer, optional, intent(in) :: chunks(5)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(5) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_5: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_5: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_5: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_integer_5: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 5
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 5 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 5
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_integer_5: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_5

  !  \brief write a 6 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_integer_6(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    integer, intent(in) :: array(:,:,:,:,:,:)         ! data to be written
    integer, optional, intent(in) :: chunks(6)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(6) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_integer_6: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_6: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_integer_6: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_integer_6: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 6
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 6 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 6
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_integer_6: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_integer_6


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_dataset_real--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief write a scalar to a hdf5 file
  subroutine hdf_write_dataset_real_0(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array                     ! data to be written
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

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_0: "//trim(dset_name)
    end if

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

    ! set rank and dims
    rank = 1
    dimsm = (/0/)
    dimsf = (/mpi_nrank/)

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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, file_space_id, dset_id, hdferror)

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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',     axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_0

  !  \brief write a 1 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_1(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:)                  ! data to be written
    integer, optional, intent(in) :: chunks(1)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(1) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_1: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_real_1: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_real_1: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_real_1: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 1
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 1 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 1
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_real_1: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_1

  !  \brief write a 2 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_2(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:,:)                ! data to be written
    integer, optional, intent(in) :: chunks(2)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(2) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_2: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_real_2: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_real_2: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_real_2: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 2
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 2 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 2
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_real_2: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_2

  !  \brief write a 3 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_3(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:,:,:)              ! data to be written
    integer, optional, intent(in) :: chunks(3)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(3) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_3: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_real_3: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_real_3: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_real_3: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 3
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 3 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 3
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_real_3: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_3

  !  \brief write a 4 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_4(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:,:,:,:)            ! data to be written
    integer, optional, intent(in) :: chunks(4)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(4) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_4: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_real_4: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_real_4: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_real_4: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 4
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 4 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 4
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_real_4: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_4

  !  \brief write a 5 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_5(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:,:,:,:,:)          ! data to be written
    integer, optional, intent(in) :: chunks(5)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(5) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_5: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_real_5: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_real_5: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_real_5: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 5
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 5 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 5
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_real_5: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_5

  !  \brief write a 6 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_real_6(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(sp), intent(in) :: array(:,:,:,:,:,:)        ! data to be written
    integer, optional, intent(in) :: chunks(6)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(6) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_real_6: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_real_6: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_real_6: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_real_6: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 6
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 6 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 6
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_real_6: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_real_6


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_dataset_double--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief write a scalar to a hdf5 file
  subroutine hdf_write_dataset_double_0(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array                     ! data to be written
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

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_0: "//trim(dset_name)
    end if

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

    ! set rank and dims
    rank = 1
    dimsm = (/0/)
    dimsf = (/mpi_nrank/)

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
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, file_space_id, dset_id, hdferror)

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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',     axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_0

  !  \brief write a 1 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_1(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:)                  ! data to be written
    integer, optional, intent(in) :: chunks(1)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(1) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_1: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_double_1: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_double_1: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_double_1: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 1
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 1 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 1
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_double_1: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_1

  !  \brief write a 2 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_2(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:,:)                ! data to be written
    integer, optional, intent(in) :: chunks(2)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(2) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_2: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_double_2: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_double_2: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_double_2: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 2
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 2 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 2
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_double_2: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_2

  !  \brief write a 3 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_3(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:,:,:)              ! data to be written
    integer, optional, intent(in) :: chunks(3)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(3) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_3: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_double_3: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_double_3: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_double_3: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 3
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 3 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 3
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_double_3: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_3

  !  \brief write a 4 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_4(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:,:,:,:)            ! data to be written
    integer, optional, intent(in) :: chunks(4)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(4) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_4: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_double_4: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_double_4: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_double_4: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 4
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 4 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 4
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_double_4: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_4

  !  \brief write a 5 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_5(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:,:,:,:,:)          ! data to be written
    integer, optional, intent(in) :: chunks(5)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(5) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_5: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_double_5: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_double_5: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_double_5: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 5
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 5 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 5
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_double_5: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_5

  !  \brief write a 6 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_double_6(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    real(dp), intent(in) :: array(:,:,:,:,:,:)        ! data to be written
    integer, optional, intent(in) :: chunks(6)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(6) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_double_6: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_double_6: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_double_6: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_double_6: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 6
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 6 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 6
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_double_6: "//trim(dset_name), &
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
    if (processor_write == -1) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_collective)
    else if (mpi_irank == processor_write) then
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsf, hdferror, &
                      file_space_id=file_space_id,               &
                      mem_space_id=mem_space_id,                 &
                      xfer_prp=dplist_independent)
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_double_6


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_dataset_character--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief write a scalar to a hdf5 file
  subroutine hdf_write_dataset_character_0(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    character(len=*), intent(in) :: array             ! data to be written
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
    integer(HSIZE_T) :: length, length_glob
    integer(HID_T) :: dtype_id
    character(len=:),allocatable :: buffer

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_character_0: "//trim(dset_name)
    end if

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

    ! set rank and dims
    rank = 1
    dimsm = (/0/)
    dimsf = (/mpi_nrank/)
    length=len(array)
    length_glob = length
    if (processor_write .ne. -1) then
      call MPI_Bcast(length, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    else
      call MPI_Allreduce(length, length_glob, 1, mpi_hsize_t, MPI_MAX, mpi_comm, mpi_ierr)
      if (length .ne. length_glob) then
        allocate(character(len=length_glob)::buffer)
        buffer = array
      end if
    end if
    call h5tcopy_f(H5T_FORTRAN_S1, dtype_id, hdferror)
    call h5tset_size_f(dtype_id, length_glob, hdferror)

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
    call h5dcreate_f(loc_id, dset_name, dtype_id, file_space_id, dset_id, hdferror)

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
    if (length .eq. length_glob) then
      if (processor_write == -1) then
        call h5dwrite_f(dset_id, dtype_id, array, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_collective)
      else if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, dtype_id, array, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_independent)
      end if
    else
      if (processor_write == -1) then
        call h5dwrite_f(dset_id, dtype_id, buffer, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_collective)
      else if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, dtype_id, buffer, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_independent)
      end if
    end if

    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',     axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if
    call hdf_write_attribute(dset_id, '', 'char_length', int(length_glob, kind=4))

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)
    call h5tclose_f(dtype_id, hdferror)
    if (length .ne. length_glob .and. processor_write .eq. -1) then
      deallocate(buffer)
    end if

  end subroutine hdf_write_dataset_character_0

  !  \brief write a 1 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_character_1(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    character(len=*), intent(in) :: array(:)          ! data to be written
    integer, optional, intent(in) :: chunks(1)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(1) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer(HSIZE_T) :: length, length_glob
    integer(HID_T) :: dtype_id
    integer :: i_1
    character(len=:), dimension(:), allocatable :: buffer

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_character_1: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_character_1: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_character_1: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_character_1: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 1
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm
    length=len(array(1))
    length_glob = length
    if (processor_write .ne. -1) then
      call MPI_Bcast(length, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    else
      call MPI_Allreduce(length, length_glob, 1, mpi_hsize_t, MPI_MAX, mpi_comm, mpi_ierr)
      if (length .ne. length_glob) then
        allocate(character(len=length_glob)::buffer(dimsm(1)))
        do i_1 = 1, int(dimsm(1))
          buffer(i_1) = array(i_1)
        end do
      end if
    end if
    call h5tcopy_f(H5T_FORTRAN_S1, dtype_id, hdferror)
    call h5tset_size_f(dtype_id, length_glob, hdferror)

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
    else if (axis_write > 1 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 1
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, dtype_id, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_character_1: "//trim(dset_name), &
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
    if (length .eq. length_glob) then
      if (processor_write == -1) then
        call h5dwrite_f(dset_id, dtype_id, array, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_collective)
      else if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, dtype_id, array, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_independent)
      end if
    else
      if (processor_write == -1) then
        call h5dwrite_f(dset_id, dtype_id, buffer, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_collective)
      else if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, dtype_id, buffer, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_independent)
      end if
    end if

    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if
    call hdf_write_attribute(dset_id, '', 'char_length', int(length_glob, kind=4))

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)
    call h5tclose_f(dtype_id, hdferror)
    if (length .ne. length_glob .and. processor_write .eq. -1) then
      deallocate(buffer)
    end if

  end subroutine hdf_write_dataset_character_1

  !  \brief write a 2 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_character_2(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    character(len=*), intent(in) :: array(:,:)        ! data to be written
    integer, optional, intent(in) :: chunks(2)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(2) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer(HSIZE_T) :: length, length_glob
    integer(HID_T) :: dtype_id
    integer :: i_1, i_2
    character(len=:), dimension(:,:), allocatable :: buffer

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_character_2: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_character_2: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_character_2: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_character_2: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 2
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm
    length=len(array(1,1))
    length_glob = length
    if (processor_write .ne. -1) then
      call MPI_Bcast(length, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    else
      call MPI_Allreduce(length, length_glob, 1, mpi_hsize_t, MPI_MAX, mpi_comm, mpi_ierr)
      if (length .ne. length_glob) then
        allocate(character(len=length_glob)::buffer(dimsm(1),dimsm(2)))
        do i_1 = 1, int(dimsm(1))
          do i_2 = 1, int(dimsm(2))
            buffer(i_1,i_2) = array(i_1,i_2)
          end do
        end do
      end if
    end if
    call h5tcopy_f(H5T_FORTRAN_S1, dtype_id, hdferror)
    call h5tset_size_f(dtype_id, length_glob, hdferror)

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
    else if (axis_write > 2 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 2
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, dtype_id, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_character_2: "//trim(dset_name), &
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
    if (length .eq. length_glob) then
      if (processor_write == -1) then
        call h5dwrite_f(dset_id, dtype_id, array, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_collective)
      else if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, dtype_id, array, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_independent)
      end if
    else
      if (processor_write == -1) then
        call h5dwrite_f(dset_id, dtype_id, buffer, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_collective)
      else if (mpi_irank == processor_write) then
        call h5dwrite_f(dset_id, dtype_id, buffer, dimsf, hdferror, &
                        file_space_id=file_space_id,               &
                        mem_space_id=mem_space_id,                 &
                        xfer_prp=dplist_independent)
      end if
    end if

    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if
    call hdf_write_attribute(dset_id, '', 'char_length', int(length_glob, kind=4))

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)
    call h5tclose_f(dtype_id, hdferror)
    if (length .ne. length_glob .and. processor_write .eq. -1) then
      deallocate(buffer)
    end if

  end subroutine hdf_write_dataset_character_2


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_write_dataset_complex_double--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief write a scalar to a hdf5 file
  subroutine hdf_write_dataset_complex_double_0(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array                  ! data to be written
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

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_0: "//trim(dset_name)
    end if

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

    ! set rank and dims
    rank = 1
    dimsm = (/0/)
    dimsf = (/mpi_nrank/)

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
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, file_space_id, dset_id, hdferror)

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
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',     axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_0

  !  \brief write a 1 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_1(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:)               ! data to be written
    integer, optional, intent(in) :: chunks(1)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(1) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_1: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_1: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_1: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_complex_double_1: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 1
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 1 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 1
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_complex_double_1: "//trim(dset_name), &
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
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_1

  !  \brief write a 2 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_2(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:,:)             ! data to be written
    integer, optional, intent(in) :: chunks(2)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(2) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_2: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_2: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_2: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_complex_double_2: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 2
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 2 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 2
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_complex_double_2: "//trim(dset_name), &
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
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_2

  !  \brief write a 3 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_3(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:,:,:)           ! data to be written
    integer, optional, intent(in) :: chunks(3)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(3) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_3: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_3: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_3: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_complex_double_3: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 3
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 3 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 3
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_complex_double_3: "//trim(dset_name), &
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
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_3

  !  \brief write a 4 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_4(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:,:,:,:)         ! data to be written
    integer, optional, intent(in) :: chunks(4)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(4) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_4: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_4: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_4: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_complex_double_4: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 4
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 4 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 4
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_complex_double_4: "//trim(dset_name), &
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
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_4

  !  \brief write a 5 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_5(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:,:,:,:,:)       ! data to be written
    integer, optional, intent(in) :: chunks(5)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(5) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_5: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_5: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_5: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_complex_double_5: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 5
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 5 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 5
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_complex_double_5: "//trim(dset_name), &
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
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_5

  !  \brief write a 6 - dimension array to a hdf5 file
  subroutine hdf_write_dataset_complex_double_6(loc_id, dset_name, array, chunks, filter, processor, axis)

    integer(HID_T), intent(in) :: loc_id              ! local id in file
    character(len=*), intent(in) :: dset_name         ! name of dataset
    complex(dp), intent(in) :: array(:,:,:,:,:,:)     ! data to be written
    integer, optional, intent(in) :: chunks(6)        ! chunk size for dataset
    character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
    integer, optional, intent(in) :: processor        ! processor that provides the data. (-1: collective write, i.e. same data on each processor_write)
    integer, optional, intent(in) :: axis             ! axis used to stack the data among processors

    integer :: rank, ii, jj
    integer(HSIZE_T), dimension(6) :: dimsf, dimsm, offset
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    ! character(len=32) :: filter_case
    integer :: hdferror, processor_write, axis_write, status(MPI_STATUS_SIZE)
    integer(HSIZE_T) :: offset_end
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_6: "//trim(dset_name)
    end if

    ! set filter
    if (present(filter)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_6: warning filter not used"
    end if

    ! set chunk (if needed)
    if (present(chunks)) then
      write (*, '(A)') "--->hdf_write_dataset_complex_double_6: warning chunks not used"
    end if

    ! set processor_write
    processor_write = -1
    if (present(processor)) then
      if ((processor .ge. 0) .and. (processor .lt. mpi_nrank)) then
        processor_write = processor
      else
        write (*, '(A, I2)') "--->hdf_write_dataset_complex_double_6: illegal processor = ", processor
      end if
    end if

    ! set rank and dims
    rank = 6
    dimsm = shape(array, KIND=HSIZE_T)
    if (processor_write .ne. -1) then
      call MPI_Bcast(dimsm, rank, mpi_hsize_t, processor_write, mpi_comm, mpi_ierr)
    end if
    dimsf = dimsm

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
    else if (axis_write > 6 .or. axis_write == 0) then
      write(*, '(A, I3, A, I3)') "Axis = ", axis_write, " is larger than array dim=", 6
      stop
    else
      call MPI_Allreduce(dimsm(axis_write), dimsf(axis_write), 1, mpi_hsize_t, MPI_SUM, mpi_comm, mpi_ierr)
      call h5screate_simple_f(rank, dimsf, file_space_id, hdferror)
    end if
    !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

    ! create dataset
    call h5dcreate_f(loc_id, dset_name, complexd_type_id, file_space_id, dset_id, hdferror)
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
                  write (*, '(A, A, I3, A)') "--->hdf_write_dataset_complex_double_6: "//trim(dset_name), &
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
    if (processor_write == -1) then
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
    end if
    !write(*,'(A20,I0)') "h5dwrite: ", hdferror

    ! write attributes
    call hdf_write_attribute(dset_id, '', 'processor', processor_write)
    call hdf_write_attribute(dset_id, '', 'axis_write',      axis_write)
    if (axis_write .ne. -1) then
      call hdf_write_attribute(dset_id, '', 'offset',    offset_glob)
      call hdf_write_attribute(dset_id, '', 'count',     count_glob)
      deallocate(offset_glob, count_glob)
    end if

    ! close all id's
    if (axis_write .ne. -1) then
      call h5sclose_f(mem_space_id, hdferror)
    end if
    call h5sclose_f(file_space_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_write_dataset_complex_double_6


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_integer--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief read a scalar from a hdf5 file
  subroutine hdf_read_dataset_integer_0(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    integer, intent(out) :: array               ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(1) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 0

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_0: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (is_parallel) then
      if (mpi_nrank .ne. mpi_nrank_old) then
        write(*, '(A)') "hdf_read_dataset_integer_0("//trim(dset_name)// &
                        "): different number of processors"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 1, mpi_nrank
        offset_glob(ii) = ii -1
        count_glob(ii)  = 1
      end do

      count_local = 1
    end if
    dimsm = (/1/)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_f(H5S_SCALAR_F, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                     array, &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_integer_0

  !  \brief read a 1 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_integer_1(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    integer, intent(out) :: array(:)            ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(1) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 1

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_1: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_integer_1("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_integer_1("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_integer_1("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_integer_1("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_integer_1("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_integer_1 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_integer_1 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                     array(1:count_local(1)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_integer_1

  !  \brief read a 2 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_integer_2(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    integer, intent(out) :: array(:,:)          ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(2) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 2

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_2: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_integer_2("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_integer_2("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_integer_2("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_integer_2("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_integer_2("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_integer_2 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_integer_2 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                     array(1:count_local(1),1:count_local(2)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_integer_2

  !  \brief read a 3 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_integer_3(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    integer, intent(out) :: array(:,:,:)        ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(3) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 3

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_3: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_integer_3("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_integer_3("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_integer_3("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_integer_3("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_integer_3("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_integer_3 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_integer_3 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_integer_3

  !  \brief read a 4 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_integer_4(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    integer, intent(out) :: array(:,:,:,:)      ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(4) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 4

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_4: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_integer_4("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_integer_4("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_integer_4("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_integer_4("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_integer_4("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_integer_4 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_integer_4 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_integer_4

  !  \brief read a 5 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_integer_5(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    integer, intent(out) :: array(:,:,:,:,:)    ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(5) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 5

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_5: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_integer_5("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_integer_5("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_integer_5("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_integer_5("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_integer_5("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_integer_5 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_integer_5 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4),1:count_local(5)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_integer_5

  !  \brief read a 6 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_integer_6(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    integer, intent(out) :: array(:,:,:,:,:,:)  ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(6) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 6

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_integer_6: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_integer_6("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_integer_6("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_integer_6("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_integer_6("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_integer_6("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_integer_6 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_integer_6 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4),1:count_local(5),1:count_local(6)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_integer_6


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_real--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief read a scalar from a hdf5 file
  subroutine hdf_read_dataset_real_0(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(sp), intent(out) :: array              ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(1) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 0

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_0: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (is_parallel) then
      if (mpi_nrank .ne. mpi_nrank_old) then
        write(*, '(A)') "hdf_read_dataset_real_0("//trim(dset_name)// &
                        "): different number of processors"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 1, mpi_nrank
        offset_glob(ii) = ii -1
        count_glob(ii)  = 1
      end do

      count_local = 1
    end if
    dimsm = (/1/)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_f(H5S_SCALAR_F, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_REAL, &
                     array, &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_real_0

  !  \brief read a 1 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_real_1(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(sp), intent(out) :: array(:)           ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(1) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 1

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_1: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_real_1("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_real_1("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_real_1("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_real_1("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_real_1("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_real_1 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_real_1 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_REAL, &
                     array(1:count_local(1)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_real_1

  !  \brief read a 2 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_real_2(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(sp), intent(out) :: array(:,:)         ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(2) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 2

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_2: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_real_2("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_real_2("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_real_2("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_real_2("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_real_2("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_real_2 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_real_2 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_REAL, &
                     array(1:count_local(1),1:count_local(2)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_real_2

  !  \brief read a 3 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_real_3(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(sp), intent(out) :: array(:,:,:)       ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(3) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 3

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_3: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_real_3("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_real_3("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_real_3("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_real_3("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_real_3("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_real_3 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_real_3 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_REAL, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_real_3

  !  \brief read a 4 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_real_4(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(sp), intent(out) :: array(:,:,:,:)     ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(4) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 4

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_4: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_real_4("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_real_4("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_real_4("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_real_4("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_real_4("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_real_4 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_real_4 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_REAL, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_real_4

  !  \brief read a 5 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_real_5(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(sp), intent(out) :: array(:,:,:,:,:)   ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(5) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 5

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_5: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_real_5("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_real_5("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_real_5("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_real_5("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_real_5("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_real_5 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_real_5 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_REAL, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4),1:count_local(5)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_real_5

  !  \brief read a 6 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_real_6(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(sp), intent(out) :: array(:,:,:,:,:,:) ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(6) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 6

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_real_6: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_real_6("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_real_6("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_real_6("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_real_6("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_real_6("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_real_6 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_real_6 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_REAL, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4),1:count_local(5),1:count_local(6)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_real_6


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_double--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief read a scalar from a hdf5 file
  subroutine hdf_read_dataset_double_0(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(dp), intent(out) :: array              ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(1) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 0

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_0: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (is_parallel) then
      if (mpi_nrank .ne. mpi_nrank_old) then
        write(*, '(A)') "hdf_read_dataset_double_0("//trim(dset_name)// &
                        "): different number of processors"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 1, mpi_nrank
        offset_glob(ii) = ii -1
        count_glob(ii)  = 1
      end do

      count_local = 1
    end if
    dimsm = (/1/)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_f(H5S_SCALAR_F, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, &
                     array, &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_double_0

  !  \brief read a 1 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_double_1(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(dp), intent(out) :: array(:)           ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(1) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 1

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_1: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_double_1("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_double_1("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_double_1("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_double_1("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_double_1("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_double_1 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_double_1 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, &
                     array(1:count_local(1)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_double_1

  !  \brief read a 2 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_double_2(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(dp), intent(out) :: array(:,:)         ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(2) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 2

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_2: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_double_2("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_double_2("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_double_2("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_double_2("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_double_2("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_double_2 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_double_2 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, &
                     array(1:count_local(1),1:count_local(2)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_double_2

  !  \brief read a 3 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_double_3(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(dp), intent(out) :: array(:,:,:)       ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(3) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 3

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_3: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_double_3("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_double_3("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_double_3("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_double_3("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_double_3("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_double_3 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_double_3 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_double_3

  !  \brief read a 4 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_double_4(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(dp), intent(out) :: array(:,:,:,:)     ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(4) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 4

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_4: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_double_4("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_double_4("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_double_4("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_double_4("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_double_4("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_double_4 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_double_4 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_double_4

  !  \brief read a 5 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_double_5(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(dp), intent(out) :: array(:,:,:,:,:)   ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(5) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 5

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_5: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_double_5("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_double_5("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_double_5("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_double_5("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_double_5("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_double_5 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_double_5 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4),1:count_local(5)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_double_5

  !  \brief read a 6 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_double_6(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    real(dp), intent(out) :: array(:,:,:,:,:,:) ! data to be read

    integer :: rank
    integer(HSIZE_T),dimension(6) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 6

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_double_6: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_double_6("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_double_6("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_double_6("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_double_6("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_double_6("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_double_6 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_double_6 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, &
                     array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4),1:count_local(5),1:count_local(6)), &
                     count_local, hdferror, &
                     mem_space_id=mem_space_id,                 &
                     file_space_id=file_space_id,               &
                     xfer_prp=dplist_collective)
    end if

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_double_6


  !!----------------------------------------------------------------------------------------
  !!--------------------------------hdf_read_dataset_complex_double--------------------------------
  !!----------------------------------------------------------------------------------------

  !  \brief read a scalar from a hdf5 file
  subroutine hdf_read_dataset_complex_double_0(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    complex(dp), intent(out) :: array           ! data to be read
    real(dp) :: buffer(2)                       ! buffer to save real and imag part
    integer :: rank
    integer(HSIZE_T),dimension(1) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 0

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_0: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (is_parallel) then
      if (mpi_nrank .ne. mpi_nrank_old) then
        write(*, '(A)') "hdf_read_dataset_complex_double_0("//trim(dset_name)// &
                        "): different number of processors"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 1, mpi_nrank
        offset_glob(ii) = ii -1
        count_glob(ii)  = 1
      end do

      count_local = 1
    end if
    dimsm = (/1/)

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    if (.not. is_parallel) then
      call h5dread_f(dset_id, complexd_field_id(1), buffer(1), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(2), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_f(H5S_SCALAR_F, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, complexd_field_id(1), buffer(1), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(2), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
    end if
    array = &
      cmplx(buffer(1), buffer(2), kind=dp)

    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_complex_double_0

  !  \brief read a 1 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_complex_double_1(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    complex(dp), intent(out) :: array(:)        ! data to be read
    real(dp), allocatable, dimension(:,:) :: buffer ! buffer to save real and imag part
    integer :: rank
    integer(HSIZE_T),dimension(1) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 1

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_1: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_complex_double_1("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_complex_double_1("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_complex_double_1("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_complex_double_1("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_complex_double_1("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_complex_double_1 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_complex_double_1 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    allocate (buffer(count_local(1), 2))
    if (.not. is_parallel) then
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,1), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,2), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,1), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,2), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
    end if
    array(1:count_local(1)) = &
      cmplx(buffer(:,1), buffer(:,2), kind=dp)

    deallocate (buffer)
    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_complex_double_1

  !  \brief read a 2 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_complex_double_2(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    complex(dp), intent(out) :: array(:,:)      ! data to be read
    real(dp), allocatable, dimension(:,:,:) :: buffer ! buffer to save real and imag part
    integer :: rank
    integer(HSIZE_T),dimension(2) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 2

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_2: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_complex_double_2("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_complex_double_2("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_complex_double_2("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_complex_double_2("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_complex_double_2("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_complex_double_2 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_complex_double_2 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    allocate (buffer(count_local(1), count_local(2), 2))
    if (.not. is_parallel) then
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,:,1), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,:,2), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,:,1), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,:,2), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
    end if
    array(1:count_local(1),1:count_local(2)) = &
      cmplx(buffer(:,:,1), buffer(:,:,2), kind=dp)

    deallocate (buffer)
    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_complex_double_2

  !  \brief read a 3 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_complex_double_3(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    complex(dp), intent(out) :: array(:,:,:)    ! data to be read
    real(dp), allocatable, dimension(:,:,:,:) :: buffer ! buffer to save real and imag part
    integer :: rank
    integer(HSIZE_T),dimension(3) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 3

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_3: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_complex_double_3("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_complex_double_3("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_complex_double_3("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_complex_double_3("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_complex_double_3("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_complex_double_3 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_complex_double_3 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    allocate (buffer(count_local(1), count_local(2), count_local(3), 2))
    if (.not. is_parallel) then
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,:,:,1), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,:,:,2), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,:,:,1), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,:,:,2), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
    end if
    array(1:count_local(1),1:count_local(2),1:count_local(3)) = &
      cmplx(buffer(:,:,:,1), buffer(:,:,:,2), kind=dp)

    deallocate (buffer)
    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_complex_double_3

  !  \brief read a 4 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_complex_double_4(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    complex(dp), intent(out) :: array(:,:,:,:)  ! data to be read
    real(dp), allocatable, dimension(:,:,:,:,:) :: buffer ! buffer to save real and imag part
    integer :: rank
    integer(HSIZE_T),dimension(4) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 4

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_4: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_complex_double_4("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_complex_double_4("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_complex_double_4("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_complex_double_4("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_complex_double_4("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_complex_double_4 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_complex_double_4 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    allocate (buffer(count_local(1), count_local(2), count_local(3), count_local(4), 2))
    if (.not. is_parallel) then
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,:,:,:,1), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,:,:,:,2), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,:,:,:,1), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,:,:,:,2), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
    end if
    array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4)) = &
      cmplx(buffer(:,:,:,:,1), buffer(:,:,:,:,2), kind=dp)

    deallocate (buffer)
    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_complex_double_4

  !  \brief read a 5 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_complex_double_5(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    complex(dp), intent(out) :: array(:,:,:,:,:)! data to be read
    real(dp), allocatable, dimension(:,:,:,:,:,:) :: buffer ! buffer to save real and imag part
    integer :: rank
    integer(HSIZE_T),dimension(5) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 5

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_5: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_complex_double_5("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_complex_double_5("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_complex_double_5("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_complex_double_5("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_complex_double_5("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_complex_double_5 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_complex_double_5 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    allocate (buffer(count_local(1), count_local(2), count_local(3), count_local(4), count_local(5), 2))
    if (.not. is_parallel) then
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,:,:,:,:,1), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,:,:,:,:,2), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,:,:,:,:,1), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,:,:,:,:,2), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
    end if
    array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4),1:count_local(5)) = &
      cmplx(buffer(:,:,:,:,:,1), buffer(:,:,:,:,:,2), kind=dp)

    deallocate (buffer)
    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_complex_double_5

  !  \brief read a 6 - dimension array from a hdf5 file
  subroutine hdf_read_dataset_complex_double_6(loc_id, dset_name, array, offset, non_parallel)

    integer(HID_T), intent(in) :: loc_id        ! local id in file
    character(len=*), intent(in) :: dset_name   ! name of dataset
    integer, optional, intent(in) :: offset(:)  ! offset of loading
    logical, optional, intent(in) :: non_parallel         ! load the data as it is
    complex(dp), intent(out) :: array(:,:,:,:,:,:)! data to be read
    real(dp), allocatable, dimension(:,:,:,:,:,:,:) :: buffer ! buffer to save real and imag part
    integer :: rank
    integer(HSIZE_T),dimension(6) :: dimsf, dimsm, offset_local, count_local
    integer(HID_T) :: dset_id, file_space_id, mem_space_id
    integer :: processor_write, axis_write
    integer(HSIZE_T), allocatable :: offset_glob(:), count_glob(:)
    integer :: hdferror, ii, jj
    logical :: is_parallel

    rank = 6

    if (hdf_print_messages) then
      write (*, '(A)') "--->hdf_read_dataset_complex_double_6: "//trim(dset_name)
    end if

    ! get file_space dimension
    call hdf_get_dims(loc_id, dset_name, dimsf)
    count_local = dimsf

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

    if (present(non_parallel) .and. rank > 0) then
      if (non_parallel) is_parallel = .false.
    end if

    ! allocate offset array
    if (is_parallel) then
      allocate(offset_glob(mpi_nrank), count_glob(mpi_nrank))
    end if

    ! syntax check of offset and set offset_glob / count_glob / count_local
    if (present(offset)) then
      if (.not. is_parallel) then
        write(*,'(A)') "hdf_read_dataset_complex_double_6("//trim(dset_name)//&
                       "): usless offset"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (size(offset) .ne. mpi_nrank) then
        write(*,'(A)') "hdf_read_dataset_complex_double_6("//trim(dset_name)//&
                       "): size of offset different than number of cores"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      if (offset(1) .ne. 0) then
        write(*, '(A)') "hdf_read_dataset_complex_double_6("//trim(dset_name)//&
                        "): offset(1) needs to be 0"
        call MPI_Abort(mpi_comm, mpi_ierr)
      end if

      do ii = 2, mpi_nrank
        if (offset(ii) < offset(ii-1) .or. offset(ii) > dimsf(axis_write)) then
          write(*,'(A)') "hdf_read_dataset_complex_double_6("//trim(dset_name)//&
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
        write(*, '(A)') "hdf_read_dataset_complex_double_6("//trim(dset_name)// &
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
          write(*, '(A)') "hdf_read_dataset_complex_double_6 ("//trim(dset_name)// &
                          "): array size is wrong"
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    else
      count_local(axis_write) = count_glob(mpi_irank+1)
      do ii = 1, rank
        jj = count_local(ii)
        if ((dimsm(ii) .lt. jj) .or. (ii .ne. axis_write .and. dimsm(ii) .ne. jj)) then
          write(*, '(A, I2, A, I2, A, I8)') "hdf_read_dataset_complex_double_6 ("//trim(dset_name)// &
                          "): array size is wrong axis=", ii, " memory siz=", dimsm(ii), " file size=", jj
          call MPI_Abort(mpi_comm, mpi_ierr)
        end if
      end do
    end if

    ! open dataset
    call h5dopen_f(loc_id, dset_name, dset_id, hdferror)

    ! read dataset
    allocate (buffer(count_local(1), count_local(2), count_local(3), count_local(4), count_local(5), count_local(6), 2))
    if (.not. is_parallel) then
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,:,:,:,:,:,1), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,:,:,:,:,:,2), &
                     dimsm, hdferror, xfer_prp=dplist_collective)
    else
      call h5screate_simple_f(rank, count_local, mem_space_id, hdferror)
      offset_local = 0
      offset_local(axis_write) = offset_glob(mpi_irank+1)
      call h5dget_space_f(dset_id, file_space_id, hdferror)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset_local, count_local, hdferror)
      call h5dread_f(dset_id, complexd_field_id(1), buffer(:,:,:,:,:,:,1), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
      call h5dread_f(dset_id, complexd_field_id(2), buffer(:,:,:,:,:,:,2), count_local, hdferror, &
                     mem_space_id=mem_space_id,    &
                     file_space_id=file_space_id,  &
                     xfer_prp=dplist_collective)
    end if
    array(1:count_local(1),1:count_local(2),1:count_local(3),1:count_local(4),1:count_local(5),1:count_local(6)) = &
      cmplx(buffer(:,:,:,:,:,:,1), buffer(:,:,:,:,:,:,2), kind=dp)

    deallocate (buffer)
    if (is_parallel) then
      deallocate(offset_glob, count_glob)
      call h5sclose_f(mem_space_id,  hdferror)
      call h5sclose_f(file_space_id, hdferror)
    end if

    ! close all id's
    call h5dclose_f(dset_id, hdferror)

  end subroutine hdf_read_dataset_complex_double_6

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
    call hdf_get_dims(loc_id, dset_name, dimsf)

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

    if (is_parallel .and. dimsf(1) .ne. mpi_nrank) then
      write(*, '(A)') "hdf_read_dataset_character_0("//trim(dset_name)// &
                        "): inconsistent number of elements v.s. number of processors"
      call MPI_Abort(mpi_comm, mpi_ierr)
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
        jj = int(dimsf(ii))
        if (ii == axis_write) jj = int(count_glob(mpi_irank+1))
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
      do ii = 1, int(dimsm(1))
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
        jj = int(dimsf(ii))
        if (ii == axis_write) jj = int(count_glob(mpi_irank+1))
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
      do ii = 1, int(dimsm(1))
        do jj = 1, int(dimsm(2))
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


  end subroutine hdf_read_dataset_character_2  !!---------------------------------------------------------------------------------------
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
