! ----------------------------------------------
! Function template for hdf_read_dataset_XXX_YYY
! You need to define FUNCP, PDTYPE, PDRANK and PDH5DTYPE before include
! ----------------------------------------------
#define GLUE2_HELPER(a,b)  a##b
#define GLUE2(a,b)  GLUE2_HELPER(a,b)
#define GLUE3_HELPER(a,b,c)  a##b##c
#define GLUE3(a,b,c) GLUE3_HELPER(a,b,c)
#define GLUE4_HELPER(a,b,c,d)  a##b##c##d
#define GLUE4(a,b,c,d) GLUE4_HELPER(a,b,c,d)
#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)

#define FUNC          GLUE3(FUNCP, _, PDRANK)

subroutine FUNC(loc_id, dset_name, array, chunks, filter)

  integer(HID_T), intent(in) :: loc_id        ! local id in file
  character(len=*), intent(in) :: dset_name   ! name of dataset
#if PDRANK == 0
  PDTYPE, intent(in) :: array                  ! data to be written
#elif PDRANK == 1
  PDTYPE, intent(in) :: array(:)               ! data to be written
#elif PDRANK == 2
  PDTYPE, intent(in) :: array(:,:)             ! data to be written
#elif PDRANK == 3
  PDTYPE, intent(in) :: array(:,:,:)           ! data to be written
#elif PDRANK == 4
  PDTYPE, intent(in) :: array(:,:,:,:)         ! data to be written
#elif PDRANK == 5
  PDTYPE, intent(in) :: array(:,:,:,:,:)       ! data to be written
#elif PDRANK == 6
  PDTYPE, intent(in) :: array(:,:,:,:,:,:)     ! data to be written
#endif

  integer :: rank
#if PDRANK > 0
  integer, optional, intent(in) :: chunks(PDRANK)        ! chunk size for dataset
  character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle')
  integer(SIZE_T) :: dims(PDRANK), cdims(PDRANK)
#else
  integer, optional, intent(in) :: chunks           ! chunk size for dataset
  character(len=*), optional, intent(in) :: filter  ! filter to use ('none', 'szip', 'gzip', 'gzip+shuffle'
  integer(SIZE_T) :: dims(1)
#endif

  integer(HID_T) :: dset_id, dspace_id, plist_id
  character(len=32) :: filter_case
  integer :: hdferror

  if (hdf_print_messages) then
      write(*,'(4A)') "--->", STRINGIFY(FUNC), ": " , trim(dset_name)
  end if

  ! set rank and dims
  rank = PDRANK
#if PDRANK == 0
  dims = (/ 0 /)

  ! create dataspace
  call h5screate_f(H5S_SCALAR_F, dspace_id, hdferror)

  ! create dataset
  call h5dcreate_f(loc_id, dset_name, PDH5DTYPE, dspace_id, dset_id, hdferror)

  !write(*,'(A20,I0)') "h5dcreate: ", hdferror
#else
  dims = shape(array, KIND=HID_T)
  !
  if (present(filter)) then
    filter_case = filter
  else
    filter_case = hdf_default_filter
  endif

  ! set chunk (if needed)
  if (present(chunks)) then
    cdims = int(chunks, SIZE_T)
  else
    cdims = 0
  endif

  ! create and set property list
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, hdferror)
  call hdf_set_property_list(plist_id, rank, dims, cdims, filter_case)

  ! create dataspace
  call h5screate_simple_f(rank, dims, dspace_id, hdferror)
  !write(*,'(A20,I0)') "h5screate_simple: ", hdferror

  ! create dataset
  call h5dcreate_f(loc_id, dset_name, PDH5DTYPE, dspace_id, dset_id, hdferror, dcpl_id=plist_id)
  !write(*,'(A20,I0)') "h5dcreate: ", hdferror
#endif


  ! write dataset
  call h5dwrite_f(dset_id, PDH5DTYPE, array, dims, hdferror)
  !write(*,'(A20,I0)') "h5dwrite: ", hdferror

  ! close all ids
  call h5sclose_f(dspace_id, hdferror)
#if PDRANK > 0
  call h5pclose_f(plist_id, hdferror)
#endif
  call h5dclose_f(dset_id, hdferror)

end subroutine FUNC