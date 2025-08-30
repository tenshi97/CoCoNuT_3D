#undef DEBUG
!> Module with generic functions for HDF5 I/O of Fortran data with metadata
module dataformat
  use precision
  use hdf5

  use mo_mpi
  use abort
  use setup_c

  implicit none
  public :: init_io, uninit_io, &
            create_data_file, open_data_file, close_data_file, flush_data_file, &
            writeq, readq, &
            hfile, pfile

  type :: datafile
    integer(kind=hid_t) :: fd
    integer(kind=hid_t) :: current_group
  end type

  type(datafile) :: hfile, pfile

! generated interfaces with perl magic
!PERL-START
!PERL my @directions = ("write", "read");
!PERL my @types = ("f8", "i4", "string", "logical");
!PERL foreach my $dir (@directions) {
  interface @[[$dir]]q
!PERL     foreach my $type (@types) {
    module procedure @[[$dir]]_@[[$type]]_scalar
!PERL         for (my $dim = 1; $dim < 8; $dim++) {
    module procedure @[[$dir]]_@[[$type]]_array_@[[$dim]]d
!PERL         }

!PERL     }
  end interface

!PERL }
!PERL-END

contains

  !> This has to be called before any output can take place,
  !> initializes the HDF5 library.
  !>
  !> \author lorenz
  !>
  subroutine init_io()
    implicit none
    integer :: error

    ! start hdf5
    call h5open_f(error)
    abort_if(error .lt. 0)

  end subroutine init_io

  !> This has to be called for properly closing the output infrastructure,
  !> finalizes the HDF5 library
  subroutine uninit_io()
    implicit none
    integer :: error

    call h5close_f(error)
    ! Too late to care about errors, so no "abort_if(error .lt. 0)"
  end subroutine uninit_io

  !> Create a new HDF5 output file, possibly MPI aware, for writing
  !>
  !> \author lorenz
  !>
  !> \param f           A type(datafile) structure
  !> \param filename    The filename for the new file
  !>
  subroutine create_data_file(f, filename)
    implicit none
    type(datafile), intent(out) :: f
    character(*), intent(in) :: filename

    integer(kind=hid_t) :: plist_id
    integer :: error, info

    f%fd = 0
    f%current_group = 0

#ifdef MPI_HYDRO
#ifndef ITASCA
    info = MPI_INFO_NULL
#else /* ITASCA */
    ! The standard settings for MPI-IO somehow do not work on
    ! the Itacsa cluster at MSI. You may try these setting
    ! when you encounter an error message like "File locking
    !failed in ADIOI_Set_lock" on other machines as well.
    call MPI_Info_create(info, error)

    ! Disables ROMIO's data-sieving */
    call MPI_Info_set(info, "romio_ds_read", "disable", error)
    call MPI_Info_set(info, "romio_ds_write", "disable", error)

    ! Enable ROMIO's collective buffering */
    call MPI_Info_set(info, "romio_cb_read", "enable", error)
    call MPI_Info_set(info, "romio_cb_write", "enable", error)
#endif /* ITASCA */
#endif /* MPI_HYDRO */

    ! Set up file access property list with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    abort_if(error .lt. 0)

#ifdef MPI_HYDRO
    ! mpi magic
      call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)
    abort_if(error .lt. 0)
#endif /* MPI_HYDRO */

    ! create a new file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, f%fd, error, access_prp=plist_id)
    abort_if(error .lt. 0)

    ! close the plist again, not needed anymore
    call h5pclose_f(plist_id, error)
    abort_if(error .lt. 0)

    call h5gopen_f(f%fd, "/", f%current_group, error)
    abort_if(error .lt. 0)

  end subroutine create_data_file

  !> Open an existing HDF5 input file, possibly MPI aware, for reading
  !>
  !> \author lorenz
  !>
  !> \param f           A type(datafile) structure
  !> \param filename    The filename for the new file
  !>
  subroutine open_data_file(f, filename)

    implicit none
    type(datafile), intent(out) :: f
    character(*), intent(in) :: filename

    integer(kind=hid_t) :: plist_id
    integer :: error

    f%fd = 0
    f%current_group = 0

    ! Set up file access property list with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    abort_if(error .lt. 0)

#ifdef MPI_HYDRO
      ! mpi magic
      call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
      abort_if(error .lt. 0)
#endif /* MPI_HYDRO */

    ! open a new file
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, f%fd, error, access_prp=plist_id)
    abort_if(error .lt. 0, "filename: " // filename)

    ! close the plist again, not needed anymore
    call h5pclose_f(plist_id, error)
    abort_if(error .lt. 0)

    call h5gopen_f(f%fd, "/", f%current_group, error)
    abort_if(error .lt. 0)

  end subroutine open_data_file

  !> Close an opened datafile
  !>
  !> \author lorenz
  !>
  !> \param f           A type(datafile) structure
  !>
  subroutine close_data_file(f)
    implicit none
    type(datafile), intent(in) :: f
    integer :: error

    ! close current time step
    call h5gclose_f(f%current_group, error)
    abort_if(error .lt. 0)

    call h5fclose_f(f%fd, error)
    abort_if(error .lt. 0)
  end subroutine close_data_file

  !> Flush an open datafile
  !>
  !> \author lorenz
  !>
  !> \param f           A type(datafile) structure
  !>
  subroutine flush_data_file(f)
    implicit none
    type(datafile), intent(in) :: f
    integer :: error

    call h5fflush_f(f%fd, H5F_SCOPE_GLOBAL_F, error)
    abort_if(error .lt. 0)
  end subroutine flush_data_file

  !> Check whether on object exists in a file
  !>
  !> \author lorenz
  !>
  !> \param     f       A type(datafile) structure
  !> \param     name    The name of the object
  !>
  !> \returns   exists  logical, true if there is such an object in the file
  function exists_in_file(f, name) result(exists)
    implicit none
    type(datafile), intent(in)   :: f
    character(len=*), intent(in) :: name
    logical :: exists
    integer :: error

    call h5lexists_f(f%current_group, name, exists, error)
    abort_if(error .lt. 0)
  end function exists_in_file

! The low level functions to write the output atoms
! are generic to every kind of variable,
! therefore only written once and reused
! with the help of perl magic;

!PERL-START
#ifdef HDF5_1_10  
!PERL   my @types = (["f8", "real(kind=8)",    "H5T_NATIVE_REAL_C_DOUBLE"   ],
!PERL                ["i4", "integer(kind=4)", "H5T_NATIVE_INTEGER_KIND(3)"],
!PERL                ["string", "character(len=*)", "dtype_id"],
!PERL                ["logical", "logical", "H5T_NATIVE_INTEGER_KIND(3)"]);
#else /* NEW_HDF5 */
!PERL   my @types = (["f8", "real(kind=8)",    "H5T_NATIVE_REAL_8"   ],
!PERL                ["i4", "integer(kind=4)", "H5T_NATIVE_INTEGER_4"],
!PERL                ["string", "character(len=*)", "dtype_id"],
!PERL                ["logical", "logical", "H5T_NATIVE_INTEGER_4"]);
#endif /* NEW_HDF5 */  
!PERL   foreach my $type (@types) {
!PERL       my $id = $$type[0];
!PERL       my $ftype = $$type[1];
!PERL       my $h5type = $$type[2];

  !> Attach an @[[$h5type]] (corresponding to fortran @[[$ftype]]) attribute to a dataset
  !>
  !> \param dset_id     dataset specifier
  !> \param name        name of the attribute
  !> \param value       content of the Attribute
  !>
  subroutine write_@[[$id]]_attribute(dset_id, name, value)
    use abort
    implicit none
    integer(kind=hid_t), intent(in) :: dset_id
    character(*), intent(in)        :: name
    @[[$ftype]], intent(in)         :: value
!PERL if ($id eq "logical") {
    integer(kind=4)                 :: value_int
!PERL }

    integer :: error
    integer(kind=hid_t) :: dspace_id, attr_id
!PERL if ($h5type eq "dtype_id") {
    integer(kind=hid_t) :: dtype_id
!PERL }

    ! Create scalar dataspace for the string attributes
    call h5screate_f(H5S_SCALAR_F, dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

!PERL if ($id eq "string") {
    ! Create the datatype for this attribute
    call h5tcopy_f(H5T_NATIVE_CHARACTER, dtype_id, error)
    abort_if(error .lt. 0)
    call h5tset_size_f(dtype_id, int(len(trim(value)), kind=size_t), error)
    abort_if(error .lt. 0)
!PERL }

    ! Create the attribute
    call h5acreate_f(dset_id, name, @[[$h5type]], dspace_id, attr_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    ! Write it
!PERL if ($id eq "logical") {
    if (value) then
      value_int = 1
    else
      value_int = 0
    end if
    call h5awrite_f(attr_id, @[[$h5type]], value_int, (/1_hsize_t/), error)
!PERL } else {
    call h5awrite_f(attr_id, @[[$h5type]], value, (/1_hsize_t/), error)
!PERL }
    abort_if(error .lt. 0, "data field: " // name)

    ! Release allocated objects

    call h5aclose_f(attr_id, error)
    abort_if(error .lt. 0, "data field: " // name)

!PERL if ($id eq "string") {
    call h5tclose_f(dtype_id, error)
    abort_if(error .lt. 0)

!PERL }
    call h5sclose_f(dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)
  end subroutine write_@[[$id]]_attribute

  !> Read a @[[$h5type]] (corresponding to fortran @[[$ftype]]) attribute of a dataset
  !>
  !> \param     f           an opened data file
  !> \param     dset_name   dataset name
  !> \param     name        attribute name
  !>
  !> \param     value       contains the contents of the attribute after return
  !>
  subroutine read_@[[$id]]_attribute(f, dset_name, name, value)
    implicit none
    type(datafile), intent(in)    :: f
    character(len=*), intent(in)  :: dset_name, name
    @[[$ftype]], intent(out)      :: value
!PERL if ($id eq "logical") {
    integer(kind=4)               :: value_int
!PERL }

    integer :: error
    integer(kind=hid_t)  :: dset_id, attr_id
!PERL if ($h5type eq "dtype_id") {
    integer(kind=hid_t)  :: dtype_id
    integer(kind=size_t) :: attribute_string_length
!PERL }

    ! Open the dataset
    call h5dopen_f(f%fd, dset_name, dset_id, error)
    abort_if(error .lt. 0, "data field: " // dset_name)

    ! Open the attribute
    call h5aopen_f(dset_id, name, attr_id, error)
    abort_if(error .lt. 0)

    ! Release the dataset
    call h5dclose_f(dset_id, error)
    abort_if(error .lt. 0)

!PERL if ($id eq "string") {
    ! Get its datatype
    call h5aget_type_f(attr_id, dtype_id, error)
    abort_if(error .lt. 0)

    ! string length
    call h5tget_size_f(dtype_id, attribute_string_length, error)
    abort_if(error .lt. 0)
    abort_if(attribute_string_length .gt. len(value))

    ! Release the attributes native datatype
    call h5tclose_f(dtype_id, error)
    abort_if(error .lt. 0)

    ! Create the datatype for the string buffer
    ! (the "memory datatype" in HDF5 lingo)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, dtype_id, error)
    abort_if(error .lt. 0)
    call h5tset_size_f(dtype_id, int(attribute_string_length, kind=size_t), error)
    abort_if(error .lt. 0)

    value = ""

!PERL }

    ! Read the attribute
!PERL if ($id eq "logical") {
    call h5aread_f(attr_id, @[[$h5type]], value_int, (/ 1_hsize_t /), error)
    value = value_int /= 0
!PERL } else {
    call h5aread_f(attr_id, @[[$h5type]], value, (/ 1_hsize_t /), error)
!PERL }
    abort_if(error .lt. 0)

    ! Release the allocated objects
    call h5aclose_f(attr_id, error)
    abort_if(error .lt. 0)

!PERL if ($id eq "string") {
    call h5tclose_f(dtype_id, error)
    abort_if(error .lt. 0)

!PERL }
  end subroutine read_@[[$id]]_attribute

!PERL   for (my $rank = 1; $rank < 8; $rank++) {
!PERL       my $gridlist   = join(", ", map {"grid" . $_} (1..$rank));
!PERL       my $dimensions = join(",", ((":") x $rank));
!PERL       my $data_shape = join(", &\n      ", map {"size(data, dim=" . $_ . ")"} (1..$rank));

  !> Instance subroutine for writing @[[$rank]] dimensional @[[$ftype]] arrays,
  !> should not be used explicitly but rather via writeq
  !>
  !> \author lorenz
  !>
  !> \param f           datafile structure describing the file and group to write to
  !> \param name        identifier for the quantity
  !> \param data        the actual array
  !> \param desc        textual description
  !> \param unit        parseable unit, optional
!PERL           for (my $i = 1; $i <= $rank; $i++) {
  !> \param grid@[[$i]]       identifier of the quantity's grid in dimension @[[$i]]
!PERL           }
  subroutine write_@[[$id]]_array_@[[$rank]]d(f, name, data, desc, unit, &
                      @[[$gridlist]], &
                      localextents, globalextents)
    implicit none
    @[[$ftype]], intent(in), dimension(@[[$dimensions]]) :: data
!PERL if ($id eq "logical") {
    integer(kind=4), dimension(@[[$dimensions]]) :: data_int(@[[$data_shape]])
!PERL }
    character(*), intent(in)                 :: @[[$gridlist]]
    integer(kind=ik), parameter              :: rank = @[[$rank]]
    type(datafile), intent(in)               :: f
    character(*), intent(in)                 :: name, desc
    character(*), intent(in), optional       :: unit
    integer(kind=ik), intent(in), optional   :: localextents(2 * rank), globalextents(2 * rank)

    logical :: write_parallel
    integer :: error, istat

    integer(kind=hsize_t)       :: globalshape(rank), localshape(rank), offset(rank)
    integer(kind=hid_t)         :: dspace_id, dset_id, memspace, plist_id

!PERL if ($h5type eq "dtype_id") {
    integer(kind=hid_t) :: dtype_id
!PERL }

    if (present(localextents) .and. present(globalextents)) then
      globalshape = globalextents(2::2) - globalextents(1::2) + 1
      localshape = localextents(2::2) - localextents(1::2) + 1
      offset = localextents(1::2) - globalextents(1::2)
      if (all(globalextents .eq. localextents)) then
        write_parallel = .false.
      else
        write_parallel = .true.
      endif
    elseif (.not. present(localextents) .and. .not. present(globalextents)) then
      globalshape = shape(data)
      localshape = shape(data)
      offset(:) = 0
      write_parallel = .false.
    else
      raise_abort("You have to specify localextents and globalextents or none! Data field: " // name)
    end if

#ifdef DEBUG
    write(*,'(a)') "write_@[[$id]]_array_@[[$rank]]d"
    write(*,*) "name = ", trim(name)
    write(*,*) "desc = ", trim(desc)
    if (present(unit)) then
      write(*,*) "unit = ", trim(unit)
    endif
!PERL           for (my $i = 1; $i <= $rank; $i++) {
    write(*,*) "grid@[[$i]] = ", trim(grid@[[$i]])
!PERL           }
    if (present(localextents)) then
    write(*,*) "global lbound = ", globalextents(1::2)
    write(*,*) "global ubound = ", globalextents(2::2)
    write(*,*) " local lbound = ", localextents(1::2)
    write(*,*) " local ubound = ", localextents(2::2)
    write(*,*) "    parallel? = ", .not. (globalextents .eq. localextents)
    endif
    write(*,*) " global shape = ", globalshape
    write(*,*) "  local shape = ", localshape
    write(*,*) "   data shape = ", shape(data)
    write(*,*) "       offset = ", offset
#endif

    abort_if(write_parallel .and. .not.(use_mpi), "data field: " // name)

    ! The array to write has to fit its specified shape
    abort_if(.not. all(shape(data) .eq. localshape), "data field: " // name)

    call h5screate_simple_f(rank, globalshape, dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

!PERL if ($id eq "string") {
    ! create the "memory datatype" for this string array
    call h5tcopy_f(H5T_NATIVE_CHARACTER, dtype_id, error)
    abort_if(error .lt. 0)
    call h5tset_size_f(dtype_id, int(len(data(@[[join(",", (("1") x $rank))]])), kind=size_t), error)
    abort_if(error .lt. 0)

!PERL }
    ! Create the dataset
    call h5dcreate_f(f%current_group, name, @[[$h5type]], dspace_id, dset_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    ! Release the dataspace
    call h5sclose_f(dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    ! define our dataset in memory
    call h5screate_simple_f(rank, localshape, memspace, error)
    abort_if(error .lt. 0, "data field: " // name)

    if ((.not. write_parallel) .and. (myproc .ne. 0)) then
      ! but shrink it to nothing in case we don't want to write
      call h5sselect_none_f(memspace, error)
      abort_if(error .lt. 0, "data field: " // name)
    endif

    ! Select hyperslab in the file.
    call h5dget_space_f(dset_id, dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, localshape, error)
    abort_if(error .lt. 0, "data field: " // name)

    if ((.not. write_parallel) .and. (myproc .ne. 0)) then
#ifdef DEBUG
      write(*,*) "Empty write in task ", myproc
#endif
      call h5sselect_none_f(dspace_id, error)
      abort_if(error .lt. 0, "data field: " // name)
    endif

    ! Create property list for collective dataset write.
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    abort_if(error .lt. 0, "data field: " // name)

#ifdef MPI_HYDRO
      !> \todo Investigate difference between H5FD_MPIO_INDEPENDENT_F and H5FD_MPIO_COLLECTIVE_F.
      !>       H5FD_MPIO_COLLECTIVE_F sounds faster (?), but for large arrays this can crash
      !>       due to memory requirements.
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
      abort_if(error .lt. 0, "data field: " // name)
#endif

    ! Do the write
!PERL if ($id eq "logical") {
    where (data)
      data_int = 1
    elsewhere
      data_int = 0
    endwhere
    call h5dwrite_f(dset_id, @[[$h5type]], data_int, localshape, error, memspace, dspace_id, plist_id)
!PERL } else {
    call h5dwrite_f(dset_id, @[[$h5type]], data, localshape, error, memspace, dspace_id, plist_id)
!PERL }
    abort_if(error .lt. 0, "data field: " // name)

    call write_string_attribute(dset_id, "description", desc)
    if (present(unit)) then
      call write_string_attribute(dset_id, "unit", unit)
    endif

!PERL           for (my $i = 1; $i <= $rank; $i++) {
    call write_string_attribute(dset_id, "grid@[[$i]]", grid@[[$i]])
!PERL           }

    ! Release allocated objects
    !> \todo check for memory leaks
    call h5dclose_f(dset_id, error)
    abort_if(error .lt. 0, "data field: " // name)

!PERL if ($id eq "string") {
    call h5tclose_f(dtype_id, error)
    abort_if(error .lt. 0, "data field: " // name)

!PERL }
    call h5sclose_f(dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    call h5sclose_f(memspace, error)
    abort_if(error .lt. 0, "data field: " // name)

    call h5pclose_f(plist_id, error)
    abort_if(error .lt. 0, "data field: " // name)

#ifdef DEBUG
    write(*,*)
#endif
  end subroutine write_@[[$id]]_array_@[[$rank]]d

!PERL   }

  !> Instance subroutine for writing a single @[[$ftype]] number,
  !> should not be used explicitly but rather via writeq
  !>
  !> \author lorenz
  !>
  !> \param f           datafile structure describing the file and group to write to
  !> \param name        identifier for the quantity
  !> \param data        the actual scalar
  !> \param desc        textual description
  !> \param unit        parseable unit, optional
  !>
  subroutine write_@[[$id]]_scalar(f, name, data, desc, unit)
    implicit none

    @[[$ftype]], intent(in)            :: data
!PERL if ($id eq "logical") {
    integer(kind=4)                    :: data_int
!PERL }
    type(datafile), intent(in)         :: f
    character(*), intent(in)           :: name, desc
    character(*), intent(in), optional :: unit

    integer :: error
    integer(kind=hid_t)         :: dspace_id, dset_id, plist_id
!PERL if ($h5type eq "dtype_id") {
    integer(kind=hid_t) :: dtype_id
!PERL }


#ifdef DEBUG
    write(*,'(a)') "write_@[[$id]]_scalar"
    write(*,*) "name = ", trim(name)
    write(*,*) "desc = ", trim(desc)
    if (present(unit)) then
      write(*,*) "unit = ", trim(unit)
    endif
    write(*,*)
#endif

    ! Create the dataspace
    call h5screate_f(H5S_SCALAR_F, dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

!PERL if ($id eq "string") {
    ! create the "memory datatype" for this string array
    call h5tcopy_f(H5T_NATIVE_CHARACTER, dtype_id, error)
    abort_if(error .lt. 0)
    call h5tset_size_f(dtype_id, int(len(data), kind=size_t), error)
    abort_if(error .lt. 0)

!PERL }
    ! Create the dataset
    call h5dcreate_f(f%current_group, name, @[[$h5type]], dspace_id, dset_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    if (myproc .ne. 0) then
      ! but shrink the dataspace to nothing in case we don't want to write
      call h5sselect_none_f(dspace_id, error)
      abort_if(error .lt. 0, "data field: " // name)
#ifdef DEBUG
      write(*,*) "Empty write in task ", myproc
#endif
    endif

    ! Create property list for collective dataset write.
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    abort_if(error .lt. 0, "data field: " // name)

#ifdef MPI_HYDRO
      !> \todo Investigate difference between H5FD_MPIO_INDEPENDENT_F and H5FD_MPIO_COLLECTIVE_F.
      !>       H5FD_MPIO_COLLECTIVE_F sounds faster (?), but for large arrays this can crash
      !>       due to memory requirements.
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
      abort_if(error .lt. 0, "data field: " // name)
#endif

    ! Do the write
!PERL if ($id eq "logical") {
    if (data) then
      data_int = 1
    else
      data_int = 0
    endif
    call h5dwrite_f(dset_id, @[[$h5type]], data_int, int((/0/), kind=hsize_t), error, xfer_prp=plist_id)
!PERL } else {
    call h5dwrite_f(dset_id, @[[$h5type]], data, int((/0/), kind=hsize_t), error, xfer_prp=plist_id)
!PERL }
    abort_if(error .lt. 0, "data field: " // name)

    call write_string_attribute(dset_id, "description", desc)
    if (present(unit)) then
      call write_string_attribute(dset_id, "unit", unit)
    endif

    ! Release allocated objects
    !> \todo check for memory leaks
    call h5dclose_f(dset_id, error)
    abort_if(error .lt. 0, "data field: " // name)

!PERL if ($id eq "string") {
    call h5tclose_f(dtype_id, error)
    abort_if(error .lt. 0, "data field: " // name)

!PERL }
    call h5sclose_f(dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    call h5pclose_f(plist_id, error)
    abort_if(error .lt. 0, "data field: " // name)

  end subroutine write_@[[$id]]_scalar

!PERL   for (my $rank = 1; $rank < 8; $rank++) {
!PERL       my $gridlist   = join(", ", map {"grid" . $_} (1..$rank));
!PERL       my $dimensions = join(",", ((":") x $rank));
!PERL       my $data_shape = join(", &\n      ", map {"size(data, dim=" . $_ . ")"} (1..$rank));

  !> Instance subroutine for reading @[[$rank]] dimensional @[[$ftype]] arrays,
  !> should not be used explicitly but rather via readq
  !>
  !> \author lorenz
  !>
  !> \param f           datafile structure describing the file and group to read from
  !> \param name        identifier for the quantity
  !> \param data        the actual array
  !> \param localextents  optional array specifying the locally accesible
  !>                      part of data, in a (/ lower_bound, upper_bound /) fashion
  !> \param globalextents likewise optional array specifying the global bounds of the quantity
  !>
  subroutine read_@[[$id]]_array_@[[$rank]]d(f, name, data, localextents, globalextents, allow_partial_read)
    implicit none

    @[[$ftype]], intent(out), dimension(@[[$dimensions]]) :: data
!PERL if ($id eq "logical") {
    integer(kind=4), dimension(@[[$dimensions]]) :: data_int(@[[$data_shape]])
!PERL }
    integer(kind=ik), parameter               :: rank = @[[$rank]]
    type(datafile), intent(in)                :: f
    character(*), intent(in)                  :: name
    integer(kind=ik), intent(in), optional    :: localextents(2 * rank), globalextents(2 * rank)
    logical, intent(in), optional             :: allow_partial_read

    logical :: dataspace_is_simple, datatype_is_equal, partial_read_ok
    integer :: error, read_rank
    integer(kind=hsize_t)       :: globalshape(rank), localshape(rank), offset(rank)
    integer(kind=hsize_t)       :: read_shape(rank), max_dims(rank)
    integer(kind=hid_t)         :: dspace_id, dset_id, memspace, plist_id
!PERL if ($h5type eq "dtype_id") {
    integer(kind=hid_t)         :: dtype_id
!PERL }
!PERL if ($id eq "string") {
    integer(kind=size_t)        :: string_length
!PERL }

    if (present(localextents) .and. present(globalextents)) then
      globalshape = globalextents(2::2) - globalextents(1::2) + 1
      localshape = localextents(2::2) - localextents(1::2) + 1
      offset = localextents(1::2) - globalextents(1::2)
    elseif (.not. present(localextents) .and. .not. present(globalextents)) then
      globalshape = shape(data)
      localshape = shape(data)
      offset(:) = 0
    else
      raise_abort("You have to specify localextents and globalextents or none! Data field: " // name)
    end if

    ! Region to read has to fit precisely into data
    abort_if(.not. all(shape(data) .eq. localshape), "data field: " // name)

#ifdef DEBUG
    write(*,'(a)') "read_@[[$id]]_array_@[[$rank]]d"
    write(*,*) "name = ", trim(name)
    write(*,*) "rank = ", rank
    if (present(localextents)) then
    write(*,*) "global lbound = ", globalextents(1::2)
    write(*,*) "global ubound = ", globalextents(2::2)
    write(*,*) " local lbound = ", localextents(1::2)
    write(*,*) " local ubound = ", localextents(2::2)
    endif
    write(*,*) " global shape = ", globalshape
    write(*,*) "  local shape = ", localshape
    write(*,*) "       offset = ", offset
#endif

    ! Open the dataset
    call h5dopen_f(f%fd, name, dset_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    ! Get its dataspace
    call h5dget_space_f(dset_id, dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    ! Make sanity checks about the shape
    call h5sis_simple_f(dspace_id, dataspace_is_simple, error)
    abort_if(error .lt. 0, "data field: " // name)
    abort_if(.not. dataspace_is_simple, "data field: " // name)

    call h5sget_simple_extent_ndims_f(dspace_id, read_rank, error)
    abort_if(error .lt. 0, "data field: " // name)
    abort_if(read_rank .ne. rank, "data field: " // name)

    call h5sget_simple_extent_dims_f(dspace_id, read_shape, max_dims, error)
#ifdef DEBUG
    write(*,*) "   read shape = ", read_shape
    write(*,*) "     max dims = ", max_dims
#endif
    abort_if(error .lt. 0, "data field: " // name)
    if (present(allow_partial_read)) then
      partial_read_ok = allow_partial_read
    else
      partial_read_ok = .false.
    endif
    if (.not. partial_read_ok) then
      abort_if(.not. all(read_shape .eq. globalshape), "data field: " // name)
      abort_if(.not. all(max_dims .eq. globalshape), "data field: " // name)
    endif

    ! Release the dataspace
    call h5sclose_f(dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    ! Select hyperslab in the file.
    call h5screate_simple_f(rank, localshape, memspace, error)
    abort_if(error .lt. 0, "data field: " // name)

    call h5dget_space_f(dset_id, dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, localshape, error)
    abort_if(error .lt. 0, "data field: " // name)

    ! Create property list for collective dataset write.
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    abort_if(error .lt. 0, "data field: " // name)

#ifdef MPI_HYDRO
      !> \todo Investigate difference between H5FD_MPIO_INDEPENDENT_F and H5FD_MPIO_COLLECTIVE_F.
      !>       H5FD_MPIO_COLLECTIVE_F sounds faster (?), but for large arrays this can crash
      !>       due to memory requirements.
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
      abort_if(error .lt. 0, "data field: " // name)
#endif

!PERL if ($id eq "string") {
    ! Get its datatype
    call h5dget_type_f(dset_id, dtype_id, error)
    abort_if(error .lt. 0)

    ! string length
    call h5tget_size_f(dtype_id, string_length, error)
    abort_if(error .lt. 0)
    abort_if(string_length .gt. len(data(@[[join(",", (("1") x $rank))]])))

    ! Release native datatype
    call h5tclose_f(dtype_id, error)
    abort_if(error .lt. 0)

    ! Create the datatype for the string buffer
    ! (the "memory datatype" in HDF5 lingo)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, dtype_id, error)
    abort_if(error .lt. 0)
    call h5tset_size_f(dtype_id, int(string_length, kind=size_t), error)
    abort_if(error .lt. 0)

    ! initilize with spaces in case the scalar length is larger than in the file
    data = ""
!PERL }

    ! Do the read
!PERL if ($id eq "logical") {
    call h5dread_f(dset_id, @[[$h5type]], data_int, localshape, error, memspace, dspace_id, plist_id)
    data = data_int /= 0
!PERL } else {
    call h5dread_f(dset_id, @[[$h5type]], data, localshape, error, memspace, dspace_id, plist_id)
!PERL }
    abort_if(error .lt. 0, "data field: " // name)

    ! Release allocated objects
    !> \todo check for memory leaks
    call h5dclose_f(dset_id, error)
    abort_if(error .lt. 0, "data field: " // name)

!PERL if ($id eq "string") {
    call h5tclose_f(dtype_id, error)
    abort_if(error .lt. 0, "data field: " // name)

!PERL }
    call h5sclose_f(dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    call h5sclose_f(memspace, error)
    abort_if(error .lt. 0, "data field: " // name)

    call h5pclose_f(plist_id, error)
    abort_if(error .lt. 0, "data field: " // name)

#ifdef DEBUG
    write(*,*)
#endif

  end subroutine read_@[[$id]]_array_@[[$rank]]d

!PERL   }

  !> Instance subroutine for reading a single @[[$ftype]] number,
  !> should not be used explicitly but rather via readq
  !>
  !> \author lorenz
  !>
  !> \param f           datafile structure describing the file and group to read from
  !> \param name        identifier for the quantity
  !> \param data        the actual scalar
  !>
  subroutine read_@[[$id]]_scalar(f, name, data)
    implicit none

    type(datafile), intent(in)           :: f
    @[[$ftype]], intent(out)             :: data
!PERL if ($id eq "logical") {
    integer(kind=4)                      :: data_int
!PERL }
    character(*), intent(in)             :: name
    integer :: error, read_rank
    integer(kind=hid_t)         :: dspace_id, dset_id, plist_id
    logical :: dataspace_is_simple
!PERL if ($h5type eq "dtype_id") {
    integer(kind=hid_t)                  :: dtype_id
!PERL }

#ifdef DEBUG
    write(*,'(a)') "read_@[[$id]]_scalar"
    write(*,*) "name = ", trim(name)
    write(*,*)
#endif

    ! Open the dataset
    call h5dopen_f(f%fd, name, dset_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    ! Get its dataspace
    call h5dget_space_f(dset_id, dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    ! Make sanity checks about the shape
    call h5sis_simple_f(dspace_id, dataspace_is_simple, error)
    abort_if(error .lt. 0, "data field: " // name)
    abort_if(.not. dataspace_is_simple, "data field: " // name)

    call h5sget_simple_extent_ndims_f(dspace_id, read_rank, error)
    abort_if(error .lt. 0, "data field: " // name)
    abort_if(read_rank .ne. 0, "data field: " // name)

    ! Create property list for collective dataset write.
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    abort_if(error .lt. 0, "data field: " // name)

#ifdef MPI_HYDRO
      !> \todo Investigate difference between H5FD_MPIO_INDEPENDENT_F and H5FD_MPIO_COLLECTIVE_F.
      !>       H5FD_MPIO_COLLECTIVE_F sounds faster (?), but for large arrays this can crash
      !>       due to memory requirements.
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
      abort_if(error .lt. 0, "data field: " // name)
#endif

!PERL if ($id eq "string") {
    ! Create the datatype for the string buffer
    ! (the "memory datatype" in HDF5 lingo)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, dtype_id, error)
    abort_if(error .lt. 0)
    call h5tset_size_f(dtype_id, int(len(data), kind=size_t), error)
    abort_if(error .lt. 0)

!PERL }

    ! Do the read
!PERL if ($id eq "logical") {
    call h5dread_f(dset_id, @[[$h5type]], data_int, int((/0/), kind=hsize_t), error, xfer_prp=plist_id)
    data = data_int /= 0
!PERL } else {
    call h5dread_f(dset_id, @[[$h5type]], data, int((/0/), kind=hsize_t), error, xfer_prp=plist_id)
!PERL }
    abort_if(error .lt. 0, "data field: " // name)

    ! Release allocated objects
    !> \todo check for memory leaks
    call h5dclose_f(dset_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    call h5sclose_f(dspace_id, error)
    abort_if(error .lt. 0, "data field: " // name)

    call h5pclose_f(plist_id, error)
    abort_if(error .lt. 0, "data field: " // name)

  end subroutine read_@[[$id]]_scalar

!PERL-END }

end module dataformat
