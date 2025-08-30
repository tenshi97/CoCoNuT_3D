#ifdef WRITE_BINARY_OUTPUT
!> Module with code specific functions for the HDF5 I/O
module dataformat_vertex

  use dataformat
  implicit none

  public

  contains

  subroutine open_vertex_file(f, filename)
    implicit none
    type(datafile), intent(out) :: f
    character(*), intent(in) :: filename

    call open_data_file(f, filename)

  end subroutine open_vertex_file

  subroutine create_vertex_file(f, filename, description)
    use configure
    use iso_c_binding
    implicit none
    type(datafile), intent(out) :: f
    character(*), intent(in) :: filename, description
    integer(kind=hid_t), target :: group_id
    integer :: error

    call create_data_file(f, filename)

#if 0
    ! write some information in a sub-group "info"
    call h5gcreate_f(f%fd, "info", group_id, error)
    abort_if(error .lt. 0)

    call write_string_attribute(group_id, "description", description)
    call version_write_hdf5(c_loc(group_id))

    call h5gclose_f(group_id, error)

    ! write the current configuration values in a sub-group "config"
    ! save the current group
    group_id = f%current_group

    call h5gcreate_f(f%fd, "config", f%current_group, error)
    abort_if(error .lt. 0)

    call write_config(f)

    call h5gclose_f(f%current_group, error)
    abort_if(error .lt. 0)

    ! restore
    f%current_group = group_id
#endif

  end subroutine create_vertex_file

  !> Call this before the data of any time step,
  !> creates a new sub-group in the HDF5 tree and
  !> writes the field "time" with the current physical time in it.
  !>
  !> \author lorenz
  !>
  !> \param f   A type(datafile) structure
  !>
  subroutine write_timestep_marker(f)
    use gfloat_hy, only: time
    use intgrs_hy, only : nstep
    type(datafile), intent(inout) :: f
    character(len=4 + 1 + 12) :: groupname
    integer :: error, pos
    integer(kind=hsize_t), parameter :: one = 1

    call h5gclose_f(f%current_group, error)
    abort_if(error .lt. 0)

    write(groupname,'(''Step#'',i012)') nstep

    pos = scan(groupname, " ")
    do while(pos > 0)
      groupname(pos:pos) = "0"
      pos = scan(groupname, " ")
    end do

    ! create group
    call h5gcreate_f(f%fd, groupname, f%current_group, error)
    abort_if(error .lt. 0)

    call writeq(f, "time", time, "physical time", "s")
  end subroutine write_timestep_marker


  !> Handy routine to write a typical radiation quantity.
  !> Based on the shape of q, the grids are inferred.
  !>
  !> \param q 5D radiation quantity
  !> \param name identifier for this quantity
  !> \param desc textual description
  !> \param unit parseable unit
  !>
  subroutine write_5d_radiation_quantity(name, q, desc, unit)
    use precision
    use configure
    use abort

    implicit none

    real(kind=rk), intent(in) :: q(:,:,:,:,:)
    character(*), intent(in) :: name, desc, unit
    character(11) :: egrid
    character(12) :: rgrid

    ! transport domain, start, end
    !> \todo these should be replaced when the mpi-cleanup branch is merged back
    integer(kind=ik) :: ys, ye
    integer(kind=ik) :: zs, ze

    ! transport domain
#ifdef MPI_HYDRO
    ys = nymoms
    ye = nymome
    zs = nzmoms
    ze = nzmome
#else /* no MPI */
    ys = config%nystrt
    ye = config%nymom
    zs = 1
    ze = config%nztra
#endif

    if(size(q, dim=2) .eq. config%iemax) then
      egrid = 'energy:emid'
    elseif (size(q, dim=2) .eq. config%iemax+1) then
      egrid = 'energy:ener'
    else
      write(0,*) "write_5d_radiation_quantity(): size(", name, ", dim=2) = ", &
                  size(q, dim=2), ", config%iemax = ", config%iemax
      raise_abort("write_5d_radiation_quantity(): Invalid extents of radiation quantity energy dimension")
    endif

    if(size(q, dim=1) .eq. config%imaxp + 3) then
      rgrid = 'radius:ralag'
    elseif (size(q, dim=1) .eq. config%imaxp + 2) then
      rgrid = 'radius:rqlag'
    else
      write(0,*) "write_5d_radiation_quantity(): size(", name, ", dim=1) = ", &
                  size(q, dim=1), ", imaxp1 = ", config%imaxp + 1
      raise_abort("write_5d_radiation_quantity(): Invalid extents of radiation quantity radial dimension")
    endif

    if (config%nystrt .eq. 0) then
      ! First index of y-dimension contains angular average, output two separate quantities
      ! no MPI in this case

      call writeq(pfile, name // '_ave', q(:,:,:,1,1), desc, unit, &
        rgrid, egrid, 'species:'//specieslist())

      if (ubound(q, dim=4) > 1) then
        call writeq(pfile, name, q(:,:,:,2:,:), desc, unit, &
          rgrid, egrid, 'species:'//specieslist(), 'theta:yzn', 'phi:zzn')
      endif

    else
      ! No angular averages present, output as is
      ! MPI possible

      call writeq(pfile, name, q, desc, unit, &
        rgrid, egrid, 'species:'//specieslist(), 'theta:yzn', 'phi:zzn', &
        (/1,size(q, dim=1), 1,size(q, dim=2), 1,config%isma, ys,ye, zs,ze/), &
        (/1,size(q, dim=1), 1,size(q, dim=2), 1,config%isma, 1,config%nymom, 1,config%nztra/))
    endif

  end subroutine write_5d_radiation_quantity

  ! Produce a string "frange(i,j,k)"
  function frange(i,j,k) result(string)
    integer, intent(in)           :: i
    integer, intent(in), optional :: j,k
    character(128)                  :: string

    if(present(j) .and. present(k)) then
      write(string,'("frange(", i5, ", ", i5, ", ", i5, ")")') i,j,k
    elseif (present(j)) then
      write(string,'("frange(", i5, ", ", i5, ")")') i,j
    else
      write(string,'("frange(", i5, ")")') i
    endif
  end function frange

  ! Produce a string list out of an array of strings
  function stringlist(a) result(string)
    character(*), intent(in) :: a(:)
    character(64)            :: format
    character(1024)          :: string
    write(format, '(a,i3,a,a,a)') "('[""',", size(a) - 1, "(a, '"", ""'), a, '""]')"
    write(string, format) a
  end function stringlist

  ! Produce a list of the neutrino species
  function specieslist() result(string)
    use configure

    implicit none

    character(32), save :: one = '["nu_e"]', two = '["nu_e", "nubar_e"]', three = '["nu_e", "nubar_e", "nu_x"]'
    character(32) :: string

    select case(config%isma)
      case(1)
        string = one
      case(2)
        string = two
      case(3)
        string = three
      case default
    end select

  end function specieslist

end module dataformat_vertex
#endif
