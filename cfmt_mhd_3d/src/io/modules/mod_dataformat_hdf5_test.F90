#ifdef WRITE_BINARY_OUTPUT
!make NO_CONFIG=1 EXTRA_CPPFLAGS="-DUNIT_TESTS -DMPI_HYDRO" DEBUG=HIGH ENABLE=MPI
program test_hdf5_io
  use dataformat
  use abort
  implicit none
  type(datafile) :: f
  character(len=200) :: string
  real(kind=rk), allocatable, dimension(:)   :: a, a2
  real(kind=rk), allocatable, dimension(:,:) :: b, b2
  integer :: i, j, error
  integer, parameter :: N = 1024

  allocate(a(N),a2(N))
  allocate(b(N,N),b2(N,N))

  a = [(real(i, kind=rk), i = 1, N)]
  b = reshape([((real(a(i)**2, kind=rk), i = 1, N), j = 1,N)], [N, N])

  a2 = 0.0_rk
  b2 = 0.0_rk

  call mpi_init(error)
  call mpi_comm_size(MPI_COMM_WORLD,nprocs,error)
  call mpi_comm_rank(MPI_COMM_WORLD,myproc,error)
  use_mpi = .true.
  call init_io

  ! create some data
  write(*,*) "Test file creation.."
  call create_data_file(f, "test.h5")

  ! shared writes
  write(string,'("range(", i5, ")")') size(a)
  call writeq(f, "a", a, "some quantity", "cm", "zone:"//trim(string))
  call writeq(f, "b", b, "some other quantity", "cm**2", "x:a", "y:a")

  ! a collective write
  ! assume b is part of a larger global array
  write(string,'("range(", i5, ")")') size(b, dim=2) * nprocs
  call writeq(f, "global", b, "a global array", "cm**2", "x:a", "y:"//trim(string), &
      (/ 1, N, N * myproc + 1, N * (myproc + 1) /), &
      (/ 1, N,              1,       N * nprocs /))


  call close_data_file(f)

  ! read back some data
  write(*,*) "Test file reading.."
  call open_data_file(f, "test.h5")

  call readq(f, "b", b2)
  abort_if(.not. all(b == b2))

  call read_string_attribute(f, "b", "grid1", string)
  call readq(f, string(index(string, ":") + 1:), a2)
  abort_if(.not. all(a == a2))

  call close_data_file(f)

  ! test an invalid read
  write(*,*) "Test reading a non-existing object (expect many errors following)"
  call expect_abort()
  call readq(f, "some_invalid_thing", a2)
  abort_if(.not. got_abort())
  write(*,*) "... success!"

  call uninit_io
  call mpi_finalize(error)

end program test_hdf5_io
#endif
