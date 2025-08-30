module testdumper
  use precision

  implicit none

  public dump_
  private

  save
  integer(kind=ik) :: counter = 0


  contains

  function check() result(res)
    use iso_c_binding
    use environment

    implicit none

    logical :: res
    character(c_char), pointer :: val(:)

    val => getenv("CHECK")

    if (associated(val) .and. size(val) .gt. 0 .and. val(1) .eq. "1") then
      res = .true.
    else
      res = .false.
    endif
  end function

  subroutine dump_(name, location, array)
    use abort
    use error_screen
    use precision
    implicit none
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: location
    real(kind=rk), intent(in)    :: array(:)

    real(kind=rk)                :: read_array(size(array))
    character(len=len(name))     :: name_read
    integer(kind=ik)             :: counter_read, n
    integer :: fd
#ifdef WITH_OPENMP
    logical :: omp_in_parallel
    logical :: omp_in_parallel_read

    integer :: omp_get_thread_num, omp_get_num_threads
    integer :: omp_get_thread_num_read, omp_get_num_threads_read

    if (omp_in_parallel()) then
      fd = 1001 + omp_get_thread_num()
    else
      fd = 1000
    endif
#else
    fd = 1000
#endif

    counter = counter + 1
    if (.not. check()) then
      n = size(array)
      write(*,'(a,a,a,i0,a,i0)') "Dumping ",name,", size = ", n, ", counter = ", counter
#ifndef WITH_OPENMP
      write(fd) name, counter, n
#else
      write(fd) name, counter, n, omp_in_parallel(), omp_get_thread_num(), omp_get_num_threads()
#endif
      write(fd) array

    else
      write(*,'(a,a,a,i0,a,i0)') "Reading ",name,", size = ", size(array), ", counter = ", counter

#ifndef WITH_OPENMP
      read(fd) name_read, counter_read, n
#else
      read(fd) name_read, counter_read, n, omp_in_parallel_read, omp_get_thread_num_read, omp_get_num_threads_read
      abort_if(omp_in_parallel_read .neqv. omp_in_parallel(), name // ", " // location)
      abort_if(omp_get_thread_num_read .ne. omp_get_thread_num(), name // ", " // location)
      abort_if(omp_get_num_threads_read .ne. omp_get_num_threads(), name // ", " // location)
#endif
      abort_if(name_read .ne. name, name // ", " // location)
      abort_if(counter_read .ne. counter, name // ", " // location)
      abort_if(n .ne. size(array), name // ", " // location)

      read(fd) read_array
      if (.not. all(array .eq. read_array)) then
        call show_error_screen(size(array), array, read_array, array - read_array)
        raise_abort("Read-in array '"//name //"' differs from reference value in " // location)
      endif

    endif
  end subroutine dump_

end module testdumper
