module abort
  implicit none

  public abort_if_x, raise_abort_x
#ifdef UNCLEAN_EXIT_FILE
  public create_unclean_exit_file
#endif

#ifdef UNIT_TESTS
  public expect_abort, got_abort
  private catch_abort, raised_aborts
#endif

  save

#ifdef UNIT_TESTS
  logical :: catch_abort = .false.
  integer :: raised_aborts = 0
#endif

  contains

#ifdef UNIT_TESTS

  !> For unit test, we also want to test whether some
  !> piece of code correctly throws an abort, by calling
  !> this routine, a flag is set that the next abort is not
  !> fatal, but expected
  subroutine expect_abort()
    implicit none
    catch_abort = .true.
    raised_aborts = 0
  end subroutine

  !> After setting expect_abort(), this checks if an abort
  !> was raised in the meantime and clears the flag to
  !> ignore aborts
  function got_abort() result(got)
    implicit none
    logical :: got
    catch_abort = .false.
    got = raised_aborts > 0
  end function

#endif /* UNIT_TESTS */


  !> Abort execution immediately if condition is true,
  !> without saving any files. Use not directly, but via
  !> pre- and postprocessor magic as "abort_if" without the
  !> parameters location and condition_stringed which will
  !> be inserted automatically
  !>
  !> \param location            "file.f90:linenumber"
  !> \param condition           if true exit program otherwise continue normally
  !> \param condition_stringed  stringified version of condition
  !> \param additional_message  optional description of the problem and location
  subroutine abort_if_x(location, condition, condition_stringed, additional_message)
#ifdef BACKTRACE_ON_ABORT
    use debug_c, only : print_backtrace
#endif

    implicit none

    logical, intent(in) :: condition
    character*(*), intent(in) :: location, condition_stringed
    character*(*), intent(in), optional :: additional_message

    if (.not. condition) return

#ifdef UNIT_TESTS
    if (catch_abort) then
      raised_aborts = raised_aborts + 1
      return
    endif
#endif


    print *,"####"
    print *,"##"
    print *,"## ABORT at ", trim(location)
    print *,"##"
    print *,"## Cause of abort: ", condition_stringed
    if (present(additional_message)) then
      print *,"## ", additional_message
    endif
    print *,"##"
    print *,"####"

#ifdef UNCLEAN_EXIT_FILE
    if (present(additional_message)) then
      call create_unclean_exit_file("abort", location, condition_stringed // ", " // additional_message)
    else
      call create_unclean_exit_file("abort", location, condition_stringed)
    endif
#endif

#ifdef BACKTRACE_ON_ABORT
    call print_backtrace
#endif

    ! does not return
    stop

  end subroutine abort_if_x

  !> Abort execution immediately without saving any files.
  !> Use not directly, but via pre- and postprocessor magic
  !> as "raise_abort" without the parameters location and
  !> condition_stringed, which will be inserted automatically
  !>
  !> \param location            "file.f90:linenumber"
  !> \param message             description of the problem and location
  subroutine raise_abort_x(location, message)
#ifdef BACKTRACE_ON_ABORT
    use debug_c, only : print_backtrace
#endif
    implicit none

    character*(*), intent(in) :: location
    character*(*), intent(in), optional :: message

#ifdef UNIT_TESTS
    if (catch_abort) then
      raised_aborts = raised_aborts + 1
      return
    endif
#endif

    print *,"####"
    print *,"##"
    print *,"## ABORT at ", trim(location)
    if (present(message)) then
      print *,"##"
      print *,"## ", message
    endif
    print *,"##"
    print *,"####"

#ifdef UNCLEAN_EXIT_FILE
    call create_unclean_exit_file("abort", location, message)
#endif
#ifdef BACKTRACE_ON_ABORT
    call print_backtrace
#endif

    ! does not return
    stop

  end subroutine raise_abort_x

#ifdef UNCLEAN_EXIT_FILE
  subroutine create_unclean_exit_file(type, location, message)
    character*(*), intent(in) :: type, location, message
    open(unit=99, file="vertex.unclean_exit")
    write(99,*) type
    write(99,*) location
    write(99,*) message
    close(99)
  end subroutine
#endif

end module
