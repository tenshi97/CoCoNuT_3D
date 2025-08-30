!> \verbatim
!> This module provides the functions error and show_error_screen,
!> which are used to print error information and to halt the execution
!> if necessary
!>
!>  Author: A. Marek, MPA, March 2009
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
module error
  use precision

  public :: error_if_x, raise_error_x

  contains

  subroutine save_files(no_savare)
    use savare_overload
#ifdef WRITE_BINARY_OUTPUT
    use output_hydro, only : write_output_files, close_files, set_filenames
    use restart, only : restrt
#endif

    use configure
    implicit none
    logical, intent(in) :: no_savare
    real(kind=rk)       :: savare_self(2), savare_children(2)
#ifdef WRITE_BINARY_OUTPUT
    if (.not.(config%use_deactivate_output)) then
       ! write output files
       call write_output_files
       call close_files

    ! increment suffix
       call set_filenames(99,99)

       if (no_savare) then
          print *,   ' W A R N I N G: savare has not been called => restart incomplete, use only for diagnosis'
          write(0,*) ' W A R N I N G: savare has not been called => restart incomplete, use only for diagnosis'
       else
          ! recall beginning of hydro-sequence
          call savare(0_ik, savare_self, savare_children)
       end if

       ! then write out restart file
       call restrt(0_ik)
    
    else
       write(0,*) "Not writing any debug output, running with config%use_deactivate_output"
    endif
#endif

  end subroutine

  !> Abort execution after saving output and restart files for diagnosis
  !> if parameter condition is true
  !>
  !> Use not directly, but via pre- and postprocessor magic as "error_if"
  !> without the parameters location and condition_stringed, which will
  !> be inserted automatically
  !>
  !> \param location            "file.f90:linenumber"
  !> \param condition           if true exit program otherwise continue normally
  !> \param condition_stringed  stringified version of condition
  !> \param additional_message  optional description of the problem and location
  !> \param no_savare           if true, do not call savare before writing output,
  !>                            therefore producing output unsuitable for restarting
  !>                            but possibly more useful for debugging
  !>
  !> \todo this function should be able to inline
  !>
  subroutine error_if_x(location, condition, condition_stringed, &
                        additional_message, no_savare)
#ifdef UNCLEAN_EXIT_FILE
    use abort, only : create_unclean_exit_file
#endif
#ifdef BACKTRACE_ON_ABORT
    use debug_c, only : print_backtrace
#endif

    implicit none

    logical, intent(in)                 :: condition
    character*(*), intent(in)           :: condition_stringed, location
    character*(*), intent(in), optional :: additional_message
    logical, intent(in), optional       :: no_savare

    logical :: no_savare_actual

    if (.not. condition) return

    if (present(no_savare)) then
      no_savare_actual = no_savare
    else
      no_savare_actual = .false.
    endif

    print *,"####"
    print *,"##"
    print *,"## ERROR at ", trim(location)
    print *,"##"
    print *,"## Cause of exit: ", condition_stringed
    if (present(additional_message)) then
    print *,"## ", additional_message
    endif
    print *,"##"
    print *,"####"

    call save_files(no_savare_actual)

#ifdef UNCLEAN_EXIT_FILE
    if (present(additional_message)) then
      call create_unclean_exit_file("error", location, condition_stringed // ", " // additional_message)
    else
      call create_unclean_exit_file("error", location, condition_stringed)
    endif
#endif

#ifdef BACKTRACE_ON_ABORT
    call print_backtrace
#endif

    ! does not return
    stop

  end subroutine error_if_x

  !> Abort execution after saving output and restart files for diagnosis
  !>
  !> Use not directly, but via pre- and postprocessor magic as "raise_error"
  !> without the parameter location, which will be inserted automatically
  !>
  !> \param location            "file.f90:linenumber"
  !> \param message  optional description of the problem and location
  !> \param no_savare           if true, do not call savare before writing output,
  !>                            therefore producing output unsuitable for restarting
  !>                            but possibly more useful for debugging
  !>
  !> \todo this function should be able to inline
  !>
  subroutine raise_error_x(location, message, no_savare)
#ifdef BACKTRACE_ON_ABORT
    use debug_c, only : print_backtrace
#endif
#ifdef UNCLEAN_EXIT_FILE
    use abort, only : create_unclean_exit_file
#endif
    implicit none

    character*(*), intent(in)     :: location
    character*(*), intent(in)     :: message
    logical, intent(in), optional :: no_savare

    logical :: no_savare_actual

    if (present(no_savare)) then
      no_savare_actual = no_savare
    else
      no_savare_actual = .false.
    endif

    print *,"####"
    print *,"##"
    print *,"## ERROR at ", trim(location)
    print *,"##"
    print *,"## ", message
    print *,"##"
    print *,"####"

    call save_files(no_savare_actual)

#ifdef UNCLEAN_EXIT_FILE
    call create_unclean_exit_file("error", location, message)
#endif
#ifdef BACKTRACE_ON_ABORT
    call print_backtrace
#endif

    ! does not return
    stop

  end subroutine raise_error_x

end module error

!>
!> \verbatim
!> This module provides the variables which are necessary for
!> determing with a security file whether code is allowed to
!> start or not
!>
!>  Author: A. Marek, MPA, March 2009
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
module io_security
  use precision

  implicit none
! LOCAL variables that are not in modules

  save

  integer(kind=ik), parameter :: id_run = 1000829
  character(LEN=15)           :: suf_sec

end module io_security
!>
!> \verbatim
!> This subroutine determines (if a security file is necessary)
!> whether the code is allowed to run or not. It will also write
!> the appropiate security file
!>
!>  Author: A. Marek, MPA, March 2009
!> \endverbatim
!>
!> \param action  string that determines whether to check or write
!>                the security file
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
subroutine io_security_file(action)
  use precision

  use io_security
  use abort

  implicit none
  character*(*), intent(in) :: action

  integer(kind=ik)             :: id_rd

  if (action .eq. "read") then
     read(*,'(a)') suf_sec

     open(10,file = 'security.'//trim(suf_sec),form = 'formatted', &
          status = 'old')
     read(10,'(i7)') id_rd
     close(10, status='delete')

     if(id_rd .ne. id_run) then
        print *, "io_read_security(): no correct security file was found!"
        stop
     endif

     print *, 'security file: security.'// trim(suf_sec),' read; id: ',id_rd

  endif

  if (action .eq. "write") then
     open(10,file = 'security.'//trim(suf_sec),form = 'formatted')
     write(10,'(i7)') id_run
     close(10)
  endif

end subroutine io_security_file
