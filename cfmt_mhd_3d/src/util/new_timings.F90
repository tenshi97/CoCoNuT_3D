module new_timings
  use precision
  use setup_c, only : ms_since_epoch
  use, intrinsic :: iso_c_binding, only : C_INT64_T
  implicit none
  save

  private
  public start_timer, end_timer

  integer, parameter :: name_length = 32

  type timing_entry_t
    character(len=name_length) :: name
    integer(kind=C_INT64_T) :: ms
    logical :: open
    type(timing_entry_t), pointer :: firstChild
    type(timing_entry_t), pointer :: lastChild
    type(timing_entry_t), pointer :: parent
    type(timing_entry_t), pointer :: nextSibling
  end type

  type(timing_entry_t), pointer :: current_timer => NULL()

  character(len=name_length), parameter :: own = "(own)"

  contains

  recursive subroutine deallocate_timer(timing_entry)
    type(timing_entry_t), pointer :: timing_entry
    type(timing_entry_t), pointer :: nextSibling
    if (associated(timing_entry%firstChild)) then
      call deallocate_timer(timing_entry%firstChild)
    endif
    nextSibling => timing_entry%nextSibling
    deallocate(timing_entry)
    nullify(timing_entry)
    if (associated(nextSibling)) then
      call deallocate_timer(nextSibling)
    endif
  end subroutine

  subroutine start_timer(name)
    character(len=*) :: name
    type(timing_entry_t), pointer :: timing_entry

    ! create new timing_entry_t node
    allocate(timing_entry)
    timing_entry%name = name
    timing_entry%ms = ms_since_epoch()
    timing_entry%open = .true.

    nullify(timing_entry%firstChild)
    nullify(timing_entry%lastChild)
    nullify(timing_entry%nextSibling)

    if (.not. associated(current_timer)) then
      nullify(timing_entry%parent)
      current_timer => timing_entry
    else
      timing_entry%parent => current_timer

      ! sort into tree
      if (current_timer%open) then
        ! one level below current_timer
        if (.not. associated(current_timer%lastChild)) then
          ! no childs yet
          current_timer%firstChild => timing_entry
          current_timer%lastChild => timing_entry
          current_timer => timing_entry
        else
          ! other childs
          current_timer%lastChild%nextSibling => timing_entry
          current_timer%lastChild => timing_entry
          current_timer => timing_entry
        endif
      else
        ! same level as current_timer
        current_timer%nextSibling => timing_entry
        current_timer%parent%lastChild => timing_entry
        current_timer => timing_entry
      endif
    endif
  end subroutine

  subroutine end_timer(name)
    use abort
    character(len=*), intent(in) :: name
    integer :: i
    if (current_timer%name(1:len(name)) .ne. name) then
      raise_abort("Expected end_timer(""" // trim(current_timer%name)  // """), but got end_timer(""" // trim(name) // """)")
    endif
    do i = len(name) + 2, name_length
      current_timer%name(i:i) = "."
    end do
    current_timer%ms = ms_since_epoch() - current_timer%ms
    current_timer%open = .false.

    ! Climb up to parent
    if (associated(current_timer%parent)) then
      current_timer => current_timer%parent
    else
      write(*,'(a)') "  . Timings"
      call print_timings(current_timer, 0)
      write(*,*)
      call deallocate_timer(current_timer)
    endif
  end subroutine

  recursive subroutine print_timings(entry, indent_level)
    type(timing_entry_t), intent(in), pointer :: entry
    integer, intent(in) :: indent_level

    character(len=64) :: format_spec

    write(format_spec,'("(",i0,"x,""|_ "",a",i0,",2x,f10.6)")') indent_level * 2 + 2, len(entry%name)
    write(*,format_spec) entry%name, real(entry%ms) / 1e6

    if (associated(entry%firstChild)) then
      write(format_spec,'("(",i0,"x,""|_ "",a",i0,",2x,f10.6)")') (indent_level + 1) * 2 + 2, len(entry%name)
      write(*,format_spec) own, real(entry%ms - ms_of_children(entry)) / 1e6
    endif

    if (associated(entry%firstChild)) then
      call print_timings(entry%firstChild, indent_level + 1)
    endif
    if (associated(entry%nextSibling)) then
      call print_timings(entry%nextSibling, indent_level)
    endif
  end subroutine

  function ms_of_children(entry) result(sum_time)
    type(timing_entry_t), intent(in), pointer :: entry
    type(timing_entry_t), pointer :: cur_entry
    integer(kind=C_INT64_T) :: sum_time
    sum_time = 0
    cur_entry => entry%firstChild
    do while (associated(cur_entry))
      sum_time = sum_time + cur_entry%ms
      cur_entry => cur_entry%nextSibling
    enddo
  end function

end module

#ifdef PROGRAM_test_timings
!make NO_CONFIG=1 DEBUG=HIGH
program test_timings
  use new_timings
  implicit none
  integer :: i

  do i = 1,3
    call start_timer("cycle")
    call a()
    call a()
    call end_timer("cycle")
  end do

  contains

  subroutine a()
    call start_timer("a")
    call b
    call c
    call end_timer("a")
  end subroutine

  subroutine b()
    call start_timer("b")
    call c
    call c
    call end_timer("b")
  end subroutine

  subroutine c()
    call start_timer("c")
    call end_timer("c")
  end subroutine
end program
#endif
