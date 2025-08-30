module stringutils
  implicit none

  contains

  function endswith(string, suffix) result(test)
    implicit none
    character(len=*), intent(in) :: string, suffix
    logical :: test
    integer :: l1, l2

    l1 = len_trim(string)
    l2 = len_trim(suffix)

    if (l2 > l1) then
      test = .false.
    else
      test = string(l1-l2+1:l1) == suffix
    endif
  end function

end module

#ifdef PROGRAM_test_stringutils
!make NO_CONFIG=1
program test_stringutils
  use stringutils
  use abort

  implicit none
  character(len=20) :: foo

  foo = "lorem ipsum"

  if (.not. endswith(foo, "psum")) then
    raise_abort("error")
  end if
end program
#endif
