module environment
  implicit none
  public :: getenv
  private

  interface
    function getenv_c(name) result(val) bind(C, name="getenv")
      use iso_c_binding
      implicit none
      character(kind=c_char), intent(in) :: name
      type(c_ptr) :: val
    end function
  end interface

  interface
    function strlen(string) result(len) bind(C, name="strlen")
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: string
      integer(c_size_t) :: len
    end function
  end interface

  contains

  function cstring(cptr) result(string)
    use iso_c_binding
    implicit none
    type(c_ptr), intent(in) :: cptr
    character(kind=c_char), pointer :: string(:)
    if (.not. c_associated(cptr)) then
      nullify(string)
    else
      call c_f_pointer(cptr, string, (/ int(strlen(cptr)) /))
    endif
  end function

  function getenv(name) result(val)
    use iso_c_binding
    implicit none
    character(*), intent(in) :: name
    character(len=len(name) + 1), target :: name_terminated
    character(kind=c_char), pointer :: val(:)

    name_terminated = name // C_NULL_CHAR

    val => cstring(getenv_c(name_terminated(1:1)))
  end function

end module environment
