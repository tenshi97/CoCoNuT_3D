module debug_c
  use precision

  implicit none
  public

#ifdef BACKTRACE_ON_ABORT
  interface
    subroutine print_backtrace() bind(C, name="print_backtrace")
      use iso_c_binding
      implicit none
    end subroutine
  end interface
#endif

end module debug_c
