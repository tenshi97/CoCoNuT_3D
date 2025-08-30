module function_pointers

  implicit none
  save

  type function_arguments_t
#ifdef INTEL_COMPILER
    ! ifort does not support empty types yet
    integer :: x__dymmy_entry
#endif
  end type

end module
