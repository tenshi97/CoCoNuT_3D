module mod_newton_raphson
  use precision
  use function_pointers
  implicit none

  public :: newton_raphson

  contains

  function newton_raphson(f, x_guess, y, rtol, atol, maxit, additional_f_args, error) result(x)
    use function_pointers
    use numutils, only : are_close
    use abort

    implicit none

    interface
      function f(arg, additional_args, error) result(res)
        use precision
        use function_pointers, only : function_arguments_t
        implicit none
        real(kind=rk), intent(in)             :: arg(:)
        real(kind=rk)                         :: res(size(arg))
        class(function_arguments_t), intent(in), optional, target :: additional_args
        logical, optional, intent(out)        :: error
      end function
    end interface

    real(kind=rk),               intent(in)  :: x_guess(:)
    real(kind=rk),               intent(in)  :: y(size(x_guess))
    real(kind=rk)                            :: x(size(x_guess))
    real(kind=rk),               intent(in)  :: rtol, atol
    integer(kind=ik),            intent(in)  :: maxit

    class(function_arguments_t), intent(in), optional, target :: additional_f_args
    logical, optional,           intent(out) :: error

    real(kind=rk), dimension(size(x_guess)) :: m, yi
    integer :: i
    logical :: f_error

    x = x_guess

    do i = 1, maxit
      if (present(additional_f_args)) then
        yi = f(x, additional_f_args, f_error)
      else
        yi = f(x, error=f_error)
      endif
      if (f_error) then
        if (present(error)) then
          error = .true.
          return
        else
          raise_abort("Error evaluating function")
        endif
      endif
      if (present(additional_f_args)) then
        m = (f(x*1.001, additional_f_args, f_error) - yi) / (0.001 * x)
      else
        m = (f(x*1.001, error=f_error) - yi) / (0.001 * x)
      endif
      if (f_error) then
        if (present(error)) then
          error = .true.
          return
        else
          raise_abort("Error evaluating function")
        endif
      endif
      x = x - (yi - y) / m
      if (all(are_close(yi, y, rtol, atol))) then
        if (present(error)) then
          error = .false.
        endif
        return
      endif
    end do

    write(6,'(a,i4,a)') "Newton Raphson iteration did not converge after ", maxit, " iterations"
    write(*,'(a)') "[  i  x(i)             yi(i)            y(i)             abs(yi(i) - y(i))  ]"
    do i  = 1, size(x)
      if (.not. are_close(yi(i), y(i), rtol, atol)) then
        write(6,'(i4,1x,5(es16.9,1x))') i, x(i), yi(i), y(i), abs(yi(i) - y(i))
      endif
    enddo
    write(6,*)

    if (present(error)) then
      error = .true.
    else
      raise_abort("Newton Raphson iteration failed")
    endif
  end function

end module mod_newton_raphson

#ifdef PROGRAM_test_newton_raphson

module test_func
  use function_pointers, only : function_arguments_t

  implicit none

  contains

  function testfunction(arg, additional_args, error) result(res)
    use precision
    use abort
    implicit none
    real(kind=rk), intent(in) :: arg(:)
    real(kind=rk) :: res(size(arg))

    class(function_arguments_t), intent(in), optional, target :: additional_args
    logical, intent(out), optional :: error

    abort_if(present(additional_args))

    res = sin(arg)
    if (present(error)) then
      error = .false.
    endif
  end function
end module test_func

!make NO_CONFIG=1 EXTRA_CPPFLAGS=-DUNIT_TESTS
program test_newton_raphson
  use precision
  use mod_newton_raphson
  use test_func
  use abort
  use numutils, only : are_close
  real(kind=rk), dimension(3) :: x_guess, x
  real(kind=rk), dimension(3) :: y
  integer :: i
  logical :: error

  x_guess = [0.1_rk, 0.1_rk, 0.1_rk]
  y       = [0.2_rk, 0.3_rk, 0.4_rk]

  x = newton_raphson(testfunction, x_guess, y, 1e-16_rk, 1e-10_rk, 100, error=error)
  if (.not. error) then
    write(*, '(''['',3(a12,1x),'']'')') "x", "y", "f(x)"
    do i = 1, size(x)
      write(*,'('' '',3(f12.8,1x),'' '')') x(i), y(i), testfunction(x(i:i))
    end do
    abort_if(.not. all( are_close(testfunction(x), y, 1e-16_rk, 1e-10_rk)))
  else
    raise_abort("newton_raphson() failed!")
  endif

end program test_newton_raphson
#endif
