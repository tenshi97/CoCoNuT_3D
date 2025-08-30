#undef DEBUG

#if 0
module secant_method
  use precision

  implicit none
  save

  contains

  !> \param    f                        function to invert
  !> \param    x0_in                    value with f(x0_in, additional_f_args) < y
  !> \param    x1_in                    value with f(x1_in, additional_f_args) > y
  !> \param    y                        desired function value
  !> \param    additional_f_args        additional function arguments
  !> \param    eps                      allowed relative error for solution
  !> \param    i_max                    maximum number of iterations
  !> \param    eval_error               If supplied, do not abort on error but return .false. in this flag
  function secant(f,x0_in,x1_in,y,additional_f_args,eps,i_max,eval_error) result(x)
    use abort
#if 0
    use function_pointers, only : function_arguments_t
    implicit none

    ! arguments
    interface
      function f(arg, additional_args, error) result(res)
        use precision
        use function_pointers, only : function_arguments_t
        implicit none
        real(kind=rk), intent(in)                       :: arg
        class(function_arguments_t), intent(in), target :: additional_args
        real(kind=rk)                                   :: res
        logical, optional, intent(out)                  :: error
      end function
    end interface

    real(kind=rk), intent(in)                       :: x0_in, x1_in, y
    class(function_arguments_t), intent(in), target :: additional_f_args
    real(kind=rk), intent(in)                       :: eps
    integer(kind=ik), intent(in)                    :: i_max
    logical, optional, intent(out)                  :: eval_error

    ! result
    real(kind=rk) :: x

    ! local variables
    integer(kind=ik) :: i
    real(kind=rk)    :: x0, x1, x2, delta_x, delta_f, y0, y1, y2
    logical          :: f_error

    x0 = x0_in
    x1 = x1_in

    if (present(eval_error)) then
      eval_error = .false.
    endif

#ifdef DEBUG
    write(*,*)
    write(*,'(a,es12.5)') "Start of secant method for y = ", y
    debug(x0_in)
    debug(x1_in)
#endif

    y0 = f(x0, additional_f_args, f_error)
    if (f_error) then
      if (present(eval_error)) then
        eval_error = f_error
        return
      else
        raise_abort("Error evaluating function")
      endif
    endif

    y1 = f(x1, additional_f_args, f_error)
    if (f_error) then
      if (present(eval_error)) then
        eval_error = f_error
        return
      else
        raise_abort("Error evaluating function")
      endif
    endif

    i=0

    if (y0 > y .or. y1 < y) then
       write(*,*) "Unsuitable initial values: y0 = ", y0, ", y1 = ", y1
       raise_abort("Unsuitable initial values!")
    endif

    do while (min(abs(y0 - y)/y, abs(y1 - y)/y) > eps)

       delta_x = x1 - x0
       delta_f = y1 - y0

       ! check for division by zero in not converged zones
       if (abs(delta_f/delta_x) <= 1e-300) then
         write(*,*) "Cannot find solution, singular point at iteration", i
         write(*,*) "[x0, x1 - x0, y0, y1 - y0]"
         write(*,'(i4,1x,4(1x,es12.5))') x0, x1 - x0, y0, y1 - y0
         write(*,*)
         raise_abort("Cannot find solution, df_dx -> 0")
       end if

       ! get point at the intersection of the
       ! lines (x0, y0)--(x1,y1) and (x0, y)--(x1, y)
       x2 = x1 - (y1 - y) * delta_x/delta_f
       y2 = f(x2, additional_f_args, f_error)
       if (f_error) then
         if (present(eval_error)) then
           eval_error = f_error
           return
         else
           raise_abort("Error evaluating function")
         endif
       endif

       ! determine new enclosing values
       if (y2 > y) then
         ! replace left value
         x1 = x2
       else
         ! replace right value
         x0 = x2
       endif

       y0 = f(x0, additional_f_args, f_error)
       if (f_error) then
         if (present(eval_error)) then
           eval_error = f_error
           return
         else
           raise_abort("Error evaluating function")
         endif
       endif

       y1 = f(x1, additional_f_args, f_error)
       if (f_error) then
         if (present(eval_error)) then
           eval_error = f_error
           return
         else
           raise_abort("Error evaluating function")
         endif
       endif

#ifdef DEBUG
       write(*,*)
       write(*,*) "Values after iteration ", i
       write(*,'(a,es12.5)') "# y = ", y
       write(*,'(2(a,es12.5))') "# x0 = ", x0, ", x1 = ", x1
       write(*,'(2(a,es12.5))') "# y0 - y = ", y0 - y, ", y1 - y = ", y1 - y
       write(*,'(2(a,es12.5))') "# |y0 - y| / y = ", abs(y0 - y)/y, ", |y1 - y|/y = ", abs(y1 - y)/y
       debug(min(abs(y0 - y)/y, abs(y1 - y)/y))
       write(*,*)
#endif

       i = i + 1

       if (i > i_max) then
          write(*,'(a,i0,a)') "## ERROR: Maximum number (", i_max, ") of iterations reached!"
          write(*,'(2(a,es12.5))') "## x0 = ", x0, ", x1 = ", x1
          write(*,'(2(a,es12.5))') "## y0 = ", y0, ", y1 = ", y1
          raise_abort("Maximum number of iterations reached!")
       endif

    enddo

    ! choose better solution
    if (abs(y0 - y) < abs(y1 - y)) then
      x = x0
    else
      x = x1
    endif

#endif

  end function secant

end module secant_method
#endif
