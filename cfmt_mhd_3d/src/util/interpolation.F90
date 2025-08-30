module interpolation
  use precision
  implicit none
  public
  save

  !> The derived type for an 1D interpolation weight
  type weight_t
    integer(kind=ik) :: n !< The index in the original array xp for this value,
                          !! that is, the the largest n for which xp(n) < value
    real(kind=rk)    :: w !< The weight between xp(n) and xp(n+1), 0 <= w < 1.
                          !! value == xp(n) * (1 - w) + w * xp(n + 1)
  end type

  contains

  !> 1D linear interpolation
  !>
  !> \param     w       The integration weights
  !> \param     yp      The function values at the grid points
  function interpol_1d(w, yp) result(y)
    type(weight_t), intent(in) :: w(:)
    real(kind=rk), intent(in) :: yp(:)
    real(kind=rk)             :: y(size(w))
    integer :: i
    do i = 1, size(w)
      y(i) = yp(w(i)%n) + w(i)%w * (yp(w(i)%n + 1) - yp(w(i)%n))
    end do
  end function


  !> 1D linear interpolation on an irregular grid
  !>
  !> This is a convenience function which also calculates the integration
  !> weights
  !>
  !> \param     x       The values at which to interpolate
  !> \param     xp      The values of the grid, must be monotonically ascending!
  !> \param     yp      The function values at the grid points
  !>
  !> Throws an abort() should x not be in [xp(1), xp(size(xp))]
  function interpol_1d_irregular(x, xp, yp) result(y)
    use abort
    real(kind=rk), intent(in) :: x(:)
    real(kind=rk), intent(in) :: xp(:)
    real(kind=rk), intent(in) :: yp(size(xp))
    real(kind=rk) :: y(size(x))

    y = interpol_1d(interpol_1d_irregular_weights(x, xp), yp)
  end function

  !> 1D linear interpolation on a regular grid
  !>
  !> This is a convenience function which also calculates the integration
  !> weights
  !>
  !> \param     x       The values at which to interpolate
  !> \param     xp      The values of the grid, must be monotonically ascending!
  !> \param     yp      The function values at the grid points
  !>
  !> Throws an abort() should x not be in [xp(1), xp(size(xp))]
  function interpol_1d_regular(x, xp, yp) result(y)
    use abort
    real(kind=rk), intent(in) :: x(:)
    real(kind=rk), intent(in) :: xp(:)
    real(kind=rk), intent(in) :: yp(size(xp))
    real(kind=rk) :: y(size(x))

    y = interpol_1d(interpol_1d_regular_weights(x, xp), yp)
  end function

  !> Get interpolation weights for 1D linear interpolation on an irregular grid
  !>
  !> \param     x       The values at which to interpolate
  !> \param     xp      The values of the grid, must be monotonically ascending!
  !>
  !> \returns   w       The integration weights, a type(weight_t) array
  !>
  !> Throws an abort() should x not be in [xp(1), xp(size(xp))]
  function interpol_1d_irregular_weights(x, xp) result(w)
    use abort
    real(kind=rk), intent(in) :: x(:)
    real(kind=rk), intent(in) :: xp(:)
    type(weight_t) :: w(size(x))

    integer :: i, n
    logical :: not_found

    character(200) :: errormessage

    do i = 1, size(x)
      n = 1
      not_found = .true.
      do while (n < size(xp) .and. not_found)
        if (xp(n) <= x(i) .and. x(i) < xp(n + 1)) then
          w(i)%n = n
          w(i)%w = (x(i) - xp(n)) / (xp(n+1) - xp(n))
          not_found = .false.
        endif
        n = n + 1
      end do
      if (not_found) then
        write(errormessage,'(''Requested value '',es11.5,'' is outside the&
          & interpolation domain of ['',es11.5,'', '',es11.5,'']'')') x(i), xp(1), xp(size(xp))
        raise_abort(errormessage)
#ifdef UNIT_TESTS
        ! when doing a unit test, raise_abort() potentially returns!
        w(i)%n = 1
        w(i)%w = 0.0_rk
#endif
      endif
    end do
  end function

  !> Get interpolation weights for 1D linear interpolation on a regular grid
  !>
  !> \param     x       The values at which to interpolate
  !> \param     xp      The values of the regular grid
  !>                    For all i: xp(i + 1) - xp(i) == dx = const
  !>
  !> \returns   w       The integration weights, a type(weight_t) array
  !>
  !> Throws an abort() should x not be in [xp(1), xp(size(xp))]
  function interpol_1d_regular_weights(x, xp) result(w)
    use abort
    real(kind=rk), intent(in) :: x(:)
    real(kind=rk), intent(in) :: xp(:)
    type(weight_t) :: w(size(x))

    integer :: i, n
    logical :: not_found

    character(200) :: errormessage

    do i = 1, size(x)
      n = floor( 1.0_rk + real(size(xp) - 1,kind=rk) * ( x(i) - xp(1)) / (xp(size(xp))-x(1)) )
      if (n < 1 .or. n > size(xp)) then
        write(errormessage,'(''Requested value '',es11.5,'' is outside the&
          & interpolation domain of ['',es11.5,'', '',es11.5,'']'')') x(i), xp(1), xp(size(xp))
        raise_abort(errormessage)
#ifdef UNIT_TESTS
        ! when doing a unit test, raise_abort() potentially returns!
        w(i)%n = 1
        w(i)%w = 0.0_rk
#endif
      endif
      w(i)%n = n
      w(i)%w = (x(i) - xp(n)) / (xp(n+1) - xp(n))
    end do
  end function

end module interpolation

#ifdef PROGRAM_test_interpolation
!make NO_CONFIG=1 EXTRA_CPPFLAGS=-DUNIT_TESTS
program test_interpolation
  use precision
  use interpolation
  use numutils, only : are_close
  use abort

  implicit none
  real(kind=rk) :: xp(5), yp(5)
  real(kind=rk) :: x(11), y_actual(11), y_interpolated(11)
  integer :: i, j

  xp = [0.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk]
  yp = stepwise_linear_test_function(xp)
  x = [(i * maxval(xp)/(size(x) - 1),i=0,10)]
  x(size(x)) = x(size(x)) * (1.0_rk - epsilon(1.0_rk))
  
  write(*,'(a)') "Test the irregular interpolation routine"
  y_actual = stepwise_linear_test_function(x)
  y_interpolated = interpol_1d_irregular(x, xp, yp)
  write(*,'(''['',3(a12,x),'']'')') "x", "f(x)", "interpol."
  do i = 1, size(x)
    write(*,'('' '',3(f12.3,x),'' '')') x(i), y_actual(i), y_interpolated(i)
  end do
  write(*,*)

  abort_if(.not. all(are_close(y_actual, y_interpolated, 1e-15, 1e-15)))

  ! xp is in fact regular.
  write(*,'(a)') "Test the regular interpolation routine"
  y_interpolated = interpol_1d_regular(x, xp, yp)
  write(*,'(''['',3(a12,x),'']'')') "x", "f(x)", "interpol."
  do i = 1, size(x)
    write(*,'('' '',3(f12.3,x),'' '')') x(i), y_actual(i), y_interpolated(i)
  end do
  write(*,*)

  abort_if(.not. all(are_close(y_actual, y_interpolated, 1e-15, 1e-15)))

  write(*,'(a)') "Test out-of-bound handling"
  call expect_abort()
  y_interpolated = interpol_1d_irregular(x, xp(1:4), yp(1:4))
  abort_if(.not. got_abort())
  write(*,*) "Correctly raised an abort()!"

  contains

  function stepwise_linear_test_function(x) result(y)
    implicit none
    real(kind=rk), intent(in) :: x(:)
    real(kind=rk) :: y(size(x))
    y = floor(x)**2 + (floor(x+1)**2 - floor(x)**2) * (x - floor(x))
  end function

end program
#endif
