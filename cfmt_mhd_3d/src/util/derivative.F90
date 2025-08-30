module mod_derivative
  use precision
  implicit none

  contains

  function diff(f, x) result(df)
    implicit none
    real(kind=rk), intent(in) :: f(:)
    real(kind=rk), intent(in) :: x(size(f))
    real(kind=rk) :: df(size(f))

    real(kind=rk), dimension(size(f)) :: x12, x01, x02

    x12 = x - cshift(x, +1)             ! x1 - x2
    x01 = cshift(x, -1) - x             ! x0 - x1
    x02 = cshift(x, -1) - cshift(x, +1) ! x0 - x2

    df = cshift(f, -1) * (x12 / (x01*x02)) + f * (1.0/x12 - 1.0/x01) - cshift(f, +1) * (x01 / (x02 * x12))

    df(1) =  f(1)  *  (x01(2)+x02(2))/(x01(2)*x02(2)) - &
             f(2)  *           x02(2)/(x01(2)*x12(2)) + &
             f(3)  *           x01(2)/(x02(2)*x12(2))

    df(size(f)) = -f(size(f)-2) *           x12(size(f)-1)/(x01(size(f)-1)*x02(size(f)-1)) + &
                   f(size(f)-1) *           x02(size(f)-1)/(x01(size(f)-1)*x12(size(f)-1)) - &
                   f(size(f)) * (x02(size(f)-1)+x12(size(f)-1))/(x02(size(f)-1)*x12(size(f)-1))

  end function

end module

#ifdef PROGRAM_test_diff
!make NO_CONFIG=1
program test_diff
  use mod_derivative
  use precision
  use abort
  integer, parameter :: N = 20
  real(kind=rk), parameter :: pi = acos(-1.0_rk)
  real(kind=rk), dimension(0:N) :: df, f, x
  real(kind=rk) :: h, eps, w
  integer :: i, o

  do o = 0, 4
    w = 1.0_rk * pi * 10.0**(-o)
    x = [ (real(i, kind=rk) / (N - 1), i = -N/2, N/2) ] * w
    h = x(1) - x(0)
    f = sin(x)
    df = diff(f, x)

    write(*,'(6(x,a12))') "x", "f", "dfdx exact", "numerical", "eps/h**2", "eps"
    do i = 0, N
      eps = abs((df(i) - cos(x(i))) / cos(x(i)))
      write(*,'(5(x,f12.8),x,es12.5)') x(i), f(i), cos(x(i)), df(i), eps/h**2, eps
      if (i > 0 .and. i < N) then
        abort_if(eps/h**2 > 0.2)
      else
        abort_if(eps/h**2 > 0.4)
      endif
    end do
    write(*,*)
  end do
end program
#endif
