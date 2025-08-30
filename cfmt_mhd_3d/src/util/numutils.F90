module numutils

  use precision

  implicit none
  save

  interface are_close
    module procedure are_close_scalar
    module procedure are_close_array
  end interface

  contains

  function are_close_scalar(a, b, rtol, atol) result(ret)
    implicit none
    real(kind=rk), intent(in) :: a, b
    real(kind=rk), intent(in) :: rtol, atol
    logical :: ret

    ret = abs(a - b) <= atol + rtol * max(abs(a),abs(b))
  end function

  function are_close_array(a, b, rtol, atol) result(ret)
    implicit none
    real(kind=rk), intent(in) :: a(:), b(size(a))
    real(kind=rk), intent(in) :: rtol, atol
    logical :: ret(size(a))

    ret = abs(a - b) <= atol + rtol * max(abs(a),abs(b))
  end function

end module
