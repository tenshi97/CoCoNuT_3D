module state
  use precision

  implicit none

  type hydro_states
    real(kind=rk) :: dt

  end type hydro_states

  type(hydro_states) :: hydro
  type(hydro_states) :: transport
end module state
