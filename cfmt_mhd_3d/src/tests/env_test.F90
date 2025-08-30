!make NO_CONFIG=1 DEBUG=HIGH
program env_test
  use iso_c_binding
  use environment
  use abort

  implicit none
  character(c_char), pointer :: val(:)

  val => getenv("PATH")

  abort_if(.not. associated(val), "Variable PATH undefined!")

  write(*,*) "PATH = '", val, "'"
end program
