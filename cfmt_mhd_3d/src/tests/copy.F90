!make CONFIG=make.inc.vertex
program copy
  use precision

  use machcons
#ifdef WRITE_BINARY_OUTPUT
  use dataformat, only : uninit_io
#endif
  use init_mpi_mod, only : init_mpi

  use initialize_array_mod, only : initialize_arrays
  use initialize_run_mod,   only : initialize_run

#ifdef WRITE_BINARY_OUTPUT
  use restart_model_mod,    only : restart_from_model

  use output_hydro, only : set_filenames
  use restart, only : restrt
#endif
  use parameters
  use configure

  use state
  use hydro_init_mod

  implicit none

  integer(kind=ik)            :: nbegin
  real(kind=rk)               :: dtt

  call set_machcons
  call init_mpi

  call initialize_run
  call initialize_arrays

  call read_parameter_files("model_init")

  call initialize_hydro(1)
  call initialize_arrays
  call initialize_hydro(2)

  call restart_from_model(config%irstrt, nbegin, dtt)

  call set_filenames(99,99)
  call restrt(0)

  call uninit_io

end program copy
