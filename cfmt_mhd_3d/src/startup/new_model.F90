module new_model_mod

  implicit none

  contains

  subroutine construct_new_model(code_start_time, nbegin)
  use precision
  use abort
  use error
  use mo_mpi
  use cputim
  use setup_c
  use gfloat_hy, only : time

  use intgrs_hy, only : nstep,  nresum
!  use param_rt,  only : p_ntr, irst_ra
  use charac, only : file_number, filpos
#ifdef WRITE_BINARY_OUTPUT
  use output_hydro, only : set_filenames, open_files, close_files, &
                           write_output_files

  use restart, only : restrt
#endif

#ifndef NOMOTO_MODEL
  use read_progenitor_model, only : read_inimod, init, progenitor_t, deallocate_progenitor
#else
  use read_progenitor_model, only : read_onemg, progenitor_t
#endif
  use grid_mod, only : grid

#if !defined(PROGRAM_remap) && !defined(NOTRA)
  use timeinfo, only : print_timestep_evolution
#endif
  use mpi_domains

  use initial_model


#ifdef CFC_TRANSPORT
  use coconut
  use gr_alloc
#endif

  use print_stdout_mod

  use hydro_areas_mod

  use configure
  use state
  implicit none

  integer(kind=ik) :: ier
  integer(kind=ik), intent(out) :: nbegin
  real(kind=rk), intent(in)     :: code_start_time(2)

  real(kind=rk), dimension(2)   :: transport_init_start_time, &
                                   transport_init_stop_time
  type(progenitor_t) :: progenitor


  if (myproc .eq. 0) then
     if (create_directories() .ne. 1) then
        raise_abort("vertex(): Unable to create directories")
     endif
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ier)
  time   = 0._rk

  hydro%dt     = config%dtini
  nstep  = 0
  nbegin = 0
  areas%nhystp = 0

  if(file_number .ne. 0) then
     raise_abort("vertex(): check input file")
  endif

  filpos = 'rewind'

  call printit_taskX(0," ")
  call printit_taskX(0,"Constructing new initial model")
  call printit_taskX(0," ")

#ifndef NOMOTO_MODEL
  progenitor = read_inimod(config%progenitor_name)
  call grid()
  call init(progenitor)
  call deallocate_progenitor(progenitor)
#else
  call grid
  call read_onemg
#endif /* NOMOTO_MODEL */

#ifdef CFC_TRANSPORT
  call init_coconut(config%irstrt)
#endif

!-PNS         tmatot(0) = pmass

  if (nstep .ne. 0)   then
     nbegin = nstep
     config%nend   = nstep + config%nend
  endif

#ifdef WRITE_BINARY_OUTPUT
  if (.not.(config%use_deactivate_output)) then
     call set_filenames(config%irstrt,config%irst_ra)
     call restrt(0)
  endif
#endif

  if (config%use_print_mpiDomain) then
     ! print MPI domain decomposition and thread affinity
     call print_mpi_setup
  endif

  call initial_model_energy_summary
#ifdef WRITE_BINARY_OUTPUT
  if (.not.(config%use_deactivate_output)) then
     call open_files
     call write_output_files
  endif
#endif

#ifndef NOTRA
  if (config%p_ntr .ne. 0) call print_timestep_evolution
#endif

  end subroutine construct_new_model


end module new_model_mod
