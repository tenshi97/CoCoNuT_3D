

module initialize_run_mod

  use configure
  implicit none

  contains

  subroutine initialize_run
  use mpi_vertex
  use setup_c
  use hydro_memory
#ifndef DEBUG_TIMINGS
  use cputim
#endif

#ifdef CHECK_THREAD_AFFINITY
  use thread_affinity
#endif

#ifdef CFC_TRANSPORT
  use gr_alloc, only : allocate_coconut_arrays
#endif

#ifdef CHECK_MEMORY
  use meminfo
#endif
  implicit none

#ifdef CHECK_MEMORY
   call init_memory_control
#endif 

   if (config%use_stdout_redirect) then
      ! re direct output and stderr of myproc != 0 to file stdout_myrank.txt and
      ! stderr_myrank.txt 
      if (use_mpi) call redirect_stdout(myproc)
   endif

#ifdef LINUX
!      call unlimit_stack
#endif

#ifndef DEBUG_TIMINGS

! next we want to be able to compute timings
      call allocate_cputim
#endif

      ! now we can allocate the arrays from module variables both in
      ! the case of MPI/OpenMP or only OpenMP mode
#ifdef CHECK_MEMORY
      meminfo_flag=.true.
#endif
      call allocate_hydro_memory


#ifdef CFC_TRANSPORT
      call allocate_coconut_arrays
#endif

#ifdef CHECK_MEMORY
      meminfo_flag=.false.
#endif

#ifdef CHECK_THREAD_AFFINITY
      call init_thread_affinity
#endif

    end subroutine initialize_run

end module initialize_run_mod
