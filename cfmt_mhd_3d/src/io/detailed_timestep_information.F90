module timestep_information_mod

  implicit none

  contains

  subroutine print_detailed_timestep_info(ncycl, nbegin)
    use precision
    use cputim
#ifdef CHECK_MEMORY
   use meminfo
#endif

    implicit none

    integer(kind=ik), intent(in) :: ncycl, nbegin

    if (ncycl.le.nbegin+6) then 
#ifndef DEBUG_TIMINGS
       cputime_flag=.true.
#endif
#ifdef CHECK_MEMORY
       meminfo_flag=.true.
#ifdef SGE_QUEUE
       call query_sge_memory
#endif
       
#endif /* CHECK_MEMORY */
    endif
         
    if ( mod(ncycl,100) .eq. 0) then
#ifndef DEBUG_TIMINGS
       cputime_flag=.true.
#endif
#ifdef CHECK_MEMORY
       meminfo_flag=.true. ! print memory info every 100 cycles

#ifdef SGE_QUEUE
       call query_sge_memory
#endif

#endif /* CHECK_MEMORY */
    endif


  end subroutine print_detailed_timestep_info

end module timestep_information_mod
