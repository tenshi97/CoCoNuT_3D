module set_cparam_mod
  use configure

  implicit none

  contains

  subroutine set_cparam(next_rho_c, starttime)
    use precision
    use intgrs_hy, only : nout1, nrst, nresum
    use gfloat_hy, only : trst, tout1
#if !defined(PROGRAM_remap)
    use coredensity
#endif
    use gfloat_hy, only : time
    use print_stdout_mod
    implicit none

    real(kind=rk), intent(out) :: next_rho_c, starttime
    real(kind=rk), parameter :: next_rho_c_factor = 1.25892541179416721_rk

    if (config%collapse_time_resolve .eq. 1) then
       next_rho_c = core_density() * next_rho_c_factor
    else
       next_rho_c = 0._rk
    endif

      nout1   = 0
      tout1   = 0._rk
      nrst    = 0
      trst    = 0._rk
      nresum  = 0

      starttime = time



#ifdef DEBUG_TIMINGS
      call printit_taskX(0," ")
      call printit_taskX(0,"===================================================")
      call printit_taskX(0," Attention you switched off all timing information")
      call printit_taskX(0," code does not know anything about the time it is ")
      call printit_taskX(0," already running")
      call printit_taskX(0,"===================================================")
      call printit_taskX(0," ")
 
#endif
  end subroutine set_cparam
end module set_cparam_mod
