module benchmark_mod
  use precision
  use print_stdout_mod
  implicit none


#if defined(BENCHMARK) && !defined(DEBUG_TIMINGS)
! special variables for benchmark timing
      real(kind=rk)              :: time_initialize_run_start(2),         &
                                    time_initialize_run_stop(2),          &
                                    time_initialize_arrays_start(2),      &
                                    time_initialize_arrays_stop(2),       &
                                    time_readparam_start(2),              &
                                    time_readparam_stop(2),               &
                                    time_restart_from_new_model_start(2), &
                                    time_restart_from_new_model_stop(2),  &
                                    time_set_counting_parameters_start(2),&
                                    time_set_counting_parameters_stop(2), &
                                    time_construct_new_model_start(2),    &
                                    time_construct_new_model_stop(2),     &
                                    time_write_final_output(2),           &
                                    time_initialize_hydro_start(2),       &
                                    time_initialize_hydro_stop(2)

#endif /* BENCHMARK */



  contains

    subroutine print_benchmark_summary(irstrt, nbegin, ncycl, starttime, stoptime, firstcycletime, secondcycletime)

      use precision
      use configure
      implicit none

      integer(kind=ik), intent(in) :: irstrt, nbegin, ncycl
      real(kind=rk), intent(in)    :: starttime(:), stoptime(:), firstcycletime(:), secondcycletime(:)
      call printit_taskX(0," ")
      call printit_taskX(0," Code runtime:       ",stoptime(2)-starttime(2))
      call printit_taskX(0," Cycles:             ",(ncycl-nbegin-1))
      call printit_taskX(0," Now at cycle:       ",ncycl)
      call printit_taskX(0," Started with cycle: ",nbegin)

      ! special benchmakr timing output in Wallclock time
      call printit_taskX(0," ")
      call printit_taskX(0,"=================================================")
      call printit_taskX(0," Benchmark times")
      call printit_taskX(0," ")
      call printit_taskX(0," Time until calculation started (incl. some ")
      call printit_taskX(0," setup:                   ",starttime(2))
      call printit_taskX(0,"    - time initialize run ",time_initialize_run_stop(2) - time_initialize_run_start(2))
      call printit_taskX(0,"    - initialize arrays   ",time_initialize_arrays_stop(2) - time_initialize_arrays_start(2))
      call printit_taskX(0,"    - read param          ",time_readparam_stop(2)  - time_readparam_start(2))
      call printit_taskX(0,"       incl. initialize_hydro(1), initialize_arrays")

      call printit_taskX(0,"    - initialize_hydro(2) ",time_initialize_hydro_stop(2) - time_initialize_hydro_start(2))
      if (irstrt .eq. 0) then
      call printit_taskX(0,"    - new model           ",time_construct_new_model_stop(2) - time_construct_new_model_start(2))
      else
      call printit_taskX(0,"    - restart model       ",time_restart_from_new_model_stop(2) - time_restart_from_new_model_start(2))
      endif
      call printit_taskX(0,"    - set counting params ",time_set_counting_parameters_stop(2) - time_set_counting_parameters_start(2))
      call printit_taskX(0," Time for first cycle:  ",firstcycletime(2) - &
                                                       starttime(2) )

      call printit_taskX(0," Time for second cycle: ",secondcycletime(2))
      if (.not.(config%use_deactivate_output)) then
         call printit_taskX(0," Time for final ouput:  ",time_write_final_output(2))
      endif
      call printit_taskX(0,"=================================================")      


    end subroutine print_benchmark_summary



end  module benchmark_mod
