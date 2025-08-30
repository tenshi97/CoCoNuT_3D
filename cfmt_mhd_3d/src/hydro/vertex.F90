!> \mainpage
!> The VERTEX-Supernova Code for simulations of core collapse supernovae
!>
!> \par provides the main program vertex
!>
!> \author M. Rampp, R. Buras, A. Marek, B. Mueller, L. Huedepohl, and F. Hanke
!>
!> \detail
!>    The multi grid/ time stepping treatment, neutrino transport,
!>    EOSLS-interface was developed by
!>
!>    - W. Keil,         MPI f. Astrophysik,     Garching
!>    - M. Rampp,        MPI f. Astrophysik,     Garching
!>
!>   the 2D-neutrino transport was developed by
!>
!>   - R. Buras         MPI f. Astrophysik,     Garching
!>   - M. Rampp         MPI f. Astrophysik,     Garching
!>
!>   all changes / improvements / additional featues since 2004 by
!>
!>   - R. Buras         MPI f. Astrophysik,     Garching
!>   - A. Marek         MPI f. Astrophysik,     Garching
!>   - B. Mueller       MPI f. Astrophysik,     Garching
!>   - L. Huedepohl     MPI f. Astrophysik,     Garching
!>   - F. Hanke         MPI f. Astrophysik,     Garching
!>
!>   Code optimisations for NEC SX systems and SGI Altix 3700
!>   - K. Benkert       HLRS                    Stuttgart
!>   - R. Fischer       NEC, Germany
!>   - R, Vogelsang     SGI, Germany!>
!>   References:
!>              Rampp M., and Janka H.-Th.,
!>              "Radiation hydrodynamics with neutrinos. Variable
!>               Eddington factor method for core-collapse supernova
!>               simulations", A&A, 396, 2002
!>
!>              Buras R., Rampp M., Janka H.-Th., and Kifonidis K.,
!>              "Two-dimensional hydrodynamic core-collapse supernova
!>               simulations with spectral neutrino transport:
!>               1 Numerical method and results of a 15msol star",
!>               A&A, 447, 2006
!> \verbatim
!>   SVN - Information
!>   $Revision:$
!>   $Date:$
!> \endverbatim
!>


!> Main routine of the VERTEX-Code
!> \author originally M. Rampp
    program vertex


      use precision
      use abort
      use error

#ifndef NOTRA

#endif
!      use phycon
      use intgrs_hy, only : nstep,  nout1, nrst, nresum, &
                             igrav
      ! use arecon_hy, only : nhystp
      use charac, only    : file_number, filpos! , calculation
   !   use param_rt, only  : produc, wcrem

#ifndef DEBUG_TIMINGS
      use cputim
#endif
      use machcons
!      use hlpare_hy

      use hydro_memory, only : deallocate_hydro_memory
      use mo_mpi

!      use io_security
      use eos_sn2, only : lsrolo
      use run_consistency, only : check_run_consistency
#ifdef WRITE_BINARY_OUTPUT
      use dataformat, only : uninit_io
#endif
#ifdef CHECK_MEMORY
      use meminfo
#endif

      use resolution_checker

      use mpi_domains

      use coredensity
      use stopentropy

      use timings, only : timings_out

!      use parameter_files


      use setup_c

      use init_mpi_mod, only : init_mpi

      use rady_mod, only : rady
      use initialize_array_mod, only : initialize_arrays
!      use initialize_hydro_mod, only : initialize_hydro
      use initialize_run_mod,   only : initialize_run

      use new_model_mod,        only : construct_new_model
      use restart_model_mod,    only : restart_from_model


      use gfloat_hy, only : tout1, trst, time! , trstrt
#ifdef WRITE_BINARY_OUTPUT
      use output_hydro, only : set_filenames, open_files, close_files, &
                               write_output_files
      use restart, only : restrt
#endif

#if !defined(PROGRAM_remap) && !defined(NOTRA)
      use timeinfo, only : print_timestep_evolution
#endif

      use set_cparam_mod             , only : set_cparam
      use timestep_information_mod,    only : print_detailed_timestep_info
      use check_security_file_mod,     only : check_for_security_file

      use spherical_neutron_star_mod,  only : spherical_neutron_star

      use print_stdout_mod

#ifdef BENCHMARK
      use benchmark_mod
#endif


      use parameters
      use configure

      use state
      use hydro_init_mod
      use machcons

#ifdef CFC_TRANSPORT
      use size_cfc, only : init_size_cfc
#endif

      implicit none

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! LOCAL variables that are not in modules

      integer(kind=ik)            :: ier, nbegin, ncycl, ithreads


      ! all times will be normalised to code_start_time !
      real(kind=rk), dimension(2) :: code_start_time,                  &
                                     calculation_stop_time,            &
                                     calculation_start_time,           &
                                     cycle_start_time, cycle_end_time, &
                                     transport_init_start_time,        &
                                     transport_init_stop_time,         &
                                     write_output_start_time,          &
                                     write_output_stop_time,           &
                                     code_run_time, time_remain
#ifdef BENCHMARK
      real(kind=rk), dimension(2) :: time_first_cycle(2), time_second_cycle(2)
#endif

      real(kind=rk)               :: dtt, ti_cyc, cpur
      real(kind=rk)               :: dt_cyc

#ifdef DEBUG_TIMINGS
      real(kind=rk)               :: cpumax=1e10
#endif

      real(kind=rk) :: next_rho_c


      integer :: i,j,k,nsb
      integer :: ierr,nsb_rcv

      real(kind=rk) :: starttime


    real(kind=rk), parameter :: next_rho_c_factor = 1.25892541179416721_rk

    real(kind=rk) :: rady_self(2), rady_children(2)

#if defined(OPENMP_HYD) && (defined(OPEN_MP_1D) || defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
! OMP variables
!$    integer OMP_GET_NUM_THREADS,numthreads
#endif


#ifndef DEBUG_TIMINGS
      cputime_flag = .true.
#endif

      call set_machcons

      call init_mpi  ! initialize MPI and the VERTEX mpi values qy_s,
                     ! qy_e, qz_s, qz_e

#ifndef DEBUG_TIMINGS
      ! all time meassurements will be normalized to this
      ! time point
      call second_v(code_start_time)
#endif
      if (use_mpi) code_start_time(2)=MPI_Wtime()

#if defined(BENCHMARK) && !defined(DEBUG_TIMINGS)
      call second_v(time_initialize_run_start)
      if (use_mpi) time_initialize_run_start(2)=MPI_Wtime()
      time_initialize_run_start(:) = time_initialize_run_start(:) - code_start_time(:)
#endif
! setup up the enviroment and control functions needed
! to run the code

#if defined(CFC_TRANSPORT) || defined(CFC_TRANSPORT2)
      call init_size_cfc
#endif

      call initialize_run

#if defined(BENCHMARK) && !defined(DEBUG_TIMINGS)
      call second_v(time_initialize_run_stop)
      if (use_mpi) time_initialize_run_stop(2)=MPI_Wtime()
      time_initialize_run_stop(:) = time_initialize_run_stop(:) - code_start_time(:)
#endif


#if defined(BENCHMARK) && !defined(DEBUG_TIMINGS)
      call second_v(time_initialize_arrays_start)
      if (use_mpi) time_initialize_arrays_start(2)=MPI_Wtime()
      time_initialize_arrays_start(:) = time_initialize_arrays_start(:) - code_start_time(:)
#endif
!! set all the hydro and transport arrays to zero
!! perform an OpenMP firsttouch
!      call initialize_arrays


#if defined(BENCHMARK) && !defined(DEBUG_TIMINGS)
      call second_v(time_initialize_arrays_stop)
      if (use_mpi) time_initialize_arrays_stop(2)=MPI_Wtime()
      time_initialize_arrays_stop(:) = time_initialize_arrays_stop(:) - code_start_time(:)
#endif


!----------------------------------------------------------------------
!     read parameter file
!----------------------------------------------------------------------
#if defined(BENCHMARK) && !defined(DEBUG_TIMINGS)
     call second_v(time_readparam_start)
     if (use_mpi) time_readparam_start(2)=MPI_Wtime()
     time_readparam_start(:) = time_readparam_start(:) - code_start_time(:)
#endif

      call read_parameter_files("model_init")

      ! set a few missing hydro parameters
      call initialize_hydro(1)


! set all the hydro and transport arrays to zero
! perform an OpenMP firsttouch
      call initialize_arrays

      if (config%produc .eq. 1) then
        call io_security_file("read")
      endif

#if defined(BENCHMARK) && !defined(DEBUG_TIMINGS)
      call second_v(time_readparam_stop)
      if (use_mpi) time_readparam_stop(2)=MPI_Wtime()
      time_readparam_stop(:) = time_readparam_stop(:) - code_start_time(:)
#endif


#if defined(BENCHMARK) && !defined(DEBUG_TIMINGS)
     call second_v(time_initialize_hydro_start)
     if (use_mpi) time_initialize_hydro_start(2)=MPI_Wtime()
     time_initialize_hydro_start(:) = time_initialize_hydro_start(:) - code_start_time(:)
#endif
     ! initialize subroutines that are needed for doing
     ! a hydro simulation
     call initialize_hydro(2)

#if defined(BENCHMARK) && !defined(DEBUG_TIMINGS)
     call second_v(time_initialize_hydro_stop)
     if (use_mpi) time_initialize_hydro_stop(2)=MPI_Wtime()
     time_initialize_hydro_stop(:) = time_initialize_hydro_stop(:) - code_start_time(:)
#endif

!----------------------------------------------------------------------
!     construct new model:
!----------------------------------------------------------------------

      if(config%irstrt .eq. 0) then

#if defined(BENCHMARK)  && !defined(DEBUG_TIMINGS)
         call second_v(time_construct_new_model_start)
         if (use_mpi) time_construct_new_model_start(2)=MPI_Wtime()
         time_construct_new_model_start(:) = time_construct_new_model_start(:) - code_start_time(:)
#endif


         call construct_new_model(code_start_time, nbegin)

#if defined(BENCHMARK)  && !defined(DEBUG_TIMINGS)
         call second_v(time_construct_new_model_stop)
         if (use_mpi) time_construct_new_model_stop(2)=MPI_Wtime()
         time_construct_new_model_stop(:) = time_construct_new_model_stop(:) - code_start_time(:)
#endif


!----------------------------------------------------------------------
!     read restart file:
!----------------------------------------------------------------------
      else

#if defined(BENCHMARK)  && !defined(DEBUG_TIMINGS)
         call second_v(time_restart_from_new_model_start)
         if (use_mpi) time_restart_from_new_model_start(2)=MPI_Wtime()
         time_restart_from_new_model_start(:) = time_restart_from_new_model_start(:) - code_start_time(:)
#endif

         call restart_from_model(config%irstrt, nbegin, dtt)

#if defined(BENCHMARK)  && !defined(DEBUG_TIMINGS)
         call second_v(time_restart_from_new_model_stop)
         if (use_mpi) time_restart_from_new_model_stop(2)=MPI_Wtime()
         time_restart_from_new_model_stop(:) = time_restart_from_new_model_stop(:) - code_start_time(:)
#endif
      endif ! config%irstrt

!-----------------------------------------------------------------------
!     set up a few (counting) parameters:
!-----------------------------------------------------------------------

#if defined(BENCHMARK)  && !defined(DEBUG_TIMINGS)
         call second_v(time_set_counting_parameters_start)
         if (use_mpi) time_set_counting_parameters_start(2)=MPI_Wtime()
         time_set_counting_parameters_start(:) = time_set_counting_parameters_start(:) - code_start_time(:)
#endif

      call set_cparam(next_rho_c, starttime)

#if defined(BENCHMARK)  && !defined(DEBUG_TIMINGS)
         call second_v(time_set_counting_parameters_stop)
         if (use_mpi) time_set_counting_parameters_stop(2)=MPI_Wtime()
         time_set_counting_parameters_stop(:) = time_set_counting_parameters_stop(:) - code_start_time(:)
#endif
!c-----------------------------------------------------------------------
!c     start of cycle loop:
!c     Attention: In contrast to the original promet-version here
!c     one cycle is identical with the largest time step!
!c-----------------------------------------------------------------------
#ifndef DEBUG_TIMINGS

      ! calculation starts in this loop
      ! time until here was spent in seting up the code
      call second_v (calculation_start_time)
      if (use_mpi) calculation_start_time(2)=MPI_Wtime()
      calculation_start_time(:) = calculation_start_time(:) - code_start_time(:)
#endif


      do ncycl = nbegin+1, config%nend


#ifndef DEBUG_TIMINGS
         call initialize_subtiming_counters
#endif

#if defined(BENCHMARK) && !defined(DEBUG_TIMINGS)

         ! benchmark timings
         call printit_taskX(0,"ncycl =",ncycl,nbegin)

         if (ncycl .eq. 2) then
            ! first cycle is done
            ! store the time needed until now

            call second_v (time_first_cycle)
            if (use_mpi) time_first_cycle(2)=MPI_Wtime()
            ! time used for the first cycle (we passed through the loop already
            ! once)
            time_first_cycle(:) = time_first_cycle(:) - code_start_time(:)

            call printit_taskX(0,"first cycle done")

         endif


#endif /* BENCHMARK */


#ifdef CHECK_MEMORY
         meminfo_flag=.false. ! we do not want to check memory in every cycle
#endif
#ifndef DEBUG_TIMINGS
         ! we do not want to measure the time in every cycle
         cputime_flag =.false.
#endif

         call print_detailed_timestep_info(ncycl, nbegin)



         if (check_for_security_file() ) then

            call printit_taskX(0," ")
            call printit_taskX(0,"encountered fastexit file")
            call printit_taskX(0," ")

            goto 200
         endif

         if (time + hydro%dt .gt. config%tmax) then

            call printit_taskX(0," ")
            call printit_taskX(0,"time + dt = ", time + hydro%dt)
            call printit_taskX(0," > tmax = ", config%tmax)
            call printit_taskX(0," ")

            goto 200
         endif


         if (config%stop_bounce .gt. 0) then
            if (stop_entropy()) then
               call printit_taskX(0," ")
               call printit_taskX(0,"entropy reached 3.2")
               call printit_taskX(0," ")
               goto 200
            endif
         endif

         if (config%max_central_density .gt. 0) then
            if (core_density() .gt. config%max_central_density) then
               call printit_taskX(0," ")
               call printit_taskX(0,"central density reached ",config%max_central_density)
               call printit_taskX(0," ")
               goto 200
            endif
         endif

         if (time - starttime > config%bad_resolution_grace_period .and. config%stop_bad_resolution .gt. 0) then
            if (check_resolution()) then
               call printit_taskX(0," ")
               call printit_taskX(0,"resolution criterium > 0.6")
               call printit_taskX(0," ")
               goto 200
            endif
         endif

         if (config%spherical_neutron_star .gt. 0) then
            call spherical_neutron_star
         endif

#ifdef OPTIMAL_13_SWITCH
         call optimal_13_switch
#endif /* OPTIMAL_13_SWITCH */

#ifdef WRITE_BINARY_OUTPUT
!c-----------------------------------------------------------------------
!c     output:
!c-----------------------------------------------------------------------

         if (.not.(config%use_deactivate_output)) then
            if (config%collapse_time_resolve .eq. 1) then
               if (core_density() .gt. next_rho_c) then
                  next_rho_c = next_rho_c * next_rho_c_factor
                  call write_output_files
               endif
            endif

            ! If write_output_files() was called above, nout1 and tou1
            ! are already reset to 0, so it is save to check here again.
            if ((nout1 .ge. config%nout) .or. (tout1 .ge. config%tout)) then
               call write_output_files
            endif

            if (nrst  .ge. config%nrstrt .or. trst  .ge. config%trstrt) then
               call close_files   ! close output file (with old suffix)
               call set_filenames(99,99)   ! increment suffix
               call restrt(0)
               call open_files   ! open pltout file with new suffix
            endif
         endif
#endif


#ifndef DEBUG_TIMINGS
         cputime_flag=.true.
         call SECOND_V (cycle_start_time)           ! measure CPU time
         if (use_mpi) cycle_start_time(2)=MPI_Wtime()

         ! now the actuall cyce starts
         cycle_start_time(:) = cycle_start_time(:) - code_start_time(:)
#endif
         nstep  = ncycl               ! update time step number


         igrav = 0

!c-----------------------------------------------------------------------
!c     hydro-/ transport-step:
!c-----------------------------------------------------------------------
         ti_cyc = time                ! time at the beginning of cycle
         call rady(dt_cyc,ti_cyc, rady_self, rady_children)
         ! calculate all the grid areas
         hydro%dt     = dt_cyc              ! with multiple time stepping
                                      ! dt = l a t e s t used time step
         time   = ti_cyc              ! time at the end of the cycle
!c-----------------------------------------------------------------------
!c     save some quantities at every timestep:
!c-----------------------------------------------------------------------
#ifndef NOTRA
         if (config%p_ntr .ne. 0) call print_timestep_evolution
#endif
!c check consistency after first timestep
         if (nstep .eq. nbegin+2) then
            call check_run_consistency(config%calculation)
         endif

!c-----------------------------------------------------------------------
!c     output/restart counters have to be updated:
!c-----------------------------------------------------------------------

         nout1  = nout1  + 1
         tout1  = tout1 + hydro%dt
         nrst   = nrst  + 1
         trst   = trst  + hydro%dt

         if (time .gt. config%tmax)  then

            call printit_taskX(0," ")
            call printit_taskX(0,"time + dt = ", time + hydro%dt)
            call printit_taskX(0," > tmax = ", config%tmax)
            call printit_taskX(0," ")

           goto 200
         endif

!c-----------------------------------------------------------------------
!c     check the remaining CPU time and stop if necessary:
!c-----------------------------------------------------------------------

#ifndef DEBUG_TIMINGS
         cputime_flag = .true.
         call SECOND_V (cycle_end_time)
         if (use_mpi) cycle_end_time(2)=MPI_Wtime()

         ! cycle is done
         cycle_end_time(:) = cycle_end_time(:) - code_start_time(:)
         code_run_time(:) = cycle_end_time(:)
#endif

#ifndef DEBUG_TIMINGS
         timer%cycle_tot = cycle_end_time - cycle_start_time
#endif
#ifdef BENCHMARK
         time_second_cycle = timer%cycle_tot
#endif
         if (ncycl.le.nbegin+6) then
            call timings_out
         endif

#ifdef CRAY
         call tremain(cpur)
         time_remain(1) = cpur
         time_remain(2) = config%cpumax - code_run_time(2)
#elif NEC_COMPILER
         call tremain(cpur)
         time_remain(1) = cpur
         time_remain(2) = config%cpumax - code_run_time(2)

!c  on RZG Himiko tremain gives 0.0 for interactive jobs
!c  so let's allow for 600s instead
         if (cpur .le. EPSILON(cpur)) cpur=600._rk

#else
         ! remaining time can be easily computed from code_run_time
         ! and given time limits
         time_remain(1) = config%cpumax - code_run_time(1)
         time_remain(2) = config%cpumax - code_run_time(2)


#endif
         select case(config%wcrem)
         case(0)
#ifndef DEBUG_TIMINGS
            ! exit if the code_run_time is approaching
            ! specified wall clock time limit (and define a safety margin of
            ! the duration of 3 average cycles

            ! exit also if the remaining time is less than the time
            ! one needs to compute 5 cycles

            if ( code_run_time(1) .ge. config%cpumax-3._rk*timer%cycle_tot(1) .or. &
                 time_remain(1) .lt. 5._rk*timer%cycle_tot(1))  then

               call printit_taskX(0," ")
               call printit_taskX(0,"CPU jobtime limit reached ! remaining:", time_remain(1))
               call printit_taskX(0," ")

               goto 200
            endif
#endif
         case(1)
#ifndef DEBUG_TIMINGS
            ! exit if the code_run_time is approaching
            ! specified wall clock time limit (and define a safety margin of
            ! the duration of 3 average cycles

            ! exit also if the remaining time is less than the time
            ! one needs to compute 5 cycles

            if ( code_run_time(2) .ge. config%cpumax-3._rk*timer%cycle_tot(2) .or. &
                 time_remain(2) .lt. 5._rk*timer%cycle_tot(2))  then

               call printit_taskX(0," ")
               call printit_taskX(0,"WC  jobtime limit reached ! remaining:", time_remain(2))
               call printit_taskX(0," ")

               goto 200
            endif
#endif
         case default
           raise_abort("vertex(): no valid jobtime limit specified")
         end select

         filpos = 'append'
!c-----------------------------------------------------------------------
!c     end of cycle loop:
!c-----------------------------------------------------------------------

      enddo
!c-----------------------------------------------------------------------
!c     end of calculation:
!c-----------------------------------------------------------------------

200  continue

      call printit_taskX(0," ")
      call printit_taskX(0,"Next cycle ", ncycl)
      call printit_taskX(0,"Cycles done",ncycl - nbegin-1)
      call printit_taskX(0," ")


#ifndef DEBUG_TIMINGS
      cputime_flag = .true.
      call SECOND_V (calculation_stop_time)
      if (use_mpi) calculation_stop_time(2)=MPI_Wtime()
      calculation_stop_time(:) = calculation_stop_time(:) - code_start_time(:)
#endif


#ifdef WRITE_BINARY_OUTPUT
      if (.not.(config%use_deactivate_output)) then

#ifndef DEBUG_TIMINGS
         ! read transport file; use nevertheless the "write" variable for
         ! timing
         call second_v(write_output_start_time)
         if (use_mpi)  write_output_start_time(2)=MPI_Wtime()
         write_output_start_time(:) = write_output_start_time(:) - code_start_time(:)
#endif
         call write_output_files
         call close_files
         call set_filenames(99,99)         ! increment suffix
         call restrt(0)

#ifndef DEBUG_TIMINGS
         call second_v(write_output_stop_time)
         if (use_mpi)  write_output_stop_time(2)=MPI_Wtime()
         write_output_stop_time(:) = write_output_stop_time(:) - code_start_time(:)

         call printit_taskX(0," ")
         call printit_taskX(0," Writing of final output:         ", &
              (write_output_stop_time(:)-write_output_start_time(:)))
         call printit_taskX(0," ")

#ifdef BENCHMARK
         time_write_final_output(:) = (write_output_stop_time(:)-write_output_start_time(:))
#endif

#endif /* DEBUG_TIMINGS */

      endif
#endif /* WRITE_BINARY_OUTPUT */

!c----------------------------------------------------------------------
!c     write security file:
!c----------------------------------------------------------------------

      if (config%produc .eq. 1) then
         call io_security_file("write")
      endif

!c      cpucyc=(cyctime1-cyctime0)/real(ncycl-nbegin+1)
!c      cpuphy=(cyctime1-cyctime0)/max(trst,1e-12)
      ithreads=1
#if defined(OPENMP_HYD) && (defined(OPEN_MP_1D) || defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP  PARALLEL
!$    ithreads=OMP_GET_NUM_THREADS()
!$OMP  END PARALLEL
#endif


#ifndef DEBUG_TIMINGS

      call printit_taskX(0,"vertex> average CPU,WC time [s] per cycle:             ", &
        (calculation_stop_time-calculation_start_time)/max(real((ncycl-nbegin-1),kind=rk),1.0_rk))

      call printit_taskX(0,"vertex> total   CPU,WC time for job  :             ",(calculation_stop_time-calculation_start_time))
      call printit_taskX(0,"vertex> #threads",ithreads)
      call printit_taskX(0,"vertex> CPU/WC: ", &
        (calculation_stop_time(1)-calculation_start_time(1)) / (calculation_stop_time(2)-calculation_start_time(2)))
#endif

#if defined(BENCHMARK) && !defined(DEBUG_TIMINGS)
      call print_benchmark_summary(config%irstrt, nbegin, ncycl, calculation_start_time, calculation_stop_time, time_first_cycle, time_second_cycle)
#endif /* BENCHMARK */


      ! now we can allocate the arrays from module variables both in
      ! the case of MPI/OpenMP or only OpenMP mode

!      call deallocate_hydro_memory

#ifdef WRITE_BINARY_OUTPUT
      if (.not.(config%use_deactivate_output)) then
         call uninit_io
      endif
#endif

      if (config%jvisc .eq. 0) then
         call printit_taskX(0," ")
         call printit_taskX(0," CAUTION: Neutrino viscosity is switched off!!")
         call printit_taskX(0," This should ONLY be done ~70 ms postbounce !!")
         call printit_taskX(0," ")
      endif
      call printit_taskX(0," ")
      call printit_taskX(0,"lsrolo was used as ",lsrolo)
      call printit_taskX(0,"OK: Vertex finished correctly!")

      if (myproc .eq. 0) then
         call version_short
      endif ! myproc

      CALL MPI_FINALIZE(ier)

    end program vertex

!c=======================================================================
!c     end of vertex (main function)
!c=======================================================================
