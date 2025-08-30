module restart_model_mod

  implicit none

  contains

  subroutine restart_from_model(restart_mode, nbegin, dtt)

    use precision
    use abort
    use error
    use gfloat_hy, only : time,  tin, gravin, gamein, gamma, &
                          gamcin, xnin
    use charac, only    : filpos
    use intgrs_hy, only : nstep !, nend
!    use param_rt, only  : irst_ra
    !    use arecon_hy, only : nhystp
#ifdef WRITE_BINARY_OUTPUT
      use output_hydro, only : set_filenames, open_files, close_files, &
                               write_output_files

      use restart, only : restrt
#endif

      use mpi_domains
#ifndef NOTRA
      use bte_init
#endif

#ifdef CFC_TRANSPORT
  use coconut
  use gr_alloc
#endif

    use hydro_areas_mod
    use configure
    use state
    implicit none

    integer(kind=ik), intent(in)  :: restart_mode
    integer(kind=ik), intent(out) :: nbegin
    real(kind=rk), intent(out)    :: dtt

#ifdef WRITE_BINARY_OUTPUT
    if (.not.(config%use_deactivate_output)) then
       call set_filenames(restart_mode,config%irst_ra)
    endif
#endif
         dtt = hydro%dt

         filpos = 'append'

         if (config%irst_ra .eq. 0) then  ! from tra.par
#ifdef WRITE_BINARY_OUTPUT
            if (.not.(config%use_deactivate_output)) then
               call open_files
               call set_filenames(restart_mode,config%irst_ra)
               call restrt(1)
            endif
#endif

#ifndef NOTRA
            call init_transport
#endif
            time         = 0._rk
            hydro%dt    = config%dtini
            nstep = 0
            areas%nhystp= 0
            raise_error("promet(): irst_ra = 0")
!            call stopit('promet: irst_ra = 0',0)
         endif
#ifdef WRITE_BINARY_OUTPUT
         if (.not.(config%use_deactivate_output)) then
            call restrt(1)
         endif
#endif

#ifdef CFC_TRANSPORT
         call init_coconut(restart_mode)
#endif

         if (config%use_print_mpiDomain) then
            ! print MPI domain decomposition and thread affinity
            call print_mpi_setup
         endif
#ifdef WRITE_BINARY_OUTPUT
         if (.not.(config%use_deactivate_output)) then
            call open_files
         endif
!         call write_output_files
#endif

         nbegin = nstep
         config%nend   = nstep + config%nend

         tin    = 0._rk
         gravin = 0._rk
         gamein = gamma
         gamcin = gamma
         xnin(:) = 1._rk

  end subroutine restart_from_model


end module restart_model_mod
