module hydro_init_mod

  implicit none

  contains

    subroutine initialize_hydro(mode)
      use precision
      use abort
      use hydro_areas_mod
      use configure
      use state
      use movgrid_hy

      use eos_sn2
      use cputim
#if defined(HTCL)
      use htcl_hy
#endif
#ifdef BURN_NETWORK
      use burnrate
#endif
#ifdef TAK_RATES
      use rateselec
#endif

#ifdef WRITE_BINARY_OUTPUT
      use dataformat, only : init_io, uninit_io
#endif
      implicit none

      integer(kind=ik) , intent(in) :: mode

      select case (mode)

         case (1)
            ! first we have to set a few remaining config variables
            ! which were not directly readin in read parameter_files

            ! set a few bounary conditions

            ! can this be removed???
            config%latybc=config%bndmny
            config%latzbc=config%bndmnz

            ! initialize the timesteps

            hydro%dt         = config%dtini
            transport%dt     = 100._rk*config%dtmax

            call init_areas

            call init_movegrid

         case(2)
#ifdef WRITE_BINARY_OUTPUT
            if (.not.(config%use_deactivate_output)) then
               call init_io
            endif
#endif

#ifdef HTCL
            call init_htcl
#endif
            call init_eos


#ifdef BURN_NETWORK
            call rates
#endif
#ifdef TAK_RATES
            call get_cap
#endif
         case default
            raise_abort("Invalid mode")
      end select

    end subroutine initialize_hydro

end module hydro_init_mod
