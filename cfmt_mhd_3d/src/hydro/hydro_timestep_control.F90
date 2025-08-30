
module timestep_areas

  public

contains

  subroutine minimum_timestep_of_areas(idtmin, isdmin, dt_nex, dt_cyc, &
                                       dt_min, dt_hyd, ndt, selftime,  &
                                       childrentime)
    use precision
    use mo_mpi
    use configure
    use hydro_areas_mod
    use specfun, only : ISRMIN_V
    use state
    use cputim
    use print_stdout_mod
    use abort
    use hydro_interface

    implicit none

    integer(kind=ik), intent(inout) :: idtmin, isdmin, ndt
    real(kind=rk), intent(inout)    :: dt_nex, dt_cyc, dt_min, dt_hyd

    integer(kind=ik)                :: ndtnew, ndtadd, ndtold, idtnex
    logical                         :: restrt=.true.
    real(kind=rk), dimension(2)     :: time_communication_start, &
                                       time_communication_end

    integer(kind=ik)                :: ndt_rcv, ierr
    real(kind=rk)                   :: dtlong
    real(kind=rk), intent(out)      :: selftime(2), childrentime(2)
    real(kind=rk), dimension(2)     :: selftime_start(2)

    selftime          = 0._rk
    childrentime      = 0._rk

#ifndef DEBUG_TIMINGS
    call second_v(selftime_start)
#endif


#ifndef CFC_TRANSPORT
    if (restrt) then
       if (config%irstrt .eq. 2) then ! if restrt-parameters taken from ppm.par
                                      ! determine the time step length of both areas
                                      ! since only the one of the first area is known
                                      ! and the one of the second could be smaller
                                      ! and therefore determine the actual hydro timestep
          call check_restrt_dt
       endif
    endif
#endif /* CFC_TRANSPORT */

    idtmin = ISRMIN_V(areas%are_nu,areas%dt_cfl,1)  ! use CFL time scale
    isdmin = areas%ix_are(idtmin,10)       ! sweep direction of MIN-area
    if (restrt) then
       dt_hyd=areas%ix_are(idtmin,11)*areas%dt_are(idtmin)
       restrt=.false.
    endif

    idtnex = idtmin
    dt_nex = 1.0e33_rk
    dt_cyc = 1.0e33_rk
    !-----------------------------------------------------------------------
    !     determine new dt_min and number of small time step ndt:
    !     - increase dt_min and ndt only, if MIN and NEXT cycle starts!
    !     - areas%dt_are still contains last used time step!
    !     - ndt is always even ?
    !-----------------------------------------------------------------------

    select case (isdmin)
    case(1)
       if(areas%nhystp .le. 10) then
          dt_min=areas%dt_cfl(idtmin)
          ndt = 2
       else
          ndtold=areas%ix_are(idtmin,11)
          dtlong=min(transport%dt,dt_hyd,config%dtmaxx)
#ifdef RAND_DT_RAD
          call RANDOM_NUMBER(rand)
          dtlong=dtlong*(1._rk+2._rk*(rand-0.5_rk)*0.1_rk)
#endif
          ndtnew=2*ceiling(0.5_rk*dtlong/areas%dt_cfl(idtmin))
          ndtadd=max(2,2*nint(0.25*ndtold))
          ndt=min(ndtold+ndtadd,ndtnew)
          dt_min=dtlong/real(ndt,kind=rk)
       end if
    case(2)
       ndt    = areas%ix_are(idtmin,11)
       call printit_taskX(0,"arecon > W A R N I N G : ndt should be even here:",ndt)
    case default

       raise_abort("rady(): arecon > wrong isdmin")
    end select


    if (use_mpi) then
#ifndef DEBUG_TIMINGS
       call second_v(time_communication_start)
#endif
       ! STILL WITH IFDEF SINCE OVERLOADED ALLREDUCE ACCEPTS ONLY REALS!!!

       ! MPI-Allreduce (Minimierung) fuer dt_min


       ndt_rcv = 0
       call MPI_AllReduce(ndt, ndt_rcv, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
       ndt=ndt_rcv
#ifndef DEBUG_TIMINGS
       call second_v(time_communication_end)

       timer%transp_comm = timer%transp_comm + time_communication_end - time_communication_start

       childrentime = childrentime + time_communication_end - time_communication_start
#endif
    endif ! use_mpi

#ifndef DEBUG_TIMINGS
    call second_v(selftime)
    selftime = selftime - selftime_start
#endif
  end subroutine minimum_timestep_of_areas

end module timestep_areas
