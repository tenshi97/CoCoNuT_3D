module neutrino_mod
!== Calculate Neutrino transport
  public

  contains
  subroutine neutrino(dt_nex, nrep, error, itnum, nenum, sigma, nstep, ti_cyc, &
                      labb, loutp, nray, ndt, dt_min, selftime, childrentime)

    use precision
    use mo_mpi
#ifndef NOTRA
    use qterms_mod
#endif
    use nusource_mod
    use cputim
    use state
    use configure
    use nutrio_hy
    use totare_hy
    use sum_fluxes
    use print_stdout_mod
    use savare_overload
    use grids
    implicit none
! LOCAL variables that are not in modules
    real(kind=rk), intent(out)   :: selftime(2), childrentime(2)
    integer(kind=ik), intent(inout) :: itnum, nenum, nstep, nray(2), ndt
    real(kind=rk), intent(inout)    :: sigma, ti_cyc, dt_min
    logical, intent(inout)          :: labb, loutp
    real(kind=rk)                :: selftime_start(2)
    real(kind=rk), dimension(2)  :: nusource_self, nusource_children, &
                                    qterms_self, qterms_children,     &
                                    savare_self, savare_children
    real(kind=rk)                :: mapra2hyd_self(2), mapra2hyd_children(2)

    logical, intent(out)         :: error
    real(kind=rk), intent(in)    :: dt_nex
    integer(kind=ik), intent(in) :: nrep

    selftime     = 0._rk
    childrentime = 0._rk
#ifndef DEBUG_TIMINGS
    call second_v(selftime_start)
#endif

    itnum = 0
    nenum = 0
    sigma = 0._rk
    error=.true.

#ifndef NOTRA
!  application (-) of the sourceterms at the OLD transp-timelevel
#ifndef NOHYDRO
    if (config%p_nbk .gt. 0) then

       call nusource(dt_nex,1, nusource_self, nusource_children)
       timer%neutrino_nusource = timer%neutrino_nusource + nusource_self
       timer%neutrino_nusource_children = timer%neutrino_nusource + nusource_children
       childrentime = childrentime + nusource_self


       sumdenu=sumdenu+sumcq(qentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
#ifdef NEW_QMO
                     +qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                     *vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
#endif /* NEW_QMO */
                     +qmytot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                     *veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e),dt_nex)

       sumdynu=sumdynu+ sumcq(qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1),dt_nex)
!#ifdef FCNC_CALC
!            sumdenu_n=sumdenu_n+
!     &           sumcq(qentot_n(1:config%qx,qy_s:qy_e,1:config%qz),dt_nex)
!            sumdynu_n=sumdynu_n+
!     &           sumcq(qyetot_n(1:config%qx,qy_s:qy_e,1:config%qz),dt_nex)
!#endif

    endif ! config%p_nbk
#endif /* NOHYDRO */

    ! TRANSPORT STEP:
    !     calculate the new neutrino quantities, in particular sourceterms
    !     (by solution of the radiative transfer equation done in qterms)
    ! TRANSPORT STEP:
    !     calculate the new neutrino quantities, in particular sourceterms
    !     (by solution of the radiative transfer equation done in qterms)

#ifdef CHECK_MEMORY
    if (meminfo_flag) call meminfo_start("total TRANSP")
#endif


   call qterms (nstep,ti_cyc,dt_nex,transport%dt,sigma,itnum,nenum,labb,.false.,loutp,nray, qterms_self, qterms_children)


   timer%neutrino_qterms = timer%neutrino_qterms + qterms_self
   timer%neutrino_qterms_children = timer%neutrino_qterms_children + qterms_children
   childrentime = childrentime + qterms_self

   !     no successful transport-timestep

   if(labb) then
      call printit_taskX(0,"arecon> error at qterms, transport timestep',' to be redone   !")


      call savare(0, savare_self, savare_children)
      timer%neutrino_savare = timer%neutrino_savare + savare_self
      timer%neutrino_savare_children = timer%neutrino_savare_children + savare_children
      childrentime = childrentime + savare_self

      !     cf Winkler&Norman: Num Radhydro p.117
      transport%dt=max(0.1_rk,0.8_rk-0.2_rk*(nrep+1))*transport%dt
!            ndt = max(int(transport%dt/dt_min),1)

      ndt=max(2*ceiling(0.5_rk*transport%dt/dt_min),2)
      dt_min=transport%dt/real(ndt,kind=rk)

#ifndef NOHYDRO
      if (config%p_nbk .gt. 0) then
         call nusource(dt_nex,99, nusource_self, nusource_children)
         timer%neutrino_nusource = timer%neutrino_nusource + nusource_self
         timer%neutrino_nusource_children = timer%neutrino_nusource_children + nusource_children
         childrentime = childrentime + nusource_self

      endif
#endif
      error=.true.
#ifndef DEBUG_TIMINGS
      call second_v(selftime)
      selftime = selftime - selftime_start
#endif

      return
   endif


#ifdef CHECK_MEMORY
   if (meminfo_flag) call meminfo_stop("total TRANSP")
#endif

   error=.false.

!  application (+) of the sourceterms at the NEW transp-timelevel

   call map_ra2hyd(.true.,mapra2hyd_self, mapra2hyd_children)
   timer%neutrino_mapra2hyd = timer%neutrino_mapra2hyd + mapra2hyd_self
   timer%neutrino_mapra2hyd_children = timer%neutrino_mapra2hyd_children + mapra2hyd_children
   childrentime = childrentime + mapra2hyd_self


#ifndef NOHYDRO
   if (config%p_nbk .gt. 0) then

      call nusource(dt_nex,2, nusource_self, nusource_children)
      timer%neutrino_nusource = timer%neutrino_nusource + nusource_self
      timer%neutrino_nusource_children = timer%neutrino_nusource_children  + nusource_children
      childrentime = childrentime + nusource_self

      sumdenu=sumdenu+ sumcq(qentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                      +qmytot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                      *veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e),dt_nex)

      sumdynu=sumdynu+sumcq(qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1),dt_nex)
#ifdef FCNC_CALC
      sumdenu_n=sumdenu_n+sumcq(qentot_n(1:config%qx,qy_s:qy_e,qz_s:qz_e),dt_nex)
      sumdynu_n=sumdynu_n+sumcq(qyetot_n(1:config%qx,qy_s:qy_e,qz_s:qz_e),dt_nex)
#endif /*FCNC_CALC*/

   endif

#endif
#endif
#ifndef DEBUG_TIMINGS
   call second_v(selftime)
   selftime = selftime - selftime_start
#endif

   return

 end subroutine neutrino
end module neutrino_mod
