!> \par driver for the hydro step
!>
!> \author A. Marek, MPA, Nov. 2009
!>
!> \detail
!> This module provides the subroutines hydro, hydro_are, detect_sh_are,
!> accel_are ... which CAN be overloaded with the COCONUT or PROMETHEUS
!> version. CAN because, ugly as the program layout is at the moment
!> some of the respective CFC-Versions are not defined here, but in
!> the respective files in the coconut subdirectory. What a mess!!
!>
!> Neverthe less IF both versions exist here, then they are overloaded.
!> IF a different version for the MPI_HYDRO exists, then they are also
!> overloaded. Somebody should clean this up!! At the moment this is
!> the easiest and "cleanest" way, without changing all the code
!>
!> Here a list of the actuall overloaded functions is given and a pointer
!> to the location of the procedures that are not overloaded
!>
!>  interface-name    real-name    compiled when?
!>  hydro             hydro_CFC    CFC_TRANSPORT
!>  hydro             hydro_PROM   if not CFC_TRANSPORT
!>  savare            savare_CFC   CFC_TRANSPORT
!>  savare            savare_PROM  if not CFC_TRANSPORT
!>
!> Here is a list of the subroutines which are only defined in hydro_area_functions.F90
!> if the PROMETHEUS version is used. This either means that these
!> subroutines do NOT exist in the COCONUT version, or are defined
!> in the subdirectory coconut. Clean up!!
!>
!> subroutine-name                 compiled when?          location
!> hydro_are                       if not CFC_TRANSPORT    hydro_area_functions.F90
!> detect_sh_are                   if not CFC_TRANSPORT    hydro_area_functions.F90
!> accel_are                       if not CFC_TRANSPORT    hydro_area_functions.F90
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!> \endverbatim
!>
module hydro_interface

  use precision



  public :: do_hydro_steps


  interface hydro_substep

#ifdef CFC_TRANSPORT
     module procedure hydro_CFC
#else
     module procedure hydro_PROM
#endif


  end interface


contains

#ifdef CFC_TRANSPORT

!> \par CFC-Version - driver for hydro step
!>
!> \author B. Mueller
!>
!> \param  config%dt_min
!> \param  ti_cyc
!> \param  dt_nex
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!> \endverbatim
!==============================================================================

 SUBROUTINE hydro_CFC(dt_min, ti_cyc, dt_nex, selftime, childrentime)

!==============================================================================

!   USE arecon_hy, ONLY: are_nu,dt_are,dt_cfl,ti_are
   USE gfloat_hy
   USE ioflx, ONLY: sumioe,sumion
   USE massio_hy
   USE nutrio_hy, ONLY: qentot,qyetot,sumdenu,sumdynu

   USE totgrq_hy, ONLY: sumdegr
   USE precision

   USE size_cfc, ONLY: n_s,n_e,o_s,o_e

   USE sum_fluxes
   USE coconut, only : coconut_hydro
   USE nusource_mod

   use hydro_areas_mod
   use cputim
   use configure
   IMPLICIT NONE

   real(kind=rk), intent(out) :: selftime(2), childrentime(2)
   real(kind=rk)              :: selftime_start(2)
   real, intent(in) :: dt_min, ti_cyc, dt_nex
   integer :: ia,i,j,k,in
   real, dimension(size(areas%dt_cfl)) :: dtnew
   integer, save :: isdtot=1
   integer :: ipos(0:config%qx)
   real(kind=rk) :: nusource_self(2), nusource_children(2)

  selftime     = 0._rk
  childrentime = 0._rk
#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

   do ia = 1, areas%are_nu         ! all areas with small time step
      areas%dt_are(ia) = dt_min
      areas%ti_are(ia) = ti_cyc + dt_nex ! time at the end of this step
   end do

   call coconut_hydro(dt_min,dtnew(1))      ! calculate area
   dtnew(1:areas%are_nu)=dtnew(1)

   areas%dt_cfl(1:areas%are_nu) = dtnew(1:areas%are_nu)

!-- time-integrate the neutrino sourceterms
   call nusource(dt_min,0, nusource_self, nusource_children)

   sumdenu=sumdenu+sumcq(qentot(1:config%qx,n_s:n_e,o_s:o_e),dt_min)
   sumdynu=sumdynu+sumcq(qyetot(1:config%qx,n_s:n_e,o_s:o_e,1),dt_min)
   sumdegr=0.0_rk

!#ifdef FCNC_CALC
!      sumdenu_n=sumdenu_n+
!     &           sumcq(qentot_n(1:config%qx,1:config%qy,1:config%qz),dt_min)
!      sumdynu_n=sumdynu_n+
!     &           sumcq(qyetot_n(1:config%qx,1:config%qy,1:config%qz),dt_min)
!#endif

   sumioe=sumioe+sumcqflux(eflxtot(2,n_s:n_e,o_s:o_e,areas%are_nu),dt_min)
   sumion=sumion+sumcqflux(xnflxtot(2,n_s:n_e,o_s:o_e,config%qn,areas%are_nu),dt_min)

#ifndef DEBUG_TIMINGS
 call second_v(selftime)
 selftime = selftime - selftime_start
#endif

   return

 END SUBROUTINE hydro_CFC

!==============================================================================



#else /* CFC_TRANSPORT2 */
!----------------------------------------------------------------------

!> \par PROMETHEUS-Version - driver for hydro step
!>
!> \author M. Rampp
!>
!> \param  dt_min
!> \param  ti_cyc
!> \param  dt_nex
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!> \endverbatim
!>
 subroutine hydro_PROM(dt_min, ti_cyc, dt_nex, selftime, childrentime)
   use precision
   use abort

   use filare_overload

!   use arecon_hy
   use totare_hy

   use sum_fluxes
   use intgrs_hy, only : igrav
   use nutrio_hy, only : sumdenu, sumdynu, qentot, qyetot, qmytot
   use totgrq_hy, only : sumdegr
   use massio_hy, only : eflxtot, dflxtot, xnflxtot
   use ioflx, only : sumioe, sumion
   use fluxcorr_overload

!== Calculate Hydro
   use gfloat_hy, ONLY : vlfrac
   use nucparam

   use mo_mpi

   use poisson_solver
   use interpol_eul
   use setinf_mod
   use hydro_area_functions
   use nusource_mod
   use tstep_mod

   use hydro_areas_mod
   use configure

   use cputim
   implicit none
! LOCAL variables that are not in modules

   real(kind=rk), intent(out)                     :: selftime(2), childrentime(2)
   real(kind=rk)                                  :: selftime_start(2),       &
                                                     setinf_self(2),          &
                                                     setinf_children(2),      &
                                                     detectshare_self(2),     &
                                                     detectshare_children(2), &
                                                     hydroare_self(2),        &
                                                     hydroare_children(2),    &
                                                     fluxcor_self(2),         &
                                                     fluxcor_children(2),     &
                                                     filare_self(2),          &
                                                     filare_children(2),      &
                                                     tstep_self(2),           &
                                                     tstep_children(2),       &
                                                     poisson_self(2),         &
                                                     poisson_children(2),     &
                                                     accelare_self(2),        &
                                                     accelare_children(2),    &
                                                     nusource_self(2),        &
                                                     nusource_children(2)

   real(kind=rk), intent(in)                      :: dt_min, ti_cyc, dt_nex
   integer(kind=ik)                               :: ia,i,j,k,in,ic,i_k,ij,irad
   real(kind=rk), dimension(size(areas%dt_cfl))   :: dtnew
   integer(kind=ik), save                         :: isdtot=1
   integer(kind=ik)                               :: ipos(0:config%qx)

   integer(kind=ik)                               :: ii1,ij1,ik1

   real(kind=rk)                                  :: ermflxtot(1:config%qy,1:config%qz)
   real(kind=rk)                                  :: tim1(2), tim2(2)


   selftime     = 0._rk
   childrentime = 0._rk

#ifndef DEBUG_TIMINGS
   call second_v(selftime_start)
#endif

!-- hydro without sourceterms for all areas
   dtnew(1:areas%are_nu)=9999._rk

   do ik1 = qz_s, qz_e
      do ij1 = qy_s, qy_e

         do irad = 1, config%qx
            dnsold(irad,ij1,ik1)=dentot(irad,ij1,ik1)
            vexold(irad,ij1,ik1)=vextot(irad,ij1,ik1)
            veyold(irad,ij1,ik1)=veytot(irad,ij1,ik1)
            vezold(irad,ij1,ik1)=veztot(irad,ij1,ik1)
         end do
      end do
   end do

   call setinf(setinf_self, setinf_children)
   childrentime = childrentime + setinf_self
   timer%hydro_setinf(:,1)=timer%hydro_setinf(:,1)+setinf_self
   timer%hydro_setinf_children(:,1)=timer%hydro_setinf_children(:,1)+setinf_children

!- multi-D shock detection for all areas
   do ia = 1, areas%are_nu
      areas%dt_are(ia) = dt_min
      call detect_sh_are(ia, detectshare_self, detectshare_children)
      childrentime = childrentime + detectshare_self
      timer%hydro_detectshare(:)=timer%hydro_detectshare(:)+ detectshare_self
      timer%hydro_detectshare_children(:)=timer%hydro_detectshare_children(:)+ detectshare_children
   enddo


   select case(isdtot)
   case (1)  ! forward-cycle
!-- x-sweep for all areas
      do ia = 1, areas%are_nu      ! all areas with small time step
         areas%dt_are(ia) = dt_min
         areas%ti_are(ia) = ti_cyc + dt_nex ! time at the end of this step
         call hydro_are(ia,1, hydroare_self, hydroare_children) ! calculate area
         timer%hydro_hydroare(:,1)=timer%hydro_hydroare(:,1)+hydroare_self
         childrentime = childrentime +  hydroare_self
         timer%hydro_hydroare_children(:,1)=timer%hydro_hydroare_children(:,1)+hydroare_self

! BEWARE: Here the overloaded subroutine filare is called. The overloading
! is done in module filare_xxx and is either the purely OPENMP filare_OPENMP
! or the hybrid MPI/OPENMP filare_MPI_HYDRO version

         call filare(ia, filare_self, filare_children)     ! fill area
         timer%hydro_filare(:,1)=timer%hydro_filare(:,1)+filare_self
         childrentime = childrentime + filare_self
         timer%hydro_filare_children(:,1)=timer%hydro_filare_children(:,1)+filare_children

         call tstep(dtnew(ia), tstep_self, tstep_children)
         timer%hydro_tstep(:,1)=timer%hydro_tstep(:,1)+tstep_self
         childrentime = childrentime + tstep_self
         timer%hydro_tstep_children(:,1)=timer%hydro_tstep_children(:,1)+tstep_children
      enddo

!-- flux correction for interfaces between areas
      call fluxcorr(fluxcor_self, fluxcor_children)
      timer%hydro_fluxcor(:,1)=timer%hydro_fluxcor(:,1)+fluxcor_self
      childrentime = childrentime + fluxcor_children
      timer%hydro_fluxcor_children(:,1)=timer%hydro_fluxcor_children(:,1)+fluxcor_children

!-- update boundaries of interfaces between areas
      call setinf(setinf_self, setinf_children)
      childrentime = childrentime + setinf_self
      timer%hydro_setinf(:,1)=timer%hydro_setinf(:,1)+setinf_self
      timer%hydro_setinf_children(:,1)=timer%hydro_setinf_children(:,1)+setinf_children


      if (config%nsdim .ge. 2)  then
!-- y-sweep for all areas
         do ia = 1, areas%are_nu   ! all areas with small time step
            areas%dt_are(ia) = dt_min
            areas%ti_are(ia) = ti_cyc + dt_nex

            call hydro_are(ia,2, hydroare_self, hydroare_children) ! calculate area
            timer%hydro_hydroare(:,2)=timer%hydro_hydroare(:,2)+hydroare_self
            childrentime = childrentime +  hydroare_children
            timer%hydro_hydroare_children(:,2)=timer%hydro_hydroare_children(:,2)+hydroare_self

! BEWARE: Here the overloaded subroutine filare is called. The overloading
! is done in module filare_xxx and is either the purely OPENMP filare_OPENMP
! or the hybrid MPI/OPENMP filare_MPI_HYDRO version

            call filare(ia, filare_self, filare_children)     ! fill area
            timer%hydro_filare(:,1)=timer%hydro_filare(:,1)+filare_self
            childrentime = childrentime + filare_self
            timer%hydro_filare_children(:,1)=timer%hydro_filare_children(:,1)+filare_children

            call tstep(dtnew(ia), tstep_self, tstep_children)
            timer%hydro_tstep(:,1)=timer%hydro_tstep(:,1)+tstep_self
            childrentime = childrentime + tstep_self
            timer%hydro_tstep_children(:,1)=timer%hydro_tstep_children(:,1)+tstep_children

         enddo
      endif

      if (config%nsdim .ge. 3)  then
!-- z-sweep for all areas
         do ia = 1, areas%are_nu   ! all areas with small time step
            areas%dt_are(ia) = dt_min
            areas%ti_are(ia) = ti_cyc + dt_nex

            call hydro_are(ia,3, hydroare_self, hydroare_children) ! calculate area
            timer%hydro_hydroare(:,3)=timer%hydro_hydroare(:,3)+hydroare_self
            childrentime = childrentime +  hydroare_children
            timer%hydro_hydroare_children(:,3)=timer%hydro_hydroare_children(:,3)+hydroare_self

            call filare(ia, filare_self, filare_children)     ! fill area
            timer%hydro_filare(:,1)=timer%hydro_filare(:,1)+filare_self
            childrentime = childrentime + filare_self
            timer%hydro_filare_children(:,1)=timer%hydro_filare_children(:,1)+filare_children

            call tstep(dtnew(ia), tstep_self, tstep_children)
            timer%hydro_tstep(:,1)=timer%hydro_tstep(:,1)+tstep_self
            childrentime = childrentime + tstep_self
            timer%hydro_tstep_children(:,1)=timer%hydro_tstep_children(:,1)+tstep_children

         enddo
      endif

   case (2)  ! backward-cycle

      if (config%nsdim .ge. 3)  then
!-- z-sweep for all areas
         do ia = 1, areas%are_nu   ! all areas with small time step
            areas%dt_are(ia) = dt_min
            areas%ti_are(ia) = ti_cyc + dt_nex

            call hydro_are(ia,3, hydroare_self, hydroare_children) ! calculate area
            timer%hydro_hydroare(:,3)=timer%hydro_hydroare(:,3)+hydroare_self
            childrentime = childrentime +  hydroare_children
            timer%hydro_hydroare_children(:,3)=timer%hydro_hydroare_children(:,3)+hydroare_self

            call filare(ia, filare_self, filare_children)     ! fill area
            timer%hydro_filare(:,2)=timer%hydro_filare(:,2)+filare_self
            childrentime = childrentime + filare_self
            timer%hydro_filare_children(:,2)=timer%hydro_filare_children(:,2)+filare_children


            call tstep(dtnew(ia), tstep_self, tstep_children)
            timer%hydro_tstep(:,2)=timer%hydro_tstep(:,2)+tstep_self
            childrentime = childrentime + tstep_self
            timer%hydro_tstep_children(:,2)=timer%hydro_tstep_children(:,2)+tstep_children

         enddo
      endif

      if (config%nsdim .ge. 2)  then
!-- y-sweep for all areas
         do ia = 1, areas%are_nu   ! all areas with small time step
            areas%dt_are(ia) = dt_min
            areas%ti_are(ia) = ti_cyc + dt_nex

            call hydro_are(ia,2, hydroare_self, hydroare_children) ! calculate area
            timer%hydro_hydroare(:,2)=timer%hydro_hydroare(:,2)+hydroare_self
            childrentime = childrentime +  hydroare_children
            timer%hydro_hydroare_children(:,2)=timer%hydro_hydroare_children(:,2)+hydroare_self

            call filare(ia, filare_self, filare_children)     ! fill area
            timer%hydro_filare(:,2)=timer%hydro_filare(:,2)+filare_self
            childrentime = childrentime + filare_self
            timer%hydro_filare_children(:,2)=timer%hydro_filare_children(:,2)+filare_children

            call tstep(dtnew(ia), tstep_self, tstep_children)
            timer%hydro_tstep(:,2)=timer%hydro_tstep(:,2)+tstep_self
            childrentime = childrentime + tstep_self
            timer%hydro_tstep_children(:,2)=timer%hydro_tstep_children(:,2)+tstep_children

         enddo
      endif

!-- update boundaries of interfaces between areas
      call setinf(setinf_self, setinf_children)
      childrentime = childrentime + setinf_self
      timer%hydro_setinf(:,2)=timer%hydro_setinf(:,2)+setinf_self
      timer%hydro_setinf_children(:,2)=timer%hydro_setinf_children(:,2)+setinf_children

!-- x-sweep for all areas
      do ia = 1, areas%are_nu      ! all areas with small time step
         areas%dt_are(ia) = dt_min
         areas%ti_are(ia) = ti_cyc + dt_nex ! time at the end of this step

         call hydro_are(ia,1, hydroare_self, hydroare_children) ! calculate area
         timer%hydro_hydroare(:,1)=timer%hydro_hydroare(:,1)+hydroare_self
         childrentime = childrentime +  hydroare_children
         timer%hydro_hydroare_children(:,1)=timer%hydro_hydroare_children(:,1)+hydroare_self

         call filare(ia, filare_self, filare_children)     ! fill area
         timer%hydro_filare(:,2)=timer%hydro_filare(:,2)+filare_self
         childrentime = childrentime + filare_self
         timer%hydro_filare_children(:,2)=timer%hydro_filare_children(:,2)+filare_children

         call tstep(dtnew(ia), tstep_self, tstep_children)
         timer%hydro_tstep(:,2)=timer%hydro_tstep(:,2)+tstep_self
         childrentime = childrentime + tstep_self
         timer%hydro_tstep_children(:,2)=timer%hydro_tstep_children(:,2)+tstep_children

      enddo

!-- flux correction for x-interfaces between areas
      call fluxcorr(fluxcor_self, fluxcor_children)
      timer%hydro_fluxcor(:,2)=timer%hydro_fluxcor(:,2)+fluxcor_self
      childrentime = childrentime + fluxcor_children
      timer%hydro_fluxcor_children(:,2)=timer%hydro_fluxcor_children(:,2)+fluxcor_children

   case default
      raise_abort("hydro(): nocase")
   end select

   areas%dt_cfl(1:areas%are_nu) = dtnew(1:areas%are_nu)

!-- accel and eos3d for all areas

!     for non Eulerian grids gpoold is interpolated to the
!      positions where gpotot is evaluated
   if (maxval(abs(ugrtot(:))).gt.0._rk) then
      xzrnew(0)       = xzltot(1)
      xzrnew(1:config%qx) = xzrtot(1:config%qx)


      do k= qz_s, qz_e
         do j = qy_s-1, qy_e

            call sort_vec(xzrold(0:),config%qx+1,1,xzrnew(0:),ipos(0:),config%qx+1,1)
            ipos(:)=ipos(:)-1
            call monintp(xzrold(0:),gpotot(0:,j,k),config%qx+1,1,ipos(0:), &
                         xzrnew(0:),gpoold(0:,j,k),config%qx+1,1,99,99)
         enddo
      enddo
   else

      gpoold(0:config%qx,qy_s-1:qy_e,qz_s:qz_e) = gpotot(0:config%qx,qy_s-1:qy_e,qz_s:qz_e)

   endif

   igrav=0
!   call poisson(poisson_self, poisson_children)
   timer%hydro_poisson(:)=timer%hydro_poisson(:)+poisson_self
   childrentime = childrentime + poisson_self
   timer%hydro_poisson_children(:)=timer%hydro_poisson_children(:)+poisson_children


   do ia = 1, areas%are_nu

      call accel_are(ia, accelare_self, accelare_children)    ! 'accel' area
      timer%hydro_accelare(:)=timer%hydro_accelare(:)+accelare_self
      childrentime = childrentime + accelare_self
      timer%hydro_accelare_children(:)=timer%hydro_accelare_children(:)+accelare_children

    ! fill area
      call filare(ia, filare_self, filare_children)     ! fill area
      timer%hydro_filare(:,2)=timer%hydro_filare(:,2)+filare_self
      childrentime = childrentime + filare_self
      timer%hydro_filare_children(:,2)=timer%hydro_filare_children(:,2)+filare_children

   enddo


   do ik1 = qz_s, qz_e
      do ij1 = qy_s, qy_e

         do irad = 1, config%qx
            acxtot(irad,ij1,ik1) = (vextot(irad,ij1,ik1) - vexold(irad,ij1,ik1))/dt_min
         end do
      end do
   end do

!-- change sweep sequence
   if (config%nsdim .ge. 2 .and. config%qy .ge. 4)  then
      if(isdtot .eq. 1) then
         isdtot = 2
      else
         isdtot = 1
      end if
   else
      isdtot = 1
   endif

   areas%ix_are(:,10) = isdtot

!-- time-integrate the neutrino sourceterms
   call nusource(dt_min,0, nusource_self, nusource_children)
   timer%hydro_nusource(:)=timer%hydro_nusource(:)+nusource_self
   childrentime = childrentime + nusource_self
   timer%hydro_nusource_children(:)=timer%hydro_nusource_children(:)+nusource_children



   sumdenu=sumdenu+sumcq(qentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                        +qmytot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                     *veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e),dt_min)
   sumdynu=sumdynu+sumcq(qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1),dt_min)
   sumdegr=sumdegr+sumcq(qgrtot(1:config%qx,qy_s:qy_e,qz_s:qz_e),dt_min)


   if (config%restmass_version .gt. 0) then

      ermflxtot(qy_s:qy_e,qz_s:qz_e) = eflxtot(2,qy_s:qy_e,qz_s:qz_e,areas%are_nu) &
                                 -dflxtot(2,qy_s:qy_e,qz_s:qz_e,areas%are_nu)*moffs
   endif

   if (config%restmass_version .eq. 1) then
      do j = qy_s, qy_e
         do k= qz_s, qz_e
            do in=1,config%qn-1
               ic=in
!#if LOW_EOS_NSE_FLAG != 22
               if (ic.gt.n_he4) ic=n_rep
!#else
!               if (ic.gt.n_he4 .or. ic .eq. n_he3 .or. ic .eq. n_tri
!     *             .or. ic .eq. n_d ) ic=n_rep
!#endif
               ermflxtot(j,k)=ermflxtot(j,k)+xnflxtot(2,j,k,in,areas%are_nu)*mbar(ic)
            enddo
         enddo
      enddo
   endif

   if (config%restmass_version .eq. 2) then
      do j = qy_s, qy_e
         do k= qz_s, qz_e
            do in=1,config%qn-1
               ic=in
               ermflxtot(j,k)=ermflxtot(j,k)+xnflxtot(2,j,k,in,areas%are_nu)*mbar(ic)
            enddo
         enddo
      enddo
   endif


   if (config%restmass_version .eq. 3) then
      do j = qy_s, qy_e
         do k= qz_s, qz_e
            do in=1,config%qn
               ic=in
               ermflxtot(j,k)=ermflxtot(j,k)+xnflxtot(2,j,k,in,areas%are_nu)*mbar(ic)
            enddo
         enddo
      enddo
   endif

  if (config%restmass_version .gt. 0) then
   sumioe=sumioe+ &
            sumcqflux(ermflxtot(qy_s:qy_e,qz_s:qz_e),dt_min)
  else
     sumioe=sumioe+ &
          sumcqflux(eflxtot(2,qy_s:qy_e,qz_s:qz_e,areas%are_nu),dt_min)
  endif

   sumion=sumion+ &

        sumcqflux(xnflxtot(2,qy_s:qy_e,qz_s:qz_e,config%qn,areas%are_nu),dt_min)


#ifndef DEBUG_TIMINGS
   call second_v(selftime)
   selftime = selftime - selftime_start
#endif

   return

end subroutine hydro_PROM

#endif /* CFC_TRANSPORT */

#ifndef CFC_TRANSPORT
#undef DEBUG
subroutine check_restrt_dt

   use totare_hy
   use bndinf_hy
   use hlpare_hy
!   use arecon_hy
   use gfloat_hy
   use intgrs_hy
   use vnew_hy
   use mesh_hy
   use mo_mpi
   use tstep_mod

   use hydro_areas_mod
   use configure
   use state
   implicit none
! LOCAL variables that are not in modules

   integer(kind=ik)                             :: ia, iarea, isd, i, ii
   real(kind=rk), dimension(size(areas%dt_cfl)) :: dtnew
   real(kind=rk)                                :: dt_old
   integer(kind=ik)                             :: j_end, k_end
   real(kind=rk)                                :: tstep_self(2), tstep_children(2)
   dtnew(1:areas%are_nu)=9999._rk

!> determine tstep for all areas
   do ia = 1, areas%are_nu      ! all areas with small time step
      iarea = ia
      are_id = iarea
      ixi    = areas%ix_are(iarea, 1)
      ixf    = areas%ix_are(iarea, 2)
      iox    = areas%ix_are(iarea, 3)
      iyi    = areas%ix_are(iarea, 4)
      iyf    = areas%ix_are(iarea, 5)
      ioy    = areas%ix_are(iarea, 6)
      izi    = areas%ix_are(iarea, 7)
      izf    = areas%ix_are(iarea, 8)
      ioz    = areas%ix_are(iarea, 9)
      isd    = areas%ix_are(iarea,10)
      !      ix_are(iarea,11) = 1

      hydro%dt     = areas%dt_are(iarea)
      dt_old       = hydro%dt

      areas%nz = izf - izi
      areas%ny = iyf - iyi
      areas%nx = ixf - ixi

      areas%nz = areas%nz/ioz + 1
      areas%ny = areas%ny/ioy + 1
      areas%nx = areas%nx/iox + 1

      if (ioy .eq. 1 .and. ioz .eq. 1) then
        j_end=qy_e
        k_end=qz_e
      else if (ioy .eq. config%qy .and. ioz .eq. config%qz) then
        j_end=qy_s
        k_end=qz_s
      end if

      velx  (1:areas%nx,qy_s:j_end,qz_s:k_end) = vextot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
      vely  (1:areas%nx,qy_s:j_end,qz_s:k_end) = veytot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
      velz  (1:areas%nx,qy_s:j_end,qz_s:k_end) = veztot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
      densty(1:areas%nx,qy_s:j_end,qz_s:k_end) = dentot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
      energy(1:areas%nx,qy_s:j_end,qz_s:k_end) = enetot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
      press (1:areas%nx,qy_s:j_end,qz_s:k_end) = pretot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
      gammae(1:areas%nx,qy_s:j_end,qz_s:k_end) = gaetot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
      gammac(1:areas%nx,qy_s:j_end,qz_s:k_end) = gactot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
      stot  (1:areas%nx,qy_s:j_end,qz_s:k_end) = stotot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
      temp  (1:areas%nx,qy_s:j_end,qz_s:k_end) = temtot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end)
      xnuc  (1:areas%nx,qy_s:j_end,qz_s:k_end,1:config%qn) = xnutot(ixi:ixf:iox,qy_s:j_end,qz_s:k_end,1:config%qn)

!-----------------------------------------------------------------------
!     Attention: For the z-direction igeomy = 5 (spherical, theta, phi) is
!     assumed! Since the resolution (ioz) may be greater than 1
!     the center of zone, the right boundaries of the zone and the
!     volume factors have to be recalculated!!!!
!-----------------------------------------------------------------------
      if (ioz .eq. 1) then
        zznl(1:areas%nz) = zzltot(izi:izf:ioz)
        zzn(1:areas%nz)  = zzntot(izi:izf:ioz)
        zznr(1:areas%nz) = zzrtot(izi:izf:ioz)
        dvz(1:areas%nz)  = dvztot(izi:izf:ioz)
      else
#ifdef DEBUG
        write(*,*) 'detect_sh_are:izi, nz,ioz ',izi, nz,ioz
#endif
        ii =  izi
        do i = 1, areas%nz
          zznl(i) = zzltot(ii)
          ii = ii + ioz
          if (ii .gt. config%qz+4) ii=1 !quick fix: ii could otherwise
          !easily exceed q if ioz=qz
        enddo

        do  i = 1, areas%nz
          zznr(i) = zznl(i+1)
        enddo

        do i = 1, areas%nz
          dvz(i) = zznr(i) - zznl(i) !igeomz = 5 is assumed!
          zzn(i) = 0.5_rk * (zznr(i) + zznl(i))
        enddo
      endif

!-----------------------------------------------------------------------
!     Attention: For the y-direction igeomy = 4 (spherical, theta) is
!     assumed! Since the resolution (ioy) may be greater than 1
!     the center of zone, the right boundaries of the zone and the
!     volume factors have to be recalculated!!!!
!-----------------------------------------------------------------------
      if(ioy .eq. 1) then
        yznl(1:areas%ny)   = yzltot(iyi:iyf:ioy)
        yzn(1:areas%ny)    = yzntot(iyi:iyf:ioy)
        yznr(1:areas%ny)   = yzrtot(iyi:iyf:ioy)
        dvy(1:areas%ny)    = dvytot(iyi:iyf:ioy)
      else
#ifdef DEBUG
        write(*,*) 'detect_sh_are:iyi, ny,ioy ',iyi, ny,ioy
#endif
        ii =  iyi
        do i = 1, areas%ny + 8
          yznl(i) = yzltot(ii)
          ii = ii + ioy
          if (ii .gt. config%qy+4) ii=1 !quick fix: ii could otherwise
          !easily exceed q if ioy=qy
        enddo

        do  i = 1, areas%ny + 7
          yznr(i) = yznl(i+1)
        enddo

        do i = 1, areas%ny + 7
          dvy(i) = cos(yznl(i)) - cos(yznr(i)) !igeomy = 4 is assumed!
          yzn(i) = 0.5_rk * (yznr(i) + yznl(i))
        enddo
      endif

!-----------------------------------------------------------------------
!     Attention: For the x-direction the resolution (iox) must be 1!
!     Otherwise the center of zone and boundaries of the zone are
!     different!!!
!-----------------------------------------------------------------------
      xznl(1:areas%nx)   = xzltot(ixi:ixf:iox)
      xzn(1:areas%nx)    = xzntot(ixi:ixf:iox)
      xznr(1:areas%nx)   = xzrtot(ixi:ixf:iox)

!> calculate now the timestep for this area
      call tstep(dtnew(ia), tstep_self, tstep_children)
   enddo

   areas%dt_cfl(1:areas%are_nu) = dtnew(1:areas%are_nu)

end subroutine check_restrt_dt

#endif /* CFC_TRANSPORT */


subroutine do_hydro_steps(dt_min, dt_nex, ti_cyc, idtmin, selftime, childrentime)

  use precision
  use configure
  use hydro_areas_mod
  use mo_mpi
  use cputim
  use grids
  use grid_mod

  implicit none
  real(kind=rk), intent(out)   :: dt_min
  real(kind=rk), intent(inout) :: dt_nex, ti_cyc
  integer(kind=ik), intent(in) :: idtmin
  real(kind=rk), intent(out)   :: selftime(2), childrentime(2)
  real(kind=rk), dimension(2)  :: selftime_start(2)

  integer(kind=ik)             :: ii, ierr

  real(kind=rk), dimension(2)  :: time_communication_start, &
                                  time_communication_end,   &
                                  grdvel_self,              &
                                  grdvel_children,          &
                                  mapra2hyd_self,           &
                                  mapra2hyd_children,       &
                                  hydrostep_self,           &
                                  hydrostep_children,       &
                                  lagrange_self,            &
                                  lagrange_children


  real(kind=rk)                :: dt_min_rcv

  selftime          = 0._rk
  childrentime      = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif


  do  ii = 1, areas%ix_are(idtmin,11)    ! number of small timesteps
     !         call setinf
     dt_min = min(areas%dt_cfl(idtmin), areas%dt_are(idtmin))
     ! - it is possible that dt becomes shorter due to CFL or TRANSP

     if (use_mpi) then
#ifndef DEBUG_TIMINGS
        call second_v(time_communication_start)
#endif

        !        MPI-Allreduce (Minimierung) fuer dt_min

        dt_min_rcv = 0._rk
        call MPI_AllReduce(dt_min, dt_min_rcv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
        dt_min= dt_min_rcv
#ifndef DEBUG_TIMINGS
        call second_v(time_communication_end)
        timer%transp_comm = timer%transp_comm + (time_communication_start-time_communication_end)
        childrentime = childrentime + time_communication_end - time_communication_start
        timer%hydro_total(:) = timer%hydro_total(:) + time_communication_end - time_communication_start


#endif
     endif ! use_mpi

     dt_nex = dt_nex + dt_min
     ! calculate therefore real dt_nex
#ifndef NOHYDRO
#ifndef CFC_TRANSPORT
     call grdvel(grdvel_self, grdvel_children)
     childrentime = childrentime + grdvel_self
     timer%hydro_total(:) = timer%hydro_total(:) + grdvel_self
#endif /*CFC_TRANSPORT*/

#ifndef NOTRA
     if (config%p_nbk .gt. 0) then
        call map_ra2hyd(.true.,mapra2hyd_self, mapra2hyd_children)
        childrentime = childrentime + mapra2hyd_self
        timer%hydro_total(:) = timer%hydro_total(:) + mapra2hyd_self
     endif

#endif /* NOTRA */

     call hydro_substep(dt_min,ti_cyc, dt_nex, hydrostep_self, hydrostep_children)

     timer%rady_hydrostep(:)          = timer%rady_hydrostep(:)          &
          + hydrostep_self
     timer%rady_hydrostep_children(:) = timer%rady_hydrostep_children(:) &
          + hydrostep_children

     childrentime = childrentime + hydrostep_self
     timer%hydro_total(:) = timer%hydro_total(:) + hydrostep_self
#endif /*NOHYDRO*/


     !     advance the (lagrangian) transp-grid
#ifndef NOTRA
     if(config%p_ntr .gt. 0 .and. config%ieul .eq. 0) then
        call lagrange(dt_min, lagrange_self, lagrange_children)
        childrentime = childrentime + lagrange_self
        timer%hydro_total(:) = timer%hydro_total(:) + lagrange_self
     endif
#endif /*NOTRA*/

  enddo ! "ndt"-loop

#ifndef DEBUG_TIMINGS
  call second_v(selftime)
  selftime = selftime - selftime_start
#endif

end subroutine do_hydro_steps

end module hydro_interface
