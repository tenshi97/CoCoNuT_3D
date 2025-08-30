#undef DIAGNOSE_2D

module grids

implicit none

contains

#if !defined(NOTRA) && !defined(PROGRAM_remap)

!>
!> \verbatim
!> advance the lagrangean grid with radial velocity be*c
!>
!>  Author: M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
 subroutine lagrange(dt, selftime, childrentime)

   use precision
   use abort
   use radial_grid_rt
   use param_rt
   use multigrid_rt
   use backquants_rt
   use phycon
   use mo_mpi
   use intgrs_hy
   use himodyn_mod
   use cputim

   use configure
   implicit none

   real(kind=rk), intent(out)  :: selftime(2), childrentime(2)
   real(kind=rk), dimension(2) :: selftime_start(2), maphyd2ra_self(2), &
                                  maphyd2ra_children(2)

   integer(kind=ik)            :: i, k, j
   real(kind=rk)               :: pid4, dt

   selftime      = 0._rk
   childrentime = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

   pid4 = 0.25_rk/pc_pi

   do k= nzmoms, nzmome
      do j= nymoms, nymome

         if (config%use_multid_collapse) then
            do i=-1,config%imaxp +1
               r(i)=ralag_3d(i,j,k)
            enddo
         else
            do i=-1,config%imaxp +1
               r(i)=ralag_1d(i)
            enddo

         endif !multid_collapse

         call himo_r(0)

         call map_hyd2ra(j,k,config%i_grtr,0, maphyd2ra_self, maphyd2ra_children)
         childrentime = childrentime + maphyd2ra_self

         if (config%use_multid_collapse) then
            do i=0,config%imaxp +1
               ralag_3d(i,j,k)    = ralag_3d(i,j,k) + dt * be(i) * pc_cl
            enddo
            if (ralag_3d(-1,j,k) .le. 0.0_rk) then 
               write (*,*) "Task ",myproc
               write (*,*) 'lagrange> i=-1 zone:',i,j,k,ralag_3d(-1,j,k)
               !     call regrid
               raise_abort("lagrange(): grids not ordered")
               !            call stopit('*** grids not ordered',0)
            endif
         else
            do i=0,config%imaxp +1
               ralag_1d(i)    = ralag_1d(i) + dt * be(i) * pc_cl
            enddo
            if (ralag_1d(-1) .le. 0.0_rk) then 
               write (*,*) "Task ",myproc
               write (*,*) 'lagrange> i=-1 zone:',i,ralag_1d(-1)
               !     call regrid
               raise_abort("lagrange(): grids not ordered")
               !            call stopit('*** grids not ordered',0)
            endif
         endif ! multid_collapse

         if (config%use_multid_collapse) then
            do i=0,config%imaxp +1
               if (ralag_3d(i,j,k) .le. (1._rk+1e-7_rk)*ralag_3d(i-1,j,k)) then 
                  write (*,*) "Task ",myproc
                  write (*,*) 'lagrange> zones cross:',i,ralag_3d(i,j,k), &
                        ralag_3d(i-1,j,k)
               !     call regrid
                  raise_abort("lagrange(): grids not ordered")
               !               call stopit('*** grids not ordered',0)
               endif
            enddo
         else
            do i=0,config%imaxp +1
               if (ralag_1d(i) .le. (1._rk+1e-7_rk)*ralag_1d(i-1)) then 
                  write (*,*) "Task ",myproc
                  write (*,*) 'lagrange> zones cross:',i,ralag_1d(i), &
                        ralag_1d(i-1)
               !     call regrid
                  raise_abort("lagrange(): grids not ordered")
               !               call stopit('*** grids not ordered',0)
               endif
            enddo
         endif ! multid_collapse

      enddo
   enddo

#ifndef DEBUG_TIMINGS
   call second_v(selftime)
   selftime = selftime - selftime_start
#endif
   return 
 end subroutine lagrange
#endif /* not PROGRAM_remap or NOTRA */

#ifndef NOTRA 
! ---------------------------------------------------------------------
!>
!> \verbatim
!> maps (per linear interpolation) the eddington-factors 
!>      living on the ralag(:,j0,:),rqlag-Grid onto the transp-grid of
!>      the current angular ray j(,k)
!>
!>  Author: M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
 subroutine map_eddi(j,k,j0,selftime, childrentime)

   use precision
#ifndef NOTRA      
   use radial_grid_rt
!   use param_rt ! forcheck
   use multigrid_rt

   use averagrid_rt
   use radfield_rt
!   use matkhelp_rt ! forcheck
   use boundary_rt
#endif
   use interpol_eul

   use abort
   use cputim
   use configure
   implicit none
! LOCAL variables that are not in modules
   real(kind=rk), intent(out)    :: selftime(2), childrentime(2)
   real(kind=rk)                 :: selftime_start(2)
   integer(kind=ik)              :: i, iposq, ie, ipos, is
   integer(kind=ik), intent (in) :: j, k, j0

   integer(kind=ik)              ::  iri(-1:config%imaxp +1)
   real(kind=rk)                 :: scr, scrtchq, scrtch

   selftime     = 0._rk
   childrentime = 0._rk
#ifndef DEBUG_TIMINGS
   call second_v(selftime_start)
#endif

   if (j0 .ne. 0) then
      raise_abort("map_eddi(): something wrong")
   endif

!      call sort_grid(ralag(-1,j0,k),config%imaxp +1+2,r(-1),iri(-1),config%imaxp +1+2)
   if (config%use_multid_collapse) then
      call sort_vec(ralag_3d(-1:,j0,k),config%imaxp +1+2,1,r(-1:),iri(-1:),config%imaxp +1+2,1)
   else
      call sort_vec(ralag_1d(-1:),config%imaxp +1+2,1,r(-1:),iri(-1:),config%imaxp +1+2,1)
   endif !multid_collapse

   do is=1,config%isma
      do ie=1,config%iemax

         do i=-1,config%imaxp +1
               !     call LOCATE(ralag(-1,j0,k),config%imaxp +1+2,r(i),ipos)
               !     ipos=ipos-2
            ipos = iri(i)-2

            if ((ipos .lt. -1) .or. (ipos .ge. config%imaxp +1)) then
               scrtch = 1.0_rk
               ipos = min(ipos,config%imaxp) !extrapolation 
               ipos = max(ipos,-1) !extrapolation 
            else
               if (config%use_multid_collapse) then
               scrtch  = (r (i)-ralag_3d(ipos,j0,k))/       &
                    (ralag_3d(ipos+1,j0,k)-ralag_3d(ipos,j0,k))
               else
               scrtch  = (r (i)-ralag_1d(ipos))/       &
                    (ralag_1d(ipos+1)-ralag_1d(ipos))
               endif ! multid_collapse

            endif 

            if (i .ge. 0) then
               if (rq(i) .ge. rqlag(ipos+1)) then
                  iposq=ipos+1
               else
                  iposq=ipos
               endif
               if ((iposq .lt. 0) .or. (iposq .ge. config%imaxp +1)) then
                  scrtchq = 1.0_rk
                  iposq = min(ipos,config%imaxp) !extrapolation 
                  iposq = max(ipos,0) !extrapolation 
               else
                  scrtchq = (rq(i)-rqlag(iposq))/   &
                       (rqlag(iposq+1)-rqlag(iposq))
               endif
               xjit(i,ie,is) = xjlag(iposq,ie,is) + &
                    scrtchq *  (xjlag(iposq+1,ie,is) - xjlag(iposq,ie,is))
               fe(i,ie,is)   = felag(iposq,ie,is) + &
                    scrtchq *  (felag(iposq+1,ie,is) - felag(iposq,ie,is))
            endif
            xhit(i,ie,is) = xhlag(ipos,ie,is)+ &
                 scrtch * (xhlag(ipos+1,ie,is  )-xhlag(ipos,ie,is  ))
            ge(i,ie,is) = gelag(ipos,ie,is)+   &
                 scrtch * (gelag(ipos+1,ie,is  )-gelag(ipos,ie,is  ))
         enddo

         do i=0,config%imaxp
            fq(i,ie,is) = polrq1(i) * fe(i,ie,is) + &
                 polrq(i) * fe(i+1,ie,is)
         enddo

         if (config%use_multid_collapse) then
            scr = (r(config%imaxp +1)-rqlag(iposq+1))/(ralag_3d(ipos+1,j0,k)-rqlag(iposq+1))
         else
            scr = (r(config%imaxp +1)-rqlag(iposq+1))/(ralag_1d(ipos+1)-rqlag(iposq+1))
         endif ! multid_collapse
         fq(config%imaxp +1,ie,is) = felag(iposq+1,ie,is)+ &
              min(scr,1.0_rk)*(fql_rd(ie,is)-felag(iposq+1,ie,is))
         xm1rd(ie,is) = xjlag(iposq+1,ie,is)+     &
              min(scr,1.0_rk)*(xjlag(iposq+2,ie,is)-xjlag(iposq+1,ie,is))
         xm1rd(ie,is) = xhit(config%imaxp +1,ie,is)/xm1rd(ie,is)
         xm1rd(ie,is) = sign(min(abs(xm1rd(ie,is)),1.0_rk),xm1rd(ie,is))
         xhinnen(ie,is) = xhit(-1,ie,is)
      enddo
   enddo


      !      write(*,*) maxval(rq),maxval(r)
      !      write(*,*) xm1rd,maxval(fq),maxval(fe),maxval(ge)

#ifndef DEBUG_TIMINGS
   call second_v(selftime)
   selftime = selftime - selftime_start
#endif
   return
 end subroutine map_eddi
#endif /* NOTRA */


#if !defined(NOTRA) && !defined(PROGRAM_remap)
!> \verbatim
!>
!>
!>     task:      
!>           - map the neutrino SOURCETERMS and local energy-density, 
!>              flux, number-density of neutrinos
!>              from the transport- to the hydro-grid 
!>
!>     note: - MUST NOT AFFECT the uncovered regions of the hydro-grid
!>           - assumes those to be filled elsewhere
!>           
!>           - assumes Eulerian FRAME quantities input
!>           - assumes that the Neutrino-Transport is done in [MeV] 
!>                      and Hydro in [CGS]     
!>
!>     input:   nxlb(j),nxub(j) bounds of the covered regions
!>              selag,sylag,smlag    sourceterms in the Eulerian FRAME
!>              enulag,dnulag,fnulag,pnulag ( e-d, e-f, n-d, n-p ) 
!>              in the Eulerian FRAME
!>
!>     output:  qentot: sourceterm for energy-eq     [erg/cm^3/s]
!>              qyetot: sourceterm for Y_e-eq        [g/cm^3/s]
!>              qmotot: sourceterm for mom-eq        [erg/cm^4]
!>
!>              enutot: Energy-density (Inertial frame)  [erg/cm^3]
!>              fnutot: Energy-Flux    (Inertial frame)  [erg/cm^2/s]
!>              pnutot: Pressure       (Inertial frame)  [erg/cm^3]
!>
!>              dnutot: Number-density (Inertial frame)  [1  /cm^3]
!>              gnutot: Number-flux    (Inertial frame)  [1  /cm^2/s]
!>
!>  Author: M. Rampp and F. Hanke
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim

 subroutine map_ra2hyd(do_mpi_comm, selftime, childrentime)
   use precision
   use abort
   use phycon
   use intgrs_hy
   use nutrio_hy
   use radial_grid_rt
   use multigrid_rt
   use stressten_rt
   use sourceterms_rt
   use overlap
   use param_rt
   use neutrinotypes
   use theta_grid_rt
   use totare_hy
   use interpol_eul
   use himodyn_mod

#ifdef CFC_TRANSPORT2
   use metric_cfc, ONLY: sqrt_gamma,alpha,phi
#endif

   use mo_mpi
   use cputim



   use configure
! LOCAL variables that are not in modules
   implicit none


   real(kind=rk), intent(out) :: selftime(2), childrentime(2)
   real(kind=rk)              :: selftime_start(2),           &
                                 time_communication_start(2), &
                                 time_communication_end(2)

   logical, intent(in)        :: do_mpi_comm

   real(kind=rk), dimension(config%qx,config%qz,config%isma) :: sbuf, rbuf

   integer(kind=ik)           :: ierr, dest, src, mpistat(MPI_STATUS_SIZE), &
                                sendcount

   integer(kind=ik)           :: i, jm1, jp1, j, is, ndim, nstrt, k, jmi, kmi
   integer(kind=ik)           ::  iri(config%qx+1)
   real(kind=rk)              ::  xl(config%qx+1),rrq(-4:config%qx+4)


#ifdef ENHANCEQ
   real(kind=rk)              :: beta, tfad=20e-3, tzero=QENHTZER, beta0=QFACENH
   logical                    :: gain
#endif



   selftime     = 0._rk
   childrentime = 0._rk

#ifndef DEBUG_TIMINGS
   call second_v(selftime_start)
#endif

   do i=1,config%qx
      xl(i)=xzltot(i)
   enddo
   xl(config%qx+1)=xzrtot(config%qx)

   do k=nzmoms,nzmome

      if ((use_mpi) .and. (config%use_spherical_eddington_factor)) then
         kmi = qz_s

      else
         kmi = kmin(k)
      endif

      do j=nymoms,nymome

         if ((use_mpi) .and. (config%use_spherical_eddington_factor)) then

            jmi = qy_s
         else
            jmi = jmin(j)
         endif

         nstrt= nxlb(j)
         ndim = nxub(j)-nstrt+1

!- (re)construct (extended) radial TRANSP-grid 
         if (config%use_multid_collapse) then
            r (-1:config%imaxp +1)   = ralag_3d(:,j,k)
         else
            r (-1:config%imaxp +1)   = ralag_1d(:)
         endif !multid_collapse
         call himo_r(0)

         do i=0,config%imaxp +1
            rrq(i)=rq(i)
         enddo
         rrq(-1)=-rrq(0)
         rrq(config%imaxp +1+1)=2._rk*r(config%imaxp +1)-rrq(config%imaxp +1)

! -interpolation of the sourceterms
         call sort_vec(r(-1:),config%imaxp +1+2,1,xl(nstrt:),iri(nstrt:),ndim+1,1)


#ifdef CFC_TRANSPORT
!     Die interpolierten Groessen selag, usw. muessen hier bereits mit den
!     entsprechenden Metrikfaktoren multipliziert worden sein (siehe
!     sourcet und enemom in qterms.F).
#endif
         call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,                &
              xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),ndim, &
              selag(-4:,j,k),qentot(nstrt:,jmi,kmi) )
         call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,                &
              xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),ndim, &
              sylag(-4:,j,k,1),qyetot(nstrt:,jmi,kmi,1) )
         call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,                &
              xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),ndim,&
              sylag(-4:,j,k,2),qyetot(nstrt:,jmi,kmi,2) )
         ! more sourcterms
         call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,                &
              xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),ndim,&
              sylag(-4:,j,k,3),qyetot(nstrt:,jmi,kmi,3) )
         call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,                &
              xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),ndim,&
              sylag(-4:,j,k,4),qyetot(nstrt:,jmi,kmi,4) )
         call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,                &
              xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),ndim,&
              sylag(-4:,j,k,5),qyetot(nstrt:,jmi,kmi,5) )
         call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,                &
              xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),ndim,&
              smlag(-4:,j,k),qmotot(nstrt:,jmi,kmi) )
#ifdef FCNC_CALC
         call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,                &
              xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),ndim,&
              selag_n(-4,j,k),qentot_n(nstrt,jmi,kmi) )
         call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,                &
              xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),ndim,&
              sylag_n(-4,j,k),qyetot_n(nstrt,jmi,kmi) )
         call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,                &
              xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),ndim,&
              smlag_n(-4,j,k),qmotot_n(nstrt,jmi,kmi) )
#endif


! -interpolation of the stress-energy-terms (for output only)
         do is=1,config%isma
            call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,             &
                 xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),  &
                 ndim,enulag(-4:,is,j,k),enutot(nstrt:,jmi,kmi,is) )
            call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,             &
                 xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),  &
                 ndim,fnulag(-4:,is,j,k),fnutot(nstrt:,jmi,kmi,is) )
            call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,             &
                 xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),  &
                 ndim,pnulag(-4:,is,j,k),pnutot(nstrt:,jmi,kmi,is) )
            call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,             &
                 xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),  &
                 ndim,dnulag(-4:,is,j,k),dnutot(nstrt:,jmi,kmi,is) )
            call cintp (r(-1:),rrq(-4:),r(0:),config%imaxp +1+1,             &
                 xzltot(nstrt:),xzrtot(nstrt:),iri(nstrt:),  &
                 ndim,gnulag(-4:,is,j,k),gnutot(nstrt:,jmi,kmi,is) )
         enddo

      enddo ! j-loop

   enddo ! k-loop



   if ((config%nsdim .ge. 1).and.((config%nymom .ne. config%qy) .or. (config%nztra .ne. config%qz))) then

      if (config%nymom .ne. config%qy .or. config%nztra .ne. config%qz) call cintp_s
   endif

   gnutot(config%qx,qy_s,qz_s,1:config%isma)=gnulag(config%imaxp +1,1:config%isma,nymoms,nzmoms)* &
        r(config%imaxp +1)**2/xzrtot(config%qx)**2
   fnutot(config%qx,qy_s,qz_s,1:config%isma)=fnulag(config%imaxp +1,1:config%isma,nymoms,nzmoms)* &
        r(config%imaxp +1)**2/xzrtot(config%qx)**2



!-- convert to CGS-units 

   do k = qz_s,qz_e
      do j=qy_s, qy_e

#ifdef ENHANCEQ
! Enhance the neutrino heating by factor beta in gain region
         gain=.true.
         do i=config%qx,1,-1
            beta=min(beta0,1.0_rk+max(0.e9_rk,(beta0-1._rk)*(time-tzero)/tfad))    

            if ( (gain) .and. (qentot(i,j,k) .ge. 0.0_rk) ) &
                 qentot(i,j,k)=beta*qentot(i,j,k)

            if ( qentot(i,j,k)-qmotot(i,j,k)*vextot(i,j,k) .lt. 0.0_rk) &
                 gain=.false.
         enddo
#endif

         do i=1,config%qx
            qentot(i,j,k)=qentot(i,j,k) * pc_meverg
            qmotot(i,j,k)=qmotot(i,j,k) * pc_meverg
            !              qyetot(i,j,k)=qyetot(i,j,k)
#ifdef FCNC_CALC
            qentot_n(i,j,k)=qentot_n(i,j,k) * pc_meverg
            qmotot_n(i,j,k)=qmotot_n(i,j,k) * pc_meverg
#endif
            do is=1,config%isma
!                  dnutot(i,j,k,is)=dnutot(i,j,k,is)            
               enutot(i,j,k,is)=enutot(i,j,k,is) * pc_meverg
               pnutot(i,j,k,is)=pnutot(i,j,k,is) &
#ifdef CFC_TRANSPORT2
                    / sqrt_gamma(i,j,k) &
#endif
                                                    * pc_meverg
               fnutot(i,j,k,is)=fnutot(i,j,k,is) * pc_meverg
            enddo
#ifdef ORIGIN
            qmotot(1,j,k)=0.0
#endif         
         enddo
      enddo
   enddo ! k-loop 2

!-- lateral neutrino acceleration (in diffusion approximation) 
   if (config%nsdim .eq. 2) then
   if (use_mpi) then
      if (do_mpi_comm .and. nprocs .gt. 1) then
         ! send left
#if !(defined(DEBUG_TIMINGS)) && !(defined(PROGRAM_remap))
         call second_v(time_communication_start)
#endif
         sbuf(:,:,:) = pnutot(:,qy_s,:,:)
         call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr)
         sendcount = config%qx*config%qz*config%isma
         call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep1, &
                      rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
                      MPI_COMM_WORLD, mpistat, ierr)
         pnutot(:,qy_e+1,:,:) = rbuf(:,:,:)

         ! send right
         sbuf(:,:,:) = pnutot(:,qy_e,:,:)
         call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
         sendcount = config%qx*config%qz*config%isma
         call MPI_SendRecv(sbuf, sendcount, MPI_DOUBLE_PRECISION, dest, tag_sweep2, &
                      rbuf, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep2,  &
                      MPI_COMM_WORLD, mpistat, ierr)
         pnutot(:,qy_s-1,:,:) = rbuf(:,:,:) 

#if !(defined(DEBUG_TIMINGS)) && !(defined(PROGRAM_remap))
         call second_v(time_communication_end)
         timer%transp_comm = timer%transp_comm + (time_communication_end-time_communication_start)
         childrentime = childrentime + (time_communication_end-time_communication_start)
#endif

      endif
   endif ! use-mpi
endif ! config%nsdim

      if (config%nymom .gt. 1 .and. config%nztra .eq. 1) then
        k=1

         do j=qy_s,qy_e

            select case (config%latybc)
            case (1)
               jp1=j+1
               if (jp1 .eq. config%nymom+1) jp1=config%nymom
               jm1=j-1
               if (jm1 .eq. 0) jm1=1
            case (4)
               jp1=j+1
               if (jp1 .eq. config%nymom+1) jp1=1
               jm1=j-1
               if (jm1 .eq. 0) jm1=config%nymom
            case default
               write (*,*) "Task ",myproc
               write (*,*) config%latybc
               raise_abort("mpa_ra2hyd(): lateral boundary condition not implemented")
            end select
            do i=1,config%qx
               if (dentot(i,j,k).ge.1.e12_rk) then
                  qmytot(i,j,k)=-SUM( neu_qw(1:config%isma)* &
                       (pnutot(i,jp1,k,1:config%isma)-pnutot(i,jm1,k,1:config%isma))) &
                       /(xzntot(i)*(th(jp1)-th(jm1)))
#ifdef CFC_TRANSPORT2
                  qmytot(i,j,k)=qmytot(i,j,k)* &
                       alpha(i,j,k)/phi(i,j,k)**2* &
                       sqrt_gamma(i,j,k)         !component in orthonormal tetrad
#endif
               else
                  qmytot(i,j,k)=0.0_rk
                  endif
            enddo
         enddo
         qmytot(1,:,:)=0.0_rk
      else
         qmytot(:,:,:)=0.0_rk
      endif

!     write(*,*) 'CHECKGRID_1 C',SUM(qentot),SUM(qmotot),SUM(qyetot)
!     raise_abort("mpa_ra2hyd(): ra2hyd ende")

#ifndef DEBUG_TIMINGS
      call second_v(selftime)
      selftime = selftime - selftime_start
#endif      
      return


 end subroutine map_ra2hyd
#endif /* not PROGRAM_remap or NOTRA */

#if !defined(NOTRA) && !defined(PROGRAM_remap)
!> \verbatim
!>
!> Purpose: - Map(radially) and average(in theta and in phi for 3d) from  the hydro-
!>            onto the transport-grid
!>              one angular zone (j>0) for the transport may cover several
!>                angular cells of the hydro-grid "single angular average"
!>              if called with j=0 the average over all angular cells 
!>                of the hydro-grid "total angular average" is taken
!>
!>          - Conserves Mass, Energy and Momentum locally
!>            
!>          - Get radial overlap of the grids
!>
!> Input: thphmean(:,:,1): weights for the single angular averages (>0)
!>        thphmean(:,:,0): weights for the total  angular average (j=0)
!>        
!>        mode: =1 map all local thermodyn quantities
!>              =0 map only density and velx (used for call by lagrange) 
!>
!> Output: nxlb(j),nxub(j): index of bounds of the lower/upper part of 
!>                                uncovered region, ie.
!>
!>           Hydro         Transp            Hydro
!>
!>        1....nxlb  |  -1...config%imaxp +1  |   nxub...nx    
!>                        of the hydro-grid  
!>        novl(j): total number of cells the sourceterms are evaluated at
!>
!> Note : works on r-grid => call copysec before
!>
!>  Author: M. Rampp and F. Hanke
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>

 subroutine map_hyd2ra(j,k,i_gr,mode, selftime, childrentime)


   use precision
   use phycon
#ifndef NOTRA

!      use nutrio_hy  ! forcheck
   use totgrq_hy
   use lotdqua_rt
   use radial_grid_rt
   use multigrid_rt
   use overlap
   use backquants_rt
   use time_rt
   use vnew_hy , only: ishck

#endif /* NOTRA */

#ifdef CFC_TRANSPORT
   use size_cfc, ONLY: m,n,o
   use metric_cfc
   use sources_cfc
   use rshlat_hy
#endif
   use mo_mpi


!      use mesh_hy
!      use vnew_hy
   use totare_hy
   use intgrs_hy
   use cputim 

   use interpol_eul


   use configure
! LOCAL variables that are not in modules

   implicit none
   real(kind=rk), intent(out)    :: selftime(2), childrentime(2)
   real(kind=rk)                 :: selftime_start(2)

   integer(kind=ik)              :: i, jm, km, kk, ii, ir, il, i_nuc, jat
   integer(kind=ik)              ::  iri(-1:config%imaxp +1+1)
   real(kind=rk)                 ::  xznave(-3:config%qx+4), &
                                     acave(-3:config%qx+4),  &
                                     epave(-3:config%qx+4),  &
                                     gaave(-3:config%qx+4),  &
                                     goave(-3:config%qx+4),  &
                                     rhave(-3:config%qx+4),  &
                                     etave(-3:config%qx+4),  &
                                     teave(-3:config%qx+4),  &
                                     vxave(-3:config%qx+4),  &
                                     vyave(-3:config%qx+4),  &
                                     vzave(-3:config%qx+4),  &
                                     xnave(-3:config%qx+4, config%qn)

   real(kind=rk)                 :: rhave_buf(-3:config%qx+4), &
                                    vxave_buf(-3:config%qx+4), &
                                    vyave_buf(-3:config%qx+4), &
                                    vzave_buf(-3:config%qx+4), &
                                    etave_buf(-3:config%qx+4), &
                                    teave_buf(-3:config%qx+4), &
                                    acave_buf(-3:config%qx+4)

#ifdef CFC_TRANSPORT
   real(kind=rk)                 :: phave(-3:config%qx+4),     &
                                    alave(-3:config%qx+4),     &
                                    wlave(-3:config%qx+4),     &
                                    brave(-3:config%qx+4),     &
                                    btave(-3:config%qx+4),     &
                                    ptave(-3:config%qx+4),     &
                                    rlave(-3:config%qx+4)
#endif

   real(kind=rk)                 :: vxq(0:config%imaxp +1+1), &
                                    vyq(0:config%imaxp +1+1), &
                                    vzq(0:config%imaxp +1+1)
   real(kind=rk)                 :: xl(config%qx+1)
   real(kind=rk)                 :: crgtm, clftm, crgtv, clftv, err
   integer(kind=ik) ,intent (in) :: j, k, i_gr, mode

   real(kind=rk) ,parameter      :: dpi4 = 0.25_rk/pc_pi,     &
                                   !,div1=1._rk/pc_mb, &
#ifdef CFC_TRANSPORT
                                    div2 =   1.0_rk
#else
                                    div2 =   1.0_rk/pc_cl
#endif


    integer(kind=ik)             :: ierr

    real(kind=rk)                :: buf(config%qx), buf_buf(config%qx)

    real(kind=rk)                :: tim2(2), tim1(2)

    selftime     = 0._rk
    childrentime = 0._rk
#ifndef DEBUG_TIMINGS
    call second_v(selftime_start)
#endif

   jat=1
   if (j .eq. 0) jat=0

#ifdef CFC_TRANSPORT
   call init_tot_arrays(jmin(j),jmax(j),kmin(k),kmax(k),mode)
#endif

! --- average hydro-quantities for each fixed radius over the angle
   rhave(:) = 0._rk
   vxave(:) = 0._rk
   vyave(:) = 0._rk

   if (mode .ne. 0) then
      etave(:)   = 0._rk
      teave(:)   = 0._rk
      vzave(:)   = 0._rk
      acave(:)   = 0._rk
      xnave(:,:) = 0._rk
      epave(:)   = 0._rk
      gaave(:)   = 0._rk
      goave(:)   = 0._rk
   endif

#ifdef CFC_TRANSPORT
   phave(:) = 0.0_rk
   alave(:) = 0.0_rk
   wlave(:) = 0.0_rk
   brave(:) = 0.0_rk
   btave(:) = 0.0_rk
   ptave(:) = 0.0_rk
   rlave(:) = 0.0_rk
#endif

#ifndef CFC_TRANSPORT
   if ((use_mpi) .and. (config%use_spherical_eddington_factor)) then
      do km=qz_s,qz_e
         do jm=qy_s,qy_e
            do i=1,config%qx 
            rhave(i) = rhave(i) + thphmean(jm,km,jat)*dentot(i,jm,km)
            vxave(i) = vxave(i) + thphmean(jm,km,jat)*vextot(i,jm,km) &
                                                *dentot(i,jm,km)
            vyave(i) = vyave(i) + thphmean(jm,km,jat)*veytot(i,jm,km) &
                                                *dentot(i,jm,km)
         enddo
         if (mode .ne. 0) then
            do i=1,config%qx 
#ifdef CFC_TRANSPORT2
               etave(i) = etave(i) + thphmean(i,jm,k,jat)*enetot(i,jm,k) &
                    *dentot(i,jm,k)
               vzave(i) = vzave(i) + thphmean(i,jm,k,jat)*veztot(i,jm,k) &
                    *dentot(i,jm,k)
               teave(i) = teave(i) + thphmean(i,jm,k,jat)*temtot(i,jm,k)
               acave(i) = acave(i) + thphmean(i,jm,k,jat)*acxtot(i,jm,k)
#else /* CFC_TRANSPORT2 */
               etave(i) = etave(i) + thphmean(jm,km,jat)*enetot(i,jm,km) &
                    *dentot(i,jm,km)
               vzave(i) = vzave(i) + thphmean(jm,km,jat)*veztot(i,jm,km) &
                    *dentot(i,jm,km)
               teave(i) = teave(i) + thphmean(jm,km,jat)*temtot(i,jm,km)
               acave(i) = acave(i) + thphmean(jm,km,jat)*acxtot(i,jm,km)
#endif /* CFC_TRANSPORT2 */
            enddo
         endif
      enddo
   enddo
   else ! mpi not set
      do km=kmin(k),kmax(k)
         do jm=jmin(j),jmax(j)
            do i=1,config%qx 
            rhave(i) = rhave(i) + thphmean(jm,km,jat)*dentot(i,jm,km)
            vxave(i) = vxave(i) + thphmean(jm,km,jat)*vextot(i,jm,km) &
                                                *dentot(i,jm,km)
            vyave(i) = vyave(i) + thphmean(jm,km,jat)*veytot(i,jm,km) &
                                                *dentot(i,jm,km)
         enddo
         if (mode .ne. 0) then
            do i=1,config%qx 
#ifdef CFC_TRANSPORT2
               etave(i) = etave(i) + thphmean(i,jm,k,jat)*enetot(i,jm,k) &
                    *dentot(i,jm,k)
               vzave(i) = vzave(i) + thphmean(i,jm,k,jat)*veztot(i,jm,k) &
                    *dentot(i,jm,k)
               teave(i) = teave(i) + thphmean(i,jm,k,jat)*temtot(i,jm,k)
               acave(i) = acave(i) + thphmean(i,jm,k,jat)*acxtot(i,jm,k)
#else /* CFC_TRANSPORT2 */
               etave(i) = etave(i) + thphmean(jm,km,jat)*enetot(i,jm,km) &
                    *dentot(i,jm,km)
               vzave(i) = vzave(i) + thphmean(jm,km,jat)*veztot(i,jm,km) &
                    *dentot(i,jm,km)
               teave(i) = teave(i) + thphmean(jm,km,jat)*temtot(i,jm,km)
               acave(i) = acave(i) + thphmean(jm,km,jat)*acxtot(i,jm,km)
#endif /* CFC_TRANSPORT2 */
            enddo
         endif
      enddo
   enddo
   endif
#else /* CFC_TRANSPORT */
   if ((use_mpi) .and. (config%use_spherical_eddington_factor)) then
      do jm=qy_s,qy_e
         do i=1,config%qx 
            rhave(i) = rhave(i) + thphmean(i,jm,k,jat)*dentot(i,jm,k)
            vxave(i) = vxave(i) + thphmean(i,jm,k,jat)*vextot(i,jm,k) &
                 *dentot(i,jm,k)
            vyave(i) = vyave(i) + thphmean(i,jm,k,jat)*veytot(i,jm,k) &
                 *dentot(i,jm,k)
            phave(i) = phave(i) + thphmean(i,jm,k,jat)*phi(i,jm,k)
            alave(i) = alave(i) + thphmean(i,jm,k,jat)*alpha(i,jm,k)
            wlave(i) = wlave(i) + thphmean(i,jm,k,jat)*wltot(i,jm,k)
            brave(i) = brave(i) + thphmean(i,jm,k,jat)*beta_up_1(i,jm,k)
            btave(i) = btave(i) + thphmean(i,jm,k,jat)*beta_up_2(i,jm,k)
            ptave(i) = ptave(i) + thphmean(i,jm,k,jat)*beta_up_k_k(i,jm,k)
            if (config%qy .gt. 1) then
               rlave(i) = rlave(i) + thphmean(i,jm,k,jat)*rshlat(i,jm,k)
            endif
         enddo
         if (mode .ne. 0) then
            do i=1,config%qx 
#ifdef CFC_TRANSPORT2
               etave(i) = etave(i) + thphmean(i,jm,k,jat)*enetot(i,jm,k) &
                    *dentot(i,jm,k)
               vzave(i) = vzave(i) + thphmean(i,jm,k,jat)*veztot(i,jm,k) &
                    *dentot(i,jm,k)
               teave(i) = teave(i) + thphmean(i,jm,k,jat)*temtot(i,jm,k)
               acave(i) = acave(i) + thphmean(i,jm,k,jat)*acxtot(i,jm,k)
#else /* CFC_TRANSPORT2 */
               etave(i) = etave(i) + thphmean(jm,km,jat)*enetot(i,jm,km) &
                    *dentot(i,jm,km)
               vzave(i) = vzave(i) + thphmean(jm,km,jat)*veztot(i,jm,km) &
                    *dentot(i,jm,km)
               teave(i) = teave(i) + thphmean(jm,km,jat)*temtot(i,jm,km)
               acave(i) = acave(i) + thphmean(jm,km,jat)*acxtot(i,jm,km)
#endif /* CFC_TRANSPORT2 */
            enddo
         endif
      enddo
   else ! mpi not set
      do jm=jmin(j),jmax(j)

         do i=1,config%qx 

            rhave(i) = rhave(i) + thphmean(i,jm,k,jat)*dentot(i,jm,k)
            vxave(i) = vxave(i) + thphmean(i,jm,k,jat)*vextot(i,jm,k) &
                 *dentot(i,jm,k)
            vyave(i) = vyave(i) + thphmean(i,jm,k,jat)*veytot(i,jm,k) &
                 *dentot(i,jm,k)
            phave(i) = phave(i) + thphmean(i,jm,k,jat)*phi(i,jm,k)
            alave(i) = alave(i) + thphmean(i,jm,k,jat)*alpha(i,jm,k)
            wlave(i) = wlave(i) + thphmean(i,jm,k,jat)*wltot(i,jm,k)
            brave(i) = brave(i) + thphmean(i,jm,k,jat)*beta_up_1(i,jm,k)
            btave(i) = btave(i) + thphmean(i,jm,k,jat)*beta_up_2(i,jm,k)
            ptave(i) = ptave(i) + thphmean(i,jm,k,jat)*beta_up_k_k(i,jm,k)
            if (config%qy .gt. 1) then
               rlave(i) = rlave(i) + thphmean(i,jm,k,jat)*rshlat(i,jm,k)
            endif

         enddo
         if (mode .ne. 0) then
            do i=1,config%qx 
#ifdef CFC_TRANSPORT2
               etave(i) = etave(i) + thphmean(i,jm,k,jat)*enetot(i,jm,k) &
                    *dentot(i,jm,k)
               vzave(i) = vzave(i) + thphmean(i,jm,k,jat)*veztot(i,jm,k) &
                    *dentot(i,jm,k)
               teave(i) = teave(i) + thphmean(i,jm,k,jat)*temtot(i,jm,k)
               acave(i) = acave(i) + thphmean(i,jm,k,jat)*acxtot(i,jm,k)
#else /* CFC_TRANSPORT2 */
               etave(i) = etave(i) + thphmean(jm,km,jat)*enetot(i,jm,km) &
                    *dentot(i,jm,km)
               vzave(i) = vzave(i) + thphmean(jm,km,jat)*veztot(i,jm,km) &
                    *dentot(i,jm,km)
               teave(i) = teave(i) + thphmean(jm,km,jat)*temtot(i,jm,km)
               acave(i) = acave(i) + thphmean(jm,km,jat)*acxtot(i,jm,km)
#endif /* CFC_TRANSPORT2 */
            enddo
         endif
      enddo
   endif

#endif /* CFC_TRANSPORT */

   if (config%use_spherical_eddington_factor) then

   if (use_mpi) then
#if !(defined(DEBUG_TIMINGS)) && !(defined(PROGRAM_remap))
      call second_v(tim1)
#endif

      call MPI_AllReduce(rhave(1:config%qx), rhave_buf(1:config%qx), config%qx, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      rhave = rhave_buf

      call MPI_AllReduce(vxave(1:config%qx), vxave_buf(1:config%qx), config%qx, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      vxave = vxave_buf

      call MPI_AllReduce(vyave(1:config%qx), vyave_buf(1:config%qx), config%qx, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      vyave = vyave_buf

      if (mode .ne. 0) then
         call MPI_AllReduce(etave(1:config%qx), etave_buf(1:config%qx), config%qx, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, MPI_COMM_WORLD, ierr)
         etave = etave_buf
         call MPI_AllReduce(teave(1:config%qx), teave_buf(1:config%qx), config%qx, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, MPI_COMM_WORLD, ierr)
         teave = teave_buf

         call MPI_AllReduce(vzave(1:config%qx), vzave_buf(1:config%qx), config%qx, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, MPI_COMM_WORLD, ierr)
         vzave = vzave_buf

         call MPI_AllReduce(acave(1:config%qx), acave_buf(1:config%qx), config%qx, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, MPI_COMM_WORLD, ierr)
         acave = acave_buf
#if !(defined(DEBUG_TIMINGS)) && !(defined(PROGRAM_remap))
         call second_v(tim2)

         timer%transp_comm = timer%transp_comm + (tim2-tim1)
#endif
      endif
   endif ! use_mpi
endif

   if (mode .ne. 0) then
      if (i_gr .ne. 0) then
#ifndef CFC_TRANSPORT2
         do i=1,config%qx 
            acave(i) = grptot(i)
            epave(i) = 0.5_rk*(ephtot(i-1)+ephtot(i))
            gaave(i) = 0.5_rk*(gamtot(i-1)+gamtot(i))
            goave(i) = 0.5_rk*(gamold(i-1)+gamold(i))
         enddo
#else /* CFC_TRANSPORT2 */
         epave(:)=1.0_rk
         do jm = jmin(j), jmax(j)
            do i=1,config%qx
               epave(i)=epave(i)+thphmean(i,jm,k,jat)*phi_potential(i,j,k)
               acave(i)=acave(i)+thphmean(i,jm,k,jat)*phi_potential(i,j,k)
            enddo
         enddo
#endif /* CFC_TRANSPORT2 */
      endif
#ifndef CFC_TRANSPORT
      if ((use_mpi) .and. (config%use_spherical_eddington_factor)) then
         do i_nuc=1,config%qn
            do km=qz_s,qz_e
               do jm=qy_s,qy_e
                  do i=1,config%qx
#ifdef CFC_TRANSPORT2
                     xnave(i,i_nuc)=xnave(i,i_nuc)+thphmean(i,jm,k,jat)* &
                       xnutot(i,jm,k,i_nuc)*dentot(i,jm,k)
#else /* CFC_TRANSPORT2 */
                     xnave(i,i_nuc)=xnave(i,i_nuc)+thphmean(jm,km,jat)* &
                       xnutot(i,jm,km,i_nuc)*dentot(i,jm,km)
#endif /* CFC_TRANSPORT2 */
                  enddo
               enddo
            enddo
         enddo
      else ! (use_mpi) .and. (config%use_spherical_eddington_factor)
         do i_nuc=1,config%qn
            do km=kmin(k),kmax(k)
               do jm=jmin(j),jmax(j)
                  do i=1,config%qx
#ifdef CFC_TRANSPORT2
                  xnave(i,i_nuc)=xnave(i,i_nuc)+thphmean(i,jm,k,jat)* &
                       xnutot(i,jm,k,i_nuc)*dentot(i,jm,k)
#else /* CFC_TRANSPORT2 */
                  xnave(i,i_nuc)=xnave(i,i_nuc)+thphmean(jm,km,jat)* &
                       xnutot(i,jm,km,i_nuc)*dentot(i,jm,km)
#endif /* CFC_TRANSPORT2 */
               enddo
            enddo
         enddo
      enddo
   endif ! USE_MPI

#else /* CFC_TRANSPORT */
   if ((use_mpi) .and. (config%use_spherical_eddington_factor)) then
      do i_nuc=1,config%qn
         do jm=qy_s,qy_e
            do i=1,config%qx
#ifdef CFC_TRANSPORT2
               xnave(i,i_nuc)=xnave(i,i_nuc)+thphmean(i,jm,k,jat)* &
                              xnutot(i,jm,k,i_nuc)*dentot(i,jm,k)
#else /* CFC_TRANSPORT2 */
               xnave(i,i_nuc)=xnave(i,i_nuc)+thphmean(jm,km,jat)* &
                              xnutot(i,jm,km,i_nuc)*dentot(i,jm,km)
#endif /* CFC_TRANSPORT2 */
            enddo
         enddo
      enddo
   else ! ((use_mpi) .and. (config%use_spherical_eddington_factor))
      do i_nuc=1,config%qn
         do jm=jmin(j),jmax(j)

            do i=1,config%qx
#ifdef CFC_TRANSPORT2
               xnave(i,i_nuc)=xnave(i,i_nuc)+thphmean(i,jm,k,jat)* &
                              xnutot(i,jm,k,i_nuc)*dentot(i,jm,k)
#else /* CFC_TRANSPORT2 */
               xnave(i,i_nuc)=xnave(i,i_nuc)+thphmean(jm,km,jat)* &
                              xnutot(i,jm,km,i_nuc)*dentot(i,jm,km)
#endif /* CFC_TRANSPORT2 */
            enddo
         enddo
      enddo
   endif

#endif /* CFC_TRANSPORT */

endif


   if ((use_mpi) .and. (config%use_spherical_eddington_factor)) then

#if !(defined(DEBUG_TIMINGS)) && !(defined(PROGRAM_remap))
         call second_v(tim1)
#endif
         do i_nuc=1,config%qn
            buf(:) = xnave(1:config%qx,i_nuc)

            call MPI_AllReduce(buf, buf_buf, config%qx, MPI_DOUBLE_PRECISION, &
                               MPI_SUM, MPI_COMM_WORLD, ierr)
            buf(:)=buf_buf(:)
            xnave(1:config%qx,i_nuc) = buf(:)
         enddo
#ifndef DEBUG_TIMINGS
         call second_v(tim2)

         timer%transp_comm = timer%transp_comm + (tim2-tim1)
#endif
      endif ! use_mpi


      do i=1,config%qx
         xznave(i)=xzntot(i)
      enddo

#if defined(CFC_TRANSPORT) && !defined(ORIGIN)
!     Additional boundary condition to resolve
!     problems with the coordinate singularity at r=0
!     in CoCoNuT
      vxave(1) = vxave(2)*xzntot(1)/xzntot(2)
#endif

      ! boundary conditions  (bottom)
      xznave(0)=-xznave(1)
      rhave(0) = rhave(1)
      vxave(0) =-vxave(1)
      vyave(0) = vyave(1)
      vzave(0) = vzave(1)
      teave(0) = teave(1)
      etave(0) = etave(1)
      acave(0) =-acave(1)
      xnave(0,:)=xnave(1,:)
      epave(0) = epave(1)
      gaave(0) = 2._rk-gaave(1)
      goave(0) = 2._rk-goave(1)
#ifdef CFC_TRANSPORT2
      alave(0) = alave(1)
      phave(0) = phave(1)
      wlave(0) = wlave(1)
      brave(0) =-brave(1)
      btave(0) = btave(1)
      ptave(0) = ptave(1)

      if (config%nsdim .ge. 2) then
         rlave(0) = rlave(1)
      endif
#endif /* CFC_TRANSPORT2 */

! boundary conditions  (top)
      xznave(config%qx+1) = 2*xznave(config%qx)-xznave(config%qx-1)
      rhave(config%qx+1) = rhave(config%qx)
      vxave(config%qx+1) = vxave(config%qx)
      vyave(config%qx+1) = vyave(config%qx)
      vzave(config%qx+1) = vzave(config%qx)
      teave(config%qx+1) = teave(config%qx)
      etave(config%qx+1) = etave(config%qx)
      acave(config%qx+1) = acave(config%qx)
      xnave(config%qx+1,:)=xnave(config%qx,:)
      epave(config%qx+1) = epave(config%qx)
      gaave(config%qx+1) = gaave(config%qx)
      goave(config%qx+1) = goave(config%qx)
#ifdef CFC_TRANSPORT2
      alave(config%qx+1) = alave(config%qx)
      phave(config%qx+1) = phave(config%qx)
      wlave(config%qx+1) = wlave(config%qx)
      brave(config%qx+1) = brave(config%qx)
      btave(config%qx+1) = btave(config%qx)
      ptave(config%qx+1) = ptave(config%qx)
      if (config%nsdim .ge. 2) then
         rlave(config%qx+1) = rlave(config%qx)
      endif
#endif /* CFC_TRANSPORT2 */




! --- do the radial interpolation conservatively 
      xl(1:config%qx)=xzltot(1:config%qx)
      xl(config%qx+1)=xzrtot(config%qx)

      call sort_vec(xl(1:),config%qx+1,1,r(-1:),iri(-1:),config%imaxp +1+2,1)

      nxub(j)=min(iri(config%imaxp +1),config%qx+1)
      nxlb(j)=max(iri(-1),0)
      if (r(config%imaxp +1).eq.xl(nxub(j))) nxub(j)=nxub(j)-1

! ... rho
      call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                  r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                  rhave(-3),rhq(0) )

! ... v_x
      call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                  r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                  vxave(-3),vxq(0) )
! ... v_y
         call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                     r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                     vyave(-3),vyq(0) )
      if (mode .ne. 0) then
! ... v_z
         call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                     r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                     vzave(-3),vzq(0) )
! ... E_tot
         call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                     r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                     etave(-3),eiq(0) )
! ... Y_e and X_{i} (baryonic mass fractions)
         do i_nuc=1,config%qn
            call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                        r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                        xnave(-3,i_nuc),xnq(0,i_nuc) )
            do i=0,config%imaxp +1
               xnq(i,i_nuc)= xnq(i,i_nuc)/rhq(i)
            enddo
         enddo

! ... T (is required as a guess-value for eos)
         call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                     r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                     teave(-3),tmq(0) )
! ... a 
         call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                     r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                     acave(-3),acq(0) )

!         write(*,*) '--> hyd2ra energy:',
!     &        SUM( etave(1:config%qx)*
!     &        (xzrtot(1:config%qx)**3-xzltot(1:config%qx)**3) ),
!     &        SUM( eiq(0:config%imaxp +1)* 
!     &        (r(0:config%imaxp +1)**3-r(-1:config%imaxp +1-1)**3) )
!
!         write(*,*) '--> hyd2ra momentum:',
!     &        SUM( vxave(1:config%qx)*
!     &        (xzrtot(1:config%qx)**3-xzltot(1:config%qx)**3) ),
!     &        SUM( vxq(0:config%imaxp +1)*
!     &        (r(0:config%imaxp +1)**3-r(-1:config%imaxp +1-1)**3) )
!         stop

      endif
      do i=0,config%imaxp +1
         vxq(i) = vxq(i)/rhq(i)
         vyq(i) = vyq(i)/rhq(i)
         beq(i) = vxq(i) * div2 
         beyq(i) = vyq(i) * div2 
         if (mode .ne. 0) then
            vzq(i) = vzq(i)/rhq(i)
#ifndef CFC_TRANSPORT
            eiq(i) = eiq(i) - 0.5_rk*rhq(i)*(vxq(i)**2 + vyq(i)**2 + vzq(i)**2 )
            acq(i) = acq(i) * div2 * div2 
#endif
         endif
      enddo

#ifdef CFC_TRANSPORT2
!     Interpolate deviation from 1 of lapse function, conformal factor and
!     Lorentz factor, otherwise cintp gives unreasonable results at the
!     grid boundary:
      alave(:)=alave(:)-1.0_rk
      phave(:)=phave(:)-1.0_rk
      wlave(:)=wlave(:)-1.0_rk
      call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                  r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                  alave(-3),alpq(0) )
      call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                  r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                  phave(-3),phigrq(0) )
      call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                  r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                  wlave(-3),wlq(0) )
      call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                  r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                  brave(-3),bshq(0) )
      call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                  r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                  btave(-3),bshthq(0) )
      call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                  r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                  ptave(-3),dlgphidtq(0) )
#ifdef LATRS_HYD
      if (config%nsdim .ge. 2) then
         call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                  r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                  rlave(-3),latrsq(0) )
      endif
#endif /* LATRS_HYD */

      alpq(:)=alpq(:)+1.0_rk
      phigrq(:)=phigrq(:)+1.0_rk
      wlq(:)=wlq(:)+1.0_rk
      dlgphidtq(:)=dlgphidtq(:)/6.0_rk

!     Phi^6 und W auf Transport-Gitter bekannt, jetzt Metrikfaktor aus
!     rhq und eiq entfernen:
      rhq(:)=rhq(:)/phigrq(:)**6/wlq(:)
      eiq(:)=eiq(:)/phigrq(:)**6/wlq(:)

#else  /* CFC_TRANSPORT2 */
      if (mode .ne. 0) then

         if (i_gr .eq. 0) then
            do i=0,config%imaxp +1
               exq(i)    = 1._rk
               gamq(i)   = 1._rk
               dtlgaq(i) = 0._rk
            enddo
         else

! ... e^phi
            call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                        r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                        epave(-3),exq(0) )
! ... Gamma
            call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                        r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                        gaave(-3),gamq(0) )
! ... Gamma_old
            call cintp (xzltot(1:),xznave(-3:),xzrtot(1:),config%qx, &
                        r(-1:),r(0),iri(-1),  config%imaxp +1+1,        &
                        goave(-3),dtlgaq(0) )
            do i=0,config%imaxp +1
               dtlgaq(i)=(1.0_rk-dtlgaq(i)/gamq(i)) * dxihm1
            enddo
         endif
      endif
#endif /* CFC_TRANSPORT2 */

#if defined(NOMOTO_MODEL) || defined(INITIAL_HYDRO_GRID_FILE)
      do i=1,config%imaxp +1
         if (rq(i).gt.6.0e7_rk) then
            beq(i)=min(beq(i), &
                 0.2_rk*(6.0e7_rk/rq(i))**3)
            wlq(i)=1.0_rk/sqrt(1.0_rk-beq(i)**2)
         end if
      end do
#endif /* NOMOTO_MODEL || (INITIAL_HYDRO_GRID_FILE)*/


#ifdef NOVELO_TR
      beq=0._rk   
      beyq=0._rk   
      acq=0._rk   
#endif


      be(-1) = 0._rk
      be(config%imaxp +1) = beq(config%imaxp +1)
      bey(-1) = 0._rk
      bey(config%imaxp +1) = beyq(config%imaxp +1)
      if (mode .ne. 0) then
         acc(-1)= 0._rk
         ex(-1) = exq(0)
         gam(-1) = gamq(0)
         dtlga(-1) = dtlgaq(0)
         acc(config%imaxp +1) = acq(config%imaxp +1)
         ex(config%imaxp +1) = exq(config%imaxp +1)
         gam(config%imaxp +1) = gamq(config%imaxp +1)
         dtlga(config%imaxp +1) = dtlgaq(config%imaxp +1)
      endif

#ifdef CFC_TRANSPORT2
      alp  (-1)=alpq  (0)
      phigr(-1)=phigrq(0)
      wl   (-1)=wlq   (0)
      bsh  (-1)=0.0_rk
      bshth(-1)=0.0_rk
      dlgphidt(-1)=dlgphidtq(0)
      alp  (config%imaxp +1)=alpq  (config%imaxp +1)
      phigr(config%imaxp +1)=phigrq(config%imaxp +1)
      wl   (config%imaxp +1)=wlq   (config%imaxp +1)
      bsh  (config%imaxp +1)=bsh   (config%imaxp +1)
      bshth(config%imaxp +1)=bshthq(config%imaxp +1)
      dlgphidt(config%imaxp +1)=dlgphidtq(config%imaxp +1)
#ifdef LATRS_HYD
      if (config%nsdim .ge. 2) then
         latrs(-1)    =latrsq(1)
         latrs(config%imaxp +1)=latrsq(config%imaxp +1)
      endif
#endif
#endif /*CFC_TRANSPORT2*/

      do i=0,config%imaxp
         clftv = 0.5_rk * dv(i)  *dpi4  
         crgtv = 0.5_rk * dv(i+1)*dpi4
         clftm = rhq(i)  *clftv
         crgtm = rhq(i+1)*crgtv
         be   (i) = ( beq(i)*clftm + beq(i+1)*crgtm ) /   &
                    (clftm + crgtm)
         bey  (i) = ( beyq(i)*clftm + beyq(i+1)*crgtm ) / &
                    (clftm + crgtm)
         if (mode .ne. 0) then
            ex   (i) = ( exq(i)*clftv + exq(i+1)*crgtv ) / &
                       (clftv + crgtv)
            gam  (i) = ( gamq(i)*clftv + gamq(i+1)*crgtv ) / &
                       (clftv + crgtv)
            acc  (i) = ( acq(i)*clftm + acq(i+1)*crgtm ) /  &
                       (clftm + crgtm)
            dtlga(i) = ( dtlgaq(i)*clftm + dtlgaq(i+1)*crgtm ) / &
                       (clftm + crgtm)
         endif
#ifdef CFC_TRANSPORT2
         alp   (i) = (  alpq(i)*clftv +  alpq(i+1)*crgtv ) / &
                     (clftv + crgtv)
         phigr (i) = (phigrq(i)*clftv +phigrq(i+1)*crgtv ) / &
                     (clftv + crgtv)
         wl    (i) = (   wlq(i)*clftv +   wlq(i+1)*crgtv ) / &
                     (clftv + crgtv)
         bsh   (i) = (  bshq(i)*clftv +  bshq(i+1)*crgtv ) / &
                     (clftv + crgtv)
         bshth (i) = (bshthq(i)*clftv +bshthq(i+1)*crgtv ) / &
                     (clftv + crgtv)
#ifdef LATRS_HYD
         latrs (i) = (latrsq(i)*clftv +latrsq(i+1)*crgtv ) / &
                     (clftv + crgtv)
#endif
#endif /* CFC_TRANSPORT2 */

      enddo


#ifdef CFC_TRANSPORT2
      detg (:)=phigr (:)**6
      detgq(:)=phigrq(:)**6

      select case (mode)
      case (2)
!     Lasse wl,wlq,detg,detgq unveraendert und lesse alte Werte aus den lag-Feldern:
         wlalt   (:)=wllag   (:,j,k)
         wlaltq  (:)=wlqlag  (:,j,k)
         detgalt (:)=detglag (:,j,k)
         detgaltq(:)=detgqlag(:,j,k)
      case (99)
!     Initialisiere wlalt und wlaltq, etc. nach Neustart (nur zu Tesztwecken),
!     konsistente Behandlung in restrt noetig!
         wlalt   (:)    =wl    (:)
         wlaltq  (:)    =wlq   (:)
         wllag   (:,j,k)=wl    (:)
         wlqlag  (:,j,k)=wlq   (:)
         detgalt (:)    =detg  (:)
         detgaltq(:)    =detgq (:)
         detglag (:,j,k)=detg  (:)
         detgqlag(:,j,k)=detgq (:)
         philag  (:,j,k)=phigr (:)
         phiqlag (:,j,k)=phigrq(:)
         beold   (:,j,k)=be    (:)
         beqold  (:,j,k)=beq   (:)
      case default
!     Verwende Lorentz- und Metrikfaktoren aus letztem Zeitschritt (wichtig fuer advec_lat):
         wl      (:)=wllag   (:,j,k)
         wlq     (:)=wlqlag  (:,j,k)
         detg    (:)=detglag (:,j,k)
         detgq   (:)=detgqlag(:,j,k)
      end select

      dwdt    (-1:config%imaxp +1)=(wl(-1:config%imaxp +1)-wlalt(-1:config%imaxp +1))*dxihm1
      dwdtq    (0:config%imaxp +1)=(wlq(0:config%imaxp +1)-wlaltq(0:config%imaxp +1))*dxihm1
      acc      (0:config%imaxp +1)=(be (0:config%imaxp +1)-beold (0:config%imaxp +1,j,k))*dxihm1*wl (0:config%imaxp +1)
      acc      (-1)      =0.0_rk
      acq      (0:config%imaxp +1)=(beq(0:config%imaxp +1)-beqold(0:config%imaxp +1,j,k))*dxihm1*wlq(0:config%imaxp +1)

      bealt(-1:config%imaxp +1)=beold(-1:config%imaxp +1,j,k)
      beqalt(0:config%imaxp +1)=beqold(0:config%imaxp +1,j,k)
#endif /* CFC_TRANSPORT2 */


      if (mode .ne. 0) then
         novl(j)=config%imaxp +1



! -- enforce X>0, SUM(X)=1
         do i_nuc=1,config%qn
            do i = 0, novl(j)
               xnq(i,i_nuc)=max(xnq(i,i_nuc),0.0_rk)
            enddo
         enddo

         do i = 0, novl(j)
            kk=SUM(MAXLOC(xnq(i,1:config%qn-1)))
            err=SUM(xnq(i,1:config%qn-1))-xnq(i,kk)
            xnq(i,kk)=1._rk-err
         enddo



         do i = 0,config%imaxp +1
            yeqold(i)=xnq(i,config%qn)
            eiqold(i)=eiq(i)
         enddo

      endif

      if (config%jvisc .eq. 3) then
         ish_ra(:) = 0.
         do km= kmin(k), kmax(k)
            do jm = jmin(j), jmax(j)
               do i = 0, config%imaxp +1
                  il = iri(i-1)
                  ir = iri(i)
                  if (xl(ir).eq.r(i)) ir=ir-1
                  do ii=il,ir
                     if (ishck(ii,jm,km).eq.1) ish_ra(i) = 1
                  enddo
               enddo
            enddo
         enddo
      endif

#ifndef DEBUG_TIMINGS
      call second_v(selftime)
      selftime = selftime - selftime_start
#endif
      return

end subroutine map_hyd2ra
#endif /* not PROGRAM_remap or NOTRA */

! -------------------------------------------------------------------
!>
!> \verbatim
!> calc. the slopes dy for a piecewise linear representation for y(x) 
!>     monotonicity-preserving
!>     Ref.: e.g. Ruffert A&A~265, 82 (1992)
!> 
!>  Author: M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
subroutine slope(x,y,nx1,dy)
  use precision

! -- 
!
! LOCAL variables that are not in modules
  implicit none

  integer(kind=ik) ,intent(in) :: nx1
  real(kind=rk) ,intent(in)    :: x(-3:nx1+4), y(-3:nx1+4)
  real(kind=rk) ,intent(out)   :: dy(nx1)
  real(kind=rk)                :: ddr, ddl
  integer(kind=ik)             :: i

  do i=1,nx1
     ddl =( y(i)   - y(i-1) ) / (x(i)  - x(i-1))
     ddr =( y(i+1) - y(i)   ) / (x(i+1)- x(i))
     dy(i) = 2._rk*max(ddl*ddr,0._rk)/(ddl+ddr+1.e-100_rk)
  enddo

end subroutine slope
!>
!> \verbatim
!>
!> Purpose:
!>       according to a monotonic piecewise linear representation 
!>          computed by "slope" performs an interpolation for the 
!>          density yc of the conserved quantity yc*yd
!>      
!>       for a spherical radial coordinate only !>!>!>
!>
!> Input: yc: density (e.g. Y_e) of the conserved quantity (e.g. n_e)
!>       
!>
!>      Note: due to the r^2 weighting in the integral over a cell,
!>              the cell average yd is different from the
!>              value in the center of the cell rh0  
!>
!>
!>      locally conserves yd 
!>      extrapolates if necessary 
!>
!>  Author: M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
subroutine cintp(xl,xm,xr,nx1,rl,rr,ipos,nr,yd,ycout)


! LOCAL variables that are not in modules
  use precision

  implicit none
  integer(kind=ik) :: ii,ilfe,j,i,ilef,irig,iii
  real(kind=rk)    :: summ
  integer(kind=ik) ,intent(in) :: nx1,nr,ipos(nr+1)
  real(kind=rk) ,intent(in) :: xl(nx1),xm(-3:nx1+4),xr(nx1), &
       rr(nr),rl(nr),yd(-3:nx1+4)
  real(kind=rk) ,intent(out) :: ycout(nr)

  real(kind=rk) dv(nx1),rh0(nx1),scrt(nx1),dd(nx1),dmas(nx1),  &
       scrt1(nr+1),scrt4(nr+1),drvol(nr),yfr(nr+1)

  real(kind=rk) xleft(nx1+1),rleft(nr+1)

  real(kind=rk) ,parameter :: drittl = 1._rk/3._rk

  xleft(1:nx1)=xl(1:nx1)
  xleft(nx1+1)=xr(nx1)
  rleft(1:nr)=rl(1:nr)
  rleft(nr+1)=rr(nr)

  call slope(xm(-3),yd(-3),nx1,dd(1))
  do i=1,nx1
     scrt (i) =         xr(i)**4.-xl(i)**4.
     dv(i) = (xr(i)**3.-xl(i)**3.)*drittl
     rh0(i) = yd(i)-0.25_rk*dd(i)*scrt(i)/dv(i)+dd(i)*xm(i)
  enddo

  do i=1,nx1
     dmas(i)=yd(i)*dv(i)
  enddo
  do j=1,nr
     drvol(j)=drittl*(rr(j)**3-rl(j)**3)  
  enddo

  do j=1,nr+1
     iii=ipos(j)
     scrt1(j) = drittl*(rleft(j)**3-xleft(iii)**3)
     scrt4(j) = 0.25_rk * (rleft(j)**4-xleft(iii)**4)
  enddo

  do j=1,nr+1
     iii=min(ipos(j),nx1)
     if (iii .eq. 0) then
        yfr(j) = yd(1)*drittl*(rl(j)**3-xl(1)**3)
     else
        yfr(j) = scrt1(j)*( rh0(iii) - dd(iii)*xm(iii) ) + scrt4(j)*dd(iii)
     endif
  enddo


  irig=ipos(1)
  do j=1,nr
     ilef=irig
     irig=ipos(j+1)
     summ = 0._rk
     do ii=ilef,irig-1
        summ=summ+dmas(ii)
     enddo
     ycout(j) = (summ+yfr(j+1)-yfr(j))/drvol(j)
  enddo

!      write(*,'(3I4)') nx1,nr,ipos(nr+1)
!      write(*,*)
!      do j=1,nr
!         iii=ipos(j)
!         write(*,'(2I4,7f10.3)') j,iii,xl(iii),rl(j),rr(j),xr(iii),
!     &        yfr(j),yfr(j+1),scrt1(j)
!      enddo
!      stop

end subroutine cintp

#if !defined(NOTRA) && !defined(PROGRAM_remap)
subroutine cintp_s

  use mo_mpi

  use phycon

  use nutrio_hy
  use totare_hy

!  use intgrs_hy
#ifndef NOTRA
  use multigrid_rt
#endif
  use param_rt

  use precision
  use abort
#if !(defined(DEBUG_TIMINGS)) && !(defined(PROGRAM_remap))
  use cputim
#endif

  use configure
! LOCAL variables that are not in modules
  implicit none

  integer(kind=ik) :: j,jmi,jma,jmil,jmir,is,jm
  integer(kind=ik) :: k,kmi,kma,kmil,kmir,km,km1,kp1

  real(kind=rk) :: jav(1:config%nymom), &
                   cyl(1:config%qy), cyr(1:config%qy), cyn(1:config%qy)

  real(kind=rk) :: kav(1:config%nztra), &
                   zl(1:config%qz), zr(1:config%qz), zn(1:config%qz)

  real(kind=rk) :: dphi

  real(kind=rk) :: ddl(1:config%qx,1:3+5*config%isma),ddr(1:config%qx,1:3+5*config%isma), &
                   dydx(1:config%qx,nymoms:nymome,nzmoms:nzmome,1:3+5*config%isma), &
                   dzdx(1:config%qx,qy_s:qy_e,nzmoms:nzmome,1:3+5*config%isma)

  real(kind=rk) :: sbuf_y(1:config%qx,nzmoms:nzmome,1:3+5*config%isma), &
                   rbuf_y(1:config%qx,nzmoms:nzmome,1:3+5*config%isma), &
                   sbuf_z(1:config%qx,qy_s:qy_e,1:3+5*config%isma), &
                   rbuf_z(1:config%qx,qy_s:qy_e,1:3+5*config%isma), &
                   ly_gridpoint(1:config%qx,nzmoms:nzmome,1:3+5*config%isma), &
                   ry_gridpoint(1:config%qx,nzmoms:nzmome,1:3+5*config%isma), &
                   lz_gridpoint(1:config%qx,qy_s:qy_e,1:3+5*config%isma), &
                   rz_gridpoint(1:config%qx,qy_s:qy_e,1:3+5*config%isma)

  integer(kind=ik) :: ierr, dest, src, mpistat(MPI_STATUS_SIZE), sendcount

  real(kind=rk) :: tim1(2), tim2(2)

!> introduce temporary arrays
  qeninp(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qentot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
  qyeinp(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1)
  qmoinp(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
  enuinp(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma) = enutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma)
  pnuinp(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma) = pnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma)
  fnuinp(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma) = fnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma)
  dnuinp(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma) = dnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma)
  gnuinp(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma) = gnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma)

!> determine grid points in theta direction
  do j=1,config%qy
    cyl(j)=cos(yzltot(j))
    cyr(j)=cos(yzrtot(j))
    cyn(j)=0.5_rk*(cyl(j)+cyr(j))
  enddo

  do j=config%nystrt,config%nymom
      jmi=jmin(j)
      jma=jmax(j)
      jav(j)=0.5_rk*( cyl(jmi)+cyr(jma) )
  enddo

!> to calculate monotonic slopes we may need coarse grid points from
!> both adjacent domains in theta-direction
!> therefore communicate the upper / lower grid point in theta-direction
  if (use_mpi) then

#if !(defined(DEBUG_TIMINGS)) && !(defined(PROGRAM_remap))
    call second_v(tim1)
#endif

    jmil = jmin(nymome)
    do k=nzmoms,nzmome
      kmi=kmin(k)
      sbuf_y(1:config%qx,k,1) = qeninp(1:config%qx,jmil,kmi)
      sbuf_y(1:config%qx,k,2) = qyeinp(1:config%qx,jmil,kmi)
      sbuf_y(1:config%qx,k,3) = qmoinp(1:config%qx,jmil,kmi)
      do is=1,config%isma
        sbuf_y(1:config%qx,k,3+is) = enuinp(1:config%qx,jmil,kmi,is)
        sbuf_y(1:config%qx,k,3+config%isma+is) = pnuinp(1:config%qx,jmil,kmi,is)
        sbuf_y(1:config%qx,k,3+2*config%isma+is) = fnuinp(1:config%qx,jmil,kmi,is)
        sbuf_y(1:config%qx,k,3+3*config%isma+is) = dnuinp(1:config%qx,jmil,kmi,is)
        sbuf_y(1:config%qx,k,3+4*config%isma+is) = gnuinp(1:config%qx,jmil,kmi,is)
      enddo
    enddo
    call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
    sendcount = config%qx*nzmom_proc*(3+5*config%isma)
    call MPI_SendRecv(sbuf_y, sendcount, MPI_DOUBLE_PRECISION, &
          dest, tag_sweep1, &
          rbuf_y, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
          MPI_COMM_WORLD, mpistat, ierr)
    ly_gridpoint(1:config%qx,nzmoms:nzmome,1) = rbuf_y(1:config%qx,nzmoms:nzmome,1)
    ly_gridpoint(1:config%qx,nzmoms:nzmome,2) = rbuf_y(1:config%qx,nzmoms:nzmome,2)
    ly_gridpoint(1:config%qx,nzmoms:nzmome,3) = rbuf_y(1:config%qx,nzmoms:nzmome,3)
    do j=0,4
      do is=1,config%isma
        ly_gridpoint(1:config%qx,nzmoms:nzmome,3+is+j*config%isma) = &
                rbuf_y(1:config%qx,nzmoms:nzmome,3+is+j*config%isma)
      enddo
    enddo

    jmir = jmin(nymoms)

    do k=nzmoms,nzmome
      kmi=kmin(k)
      sbuf_y(1:config%qx,k,1) = qeninp(1:config%qx,jmir,kmi)
      sbuf_y(1:config%qx,k,2) = qyeinp(1:config%qx,jmir,kmi)
      sbuf_y(1:config%qx,k,3) = qmoinp(1:config%qx,jmir,kmi)
      do is=1,config%isma
        sbuf_y(1:config%qx,k,3+is) = enuinp(1:config%qx,jmir,kmi,is)
        sbuf_y(1:config%qx,k,3+config%isma+is) = pnuinp(1:config%qx,jmir,kmi,is)
        sbuf_y(1:config%qx,k,3+2*config%isma+is) = fnuinp(1:config%qx,jmir,kmi,is)
        sbuf_y(1:config%qx,k,3+3*config%isma+is) = dnuinp(1:config%qx,jmir,kmi,is)
        sbuf_y(1:config%qx,k,3+4*config%isma+is) = gnuinp(1:config%qx,jmir,kmi,is)
      enddo
    enddo
    call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr)
    sendcount = config%qx*nzmom_proc*(3+5*config%isma)
    call MPI_SendRecv(sbuf_y, sendcount, MPI_DOUBLE_PRECISION, &
          dest, tag_sweep1, &
          rbuf_y, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
          MPI_COMM_WORLD, mpistat, ierr)
    ry_gridpoint(1:config%qx,nzmoms:nzmome,1) = rbuf_y(1:config%qx,nzmoms:nzmome,1)
    ry_gridpoint(1:config%qx,nzmoms:nzmome,2) = rbuf_y(1:config%qx,nzmoms:nzmome,2)
    ry_gridpoint(1:config%qx,nzmoms:nzmome,3) = rbuf_y(1:config%qx,nzmoms:nzmome,3)
    do j=0,4
      do is=1,config%isma
        ry_gridpoint(1:config%qx,nzmoms:nzmome,3+is+j*config%isma) = &
                rbuf_y(1:config%qx,nzmoms:nzmome,3+is+j*config%isma)
      enddo
    enddo
#if !(defined(DEBUG_TIMINGS)) && !(defined(PROGRAM_remap))
    call second_v(tim2)

    timer%hydro_comm = timer%hydro_comm + (tim2-tim1)
#endif
  endif

!> at first calculate the slopes dydx in theta direction
!> for coarse grid points in theta and phi direction
  do k=nzmoms,nzmome
    kmi= kmin(k)
    !> reflecting boundary in y-direction
    if (config%latybc .eq. 1) then
      do j=nymoms,nymome
        if (j.ne.1 .and. j.ne.config%nymom) then
          jmil=jmin(j-1)
          jmi= jmin(j)
          jmir=jmin(j+1)

          if (j.eq.nymoms .and. use_mpi) then
            ddl(1:config%qx,1) =( qeninp(1:config%qx,jmi,kmi)  -  ly_gridpoint(1:config%qx,k,1) ) / (jav(j)  - jav(j-1))
            ddl(1:config%qx,2) =( qyeinp(1:config%qx,jmi,kmi)  - ly_gridpoint(1:config%qx,k,2) ) / (jav(j)  - jav(j-1))
            ddl(1:config%qx,3) =( qmoinp(1:config%qx,jmi,kmi)  - ly_gridpoint(1:config%qx,k,3) ) / (jav(j)  - jav(j-1))
            do is=1,config%isma
              ddl(1:config%qx,3+is)        = &
                ( enuinp(1:config%qx,jmi,kmi,is)  - ly_gridpoint(1:config%qx,k,3+is) ) / (jav(j)  - jav(j-1))
              ddl(1:config%qx,3+config%isma+is)   = &
                ( pnuinp(1:config%qx,jmi,kmi,is)  - ly_gridpoint(1:config%qx,k,3+config%isma+is) ) / (jav(j)  - jav(j-1))
              ddl(1:config%qx,3+2*config%isma+is) = &
                ( fnuinp(1:config%qx,jmi,kmi,is)  - ly_gridpoint(1:config%qx,k,3+2*config%isma+is) ) / (jav(j)  - jav(j-1))
              ddl(1:config%qx,3+3*config%isma+is) = &
                ( dnuinp(1:config%qx,jmi,kmi,is)  - ly_gridpoint(1:config%qx,k,3+3*config%isma+is) ) / (jav(j)  - jav(j-1))
              ddl(1:config%qx,3+4*config%isma+is) = &
                ( gnuinp(1:config%qx,jmi,kmi,is)  - ly_gridpoint(1:config%qx,k,3+3*config%isma+is) ) / (jav(j)  - jav(j-1))
            enddo
          else
            ddl(1:config%qx,1) =( qeninp(1:config%qx,jmi,kmi)  - qeninp(1:config%qx,jmil,kmi) ) / (jav(j)  - jav(j-1))
            ddl(1:config%qx,2) =( qyeinp(1:config%qx,jmi,kmi)  - qyeinp(1:config%qx,jmil,kmi) ) / (jav(j)  - jav(j-1))
            ddl(1:config%qx,3) =( qmoinp(1:config%qx,jmi,kmi)  - qmoinp(1:config%qx,jmil,kmi) ) / (jav(j)  - jav(j-1))
            do is=1,config%isma
              ddl(1:config%qx,3+is)        = & 
                ( enuinp(1:config%qx,jmi,kmi,is)  - enuinp(1:config%qx,jmil,kmi,is) ) / (jav(j)  - jav(j-1))
              ddl(1:config%qx,3+config%isma+is)   = &
                ( pnuinp(1:config%qx,jmi,kmi,is)  - pnuinp(1:config%qx,jmil,kmi,is) ) / (jav(j)  - jav(j-1))
              ddl(1:config%qx,3+2*config%isma+is) = &
                ( fnuinp(1:config%qx,jmi,kmi,is)  - fnuinp(1:config%qx,jmil,kmi,is) ) / (jav(j)  - jav(j-1))
              ddl(1:config%qx,3+3*config%isma+is) = &
                ( dnuinp(1:config%qx,jmi,kmi,is)  - dnuinp(1:config%qx,jmil,kmi,is) ) / (jav(j)  - jav(j-1))
              ddl(1:config%qx,3+4*config%isma+is) = &
                ( gnuinp(1:config%qx,jmi,kmi,is)  - gnuinp(1:config%qx,jmil,kmi,is) ) / (jav(j)  - jav(j-1))
            enddo
          endif

          if (j.eq.nymome .and. use_mpi) then
            ddr(1:config%qx,1) =( ry_gridpoint(1:config%qx,k,1)  - qeninp(1:config%qx,jmi,kmi) ) / (jav(j+1)  - jav(j))
            ddr(1:config%qx,2) =( ry_gridpoint(1:config%qx,k,2)  - qyeinp(1:config%qx,jmi,kmi) ) / (jav(j+1)  - jav(j))
            ddr(1:config%qx,3) =( ry_gridpoint(1:config%qx,k,3)  - qmoinp(1:config%qx,jmi,kmi) ) / (jav(j+1)  - jav(j))
            do is=1,config%isma
              ddr(1:config%qx,3+is)        = &
                ( ry_gridpoint(1:config%qx,k,3+is)  - enuinp(1:config%qx,jmi,kmi,is) ) / (jav(j+1)  - jav(j))
              ddr(1:config%qx,3+config%isma+is)   = &
                ( ry_gridpoint(1:config%qx,k,3+config%isma+is)  - pnuinp(1:config%qx,jmi,kmi,is) ) / (jav(j+1)  - jav(j))
              ddr(1:config%qx,3+2*config%isma+is) = &
                ( ry_gridpoint(1:config%qx,k,3+2*config%isma+is)  - fnuinp(1:config%qx,jmi,kmi,is) ) / (jav(j+1)  - jav(j))
              ddr(1:config%qx,3+3*config%isma+is) = &
                ( ry_gridpoint(1:config%qx,k,3+3*config%isma+is)  - dnuinp(1:config%qx,jmi,kmi,is) ) / (jav(j+1)  - jav(j))
              ddr(1:config%qx,3+4*config%isma+is) = &
                ( ry_gridpoint(1:config%qx,k,3+4*config%isma+is)  - gnuinp(1:config%qx,jmi,kmi,is) ) / (jav(j+1)  - jav(j))
            enddo
          else
            ddr(1:config%qx,1) =( qeninp(1:config%qx,jmir,kmi)  - qeninp(1:config%qx,jmi,kmi) ) / (jav(j+1)  - jav(j))
            ddr(1:config%qx,2) =( qyeinp(1:config%qx,jmir,kmi)  - qyeinp(1:config%qx,jmi,kmi) ) / (jav(j+1)  - jav(j))
            ddr(1:config%qx,3) =( qmoinp(1:config%qx,jmir,kmi)  - qmoinp(1:config%qx,jmi,kmi) ) / (jav(j+1)  - jav(j))
            do is=1,config%isma
              ddr(1:config%qx,3+is)        = &
                ( enuinp(1:config%qx,jmir,kmi,is)  - enuinp(1:config%qx,jmi,kmi,is) ) / (jav(j+1)  - jav(j))
              ddr(1:config%qx,3+config%isma+is)   = & 
                ( pnuinp(1:config%qx,jmir,kmi,is)  - pnuinp(1:config%qx,jmi,kmi,is) ) / (jav(j+1)  - jav(j))
              ddr(1:config%qx,3+2*config%isma+is) = &
                ( fnuinp(1:config%qx,jmir,kmi,is)  - fnuinp(1:config%qx,jmi,kmi,is) ) / (jav(j+1)  - jav(j))
              ddr(1:config%qx,3+3*config%isma+is) = &
                ( dnuinp(1:config%qx,jmir,kmi,is)  - dnuinp(1:config%qx,jmi,kmi,is) ) / (jav(j+1)  - jav(j))
              ddr(1:config%qx,3+4*config%isma+is) = &
                ( gnuinp(1:config%qx,jmir,kmi,is)  - gnuinp(1:config%qx,jmi,kmi,is) ) / (jav(j+1)  - jav(j))
            enddo
          endif

          dydx(1:config%qx,j,k,1) = &
                2._rk*max(ddl(1:config%qx,1)*ddr(1:config%qx,1),0._rk)/( ddl(1:config%qx,1)+ddr(1:config%qx,1)+1.e-100_rk)
          dydx(1:config%qx,j,k,2) = &
                2._rk*max(ddl(1:config%qx,2)*ddr(1:config%qx,2),0._rk)/( ddl(1:config%qx,2)+ddr(1:config%qx,2)+1.e-100_rk)
          dydx(1:config%qx,j,k,3) = &
                2._rk*max(ddl(1:config%qx,3)*ddr(1:config%qx,3),0._rk)/( ddl(1:config%qx,3)+ddr(1:config%qx,3)+1.e-100_rk)
          do is=1,config%isma
            dydx(1:config%qx,j,k,3+is) = 2._rk*max(ddl(1:config%qx,3+is)*ddr(1:config%qx,3+is),0._rk)/ &
              ( ddl(1:config%qx,3+is)+ddr(1:config%qx,3+is)+1.e-100_rk)
            dydx(1:config%qx,j,k,3+config%isma+is) = & 
              2._rk*max(ddl(1:config%qx,3+config%isma+is)*ddr(1:config%qx,3+config%isma+is),0._rk)/ &
              ( ddl(1:config%qx,3+config%isma+is)+ddr(1:config%qx,3+config%isma+is)+1.e-100_rk)
            dydx(1:config%qx,j,k,3+2*config%isma+is) = & 
              2._rk*max(ddl(1:config%qx,3+2*config%isma+is)*ddr(1:config%qx,3+2*config%isma+is),0._rk)/ &
              ( ddl(1:config%qx,3+2*config%isma+is)+ddr(1:config%qx,3+2*config%isma+is)+1.e-100_rk)
            dydx(1:config%qx,j,k,3+3*config%isma+is) = &
              2._rk*max(ddl(1:config%qx,3+3*config%isma+is)*ddr(1:config%qx,3+3*config%isma+is),0._rk)/ &
              ( ddl(1:config%qx,3+3*config%isma+is)+ddr(1:config%qx,3+3*config%isma+is)+1.e-100_rk)
            dydx(1:config%qx,j,k,3+4*config%isma+is) = &
              2._rk*max(ddl(1:config%qx,3+4*config%isma+is)*ddr(1:config%qx,3+4*config%isma+is),0._rk)/ &
              ( ddl(1:config%qx,3+4*config%isma+is)+ddr(1:config%qx,3+4*config%isma+is)+1.e-100_rk)
          enddo
        endif
      enddo ! y-loop

      if (nymoms .eq. 1) dydx(1:config%qx,1,k,1:3+5*config%isma)  = 0._rk
      if (nymome .eq. config%nymom) dydx(1:config%qx,config%nymom,k,1:3+5*config%isma) = 0._rk

    else
      raise_abort("cintp_s(): latybc not equal 1")
    endif
  enddo ! k-loop
!> now interpolate at every coarse k-gridpoint in y(theta)-direction
  do k=nzmoms,nzmome
    kmi=kmin(k)
    do j=nymoms,nymome
      jmi=jmin(j)
      jma=jmax(j)
      do jm=jmi,jma
        qentot(1:config%qx,jm,kmi) = dydx(1:config%qx,j,k,1)*(cyn(jm)-jav(j))+qeninp(1:config%qx,jmi,kmi)
        qyetot(1:config%qx,jm,kmi,1) = dydx(1:config%qx,j,k,2)*(cyn(jm)-jav(j))+qyeinp(1:config%qx,jmi,kmi)
        qmotot(1:config%qx,jm,kmi) = dydx(1:config%qx,j,k,3)*(cyn(jm)-jav(j))+qmoinp(1:config%qx,jmi,kmi)
        do is=1,config%isma
          enutot(1:config%qx,jm,kmi,is) = &
                dydx(1:config%qx,j,k,3+is)*(cyn(jm)-jav(j))+enuinp(1:config%qx,jmi,kmi,is)
          pnutot(1:config%qx,jm,kmi,is) = &
                dydx(1:config%qx,j,k,3+config%isma+is)*(cyn(jm)-jav(j))+pnuinp(1:config%qx,jmi,kmi,is)
          fnutot(1:config%qx,jm,kmi,is) = &
                dydx(1:config%qx,j,k,3+2*config%isma+is)*(cyn(jm)-jav(j))+fnuinp(1:config%qx,jmi,kmi,is)
          dnutot(1:config%qx,jm,kmi,is) = &
                dydx(1:config%qx,j,k,3+3*config%isma+is)*(cyn(jm)-jav(j))+dnuinp(1:config%qx,jmi,kmi,is)
          gnutot(1:config%qx,jm,kmi,is) = &
                dydx(1:config%qx,j,k,3+4*config%isma+is)*(cyn(jm)-jav(j))+gnuinp(1:config%qx,jmi,kmi,is)
        enddo
      enddo
    enddo ! y-loop
  enddo ! z-loop

!> the following is only necassary for the 3d case
  if (config%nsdim .eq. 3) then
!> switch again to the temporary arrays
  qeninp(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qentot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
  qyeinp(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1)
  qmoinp(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
  enuinp(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma) = enutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma)
  pnuinp(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma) = pnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma)
  fnuinp(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma) = fnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma)
  dnuinp(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma) = dnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma)
  gnuinp(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma) = gnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%isma)

!> in phi direction the grid is equidistant
  dphi = 2._rk*pc_pi/config%nztra

!> determine grid points in phi direction
  do k=1,config%qz
    zl(k) = zzltot(k)
    zr(k) = zzrtot(k)
    zn(k)=0.5_rk*(zl(k)+zr(k))
  enddo

  do k=1,config%nztra
    kmi=kmin(k)
    kma=kmax(k)
    kav(k)=0.5_rk*( zl(kmi)+zr(kma) )
  enddo

!> to calculate monotonic slopes we may need coarse grid points from
!> both adjacent domains in phi-direction
!> therefore communicate the upper / lower grid point in phi-direction
!> for every fine theta-gridpoint
  if (use_mpi) then
#if !(defined(DEBUG_TIMINGS)) && !(defined(PROGRAM_remap))
    call second_v(tim1)
#endif
    kmil = kmin(nzmome)

    do jm=qy_s,qy_e
      sbuf_z(1:config%qx,jm,1) = qeninp(1:config%qx,jm,kmil)
      sbuf_z(1:config%qx,jm,2) = qyeinp(1:config%qx,jm,kmil)
      sbuf_z(1:config%qx,jm,3) = qmoinp(1:config%qx,jm,kmil)
      do is=1,config%isma
        sbuf_z(1:config%qx,jm,3+is) = enuinp(1:config%qx,jm,kmil,is)
        sbuf_z(1:config%qx,jm,3+config%isma+is) = pnuinp(1:config%qx,jm,kmil,is)
        sbuf_z(1:config%qx,jm,3+2*config%isma+is) = fnuinp(1:config%qx,jm,kmil,is)
        sbuf_z(1:config%qx,jm,3+3*config%isma+is) = dnuinp(1:config%qx,jm,kmil,is)
        sbuf_z(1:config%qx,jm,3+4*config%isma+is) = gnuinp(1:config%qx,jm,kmil,is)
      enddo
    enddo
    call MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
    sendcount = config%qx*qy_proc*(3+5*config%isma)
    call MPI_SendRecv(sbuf_z, sendcount, MPI_DOUBLE_PRECISION, &
          dest, tag_sweep1, &
          rbuf_z, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
          MPI_COMM_WORLD, mpistat, ierr)
    lz_gridpoint(1:config%qx,qy_s:qy_e,1) = rbuf_z(1:config%qx,qy_s:qy_e,1)
    lz_gridpoint(1:config%qx,qy_s:qy_e,2) = rbuf_z(1:config%qx,qy_s:qy_e,2)
    lz_gridpoint(1:config%qx,qy_s:qy_e,3) = rbuf_z(1:config%qx,qy_s:qy_e,3)
    do j=0,4
      do is=1,config%isma
        lz_gridpoint(1:config%qx,qy_s:qy_e,3+is+j*config%isma) = rbuf_z(1:config%qx,qy_s:qy_e,3+is+j*config%isma)
      enddo
    enddo

    kmir = kmin(nzmoms)

    do jm=qy_s,qy_e
      sbuf_z(1:config%qx,jm,1) = qeninp(1:config%qx,jm,kmir)
      sbuf_z(1:config%qx,jm,2) = qyeinp(1:config%qx,jm,kmir)
      sbuf_z(1:config%qx,jm,3) = qmoinp(1:config%qx,jm,kmir)
      do is=1,config%isma
        sbuf_z(1:config%qx,jm,3+is) = enuinp(1:config%qx,jm,kmir,is)
        sbuf_z(1:config%qx,jm,3+config%isma+is) = pnuinp(1:config%qx,jm,kmir,is)
        sbuf_z(1:config%qx,jm,3+2*config%isma+is) = fnuinp(1:config%qx,jm,kmir,is)
        sbuf_z(1:config%qx,jm,3+3*config%isma+is) = dnuinp(1:config%qx,jm,kmir,is)
        sbuf_z(1:config%qx,jm,3+4*config%isma+is) = gnuinp(1:config%qx,jm,kmir,is)
      enddo
    enddo
    call MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
    sendcount = config%qx*qy_proc*(3+5*config%isma)
    call MPI_SendRecv(sbuf_z, sendcount, MPI_DOUBLE_PRECISION, &
          dest, tag_sweep1, &
          rbuf_z, sendcount, MPI_DOUBLE_PRECISION, src, tag_sweep1,  &
          MPI_COMM_WORLD, mpistat, ierr)
    rz_gridpoint(1:config%qx,qy_s:qy_e,1) = rbuf_z(1:config%qx,qy_s:qy_e,1)
    rz_gridpoint(1:config%qx,qy_s:qy_e,2) = rbuf_z(1:config%qx,qy_s:qy_e,2)
    rz_gridpoint(1:config%qx,qy_s:qy_e,3) = rbuf_z(1:config%qx,qy_s:qy_e,3)
    do j=0,4
      do is=1,config%isma
        rz_gridpoint(1:config%qx,qy_s:qy_e,3+is+j*config%isma) = rbuf_z(1:config%qx,qy_s:qy_e,3+is+j*config%isma)
      enddo
    enddo
#if !(defined(DEBUG_TIMINGS)) && !(defined(PROGRAM_remap))
    call second_v(tim2)

    timer%hydro_comm = timer%hydro_comm + (tim2-tim1)
#endif

 endif

!> now calculate the slopes dzdx in phi direction
!> for coarse grid points in phi and for fine grid points in theta direction
  do jm=qy_s,qy_e
    !> periodic boundary in z-direction
    if (config%latzbc .eq. 4) then
      do k=nzmoms,nzmome
        km1=k-1
        kp1=k+1
        if(k .eq. 1) km1=config%nztra
        if(k .eq. config%nztra) kp1=1
        kmil=kmin(km1)
        kmi= kmin(k)
        kmir=kmin(kp1)

        if (k.eq.nzmoms .and. use_mpi) then
          ddl(1:config%qx,1) =( qeninp(1:config%qx,jm,kmi)  -  lz_gridpoint(1:config%qx,jm,1) ) / dphi
          ddl(1:config%qx,2) =( qyeinp(1:config%qx,jm,kmi)  - lz_gridpoint(1:config%qx,jm,2) ) / dphi
          ddl(1:config%qx,3) =( qmoinp(1:config%qx,jm,kmi)  - lz_gridpoint(1:config%qx,jm,3) ) / dphi
          do is=1,config%isma
            ddl(1:config%qx,3+is)        = &
                ( enuinp(1:config%qx,jm,kmi,is)  - lz_gridpoint(1:config%qx,jm,3+is) ) / dphi
            ddl(1:config%qx,3+config%isma+is)   = &
                ( pnuinp(1:config%qx,jm,kmi,is)  - lz_gridpoint(1:config%qx,jm,3+config%isma+is) ) / dphi
            ddl(1:config%qx,3+2*config%isma+is) = &
                ( fnuinp(1:config%qx,jm,kmi,is)  - lz_gridpoint(1:config%qx,jm,3+2*config%isma+is) ) / dphi
            ddl(1:config%qx,3+3*config%isma+is) = &
                ( dnuinp(1:config%qx,jm,kmi,is)  - lz_gridpoint(1:config%qx,jm,3+3*config%isma+is) ) / dphi
            ddl(1:config%qx,3+4*config%isma+is) = &
                ( gnuinp(1:config%qx,jm,kmi,is)  - lz_gridpoint(1:config%qx,jm,3+4*config%isma+is) ) / dphi
          enddo
        else
          ddl(1:config%qx,1) =( qeninp(1:config%qx,jm,kmi)  - qeninp(1:config%qx,jm,kmil)  ) / dphi
          ddl(1:config%qx,2) =( qyeinp(1:config%qx,jm,kmi)  - qyeinp(1:config%qx,jm,kmil) ) / dphi
          ddl(1:config%qx,3) =( qmoinp(1:config%qx,jm,kmi)  - qmoinp(1:config%qx,jm,kmil) ) / dphi
          do is=1,config%isma
            ddl(1:config%qx,3+is)        =( enuinp(1:config%qx,jm,kmi,is)  - enuinp(1:config%qx,jm,kmil,is) ) / dphi
            ddl(1:config%qx,3+config%isma+is)   =( pnuinp(1:config%qx,jm,kmi,is)  - pnuinp(1:config%qx,jm,kmil,is) ) / dphi
            ddl(1:config%qx,3+2*config%isma+is) =( fnuinp(1:config%qx,jm,kmi,is)  - fnuinp(1:config%qx,jm,kmil,is) ) / dphi
            ddl(1:config%qx,3+3*config%isma+is) =( dnuinp(1:config%qx,jm,kmi,is)  - dnuinp(1:config%qx,jm,kmil,is) ) / dphi
            ddl(1:config%qx,3+4*config%isma+is) =( gnuinp(1:config%qx,jm,kmi,is)  - gnuinp(1:config%qx,jm,kmil,is) ) / dphi
          enddo
        endif

        if (k.eq.nzmome .and. use_mpi) then
          ddr(1:config%qx,1) =( rz_gridpoint(1:config%qx,jm,1) - qeninp(1:config%qx,jm,kmi) ) / dphi
          ddr(1:config%qx,2) =( rz_gridpoint(1:config%qx,jm,2) - qyeinp(1:config%qx,jm,kmi) ) / dphi
          ddr(1:config%qx,3) =( rz_gridpoint(1:config%qx,jm,3) - qmoinp(1:config%qx,jm,kmi) ) / dphi
          do is=1,config%isma
            ddr(1:config%qx,3+is)        = &
                ( rz_gridpoint(1:config%qx,jm,3+is) - enuinp(1:config%qx,jm,kmi,is) ) / dphi
            ddr(1:config%qx,3+config%isma+is)   = &
                ( rz_gridpoint(1:config%qx,jm,3+config%isma+is) - pnuinp(1:config%qx,jm,kmi,is) ) / dphi
            ddr(1:config%qx,3+2*config%isma+is) = &
                ( rz_gridpoint(1:config%qx,jm,3+2*config%isma+is) - fnuinp(1:config%qx,jm,kmi,is) ) / dphi
            ddr(1:config%qx,3+3*config%isma+is) = &
                ( rz_gridpoint(1:config%qx,jm,3+3*config%isma+is) - dnuinp(1:config%qx,jm,kmi,is) ) / dphi
            ddr(1:config%qx,3+4*config%isma+is) =( rz_gridpoint(1:config%qx,jm,3+4*config%isma+is) - gnuinp(1:config%qx,jm,kmi,is) ) / dphi
          enddo
        else
          ddr(1:config%qx,1) =( qeninp(1:config%qx,jm,kmir) - qeninp(1:config%qx,jm,kmi) ) / dphi
          ddr(1:config%qx,2) =( qyeinp(1:config%qx,jm,kmir) - qyeinp(1:config%qx,jm,kmi) ) / dphi
          ddr(1:config%qx,3) =( qmoinp(1:config%qx,jm,kmir) - qmoinp(1:config%qx,jm,kmi) ) / dphi
          do is=1,config%isma
            ddr(1:config%qx,3+is)        = &
                ( enuinp(1:config%qx,jm,kmir,is) - enuinp(1:config%qx,jm,kmi,is) ) / dphi
            ddr(1:config%qx,3+config%isma+is)   = &
                ( pnuinp(1:config%qx,jm,kmir,is) - pnuinp(1:config%qx,jm,kmi,is) ) / dphi
            ddr(1:config%qx,3+2*config%isma+is) = &
                ( fnuinp(1:config%qx,jm,kmir,is) - fnuinp(1:config%qx,jm,kmi,is) ) / dphi
            ddr(1:config%qx,3+3*config%isma+is) = &
                ( dnuinp(1:config%qx,jm,kmir,is) - dnuinp(1:config%qx,jm,kmi,is) ) / dphi
            ddr(1:config%qx,3+4*config%isma+is) =( gnuinp(1:config%qx,jm,kmir,is) - gnuinp(1:config%qx,jm,kmi,is) ) / dphi
          enddo
        endif

        dzdx(1:config%qx,jm,k,1) = 2._rk*max(ddl(1:config%qx,1)*ddr(1:config%qx,1),0._rk)/ &
                ( ddl(1:config%qx,1)+ddr(1:config%qx,1)+1.e-100_rk)
        dzdx(1:config%qx,jm,k,2) = 2._rk*max(ddl(1:config%qx,2)*ddr(1:config%qx,2),0._rk)/ &
                ( ddl(1:config%qx,2)+ddr(1:config%qx,2)+1.e-100_rk)
        dzdx(1:config%qx,jm,k,3) = 2._rk*max(ddl(1:config%qx,3)*ddr(1:config%qx,3),0._rk)/ &
                ( ddl(1:config%qx,3)+ddr(1:config%qx,3)+1.e-100_rk)
        do is=1,config%isma
          dzdx(1:config%qx,jm,k,3+is) = 2._rk*max(ddl(1:config%qx,3+is)*ddr(1:config%qx,3+is),0._rk)/ &
            ( ddl(1:config%qx,3+is)+ddr(1:config%qx,3+is)+1.e-100_rk)
          dzdx(1:config%qx,jm,k,3+config%isma+is) = &
            2._rk*max(ddl(1:config%qx,3+config%isma+is)*ddr(1:config%qx,3+config%isma+is),0._rk)/ &
            ( ddl(1:config%qx,3+config%isma+is)+ddr(1:config%qx,3+config%isma+is)+1.e-100_rk)
          dzdx(1:config%qx,jm,k,3+2*config%isma+is) = &
            2._rk*max(ddl(1:config%qx,3+2*config%isma+is)*ddr(1:config%qx,3+2*config%isma+is),0._rk)/ &
            ( ddl(1:config%qx,3+2*config%isma+is)+ddr(1:config%qx,3+2*config%isma+is)+1.e-100_rk)
          dzdx(1:config%qx,jm,k,3+3*config%isma+is) = &
            2._rk*max(ddl(1:config%qx,3+3*config%isma+is)*ddr(1:config%qx,3+3*config%isma+is),0._rk)/ &
            ( ddl(1:config%qx,3+3*config%isma+is)+ddr(1:config%qx,3+3*config%isma+is)+1.e-100_rk)
          dzdx(1:config%qx,jm,k,3+4*config%isma+is) = &
            2._rk*max(ddl(1:config%qx,3+4*config%isma+is)*ddr(1:config%qx,3+4*config%isma+is),0._rk)/ &
            ( ddl(1:config%qx,3+4*config%isma+is)+ddr(1:config%qx,3+4*config%isma+is)+1.e-100_rk)
        enddo
      enddo ! k-loop
    else
      raise_abort("cintp_s(): latzbc not equal 4")
    endif
  enddo ! jm-loop

!> finally interpolate at every fine y-gridpoint in z(phi)-direction
  do jm=qy_s,qy_e
    do k=nzmoms,nzmome
      kmi=kmin(k)
      kma=kmax(k)
      do km=kmi,kma
        qentot(1:config%qx,jm,km) = &
                dzdx(1:config%qx,jm,k,1)*(zn(km)-kav(k)) + qeninp(1:config%qx,jm,kmi)
        qyetot(1:config%qx,jm,km,1) = &
                dzdx(1:config%qx,jm,k,2)*(zn(km)-kav(k))+qyeinp(1:config%qx,jm,kmi)
        qmotot(1:config%qx,jm,km) = &
                dzdx(1:config%qx,jm,k,3)*(zn(km)-kav(k))+qmoinp(1:config%qx,jm,kmi)
        do is=1,config%isma
          enutot(1:config%qx,jm,km,is) = &
                dzdx(1:config%qx,jm,k,3+is)*(zn(km)-kav(k))+enuinp(1:config%qx,jm,kmi,is)
          pnutot(1:config%qx,jm,km,is) = &
                dzdx(1:config%qx,jm,k,3+config%isma+is)*(zn(km)-kav(k))+pnuinp(1:config%qx,jm,kmi,is)
          fnutot(1:config%qx,jm,km,is) = &
                dzdx(1:config%qx,jm,k,3+2*config%isma+is)*(zn(km)-kav(k))+fnuinp(1:config%qx,jm,kmi,is)
          dnutot(1:config%qx,jm,km,is) = &
                dzdx(1:config%qx,jm,k,3+3*config%isma+is)*(zn(km)-kav(k))+dnuinp(1:config%qx,jm,kmi,is)
          gnutot(1:config%qx,jm,km,is) = &
                dzdx(1:config%qx,jm,k,3+4*config%isma+is)*(zn(km)-kav(k))+gnuinp(1:config%qx,jm,kmi,is)
        enddo
      enddo ! km-loop
    enddo ! k-loop
  enddo ! jm-loop
endif ! config%nsdim
  return
end subroutine cintp_s
#endif /* not PROGRAM_remap or NOTRA */

#ifndef NOTRA
subroutine extend_grid(j,k)

  use precision

#ifndef NOTRA

  use radial_grid_rt
!      use multigrid_rt   ! forcheck
  use overlap
#endif


!      use mesh_hy
  use totare_hy 
  use intgrs_hy
  use hydro_areas_mod
  use configure
! LOCAL variables that are not in modules

  implicit none
  integer(kind=ik) :: i
  integer(kind=ik) ,intent(in) :: j,k

!      if (nxlb(j) .ge. 1) then
!         if (r(-1).eq.xzltot(nxlb(j))) then
!            rr(-2)       = xzltot(nxlb(j)-1)
!         else
!            rr(-2)       = xzltot(nxlb(j))
!         endif
!      else
!         rr(-2)       = zero
!      endif
  rr(-2)=1.e45_rk

  do i=-1,config%imaxp +1
!     if (config%use_multid_collapse) then
!         rr(i)=ralag(i,j,k)
!     else
!         rr(i)=ralag(i)
!     endif
     rr(i)=r(i)
  enddo
  if (nxub(j) .lt. areas%nx) then
     if (r(config%imaxp +1) .eq. xzrtot(nxub(j))) then
        rr(config%imaxp +1+1) = xzrtot(nxub(j)+1)
     else
        rr(config%imaxp +1+1) = xzrtot(nxub(j))
     endif
  else
     if (r(config%imaxp +1) .ge. xzrtot(areas%nx)) then
        rr(config%imaxp +1+1) = 2._rk*rr(config%imaxp +1)-rr(config%imaxp)
     else
        rr(config%imaxp +1+1) = xzrtot(areas%nx)
     endif
     !        rr(config%imaxp +1+1) = max(xzrtot(nx),2._rk*rr(config%imaxp +1)-rr(imaxp))
  endif

!      write(76,*) rr(config%imaxp +1),rr(config%imaxp +1+1),xzrtot(nxub(j)),nxub(j)


  return
end subroutine extend_grid
#endif /* NOTRA */


#ifdef CFC_TRANSPORT
!=======================================================================

SUBROUTINE init_tot_arrays(jmi,jma,kmi,kma,mode)

!=======================================================================
! Routine uebertraegt Hydro-Arrays aus CoCoNuT in die ...tot-Arrays
! in PROMETHEUS. Vorsicht: Geschwindigkeiten in PROMETHEUS/VERTEX
! werden anders als in CoCoNuT immer in einer orthonormalen
! Tetradenbasis angegeben.
!-----------------------------------------------------------------------

  USE precision

#ifndef NOTRA
#endif
  USE totare_hy
  USE phycon
  USE nutrio_hy

  USE size_cfc, ONLY: m,n,o
  USE parameters_cfc
  USE grid_cfc
  USE perm_aux_cfc
  USE hydro_primitives_cfc
  USE fluxes_cfc
  USE conserved_cfc
  USE metric_cfc
  USE shock_cfc
  use cputim
  use configure
  IMPLICIT NONE


  integer(kind=ik), intent(in) :: jmi, jma, kmi, kma, mode
  real(kind=rk)                ::  paralleltime_start(2), paralleltime_end(2)
  integer(kind=ik)             :: i,j,k,jk,n_zones,o_zones

  paralleltime_start(:) = 0._rk
  paralleltime_end(:)   = 0_rk

  ! Folgende Schleife im Prinzip parallellisierbar, aber Vorsicht
  ! wegen Scoping der Variablen. Nur empfehlenswert und noetig
  ! fuer 1D-Eddington-Faktor.

  n_zones = jma - jmi + 1
  o_zones = kma - kmi + 1

#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
  call second_v(paralleltime_start)
#endif
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j,k,jk) &
!$OMP SHARED(mode,jmi,jma,kmi,kma)
#endif
  do jk = 1, n_zones * o_zones

     k = int((jk + n_zones - 1) / n_zones )
     j = (jmi - 1) + (jk - (k - 1) * n_zones)
     k = k + kmi - 1

     if (mode.ne.-1) then
        do i=1,config%qx
           sqrt_gamma(i,j,k)=phi(i,j,k)**6 
           dentot(i,j,k)=rho(i,j,k)*rho_geom_factor_inv* &
                sqrt_gamma(i,j,k)*w(i,j,k)
           vextot(i,j,k)=v_1(i,j,k)/phi(i,j,k)**2
!           veytot(i,j,k)=0.5*(vav_theta(i,j-1,k)+vav_theta(i,j,k))/ &
!                (phi(i,j,k)**2*r(i))
           veytot(i,j,k)=v_2(i,j,k)/(phi(i,j,k)**2*r(i))
           veztot(i,j,k)=v_3(i,j,k)/(phi(i,j,k)**2*r(i)*sinus_theta(j))
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
           xnutot(i,j,k,1:config%qn)=xnnu(i,j,k,1:config%qn)
           wltot (i,j,k)=w  (i,j,k)
           gpotot(i,j,k)=phi(i,j,k)
        end do
        if (mode.ne.0) then
           do i=1,config%qx
              acxtot(i,j,k)=alpha(i,j,k) !temporaer, sonst acc(i,j,k)
              temtot(i,j,k)=t  (i,j,k)
              enetot(i,j,k)=eps(i,j,k)*pc_cl**2
           end do
        end if

     else                   !copy all the arrays needed in pltout
        do i=1,config%qx
           dentot(i,j,k)=rho    (i,j,k)*pc_geog
           xnutot(i,j,k,1:config%qn)=xnnu(i,j,k,1:config%qn)
           wltot(i,j,k) =w      (i,j,k)
           acxtot(i,j,k)=0.0_rk !acc    (i,j,k)*pc_cl**2
           temtot(i,j,k)=t      (i,j,k)
           enetot(i,j,k)=tau_hat(i,j,k)*pc_geoe/d_cap_hat(i,j,k)
           stotot(i,j,k)=entropy(i,j,k)
           gactot(i,j,k)=gamm   (i,j,k)
           pretot(i,j,k)=p      (i,j,k)*pc_geoe
           vextot(i,j,k)=v_1    (i,j,k)*pc_cl/phi(i,j,k)**2
           veytot(i,j,k)=v_2    (i,j,k)*pc_cl/(phi(i,j,k)**2*r(i))
           veztot(i,j,k)=v_3    (i,j,k)*pc_cl/(phi(i,j,k)**2*r(i)*sinus_theta(j))
#ifdef CFC_TRANSPORT2
           gpotot(i,j,k)=phi    (i,j,k)
#else
           gpotot(i,j,k)=phi_potential(i,j,k)
#endif
           cpotot(i,j,k,1)=cpot(i,j,k,1)
           cpotot(i,j,k,2)=cpot(i,j,k,2)
           cpotot(i,j,k,3)=cpot(i,j,k,3)
           cpotot(i,j,k,4)=cpot(i,j,k,4)
        end do
        where(lshock(:,j,k))
           ishtot(:,j,k)=1
        elsewhere
           ishtot(:,j,k)=0
        end where
     endif
  end do
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
  call second_v(paralleltime_end)
#endif
#endif
  timer%omp_par = timer%omp_par + paralleltime_end - paralleltime_start

  return

END SUBROUTINE init_tot_arrays

!=======================================================================


#ifndef NOTRA
!=======================================================================

SUBROUTINE map_nutra2cfc

!=======================================================================

  USE precision

  USE totare_hy
  USE nutrio_hy
  USE multigrid_rt
  USE phycon

  USE size_cfc, ONLY: m,n,o,n_s,n_e,o_s,o_e,n_loc,o_loc
  USE nutra_cfc
  USE metric_cfc

  use configure

  IMPLICIT NONE

  integer(kind=ik) :: i,is,j,k,jk

#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j,k,jk) &
!$OMP SHARED(enu,fnu,pnu)
#endif
  do jk = 1, n_loc * o_loc

     k = int((jk + n_loc - 1) / n_loc )
     j = (n_s - 1) + (jk - (k - 1) * n_loc)
     k = k + o_s - 1

     do i=1,config%qx
        enu(i,j,k)=sum(enutot(i,j,k,1:config%isma))*pc_egeo
        fnu(i,j,k)=sum(enutot(i,j,k,1:config%isma))*pc_egeo/ pc_cl*phi(i,j,k)**2
        if (dentot(i,j,k).gt.1.0d12) then
           pnu(i,j,k)=sum(pnutot(i,j,k,1:config%isma))*pc_egeo !/ & 
           !                   sqrt_gamma(i,j,k)
        else
           pnu(i,j,k)=0.0_rk
        end if
     end do
  end do
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
!$OMP END PARALLEL DO
#endif

  return

END SUBROUTINE map_nutra2cfc

!=======================================================================
#endif /* NOTRA */



#ifndef NOTRA
!=======================================================================

SUBROUTINE calc_lat_redsh

!=======================================================================
!
! Calculate the exact redshift coefficients in the optically thick regime
! according to optisch_dicht.nb
!
!-----------------------------------------------------------------------

 USE precision

 USE nutrio_hy
 USE rshlat_hy
 USE phycon
 USE backquants_rt

 USE size_cfc, ONLY: m,n,o,n_s,n_e,o_s,o_e,n_loc,o_loc
 USE nutra_cfc
 USE grid_cfc
 USE perm_aux_cfc
 USE hydro_primitives_cfc
 USE metric_cfc

#ifdef MPI_HYDRO
 USE mo_mpi
#endif

 use configure
 IMPLICIT NONE

 real(kind=rk) :: scr(1:config%qx,n_s-1:n_e+1,o_s-1:o_e+1)

#ifdef MPI_HYDRO
 real (kind=rk) :: sbufy (1:m+4, o_s:o_e)
 real (kind=rk) :: rbufy (1:m+4, o_s:o_e)
 real (kind=rk) :: sbufz (1:m+4, n_s:n_e)
 real (kind=rk) :: rbufz (1:m+4, n_s:n_e)

 integer (kind=ik) :: src, dest, ierr, mpistat (MPI_STATUS_SIZE)
 integer (kind=ik), parameter :: tag_ysndl=2000_ik
 integer (kind=ik), parameter :: tag_ysndr=2001_ik
 integer (kind=ik), parameter :: tag_zsndl=2002_ik
 integer (kind=ik), parameter :: tag_zsndr=2003_ik
#endif

 integer(kind=ik) :: i,j,k,jk,jp1,jm1

#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j,k,jk) &
!$OMP SHARED(rshlat,sinus_theta,w,alpha,v_2,phi,rr_inv, &
!$OMP       beta_up_2,scr,cosin_theta,cosin_theta_if)
#endif
 do jk = 1, n_loc * o_loc

    k = int((jk + n_loc - 1) / n_loc )
    j = (n_s - 1) + (jk - (k - 1) * n_loc)
    k = k + o_s - 1

    do i=1,config%qx
       scr(i,j,k)=sinus_theta(j)*w(i,j,k)* &
            (alpha(i,j,k)*v_2(i,j,k)/phi(i,j,k)**4* &
               rr_inv(i)-sqrt_gamma(i,j,k)*beta_up_2(i,j,k))
    end do
 end do

#ifdef MPI_HYDRO
 !send to left neighbour in THETA direction -----------------
 sbufy(1:m,o_s:o_e)=scr(1:m,n_s,o_s:o_e)
 call MPI_Cart_shift(cart_comm,0,-1,src,dest,ierr)
 call MPI_Sendrecv (sbufy,m*o_loc, &
      MPI_DOUBLE_PRECISION,dest,tag_ysndl,rbufy, &
      m*o_loc,MPI_DOUBLE_PRECISION, &
      src,tag_ysndl,cart_comm,mpistat,ierr)
 scr(1:m,n_e+1,o_s:o_e)=rbufy(1:m,o_s:o_e)

 !send to right neighbour in THETA direction -----------------
 sbufy(1:m,o_s:o_e)=scr(1:m,n_e,o_s:o_e)
 call MPI_Cart_shift(cart_comm,0,1,src,dest,ierr)
 call MPI_Sendrecv (sbufy,m*o_loc, &
      MPI_DOUBLE_PRECISION,dest,tag_ysndl,rbufy, &
      m*o_loc,MPI_DOUBLE_PRECISION, &
      src,tag_ysndl,cart_comm,mpistat,ierr)
 scr(1:m,n_s-1,o_s:o_e)=rbufy(1:m,o_s:o_e)
#endif /*MPI_HYDRO*/


#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j,k,jk) &
!$OMP SHARED(rshlat,sinus_theta,w,alpha,v_2,phi,rr_inv, &
!$OMP       beta_up_2,scr,cosin_theta,cosin_theta_if)
#endif
 do jk = 1, n_loc * o_loc

    k = int((jk + n_loc - 1) / n_loc )
    j = (n_s - 1) + (jk - (k - 1) * n_loc)
    k = k + o_s - 1

    if (j .ne. 1 .and. j .ne. config%qy) then

       jp1=min(j+1,config%qy)
       jm1=max(j-1,1)
       do i=1,config%qx
          rshlat(i,j,k)=-0.5*w(i,j,k)* (scr(i,jp1,k)-scr(i,jm1,k))/ &
               (cosin_theta(jp1)-cosin_theta(jm1))
       end do

    else if (j .eq. 1) then

       jp1=j+1
       do i=1,config%qx
          rshlat(i,j,k)=-0.5*w(i,j,k)* (scr(i,jp1,k))/ &
               (cosin_theta(jp1)-cosin_theta_if(0))
       end do

    else if (j .eq. config%qy) then

       jm1=j-1
       do i=1,config%qx
          rshlat(i,j,k)=-0.5*w(i,j,k)* (-scr(i,jm1,k))/ &
               (cosin_theta_if(j)-cosin_theta(jm1))
       end do

    end if
 end do


 if (config%qz .gt. 1) then
    ! TO DO: compute advection velocities in VARPHI direction

#ifdef MPI_HYDRO
    !send to left neighbour in VARPHI direction ----------------
    sbufy(1:m,n_s:n_e)=scr(1:m,n_s:n_e,o_s)
    call MPI_Cart_shift(cart_comm,1,-1,src,dest,ierr)
    call MPI_Sendrecv (sbufy,m*n_loc, &
           MPI_DOUBLE_PRECISION,dest,tag_ysndl,rbufy, &
           m*n_loc,MPI_DOUBLE_PRECISION, &
           src,tag_ysndl,cart_comm,mpistat,ierr)
    scr(1:m,n_s:n_e,o_e+1)=rbufy(1:m,n_s:n_e)

    !send to right neighbour in THETA direction ----------------
    sbufy(1:m,n_s:n_e)=scr(1:m,n_s:n_e,o_e)
    call MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
    call MPI_Sendrecv (sbufy,m*n_loc, &
           MPI_DOUBLE_PRECISION,dest,tag_ysndl,rbufy, &
           m*n_loc,MPI_DOUBLE_PRECISION, &
           src,tag_ysndl,cart_comm,mpistat,ierr)
    scr(1:m,n_s:n_e,o_s-1)=rbufy(1:m,n_s:n_e)
#endif /*MPI_HYDRO*/

    ! TO DO: compute contribution of VARPHI velocities to compression rate

 endif

 return

END SUBROUTINE calc_lat_redsh

!     =======================================================================

#endif /* NOTRA */
#endif /* CFC_TRANSPORT */

end module grids
