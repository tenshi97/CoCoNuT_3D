module eos3d_routine

  private

#if ! defined(PROGRAM_remap)
  public :: eos3d
#ifdef CFC_TRANSPORT
  public :: eos3ds_CFC
#else
  public :: eos3ds_PROM
#endif
#else
  public :: eos3ds_remaper
#endif

            
contains

#if defined(PROGRAM_remap)
subroutine eos3ds_remaper (mode, j, k, qin, return_value, ler)

  use precision
  use oldhydro
  use newhydro
!  use restart, only : use_temp_restart
  use eos_sn2
  use dimensions

  use configure
  implicit none

  integer(kind=ik), intent(in)  :: mode, j, k, qin
  real(kind=rk), intent(inout)  :: return_value(:,:,:)
  real(kind=rk)                 :: eos_self(2), eos_children(2)
  real(kind=rk)                 :: childrentime(2)
  logical, intent(out)          :: ler
  integer(kind=ik)              :: i,n

  real(kind=rk), dimension(qin) :: rho, tmp, ek, ei, p, gamc, s, ccu, cce, ccn, &
                                   ccp, dmy1
  real(kind=rk)                 :: dmy2(qin,2), xn(qin,config%qn), xn_save(qin,config%qn)


  select case (mode)
  case(1)

     do i = 1, config%qx
        rho (i) = den_old(i,j,k)
        tmp (i) = tem_old(i,j,k)
        ek  (i) = 0.5 * (vex_old(i,j,k)**2 + vey_old(i,j,k)**2 + &
                         vez_old(i,j,k)**2)  
     enddo

     do n=1,config%qn
        do i = 1, config%qx
           xn(i,n)      = xnu_old(i,j,k,n)
           xn_save(i,n) = xn(i,n)
        enddo
     enddo


     if (config%use_temp_restart) then
        call eos (rho(1:config%qx), tmp(1:config%qx), xn(1:config%qx,:),   &
                  dmy1(1:config%qx), dmy2(1:config%qx,:), ei(1:config%qx), &
                  p(1:config%qx), gamc(1:config%qx), s(1:config%qx),       &
                  ccu(1:config%qx), cce(1:config%qx), ccn(1:config%qx),    &
                  ccp(1:config%qx), eos_self, eos_children, mode=mode, nsemode=0, ler=ler)
     endif
     childrentime = childrentime + eos_self

     do n=1,config%qn
        do i = 1, config%qx
           if (xn(i,n) .ne. 0._rk) then
              if (abs(xn(i,n)-xn_save(i,n))/xn(i,n) .gt. 1e-6_rk) then
                write (*,*) "err1 ",mode,i,n,xn(i,n),xn_save(i,n),rho(i)/1e11
              endif
           endif
        enddo
     enddo
     
     do i = 1, config%qx
        return_value(i,j,k) = ei (i)/rho(i) + ek(i)
     enddo

  case(2)
     do i = 1, nxnew
        rho (i) = dennew(i,j,k)
        tmp (i) = temnew(i,j,k)
        ek  (i) = 0.5 * (vexnew(i,j,k)**2 + veynew(i,j,k)**2 + &
                         veznew(i,j,k)**2)                   
        ei  (i) = (enenew(i,j,k) - ek(i)) * rho(i)
        ei  (i) = max (ei(i), 1e10_rk)
     enddo
     
     do n=1,config%qn
        do i = 1, nxnew
           xn(i,n) = xnunew(i,j,k,n)
           xn_save(i,n)=xn(i,n)
        enddo
     enddo
     
!         do i = 1, nxnew
!            write(*,'(I5,3(1pe12.5))') i,rho (i),tmp(i),ei(i)
!         enddo
!         stop

   if (config%use_temp_restart) then
      call eos (rho(1:nxnew), tmp(1:nxnew), xn(1:nxnew,:),   &
                dmy1(1:nxnew), dmy2(1:nxnew,:), ei(1:nxnew), &
                p(1:nxnew), gamc(1:nxnew), s(1:nxnew),       &
                ccu(1:nxnew), cce(1:nxnew), ccn(1:nxnew),    &
                ccp(1:nxnew), eos_self, eos_children, mode=mode, nsemode=0, ler=ler)
      childrentime = childrentime + eos_self
   endif

   do i = 1, nxnew
      return_value(i,j,k) = tmp (i)
   enddo


   do n=1,config%qn
      do i = 1, config%qx
         if (xn(i,n) .ne. 0._rk) then
            if (abs(xn(i,n)-xn_save(i,n))/xn(i,n) .gt. 1e-6_rk) then
              write (*,*) "err2 ",mode,i,n,xn(i,n),xn_save(i,n),rho(i)/ 1e11
            endif
         endif
      enddo
   enddo
end select
  
return
end subroutine eos3ds_remaper
#endif /* PROGRAM_remap */

#if !(defined(PROGRAM_remap))
!>=======================================================================
!>      driver for eos
!>=======================================================================
!>
!>
!>     task:  copy 3D arrays to 1D slices thata re used in eos
!>
!> Author : W. Keil
!>=======================================================================
!>
!> \param  mode  eos mode 1 = tem, 2 = energy given
!> \param  err   error flag              
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
      subroutine eos3d (mode,ler, selftime, childrentime)    

      use precision

      use intgrs_hy
      use mo_mpi
      use hydro_areas_mod
      use configure
      use cputim
      implicit none
! LOCAL variables that are not in modules
      real(kind=rk), intent(out)   :: selftime(2), childrentime(2)
      real(kind=rk)                :: selftime_start(2)
      real(kind=rk)                :: eos3ds_self(2), eos3ds_children(2)
      real(kind=rk)                :: paralleltime_start(2), &
                                      paralleltime_end(2)
      integer(kind=ik)             :: jk, k, j
      integer(kind=ik), intent(in) :: mode
      logical, intent(out)         :: ler

      logical                      :: ler_2d(qy_proc,qz_proc)

      selftime           = 0._rk
      childrentime       = 0._rk
      paralleltime_start = 0._rk
      paralleltime_end   = 0._rk
#ifndef DEBUG_TIMINGS
      call second_v(selftime_start)
#endif

      ler_2d=.false.

#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
      call second_v(paralleltime_start)
#endif
!$omp  parallel do                              &
!$omp  default (none)                           &
!$omp  private( jk,k,j, eos3ds_self, eos3ds_children )               &
!$omp  shared ( config,mode,ler_2d,qy_s,qy_proc,qz_s,qz_proc, areas, childrentime) &
!$omp  schedule(static)
#endif
!      do jk = 1, qy_proc * qz_proc
      do jk = 1, min(areas%ny,qy_proc) * min(config%nz,qz_proc)
         k = int( ( jk+qy_proc-1 ) / qy_proc )
         j = jk - ( k-1 ) * qy_proc
         
#ifdef CFC_TRANSPORT
         call eos3ds_CFC(mode,qy_s+j-1,qz_s+k-1,ler_2d(j,k), eos3ds_self, eos3ds_children)
#else
         call eos3ds_PROM(mode,qy_s+j-1,qz_s+k-1,ler_2d(j,k), eos3ds_self, eos3ds_children)
#endif
         childrentime = childrentime + eos3ds_self
      enddo
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
      call second_v(paralleltime_end)
#endif
#endif
      timer%omp_par = timer%omp_par + paralleltime_end - paralleltime_start

      ler=any(ler_2d(:,:))
      if (ler) then
         write(*,*) 'eos3d> eos3ds failed in sectors (j,k):'

         do k=1,qz_proc
            do j=1,qy_proc
               if (ler_2d(j,k)) write(*,*) qy_s+j-1,qz_s+k-1
            end do
         enddo

      endif

#ifndef DEBUG_TIMINGS
      call second_v(selftime)
      selftime = selftime - selftime_start
#endif
      return
    end subroutine eos3d


#ifdef CFC_TRANSPORT
!>=======================================================================
!>    calculate equation of state for given energy, density,
!>     and composition on entire grid.  
!>
!>     CFC version
!>
!> Author : B. Mueller
!>=======================================================================
!>
!> \param  mode  eos mode 1 = tem, 2 = energy given
!> \param  j
!> \param  k
!> \param  err   error flag              
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
      SUBROUTINE eos3ds_CFC( mode,j,k,ler, selftime, childrentime )

      use precision

      use size_cfc
      use hydro_primitives_cfc

      use eos_sn2
      use nucparam
      use phycon
      use hydro_areas_mod
      use configure
      use abort
      use cputim
      IMPLICIT NONE

      integer(kind=ik), intent(in)  :: mode,j,k
      logical, intent(out)          :: ler
      real(kind=rk) , intent(out)   :: selftime(2), childrentime(2)
      real(kind=rk)                 :: selftime_start(2)
      real(kind=rk)                 :: eos_self(2), eos_children(2)

      real(kind=rk) :: den(m),tmp(m),ek(m),ei(m),pre(m),e(m),s(m), &
                       gamc(m),xn(m,config%qn),kap(m),xnusum(m)

! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=rk) :: ccu(m), cce(m), ccn(m), ccp(m)
      real(kind=rk) :: dmy1(m), dmy2(m,2) !dummys
! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer(kind=ik) :: i,i_n

      selftime     = 0._rk
      childrentime = 0._rk
#ifndef DEBUG_TIMINGS
      call second_v(selftime_start)
#endif
      ler=.true.

!     Work-Arrays anlegen und in cgs-Einheiten konvertieren:
      do i = 1, m
         den (i) = rho(i,j,k)*pc_geog
         ei  (i) = (rho(i,j,k)*eps(i,j,k))*pc_geoe
         tmp (i) = t  (i,j,k)
      end do

      do i_n=1, config%qn
         do i = 1, m
            xn(i,i_n) = xnnu(i,j,k,i_n)
         end do
      end do
      xnusum (1:m) = 0.0_rk
      do i_n=1, config%qn-1
         do i = 1, m
            xnusum(i) = xnusum (i) + xnnu(i,j,k,i_n)
         end do
      end do
      do i_n=1, config%qn-1
         do i = 1, m
            xn(i,i_n) = xnnu(i,j,k,i_n) / xnusum(i)
         end do
      end do

      
      call eos( den(1:m),tmp(1:m),xn(1:m,:),dmy1(1:m),        &
                dmy2(1:m,:),ei(1:m),pre(1:m),gamc(1:m),s(1:m),&
                ccu(1:m),cce(1:m),ccn(1:m),ccp(1:m), eos_self, eos_children,          &
                mode=mode,nsemode=0,ler=ler)
      childrentime = childrentime + eos_self

      if (ler) then
         write (*,*) 'ERROR: eos3ds failed in sector',j,k
         write (*,*) 'mode =',mode
         write (*,*) 'den - tmp - y_e - ei'
         do i=1,m
            write(*,'(I6,4D16.7)') i,den(i),tmp(i),xn(i,config%qn),ei(i)
         end do
         raise_abort("ERROR: eos3ds failed")
      end if


! --  Zurueckschreiben und in geometrische Einheiten konvertieren:
      do i = 1, m
         if (mode .eq. 1 .or. mode .eq. 4) then
            eps (i,j,k) = ei (i)*pc_egeo/rho(i,j,k)
         end if
         p      (i,j,k) = pre(i)*pc_egeo
         t      (i,j,k) = tmp(i)
         entropy(i,j,k) = s  (i)
         gamm   (i,j,k) = gamc(i)
         xnnu   (i,j,k,:)= xn(i,:)
!         gammac (i,j,k) = gamc(i)
!         gammae (i,j,k) = p(i) / ei(i) + 1.0
!         ! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         cpot (i,j,k,1)  = ccu(i)
         cpot (i,j,k,2)  = cce(i)
         cpot (i,j,k,3)  = ccn(i)
         cpot (i,j,k,4)  = ccp(i)
!         ! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Zusaetzliche Groessen fuer den CFC-Code:
#ifdef CFC_TRANSPORT2
         h           (i,j,k) = 1.0+eps(i,j,k)+p(i,j,k)/rho(i,j,k)
         c_sound_squared(i,j,k) = gamc(i)*p(i,j,k)/ (rho(i,j,k)*h(i,j,k))
#else
!         v_squared (i, j, k) = &
!              v_1(i,j,k)**2*g_up_11(i,j,k)+ &
!              v_2(i,j,k)**2*g_up_22(i,j,k)+ &
!              v_3(i,j,k)**2*g_up_33(i,j,k)
         h           (i,j,k) = eps(i,j,k)+p(i,j,k)/rho(i,j,k)+ &
              0.5_rk*v_squared(i,j,k)
         c_sound_squared(i,j,k) = gamc(i)*p(i,j,k)/rho(i,j,k)
#endif /* CFC_TRANSPORT2 */
      end do

      ler=.false.
#ifndef DEBUG_TIMINGS
      call second_v(selftime)
      selftime = selftime - selftime_start
#endif
      return

      END SUBROUTINE eos3ds_CFC

#else /* CFC_TRANSPORT */

! --------------------------------------------------------------
!>=======================================================================
!>    calculate equation of state for given energy, density,
!>     and composition on entire grid.  
!>
!>     PROMETHEUS version
!>
!> Author : W. Keil
!>=======================================================================
!>
!> \param  mode  eos mode 1 = tem, 2 = energy given
!> \param  j
!> \param  k
!> \param  err   error flag              
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
      subroutine eos3ds_PROM (mode, j, k, ler, selftime, childrentime)

      use precision
      use intgrs_hy
      use nutrio_hy
      use vnew_hy
      use gfloat_hy
      use hydro_hy
      use physcs_hy
      use eos_sn2
      use cputim
      use hydro_areas_mod
      use configure
! LOCAL variables that are not in modules

      implicit none
!$    integer(kind=ik) OMP_GET_THREAD_NUM
      real(kind=rk), intent(out)      :: selftime(2), childrentime(2)
      real(kind=rk)                   :: selftime_start(2)
      integer(kind=ik)                :: i, n
      integer(kind=ik), intent(in)    :: mode, j, k
      logical, intent(out)            :: ler

      real(kind=rk)                   :: eos_self(2), eos_children(2)

      real(kind=rk)                   :: ccu(config%q), cce(config%q), &
                                         ccn(config%q), ccp(config%q)
      real(kind=rk)                   :: dmy1(config%q), dmy2(config%q,2) !dummys
      selftime     = 0._rk
      childrentime = 0._rk
#ifndef DEBUG_TIMINGS
      call second_v(selftime_start)
#endif
      do i = 1, areas%nx
         rho (i) = densty(i,j,k)
         tmp (i) = temp  (i,j,k)
         ek  (i) = 0.5_rk * (velx(i,j,k)**2 + vely(i,j,k)**2 + &
                             velz(i,j,k)**2                  )
         ei  (i) = (energy(i,j,k) - ek(i)) * rho(i)
         ei  (i) = max (ei(i), config%smallp)
      enddo

      do n=1,config%qn
         do i = 1, areas%nx
            xn(i,n) = xnuc(i,j,k,n)
         enddo
      enddo


      call eos (rho(1:areas%nx),tmp(1:areas%nx),xn(1:areas%nx,:),  &
               dmy1(1:areas%nx),dmy2(1:areas%nx,:),ei(1:areas%nx), &
                  p(1:areas%nx),gamc(1:areas%nx),s(1:areas%nx),    &
                ccu(1:areas%nx),cce(1:areas%nx),ccn(1:areas%nx),   &
                ccp(1:areas%nx), eos_self, eos_children,mode=mode,nsemode=0,ler=ler)
      childrentime = childrentime + eos_self
      if (ler) then
#ifndef DEBUG_TIMINGS
         call second_v(selftime)
         selftime = selftime - selftime_start
#endif
         return
      endif
!
! --  copy back
 
      do i = 1, areas%nx
         energy(i,j,k) = ei(i) / rho(i) + ek(i)
         press (i,j,k) = p   (i)
         temp  (i,j,k) = tmp (i)
         stot  (i,j,k) = s   (i)
         gammac(i,j,k) = gamc(i)
         gammae(i,j,k) = p(i) / ei(i) + 1.0_rk
         ! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         cpo(i,j,k,1)  = ccu(i)
         cpo(i,j,k,2)  = cce(i)
         cpo(i,j,k,3)  = ccn(i)
         cpo(i,j,k,4)  = ccp(i)
         ! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo

!  the mass fractions are changed by eos in the high-density regime
      do n = 1, config%qn-1
         do i = 1, areas%nx
            xnuc(i,j,k,n) = xn(i,n)
         enddo
      enddo
#ifndef DEBUG_TIMINGS
      call second_v(selftime)
      selftime = selftime - selftime_start
#endif
     return
    end subroutine eos3ds_PROM

#endif /* CFC_TRANSPORT */
  
#endif /* PROGRAM_remap */
end module eos3d_routine
