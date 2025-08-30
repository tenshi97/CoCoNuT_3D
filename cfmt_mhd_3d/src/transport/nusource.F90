module nusource_mod

implicit none

contains

#undef DEBUG

subroutine nusource(dt_n,mode, selftime, childrentime)

!-----------------------------------------
! Autor           : Markus Rampp (MPA) 
! Modul           : $Id: rady.F,v 1.15 2006/05/11 15:44:04 rburas Exp $
! Version         : $Revision: 1.15 $
! Date            : $Date: 2006/05/11 15:44:04 $
!
!     task        :  
!
!     input:   
!
!-----------------------------------------
  use precision
  use abort
  use error
  
  use phycon
  use param_rt
  use intgrs_hy
  use totare_hy
!  use arecon_hy
  use nutrio_hy
  use ioflx
  use totgrq_hy
  use massio_hy
  use gfloat_hy
  use vnuw_hy
  use hlpare_hy
  use specfun

#ifndef NOTRA
  use qave_overload
#endif

#ifdef CFC_TRANSPORT
  use size_cfc
  use conserved_cfc, ONLY: s_1_hat,s_2_hat,s_3_hat,d_cap_xnu_hat, &
                           tau_hat,d_cap_hat !debug
  use hydro_primitives_cfc, ONLY: v_1,p !debug
  use metric_cfc, ONLY: alpha

  use grids, only : init_tot_arrays
  use cons_check
  use recover_prim_vars
#ifdef CFC_TRANSPORT2
  use multigrid_rt, ONLY: sthe
#endif 
#endif /* CFC_TRANSPORT2 */
  use mo_mpi
  use eos3d_routine, only : eos3d
  use cpyare_mod
      
  use nusource_data
 

  use hydro_areas_mod
  use cputim
  use configure
  implicit none
! LOCAL variables that are not in modules

  real(kind=rk), intent(out)    :: selftime(2), childrentime(2)
  real(kind=rk)                 :: selftime_start(2), eos3d_self(2), &
                                   eos3d_children(2)
  real(kind=rk)                 :: paralleltime_start(2), paralleltime_end(2)
#ifdef CFC_TRANSPORT2
  real(kind=rk)                 :: wm1_v2w, v2, winv
#endif

  integer(kind=ik)              :: i, j, k, jk, idtmin
!-----------------------------------------------------------------------
!     quantities of a distinct calculation area:
!-----------------------------------------------------------------------
  integer(kind=ik) ,intent(in) :: mode
  real(kind=rk), intent(in)    :: dt_n
  logical                      :: ler
  integer(kind=ik)             :: ia

  selftime           = 0._rk
  childrentime       = 0._rk
  paralleltime_start = 0._rk
  paralleltime_end   = 0._rk

  if (config%p_nbk.eq.0) return
#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

  select case (mode)
         
!     In dieser Routine: 1:config%qy->qy_s:qy_e

  case(0)


     syeold(1:config%qx,qy_s:qy_e,qz_s:qz_e)= syeold(1:config%qx,qy_s:qy_e,qz_s:qz_e) + dt_n* &
                                        qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1)
     syeold1(1:config%qx,qy_s:qy_e,qz_s:qz_e)= syeold1(1:config%qx,qy_s:qy_e,qz_s:qz_e) + dt_n* &
                                         qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,2)
     syeold2(1:config%qx,qy_s:qy_e,qz_s:qz_e)= syeold2(1:config%qx,qy_s:qy_e,qz_s:qz_e) + dt_n* &
                                         qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,3)
     syeold3(1:config%qx,qy_s:qy_e,qz_s:qz_e)= syeold3(1:config%qx,qy_s:qy_e,qz_s:qz_e) + dt_n* &
                                         qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,4)
     syeold4(1:config%qx,qy_s:qy_e,qz_s:qz_e)= syeold4(1:config%qx,qy_s:qy_e,qz_s:qz_e) + dt_n* &
                                         qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,5)
     
     smoold(1:config%qx,qy_s:qy_e,qz_s:qz_e)= smoold(1:config%qx,qy_s:qy_e,qz_s:qz_e) + dt_n* &
                                        qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     smyold(1:config%qx,qy_s:qy_e,qz_s:qz_e)= smyold(1:config%qx,qy_s:qy_e,qz_s:qz_e) + dt_n* &
                                        qmytot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     senold(1:config%qx,qy_s:qy_e,qz_s:qz_e)= senold(1:config%qx,qy_s:qy_e,qz_s:qz_e) + dt_n* &
                                        qentot(1:config%qx,qy_s:qy_e,qz_s:qz_e)

#ifndef DEBUG_TIMINGS
     call second_v(selftime)
     selftime = selftime - selftime_start
#endif
     return

  case(1)



     qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1) = -syeold(1:config%qx,qy_s:qy_e,qz_s:qz_e)/dt_n
     qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,2) = -syeold1(1:config%qx,qy_s:qy_e,qz_s:qz_e)/dt_n
     qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,3) = -syeold2(1:config%qx,qy_s:qy_e,qz_s:qz_e)/dt_n
     qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,4) = -syeold3(1:config%qx,qy_s:qy_e,qz_s:qz_e)/dt_n
     qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,5) = -syeold4(1:config%qx,qy_s:qy_e,qz_s:qz_e)/dt_n
!    qmotot  (1:config%qx,qy_s:qy_e,qz_s:qz_e) = -smoold(1:config%qx,qy_s:qy_e,qz_s:qz_e)/dt_n
     qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e)  = 0.0_rk
     qmytot(1:config%qx,qy_s:qy_e,qz_s:qz_e)  = 0.0_rk
     qentot(1:config%qx,qy_s:qy_e,qz_s:qz_e)  = -senold(1:config%qx,qy_s:qy_e,qz_s:qz_e)/dt_n


  case(2)


     qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                                         - smoold(1:config%qx,qy_s:qy_e,qz_s:qz_e)/dt_n
     qmytot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qmytot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                                          - smyold(1:config%qx,qy_s:qy_e,qz_s:qz_e)/dt_n

     syeold(1:config%qx,qy_s:qy_e,qz_s:qz_e)  = 0.0_rk
     smoold(1:config%qx,qy_s:qy_e,qz_s:qz_e)  = 0.0_rk
     smyold(1:config%qx,qy_s:qy_e,qz_s:qz_e)  = 0.0_rk
     senold(1:config%qx,qy_s:qy_e,qz_s:qz_e)  = 0.0_rk


  case(99)
     syeold(:,:,:)  = 0._rk
     syeold1(:,:,:) = 0._rk
     smoold(:,:,:)  = 0._rk
     smyold(:,:,:)  = 0._rk
     senold(:,:,:)  = 0._rk
#ifndef DEBUG_TIMINGS
     call second_v(selftime)
     selftime = selftime - selftime_start
#endif
     return
  case default
     raise_abort("nusource(): wrong mode")
  end select


!-----------------------------------------------------------------------
!     Apply sourceterms to all areas
!-----------------------------------------------------------------------

!- first calculate angular averages of sourceterms

#ifndef NOTRA
  if (config%p_ntr .ne. 0) then
     do ia = 1, areas%are_nu 
        ! Beware here we call the overloaded subroutine qave, whcih is either
        ! the purely OPENMP qave_are_OPENMP or the hybrid MPI/OPENMP 
        ! qave_are_MPI_HYDRO version. The overloading is done in module qave_overload
        call qave_are(ia)
     enddo
  endif
#endif

#ifdef CFC_TRANSPORT
!     CoCoNuT-Version:
      
  call init_tot_arrays(n_s,n_e,o_s,o_e,0)

#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
  call second_v(paralleltime_start)
#endif
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j,k,jk,v2,wm1_v2w,winv) &
!$OMP SHARED(config,s_1_hat,s_2_hat,s_3_hat,tau_hat, &
!$OMP        d_cap_xnu_hat,dt_n,vextot,veytot,veztot,qentot,qyetot, &
!$OMP        qmotot,qmytot,gpotot,xzntot)
#endif
  do jk = 1, n_loc * o_loc
     k = int((jk + n_loc - 1) / n_loc )
     j = (n_s - 1) + (jk - (k - 1) * n_loc)
     k = k + o_s - 1

     do i=1,config%qx

#ifdef CFC_TRANSPORT2
!     Relativistic version:

        if (config%nsdim .ge. 2) then
           v2=vextot(i,j,k)**2+veytot(i,j,k)**2+veztot(i,j,k)**2
           wm1_v2w=0.5_rk*(1.0+v2*(0.25_rk+v2*(0.125_rk+              &
                v2*(0.078125_rk+v2*(0.0546875_rk+v2*(0.041015625_rk + &
                v2*(0.0322265625_rk+v2*0.02618408203125_rk)))))))
           if (v2.gt.0.01_rk) wm1_v2w= (wltot(i,j,k)-1.0d0)/(v2*wltot(i,j,k))
           winv=1.0/wltot(i,j,k)
        endif

        if (config%nsdim .ge. 2) then
           s_1_hat(i,j,k)= s_1_hat(i,j,k) + dt_n*pc_egeo*            &
                (pc_cl*qmotot(i,j,k)*(winv+vextot(i,j,k)**2*wm1_v2w) &
                +vextot(i,j,k)*qentot(i,j,k))*gpotot(i,j,k)**2
        else
           s_1_hat(i,j,k)= s_1_hat(i,j,k) + dt_n*pc_egeo*            &
                (pc_cl*qmotot(i,j,k)                                 &
                +vextot(i,j,k)*qentot(i,j,k))*gpotot(i,j,k)**2
        endif

         if (config%nsdim .ge. 2) then
            s_2_hat(i,j,k)= s_2_hat(i,j,k) + dt_n*pc_egeo*                &
                  (pc_cl*qmytot(i,j,k)                                    &
                 +pc_cl*qmotot(i,j,k)*vextot(i,j,k)*veytot(i,j,k)*wm1_v2w &
                 +veytot(i,j,k)*qentot(i,j,k))*gpotot(i,j,k)**2*xzntot(i)
        
            s_3_hat(i,j,k)= s_3_hat(i,j,k)+dt_n*pc_egeo*                  &
                  (+pc_cl*qmotot(i,j,k)*                                  &
                 vextot(i,j,k)*veztot(i,j,k)*wm1_v2w                      &
                 +veztot(i,j,k)*qentot(i,j,k))*gpotot(i,j,k)**2*xzntot(i) &
                 *sthe(j)
         endif


         if (config%nsdim .ge. 2) then
            tau_hat(i,j,k)= tau_hat(i,j,k) + dt_n*pc_egeo* &
                 (             qentot(i,j,k)+              &
                 veytot(i,j,k)*qmytot(i,j,k)*pc_cl+        &
                 vextot(i,j,k)*qmotot(i,j,k)*pc_cl) 
         else
            tau_hat(i,j,k)= tau_hat(i,j,k) + dt_n*pc_egeo* &
                 (             qentot(i,j,k)+              &
                 vextot(i,j,k)*qmotot(i,j,k)*pc_cl) 
         endif


        d_cap_xnu_hat(i,j,k,config%qn)= d_cap_xnu_hat(i,j,k,config%qn) &
                                       + dt_n*pc_ggeo*qyetot(i,j,k,1)
#else /* CFC_TRANSPORT2 */
!     Newtonian version:
        s_1_hat(i,j,k)= s_1_hat(i,j,k)+ dt_n*pc_egeo*(       pc_cl*qmotot(i,j,k)) 

        if (config%nsdim .ge. 2) then
           s_2_hat(i,j,k)= s_2_hat(i,j,k)+ dt_n*pc_egeo* &
                (        pc_cl*qmytot(i,j,k))*xzntot(i)
        endif

        if (config%nsdim .ge. 2) then
           tau_hat(i,j,k)= tau_hat(i,j,k)+ dt_n*pc_egeo* &
                (             qentot(i,j,k)+  &
                veytot(i,j,k)*qmytot(i,j,k)*pc_cl+ &
                vextot(i,j,k)*qmotot(i,j,k)*pc_cl)
        else
           tau_hat(i,j,k)= tau_hat(i,j,k)+ dt_n*pc_egeo* &
                (             qentot(i,j,k)+  &
                vextot(i,j,k)*qmotot(i,j,k)*pc_cl)
        endif

              
        d_cap_xnu_hat(i,j,k,config%qn)= d_cap_xnu_hat(i,j,k,config%qn) + dt_n*pc_ggeo*qyetot(i,j,k,1)

#endif /*CFC_TRANSPORT2*/

     end do

#ifdef AXIS
     if (j .eq. 1 .or. j .eq. config%qy) then
        do i=1,config%qx
           s_2_hat (i,j,k)=0.0
        end do
     end if
#endif

  end do
#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
  call second_v(paralleltime_end)
#endif
#endif
  timer%omp_par = timer%omp_par + paralleltime_end - paralleltime_start

#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
  call second_v(paralleltime_start)
#endif
!$OMP PARALLEL
#endif
  
      call recover_primitives

#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END PARALLEL
#ifndef DEBUG_TIMINGS
  call second_v(paralleltime_end)
#endif
#endif
  timer%omp_par = timer%omp_par + paralleltime_end - paralleltime_start
  

#else /* CFC_TRANNSPORT */

!     PROMETHEUS-Version:


#ifdef DAMPNSV
      idtmin = ISRMIN_V(are_nu,dt_cfl,1)  ! use CFL time scale
      vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e)= vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                                        + dt_n*( qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e)/ &
                                        dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) ) &
                                        * 0.9**ix_are(idtmin,11)
#else
      vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e)= vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                                         + dt_n*( qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e)/ &
                                         dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) )
#endif /*DAMPNSV*/

      veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e)= veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                                        + dt_n*( qmytot(1:config%qx,qy_s:qy_e,qz_s:qz_e)/ &
                                        dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) )

      enetot(1:config%qx,qy_s:qy_e,qz_s:qz_e)= enetot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                                        + dt_n*( veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
                                        *qmytot(1:config%qx,qy_s:qy_e,qz_s:qz_e)/ &
                                        dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) )
      
!      veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = ?
!      veztot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = ?

      enetot(1:config%qx,qy_s:qy_e,qz_s:qz_e)= enetot(1:config%qx,qy_s:qy_e,qz_s:qz_e) &
#ifdef NEW_QMO
                                        + dt_n*( (qentot(1:config%qx,qy_s:qy_e,qz_s:qz_e)+ &
                                        qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e)* &
                                        vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e))/ &
#else
                                        + dt_n*( qentot(1:config%qx,qy_s:qy_e,qz_s:qz_e)/ &
#endif /*NEW_QMO*/

                                        dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) )

     xnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,config%qn)=xnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,config%qn) &
                                          + dt_n*( qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1)/ &
                                            dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) )


     call cpyare(2)       ! EOS cannot use ***tot-arrays 
     areas%nx=config%qx
     areas%ny=config%qy
     areas%nz=config%qz


     call eos3d (2,ler, eos3d_self, eos3d_children)
     if (ler) raise_error("nusource(): eos failed")
     call cpyare(99)  
#endif /*CFC_TRANSPORT*/

#ifndef DEBUG_TIMINGS
     call second_v(selftime)
     selftime = selftime - selftime_start
#endif
     return
end subroutine nusource

end module nusource_mod
