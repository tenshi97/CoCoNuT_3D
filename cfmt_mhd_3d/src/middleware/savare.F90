!> \verbatim
!> This module provides the subroutine savare which is overloaded
!> with the CFC or the PROMETHEUS version
!>
!> Here is a list of subroutines which are overloaded depending on whether
!> CFC or PROMETHEUS is uses
!>
!>  interface-name    real-name        compiled when?
!>  savare            savare_CFC       CFC_TRANSPORT
!>  savare            savare_PROM      if not CFC_TRANSPORT
!>
!>  The purpose of savare is to get or put 3D calculation area
!>
!>  Author: A. Marek, MPA, Nov. 2009
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
module savare_overload

  use precision
  use abort

#ifdef CFC_TRANSPORT
  public :: savare, savare_CFC
#else
  public :: savare, savare_PROM
#endif


  interface savare
#ifdef CFC_TRANSPORT
     module procedure savare_CFC
#else
     module procedure savare_PROM
#endif

  end interface

  contains


#ifdef CFC_TRANSPORT

!=======================================================================

!> \verbatim
!>=======================================================================
!>      CFC-Version
!>=======================================================================
!>
!>
!>     task:  get or put 3D calculation area
!>
!> Author : B. Mueller
!>=======================================================================
!> \endverbatim
!>
!> \param  isw  0: copy vexanc -> vextot  (get old values)
!>              1: copy vextot -> vexanc  (put old values)
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
subroutine savare_CFC(isw, savare_self, savare_children)

  use totgrq_hy
  use nutrio_hy
  ! use arecon_hy
  use totare_hy
  use ioflx
  use intgrs_hy

  use size_cfc
  use hydro_primitives_cfc
  use conserved_cfc
  use metric_cfc
  use parameters_cfc

  use metric
  use ancient_hy
  use ancient_cfc


  use hydro_areas_mod
  use configure
  use abort
  use cputim
  IMPLICIT NONE

  real(kind=rk)                :: savare_self(2), savare_children(2)
  real(kind=rk)                :: paralleltime_start(2), paralleltime_end(2)
  real(kind=rk)                :: selftime_start(2)
  integer(kind=ik), intent(in) :: isw

  logical                      :: ler

  real(kind=rk), save          :: sioe_anc, sion_anc, sdegr_anc, &
                                  sdenu_anc, sdynu_anc
  real(kind=rk), save          :: delta_t_anc, delta_t_old_anc,  &
                                  t_total_anc, iteration_anc

  integer                      :: i,j,k,jk

  savare_self     = 0._rk
  savare_children = 0._rk
#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif
  
  select case (isw)
  case (0)                  !restore
     
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
#ifndef DEBUG_TIMINGS
  call second_v(paralleltime_start)
#endif
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,k,jk) &
!$OMP SHARED(n_loc,o_loc,n_s,o_s,n_e,o_e) &
!$OMP DEFAULT(SHARED)

!$OMP DO &
!$OMP SCHEDULE(static)
#endif
     do jk = 1, n_loc * o_loc
        k = int((jk + n_loc - 1) / n_loc )
        j = (n_s - 1) + (jk - (k - 1) * n_loc)
        k = k + o_s - 1

!------ hydro quantities: --------------------------------------------------
        do i=1,config%qx
           rho(i,j,k) = denanc(i,j,k)
           v_1(i,j,k) = vexanc(i,j,k)
           v_2(i,j,k) = veyanc(i,j,k)
           v_3(i,j,k) = vezanc(i,j,k)
           t  (i,j,k) = temanc(i,j,k)

           eps(i,j,k) = epsanc(i,j,k)
           p  (i,j,k) = preanc(i,j,k)
           w  (i,j,k) = wloanc(i,j,k)

           xnnu(i,j,k,1:config%qn) = xnuanc(i,j,k,1:config%qn)

           qentot(i,j,k) = qenanc(i,j,k)
           qmotot(i,j,k) = qmoanc(i,j,k)
           qyetot(i,j,k,1) = qyeanc(i,j,k,1)
           qyetot(i,j,k,2) = qyeanc(i,j,k,2)
           qyetot(i,j,k,3) = qyeanc(i,j,k,3)
           qyetot(i,j,k,4) = qyeanc(i,j,k,4)
           qyetot(i,j,k,5) = qyeanc(i,j,k,5)

           d_cap_hat(i,j,k) = d_cap_hat_anc (i,j,k)
           s_1_hat  (i,j,k) = s_1_hat_anc   (i,j,k)
           s_2_hat  (i,j,k) = s_2_hat_anc   (i,j,k)
           s_3_hat  (i,j,k) = s_3_hat_anc   (i,j,k)
           tau_hat  (i,j,k) = tau_hat_anc   (i,j,k)
           d_cap_xnu_hat (i,j,k,1:config%qn) = d_cap_xnu_hat_anc (i,j,k,1:config%qn)

        end do

!------ metric quantities: -------------------------------------------------
        do i=0,config%qx+1
           phi      (i,j,k) = phi_anc      (i,j,k)
           alpha    (i,j,k) = alpha_anc    (i,j,k)
           beta_up_1(i,j,k) = beta_up_1_anc(i,j,k)
           beta_up_2(i,j,k) = beta_up_2_anc(i,j,k)
           beta_up_3(i,j,k) = beta_up_3_anc(i,j,k)
        end do

        if (config%nsdim .ge. 2) then
        if (j.eq.n_s) then
           do i=0,config%qx+1
              phi      (i,n_s-1,k) = phi_anc      (i,n_s-1,k)
              alpha    (i,n_s-1,k) = alpha_anc    (i,n_s-1,k)
              beta_up_1(i,n_s-1,k) = beta_up_1_anc(i,n_s-1,k)
              beta_up_2(i,n_s-1,k) = beta_up_2_anc(i,n_s-1,k)
              beta_up_3(i,n_s-1,k) = beta_up_3_anc(i,n_s-1,k)
           end do
        end if

        if (j.eq.n_e) then
           do i=0,config%qx+1
              phi      (i,n_e+1,k) = phi_anc      (i,n_e+1,k)
              alpha    (i,n_e+1,k) = alpha_anc    (i,n_e+1,k)
              beta_up_1(i,n_e+1,k) = beta_up_1_anc(i,n_e+1,k)
              beta_up_2(i,n_e+1,k) = beta_up_2_anc(i,n_e+1,k)
              beta_up_3(i,n_e+1,k) = beta_up_3_anc(i,n_e+1,k)
           end do
        end if
     endif

     if (config%nsdim .eq. 3) then
        if (k.eq.o_s) then
           do i=0,config%qx+1
              phi      (i,j,o_s-1) = phi_anc      (i,j,o_s-1)
              alpha    (i,j,o_s-1) = alpha_anc    (i,j,o_s-1)
              beta_up_1(i,j,o_s-1) = beta_up_1_anc(i,j,o_s-1)
              beta_up_2(i,j,o_s-1) = beta_up_2_anc(i,j,o_s-1)
              beta_up_3(i,j,o_s-1) = beta_up_3_anc(i,j,o_s-1)
           end do
        end if

        if (k.eq.o_e) then
           do i=0,config%qx+1
              phi      (i,j,o_e+1) = phi_anc      (i,j,o_e+1)
              alpha    (i,j,o_e+1) = alpha_anc    (i,j,o_e+1)
              beta_up_1(i,j,o_e+1) = beta_up_1_anc(i,j,o_e+1)
              beta_up_2(i,j,o_e+1) = beta_up_2_anc(i,j,o_e+1)
              beta_up_3(i,j,o_e+1) = beta_up_3_anc(i,j,o_e+1)
           end do
        end if
        endif

     end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END PARALLEL
#ifndef DEBUG_TIMINGS
  call second_v(paralleltime_end)
#endif
#endif
  timer%omp_par = timer%omp_par + paralleltime_end - paralleltime_start

!     grid quantities:
     xzltot(1:config%qx) = xzlanc(1:config%qx)
     xzrtot(1:config%qx) = xzranc(1:config%qx)
     xzntot(1:config%qx) = 0.5_rk*(xzlanc(1:config%qx)+xzranc(1:config%qx))
!???
     dvxtot(1:config%qx) = (xzranc(1:config%qx)**3-xzlanc(1:config%qx)**3)/3.0_rk

     yzltot(1:config%qy) = yzlanc(1:config%qy)
     yzrtot(1:config%qy) = yzranc(1:config%qy)
     yzntot(1:config%qy) = 0.5*(yzlanc(1:config%qy)+yzranc(1:config%qy))
!???
! dvytot=
     zzltot(1:config%qz) = zzlanc(1:config%qz)
     zzrtot(1:config%qz) = zzranc(1:config%qz)
     zzntot(1:config%qz) = 0.5*(zzlanc(1:config%qz)+zzranc(1:config%qz))
! ???
! dvztot=

     call calculate_metric

     areas%nhystp = nhystpanc
     areas%ix_are(:,:)=ix_areanc(:,:)

     delta_t=delta_t_anc
     delta_t_old=delta_t_old_anc
     t_total=t_total_anc
     iteration=iteration_anc

! some output Quantities
     sumioe  = sioe_anc
     sumion  = sion_anc
     sumdegr = sdegr_anc
     sumdenu = sdenu_anc
     sumdynu = sdynu_anc

  case (1)                  !save


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,k,jk) &
!$OMP SHARED(n_loc,o_loc,n_s,o_s,n_e,o_e) &
!$OMP DEFAULT(SHARED)

!$OMP DO &
!$OMP SCHEDULE(static)
#endif
     do jk = 1, n_loc * o_loc
        k = int((jk + n_loc - 1) / n_loc )
        j = (n_s - 1) + (jk - (k - 1) * n_loc)
        k = k + o_s - 1

!------ hydro quantities: --------------------------------------------------
        do i = 1, config%qx
           denanc(i,j,k) = rho(i,j,k)
           vexanc(i,j,k) = v_1(i,j,k)
           veyanc(i,j,k) = v_2(i,j,k)
           vezanc(i,j,k) = v_3(i,j,k)
           temanc(i,j,k) = t  (i,j,k)

           epsanc(i,j,k) = eps(i,j,k)
           preanc(i,j,k) = p  (i,j,k)
           wloanc (i,j,k) = w(i,j,k)

           xnuanc(i,j,k,1:config%qn) = xnnu(i,j,k,1:config%qn)

           qenanc(i,j,k) = qentot(i,j,k)
           qmoanc(i,j,k) = qmotot(i,j,k)
           qyeanc(i,j,k,1) = qyetot(i,j,k,1)
           qyeanc(i,j,k,2) = qyetot(i,j,k,2)
           qyeanc(i,j,k,3) = qyetot(i,j,k,3)
           qyeanc(i,j,k,4) = qyetot(i,j,k,4)
           qyeanc(i,j,k,5) = qyetot(i,j,k,5)
           d_cap_hat_anc (i,j,k) = d_cap_hat (i,j,k)
           s_1_hat_anc   (i,j,k) = s_1_hat   (i,j,k)
           s_2_hat_anc   (i,j,k) = s_2_hat   (i,j,k)
           s_3_hat_anc   (i,j,k) = s_3_hat   (i,j,k)
           tau_hat_anc   (i,j,k) = tau_hat   (i,j,k)
           d_cap_xnu_hat_anc (i,j,k,1:config%qn) = d_cap_xnu_hat (i,j,k,1:config%qn)

        end do
     end do


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO &
!$OMP SCHEDULE(static)
#endif
     do jk = 1, n_loc * o_loc

        k = int((jk + n_loc - 1) / n_loc )
        j = (n_s - 1) + (jk - (k - 1) * n_loc)
        k = k + o_s - 1

!------ metric quantities: -------------------------------------------------
        do i=0,config%qx+1
           phi_anc      (i,j,k) = phi      (i,j,k)
           alpha_anc    (i,j,k) = alpha    (i,j,k)
           beta_up_1_anc(i,j,k) = beta_up_1(i,j,k)
           beta_up_2_anc(i,j,k) = beta_up_2(i,j,k)
           beta_up_3_anc(i,j,k) = beta_up_3(i,j,k)
        end do

        if (config%nsdim .ge. 2) then
        if (j.eq.n_s) then
           do i=0,config%qx+1
              phi_anc      (i,n_s-1,k) = phi      (i,n_s-1,k)
              alpha_anc    (i,n_s-1,k) = alpha    (i,n_s-1,k)
              beta_up_1_anc(i,n_s-1,k) = beta_up_1(i,n_s-1,k)
              beta_up_2_anc(i,n_s-1,k) = beta_up_2(i,n_s-1,k)
              beta_up_3_anc(i,n_s-1,k) = beta_up_3(i,n_s-1,k)
           end do
        endif

        if (j.eq.n_e) then
           do i=0,config%qx+1
              phi_anc      (i,n_e+1,k) = phi      (i,n_e+1,k)
              alpha_anc    (i,n_e+1,k) = alpha    (i,n_e+1,k)
              beta_up_1_anc(i,n_e+1,k) = beta_up_1(i,n_e+1,k)
              beta_up_2_anc(i,n_e+1,k) = beta_up_2(i,n_e+1,k)
              beta_up_3_anc(i,n_e+1,k) = beta_up_3(i,n_e+1,k)
           end do
        endif
        endif

        if (config%nsdim .eq. 3) then
        if (k.eq.o_s) then
           do i=0,config%qx+1
              phi_anc      (i,j,o_s-1) = phi      (i,j,o_s-1)
              alpha_anc    (i,j,o_s-1) = alpha    (i,j,o_s-1)
              beta_up_1_anc(i,j,o_s-1) = beta_up_1(i,j,o_s-1)
              beta_up_2_anc(i,j,o_s-1) = beta_up_2(i,j,o_s-1)
              beta_up_3_anc(i,j,o_s-1) = beta_up_3(i,j,o_s-1)
           end do
        endif

        if (k.eq.o_e) then
           do i=0,config%qx+1
              phi_anc      (i,j,o_e+1) = phi      (i,j,o_e+1)
              alpha_anc    (i,j,o_e+1) = alpha    (i,j,o_e+1)
              beta_up_1_anc(i,j,o_e+1) = beta_up_1(i,j,o_e+1)
              beta_up_2_anc(i,j,o_e+1) = beta_up_2(i,j,o_e+1)
              beta_up_3_anc(i,j,o_e+1) = beta_up_3(i,j,o_e+1)
           end do
        endif
        endif

     end do !j,k-loop


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END PARALLEL
#endif

!     grid quantities:
     xzlanc(1:config%qx) = xzltot(1:config%qx)
     xzranc(1:config%qx) = xzrtot(1:config%qx)

     yzlanc(1:config%qy) = yzltot(1:config%qy)
     yzranc(1:config%qy) = yzrtot(1:config%qy)

     zzlanc(1:config%qz) = zzltot(1:config%qz)
     zzranc(1:config%qz) = zzrtot(1:config%qz)

     nhystpanc = areas%nhystp
     ix_areanc(:,:)=areas%ix_are(:,:)

     delta_t_anc=delta_t
     delta_t_old_anc=delta_t_old
     t_total_anc=t_total
     iteration_anc=iteration

! some output Quantities
     sioe_anc   = sumioe
     sion_anc   = sumion
     sdegr_anc  = sumdegr
     sdenu_anc  = sumdenu
     sdynu_anc  = sumdynu

  case default
     raise_abort('savare> mode not implemented')
  end select

#ifndef DEBUG_TIMINGS
 call second_v(savare_self)
 savare_self = savare_self - selftime_start
#endif
  
  return

end subroutine savare_CFC



#else /* CFC_TRANSPORT */

!=======================================================================
!> \verbatim
!>=======================================================================
!>      CFC-Version
!>=======================================================================
!>
!>
!>     task:  get or put 3D calculation area
!>
!> Author : M. Rampp
!>=======================================================================
!> \endverbatim
!>
!> \param  isw  0: copy vexanc -> vextot  (get old values)
!>              1: copy vextot -> vexanc  (put old values)
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
subroutine savare_PROM(isw, selftime, childrentime)

  use precision
  use abort
  use totgrq_hy
  use nutrio_hy
  ! use arecon_hy
  use ioflx
  use totare_hy
  use intgrs_hy
  use mo_mpi
  use eos3d_routine, only : eos3d
  use ancient_hy
  use cpyare_mod


  use hydro_areas_mod
  use cputim
  use configure
  implicit none
! LOCAL variables that are not in modules


  real(kind=rk), intent(out)   :: selftime(2), childrentime(2)
  real(kind=rk)                :: selftime_start(2)

  integer(kind=ik), intent(in) :: isw
  logical                      :: ler
  real(kind=rk), save          :: sioe_anc, sion_anc,   &
                                  sdegr_anc, sdenu_anc, &
                                  sdynu_anc
  real(kind=rk)                :: eos3d_self(2), eos3d_children(2)

  selftime     = 0._rk
  childrentime = 0._rk
#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

  select case (isw)
  case (0)
     dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = denanc(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = vexanc(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = veyanc(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     veztot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = vezanc(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     temtot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = temanc(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     qentot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qenanc(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qmoanc(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1) = qyeanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,1)
     qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,2) = qyeanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,2)
     qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,3) = qyeanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,3)
     qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,4) = qyeanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,4)
     qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,5) = qyeanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,5)
     xnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%qn) = xnuanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%qn)

     gpotot(0:config%qx,qy_s-1:qy_e,qz_s:qz_e) = gpoanc(0:config%qx,qy_s-1:qy_e,qz_s:qz_e)

     gamtot(0:config%qx) = gamanc(0:config%qx)


     xzltot(1:config%qx) = xzlanc(1:config%qx)
     xzrtot(1:config%qx) = xzranc(1:config%qx)
     xzntot(1:config%qx) = 0.5_rk*(xzlanc(1:config%qx)+xzranc(1:config%qx))
!???
     dvxtot(1:config%qx) = (xzranc(1:config%qx)**3-xzlanc(1:config%qx)**3) /3._rk

     yzltot(1:config%qy) = yzlanc(1:config%qy)
     yzrtot(1:config%qy) = yzranc(1:config%qy)
     yzntot(1:config%qy) = 0.5_rk*(yzlanc(1:config%qy)+yzranc(1:config%qy))
!???
!      dvytot=
     zzltot(1:config%qz) = zzlanc(1:config%qz)
     zzrtot(1:config%qz) = zzranc(1:config%qz)
     zzntot(1:config%qz) = 0.5_rk*(zzlanc(1:config%qz)+zzranc(1:config%qz))
! ???
!      dvztot=

     call cpyare(2)       ! EOS cannot use ***tot-arrays
     areas%nx=config%qx
     areas%ny=config%qy
     areas%nz=config%qz

     call eos3d (1,ler, eos3d_self, eos3d_children)
     if (ler) then
        raise_abort("savare(): eos failed")
     endif
     call cpyare(99)

     areas%nhystp = nhystpanc
     areas%ix_are(:,:)=ix_areanc(:,:)

! some output Quantities
     sumioe  = sioe_anc
     sumion  = sion_anc
     sumdegr = sdegr_anc
     sumdenu = sdenu_anc
     sumdynu = sdynu_anc

  case (1)


     denanc(1:config%qx,qy_s:qy_e,qz_s:qz_e) = dentot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     vexanc(1:config%qx,qy_s:qy_e,qz_s:qz_e) = vextot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     veyanc(1:config%qx,qy_s:qy_e,qz_s:qz_e) = veytot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     vezanc(1:config%qx,qy_s:qy_e,qz_s:qz_e) = veztot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     temanc(1:config%qx,qy_s:qy_e,qz_s:qz_e) = temtot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     qenanc(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qentot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     qmoanc(1:config%qx,qy_s:qy_e,qz_s:qz_e) = qmotot(1:config%qx,qy_s:qy_e,qz_s:qz_e)
     qyeanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,1) = qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1)
     qyeanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,2) = qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,2)
     qyeanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,3) = qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,3)
     qyeanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,4) = qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,4)
     qyeanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,5) = qyetot(1:config%qx,qy_s:qy_e,qz_s:qz_e,5)
     xnuanc(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%qn) = xnutot(1:config%qx,qy_s:qy_e,qz_s:qz_e,1:config%qn)
     gpoanc(0:config%qx,qy_s-1:qy_e,qz_s:qz_e) = gpotot(0:config%qx,qy_s-1:qy_e,qz_s:qz_e)

     gamanc(0:config%qx) = gamtot(0:config%qx)


     xzlanc(1:config%qx) = xzltot(1:config%qx)
     xzranc(1:config%qx) = xzrtot(1:config%qx)

     yzlanc(1:config%qy) = yzltot(1:config%qy)
     yzranc(1:config%qy) = yzrtot(1:config%qy)

     zzlanc(1:config%qz) = zzltot(1:config%qz)
     zzranc(1:config%qz) = zzrtot(1:config%qz)

     nhystpanc = areas%nhystp
     ix_areanc(:,:)=areas%ix_are(:,:)

! some output Quantities
     sioe_anc   = sumioe
     sion_anc   = sumion
     sdegr_anc  = sumdegr
     sdenu_anc  = sumdenu
     sdynu_anc  = sumdynu

  case default
     raise_abort("savare(): mode not implemented")

  end select

#ifndef DEBUG_TIMINGS
 call second_v(selftime)
 selftime = selftime - selftime_start
#endif

  return
end subroutine savare_PROM

#endif /*CFC_TRANSPORT*/



end module savare_overload
