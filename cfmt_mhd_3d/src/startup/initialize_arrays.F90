module initialize_array_mod

implicit none

contains

!> \par initialzes all main hydro and transport arrays
!>
!> \author unknown
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!> \endverbatim
!>
subroutine initialize_arrays
  use precision

  use ioflx 
  use dquants 
  use cnsv 
!#if defined(SGI) || defined(OPA)

#ifndef NOTRA

!  use radial_grid_rt
!  use energy_grid_rt
!  use radfield_rt
  use param_rt
  use multigrid_rt
!  use backquants_rt
  use lagradq_rt
  use lagrold_rt
!  use tanray_grids_rt
!  use matkhelp_rt
!  use matkoeff_rt
!  use lotdqua_rt
  use boundary_rt
!  use neskernel_rt
!  use pairkernel_rt
  use sourceterms_rt
  use time_rt
!  use theta_grid_rt
!  use speccut_rt
  use lagback_rt 
!  use neutrinotypes

#endif /* NOTRA */

!  use phycon

!  use charac
!  use intgrs_hy
!  use totgrq_hy
!  use arecon_hy
  use nutrio_hy
!  use marker_hy
!  use massio_hy
!  use revsho_hy
!  use gfloat_hy

  use totare_hy

  use vnew_hy
!#endif /* SGI or OPA */
  use nusource_data

  use vnuw_hy

  use mo_mpi

#ifdef CHECK_THREAD_AFFINITY
  use thread_affinity
#endif


  use configure
  IMPLICIT NONE

  integer(kind=ik) :: i,j,k,l,jk,n1,n2,n3,m

#ifdef FCNC_CALC
  real(kind=rk) GTMR_delta

  open(15,file='GTMR_delta')
  read(15,*) GTMR_delta
  read(15,*) fcnc_fac
  close(15)

  fcnc_delta(1)=GTMR_delta
  fcnc_delta(2)=GTMR_delta
  fcnc_delta(3)=0.
  
  if (fcnc_fac.eq.2.) then
     fcnc_delta(1:2) = 0.
     fcnc_delta(3)   = GTMR_delta
     fcnc_fac = 0.
  endif
  if (myproc .eq. 0) then
     
     write (*,*) "Task ",myproc,' GTMR, mode ',fcnc_fac,', delta = ',fcnc_delta(1:3)
  endif


  if (isma.eq.1) fcnc_fac = 0.
#endif /* FCNC_CALC */

     do m = 1, 5
#if defined(OPENMP_HYD) && (defined(OPEN_MP_3D))
!$omp parallel do private(j,k)
#endif

        do k = qz_s, qz_e

#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D))
!$omp parallel do private(j)
#endif

           do j = qy_s, qy_e 
              qye (1:config%qx,j,k,m) = 0._rk
           enddo
        enddo
     enddo

#if defined(OPENMP_HYD) && (defined(OPEN_MP_3D))
!$omp parallel do private(j,k)
#endif
     do k = qz_s, qz_e
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D))
!$omp parallel do private(j)
#endif
        
        do j = qy_s, qy_e 
           qen (1:config%qx,j,k) = 0._rk
           qmo (1:config%qx,j,k) = 0._rk
           qmy (1:config%qx,j,k) = 0._rk
        enddo
     enddo

#ifndef NOTRA

  if (config%p_ntr .ne. 0) then

     delyem_o=0.0_rk
     delden_o=0.0_rk
     deltem_o=0.0_rk
     deleni_o=0.0_rk
     sigma_o=0.0_rk
  
     sumioe=0.0_rk
     sumion=0.0_rk
     sumiom=0.0_rk
     sumioe_n=0.0_rk
     sumion_n=0.0_rk
     
     slep(:)=0.0_rk
     dtot(:)=0.0_rk !?
     
     eges_old=0.0_rk
     eges_old2=0.0_rk
     nges_old=0.0_rk
     nges_old2=0.0_rk
     tten_old=0.0_rk
     elto_old=0.0_rk

     
  endif ! p_ntr

#endif /* NOTRA */


!#if defined(SGI) || defined(OPA)

!*     ---------------------------------
!*     Perform a parallelization over jk


  n1 = qz_proc
  n2 = qy_proc
  n3 = config%qx
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel do shared(n1,n2,n3) private(i,j,k,jk)
#endif
  do jk = 1, n1*n2

     k  = int( ( jk+n2-1 ) / n2 )
     j  = jk - (k-1) * n2 + (qy_s - 1)
     k  = k + (qz_s-1)


     do i = 1, n3
        ishck (i,j,k) = 0
        velx  (i,j,k) = 0._rk
        vely  (i,j,k) = 0._rk
        velz  (i,j,k) = 0._rk
        densty(i,j,k) = 0._rk
        energy(i,j,k) = 0._rk
        press (i,j,k) = 0._rk
        temp  (i,j,k) = 0._rk
        gpot  (i,j,k) = 0._rk
        
        vextot(i,j,k) = 0._rk
        vextot(i,j,k) = 0._rk
        veytot(i,j,k) = 0._rk
        veztot(i,j,k) = 0._rk
        vexold(i,j,k) = 0._rk
        veyold(i,j,k) = 0._rk
        vezold(i,j,k) = 0._rk
              
        dnsold(i,j,k) = 0._rk
        acxtot(i,j,k) = 0._rk
        qgrtot(i,j,k) = 0._rk
        dentot(i,j,k) = 0._rk
        enetot(i,j,k) = 0._rk
        pretot(i,j,k) = 0._rk
        
        qyetot(i,j,k,1) = 0._rk
        qyetot(i,j,k,:) = 0._rk
        qentot(i,j,k) = 0._rk
        qmotot(i,j,k) = 0._rk
        qmytot(i,j,k) = 0._rk

        gaetot(i,j,k) = 0._rk
        gactot(i,j,k) = 0._rk
        stotot(i,j,k) = 0._rk
        temtot(i,j,k) = 0._rk
        gpotot(i,j,k) = 0._rk
        gpoold(i,j,k) = 0._rk

        syeold(i,j,k)  = 0._rk
        syeold1(i,j,k) = 0._rk
        syeold2(i,j,k) = 0._rk       
        syeold3(i,j,k) = 0._rk
        syeold4(i,j,k) = 0._rk           

        smoold(i,j,k) = 0._rk   
        smyold(i,j,k) = 0._rk   
        senold(i,j,k) = 0._rk  
#ifndef NOTRA
        qyeinp(i,j,k) = 0._rk
        qeninp(i,j,k) = 0._rk
        qmoinp(i,j,k) = 0._rk
#endif /* NOTRA */
#ifndef NOTRA
        if (config%p_ntr .ne. 0) then
           thphmass(i,j,k,0) = 0._rk
           thphmass(i,j,k,1) = 0._rk
        endif
#endif /* NOTRA */
     enddo
     dmdtiotot(1:6,j,k) = 0._rk
  enddo

  do l = 1, qc
#if defined(OPENMP_HYD) &&  (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel do shared(n1,n2,n3) private(i,j,k,jk)
#endif
     do jk = 1, n1*n2
        k  = int( ( jk+n2-1 ) / n2 )
        j  = jk - (k-1) * n2 + (qy_s - 1)
        k  = k + (qz_s -1)
        do i=1,n3
           cpotot(i,j,k,l)=0._rk
           cpo(i,j,k,l)=0._rk
        enddo
     enddo
  enddo

  do l = 1, config%isma
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel do shared(n1,n2,n3) private(i,j,k,jk)
#endif
     do jk = 1, n1*n2
        k  = int( ( jk+n2-1 ) / n2 )
        j  = jk - (k-1) * n2 + (qy_s - 1)
        k  = k + (qz_s -1)


        do i=1,n3
           enutot(i,j,k,l)=0._rk
           fnutot(i,j,k,l)=0._rk
           pnutot(i,j,k,l)=0._rk
           dnutot(i,j,k,l)=0._rk
           gnutot(i,j,k,l)=0._rk
#ifndef NOTRA
           enuinp(i,j,k,l)=0._rk
           fnuinp(i,j,k,l)=0._rk
           pnuinp(i,j,k,l)=0._rk
           dnuinp(i,j,k,l)=0._rk
           gnuinp(i,j,k,l)=0._rk
#endif /* NOTRA */
        enddo
     enddo
  enddo

  do l = 1, config%qn
#if defined(OPENMP_HYD) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel do shared(n1,n2,n3) private(i,j,k,jk)
#endif
     do jk = 1, n1*n2
        k  = int( ( jk+n2-1 ) / n2 )
        j  = jk - (k-1) * n2 + (qy_s - 1)
        k  = k + (qz_s -1)


        do i=1,n3
           xnuc  (i,j,k,l) = 0._rk
           xnutot(i,j,k,l) = 0._rk
        enddo
     enddo
  enddo

#ifndef NOTRA

  if (config%p_ntr .ne. 0) then

     do m = 0, 9
#if defined(OPENMP_TRANSPORT) && defined(OPEN_MP_3D)
!$omp parallel do private(j,k)
#endif
        do k = nzmoms, nzmome

#if defined(OPENMP_TRANSPORT) && defined(OPEN_MP_2D)
!$omp parallel do private(j)
#endif
           do j = nymoms, nymome
              wwlag(0:config%imaxp+1,config%iemax,config%isma,j,k,m) = 0._rk
           enddo
        end do
     end do
  endif

#endif /* NOTRA */

#ifndef NOTRA
 if (config%p_ntr .ne.0) then


     if (use_mpi) then

#if defined(OPENMP_TRANSPORT) && defined(OPEN_MP_3D)
!$omp parallel do private(j,k)
#endif

        do k = nzmoms, nzmome

#if defined(OPENMP_TRANSPORT) && defined(OPEN_MP_2D)
!$omp parallel do private(j)
#endif

           do j = nymoms, nymome

              if (config%use_multid_collapse) then
                 ralag_3d   (-1:config%imaxp+1,j,k) = 0._rk
                 ralagold_3d(-1:config%imaxp+1,j,k) = 0._rk
              else
                 ralag_1d   (-1:config%imaxp+1) = 0._rk
                 ralagold_1d(-1:config%imaxp+1) = 0._rk
              endif
              xlag (-1:config%imaxp*2+3,1:config%iemax,1:config%isma,j,k) = 0._rk
              wlag (-1:config%imaxp*2+3,1:config%iemax,1:config%isma,j,k) = 0._rk
              xlags(-1:config%imaxp*2+3,1:config%iemax,1:config%isma,j,k) = 0._rk
              wlags(-1:config%imaxp*2+3,1:config%iemax,1:config%isma,j,k) = 0._rk
              beythq_2d(-1:config%imaxp+1,j,k) = 0._rk
              beyqth_2d( 0:config%imaxp+1,j,k) = 0._rk
              
              selag(-4:config%imaxp+1+4,j,k) = 0._rk
              smlag(-4:config%imaxp+1+4,j,k) = 0._rk
              sylag(-4:config%imaxp+1+4,j,k,1) = 0._rk
              sylag(-4:config%imaxp+1+4,j,k,:) = 0._rk
              
              elag  (-1:config%imaxp*2+3,1:config%iemax,1:2*config%isma  ,j,k)         = 0._rk
              riplag(-1:config%imaxp+1,config%cmin:config%imaxp+1,1:config%iemax,1:config%isma,j,k)  = 0._rk
              rimlag(-1:config%imaxp+1,config%cmin:config%imaxp+1,1:config%iemax,1:config%isma,j,k)  = 0._rk
              riplags(-1:config%imaxp+1,config%cmin:config%imaxp+1,1:config%iemax,1:config%isma,j,k) = 0._rk
              rimlags(-1:config%imaxp+1,config%cmin:config%imaxp+1,1:config%iemax,1:config%isma,j,k) = 0._rk
              
           enddo
        enddo

     else ! no mpi used

#if defined(OPENMP_TRANSPORT) && defined(OPEN_MP_3D)
!$omp parallel do private(j,k)
#endif

        do k = nzmoms, nzmome

#if defined(OPENMP_TRANSPORT) && defined(OPEN_MP_2D)
!$omp parallel do private(j)
#endif

           do j = nymoms, nymome

              beythq_2d(-1:config%imaxp+1,j,k) = 0._rk
              beyqth_2d( 0:config%imaxp+1,j,k) = 0._rk
              
              selag(-4:config%imaxp+1+4,j,k) = 0._rk
              smlag(-4:config%imaxp+1+4,j,k) = 0._rk
              sylag(-4:config%imaxp+1+4,j,k,1) = 0._rk
              sylag(-4:config%imaxp+1+4,j,k,:) = 0._rk

           enddo
        enddo
#if defined(OPENMP_TRANSPORT) && defined(OPEN_MP_3D)
!$omp parallel do private(j,k)
#endif
        do k = nzmoms, nzmome

#if defined(OPENMP_TRANSPORT) && defined(OPEN_MP_2D)
!$omp parallel do private(j)
#endif
           do j = config%nystrt, config%nytra
              elag  (-1:config%imaxp*2+3,1:config%iemax,1:2*config%isma  ,j,k)         = 0._rk
              riplag(-1:config%imaxp+1,config%cmin:config%imaxp+1,1:config%iemax,1:config%isma,j,k)  = 0._rk
              rimlag(-1:config%imaxp+1,config%cmin:config%imaxp+1,1:config%iemax,1:config%isma,j,k)  = 0._rk

              riplags(-1:config%imaxp+1,config%cmin:config%imaxp+1,1:config%iemax,1:config%isma,j,k) = 0._rk
              rimlags(-1:config%imaxp+1,config%cmin:config%imaxp+1,1:config%iemax,1:config%isma,j,k) = 0._rk
              
           enddo

#if defined(OPENMP_TRANSPORT) && defined(OPEN_MP_2D)
!$omp parallel do private(j)
#endif
           do j = config%nystrt, config%nymom

              if (config%use_multid_collapse) then
                 ralag_3d   (-1:config%imaxp+1,j,k) = 0._rk
                 ralagold_3d(-1:config%imaxp+1,j,k) = 0._rk
              else
                 ralag_1d   (-1:config%imaxp+1) = 0._rk
                 ralagold_1d(-1:config%imaxp+1) = 0._rk
              endif ! multid_collapse
                 
              xlag (-1:config%imaxp*2+3,1:config%iemax,1:config%isma,j,k) = 0._rk
              wlag (-1:config%imaxp*2+3,1:config%iemax,1:config%isma,j,k) = 0._rk

              xlags(-1:config%imaxp*2+3,1:config%iemax,1:config%isma,j,k) = 0._rk
              wlags(-1:config%imaxp*2+3,1:config%iemax,1:config%isma,j,k) = 0._rk

           enddo
        enddo
     endif

  endif ! p_ntr
#endif /* NOTRA */

!#endif /* SGI */

#ifdef CHECK_THREAD_AFFINITY
  call check_thread_affinity
#endif



end subroutine initialize_arrays

end module initialize_array_mod
