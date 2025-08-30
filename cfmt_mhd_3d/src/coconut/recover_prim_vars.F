#ifdef ONEMG_EOS
#define SPIKE_DETEC
#else
#undef SPIKE_DETEC
#endif /* ONEMG_EOS */


c     ==================================================================

      MODULE recover_prim_vars

c     ==================================================================

      IMPLICIT NONE

      CONTAINS



#ifdef CFC_TRANSPORT2
c     ==================================================================

      SUBROUTINE recover_primitives

c     ==================================================================

c     Relativistic version

c     ------------------------------------------------------------------

      USE precision

      USE size_cfc

      USE parameters_cfc
      USE grid_cfc
      USE conserved_cfc
      USE hydro_primitives_cfc
      USE metric_cfc
      USE perm_aux_cfc
      USE output_hydro
      USE eos_sn2
#ifdef ONEMG_EOS
      USE burn_lim, ONLY: thigh2
#endif

      USE abort

      use configure

      IMPLICIT NONE
      
      real(kind=rk) :: p_error(m), p_old(m)
      real(kind=rk) :: f_p, df_dp
      real(kind=rk) :: eos_self(2), eos_children(2)
      real(kind=rk) :: den_tmp(m),ei_tmp(m),
     &     den_wrk(m),ei_wrk(m),tem_wrk(m),pre_wrk(m),
     &     s_wrk(m),gam_wrk(m),xnu_wrk(m,config%qn),kap_wrk(m)
      logical :: lconv(m), use_ls_eos(m), use_ls_wrk(m)

! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=rk) ::  ccu(m), cce(m), ccn(m), ccp(m)
      real(kind=rk) :: dmy1(m), dmy2(m,2) !dummys
! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CFC_MHD
      real(kind=rk) :: b0, b1, b2, b3, bsq, sqgi, bivi
#endif /* CFC_MHD */

      integer(kind=ik) :: i, j, k, jk, i_nuc, kt, recovery_iteration
      logical :: ler, converged

#ifdef SPIKE_DETEC
      real(kind=rk), dimension(m,10) :: temk_tmp
      real(kind=rk) :: medtem
      real(kind=rk), parameter :: dT = 1.15_rk
      integer(kind=ik) :: iplus,iminus
      logical, dimension(m) :: ave_tem
#endif /* SPIKE_DETEC */

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,p_error,p_old,f_p,df_dp,den_tmp,ei_tmp,
!$OMP&         den_wrk,ei_wrk,tem_wrk,pre_wrk,s_wrk,gam_wrk,
!$OMP&         xnu_wrk,kap_wrk,lconv,use_ls_eos,use_ls_wrk,
!$OMP&         ccu,cce,ccn,ccp,dmy1,dmy2,i_nuc,kt,ler,
#ifdef SPIKE_DETEC
!$OMP&         temk_tmp,medtem,iplus,iminus,ave_tem,
#endif
!$OMP&         recovery_iteration)
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         
         do i_nuc = 1, config%qn
            do i = 1, m
               xnnu (i, j, k, i_nuc) =
     &              d_cap_xnu_hat (i, j, k, i_nuc)/
     &              d_cap_hat     (i, j, k)
            end do
         end do

#if defined(CFC_MHD) && defined(PASSIVE_MHD)
c     In this case, we need not identify b^i with the field in
c     the comoving frame and can just set b^i=B^i=\hat{B}^i \Phi^{-6}
         do i = 1, m
            b_1 (i, j, k) = b_cap_1_hat (i, j, k) *
     &           sqrt_gamma_inv (i, j, k)
            b_2 (i, j, k) = b_cap_2_hat (i, j, k) *
     &           sqrt_gamma_inv (i, j, k)
            b_3 (i, j, k) = b_cap_3_hat (i, j, k) *
     &           sqrt_gamma_inv (i, j, k)
         end do
#endif /* CFC_MHD && PASSIVE_MHD */

#if 0
         if (any (xnnu(1:m,j,k,:).lt. -1.0e-13_rk)) then
            print *,'recovery error (negative mass fractions)',j,k
            do i = 1, m
               if (any(xnnu(i,j,k,:).lt.-1.0e-13_rk)) then
                  print *,i,xnnu(i,j,k,:)
               end if
            end do
         end if
         do i = 1, m
            if (abs(sum(xnnu(i,j,k,1:config%qn-1))-1.0_rk).gt. 
     &           1e-12_rk) then
               print *,'recovery error (sum)',i,j,k
               print *,i,sum(xnnu(i,j,k,1:config%qn-1))-1.0_rk
            end if
         end do
#endif
         
         lconv(1:m)=.false.
         p_error(1:m) = 1.0_rk

         rec_cyc: do recovery_iteration = 1,100
ccc   initialize p_old
            p_old(1:m)= p(1:m,j,k)
            do i = 1, m

               if (.not.lconv(i)) then
ccc   do Newton-Raphson iteration

#if defined(CFC_MHD) && !defined(PASSIVE_MHD)
                  sqgi = sqrt_gamma_inv (i, j, k)
                  bivi = sqgi / phi (i, j, k) ** 2 *
     &                 (b_cap_1_hat (i, j, k) * v_1 (i, j, k) +
     &                  b_cap_2_hat (i, j, k) * v_2 (i, j, k) *
     &                  r_inv (i) +
     &                  b_cap_3_hat (i, j, k) * v_3 (i, j, k) *
     &                  r_inv (i) * sinus_theta_inv (j))
                  b0 = bivi * w (i, j, k) / alpha (i, j, k)
c covariant components of b_i within iteration
                  b_1 (i, j, k) = b_cap_1_hat (i, j, k) /
     &                 phi (i, j, k) ** 2 /
     &                 w (i, j, k) + w (i, j, k) * bivi *
     &                 (v_1 (i, j, k) - beta_up_1 (i, j, k) *
     &                 g_11 (i, j, k) / alpha (i, j, k))
                  b_2 (i, j, k) = b_cap_2_hat (i, j, k) /
     &                 phi (i, j, k) ** 2 /
     &                 w (i, j, k) + w (i, j, k) * bivi *
     &                 (v_2 (i, j, k) - beta_up_2 (i, j, k) *
     &                 g_22 (i, j, k) / alpha (i, j, k))
                  b_3 (i, j, k) = b_cap_3_hat (i, j, k)  /
     &                 phi (i, j, k) ** 2 /
     &                 w (i, j, k) + w (i, j, k) * bivi *
     &                 (v_3 (i, j, k) - beta_up_3 (i, j, k) *
     &                 g_33 (i, j, k) / alpha (i, j, k))
                  bsq = (b_cap_1_hat (i, j, k) ** 2 +
     &                   b_cap_2_hat (i, j, k) ** 2 +
     &                   b_cap_3_hat (i, j, k) ** 3) *
     &                   (sqgi / w (i, j, k)) ** 2 +
     &                   bivi ** 2

                  v_1 (i, j, k) =
     &                 (s_1_hat (i, j, k) + alpha (i, j, k) *
     &                 b0 * b_1 (i, j, k)) / (tau_hat (i, j, k) + 
     &                 d_cap_hat (i, j, k) + (p_old (i) + bsq) *
     &                 sqrt_gamma (i, j, k))
                  v_2 (i, j, k) =
     &                 (s_2_hat (i, j, k) + alpha (i, j, k) *
     &                 b0 * b_2 (i, j, k)) / (tau_hat (i, j, k) +
     &                 d_cap_hat (i, j, k) + (p_old (i) + bsq)*
     &                 sqrt_gamma (i, j, k))
                  v_3 (i, j, k) =
     &                 (s_3_hat (i, j, k) + alpha (i, j, k) *
     &                 b0 * b_3 (i, j, k)) / (tau_hat (i, j, k) +
     &                 d_cap_hat (i, j, k) + (p_old (i) + bsq) *
     &                 sqrt_gamma (i, j, k))
                  
                  v_squared (i, j, k) =
     &                 g_up_11 (i, j, k) * v_1 (i, j, k) ** 2 +
     &                 g_up_22 (i, j, k) * v_2 (i, j, k) ** 2 +
     &                 g_up_33 (i, j, k) * v_3 (i, j, k) ** 2
                  
                  w (i, j, k) =
     &                 1.0_rk / sqrt (1.0_rk - v_squared (i, j, k))
                  
                  rho (i, j, k) = d_cap_hat (i, j, k) /
     &                 (w (i, j, k) * sqrt_gamma (i, j, k))
                  
                  eps (i, j, k) = (tau_hat (i, j, k) +
     &                 (p_old (i) * (1.0_rk - w (i, j, k) ** 2) +
     &                 bsq * (0.5_rk - w (i, j, k) ** 2)) *
     &                 sqrt_gamma (i, j, k) +
     &                 d_cap_hat (i, j, k) *
     &                 (1.0_rk - w (i, j, k)) +
     &                 alpha (i, j, k) * b0 ** 2) /
     &                 (d_cap_hat (i, j, k) * w (i, j, k))

#else
                  v_1 (i, j, k) =
     &                 s_1_hat (i, j, k) / (tau_hat (i, j, k) + 
     &                 d_cap_hat (i, j, k) + p_old (i) *
     &                 sqrt_gamma (i, j, k))
                  v_2 (i, j, k) =
     &                 s_2_hat (i, j, k) / (tau_hat (i, j, k) +
     &                 d_cap_hat (i, j, k) + p_old (i) *
     &                 sqrt_gamma (i, j, k))
                  v_3 (i, j, k) =
     &                 s_3_hat (i, j, k) / (tau_hat (i, j, k) +
     &                 d_cap_hat (i, j, k) + p_old (i) *
     &                 sqrt_gamma (i, j, k))
                  
                  v_squared (i, j, k) =
     &                 g_up_11 (i, j, k) * v_1 (i, j, k) ** 2 +
     &                 g_up_22 (i, j, k) * v_2 (i, j, k) ** 2 +
     &                 g_up_33 (i, j, k) * v_3 (i, j, k) ** 2
                  
                  w (i, j, k) =
     &                 1.0_rk / sqrt (1.0_rk - v_squared (i, j, k))
                  
                  rho (i, j, k) = d_cap_hat (i, j, k) /
     &                 (w (i, j, k) * sqrt_gamma (i, j, k))
                  
                  eps (i, j, k) = (tau_hat (i, j, k) + p_old (i) *
     &                 sqrt_gamma (i, j, k) *
     &                 (1.0_rk - w (i, j, k) ** 2) +
     &                 d_cap_hat (i, j, k) *
     &                 (1.0_rk - w (i, j, k))) /
     &                 (d_cap_hat (i, j, k) * w (i, j, k))
#endif /* CFC_MHD */
                  
                  
                  den_tmp (i) = rho(i, j, k) * rho_geom_factor_inv
                  ei_tmp  (i) = eps(i, j, k) * p_geom_factor_inv *
     &                 rho(i, j, k)
               end if
            end do
            
            if (recovery_iteration.eq.1) then
c     Assignment of cells to high- and low-density EOS outside iteration
c     to avoid trouble at the transition density
               do i = 1, m
#ifdef ONEMG_EOS
                  use_ls_eos(i)=den_tmp(i).gt.lsrolo2 .or.
     &                 (den_tmp(i).ge.lsrolo .and.
     &                 t(i,j,k).ge.thigh2)
#else
                  use_ls_eos(i)=den_tmp(i).gt.lsrolo
#endif
               end do
            end if

c     --------------------------------------------------------------------                  
c     Gather:
            kt=0
            do i=1,m
               if (.not.lconv(i)) then
                  kt=kt+1
                  den_wrk(kt)=den_tmp(i)
                  ei_wrk (kt)=ei_tmp (i)
                  tem_wrk(kt)=t(i,j,k)
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
                  xnu_wrk(kt,1:config%qn)=xnnu(i,j,k,1:config%qn)
                  use_ls_wrk(kt)=use_ls_eos(i)
               end if
            end do
            
c     Call EOS:
            call eos(den_wrk(1:kt),tem_wrk(1:kt),xnu_wrk(1:kt,:),
     &           dmy1(1:kt),dmy2(1:kt,:),ei_wrk(1:kt),pre_wrk(1:kt),
     &           gam_wrk(1:kt),s_wrk(1:kt),ccu(1:kt),cce(1:kt),
     &           ccn(1:kt),ccp(1:kt),eos_self, eos_children,mode=3,
     &           nsemode=0,ler=ler,
     &           use_ls_eos=use_ls_wrk(1:kt))
            if (ler) then
               raise_abort("EoS call failed!")
            end if
            
c     Scatter:
            kt=0
            do i=1,m
               if (.not.lconv(i)) then
                  kt=kt+1
                  t      (i,j,k)   =tem_wrk(kt)
                  p      (i,j,k)   =pre_wrk(kt)*p_geom_factor
                  entropy(i,j,k)   =s_wrk(kt)
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
                  xnnu (i,j,k,1:config%qn)=xnu_wrk(kt,1:config%qn)
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
                  d_cap_xnu_hat(i,j,k,1:config%qn)=
     &                 d_cap_hat(i,j,k)*xnu_wrk(kt,1:config%qn)
                  cpot (i,j,k,1)  = ccu(kt)
                  cpot (i,j,k,2)  = cce(kt)
                  cpot (i,j,k,3)  = ccn(kt)
                  cpot (i,j,k,4)  = ccp(kt)
                  gamm (i,j,k)     =gam_wrk(kt)
                  h      (i,j,k)   = 1.0_rk+eps(i,j,k)+
     &                 p(i,j,k)/rho(i,j,k)
                  c_sound_squared(i,j,k) =
     &                 gam_wrk(kt)*p(i,j,k)/
     &                 (rho(i,j,k)*h(i,j,k))
               end if
            end do
c     --------------------------------------------------------------------
                  
            do i=1,m
               if (.not.lconv(i)) then
                  f_p = p (i, j, k) - p_old (i)
                  
                  df_dp = v_squared (i, j, k) *
     &                 c_sound_squared (i, j, k) - 1.0_rk
                  
                  p (i, j, k) = p_old (i) - f_p / df_dp
                  
                  p_error(i) = abs (p (i, j, k) /
     &                 p_old (i) - 1.0_rk)
                  
                  if (p_error(i).le.
     &                 conversion_tolerance) then
                     lconv(i) = .true.                     
                     p_old (i) = p (i, j, k)
                  end if
                  
               end if
            enddo

            if (.not.any(.not.lconv(1:m))) exit rec_cyc
            
            if (recovery_iteration.ge.30.and.
     &           .not.any(p_error(1:m).ge.1.0e-4_rk)) then
               print*
               print*, 'Nonfatal problem in:'
               print*, '   Recovery of primitive quantities'
               print*
               print*, 'Diagnosis:'
               print*, '   No high precision reached in ' //
     &              'recovery iteration'
               print*
               print*, 'Action:'
               print*, '   Continuing'
               print*
               do i=1,m
                  if (.not.lconv(i))
     &                 print *,i,rho(i,j,k)*rho_geom_factor_inv,
     &                 w(i,j,k)
               end do
               lconv(:)  = .true.
            endif
         end do rec_cyc
         

#if defined (CFC_MHD)
c     Convert b_i to correct basis again
         do i = 1, m
            b_1 (i, j, k) = b_1 (i, j ,k) * g_up_11 (i, j, k)
            b_2 (i, j, k) = b_2 (i, j ,k) * g_up_22 (i, j, k)
            b_3 (i, j, k) = b_3 (i, j ,k) * g_up_33 (i, j, k)
         end do
#endif  /* CFC_MHD */


#ifdef SPIKE_DETEC
c     ------------------------------------------------------------------
c     temperature spikes:
         ave_tem(:) = .false.
         temk_tmp(:,:) = spread(t(1:m,j,k),dim=2,ncopies=10)
         do iminus = -2, -1
            do iplus = 1, 2
               do i = 6, m-5
                  if ((t(i,j,k) .gt. temk_tmp(i+iminus,11+iminus)*dT 
     &                 .and. t(i,j,k) .gt. temk_tmp(i+iplus,iplus)*dT)
     &                 .or. 
     &                 (temk_tmp(i+iminus,11+iminus) .gt. t(i,j,k)*dT 
     &                 .and. temk_tmp(i+iplus,iplus) .gt. t(i,j,k)*dT))
     &                 then
                     ave_tem(i) = .true.
                  end if
               end do
            end do
         end do
         
         do i = 6, m-5
            if (ave_tem(i)) then
!               print *,'Smoothing',i,t(i-2:i+2,j,k)
               medtem = (temk_tmp(i-5,1)+temk_tmp(i-4,2)+
     &              temk_tmp(i-3,3)+temk_tmp(i-2,4)+temk_tmp(i-1,5)+
     &              temk_tmp(i+1,6)+temk_tmp(i+2,7)+temk_tmp(i+3,8)+
     &              temk_tmp(i+4,9)+temk_tmp(i+5,10))/10.0d0
               t(i,j,k) = medtem
            end if
         end do
#endif



      enddo                     !j,k-loop

      call check_for_floating_exceptions_hydro
      
      return

      END SUBROUTINE recover_primitives

c     ==================================================================


#else

c     ==================================================================

      SUBROUTINE recover_primitives

c     ==================================================================

c     Newtonian version

c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
      USE conserved_cfc
      USE hydro_primitives_cfc
      USE metric_cfc
      USE perm_aux_cfc
      USE eos_sn2

      USE configure
      USE phycon!debug
      USE abort

      IMPLICIT NONE

      real(kind=rk) :: den_tmp(m),ei_tmp(m),xnu_tmp(m,config%qn)
      logical :: lconv(m) 

      real(kind=rk) :: ccu(m), cce(m), ccn(m), ccp(m)
      real(kind=rk) :: dmy1(m), dmy2(m,2) !dummys
      real(kind=rk) :: bsq
      real(kind=rk) :: eos_self(2), eos_children(2)
      
      integer(kind=ik) :: i, j, k, jk, i_nuc
      logical :: ler

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,i_nuc,den_tmp,ei_tmp,xnu_tmp,ler)
!$OMP& SCHEDULE(static)
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m
ccc   get density and velocities
            rho (i, j, k) = d_cap_hat (i, j, k)
            v_1 (i, j, k) = s_1_hat (i, j, k) / rho (i, j, k)
            v_2 (i, j, k) = s_2_hat (i, j, k) / rho (i, j, k)
            v_3 (i, j, k) = s_3_hat (i, j, k) / rho (i, j, k)
            
ccc   calculate v_squared (needed for e.o.s.)
            v_squared (i, j, k) =
     &           g_up_11 (i, j, k) * v_1 (i, j, k) ** 2 +
     &           g_up_22 (i, j, k) * v_2 (i, j, k) ** 2 +
     &           g_up_33 (i, j, k) * v_3 (i, j, k) ** 2
            
ccc   get internal energy density from tau
            eps (i, j, k) = tau_hat (i, j, k) / rho (i, j, k) -
     &           0.5_rk * v_squared (i, j, k)

#ifdef CFC_MHD
            b_1 (i, j, k) = b_cap_1_hat (i, j, k)
            b_2 (i, j, k) = b_cap_2_hat (i, j, k)
            b_3 (i, j, k) = b_cap_3_hat (i, j, k)
            bsq = b_1 (i, j, k) ** 2 +
     &           b_2 (i, j, k) ** 2 + 
     &           b_3 (i, j, k) ** 2
#ifndef CONSTRAINED
     &           + psi (i, j, k) ** 2
#endif
 
#ifndef PASSIVE_MHD
            eps (i, j, k) = eps (i, j, k) - 0.5_rk *
     &           bsq / rho (i, j, k)
c            IF (eps (i,j,k) * rho (i, j, k) .LT.
            IF (p (i, j, k) .LT.
     &           MAX (0.05_rk * (0.5_rk * bsq +
     &           0.5_rk * rho (i, j, k) * v_squared (i, j, k)),
     &           0.5_rk * 0.5_rk *  bsq)) then
c               print *,'internal energy equation',i,j,k,
c     &              eps (i, j, k), eps_aux (i, j, k) / rho (i, j, k)
               eps (i, j, k) = max(eps (i, j, k),
     &              eps_aux (i, j, k) / rho (i, j, k))
ccc               eps (i, j, k) = 
ccc     &              eps_aux (i, j, k) / rho (i, j, k)
            END IF
#endif /* PASSIVE MHD */
#endif /* CFC_MHD */

         end do
         
         do i_nuc = 1, config%qn
            do i = 1, m
               xnnu (i, j, k, i_nuc) =
     &              d_cap_xnu_hat (i, j, k, i_nuc)/
     &              d_cap_hat     (i, j, k)
            end do
         end do
         
c     copy and convert from geoemtrical units
         den_tmp(1:m) = rho(1:m,j,k)*rho_geom_factor_inv
         ei_tmp(1:m)  = eps(1:m,j,k)*p_geom_factor_inv*
     &        rho(1:m,j,k)
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
         xnu_tmp(1:m,1:config%qn)    =xnnu(1:m,j,k,1:config%qn)
         
c     call EOS
         call eos(den_tmp(1:m),t(1:m,j,k),xnu_tmp(1:m,:),
     &        dmy1(1:m),dmy2(1:m,:),ei_tmp(1:m),
     &        p(1:m,j,k),gamm(1:m,j,k),entropy(1:m,j,k),
     &        ccu(1:m),cce(1:m),ccn(1:m),ccp(1:m),
     &        eos_self, eos_children,
     &        mode=3,nsemode=0,ler=ler)
         if (ler) then
            raise_abort("EoS call failed!") 
         end if

c     copy back and convert back to geomtrical units
         p    (1:m,j,k)     = p(1:m,j,k)*p_geom_factor
         xnnu (1:m,j,k,1:config%qn)= xnu_tmp(1:m,1:config%qn)
         do i_nuc = 1, config%qn
            d_cap_xnu_hat(1:m,j,k,i_nuc)=
     &           d_cap_hat(1:m,j,k)*xnu_tmp(1:m,i_nuc)
         end do
         h    (1:m,j,k)     =
     &        eps(1:m,j,k)+
     &        p(1:m,j,k)/rho(1:m,j,k)+
     &        0.5_rk*v_squared(1:m,j,k)
         c_sound_squared(1:m,j,k) =
     &        gamm(1:m,j,k)*p(1:m,j,k)/rho(1:m,j,k)
         w    (1:m,j,k)     = 1.0_rk
         cpot (1:m,j,k,1)  = ccu(1:m)
         cpot (1:m,j,k,2)  = cce(1:m)
         cpot (1:m,j,k,3)  = ccn(1:m)
         cpot (1:m,j,k,4)  = ccp(1:m)
#if defined(CFC_MHD) &&         !defined(PASSIVE_MHD)
         do i = 1, m
c            if (eps (i,j,k) * rho (i, j, k) .lt.
c     &           3e-1_rk * bsq)
c     &           print *,'internal energy equation',i,j,k,
c     &           eps (i, j, k), eps_aux (i, j, k) / rho (i, j, k),
c     &           entropy (i, j, k)
c     if (p (i,j,k) .lt.
c     &        1e-2_rk * tau_hat (i, j, k)) then
c     Recompute tau_hat in case we solved the internal energy equation 
            tau_hat   (i, j, k) = rho (i, j, k) *
     &           (eps (i, j, k) + 0.5_rk * v_squared (i,j,k)) +
     &           0.5_rk * 
     &           (b_1 (i, j, k) ** 2 + b_2 (i, j, k) ** 2 +
     &           b_3 (i, j, k) ** 2)
#ifndef CONSTRAINED
     &           + 0.5_rk * psi (i, j, k) ** 2
#endif
c            endif
         enddo
#endif /* CFC_MHD and not PASSIVE_MHD */            
            
      enddo
      
      call check_for_floating_exceptions_hydro
      

      return

      END SUBROUTINE recover_primitives

c     ==================================================================

#endif /*CFC_TRANSPORT2*/


c     ==================================================================

      SUBROUTINE check_for_floating_exceptions_hydro

c     ==================================================================

      USE precision

      USE size_cfc
      USE hydro_primitives_cfc
      USE parameters_cfc

      USE abort

      IMPLICIT NONE
      
      integer(kind=ik) :: i, j, k, jk

      logical :: ler

      ler=.false.

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
!$OMP& SCHEDULE(static)
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m

c     check for NANs
            if (rho (i, j, k).ne.rho (i, j, k)) ler=.true.
            if (v_1 (i, j, k).ne.v_1 (i, j, k)) ler=.true.
            if (v_2 (i, j, k).ne.v_2 (i, j, k)) ler=.true.
            if (v_3 (i, j, k).ne.v_3 (i, j, k)) ler=.true.
            if (eps (i, j, k).ne.eps (i, j, k)) ler=.true.
            
c     check for INFs
            if (2.0_rk * rho (i, j, k).eq.rho (i, j, k).and.
     &           rho (i, j, k).ne.0.0_rk) ler=.true.
            if (2.0_rk * v_1 (i, j, k).eq.v_1 (i, j, k).and.
     &           v_1 (i, j, k).ne.0.0_rk) ler=.true.
            if (2.0_rk * v_2 (i, j, k).eq.v_2 (i, j, k).and.
     &           v_2 (i, j, k).ne.0.0_rk) ler=.true.
            if (2.0_rk * v_3 (i, j, k).eq.v_3 (i, j, k).and.
     &           v_3 (i, j, k).ne.0.0_rk) ler=.true.
            if (2.0_rk * eps (i, j, k).eq.eps (i, j, k).and.
     &           eps (i, j, k).ne.0.0_rk) ler=.true.

         enddo
      enddo
      
      if (.not.ler) goto 111

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
      
      print*
      print*, 'Fatal problem in:'
      print*, '   Recovery of primitive quantities'
      print*
      print*, 'Diagnosis:'
      print*, '   Floating point exception in hydro quantities'
      print*
      print*, 'Action:'
      print*, '   Stopping'
      print*
      raise_abort ("Floating point exception")
      
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif      

 111  continue
      
      return

      END SUBROUTINE check_for_floating_exceptions_hydro

c     ==================================================================



      END MODULE recover_prim_vars

c     ==================================================================
