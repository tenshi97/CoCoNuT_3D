c     ==================================================================

      MODULE gr_sources

c     ==================================================================

      IMPLICIT NONE

      CONTAINS



c     ==================================================================

      SUBROUTINE calculate_sources

c     ==================================================================

      USE precision

      USE size_cfc
      USE grid_cfc
      USE conserved_cfc
      USE hydro_primitives_cfc
      USE perm_aux_cfc
      USE sources_cfc
      USE metric_cfc
      USE parameters_cfc
      USE nutra_cfc

      use configure

      IMPLICIT NONE

      real (kind=rk) :: em_tensor (0 : 3, 0 : 3)
      real (kind=rk) :: em_tensor_excrv_contracted

      real (kind=rk) :: alpha_inv, alpha_2_inv

      real (kind=rk) :: v_up_1, v_up_2, v_up_3, v_up_1_hat, v_up_2_hat,
     &     v_up_3_hat, s_up_1_hat, s_up_2_hat, s_up_3_hat

      real (kind=rk) :: rho_h_ww
#ifdef CFC_MHD
      real (kind=rk) :: bsq, p_sav, delta_t_inv
#endif

      real (kind=rk) :: v_i_hat_times_dln_sqrt_gamma_di

      integer (kind=ik) :: i, j, k, jk, l, i_nuc

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,alpha_inv,alpha_2_inv,v_up_1,v_up_2,v_up_3,
!$OMP& v_up_1_hat,v_up_2_hat,v_up_3_hat,
!$OMP& v_i_hat_times_dln_sqrt_gamma_di,
!$OMP& rho_h_ww,s_up_1_hat,s_up_2_hat,s_up_3_hat)
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

#ifdef CFC_MHD
         delta_t_inv = 1.0_rk / delta_t !needed by every OpenMP thread
#endif
         
         do i = 1, m

            alpha_inv = 1.0_rk / alpha (i, j, k)
            alpha_2_inv = alpha_inv ** 2
            
            v_up_1 = g_up_11 (i, j, k) * v_1 (i, j, k)
            v_up_2 = g_up_22 (i, j, k) * v_2 (i, j, k)
            v_up_3 = g_up_33 (i, j, k) * v_3 (i, j, k)
            
            v_up_1_hat = v_up_1 - beta_up_1 (i, j, k) * alpha_inv
            v_up_2_hat = v_up_2 - beta_up_2 (i, j, k) * alpha_inv
            v_up_3_hat = v_up_3 - beta_up_3 (i, j, k) * alpha_inv
            
              

c     calculate energy momentum tensor
c     ==================================================================

#ifdef CFC_TRANSPORT2
            rho_h_ww = rho (i, j, k) * h (i, j, k) *
     &           w (i, j, k) ** 2
#else /* CFC_TRANSPORT2 */
            rho_h_ww = rho (i, j, k)
#endif /* CFC_TRANSPORT2 */
            
#if defined(CFC_MHD) &&         !defined(PASSIVE_MHD)
c     Eq. (22) in Moesta et al. (2013), remember that b_i are
c     contravariant components
            bsq = (b_1 (i, j, k) ** 2 + b_2 (i, j, k) ** 2 + 
     &           b_3 (i, j, k) ** 2) * phi (i, j, k) ** 4 /
     &           w (i, j, k) ** 2
#ifdef CFC_TRANSPORT2
     &           + (b_1 (i, j, k) * v_up_1 +
     &           r_inv (i) * (b_2 (i, j, k) * v_up_2 +
     &           b_3 (i, j, k) * sinus_theta_inv (j) * v_up_3)) ** 2
            rho_h_ww = rho_h_ww + bsq * w (i, j, k) ** 2
#endif /* CFC_TRANSPORT2*/
            p_sav = p (i, j, k)
            p (i, j, k) = p (i, j, k) + 0.5_rk * bsq
#endif /*CFC_MHD && !PASSIVE_MHD */
            
            s_up_1_hat = rho_h_ww * v_up_1_hat
            s_up_2_hat = rho_h_ww * v_up_2_hat
            s_up_3_hat = rho_h_ww * v_up_3_hat
            
            em_tensor (0, 0) =
     &           alpha_2_inv * (rho_h_ww - p (i, j, k))
            em_tensor (0, 1) = alpha_inv * (s_up_1_hat +
     &           alpha_inv * beta_up_1 (i, j, k) * p (i, j, k))
            em_tensor (0, 2) = alpha_inv * (s_up_2_hat +
     &           alpha_inv * beta_up_2 (i, j, k) * p (i, j, k))
            em_tensor (0, 3) = alpha_inv * (s_up_3_hat +
     &           alpha_inv * beta_up_3 (i, j, k) * p (i, j, k))
            em_tensor (1, 1) = s_up_1_hat * v_up_1_hat +
     &           (g_up_11 (i, j, k) - beta_up_1 (i, j, k) ** 2 *
     &           alpha_2_inv) * p (i, j, k)
            em_tensor (1, 2) = s_up_1_hat * v_up_2_hat -
     &           beta_up_1 (i, j, k) * beta_up_2 (i, j, k) *
     &           alpha_2_inv * p (i, j, k)
            em_tensor (1, 3) = s_up_1_hat * v_up_3_hat -
     &           beta_up_1 (i, j, k) * beta_up_3 (i, j, k) *
     &           alpha_2_inv * p (i, j, k)
            em_tensor (2, 2) = s_up_2_hat * v_up_2_hat +
     &           (g_up_22 (i, j, k) - beta_up_2 (i, j, k) ** 2 *
     &           alpha_2_inv) * p (i, j, k)
            em_tensor (2, 3) = s_up_2_hat * v_up_3_hat -
     &           beta_up_2 (i, j, k) * beta_up_3 (i, j, k) *
     &           alpha_2_inv * p (i, j, k)
            em_tensor (3, 3) = s_up_3_hat * v_up_3_hat +
     &           (g_up_33 (i, j, k) - beta_up_3 (i, j, k) ** 2 *
     &           alpha_2_inv) * p (i, j, k)

#if defined(CFC_MHD) && !defined(PASSIVE_MHD)
c     Subtract magnetic pressure again.
            p (i, j, k) = p_sav
#ifdef CFC_TRANSPORT2
            stop 'GRMHD does not work yet.'
#else
            em_tensor (1, 1) = em_tensor (1, 1) -
     &           b_1 (i, j, k) ** 2 * g_up_11 (i, j, k) 
            em_tensor (1, 2) = em_tensor (1, 2) -
     &           b_1 (i, j, k) * b_2 (i, j, k) *
     &           sqrt (g_up_11 (i, j, k) * g_up_22 (i, j, k))
            em_tensor (1, 3) = em_tensor (1, 3) -
     &           b_1 (i, j, k) * b_3 (i, j, k) *
     &           sqrt (g_up_11 (i, j, k) * g_up_33 (i, j, k))
            em_tensor (2, 2) = em_tensor (2, 2) -
     &           b_2 (i, j, k) ** 2 * g_up_22 (i, j, k)
            em_tensor (2, 3) = em_tensor (2, 3) -
     &           b_2 (i, j, k) * b_3 (i, j, k) *
     &           sqrt (g_up_22 (i, j, k) * g_up_33 (i, j, k))
            em_tensor (3, 3) = em_tensor (3, 3) -
     &           b_3 (i, j, k) ** 2 * g_up_33 (i, j, k)
#endif /* CFC_TRANSPORT2*/            
#endif /* CFC_MHD and not PASSIVE_MHD */
            
            em_tensor (1, 0) = em_tensor (0, 1)
            em_tensor (2, 0) = em_tensor (0, 2)
            em_tensor (2, 1) = em_tensor (1, 2)
            em_tensor (3, 0) = em_tensor (0, 3)
            em_tensor (3, 1) = em_tensor (1, 3)
            em_tensor (3, 2) = em_tensor (2, 3)
            
            em_tensor_excrv_contracted =
     &           em_tensor (1, 1) *
     &           excrv (i, j, k, 1, 1) +
     &           em_tensor (2, 2) *
     &           excrv (i, j, k, 2, 2) +
     &           em_tensor (3, 3) *
     &           excrv (i, j, k, 3, 3) +
     &           em_tensor (1, 2) *
     &           excrv (i, j, k, 1, 2) +
     &           em_tensor (1, 3) *
     &           excrv (i, j, k, 1, 3) +
     &           em_tensor (2, 3) *
     &           excrv (i, j, k, 2, 3) +
     &           em_tensor (2, 1) *
     &           excrv (i, j, k, 2, 1) +
     &           em_tensor (3, 1) *
     &           excrv (i, j, k, 3, 1) +
     &           em_tensor (3, 2) *
     &           excrv (i, j, k, 3, 2)



c     calculate bare sources
c     ==================================================================
               
            s_total (i, j, k, 5) =
     &           em_tensor (0, 0) *
     &           (excrv_beta (i, j, k) -
     &           beta_up_1 (i, j, k) * dalpha_dr (i, j, k) -
     &           beta_up_2 (i, j, k) * dalpha_dtheta (i, j, k) -
     &           beta_up_3 (i, j, k) * dalpha_dvarphi (i, j, k)) +
     &           em_tensor (0, 1) * (2.0_rk * beta_up_1 (i, j, k) *
     &           excrv (i, j, k, 1, 1)) +
     &           em_tensor (0, 2) * (2.0_rk * beta_up_2 (i, j, k) *
     &           excrv (i, j, k, 1, 2)) +
     &           em_tensor (0, 3) * (2.0_rk * beta_up_3 (i, j, k) *
     &           excrv (i, j, k, 1, 3)) +
     &           em_tensor_excrv_contracted

c     ------------------------------------------------------------------

            s_total (i, j, k, 1) = 0.0_rk

#if MOMENTUM_EQUATION_MODE>0
c     ------------------------------------------------------------------
c     Recalculate energy momentum tensor (without pressure terms)
            
            em_tensor (0, 0) = alpha_2_inv * rho_h_ww
            em_tensor (0, 1) = alpha_inv * s_up_1_hat
            em_tensor (0, 2) = alpha_inv * s_up_2_hat
            em_tensor (0, 3) = alpha_inv * s_up_3_hat 
            

#if MOMENTUM_EQUATION_MODE==2
            em_tensor (1, 1) = s_up_1_hat * v_up_1_hat 
            em_tensor (1, 2) = s_up_1_hat * v_up_2_hat
            em_tensor (1, 3) = s_up_1_hat * v_up_3_hat
            em_tensor (2, 2) = s_up_2_hat * v_up_2_hat * rr (i)
            em_tensor (2, 3) = s_up_2_hat * v_up_3_hat
#if 0
c     Also eliminate geometric source terms from varphi-component
c     of velocity
            em_tensor (3, 3) = s_up_3_hat * v_up_3_hat * rr (i) *
     &           ssinus_theta (j)
#else
            em_tensor (3, 3) = s_up_3_hat * v_up_3_hat
#endif
#else
            em_tensor (1, 1) = s_up_1_hat * v_up_1_hat 
            em_tensor (1, 2) = s_up_1_hat * v_up_2_hat
            em_tensor (1, 3) = s_up_1_hat * v_up_3_hat
            em_tensor (2, 2) = s_up_2_hat * v_up_2_hat 
            em_tensor (2, 3) = s_up_2_hat * v_up_3_hat
            em_tensor (3, 3) = s_up_3_hat * v_up_3_hat 
#endif /* MOMENTUM_EQUATION_MODE==2 */
#endif /* MOMENTUM_EQUATION_MODE>0 */


#if defined(CFC_MHD) && !defined(PASSIVE_MHD)
#ifdef CFC_TRANSPORT2
            stop 'GRMHD does not work yet.'
#else
c     only the diagonal components are needed:
            em_tensor (1, 1) = em_tensor (1, 1) -
     &           b_1 (i, j, k) ** 2 * g_up_11 (i, j, k) 
            em_tensor (2, 2) = em_tensor (2, 2) -
#if MOMENTUM_EQUATION_MODE==2
     &           b_2 (i, j, k) ** 2 
#else            
     &           b_2 (i, j, k) ** 2 * g_up_22 (i, j, k)
#endif            
            em_tensor (3, 3) = em_tensor (3, 3) -
     &           b_3 (i, j, k) ** 2 * g_up_33 (i, j, k)
#endif /* CFC_TRANSPORT2*/            
#endif /* CFC_MHD and not PASSIVE_MHD */

            
c     ------------------------------------------------------------------
               
            s_total (i, j, k, 2) = 0.5_rk * (
     &           em_tensor (0, 0) * dg_00_dr (i, j, k) +
#if MOMENTUM_EQUATION_MODE==2
#if 0
     &           (em_tensor (1, 1) + em_tensor (2, 2) +
     &           em_tensor (3, 3)) *
     &           dg_11_dr (i, j, k)) +
#else
     &           (em_tensor (1, 1) + em_tensor (2, 2)) *
     &           dg_11_dr (i, j, k) +
     &           em_tensor (3, 3) * dg_33_dr (i, j, k)) +
#endif
#else
     &           em_tensor (1, 1) * dg_11_dr (i, j, k) +
     &           em_tensor (2, 2) * dg_22_dr (i, j, k) +
     &           em_tensor (3, 3) * dg_33_dr (i, j, k)) +
#endif /* MOMENTUM_EQUATION_MODE == 2 */
     &           em_tensor (0, 1) * dg_01_dr (i, j, k) +
     &           em_tensor (0, 2) * dg_02_dr (i, j, k) +
     &           em_tensor (0, 3) * dg_03_dr (i, j, k)
               
c     ------------------------------------------------------------------
            
            s_total (i, j, k, 3) = 0.5_rk * (
     &           em_tensor (0, 0) * dg_00_dtheta (i, j, k) +
#if MOMENTUM_EQUATION_MODE==2
#if 0
     &           (em_tensor (1, 1) + em_tensor (2, 2) +
     &           em_tensor (3, 3)) *
     &           dg_11_dtheta (i, j, k)) + 
#else
     &           (em_tensor (1, 1) + em_tensor (2, 2)) *
     &           dg_11_dtheta (i, j, k) +
     &           em_tensor (3, 3) * dg_33_dtheta (i, j, k)) +
#endif
#else
     &           em_tensor (1, 1) * dg_11_dtheta (i, j, k) +
     &           em_tensor (2, 2) * dg_22_dtheta (i, j, k) +
     &           em_tensor (3, 3) * dg_33_dtheta (i, j, k)) +
#endif /* MOMENTUM_EQUATION_MODE==2 */
            
     &           em_tensor (0, 1) * dg_01_dtheta (i, j, k) +
     &           em_tensor (0, 2) * dg_02_dtheta (i, j, k) +
     &           em_tensor (0, 3) * dg_03_dtheta (i, j, k)
            
c     ------------------------------------------------------------------

            if (config%nsdim .eq. 3) then
               s_total (i, j, k, 4) = 0.5_rk *
     &              (em_tensor (0, 0) * dg_00_dvarphi (i, j, k) +
     &              em_tensor (1, 1) * dg_11_dvarphi (i, j, k) +
     &              em_tensor (2, 2) * dg_22_dvarphi (i, j, k) +
     &              em_tensor (3, 3) * dg_33_dvarphi (i, j, k)) +
     &              em_tensor (0, 1) * dg_01_dvarphi (i, j, k) +
     &              em_tensor (0, 2) * dg_02_dvarphi (i, j, k) +
     &              em_tensor (0, 3) * dg_03_dvarphi (i, j, k)
            
c     ------------------------------------------------------------------
            endif

#ifdef CFC_MHD
c            s_total (i, j, k, 6) = - alpha (i, j, k) *
c     &           g_up_11 (i, j, k) *
c     &           (delta_r_amp_inv (i) *
c     &           (psi (i + 1, j, k) +
c     &           (aamp (i) - 1.0_rk) * psi (i, j ,k) -
c     &           aamp (i) * psi (i - 1, j, k))) !* 0.0_rk

c            s_total (i, j, k, 7) = - 0.5_rk * alpha (i, j, k) *
c     &           g_up_22 (i, j, k) *
c     &           (psi (i, j + 1, k) -  psi (i, j - 1 ,k)) *
c     &           delta_theta_inv !* 0.0_rk

#if DIMN<3
c            s_total (i, j, k, 8) = - 0.5_rk * alpha (i, j, k) *
c     &           sqrt_gamma (i, j, k) * g_up_33 (i, j, k) *
c     &           (psi (i, j, k + 1) -  psi (i, j ,k - 1)) *
c     &           delta_varphi_inv
#else
c            s_total (i, j, k, 8) = 0.0_rk
#endif /* DIMN==3 */            

#endif /* CFC_MHD */            
         enddo
      enddo


c     add Newtonian gravity potential in the Newtonian approximation
c     ==================================================================

#ifndef CFC_TRANSPORT2
c      call calculate_newtonian_potential
      call add_newtonian_sources
#endif


c     add correct factor for exact conservation option
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,l)
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m
            s_total (i, j, k, 2) =
     &           s_total (i, j, k, 2) * sqrt_gamma (i, j, k)
            s_total (i, j, k, 3) =
     &           s_total (i, j, k, 3) * sqrt_gamma (i, j, k)
            s_total (i, j, k, 4) =
     &           s_total (i, j, k, 4) * sqrt_gamma (i, j, k)
            s_total (i, j, k, 5) =
     &           s_total (i, j, k, 5) * sqrt_gamma (i, j, k)

#ifdef ENECONS
#ifdef CFC_TRANSPORT2
            s_total (i, j, k, 5) = alpha (i, j, k) *
     &           s_total (i, j, k, 5)
#else 
            s_total(i, j, k, 5)=0.0_rk
#endif /* CFC_TRANSPORT2 */
#endif /* ENECONS */

#ifdef CFC_MHD
c            s_total (i, j, k, 6) =
c     &           s_total (i, j, k, 6) * sqrt_gamma (i, j, k)
c            s_total (i, j, k, 7) =
c     &           s_total (i, j, k, 7) * sqrt_gamma (i, j, k)
c            s_total (i, j, k, 8) =
c     &           s_total (i, j, k, 8) * sqrt_gamma (i, j, k)


c     Energy source term from divergence cleaning
c     Still unclear what this term should be in the GR case,
c     especially with ENECONS
c            s_total (i, j, k, 5) = s_total (i, j, k, 5) +
c     &           (b_1 (i, j, k) * s_total (i, j, k, 6) +
c     &           (b_2 (i, j, k) * s_total (i, j, k, 7) +
c     &            b_3 (i, j, k) * s_total (i, j, k, 8) *
c     &           sinus_theta (j)) * r (i)) * phi (i, j, k) ** 2
            
c            s_total (i, j, k, 9) =
c     &           - 0.125_rk * delta_t_inv * psi (i, j, k) 
            s_total (i, j, k, 9) =
     &           - 0.125_rk * psi (i, j, k) *
     &           min (delta_t_inv,
     &           sqrt (c_sound_squared (i, j, k) +
     &           (b_1 (i, j, k) ** 2 + b_2 (i, j, k) ** 2 +
     &           b_3 (i, j, k) ** 2) / rho (i, j, k)) /
     &           min (delta_r (i), r (i) * delta_theta))
#ifdef CFC_TRANSPORT2
#error
#endif            
#endif /* CFC_MHD */
            
         end do

#ifdef ORIGIN
         s_total (1,j,k,2)=0.0_rk
#endif
            
      end do                    !j,k-loop


      return

      END SUBROUTINE calculate_sources

c     ==================================================================



#ifndef CFC_TRANSPORT2
c     ==================================================================

      SUBROUTINE calculate_newtonian_potential

c     ==================================================================

      USE precision

      USE size_cfc
      USE grid_cfc
      USE hydro_primitives_cfc
      USE perm_aux_cfc
      USE metric_cfc
      USE modpot_cfc
      USE parameters_cfc
      USE temp_aux_cfc
      USE poisson

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk

      if (newtonian_potential_form.ne.1)
     &     call calculate_modified_potential
      
c     use expansion method to calculate Newtonian potential
c     (gives smaller gravitational wave signal for rotating neutron star in equilibrium)

      if (newtonian_potential_form.eq.9) then

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            do i = 1, m
               rho_grav (i, j, k) = 4.0_rk * pi * 
     &              (rho (i, j, k) + p (i, j, k)) *
     +              (1.0_rk + eps (i, j, k)) / ((1.0_rk - 2.0_rk *
     +              m_spherical_tov (i) * r_inv (i)) *
     +              sqrt (1.0_rk - v_2 (i, j, k) ** 2 *
     +              rr_inv (i) - v_3 (i, j, k) ** 2 *
     +              rr_inv (i) * ssinus_theta_inv (j)))
            end do
         end do

      else 

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_loc * o_loc

            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            
            do i = 1, m
               rho_grav (i, j, k) = 4.0_rk * pi * rho (i, j, k)
            end do
         end do
      end if

      call loese_poisson(rho_grav, phi_potential)

c     adjust modified Newtonian potential in r-derivatives
      if (newtonian_potential_form.ne.1) then

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINLGE
#endif
         
c     adjust Newtonian potential in r-derivatives
         do i = 0, m + 1
            delta_potential_spherical (i) =
     &           potential_spherical_tov (i) -
     &           potential_spherical_newton (i)
         enddo
         
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINLGE
#endif
         
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            
            do i = 0, m + 1
               phi_potential (i, j, k) =
     &              phi_potential (i, j, k) +
     &                 delta_potential_spherical (i)
            enddo

            if (j .eq. n_s) then
               do i = 0, m + 1
                  phi_potential (i, j - 1, k) =
     &                 phi_potential (i, j - 1, k) +
     &                 delta_potential_spherical (i)
               enddo
            end if
            if (j .eq. n_e) then
               do i = 0, m + 1
                  phi_potential (i, j + 1, k) =
     &                 phi_potential (i, j + 1, k) +
     &                 delta_potential_spherical (i)
               enddo
            end if
            if (k .eq. o_s) then
               do i = 0, m + 1
                  phi_potential (i, j, k - 1) =
     &                 phi_potential (i, j, k - 1) +
     &                 delta_potential_spherical (i)
               enddo
            end if
            if (k .eq. o_e) then
               do i = 0, m + 1
                  phi_potential (i, j, k + 1) =
     &                 phi_potential (i, j, k + 1) +
     &                 delta_potential_spherical (i)
               enddo
            end if


         enddo
      endif



c     ------------------------------------------------------------------
c     Calculate derivatives of the the gravitational potential

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
      do jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         
         do i = 1, m
            
            dphi_potential_dr (i, j, k) =
     &           delta_r_amp_inv (i) *
     &           (phi_potential (i + 1, j, k) +
     &           (aamp (i) - 1.0_rk) * phi_potential (i, j, k) -
     &           aamp (i) * phi_potential (i - 1, j, k))
               
            dphi_potential_dtheta (i, j, k) =
     &           0.5_rk * delta_theta_inv *
     &           (phi_potential (i, j + 1, k) -
     &           phi_potential (i, j - 1, k))
            
            dphi_potential_dvarphi (i, j, k) = 0.0_rk
            
         enddo
      enddo
                     
            
      return

      END SUBROUTINE calculate_newtonian_potential

c     ==================================================================




c     ==================================================================

      SUBROUTINE add_newtonian_sources

c     ==================================================================

      USE precision

      USE size_cfc
      USE hydro_primitives_cfc
      USE perm_aux_cfc
      USE metric_cfc
      USE sources_cfc
      USE parameters_cfc

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
      do jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m
            s_total (i, j, k, 2) = s_total (i, j, k, 2) -
     &           rho (i, j, k) * dphi_potential_dr (i, j, k)
            s_total (i, j, k, 3) = s_total (i, j, k, 3) -
     &           rho (i, j, k) * dphi_potential_dtheta (i, j, k)
            s_total (i, j, k, 4) = s_total (i, j, k, 4) -
     &           rho (i, j, k) * dphi_potential_dvarphi (i, j, k)
#ifndef ENECONS
            s_total (i, j, k, 5) = s_total (i, j, k, 5) -
     &           rho (i, j, k) *
     &           (v_1 (i, j, k) * g_up_11 (i, j, k) *
     &           dphi_potential_dr (i, j, k) +
     &           v_2 (i, j, k) * g_up_22 (i, j, k) *
     &           dphi_potential_dtheta (i, j, k) +
     &           v_3 (i, j, k) * g_up_33 (i, j, k) *
     &           dphi_potential_dvarphi (i, j, k))
#endif
         enddo
      enddo

      return

      END SUBROUTINE add_newtonian_sources

c     ==================================================================



c     ==================================================================

      SUBROUTINE calculate_modified_potential

c     ==================================================================

c     Modified Newtonian potential according to Marek et al. (2006), 
c     Mueller et al. (2008).

c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE grid_cfc
      USE hydro_primitives_cfc
      USE modpot_cfc
      USE parameters_cfc
      USE perm_aux_cfc
      USE nutra_cfc

      USE totare_hy, ONLY: dvxtot,dvytot,dvztot
      USE totgrq_hy, ONLY: ephtot

#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk

#ifdef MPI_HYDRO
      real (kind=rk) :: buf (m + 1, 8), buf2 (m + 1, 8)
      integer (kind=ik) :: ierr
#endif

      real (kind=rk) :: dvom, grid_factor



#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif

c     calculate angular averaged values for density, pressure, enthalpy,
c     and radial velocity at cell centers and interfaces

      grid_factor = 0.0_rk
      do k = 1, o
         do j = 1, n
            grid_factor = grid_factor + dvytot (j) * dvztot (k)
         enddo
      end do
      grid_factor = 1.0_rk / grid_factor

      rho_tangential_1_spherical (:) = 0.0_rk
      rho_tangential_2_spherical (:) = 0.0_rk
      eps_spherical (:) = 0.0_rk
      p_spherical (:) = 0.0_rk
      v_1_spherical (:) = 0.0_rk
      pnu_spherical (:) = 0.0_rk


c#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
c!$OMP DO
c!$OMP& PRIVATE(i,j,k,jk)
c#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         
         dvom = dvytot (j) * dvztot (k)

         if (newtonian_potential_form.ne.9) then

            do i = 1, m
               rho_tangential_1_spherical (i) =
     &              rho_tangential_1_spherical (i) +
     &              rho (i, j, k) * dvom
               eps_spherical (i) = eps_spherical (i) +
     &              (eps (i, j, k) * rho (i, j, k) +
     &              enu (i, j, k)) * dvom
               p_spherical (i) = p_spherical (i) +
     &              p (i, j, k) * dvom
               v_1_spherical (i) = v_1_spherical (i) +
     &              v_1 (i, j, k) * rho (i, j, k) * dvom
               pnu_spherical (i) = pnu_spherical (i) +
     &              pnu (i, j, k) * dvom
            enddo
            rho_tangential_1_spherical (m + 1) =
     &           rho_tangential_1_spherical (m + 1) +
     &           rho (m + 1, j, k) * dvom
            eps_spherical (m + 1) = eps_spherical (m + 1) +
     &           (eps (m + 1, j, k) * rho (m + 1, j, k) +
     &           enu (m, j, k)) * dvom
            p_spherical (m + 1) = p_spherical (m + 1) +
     &           p (m + 1, j, k) * dvom
            v_1_spherical (m + 1) = v_1_spherical (m + 1) +
     &           v_1 (m + 1, j, k) * rho (m + 1, j, k) * dvom
            pnu_spherical (m + 1) = pnu_spherical (m + 1) +
     &           pnu (m, j, k) * dvom
            
         else

            do i = 1, m
               v_tangential_1 = sqrt (1.0_rk - v_2 (i, j, k) ** 2 *
     &              rr_inv (i) - v_3 (i, j, k) ** 2 * rr_inv (i) *
     &              ssinus_theta_inv (j))
               v_tangential_2 = v_tangential_1 ** 3
               rho_tangential_1_spherical (i) =
     &              rho_tangential_1_spherical (i) + rho (i, j, k) /
     &              v_tangential_1 * dvom
               rho_tangential_2_spherical (i) =
     &              rho_tangential_2_spherical (i) + rho (i, j, k) / 
     &              v_tangential_2 * dvom
               eps_spherical (i) = eps_spherical (i) +
     &              (eps (i, j, k) * rho (i, j, k) + 
     &              enu (i, j, k)) * dvom
               p_spherical (i) = p_spherical (i) +
     &              p (i, j, k) * dvom
               v_1_spherical (i) = v_1_spherical (i) +
     &              v_1 (i, j, k) * rho (i, j, k) * dvom
               pnu_spherical (i) = pnu_spherical (i) +
     &              pnu (i, j, k) * dvom
            enddo
            i = m + 1
            v_tangential_1 = sqrt (1.0_rk - v_2 (i, j, k) ** 2 *
     &           rr_inv (i) - v_3 (i, j, k) ** 2 * rr_inv (i) *
     &           ssinus_theta_inv (j))
            v_tangential_2 = v_tangential_1 ** 3
            rho_tangential_1_spherical (i) =
     &           rho_tangential_1_spherical (i) + rho (i, j, k) /
     &           v_tangential_1 * dvom
            rho_tangential_2_spherical (i) =
     &           rho_tangential_2_spherical (i) + rho (i, j, k) / 
     &           v_tangential_2 * dvom
            eps_spherical (i) = eps_spherical (i) +
     &           (eps (i, j, k) * rho (i, j, k) + 
     &           enu (m, j, k)) * dvom
            p_spherical (i) = p_spherical (i) +
     &           p (i, j, k) * dvom
            v_1_spherical (i) = v_1_spherical (i) +
     &           v_1 (i, j, k) * rho (i, j, k) * dvom
            pnu_spherical (i) = pnu_spherical (i) +
     &           pnu (m, j, k) * dvom

         endif
      enddo
         
      rho_mod_spherical (:) = 0.0_rk

c#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
c!$OMP DO
c!$OMP& PRIVATE(i,j,k,jk)
c#endif
      do jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         
         dvom = dvytot (j) * dvztot (k)
         
         if (newtonian_potential_form.ne.9) then
            do i = 1, m + 1
               rho_mod_spherical (i) =
     &              rho_mod_spherical (i) + rho (i, j, k) * dvom
            enddo
         else
            do i = 1, m
               rho_mod_spherical (i) =
     &              rho_mod_spherical (i) + (rho (i, j, k) +
     &              p (i, j, k)) * (1.0_rk + eps (i, j, k)) /
     &              (1.0_rk  - 2.0_rk * m_spherical_tov (i) *
     &              r_inv (i)) / sqrt (1.0_rk - v_2 (i, j, k) ** 2 *
     &              rr_inv (i) - v_3 (i, j, k) ** 2 * rr_inv (i) *
     &              ssinus_theta (j)) * dvom
            enddo
         endif
      enddo
      
#ifdef MPI_HYDRO
      buf (1 : m + 1, 1) = rho_tangential_1_spherical (1 : m + 1)
      buf (1 : m + 1, 2) = rho_tangential_2_spherical (1 : m + 1)
      buf (1 : m + 1, 3) = eps_spherical (1 : m + 1)
      buf (1 : m + 1, 4) = p_spherical (1 : m + 1)
      buf (1 : m + 1, 5) = h_spherical (1 : m + 1)
      buf (1 : m + 1, 6) = v_1_spherical (1 : m + 1)
      buf (1 : m + 1, 7) = pnu_spherical (1 : m + 1)
      buf (1 : m + 1, 8) = rho_mod_spherical (1 : m + 1)
      call MPI_Allreduce(buf, buf2, (m + 1_ik) * 8_ik,
     &     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      rho_tangential_1_spherical (1 : m + 1) = buf2 (1 : m + 1, 1)
      rho_tangential_2_spherical (1 : m + 1) = buf2 (1 : m + 1, 2)
      eps_spherical (1 : m + 1)              = buf2 (1 : m + 1, 3)
      p_spherical (1 : m + 1)                = buf2 (1 : m + 1, 4)
      h_spherical (1 : m + 1)                = buf2 (1 : m + 1, 5)
      v_1_spherical (1 : m + 1)              = buf2 (1 : m + 1, 6)
      pnu_spherical (1 : m + 1)              = buf2 (1 : m + 1, 7)
      rho_mod_spherical (1 : m + 1)          = buf2 (1 : m + 1, 8)
#endif /* MPI_HYDRO */         

      rho_tangential_1_spherical (:) =
     &     rho_tangential_1_spherical (:) * grid_factor
      rho_tangential_2_spherical (:) =
     &     rho_tangential_2_spherical (:) * grid_factor
      eps_spherical (:) = eps_spherical (:) * grid_factor /
     &     rho_tangential_1_spherical (:)
      p_spherical (:) = p_spherical (:) * grid_factor
      if (newtonian_potential_form.eq.6.or.
     &     newtonian_potential_form.eq.7) then
         h_spherical (:) = 1.0_rk + p_spherical (:) /
     &        rho_tangential_1_spherical (:)
      else
         h_spherical (:) = 1.0_rk + eps_spherical (:) +
     &        p_spherical (:) / rho_tangential_1_spherical (:)
      endif
      v_1_spherical (:) = v_1_spherical (:) * grid_factor /
     &     rho_tangential_1_spherical (:)
      pnu_spherical (:) = pnu_spherical (:) * grid_factor
      rho_mod_spherical (:) = rho_mod_spherical (:) * grid_factor
      
      if (newtonian_potential_form.eq.3.or.
     &     newtonian_potential_form.eq.7) then
         energy_factor = 0.0_rk
      else
         energy_factor = 1.0_rk
      endif

      if (newtonian_potential_form.eq.4) then
         velocity_factor = 0.0_rk
      else
         velocity_factor = 1.0_rk
      endif



c     calculate spherical TOV mass
      
      m_spherical_tov_if (0) = 0.0_rk
      
      do i = 1, m
         
         if (newtonian_potential_form.eq.8.or.
     &        newtonian_potential_form.eq.9) then
            gamma_mass = sqrt (1.0_rk - 2.0_rk *
     &           m_spherical_tov_if (i - 1) * r_inv (i))
         else
            gamma_mass = 1.0_rk
         endif
         
         m_spherical_tov_if (i) =
     &        m_spherical_tov_if (i - 1) + 4.0_rk * pi *
     &        rho_tangential_1_spherical (i) * gamma_mass *
     &        (1.0_rk + energy_factor * eps_spherical (i)) * dvxtot (i)
         
         m_spherical_tov (i) = 0.5_rk *
     &        (m_spherical_tov_if (i - 1) +
     &        m_spherical_tov_if (i))
         
      enddo
      
      
c     calculate spherical Newtonian mass
               
      m_spherical_newton_if (0) = 0.0_rk

      do i = 1, m
         m_spherical_newton_if (i) =
     &        m_spherical_newton_if (i - 1) +
     &        4.0_rk * pi * rho_mod_spherical (i) * dvxtot (i)
      enddo
      
      
c     calculate spherical TOV potential and Newtonian potential
      
      potential_spherical_tov (m + 1) = 
     &     - m_spherical_tov_if (m) * r_inv (m + 1)
      potential_spherical_newton (m + 1) =
     &     - m_spherical_newton_if (m) * r_inv (m + 1)
      
      do i = m, 1, -1

         if (newtonian_potential_form.eq.5) then
            gamma = 1.0_rk
         else
            gamma = (1.0_rk + v_1_spherical (i) ** 2 * velocity_factor -
     &           2.0_rk * m_spherical_tov_if (i) * r_inv (i))
         endif
         
         if (newtonian_potential_form.ne.9) then
            gamma_tangential = 1.0_rk
         else
            gamma_tangential =
     &           rho_tangential_2_spherical (i) /
     &           rho_tangential_1_spherical (i)
         endif
         
         potential_spherical_tov (i) =
     &        potential_spherical_tov (i + 1) -
     &        delta_r_if (i) * m_spherical_tov_if (i) * 
     &        0.5_rk * (h_spherical (i) + h_spherical (i + 1)) /
     &        rr_if (i) / gamma * (1.0_rk + 4.0_rk * pi * rr_if (i) *
     &        r_if (i) * 
     &        0.5_rk * (p_spherical (i) + p_spherical (i + 1)) / 
     &        m_spherical_tov_if (i)) *
     &        gamma_tangential
         
         potential_spherical_newton (i) =
     &        potential_spherical_newton (i + 1) -
     &        delta_r_if (i) * m_spherical_newton_if (i) / rr_if (i)

      enddo
      potential_spherical_tov (0) = potential_spherical_tov (1) 
      potential_spherical_newton (0) = potential_spherical_newton (1)


      ephtot (m) = potential_spherical_tov (m)
      do i = m - 1, 1, -1
         ephtot (i) = ephtot (i + 1) +
     &        (potential_spherical_tov (i) -
     &        potential_spherical_tov (i + 1)) * 2.0_rk /
     &        (h_spherical (i) + h_spherical (i + 1))
      end do
      ephtot (0) = ephtot (1)
      ephtot (:) = exp (ephtot (:))

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif
      
      return

      END SUBROUTINE calculate_modified_potential

c     ==================================================================

#endif /* CFC_TRANSPORT2 */



      END MODULE gr_sources

c     ==================================================================
