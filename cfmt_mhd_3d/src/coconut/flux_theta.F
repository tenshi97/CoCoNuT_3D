#if defined (CFC_MHD) && !defined(PASSIVE_MHD)
#define MHDFLX
#endif
      
c     ==================================================================

      MODULE gr_flux_theta

c     ==================================================================

      IMPLICIT NONE

      CONTAINS


#ifdef CFC_TRANSPORT2


c     ==================================================================

      SUBROUTINE calculate_flux_theta

c     ==================================================================

c     Computes HLLC/HLLE fluxes in theta-direction (GR version)

c     ------------------------------------------------------------------

      USE precision
      
      USE size_cfc
      USE parameters_cfc
      USE conserved_cfc
      USE interface_hydro_cfc
      USE grid_cfc
      USE fluxes_cfc
      USE metric_cfc
      USE perm_aux_cfc
      USE hydro_primitives_cfc

      USE shock_cfc
      
      use configure
      
      IMPLICIT NONE

      integer (kind=ik) :: i, j, jk, k, l, ll

      real (kind=rk) :: f_num_theta_if (1 : 5, 0 : 1)

      real (kind=rk) :: w_if (0 : 1)

      real (kind=rk) :: d_cap_if (0 : 1), s_1_if (0 : 1),
     &     s_2_if (0 : 1), s_3_if (0 : 1), tau_if (0 : 1)

      real (kind=rk) :: v_up_2_if (0 : 1)
      real (kind=rk) :: v_up_2_hat_if (0 : 1)

      real (kind=rk) :: lambda_minus_if (0 : 1)
      real (kind=rk) :: lambda_plus_if (0 : 1)

      real (kind=rk) :: lambda_minus_min_if, lambda_plus_max_if,
     &     lambda_diff_if
      real (kind=rk) :: c_sound_if (0 : 1)

      real (kind=rk) :: tmp_1, tmp_2, tmp_3, rho_h_ww_if

      real (kind=rk) :: e_hll,s_2_hll,lambda_star,lambda_star_up,
     &     p_star,d_cap_star,s_1_star,s_2_star,s_3_star,
     &     tau_star,scr1
#ifdef CFC_MHD
      real (kind=rk) :: v_1_upw, v_2_upw, v_3_upw
#endif /* CFC_MHD */

      logical :: lscr1

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,l,ll)
#endif
      do jk = 1, (n_loc + 1) * o_loc

         k = int((jk + (n_loc + 1) - 1) / (n_loc + 1))
         j = (n_s - 2) + (jk - (k - 1) * (n_loc + 1))
         k = k + o_s - 1

         do i = 1, m

ccc   ==================================================================
ccc   left interface
            
            c_sound_if (0) =
     &           sqrt (c_sound_squared_if (i, j + 1, k, 0))
ccc   calculate sound speed c_s
            
            w_if (0) = 1.0_rk /
     &           sqrt (1.0_rk - v_squared_if (i, j + 1, k, 0))
ccc   calculate Lorentz factor W
            
            v_up_2_if (0) =
     &           g_up_22_theta (i, j, k) * v_2_if (i, j + 1, k, 0)
ccc   calculate interface value of v^i = g^k v_i

            v_up_2_hat_if (0) = v_up_2_if (0) -
     &           beta_up_2_theta (i, j, k) / alpha_theta (i, j, k)
            rho_h_ww_if = rho_if (i, j + 1, k, 0) *
     &           h_if (i, j + 1, k, 0) * w_if (0) ** 2
            
            d_cap_if (0) = rho_if (i, j + 1, k, 0) * w_if (0)
            s_1_if (0) = rho_h_ww_if * v_1_if (i, j + 1, k, 0)
            s_2_if (0) = rho_h_ww_if * v_2_if (i, j + 1, k, 0)
            s_3_if (0) = rho_h_ww_if * v_3_if (i, j + 1, k, 0)
            tau_if (0) = rho_h_ww_if - d_cap_if (0) -
     &           p_if (i, j + 1, k, 0)
            
ccc   ==================================================================
               
ccc   calculate auxiliary terms: tmp_1, tmp_2, tmp_3
            tmp_1 = v_up_2_if (0) *
     &           (1.0_rk - c_sound_squared_if (i, j + 1, k, 0))
ccc   tmp_1 = v^1 (1 - v_s^2)
            tmp_2 = (1.0_rk - v_squared_if (i, j + 1, k, 0)) *
     &           (g_up_22_theta (i, j, k) *
     &           (1.0_rk - v_squared_if (i, j + 1, k, 0) *
     &           c_sound_squared_if (i, j + 1, k, 0)) - 
     &           v_up_2_if (0) ** 2 *
     &           (1.0_rk - c_sound_squared_if (i, j + 1, k, 0)))
ccc   tmp_2 = (1 - v^2) (g^11 (1 - c_s^2 v^2) - v^1 v^1 (1 - c_s^2))
            tmp_3 = alpha_theta (i, j, k) / (1.0_rk - 
     &           v_squared_if (i, j + 1, k, 0) *
     &           c_sound_squared_if (i, j + 1, k, 0))
ccc   tmp_3 = alpha / (1 - c_s^2 v^2)
               
ccc   calculate eigenvalues --------------------------------------------
            lambda_minus_if (0) = tmp_3 * (tmp_1 - c_sound_if (0) * 
     &           sqrt (tmp_2)) - beta_up_2_theta (i, j, k)
ccc   lambda_-
            lambda_plus_if (0) = tmp_3 * (tmp_1 + c_sound_if (0) * 
     &           sqrt (tmp_2)) - beta_up_2_theta (i, j, k)
ccc   lambda_+



ccc   ==================================================================
ccc   right interface

            c_sound_if (1) = sqrt (c_sound_squared_if (i, j, k, 1))
ccc   calculate sound speed c_s

            w_if (1) = 1.0_rk /
     &           sqrt (1.0_rk - v_squared_if (i, j, k, 1))
ccc   calculate Lorentz factor W
               
            v_up_2_if (1) =
     &           g_up_22_theta (i, j, k) * v_2_if (i, j, k, 1)
ccc   calculate interface value of v^i = g^k v_i
               
            v_up_2_hat_if (1) = v_up_2_if (1) -
     &           beta_up_2_theta (i, j, k) / alpha_theta (i, j, k)
            rho_h_ww_if = rho_if (i, j, k, 1) *
     &           h_if (i, j, k, 1) * w_if (1) ** 2
            
            d_cap_if (1) = rho_if (i, j, k, 1) * w_if (1)
            s_1_if (1) = rho_h_ww_if * v_1_if (i, j, k, 1)
            s_2_if (1) = rho_h_ww_if * v_2_if (i, j, k, 1)
            s_3_if (1) = rho_h_ww_if * v_3_if (i, j, k, 1)
            tau_if (1) =
     &           rho_h_ww_if - d_cap_if (1) - p_if (i, j, k, 1)
               
ccc   ==================================================================
               
ccc   calculate auxiliary terms: tmp_1, tmp_2, tmp_3
            tmp_1 = v_up_2_if (1) *
     &           (1.0_rk - c_sound_squared_if (i, j, k, 1))
ccc   tmp_1 = v^1 (1 - v_s^2)
            tmp_2 = (1.0_rk - v_squared_if (i, j, k, 1)) *
     &           (g_up_22_theta (i, j, k) *
     &           (1.0_rk - v_squared_if (i, j, k, 1) *
     &           c_sound_squared_if (i, j, k, 1)) - 
     &           v_up_2_if (1) ** 2 *
     &           (1.0_rk - c_sound_squared_if (i, j, k, 1)))
ccc   tmp_2 = (1 - v^2) (g^11 (1 - c_s^2 v^2) - v^1 v^1 (1 - c_s^2))
            tmp_3 = alpha_theta (i, j, k) / (1.0_rk - 
     &           v_squared_if (i, j, k, 1) *
     &           c_sound_squared_if (i, j, k, 1))
ccc   tmp_3 = alpha / (1 - c_s^2 v^2)
            
ccc   calculate eigenvalues --------------------------------------------
            lambda_minus_if (1) = tmp_3 * (tmp_1 - c_sound_if (1) *
     &           sqrt (tmp_2)) - beta_up_2_theta (i, j, k)
ccc   lambda_-
            lambda_plus_if (1) = tmp_3 * (tmp_1 + c_sound_if (1) * 
     &           sqrt (tmp_2)) - beta_up_2_theta (i, j, k)
ccc   lambda_+



ccc   ================================================================= 
ccc   calculate maximum of left and right interface eigenvalues
ccc   and their difference
            lambda_minus_min_if = min (lambda_minus_if (0),
     &           lambda_minus_if (1), 0.0_rk)
            lambda_plus_max_if = max (lambda_plus_if (0),
     &           lambda_plus_if (1), 0.0_rk)
            
            lambda_diff_if = lambda_plus_max_if - lambda_minus_min_if



ccc   ================================================================= 
ccc   calculate flux at left and right interface            
            do l = 0, 1  
               f_num_theta_if (1, l) =
     &              d_cap_if (l) * v_up_2_hat_if (l)
               f_num_theta_if (2, l) =
     &              s_1_if (l) * v_up_2_hat_if (l)
               f_num_theta_if (3, l) =
     &              s_2_if (l) * v_up_2_hat_if (l) + 
     &              p_if (i, j + 1 - l, k, l)
               f_num_theta_if (4, l) =
     &              s_3_if (l) * v_up_2_hat_if (l)
               f_num_theta_if (5, l) =
     &              tau_if (l) * v_up_2_hat_if (l) + 
     &              p_if (i, j + 1 - l, k, l) * v_up_2_if (l)
               
            enddo



ccc   ==================================================================
ccc   calculate numerical fluxes at interface
            f_num_if (i, j, k, 1) =
     &           (lambda_plus_max_if * f_num_theta_if (1, 1) -
     &           lambda_minus_min_if * f_num_theta_if (1, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (d_cap_if (0) - d_cap_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 2) =
     &           (lambda_plus_max_if * f_num_theta_if (2, 1) -
     &           lambda_minus_min_if * f_num_theta_if (2, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_1_if (0) - s_1_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 3) =
     &           (lambda_plus_max_if * f_num_theta_if (3, 1) -
     &           lambda_minus_min_if * f_num_theta_if (3, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_2_if (0) - s_2_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 4) =
     &           (lambda_plus_max_if * f_num_theta_if (4, 1) -
     &           lambda_minus_min_if * f_num_theta_if (4, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_3_if (0) - s_3_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 5) =
     &           (lambda_plus_max_if * f_num_theta_if (5, 1) -
     &           lambda_minus_min_if * f_num_theta_if (5, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (tau_if (0) - tau_if (1))) / lambda_diff_if
               
c     calculate HLL state vector (only 2 components needed)
            e_hll = (lambda_plus_max_if *
     &           (tau_if (0) + d_cap_if (0))-
     &           lambda_minus_min_if *
     &           (tau_if (1) + d_cap_if (1))+
     &           (f_num_theta_if (5, 1) +
     &           f_num_theta_if (1, 1)) -
     &           (f_num_theta_if (5, 0) +
     &           f_num_theta_if (1, 0))) /
     &           lambda_diff_if
            s_2_hll = (lambda_plus_max_if * s_2_if (0) -
     &           lambda_minus_min_if * s_2_if (1) +
     &           f_num_theta_if (3, 1) -
     &           f_num_theta_if (3, 0)) /
     &           lambda_diff_if
c     calculate physical root of (18)
            scr1=e_hll+f_num_if(i,j,k,3)
            lambda_star=2.0_rk*s_2_hll/
     &           (scr1+sqrt(scr1**2-
     &           4.0_rk*
     &           (f_num_if(i,j,k,5)+
     &           f_num_if(i,j,k,1))*
     &           s_2_hll))
            lambda_star_up=lambda_star*g_up_22_theta(i,j,k)-
     &           beta_up_2_theta(i,j,k)/alpha_theta(i,j,k)

c     calculate state variables for intermediate state
c     (depending on the sign of lambda_star_up):
            if (lambda_star_up .lt. 0.0_rk) then
               p_star=lambda_star*
     &              (lambda_plus_max_if *
     &              (tau_if(0)+d_cap_if(0)+
     &              p_if(i,j+1,k,0))-
     &              s_2_if(0)*g_up_22_theta(i,j,k))-
     &              s_2_if(0)*(lambda_plus_max_if-
     &              v_up_2_hat_if(0))+p_if(i,j+1,k,0)
               d_cap_star=d_cap_if(0)*
     &              (lambda_plus_max_if-v_up_2_hat_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_1_star=s_1_if(0)*
     &              (lambda_plus_max_if-v_up_2_hat_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_2_star=(s_2_if(0)*
     &              (lambda_plus_max_if-v_up_2_hat_if(0))+
     &              p_star-p_if(i,j+1,k,0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_3_star=s_3_if(0)*
     &              (lambda_plus_max_if-v_up_2_hat_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               tau_star=(tau_if(0)*
     &              (lambda_plus_max_if-v_up_2_hat_if(0))+
     &              p_star*lambda_star_up-
     &              p_if(i,j+1,k,0)*v_up_2_hat_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
            else
               p_star=lambda_star*
     &              (lambda_minus_min_if *
     &              (tau_if(1)+d_cap_if(1)+
     &              p_if(i,j,k,1))-
     &              s_2_if(1)*g_up_22_theta(i,j,k))-
     &              s_2_if(1)*(lambda_minus_min_if-
     &              v_up_2_hat_if(1))+p_if(i,j,k,1)
               d_cap_star=d_cap_if(1)*
     &              (lambda_minus_min_if-v_up_2_hat_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_1_star=s_1_if(1)*
     &              (lambda_minus_min_if-v_up_2_hat_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_2_star=(s_2_if(1)*
     &              (lambda_minus_min_if-v_up_2_hat_if(1))+
     &              p_star-p_if(i,j,k,1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_3_star=s_3_if(1)*
     &              (lambda_minus_min_if-v_up_2_hat_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               tau_star=(tau_if(1)*
     &              (lambda_minus_min_if-v_up_2_hat_if(1))+
     &              p_star*lambda_star_up-
     &              p_if(i,j,k,1)*v_up_2_hat_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
            end if
            if (lambda_plus_max_if .le. 0.0_rk) then
               p_star=p_if(i,j+1,k,0)
            end if
            
            if (lambda_minus_min_if .ge. 0.0_rk) then
               p_star=p_if(i,j  ,k,1)
            end if
            
            pav_theta (i, j, k) = p_star

c     Use HLLC flux for subsonic flow
            lscr1=(lambda_plus_max_if*
     &           lambda_minus_min_if .lt. 0.0_rk) .and.
     &           .not.(lshock(i,j,k) .or. lshock(i,j+1,k))
            if (lscr1) then
               f_num_if (i, j, k, 1) =
     &              d_cap_star*lambda_star_up
               f_num_if (i, j, k, 2) =
     &              s_1_star*lambda_star_up
#if MOMENTUM_EQUATION_MODE==0
               f_num_if (i, j, k, 3) =
     &              s_2_star*lambda_star_up+p_star
#else
               f_num_if (i, j, k, 3) =
     &              s_2_star*lambda_star_up
#endif /* MOMENTUM_EQUATION_MODE==0 */
               f_num_if (i, j, k, 4) =
     &              s_3_star*lambda_star_up
               f_num_if (i, j, k, 5) =
     &              tau_star*lambda_star_up+
     &              p_star*(lambda_star_up+
     &              beta_up_2_theta(i,j,k)/
     &              alpha_theta(i,j,k))
#if MOMENTUM_EQUATION_MODE>0
            else
c     Recalculate HLLE momentum flux (without pressure contribution)
               f_num_theta_if (3, 0) =
     &              s_2_if (0) * v_up_2_hat_if (0)
               f_num_theta_if (3, 1) =
     &              s_2_if (1) * v_up_2_hat_if (1) 
               f_num_if (i, j, k, 3) =
     &              (lambda_plus_max_if * f_num_theta_if (3, 1) -
     &              lambda_minus_min_if * f_num_theta_if (3, 0) +
     &              lambda_plus_max_if * lambda_minus_min_if *
     &              (s_2_if (0) - s_2_if (1))) / lambda_diff_if
c     Recalculate effective pressure at interface
               pav_theta (i, j, k) =
     &              (lambda_plus_max_if * p_if (i, j    , k, 1) -
     &              lambda_minus_min_if * p_if (i, j + 1, k, 0)) / 
     &              lambda_diff_if
#endif /* MOMENTUM_EQUATION_MODE>0 */
            end if


            if (f_num_if (i, j, k, 1) .gt. 0.0_rk) then
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
               f_num_xnu_if (i, j, k, 1 : config%qn) = 
     &              f_num_if (i, j, k, 1) *
     &              xnu_if (i, j, k, 1, 1 : config%qn)
            else
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
               f_num_xnu_if (i, j, k, 1 : config%qn) =
     &              f_num_if (i, j, k, 1) *
     &              xnu_if (i, j + 1, k, 0, 1 : config%qn)
            end if

#undef EDDY
#ifdef EDDY
            scr1=ABS (v_1 (i, j, k) - v_1 (i, j + 1, k))
            f_num_if (i, j, k, 1) = f_num_if (i, j, k, 1) + scr1 *
     &           (d_cap_hat (i, j, k) - d_cap_hat (i, j + 1, k))
            f_num_if (i, j, k, 2) = f_num_if (i, j, k, 2) + scr1 *
     &           (s_1_hat (i, j, k) - s_1_hat (i, j + 1, k))
            f_num_if (i, j, k, 3) = f_num_if (i, j, k, 3) + scr1 *
     &           (s_2_hat (i, j, k) - s_2_hat (i, j + 1, k))
            f_num_if (i, j, k, 4) = f_num_if (i, j, k, 4) + scr1 *
     &           (s_3_hat (i, j, k) - s_3_hat (i, j + 1, k))
            f_num_if (i, j, k, 5) = f_num_if (i, j, k, 5) + scr1 *
     &           (tau_hat (i, j, k) - tau_hat (i, j + 1, k))
            f_num_xnu_if (i, j, k, 1 : config%qn) =
     &           f_num_xnu_if (i, j, k, 1 : config%qn) +
     &           (d_cap_xnu_hat (i, j, k, 1 : config%qn) -
     &           d_cap_xnu_hat (i, j + 1, k, 1 : config%qn))
#endif

#ifdef CFC_MHD
            if (lscr1) then 
               v_1_upw = s_1_star / (tau_star + d_cap_star + p_star)
               v_2_upw = lambda_star_up
               v_3_upw = s_3_star / (tau_star + d_cap_star + p_star)
            else if (lambda_plus_max_if .le. 0.0_rk) then
               v_1_upw = v_1_if (i, j + 1, k, 0)
               v_2_upw = v_up_2_hat_if (0)
               v_3_upw = v_3_if (i, j + 1, k, 0)
            else
               v_1_upw = v_1_if (i, j, k, 1)
               v_2_upw = v_up_2_hat_if (1)
               v_3_upw = v_3_if (i, j, k, 1)
            end if
c     Compute v* and transform to orthonormal components
            v_1_upw = (alpha_theta (i,j,k) * v_1_upw *
     &           g_up_11_theta (i,j,k) - beta_up_1_theta (i,j,k))
            v_2_upw = (alpha_theta (i,j,k) * v_2_upw) *
     &           r (i)
            v_3_upw = (alpha_theta (i,j,k) * v_3_upw *
     &           g_up_33_theta (i,j,k) -
     &           beta_up_3_theta (i,j,k)) *
     &           r (i) * sinus_theta_if (j)
c     Test
c$$$            if (r (i) .lt. 1e8_rk) then
c$$$               v_1_upw = -5e-10_rk * r (i)
c$$$            else
c$$$               v_1_upw = -5e-10_rk * 1e8_rk * sqrt(1e8_rk/r(i))
c$$$            end if
c$$$            v_2_upw = 0.0d0
c$$$            v_3_upw = 0.0d0
#endif /* CFC_MHD */
#ifdef PASSIVE_MHD
            if (v_2_upw .gt. 0.0_rk) then
               emf_face (i, j, k, 1, 2) =
     &              v_2_upw * b_3_if (i, j, k, 1) -
     &              v_3_upw * b_2_if (i, j, k, 1)
               emf_face (i, j, k, 2, 2) =
     &              v_1_upw * b_2_if (i, j, k, 1) -
     &              v_2_upw * b_1_if (i, j, k, 1)
            else
               emf_face (i, j, k, 1, 2) =
     &              v_2_upw * b_3_if (i, j+1, k, 0) -
     &              v_3_upw * b_2_if (i, j+1, k, 0)
               emf_face (i, j, k, 2, 2) =
     &              v_1_upw * b_2_if (i, j+1, k, 0) -
     &              v_2_upw * b_1_if (i, j+1, k, 0)
c     &              v_1_upw * sinus_theta_if (j)
c     &              0.0_rk
            end if
c            if (i.eq.32) print *,'flux_emf',
c     &           v_1_upw,b_2_if (i, j+1, k, 0),
c     &           v_2_upw,b_1_if (i, j+1, k, 0)
            emf_face (i, j, k, 1, 2) = emf_face (i, j, k , 1, 2) *
     &           phi_theta (i, j, k) ** 6
            emf_face (i, j, k, 2, 2) = emf_face (i, j, k , 2, 2) *
     &           phi_theta (i, j, k) ** 6
#endif /* PASSIVE_MHD */

         enddo
      enddo

c      print *,'emf_chk1',emf_face(32,:,1,2,2)

      return

      END SUBROUTINE calculate_flux_theta

c     ==================================================================


#else /* CFC_TRANSPORT2 */

#ifndef MHDFLX
      
c     ==================================================================

      SUBROUTINE calculate_flux_theta

c     ==================================================================

c     Computes HLLC/HLLE fluxes in theta-direction (Newtonian version)

c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE parameters_cfc
      USE interface_hydro_cfc
      USE grid_cfc
      USE fluxes_cfc
      USE metric_cfc
      USE perm_aux_cfc
      USE shock_cfc

      use configure
      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk, l, ll

      real (kind=rk) :: f_num_theta_if (1 : 5, 0 : 1)

      real (kind=rk) :: d_cap_if (0 : 1), s_1_if (0 : 1),
     &     s_2_if (0 : 1),  s_3_if (0 : 1), tau_if (0 : 1)
      real (kind=rk) :: v_up_2_if (0 : 1)

      real (kind=rk) :: lambda_minus_if (0 : 1)
      real (kind=rk) :: lambda_plus_if (0 : 1)

      real (kind=rk) :: lambda_minus_min_if, lambda_plus_max_if,
     &     lambda_diff_if
      real (kind=rk) :: c_sound_if (0 : 1)

      real (kind=rk) :: e_hll,s_2_hll,lambda_star,lambda_star_up,
     &     p_star,d_cap_star,s_1_star,s_2_star,s_3_star,
     &     tau_star,scr1
#ifdef CFC_MHD
      real (kind=rk) :: v_1_upw, v_2_upw, v_3_upw
#endif /* CFC_MHD */
      logical :: lscr1


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,l,ll)
#endif
      do jk = 1, (n_loc + 1) * o_loc

         k = int((jk + (n_loc + 1) - 1) / (n_loc + 1))
         j = (n_s - 2) + (jk - (k - 1) * (n_loc + 1))
         k = k + o_s - 1

         do i = 1, m

ccc   ==================================================================
ccc   left interface
               
            v_up_2_if (0) =
     &           g_up_22_theta (i, j, k) * v_2_if (i, j + 1, k, 0)
ccc   calculate interface value of v^i = g^k v_i
               
            d_cap_if (0) = rho_if (i, j + 1, k, 0)
            s_1_if (0) = rho_if (i, j + 1, k, 0) *
     &           v_1_if (i, j + 1, k, 0)
            s_2_if (0) = rho_if (i, j + 1, k, 0) *
     &           v_2_if (i, j + 1, k, 0)
            s_3_if (0) = rho_if (i, j + 1, k, 0) *
     &           v_3_if (i, j + 1, k, 0)
            tau_if (0) = rho_if (i, j + 1, k, 0) *
     &           h_if (i, j + 1, k, 0) - p_if (i, j + 1, k, 0)
            
            c_sound_if (0) =
     &           sqrt (c_sound_squared_if (i, j + 1, k, 0))
ccc   calculate sound speed c_s
            
ccc   calculate eigenvalues --------------------------------------------
            lambda_minus_if (0) = v_up_2_if (0) - c_sound_if (0) *
     &           sqrt (g_up_22_theta (i, j, k))
ccc   lambda_-
            lambda_plus_if (0) = v_up_2_if (0) + c_sound_if (0) *
     &           sqrt (g_up_22_theta (i, j, k))
ccc   lambda_+

ccc   ==================================================================
ccc   right interface
            
            v_up_2_if (1) =
     &           g_up_22_theta (i, j, k) * v_2_if (i, j, k, 1)
ccc   calculate interface value of v^i = g^k v_i
               
            d_cap_if (1) = rho_if (i, j, k, 1)
            s_1_if (1) = rho_if (i, j, k, 1) * v_1_if (i, j, k, 1)
            s_2_if (1) = rho_if (i, j, k, 1) * v_2_if (i, j, k, 1)
            s_3_if (1) = rho_if (i, j, k, 1) * v_3_if (i, j, k, 1)
            tau_if (1) = rho_if (i, j, k, 1) *
     &           h_if (i, j, k, 1) - p_if (i, j, k, 1)
            
            c_sound_if (1) = sqrt (c_sound_squared_if (i, j, k, 1))
ccc   calculate sound speed c_s
            
ccc   calculate eigenvalues --------------------------------------------
            lambda_minus_if (1) = v_up_2_if (1) - c_sound_if (1) *
     &           sqrt (g_up_22_theta (i, j, k))
ccc   lambda_-
            lambda_plus_if (1) = v_up_2_if (1) + c_sound_if (1) *
     &           sqrt (g_up_22_theta (i, j, k))
ccc   lambda_+

ccc   ==================================================================
ccc   calculate maximum of left and right interface eigenvalues
ccc   and their difference
            lambda_minus_min_if = min (lambda_minus_if (0),
     &           lambda_minus_if (1), 0.0_rk)
            lambda_plus_max_if = max (lambda_plus_if (0),
     &           lambda_plus_if (1), 0.0_rk)
            
            lambda_diff_if = lambda_plus_max_if - lambda_minus_min_if

ccc   ==================================================================
ccc   calculate flux at left and right interface            
            do l = 0, 1
               f_num_theta_if (1, l) = d_cap_if (l) * v_up_2_if (l)
               f_num_theta_if (2, l) = s_1_if (l) * v_up_2_if (l)
               f_num_theta_if (3, l) = s_2_if (l) * v_up_2_if (l) + 
     &              p_if (i, j + 1 - l, k, l)
               f_num_theta_if (4, l) = s_3_if (l) * v_up_2_if (l)
               f_num_theta_if (5, l) = (tau_if (l) +
     &              p_if (i, j + 1 - l, k, l)) * v_up_2_if (l)
            enddo


ccc   ==================================================================
ccc   calculate numerical fluxes at interface
            f_num_if (i, j, k, 1) =
     &           (lambda_plus_max_if * f_num_theta_if (1, 1) -
     &           lambda_minus_min_if * f_num_theta_if (1, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (d_cap_if (0) - d_cap_if (1))) / lambda_diff_if
               
            f_num_if (i, j, k, 2) =
     &           (lambda_plus_max_if * f_num_theta_if (2, 1) -
     &           lambda_minus_min_if * f_num_theta_if (2, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_1_if (0) - s_1_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 3) =
     &           (lambda_plus_max_if * f_num_theta_if (3, 1) -
     &           lambda_minus_min_if * f_num_theta_if (3, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_2_if (0) - s_2_if (1))) / lambda_diff_if
               
            f_num_if (i, j, k, 4) =
     &           (lambda_plus_max_if * f_num_theta_if (4, 1) -
     &           lambda_minus_min_if * f_num_theta_if (4, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_3_if (0) - s_3_if (1))) / lambda_diff_if
               
            f_num_if (i, j, k, 5) =
     &           (lambda_plus_max_if * f_num_theta_if (5, 1) -
     &           lambda_minus_min_if * f_num_theta_if (5, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (tau_if (0) - tau_if (1))) / lambda_diff_if
               
c     calculate HLL state vector (only 2 components needed)
            e_hll = (lambda_plus_max_if * d_cap_if (0)-
     &           lambda_minus_min_if * d_cap_if (1)+
     &           f_num_theta_if (1, 1) -
     &           f_num_theta_if (1, 0)) /
     &           lambda_diff_if
            s_2_hll = (lambda_plus_max_if * s_2_if (0) -
     &           lambda_minus_min_if * s_2_if (1) +
     &           f_num_theta_if (3, 1) -
     &           f_num_theta_if (3, 0)) /
     &           lambda_diff_if
c     calculate physical root of (18)
            lambda_star=s_2_hll/e_hll
            lambda_star_up=lambda_star*g_up_22_theta(i,j,k)

c     calculate state variables for intermediate state
c     (depending on the sign of lambda_star_up):
            if (lambda_star_up .lt. 0.0_rk) then
               p_star=p_if(i,j+1,k,0)+
     &              d_cap_if(0)*
     &              (lambda_plus_max_if-v_up_2_if(0))*
     &              (lambda_star-v_2_if(i,j+1,k,0))
               d_cap_star=d_cap_if(0)*
     &              (lambda_plus_max_if-v_up_2_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_1_star=s_1_if(0)*
     &              (lambda_plus_max_if-v_up_2_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_2_star=(s_2_if(0)*
     &              (lambda_plus_max_if-v_up_2_if(0))+
     &              p_star-p_if(i,j+1,k,0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_3_star=s_3_if(0)*
     &              (lambda_plus_max_if-v_up_2_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               tau_star=(tau_if(0)*
     &              (lambda_plus_max_if-v_up_2_if(0))+
     &              p_star*lambda_star_up-
     &              p_if(i,j+1,k,0)*v_up_2_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
            else
               p_star=p_if(i,j,k,1)+
     &              d_cap_if(1)*
     &              (lambda_minus_min_if-v_up_2_if(1))*
     &              (lambda_star_up-v_2_if(i,j,k,1))
               d_cap_star=d_cap_if(1)*
     &              (lambda_minus_min_if-v_up_2_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_1_star=s_1_if(1)*
     &              (lambda_minus_min_if-v_up_2_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_2_star=(s_2_if(1)*
     &              (lambda_minus_min_if-v_up_2_if(1))+
     &              p_star-p_if(i,j,k,1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_3_star=s_3_if(1)*
     &              (lambda_minus_min_if-v_up_2_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               tau_star=(tau_if(1)*
     &              (lambda_minus_min_if-v_up_2_if(1))+
     &              p_star*lambda_star_up-
     &              p_if(i,j,k,1)*v_up_2_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
            end if
            
            if (lambda_plus_max_if .le. 0.0_rk) then
               p_star=p_if(i,j+1,k,0)
            end if
            
            if (lambda_minus_min_if .ge. 0.0_rk) then
               p_star=p_if(i,j  ,k,1)
            end if
            
            pav_theta (i, j, k) = p_star
            
c     Use HLLC flux for subsonic flow
            lscr1=(lambda_plus_max_if*
     &           lambda_minus_min_if .lt. 0.0_rk) .and.
     &           .not.(lshock(i,j,k) .or. lshock(i,j+1,k))
            if (lscr1) then
               f_num_if (i, j, k, 1) =
     &              d_cap_star*lambda_star_up
               f_num_if (i, j, k, 2) =
     &              s_1_star*lambda_star_up
#if MOMENTUM_EQUATION_MODE==0
               f_num_if (i, j, k, 3) =
     &              s_2_star*lambda_star_up+p_star
#else
               f_num_if (i, j, k, 3) =
     &              s_2_star*lambda_star_up
#endif /* MOMENTUM_EQUATION_MODE==0 */
               f_num_if (i, j, k, 4) =
     &              s_3_star*lambda_star_up
               f_num_if (i, j, k, 5) =
     &              tau_star*lambda_star_up+
     &              p_star*lambda_star_up
#if MOMENTUM_EQUATION_MODE>0
            else
c     Recalculate HLLE momentum flux (without pressure contribution)
               f_num_theta_if (3, 0) =
     &              s_2_if (0) * v_up_2_if (0)
               f_num_theta_if (3, 1) =
     &              s_2_if (1) * v_up_2_if (1) 
               f_num_if (i, j, k, 3) =
     &              (lambda_plus_max_if * f_num_theta_if (3, 1) -
     &              lambda_minus_min_if * f_num_theta_if (3, 0) +
     &              lambda_plus_max_if * lambda_minus_min_if *
     &              (s_2_if (0) - s_2_if (1))) / lambda_diff_if
c     Recalculate effective pressure at interface
               pav_theta (i, j, k) =
     &              (lambda_plus_max_if * p_if (i, j    , k, 1) -
     &              lambda_minus_min_if * p_if (i, j + 1, k, 0)) / 
     &              lambda_diff_if
#endif /* MOMENTUM_EQUATION_MODE>0 */
            end if
            
            
            if (f_num_if (i, j, k, 1) .gt. 0.0_rk) then
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
               f_num_xnu_if (i, j, k, 1 : config%qn) = 
     &              f_num_if (i, j, k, 1) *
     &              xnu_if (i, j, k, 1, 1 : config%qn)
            else
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
               f_num_xnu_if (i, j, k, 1 : config%qn) =
     &              f_num_if (i, j, k, 1) *
     &              xnu_if (i, j + 1, k, 0, 1 : config%qn)
            end if

#ifdef CFC_MHD
            if (lscr1) then 
               v_1_upw = s_1_star / d_cap_star
               v_2_upw = lambda_star_up
               v_3_upw = s_3_star / d_cap_star
            else if (lambda_plus_max_if .le. 0.0_rk) then
               v_1_upw = v_1_if (i, j + 1, k, 0)
               v_2_upw = v_up_2_if (0)
               v_3_upw = v_3_if (i, j + 1, k, 0)
            else
               v_1_upw = v_1_if (i, j, k, 1)
               v_2_upw = v_up_2_if (1)
               v_3_upw = v_3_if (i, j, k, 1)
            end if
c     Compute v* and transform to orthonormal components
            v_1_upw = (alpha_theta (i,j,k) * v_1_upw *
     &           g_up_11_theta (i,j,k) - beta_up_1_theta (i,j,k))
            v_2_upw = (alpha_theta (i,j,k) * v_2_upw) *
     &           r (i)
            v_3_upw = (alpha_theta (i,j,k) * v_3_upw *
     &           g_up_33_theta (i,j,k) -
     &           beta_up_3_theta (i,j,k)) *
     &           r (i) * sinus_theta_if (j)

            if (v_2_upw .gt. 0.0_rk) then
               emf_face (i, j, k, 1, 2) =
     &              v_2_upw * b_3_if (i, j, k, 1) -
     &              v_3_upw * b_2_if (i, j, k, 1)
               emf_face (i, j, k, 2, 2) =
     &              v_1_upw * b_2_if (i, j, k, 1) -
     &              v_2_upw * b_1_if (i, j, k, 1)
            else
               emf_face (i, j, k, 1, 2) =
     &              v_2_upw * b_3_if (i, j+1, k, 0) -
     &              v_3_upw * b_2_if (i, j+1, k, 0)
               emf_face (i, j, k, 2, 2) =
     &              v_1_upw * b_2_if (i, j+1, k, 0) -
     &              v_2_upw * b_1_if (i, j+1, k, 0)
            end if
            emf_face (i, j, k, 1, 2) = emf_face (i, j, k , 1, 2) *
     &           phi_theta (i, j, k) ** 6
            emf_face (i, j, k, 2, 2) = emf_face (i, j, k , 2, 2) *
     &           phi_theta (i, j, k) ** 6
#endif /* CFC_MHD */

         enddo
      enddo

      return

      END SUBROUTINE calculate_flux_theta

c     ==================================================================

#else /* MHDFLX */
      
c     ==================================================================

      SUBROUTINE calculate_flux_theta

c     ==================================================================

c     Computes HLLD/HLLE fluxes in theta-direction (Newtonian MHD version)

c     ------------------------------------------------------------------
c     Different from the pure hydro version, all velocities are
c     transformed to an orthonormal frame.
c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE parameters_cfc
      USE interface_hydro_cfc, b_1_if_glob => b_1_if,
     &     b_2_if_glob => b_2_if, b_3_if_glob => b_3_if,
     &     psi_if_glob => psi_if
      USE grid_cfc
      USE fluxes_cfc
      USE metric_cfc
      USE perm_aux_cfc
      USE shock_cfc

c#ifdef CONVECTION
c      USE hydro_primitives_cfc, ONLY: rho,eps,p,gamm,
c     &     c_sound_squared,xnnu,v_1
      USE hydro_areas_mod
      USE configure
c#endif
      USE hydro_primitives_cfc

      use configure

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk, l, ll

      real (kind=rk) :: f_num_theta_if (1 : 8, 0 : 1)
      real (kind=rk) :: f_num_xnu_theta_if (1 : config%qn, 0 : 1)

      real (kind=rk) :: d_cap_if (0 : 1), s_1_if (0 : 1),
     &     s_2_if (0 : 1), s_3_if (0 : 1), tau_if (0 : 1)
      real (kind=rk) :: d_cap_xnu_if (1 : config%qn, 0 : 1)
      real (kind=rk) :: b_1_if (0 : 1), b_2_if (0 : 1),
     &     b_3_if (0 : 1), cs2_if (0 : 1), bsq_if (0 : 1),
     &     v_or_1_if (0 : 1), v_or_2_if (0 : 1), v_or_3_if (0 : 1),
     &     psi_if (0 : 1)
      
      real (kind=rk) :: lambda_minus_if (0 : 1)
      real (kind=rk) :: lambda_plus_if (0 : 1)

      real (kind=rk) :: lambda_minus_min_if, lambda_plus_max_if,
     &     lambda_diff_if, lambda_b_if
      real (kind=rk) :: c_sound_if (0 : 1)

      real (kind=rk) :: e_hll,s_2_hll,lambda_star,
     &     lambda_a_minus, lambda_a_plus,
     &     p_star,d_cap_star,v_1_star,v_2_star,v_3_star,
     &     b_1_star,b_2_star,b_3_star,psi_star,
     &     tau_star,scr1,sw1,vcln,vclnl,vclnr
      real (kind=rk) ::
     &     dcsl, dcsr, psl, psr, esl, esr,
     &     v1sl, v3sl, v1sr, v3sr,
     &     b1sl, b3sl, b1sr, b3sr,
     &     swl, swr, tmp1, tmp2
      real (kind=rk) :: slft, srght, alft, argt
      real (kind=rk) :: v_1_upw, v_2_upw, v_3_upw
      logical :: lscr1, lscr2

      integer :: iarea
      integer (kind=ik) :: ixi,ixf,iox,iyi,iyf,ioy,izi,izf,ioz,isd
      integer :: filter(1:m)


      if (config%nsdim .gt. 1) then
         do iarea=1,areas%are_nu

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
            
            if (ioy .eq.  config%qy) then
               filter (ixi:ixf) = 0
            else
               filter (ixi:ixf) = 1
            endif
         enddo
      end if
      
      
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,l,ll)
#endif
      do jk = 1, (n_loc + 1) * o_loc

         k = int((jk + (n_loc + 1) - 1) / (n_loc + 1))
         j = (n_s - 2) + (jk - (k - 1) * (n_loc + 1))
         k = k + o_s - 1

         if (j .eq. 0 .or. j .eq. n) then
            sw1 = 0.0_rk
         else
            sw1 = 1.0_rk
         endif

         do i = 1, m

            lscr2=.true.
            
ccc   ==================================================================
ccc   left interface
               
            v_or_1_if (0) =
     &           sqrt (g_up_11_theta (i, j, k)) * v_1_if (i, j+1, k, 0)
            v_or_2_if (0) =
     &           sqrt (g_up_22_theta (i, j, k)) * v_2_if (i, j+1, k, 0)
            v_or_3_if (0) =
     &           sqrt (g_up_33_theta (i, j, k)) * v_3_if (i, j+1, k, 0)
ccc   calculate interface value of v^i = g^k v_i
            
            b_1_if (0) = b_1_if_glob (i, j + 1, k, 0)
            b_2_if (0) = b_2_if_glob (i, j + 1, k, 0)
            b_3_if (0) = b_3_if_glob (i, j + 1, k, 0)
            bsq_if (0) =
     &           b_1_if (0) ** 2 + b_2_if (0) ** 2 + b_3_if (0) ** 2
            psi_if (0) = psi_if_glob (i, j + 1, k, 0)
               
            d_cap_if (0) = rho_if (i, j + 1, k, 0)
            s_1_if (0) = d_cap_if (0) * v_or_1_if (0)
            s_2_if (0) = d_cap_if (0) * v_or_2_if (0)
            s_3_if (0) = d_cap_if (0) * v_or_3_if (0)
            tau_if (0) = d_cap_if (0) *
     &           h_if (i, j + 1, k, 0) - p_if (i, j + 1, k, 0)
     &           + 0.5_rk * bsq_if (0)

            d_cap_xnu_if (1 : config%qn, 0) = d_cap_if (0) *
     &           xnu_if (i, j + 1, k, 0, 1 : config%qn)

            c_sound_if (0) =
     &           sqrt (c_sound_squared_if (i, j + 1, k, 0))! +
c     &           v_or_3_if (0) **2)
ccc   calculate sound speed c_s
            cs2_if (0) = c_sound_squared_if (i, j + 1, k, 0)! +
c     &           v_or_3_if (0) **2

ccc   calculate eigenvalues --------------------------------------------
            tmp1 = cs2_if (0) * d_cap_if (0) + bsq_if (0)
            tmp2 = sqrt ((tmp1 + 
     &           sqrt (tmp1 ** 2 - 4.0_rk * cs2_if (0) * d_cap_if (0) *
     &           0*b_2_if (0) ** 2)) / (2.0_rk * d_cap_if (0)))
            lambda_minus_if (0) = v_or_2_if (0) - tmp2
ccc   lambda_-
            lambda_plus_if  (0) = v_or_2_if (0) + tmp2
ccc   lambda_+


ccc   ==================================================================
ccc   right interface
            
            v_or_1_if (1) =
     &           sqrt (g_up_11_theta (i, j, k)) * v_1_if (i, j, k, 1)
            v_or_2_if (1) =
     &           sqrt (g_up_22_theta (i, j, k)) * v_2_if (i, j, k, 1)
            v_or_3_if (1) =
     &           sqrt (g_up_33_theta (i, j, k)) * v_3_if (i, j, k, 1)
ccc   calculate interface value of v^i = g^k v_i
            
            b_1_if (1) = b_1_if_glob (i, j , k, 1)
            b_2_if (1) = b_2_if_glob (i, j , k, 1)
            b_3_if (1) = b_3_if_glob (i, j , k, 1)
            bsq_if (1) =
     &           b_1_if (1) ** 2 + b_2_if (1) ** 2 + b_3_if (1) ** 2
            psi_if (1) = psi_if_glob (i, j, k, 1)
               
            d_cap_if (1) = rho_if (i, j, k, 1)
            s_1_if (1) = d_cap_if (1) * v_or_1_if (1)
            s_2_if (1) = d_cap_if (1) * v_or_2_if (1)
            s_3_if (1) = d_cap_if (1) * v_or_3_if (1)
            tau_if (1) = d_cap_if (1) *
     &           h_if (i, j, k, 1) - p_if (i, j, k, 1)
     &           + 0.5_rk * bsq_if (1)

            d_cap_xnu_if (1 : config%qn, 1) = rho_if (i, j, k, 1) *
     &           xnu_if (i, j, k, 1, 1 : config%qn)
            
            c_sound_if (1) = sqrt (c_sound_squared_if (i, j, k, 1))! +
c     &           v_or_3_if (1) **2)
ccc   calculate sound speed c_s
            cs2_if (1) = c_sound_squared_if (i, j, k, 1)! +
c     &           v_or_3_if (1) **2

ccc   calculate eigenvalues --------------------------------------------
            tmp1 = cs2_if (1) * d_cap_if (1) + bsq_if (1)
            tmp2 = sqrt ((tmp1 + 
     &           sqrt (tmp1 ** 2 - 4.0_rk * cs2_if (1) * d_cap_if (1) *
     &           0*b_2_if (1) ** 2)) / (2.0_rk * d_cap_if (1)))
            lambda_minus_if (1) = v_or_2_if (1) - tmp2
ccc   lambda_-
            lambda_plus_if  (1) = v_or_2_if (1) + tmp2
ccc   lambda_+

  
ccc   ==================================================================
ccc   calculate maximum of left and right interface eigenvalues
ccc   and their difference
            lambda_minus_min_if = min (lambda_minus_if (0),
     &           lambda_minus_if (1), 0.0_rk)
            lambda_plus_max_if = max (lambda_plus_if (0),
     &           lambda_plus_if (1), 0.0_rk)
#if MHD_SOLVER==HLLE
            if (filter (i) .eq. 0) then
c     We must reduce the diffusivity in the spherical core
               scr1 = r (i) / delta_t /
     &              MAX(lambda_plus_max_if, -lambda_minus_min_if)
               lambda_plus_max_if = scr1 * lambda_plus_max_if
               lambda_minus_min_if = scr1 * lambda_minus_min_if
            end if
#endif /* MHD_SOLVER==HLLE */
            
            lambda_diff_if = lambda_plus_max_if - lambda_minus_min_if
            lambda_b_if = max (abs (v_or_2_if (0)), abs (v_or_2_if (1)))


ccc   ==================================================================
c     Treat 2 * 2 subsytem first

c advective flux of psi           
            f_num_if (i, j, k, 11) = 0.25_rk *
     &          ((psi_if (1) + psi_if (0)) * 
     &           (v_or_2_if (0) + v_or_2_if (1)) -
     &           (psi_if (1) - psi_if (0)) *
     &            abs (v_or_2_if (0) - v_or_2_if (1)))

            vcln = 0.25_rk * min (r (i) * delta_theta, delta_r (i)) /
     &           delta_t !cleaning speed
c            vcln = min (vcln, 0.5_rk *
c     &          (sqrt (cs2_if (0) + bsq_if (0) / d_cap_if (0)) +
c     &           sqrt (cs2_if (1) + bsq_if (1) / d_cap_if (1))))
c            vcln = 0.0_rk
            vclnr = min(vcln, sqrt(cs2_if(0) + bsq_if(0) / d_cap_if(0)))
            vclnl = min(vcln, sqrt(cs2_if(1) + bsq_if(1) / d_cap_if(1)))
c            vcln = 0.5_rk * (vclnl + vclnl)
            vcln = max (vclnl, vclnr)
            
c            srgt = 0.5_rk * (b_2_if (0) + psi_if (0))
            slft = 0.5_rk * (b_2_if (1) + psi_if (1))
            argt = 0.5_rk * (b_2_if (0) - psi_if (0))
c            alft = 0.5_rk * (b_2_if (1) - psi_if (1))

c            if (vcln + 0.5_rk * (v_or_2_if (0) + v_or_2_if (1))
c     &           .lt. 0.0_rk) then
c               b_2_star = b_2_if (0)
c               psi_star = psi_if (0)
c            else if (- vcln + 0.5_rk * (v_or_2_if (0) + v_or_2_if (1))
c     &              .gt. 0.0_rk) then
c               b_2_star = b_2_if (1)
cc               psi_star = psi_if (0)
c            else            
            b_2_star = 0.5_rk *
     &           (b_2_if (0) + b_2_if(1) -
     &           (psi_if (0) - psi_if (1)) * sw1)
            psi_star = 0.5_rk * 
     &           (psi_if (0) + psi_if (1) -
     &           (b_2_if (0) - b_2_if (1)) * sw1)
c            psi_star =
c     &           (vclnl * slft - vclnr * argt) / vcln
c            end if
c            b_2_if (0) = b_2_star
c            b_2_if (1) = b_2_star


            
ccc   ==================================================================
ccc   calculate flux at left and right interface            
            do l = 0, 1
               f_num_theta_if (1, l) = d_cap_if (l) * v_or_2_if (l)
               f_num_theta_if (2, l) = s_1_if (l) * v_or_2_if (l) -
     &              b_2_if (l) * b_1_if (l)
               f_num_theta_if (3, l) = s_2_if (l) * v_or_2_if (l) + 
     &              p_if (i, j + 1 - l, k, l) + 0.5_rk * bsq_if (l) -
     &              b_2_if (l) ** 2
               f_num_theta_if (4, l) = s_3_if (l) * v_or_2_if (l) -
     &              b_2_if (l) * b_3_if (l)
               f_num_theta_if (5, l) = (tau_if (l) +
     &              p_if (i, j + 1 - l, k, l) + 0.5_rk * bsq_if (l)) *
     &              v_or_2_if (l) -  b_2_if (l) *
     &              (v_or_1_if (l) * b_1_if (l) +
     &              v_or_2_if (l) * b_2_if (l) +
     &              v_or_3_if (l) * b_3_if (l))
               f_num_theta_if (6, l) = - e_3_if (i, j + 1 - l, k, l)
               f_num_theta_if (7, l) =   0.0_rk
               f_num_theta_if (8, l) =   e_1_if (i, j + 1 - l, k, l)

               f_num_xnu_theta_if (1 : config%qn , l) =
     &              d_cap_xnu_if (1 : config%qn, l) *
     &              v_or_2_if (l)
            enddo



ccc   ==================================================================
ccc   calculate numerical fluxes at interface
            f_num_if (i, j, k, 1) =
     &           (lambda_plus_max_if * f_num_theta_if (1, 1) -
     &           lambda_minus_min_if * f_num_theta_if (1, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (d_cap_if (0) - d_cap_if (1))) / lambda_diff_if

            f_num_if (i, j, k, 2) =
     &           (lambda_plus_max_if * f_num_theta_if (2, 1) -
     &           lambda_minus_min_if * f_num_theta_if (2, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_1_if (0) - s_1_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 3) =
     &           (lambda_plus_max_if * f_num_theta_if (3, 1) -
     &           lambda_minus_min_if * f_num_theta_if (3, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_2_if (0) - s_2_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 4) =
     &           (lambda_plus_max_if * f_num_theta_if (4, 1) -
     &           lambda_minus_min_if * f_num_theta_if (4, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_3_if (0) - s_3_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 5) =
     &           (lambda_plus_max_if * f_num_theta_if (5, 1) -
     &           lambda_minus_min_if * f_num_theta_if (5, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (tau_if (0) - tau_if (1))) / lambda_diff_if

            f_num_if (i, j, k, 6) =
     &           (lambda_plus_max_if * f_num_theta_if (6, 1) -
     &           lambda_minus_min_if * f_num_theta_if (6, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if * filter (i) *
     &           (b_1_if (0) - b_1_if (1))) / lambda_diff_if

            f_num_if (i, j, k, 7) = psi_star * vcln

            f_num_if (i, j, k, 8) =
     &           (lambda_plus_max_if * f_num_theta_if (8, 1) -
     &           lambda_minus_min_if * f_num_theta_if (8, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if * filter (i) *
     &           (b_3_if (0) - b_3_if (1))) / lambda_diff_if

            f_num_if (i, j, k, 9) = b_2_star

            f_num_if (i, j, k, 10) =
     &           (lambda_plus_max_if * d_cap_if (1) * 
     &           eps_if (i, j, k, 1) * v_or_2_if (1) -
     &           lambda_minus_min_if * d_cap_if (0) *
     &           eps_if (i, j + 1, k, 0) * v_or_2_if (0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (d_cap_if (0) * eps_if (i, j + 1, k, 0) -
     &           d_cap_if (1) * eps_if (i, j, k, 1))) / lambda_diff_if

#ifdef CONSTRAINED
            emf_face (i, j, k, 1, 2) = ! (v \times B)_1 = v_2 B_3 - v_3 B_2
     &           (lambda_plus_max_if * f_num_theta_if (7, 1) -
     &           lambda_minus_min_if * f_num_theta_if (7, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if * filter (i) *
     &           (b_3_if (0) - b_3_if (1))) / lambda_diff_if
            
            emf_face (i, j, k, 2, 2) = ! (v \times B)_3 = v_1 B_2 - v_2 B_1
     &           - (lambda_plus_max_if * f_num_theta_if (6, 1) -
     &           lambda_minus_min_if * f_num_theta_if (6, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if * filter (i) *
     &           (b_1_if (0) - b_1_if (1))) / lambda_diff_if
#endif /* CONSTRAINED */
            
            if (f_num_if (i, j, k, 1) .gt. 0.0_rk) then
               v_2_star = f_num_if (i, j, k, 1) / d_cap_if (1)
            else
               v_2_star = f_num_if (i, j, k, 1) / d_cap_if (0)
            endif
            
c     calculate HLL state vector (only 2 components needed)
            e_hll = (lambda_plus_max_if * d_cap_if (0) -
     &           lambda_minus_min_if * d_cap_if (1) +
     &           f_num_theta_if (1, 1) - f_num_theta_if (1, 0)) /
     &           lambda_diff_if
            s_2_hll = (lambda_plus_max_if * s_2_if (0) -
     &           lambda_minus_min_if * s_2_if (1) +
     &           f_num_theta_if (3, 1) - f_num_theta_if (3, 0)) /
     &           lambda_diff_if
            lambda_star=s_2_hll/e_hll!*filter(i)

c     calculate state variables for intermediate states

c     right star state
            psr  =p_if(i,j+1,k,0)+0.5_rk*bsq_if(0)+
     &           d_cap_if(0)*
     &           (lambda_plus_max_if-v_or_2_if(0))*
     &           (lambda_star-v_or_2_if(0))
            dcsr = d_cap_if(0)*
     &           (lambda_plus_max_if-v_or_2_if(0))/
     &           (lambda_plus_max_if-lambda_star)            
            scr1 = dcsr *
     &           (lambda_plus_max_if - lambda_star) ** 2 -
     &           b_2_if (0) ** 2
            if (scr1 .ge. 1e-5_rk * psr) then
               v1sr = v_or_1_if (0) -
     &              b_2_if (0) * b_1_if (0) *
     &              (lambda_star - v_or_2_if (0)) / scr1
               v3sr = v_or_3_if (0) -
     &              b_2_if (0) * b_3_if (0) *
     &              (lambda_star - v_or_2_if (0)) / scr1
               b1sr = b_1_if (0) * (d_cap_if (0) *
     &              (lambda_plus_max_if - v_or_2_if (0)) ** 2 -
     &              b_2_if (0) ** 2) / scr1
               b3sr = b_3_if (0) * (d_cap_if (0) *
     &              (lambda_plus_max_if - v_or_2_if (0)) ** 2 -
     &              b_2_if (0) ** 2) / scr1
            else
               lscr2 = .false.
               v1sr = v_or_1_if (0)
               v3sr = v_or_3_if (0)
               b1sr = b_1_if (0)
               b3sr = b_3_if (0)               
            end if
            esr  = (tau_if(0)*
     &           (lambda_plus_max_if-v_or_2_if(0))+
     &           psr*lambda_star-
     &           (p_if(i,j+1,k,0)+0.5_rk*bsq_if(0))*v_or_2_if(0)+
     &           b_2_if (0) * (
     &           v_or_1_if (0) * b_1_if (0) +
     &           v_or_2_if (0) * b_2_if (0) +
     &           v_or_3_if (0) * b_3_if (0) -
     &           lambda_star * b_2_if (0) - v1sr * b1sr - v3sr * b3sr))/
     &           (lambda_plus_max_if-lambda_star)
            lambda_a_plus = lambda_star + abs(b_2_if(0))/sqrt (dcsr)
            
c     left star state
            psl  =p_if(i,j,k,1)+0.5_rk*bsq_if(1)+
     &           d_cap_if(1)*
     &           (lambda_minus_min_if-v_or_2_if(1))*
     &           (lambda_star-v_or_2_if(1))
            dcsl = d_cap_if(1)*
     &           (lambda_minus_min_if-v_or_2_if(1))/
     &           (lambda_minus_min_if-lambda_star)
            scr1 = dcsr *
     &           (lambda_minus_min_if - lambda_star) ** 2 -
     &           b_2_if (1) ** 2
            if (scr1 .ge. 1e-5_rk * psl) then
               v1sl = v_or_1_if (1) -
     &              b_2_if (1) * b_1_if (1) *
     &              (lambda_star - v_or_2_if (1)) / scr1
               v3sl = v_or_3_if (1) -
     &              b_2_if (1) * b_3_if (1) *
     &              (lambda_star - v_or_2_if (1)) / scr1
               b1sl = b_1_if (1) * (d_cap_if (1) *
     &              (lambda_minus_min_if - v_or_2_if (1)) ** 2 -
     &              b_2_if (1) ** 2) / scr1
               b3sl = b_3_if (1) * (d_cap_if (1) *
     &              (lambda_minus_min_if - v_or_2_if (1)) ** 2 -
     &              b_2_if (1) ** 2) / scr1
            else
               lscr2 = .false.
               v1sl = v_or_1_if (1)
               v3sl = v_or_3_if (1)
               b1sl = b_1_if (1)
               b3sl = b_3_if (1)               
            end if
            esl  = (tau_if(1)*
     &           (lambda_minus_min_if-v_or_2_if(1))+
     &           psl*lambda_star-
     &           (p_if(i,j,k,1)+0.5_rk*bsq_if(1))*v_or_2_if(1)+
     &           b_2_if (1) * (
     &           v_or_1_if (1) * b_1_if (1) +
     &           v_or_2_if (1) * b_2_if (1) +
     &           v_or_3_if (1) * b_3_if (1) -
     &           lambda_star * b_2_if (1) - v1sl * b1sl - v3sl * b3sl))/
     &           (lambda_minus_min_if-lambda_star)
            lambda_a_minus = lambda_star - abs(b_2_if(1))/sqrt (dcsl)

            swl = 0.5_rk + sign (0.5_rk,  lambda_star)
            if (lambda_star .eq. 0.0_rk) swl = 0.5_rk
            swr = 1.0_rk - swl
            p_star     = psl * swl + psr * swr
            d_cap_star = dcsl * swl + dcsr * swr
            v_1_star   = v1sl * swl + v1sr * swr
            v_3_star   = v3sl * swl + v3sr * swr
            b_1_star   = b1sl * swl + b1sr * swr
c            b_2_star   = 0.5_rk * (b_2_if (0) + b_2_if (1))
            b_3_star   = b3sl * swl + b3sr * swr
            tau_star   = esl  * swl + esr * swr
c     If  we're in one of the star states, this is already the state
c     that we need. If we're between the Alfven waves, we need
c     compute the double star state.


#if MHD_SOLVER == HLLD            
            if (lambda_a_minus * lambda_a_plus .le. 0.0_rk .and.
     &           lambda_a_minus .gt. lambda_minus_min_if .and.
     &           lambda_a_plus  .lt. lambda_plus_max_if) then
               dcsl = sqrt (dcsl)
               dcsr = sqrt (dcsr)
               scr1 = dcsl+dcsr
               v_1_star = (dcsl*v1sl+dcsr*v1sr+
     &              sign(b1sr-b1sl,b_2_star))/scr1
               v_3_star = (dcsl*v3sl+dcsr*v3sr+
     &              sign(b3sr-b3sl,b_2_star))/scr1
               b_1_star = (dcsl*b1sr+dcsr*b1sl+
     &              dcsl*dcsr*sign(v1sr-v1sl,b_2_star))/scr1
               b_3_star = (dcsl*b3sr+dcsr*b3sl+
     &              dcsl*dcsr*sign(v3sr-v3sl,b_2_star))/scr1
               tau_star = tau_star +
     &              swr  * dcsr*
     &              sign(v1sr*b1sr+v3sr*b3sr- !v_2 B_2 cancels
     &              v_1_star*b_1_star-v_3_star*b_3_star,b_2_star) -
     &              swl  * dcsl*
     &              sign(v1sl*b1sl+v3sl*b3sl- !v_2 B_2 cancels
     &              v_1_star*b_1_star-v_3_star*b_3_star,b_2_star)
            end if
#endif /* MHD_SOLVER == HLLD */
            
            pav_theta (i, j, k) = p_star

            if (min (lambda_a_minus-lambda_minus_min_if,
     &           lambda_plus_max_if-lambda_a_plus)
     &       .lt. 5e-2_rk*lambda_diff_if) lscr2 = .false.

c     Use HLLD flux for subsonic flow
            lscr1=(lambda_plus_max_if*
     &           lambda_minus_min_if .lt. 0.0_rk).and.
     &           .not.(lshock(i,j,k).or.lshock(i,j+1,k))
            if (filter (i) .eq. 0) lscr1 = .false.
#if MHD_SOLVER == HLLE
            lscr1 = .false.
#endif            
            if (lscr1.and.lscr2) then
               f_num_if (i, j, k, 1) =
     &              d_cap_star*lambda_star
               f_num_if (i, j, k, 2) =
     &              d_cap_star*v_1_star*lambda_star-b_1_star*b_2_star
#if MOMENTUM_EQUATION_MODE==0
               f_num_if (i, j, k, 3) =
     &              d_cap_star*lambda_star**2+p_star-b_2_star**2
#else
               f_num_if (i, j, k, 3) =
     &              d_cap_star*lambda_star**2-b_2_star**2
#endif /* MOMENTUM_EQUATION_MODE==0 */
               f_num_if (i, j, k, 4) =
     &              d_cap_star*v_3_star*lambda_star-b_3_star*b_2_star
               f_num_if (i, j, k, 5) =
     &              tau_star*lambda_star+
     &              p_star*lambda_star - b_2_star *
     &              (lambda_star*b_2_star+
     &              v_1_star*b_1_star+v_3_star*b_3_star)
               
               v_1_upw = v_1_star
               v_2_upw = lambda_star * filter (i)
               v_3_upw = v_3_star

               f_num_if (i, j, k, 6) =
     &              v_2_upw * b_1_star - v_1_upw * b_2_star
               f_num_if (i, j, k, 8) =
     &              v_2_upw * b_3_star - v_3_upw * b_2_star
               
#ifdef CONSTRAINED               
               emf_face (i, j, k, 1, 2) =
     &              v_2_upw * b_3_star - 0*v_3_upw * b_2_star
               emf_face (i, j, k, 2, 2) =
     &              0*v_1_upw * b_2_star - 1*v_2_upw * b_1_star
#if DIMN==2               
c     Upwind selection for v_3 * B_1 term -- works only in 2D
               if (v_2_upw .gt. 0.0_rk) then
                  emf_face (i, j, k, 1, 2) = emf_face (i, j, k, 1, 2)
     &                 - v_3 (i, j, k) * b_2 (i, j, k)
               else
                  emf_face (i, j, k, 1, 2) = emf_face (i, j, k, 1, 2)
     &                 - v_3 (i, j + 1, k) * b_2 (i, j + 1, k)
               end if
#endif /* DIMN==2*/
c     diffusive term for stabilising constraint terms
c               emf_face (i, j, k, 1, 2) = emf_face (i, j, k, 1, 2) +
c     &              2.0_rk *
c     &              max(abs (v_or_3_if (1)), abs(v_or_3_if (0))) *
c     &              (b_3_if (1) - b_3_if (0))
#endif /* CONSTRAINED */
            f_num_if (i, j, k, 6) = f_num_if (i, j, k, 6) +
     &           0.5_rk * lambda_diff_if * filter (i) *
     &           (b_1_if (1) - b_1_if (0))
            f_num_if (i, j, k, 8) = f_num_if (i, j, k, 8) +
     &           0.5_rk * lambda_diff_if * filter (i) *
     &              (b_3_if (1) - b_3_if (0))
               
               v_2_star = lambda_star
               
#if MOMENTUM_EQUATION_MODE>0
            else
c     Recalculate HLLE momentum flux (without pressure contribution)
               f_num_theta_if (3, 0) = s_2_if (0) * v_or_2_if (0) -
     &              b_2_if (0) ** 2
               f_num_theta_if (3, 1) = s_2_if (1) * v_or_2_if (1) -
     &              b_2_if (1) ** 2
               f_num_if (i, j, k, 3) =
     &              (lambda_plus_max_if * f_num_theta_if (3, 1) -
     &              lambda_minus_min_if * f_num_theta_if (3, 0) +
     &              lambda_plus_max_if * lambda_minus_min_if *
     &              (s_2_if (0) - s_2_if (1))) / lambda_diff_if
c     Recalculate effective pressure at interface
               pav_theta (i, j, k) =
     &              (lambda_plus_max_if * (p_if (i, j    , k, 1) +
     &              0.5_rk * bsq_if (1)) -
     &              lambda_minus_min_if * (p_if (i, j + 1, k, 0) +
     &              0.5_rk * bsq_if (0))) / 
     &              lambda_diff_if
#endif /* MOMENTUM_EQUATION_MODE>0 */
            end if

            if (filter (i) .eq. 0) then
               scr1 = MIN (1.0_rk, 0.5_rk * r (i) * delta_theta /
     &              delta_t /
     &              MAX(lambda_plus_max_if, -lambda_minus_min_if))
               lambda_plus_max_if = scr1 * lambda_plus_max_if
               lambda_minus_min_if = scr1 * lambda_minus_min_if
               lambda_diff_if = lambda_plus_max_if - lambda_minus_min_if

               f_num_if (i, j, k, 6) =
     &              (0*lambda_plus_max_if * f_num_theta_if (6, 1) -
     &              0*lambda_minus_min_if * f_num_theta_if (6, 0) +
     &              lambda_plus_max_if * lambda_minus_min_if *
     &              (b_1_if (0) - b_1_if (1))) / lambda_diff_if

               f_num_if (i, j, k, 8) =
     &              (0*lambda_plus_max_if * f_num_theta_if (8, 1) -
     &              0*lambda_minus_min_if * f_num_theta_if (8, 0) +
     &              lambda_plus_max_if * lambda_minus_min_if *
     &              (b_3_if (0) - b_3_if (1))) / lambda_diff_if

#ifdef CONSTRAINED
               emf_face (i, j, k, 1, 2) = ! (v \times B)_1 = v_2 B_3 - v_3 B_2
     &              (lambda_plus_max_if * f_num_theta_if (8, 1) -
     &              lambda_minus_min_if * f_num_theta_if (8, 0) +
     &              lambda_plus_max_if * lambda_minus_min_if *
     &              (b_3_if (0) - b_3_if (1))) / lambda_diff_if
               
               emf_face (i, j, k, 2, 2) = ! (v \times B)_3 = v_1 B_2 - v_2 B_1
     &              - (lambda_plus_max_if * f_num_theta_if (6, 1) -
     &              lambda_minus_min_if * f_num_theta_if (6, 0) +
     &              lambda_plus_max_if * lambda_minus_min_if *
     &              (b_1_if (0) - b_1_if (1))) / lambda_diff_if
#endif /* CONSTRAINED */               
            end if

c add energy flux in divergence cleaning mode
            tmp1 = 0.5_rk * (b_2_if (0) + b_2_if (1))
            tmp2 = 0.5_rk * (psi_if (0) + psi_if (1))
            f_num_if (i, j, k, 5) = f_num_if (i, j, k, 5) +
     &           vcln  * (tmp1 * tmp2 +
     &           0.5_rk * tmp1 * (b_2_if (0) - b_2_if (1)) +
     &           0.5_rk * tmp2 * (psi_if (0) - psi_if (1)))           
c     &           f_num_if (i, j, k, 7) * b_2_star +
c     &           f_num_if (i, j, k, 11) * psi_star
            
c     Transform to contravriant/mixed flux components in
c     spherical polar reference metric:
            f_num_if (i, j, k, 1) = f_num_if (i, j, k, 1) * r_inv (i)
            f_num_if (i, j, k, 2) = f_num_if (i, j, k, 2) * r_inv (i)
            !f_num_if (i, j, k, 3) = f_num_if (i, j, k, 3)
            f_num_if (i, j, k, 4) = f_num_if (i, j, k, 4) *
     &           sinus_theta_if (j)
            f_num_if (i, j, k, 5) = f_num_if (i, j, k, 5) * r_inv (i)
            f_num_if (i, j, k, 6) = f_num_if (i, j, k, 6) * r_inv (i)
c            f_num_if (i, j, k, 7) = f_num_if (i, j, k, 7) * rr_inv (i)
            f_num_if (i, j, k, 8) = f_num_if (i, j, k, 8) * rr_inv (i) *
     &           sinus_theta_inv_if (j)
            f_num_if (i, j, k, 9) = f_num_if (i, j, k, 9) * r_inv (i)
            f_num_if (i, j, k, 10) = f_num_if (i, j, k, 10) * r_inv (i)
            f_num_if (i, j, k, 11) = f_num_if (i, j, k, 11) * r_inv (i)

            if (f_num_if(i,j,k,1) .gt. 0.0_rk) then 
               if (lscr1) f_num_if (i, j, k, 10) = 
     &              f_num_if(i, j, k, 1) *
     &              eps_if (i, j, k, 1)
               f_num_xnu_if (i, j, k, 1 : config%qn) = 
     &              f_num_if(i, j, k, 1) *
     &              xnu_if (i, j, k, 1, 1 : config%qn)
            else
               if (lscr1) f_num_if (i, j, k, 10) = 
     &              f_num_if(i, j, k, 1) *
     &              eps_if (i, j + 1, k, 0)
               f_num_xnu_if (i, j, k, 1 : config%qn) = 
     &              f_num_if(i, j, k, 1) *
     &              xnu_if (i, j + 1, k, 0, 1 : config%qn)            
            end if
         enddo
      enddo

      return

      END SUBROUTINE calculate_flux_theta

c     ==================================================================

#endif /* MHDFLX */

      
#endif /* CFC_TRANSPORT2 */



      END MODULE gr_flux_theta

c     ==================================================================
