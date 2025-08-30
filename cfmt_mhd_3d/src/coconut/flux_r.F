#if defined (CONVECTION_1D) 
#define CONVECTION
#endif

#if defined (CFC_MHD) && !defined(PASSIVE_MHD)
#define MHDFLX
#if !(MHD_SOLVER == HLLE) && !(MHD_SOLVER == HLLC) && !(MHD_SOLVER == HLLD)
#error No valid MHD Riemann solver selected.
#endif       
#endif /* CFC_MHD and not PASSIVE_MHD */
      
c     ==================================================================

      MODULE gr_flux_r

c     ==================================================================

      IMPLICIT NONE

      CONTAINS




#ifdef CFC_TRANSPORT2

c     ==================================================================

      SUBROUTINE calculate_flux_r

c     ==================================================================

c     Computes HLLC/HLLE fluxes in r-direction (Relativistic version)

c     ------------------------------------------------------------------

      
ccc   This subroutine computes the numerical fluxes
ccc   in the radial direction

      USE precision
      USE phycon, ONLY: pc_ggeo

      USE size_cfc
      USE parameters_cfc
      USE interface_hydro_cfc
      USE grid_cfc
      USE fluxes_cfc
      USE metric_cfc
      USE perm_aux_cfc
      USE shock_cfc

#ifdef CONVECTION
      USE hydro_primitives_cfc, ONLY: rho,eps,p,gamm,
     &     c_sound_squared,xnnu,v_1
      USE hydro_areas_mod
#endif

      use configure

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk, l, ll

      real (kind=rk) :: f_num_r_if (1 : 5, 0 : 1)
      real (kind=rk) :: f_num_xnu_r_if (1 : config%qn, 0 : 1)

      real (kind=rk) :: w_if (0 : 1)

      real (kind=rk) :: d_cap_if (0 : 1), s_1_if (0 : 1),
     &     s_2_if (0 : 1), s_3_if (0 : 1), tau_if (0 : 1)
      real (kind=rk) :: d_cap_xnu_if (1 : config%qn, 0 : 1)

      real (kind=rk) :: v_up_1_if (0 : 1)
      real (kind=rk) :: v_up_1_hat_if (0 : 1)

      real (kind=rk) :: lambda_minus_if (0 : 1)
      real (kind=rk) :: lambda_plus_if (0 : 1)

      real (kind=rk) :: lambda_minus_min_if, lambda_plus_max_if,
     &     lambda_diff_if
      real (kind=rk) :: c_sound_if (0 : 1)

      real (kind=rk) :: tmp_1, tmp_2, tmp_3, rho_h_ww_if

      real (kind=rk) :: e_hll,s_1_hll,lambda_star,lambda_star_up,
     &     p_star,d_cap_star,s_1_star,s_2_star,s_3_star,
     &     tau_star,scr1

#ifdef CFC_MHD
      real (kind=rk) :: v_1_upw, v_2_upw, v_3_upw
#endif /* CFC_MHD */
      logical :: lscr1

#ifdef CONVECTION
      real (kind=rk) :: cledoux,lmix,vmix,glocal,rhotot,deltap
      integer :: i_sonic,filter(0:m),iarea
      integer (kind=ik) :: ixi,ixf,iox,iyi,iyf,ioy,izi,izf,ioz,isd
#endif


#ifdef CONVECTION
c      if (config%nsdim .eq. 1) then
      i_sonic=m-1
#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,filter)
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         do i=2,m
            if (v_1(i,j,k)**2 .gt. c_sound_squared(i,j,k))
     &           i_sonic=min(i_sonic,i)
c            if (rho(i,j,k) .lt. 5e12_rk * pc_ggeo)
c     &           i_sonic=min(i_sonic,i)
         end do
      end do
      filter(0)=0
      filter(1:i_sonic-1)=1
      filter(i_sonic:m)=0
c     endif
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
            
            if (ioy .ne.  config%qy) then
               filter (max(ixi-1,1):min(ixf+1,config%qx)) = 0
            endif
         enddo
      end if

#endif /* CONVECTION */

#if  defined(OPENMP_CFC) && ( defined(OPENMP_CFC) && defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,l,ll)
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 0, m

ccc   ==================================================================
ccc   left interface
               
            c_sound_if (0) =
     &           sqrt (c_sound_squared_if (i + 1, j, k, 0))

ccc   calculate sound speed c_s
            
            w_if (0) = 1.0_rk /
     &           sqrt (1.0_rk - v_squared_if (i + 1, j, k, 0))
ccc   calculate Lorentz factor W

            v_up_1_if (0) =
     &           g_up_11_r (i, j, k) * v_1_if (i + 1, j, k, 0)
ccc   calculate interface value of v^i = g^k v_i

            v_up_1_hat_if (0) = v_up_1_if (0) -
     &           beta_up_1_r (i, j, k) / alpha_r (i, j, k)
            rho_h_ww_if = rho_if (i + 1, j, k, 0) *
     &           h_if (i + 1, j, k, 0) * w_if (0) ** 2
            
            d_cap_if (0) = rho_if (i + 1, j, k, 0) * w_if (0)
            s_1_if (0) = rho_h_ww_if * v_1_if (i + 1, j, k, 0)
            s_2_if (0) = rho_h_ww_if * v_2_if (i + 1, j, k, 0)
            s_3_if (0) = rho_h_ww_if * v_3_if (i + 1, j, k, 0)
            tau_if (0) = rho_h_ww_if - d_cap_if (0) -
     &           p_if (i + 1, j, k, 0)

            d_cap_xnu_if (1 : config%qn, 0) = rho_if (i + 1, j, k, 0) *
     &           w_if (0) * xnu_if (i + 1, j, k, 0, 1 : config%qn)
            
ccc   ==================================================================
               
ccc   calculate auxiliary terms: tmp_1, tmp_2, tmp_3
            tmp_1 = v_up_1_if (0) *
     &           (1.0_rk - c_sound_squared_if (i + 1, j, k, 0))
ccc   tmp_1 = v^1 (1 - v_s^2)
            tmp_2 = (1.0_rk - v_squared_if (i + 1, j, k, 0)) *
     &           (g_up_11_r (i, j, k) *
     &           (1.0_rk - v_squared_if (i + 1, j, k, 0) *
     &           c_sound_squared_if (i + 1, j, k, 0)) - 
     &           v_up_1_if (0) ** 2 *
     &           (1.0_rk - c_sound_squared_if (i + 1, j, k, 0)))
ccc   tmp_2 = (1 - v^2) (g^11 (1 - c_s^2 v^2) - v^1 v^1 (1 - c_s^2))
            tmp_3 = alpha_r (i, j, k) / (1.0_rk - 
     &           v_squared_if (i + 1, j, k, 0) *
     &           c_sound_squared_if (i + 1, j, k, 0))
ccc   tmp_3 = alpha / (1 - c_s^2 v^2)
               
ccc   calculate eigenvalues --------------------------------------------
            lambda_minus_if (0) = tmp_3 * (tmp_1 - c_sound_if (0) * 
     &           sqrt (tmp_2)) - beta_up_1_r (i, j, k)
ccc   lambda_-
            lambda_plus_if (0) = tmp_3 * (tmp_1 + c_sound_if (0) * 
     &           sqrt (tmp_2)) - beta_up_1_r (i, j, k)
ccc   lambda_+



ccc   ==================================================================
ccc   right interface

            c_sound_if (1) = sqrt (c_sound_squared_if (i, j, k, 1))
ccc   calculate sound speed c_s

            w_if (1) = 1.0_rk /
     &           sqrt (1.0_rk - v_squared_if (i, j, k, 1))
ccc   calculate Lorentz factor W
               
            v_up_1_if (1) = g_up_11_r (i, j, k) * v_1_if (i, j, k, 1)
ccc   calculate interface value of v^i = g^k v_i
               
            v_up_1_hat_if (1) = v_up_1_if (1) -
     &           beta_up_1_r (i, j, k) / alpha_r (i, j, k)
            rho_h_ww_if = rho_if (i, j, k, 1) *
     &           h_if (i, j, k, 1) * w_if (1) ** 2
            
            d_cap_if (1) = rho_if (i, j, k, 1) * w_if (1)
            s_1_if (1) = rho_h_ww_if * v_1_if (i, j, k, 1)
            s_2_if (1) = rho_h_ww_if * v_2_if (i, j, k, 1)
            s_3_if (1) = rho_h_ww_if * v_3_if (i, j, k, 1)
            tau_if (1) =
     &           rho_h_ww_if - d_cap_if (1) - p_if (i, j, k, 1)

            d_cap_xnu_if (1 : config%qn, 1) = rho_if (i, j, k, 1) *
     &           w_if (1) * xnu_if (i, j, k, 1, 1 : config%qn)
               
ccc   ==================================================================
               
ccc   calculate auxiliary terms: tmp_1, tmp_2, tmp_3
            tmp_1 = v_up_1_if (1) *
     &           (1.0_rk - c_sound_squared_if (i, j, k, 1))
ccc   tmp_1 = v^1 (1 - v_s^2)
            tmp_2 = (1.0_rk - v_squared_if (i, j, k, 1)) *
     &           (g_up_11_r (i, j, k) *
     &           (1.0_rk - v_squared_if (i, j, k, 1) *
     &           c_sound_squared_if (i, j, k, 1)) - 
     &           v_up_1_if (1) ** 2 *
     &           (1.0_rk - c_sound_squared_if (i, j, k, 1)))
ccc   tmp_2 = (1 - v^2) (g^11 (1 - c_s^2 v^2) - v^1 v^1 (1 - c_s^2))
            tmp_3 = alpha_r (i, j, k) / (1.0_rk - 
     &           v_squared_if (i, j, k, 1) *
     &           c_sound_squared_if (i, j, k, 1))
ccc   tmp_3 = alpha / (1 - c_s^2 v^2)
            
ccc   calculate eigenvalues --------------------------------------------
            lambda_minus_if (1) = tmp_3 * (tmp_1 - c_sound_if (1) *
     &           sqrt (tmp_2)) - beta_up_1_r (i, j, k)
ccc   lambda_-
            lambda_plus_if (1) = tmp_3 * (tmp_1 + c_sound_if (1) * 
     &              sqrt (tmp_2)) - beta_up_1_r (i, j, k)
ccc   lambda_+


ccc   ==================================================================
ccc   calculate maximum of left and right interface eigenvalues
c     and their difference
            lambda_minus_min_if = min (lambda_minus_if (0),
     &           lambda_minus_if (1), 0.0_rk)
            lambda_plus_max_if = max (lambda_plus_if (0),
     &           lambda_plus_if (1), 0.0_rk)
            
            lambda_diff_if = lambda_plus_max_if - lambda_minus_min_if


ccc   ==================================================================
ccc   calculate flux at left and right interface            
            do l = 0, 1  
               f_num_r_if (1, l) = d_cap_if (l) * v_up_1_hat_if (l)
               f_num_r_if (2, l) = s_1_if (l) * v_up_1_hat_if (l) + 
     &              p_if (i + 1 - l, j, k, l)
               f_num_r_if (3, l) = s_2_if (l) * v_up_1_hat_if (l)
               f_num_r_if (4, l) = s_3_if (l) * v_up_1_hat_if (l)
               f_num_r_if (5, l) = tau_if (l) * v_up_1_hat_if (l) + 
     &              p_if (i + 1 - l, j, k, l) * v_up_1_if (l)

               f_num_xnu_r_if (1 : config%qn, l) =
     &              d_cap_xnu_if (1 : config%qn, l) *
     &              v_up_1_hat_if (l)
               
            enddo



ccc   ==================================================================
ccc   calculate numerical HLLE fluxes at interface
            f_num_if (i, j, k, 1) =
     &           (lambda_plus_max_if * f_num_r_if (1, 1) -
     &           lambda_minus_min_if * f_num_r_if (1, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (d_cap_if (0) - d_cap_if (1))) / lambda_diff_if
               
            f_num_if (i, j, k, 2) =
     &           (lambda_plus_max_if * f_num_r_if (2, 1) -
     &           lambda_minus_min_if * f_num_r_if (2, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_1_if (0) - s_1_if (1))) / lambda_diff_if

            f_num_if (i, j, k, 3) =
     &           (lambda_plus_max_if * f_num_r_if (3, 1) -
     &           lambda_minus_min_if * f_num_r_if (3, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_2_if (0) - s_2_if (1))) / lambda_diff_if
               
            f_num_if (i, j, k, 4) =
     &           (lambda_plus_max_if * f_num_r_if (4, 1) -
     &           lambda_minus_min_if * f_num_r_if (4, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_3_if (0) - s_3_if (1))) / lambda_diff_if
               
            f_num_if (i, j, k, 5) =
     &           (lambda_plus_max_if * f_num_r_if (5, 1) -
     &           lambda_minus_min_if * f_num_r_if (5, 0) + 
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (tau_if (0) - tau_if (1))) / lambda_diff_if

c     calculate HLL state vector (only 2 components needed)
            e_hll = (lambda_plus_max_if *
     &           (tau_if (0) + d_cap_if (0))-
     &           lambda_minus_min_if *
     &           (tau_if (1) + d_cap_if (1))+
     &           (f_num_r_if (5, 1) + f_num_r_if (1, 1)) -
     &           (f_num_r_if (5, 0) + f_num_r_if (1, 0))) /
     &           lambda_diff_if
            s_1_hll = (lambda_plus_max_if * s_1_if (0) -
     &           lambda_minus_min_if * s_1_if (1) +
     &           f_num_r_if (2, 1) - f_num_r_if (2, 0)) /
     &           lambda_diff_if
c     calculate physical root of (18)
            scr1=e_hll+f_num_if(i,j,k,2)
            lambda_star=2.0_rk*s_1_hll/
     &           (scr1+sqrt(scr1**2-
     &           4.0_rk*
     &           (f_num_if(i,j,k,5)+
     &           f_num_if(i,j,k,1))*
     &           s_1_hll))
            lambda_star_up=lambda_star*g_up_11_r(i,j,k)-
     &           beta_up_1_r(i,j,k)/alpha_r(i,j,k)
c     calculate state variables for intermediate state
c     (depending on the sign of lambda_star_up):
            if (lambda_star_up .lt. 0.0_rk) then
               p_star=lambda_star*
     &              (lambda_plus_max_if *
     &              (tau_if(0)+d_cap_if(0)+
     &              p_if(i+1,j,k,0))-
     &              s_1_if(0)*g_up_11_r(i,j,k))-
     &              s_1_if(0)*(lambda_plus_max_if-
     &              v_up_1_hat_if(0))+p_if(i+1,j,k,0)
               d_cap_star=d_cap_if(0)*
     &              (lambda_plus_max_if-v_up_1_hat_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_1_star=(s_1_if(0)*
     &              (lambda_plus_max_if-v_up_1_hat_if(0))+
     &              p_star-p_if(i+1,j,k,0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_2_star=s_2_if(0)*
     &              (lambda_plus_max_if-v_up_1_hat_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_3_star=s_3_if(0)*
     &              (lambda_plus_max_if-v_up_1_hat_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               tau_star=(tau_if(0)*
     &              (lambda_plus_max_if-v_up_1_hat_if(0))+
     &              p_star*lambda_star_up-
     &              p_if(i+1,j,k,0)*v_up_1_hat_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
            else
               p_star=lambda_star*
     &              (lambda_minus_min_if *
     &              (tau_if(1)+d_cap_if(1)+
     &              p_if(i,j,k,1))-
     &              s_1_if(1)*g_up_11_r(i,j,k))-
     &              s_1_if(1)*(lambda_minus_min_if-
     &              v_up_1_hat_if(1))+p_if(i,j,k,1)
               d_cap_star=d_cap_if(1)*
     &              (lambda_minus_min_if-v_up_1_hat_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_1_star=(s_1_if(1)*
     &              (lambda_minus_min_if-v_up_1_hat_if(1))+
     &              p_star-p_if(i,j,k,1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_2_star=s_2_if(1)*
     &              (lambda_minus_min_if-v_up_1_hat_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_3_star=s_3_if(1)*
     &              (lambda_minus_min_if-v_up_1_hat_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               tau_star=(tau_if(1)*
     &              (lambda_minus_min_if-v_up_1_hat_if(1))+
     &              p_star*lambda_star_up-
     &              p_if(i,j,k,1)*v_up_1_hat_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
            end if
            if (lambda_plus_max_if .le. 0.0_rk) then
               p_star=p_if(i+1,j,k,0)
            end if
            if (lambda_minus_min_if .ge. 0.0_rk) then
               p_star=p_if(i  ,j,k,1)
            endif

            pav_r (i, j, k) = p_star

c     Use HLLC flux for subsonic flow
            lscr1=(lambda_plus_max_if*
     &           lambda_minus_min_if .lt. 0.0_rk).and.
     &           .not.(lshock(i,j,k).or.lshock(i+1,j,k))
            if (lscr1) then
               f_num_if (i, j, k, 1) =
     &              d_cap_star*lambda_star_up
#if MOMENTUM_EQUATION_MODE==0
               f_num_if (i, j, k, 2) =
     &              s_1_star*lambda_star_up+p_star
#else
               f_num_if (i, j, k, 2) =
     &              s_1_star*lambda_star_up
#endif /* MOMENTUM_EQUATION_MODE==0 */
               f_num_if (i, j, k, 3) =
     &              s_2_star*lambda_star_up
               f_num_if (i, j, k, 4) =
     &              s_3_star*lambda_star_up
               f_num_if (i, j, k, 5) =
     &              tau_star*lambda_star_up+
     &              p_star*(lambda_star_up+
     &              beta_up_1_r(i,j,k)/
     &              alpha_r(i,j,k))
#if MOMENTUM_EQUATION_MODE>0
            else
c     Recalculate HLLE momentum flux (without pressure contribution)
               f_num_r_if (2, 0) = s_1_if (0) *
     &              v_up_1_hat_if (0) 
               f_num_r_if (2, 1) = s_1_if (1) *
     &              v_up_1_hat_if (1) 
               f_num_if (i, j, k, 2) =
     &              (lambda_plus_max_if * f_num_r_if (2, 1) -
     &              lambda_minus_min_if * f_num_r_if (2, 0) +
     &              lambda_plus_max_if * lambda_minus_min_if *
     &              (s_1_if (0) - s_1_if (1))) / lambda_diff_if
c     Recalculate effective pressure at interface
               pav_r (i, j, k) =
     &              (lambda_plus_max_if * p_if (i    , j, k, 1) -
     &              lambda_minus_min_if * p_if (i + 1, j, k, 0)) / 
     &              lambda_diff_if
#endif /* MOMENTUM_EQUATION_MODE>0 */
            end if
               
            if (f_num_if (i, j, k, 1) .gt. 0.0_rk) then
               f_num_xnu_if (i, j, k, 1 : config%qn) = 
     &              f_num_if (i, j, k, 1) *
     &              xnu_if (i, j, k, 1, 1 : config%qn)
            else
               f_num_xnu_if (i, j, k, 1 : config%qn) = 
     &              f_num_if (i, j, k, 1) *
     &              xnu_if (i + 1, j, k, 0, 1 : config%qn)
            end if

#ifdef CFC_MHD
            if (lscr1) then 
               v_1_upw = lambda_star_up
c               v_1_upw = s_1_star / (tau_star + d_cap_star + p_star)
               v_2_upw = s_2_star / (tau_star + d_cap_star + p_star)
               v_3_upw = s_3_star / (tau_star + d_cap_star + p_star)
            else if (lambda_plus_max_if .le. 0.0_rk) then
               v_1_upw = v_up_1_hat_if (0)
               v_2_upw = v_2_if (i + 1, j, k, 0)
               v_3_upw = v_3_if (i + 1, j, k, 0)
            else
               v_1_upw = v_up_1_hat_if (1)
               v_2_upw = v_2_if (i, j, k, 1)
               v_3_upw = v_3_if (i, j, k, 1)
            end if
c     Compute v* and transform to orthonormal components
            v_1_upw = alpha_r (i,j,k) * v_1_upw
            v_2_upw = (alpha_r (i,j,k) * v_2_upw * g_up_22_r (i,j,k) -
     &           beta_up_2_r (i,j,k)) *
     &           r_if (i)
            v_3_upw = (alpha_r (i,j,k) * v_3_upw * g_up_33_r (i,j,k) -
     &           beta_up_3_r (i,j,k)) *
     &           r_if (i) * sinus_theta (j)

            if (v_1_upw .gt. 0.0_rk) then
               emf_face (i, j, k, 1, 1) =
     &              v_3_upw * b_1_if (i, j, k, 1) -
     &              v_1_upw * b_3_if (i, j, k, 1) 
               emf_face (i, j, k, 2, 1) =
     &              v_1_upw * b_2_if (i, j, k, 1) -
     &              v_2_upw * b_1_if (i, j, k, 1)
            else
               emf_face (i, j, k, 1, 1) =
     &              v_3_upw * b_1_if (i+1, j, k, 0) -
     &              v_1_upw * b_3_if (i+1, j, k, 0)
               emf_face (i, j, k, 2, 1) =
     &              v_1_upw * b_2_if (i+1, j, k, 0) -
     &              v_2_upw * b_1_if (i+1, j, k, 0)
            end if
            emf_face (i, j, k, 1, 1) = emf_face (i, j, k , 1, 1) *
     &           phi_r (i, j, k) ** 6
            emf_face (i, j, k, 2, 1) = emf_face (i, j, k , 2, 1) *
     &           phi_r (i, j, k) ** 6

#endif /* CFC_MHD */



#ifdef CONVECTION

c           if (config%nsdim .eq. 1) then
c     ------------------------------------------------------------------
c     Convective energy and partial mass fluxes from mixing length
c     approach:

               rhotot=0.5_rk*
     &              (rho(i,j,k)*(1.0_rk+eps(i,j,k))+
     &              rho(i+1,j,k)*(1.0_rk+eps(i+1,j,k)))
c               deltap=max(abs((p(i+1,j,k)-p(i,j,k))),
c     &              1.0e-8_rk*p(i,j,k))
               deltap = ABS (p(i+1,j,k)-p(i,j,k))

c     Local value of Ledoux criterion:
               cledoux=max(
     &              ((rho(i+1,j,k)*(1.0_rk+eps(i+1,j,k))-
     &              rho(i,j,k)*(1.0_rk+eps(i,j,k)))-
     &              2.0_rk*(p(i+1,j,k)-p(i,j,k))/
     &              (c_sound_squared(i+1,j,k)+c_sound_squared(i,j,k)))/
     &              delta_r_if(i),
     &              0.0_rk)*filter(i)

c     Mixing length from pressure scale height:
c               lmix=0.5_rk*(p(i+1,j,k)+p(i,j,k))/deltap*delta_r_if(i)
               IF (deltap .NE. 0.0_rk) THEN
                  lmix=0.5_rk*(p(i+1,j,k)+p(i,j,k))/
     &                 deltap*delta_r_if(i)
               ELSE
                  lmix = 0.0_rk
               ENDIF

c     Local gravitational acceleration:
               glocal=abs(alpha(i+1,j,k)-alpha(i,j,k))/
     &              delta_r_if(i)

c     Typical velocity of convective fluxes:
               vmix=lmix*sqrt(0.5_rk*glocal*cledoux/
     &              rhotot)
               vmix=min(vmix*CONV_PARAM_1,0.5_rk*
     &              (c_sound_if(0)+c_sound_if(1))*
     &              min(1.0_rk,delta_r_if(i)/lmix))
     
               f_num_xnu_if(i,j,k,:)=f_num_xnu_if(i,j,k,:)-
     &              0.25_rk*(rho(i,j,k)+rho(i+1,j,k))*
     &              vmix*(xnnu(i+1,j,k,:)-xnnu(i,j,k,:))*
     &              lmix/delta_r_if(i)
c$$$               if (j.eq.1) print *,'mlt',i,f_num_if(i,j,k,5),lmix,
c$$$     &              cledoux,
c$$$     &              filter(i),
c$$$     &              0.25_rk*vmix*
c$$$     &              (rho(i,j,k)+rho(i+1,j,k))*
c$$$     &              ((eps(i+1,j,k)-eps(i,j,k))+
c$$$     &              0.5_rk*(p(i+1,j,k)+p(i,j,k))*
c$$$     &              (1.0_rk/rho(i+1,j,k)-1.0_rk/rho(i,j,k)))*
c$$$     &              lmix/delta_r_if(i)
               f_num_if(i,j,k,5)=f_num_if(i,j,k,5)-
     &              0.25_rk*vmix*
     &              (rho(i,j,k)+rho(i+1,j,k))*
     &              ((eps(i+1,j,k)-eps(i,j,k))+
     &              0.5_rk*(p(i+1,j,k)+p(i,j,k))*
     &              (1.0_rk/rho(i+1,j,k)-1.0_rk/rho(i,j,k)))*
     &              lmix/delta_r_if(i)
c     ------------------------------------------------------------------
c            endif
#endif /* CONVECTION */

          enddo
      enddo

      return

      END SUBROUTINE calculate_flux_r

c     ==================================================================

#else /* CFC_TRANSPORT2 */


#ifndef MHDFLX
c     ==================================================================

      SUBROUTINE calculate_flux_r

c     ==================================================================

c     Computes HLLC/HLLE fluxes in r-direction (Newtonian version)

c     ------------------------------------------------------------------

      
ccc   This subroutine computes the numerical fluxes
ccc   in the radial direction
      
      USE precision

      USE size_cfc
      USE parameters_cfc
      USE interface_hydro_cfc
      USE grid_cfc
      USE fluxes_cfc
      USE metric_cfc
      USE perm_aux_cfc
      USE shock_cfc

#ifdef CONVECTION
      USE hydro_primitives_cfc, ONLY: rho,eps,p,gamm,
     &     c_sound_squared,xnnu,v_1
      USE hydro_areas_mod
#endif
      USE hydro_primitives_cfc

      use configure

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk, l, ll

      real (kind=rk) :: f_num_r_if (1 : 5, 0 : 1)
      real (kind=rk) :: f_num_xnu_r_if (1 : config%qn, 0 : 1)

      real (kind=rk) :: d_cap_if (0 : 1), s_1_if (0 : 1),
     &     s_2_if (0 : 1), s_3_if (0 : 1), tau_if (0 : 1)
      real (kind=rk) :: d_cap_xnu_if (1 : config%qn, 0 : 1)
      real (kind=rk) :: v_up_1_if (0 : 1)

      real (kind=rk) :: lambda_minus_if (0 : 1)
      real (kind=rk) :: lambda_plus_if (0 : 1)

      real (kind=rk) :: lambda_minus_min_if, lambda_plus_max_if,
     &     lambda_diff_if
      real (kind=rk) :: c_sound_if (0 : 1)

      real (kind=rk) :: e_hll,s_1_hll,lambda_star,lambda_star_up,
     &     p_star,d_cap_star,s_1_star,s_2_star,s_3_star,
     &     tau_star,scr1
#ifdef CFC_MHD
      real (kind=rk) :: v_1_upw, v_2_upw, v_3_upw
#endif /* CFC_MHD */
      logical :: lscr1

#ifdef CONVECTION
      real (kind=rk) :: cledoux,lmix,vmix,glocal,rhotot,deltap
      integer :: i_sonic,filter(0:m),iarea
      integer (kind=ik) :: ixi,ixf,iox,iyi,iyf,ioy,izi,izf,ioz,isd
#endif


#ifdef CONVECTION
      i_sonic=m-1
#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,filter)
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         do i=2,m
            if (v_1(i,j,k)**2 .gt. c_sound_squared(i,j,k))
     &           i_sonic=min(i_sonic,i)
c            if (rho(i,j,k) .lt. 5e12_rk * pc_ggeo)
c     &           i_sonic=min(i_sonic,i)
         end do
      end do
      filter(0)=0
      filter(1:i_sonic-1)=1
      filter(i_sonic:m)=0
c     endif
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
            
            if (ioy .ne.  config%qy) then
               filter (max(ixi-1,1):min(ixf+1,config%qx)) = 0
            endif
         enddo
      end if
#endif /* CONVECTION */


#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,i_n)
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 0, m

ccc   ==================================================================
ccc   left interface
               
            v_up_1_if (0) =
     &           g_up_11_r (i, j, k) * v_1_if (i + 1, j, k, 0)
ccc   calculate interface value of v^i = g^k v_i
               
            d_cap_if (0) = rho_if (i + 1, j, k, 0)
            s_1_if (0) = rho_if (i + 1, j, k, 0) *
     &           v_1_if (i + 1, j, k, 0)
            s_2_if (0) = rho_if (i + 1, j, k, 0) *
     &           v_2_if (i + 1, j, k, 0)
            s_3_if (0) = rho_if (i + 1, j, k, 0) *
     &           v_3_if (i + 1, j, k, 0)
            tau_if (0) = rho_if (i + 1, j, k, 0) *
     &           h_if (i + 1, j, k, 0) - p_if (i + 1, j, k, 0)

            d_cap_xnu_if (1 : config%qn, 0) = rho_if (i + 1, j, k, 0) *
     &           xnu_if (i + 1, j, k, 0, 1 : config%qn)

            c_sound_if (0) =
     &           sqrt (c_sound_squared_if (i + 1, j, k, 0))
ccc   calculate sound speed c_s

ccc   calculate eigenvalues --------------------------------------------
            lambda_minus_if (0) = v_up_1_if (0) - c_sound_if (0) *
     &           sqrt (g_up_11_r (i, j, k))
ccc   lambda_-
            lambda_plus_if (0) = v_up_1_if (0) + c_sound_if (0) *
     &           sqrt (g_up_11_r (i, j, k))
ccc   lambda_+


ccc   ==================================================================
ccc   right interface
            
            v_up_1_if (1) = g_up_11_r (i, j, k) * v_1_if (i, j, k, 1)
ccc   calculate interface value of v^i = g^k v_i
            
            d_cap_if (1) = rho_if (i, j, k, 1)
            s_1_if (1) = rho_if (i, j, k, 1) * v_1_if (i, j, k, 1)
            s_2_if (1) = rho_if (i, j, k, 1) * v_2_if (i, j, k, 1)
            s_3_if (1) = rho_if (i, j, k, 1) * v_3_if (i, j, k, 1)
            tau_if (1) = rho_if (i, j, k, 1) *
     &           h_if (i, j, k, 1) - p_if (i, j, k, 1)

            d_cap_xnu_if (1 : config%qn, 1) = rho_if (i, j, k, 1) *
     &           xnu_if (i, j, k, 1, 1 : config%qn)
            
            c_sound_if (1) = sqrt (c_sound_squared_if (i, j, k, 1))
ccc   calculate sound speed c_s

ccc   calculate eigenvalues --------------------------------------------
            lambda_minus_if (1) = v_up_1_if (1) - c_sound_if (1) *
     &           sqrt (g_up_11_r (i, j, k))
ccc   lambda_-
            lambda_plus_if (1) = v_up_1_if (1) + c_sound_if (1) *
     &           sqrt (g_up_11_r (i, j, k))
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
               f_num_r_if (1, l) = d_cap_if (l) * v_up_1_if (l)
               f_num_r_if (2, l) = s_1_if (l) * v_up_1_if (l) + 
     &              p_if (i + 1 - l, j, k, l)
               f_num_r_if (3, l) = s_2_if (l) * v_up_1_if (l)
               f_num_r_if (4, l) = s_3_if (l) * v_up_1_if (l)
               f_num_r_if (5, l) = (tau_if (l) +
     &              p_if (i + 1 - l, j, k, l)) * v_up_1_if (l)

               f_num_xnu_r_if (1 : config%qn , l) =
     &              d_cap_xnu_if (1 : config%qn, l) *
     &              v_up_1_if (l)
            enddo



ccc   ==================================================================
ccc   calculate numerical fluxes at interface
            f_num_if (i, j, k, 1) =
     &           (lambda_plus_max_if * f_num_r_if (1, 1) -
     &           lambda_minus_min_if * f_num_r_if (1, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (d_cap_if (0) - d_cap_if (1))) / lambda_diff_if
               
            f_num_if (i, j, k, 2) =
     &           (lambda_plus_max_if * f_num_r_if (2, 1) -
     &           lambda_minus_min_if * f_num_r_if (2, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_1_if (0) - s_1_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 3) =
     &           (lambda_plus_max_if * f_num_r_if (3, 1) -
     &           lambda_minus_min_if * f_num_r_if (3, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_2_if (0) - s_2_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 4) =
     &           (lambda_plus_max_if * f_num_r_if (4, 1) -
     &           lambda_minus_min_if * f_num_r_if (4, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_3_if (0) - s_3_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 5) =
     &           (lambda_plus_max_if * f_num_r_if (5, 1) -
     &           lambda_minus_min_if * f_num_r_if (5, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (tau_if (0) - tau_if (1))) / lambda_diff_if
               
c     calculate HLL state vector (only 2 components needed)
            e_hll = (lambda_plus_max_if * d_cap_if (0) -
     &           lambda_minus_min_if * d_cap_if (1) +
     &           f_num_r_if (1, 1) - f_num_r_if (1, 0)) /
     &           lambda_diff_if
            s_1_hll = (lambda_plus_max_if * s_1_if (0) -
     &           lambda_minus_min_if * s_1_if (1) +
     &           f_num_r_if (2, 1) - f_num_r_if (2, 0)) /
     &           lambda_diff_if
c     calculate physical root of (18)
            lambda_star=s_1_hll/e_hll
            lambda_star_up=lambda_star*g_up_11_r(i,j,k)
c     calculate state variables for intermediate state
c     (depending on the sign of lambda_star_up):
            if (lambda_star_up .lt. 0.0_rk) then
               p_star=p_if(i+1,j,k,0)+
     &              d_cap_if(0)*
     &              (lambda_plus_max_if-v_up_1_if(0))*
     &              (lambda_star-v_1_if(i+1,j,k,0))
               d_cap_star=d_cap_if(0)*
     &              (lambda_plus_max_if-v_up_1_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_1_star=(s_1_if(0)*
     &              (lambda_plus_max_if-v_up_1_if(0))+
     &              p_star-p_if(i+1,j,k,0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_2_star=s_2_if(0)*
     &              (lambda_plus_max_if-v_up_1_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               s_3_star=s_3_if(0)*
     &              (lambda_plus_max_if-v_up_1_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
               tau_star=(tau_if(0)*
     &              (lambda_plus_max_if-v_up_1_if(0))+
     &              p_star*lambda_star_up-
     &              p_if(i+1,j,k,0)*v_up_1_if(0))/
     &              (lambda_plus_max_if-lambda_star_up)
            else
               p_star=p_if(i,j,k,1)+
     &              d_cap_if(1)*
     &              (lambda_minus_min_if-v_up_1_if(1))*
     &              (lambda_star-v_1_if(i,j,k,1))
               d_cap_star=d_cap_if(1)*
     &              (lambda_minus_min_if-v_up_1_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_1_star=(s_1_if(1)*
     &              (lambda_minus_min_if-v_up_1_if(1))+
     &              p_star-p_if(i,j,k,1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_2_star=s_2_if(1)*
     &              (lambda_minus_min_if-v_up_1_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               s_3_star=s_3_if(1)*
     &              (lambda_minus_min_if-v_up_1_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
               tau_star=(tau_if(1)*
     &              (lambda_minus_min_if-v_up_1_if(1))+
     &              p_star*lambda_star_up-
     &              p_if(i,j,k,1)*v_up_1_if(1))/
     &              (lambda_minus_min_if-lambda_star_up)
            end if
            if (lambda_plus_max_if .le. 0.0_rk) then
               p_star=p_if(i+1,j,k,0)
            end if
            if (lambda_minus_min_if .ge. 0.0_rk) then
               p_star=p_if(i  ,j,k,1)
            endif
            
            pav_r (i, j, k) = p_star

c     Use HLLC flux for subsonic flow
            lscr1=(lambda_plus_max_if*
     &           lambda_minus_min_if .lt. 0.0_rk).and.
     &           .not.(lshock(i,j,k).or.lshock(i+1,j,k))
#ifdef HLLE_AT_SONIC_POINT
            if (min(lambda_plus_max_if,-lambda_minus_min_if).lt.
     &           0.1_rk*(lambda_plus_max_if-lambda_minus_min_if))
     &           lscr1 = .false.
#endif /* HLLE_AT_SONIC_POINT */
            if (lscr1) then
               f_num_if (i, j, k, 1) =
     &              d_cap_star*lambda_star_up
#if MOMENTUM_EQUATION_MODE==0
               f_num_if (i, j, k, 2) =
     &              s_1_star*lambda_star_up+p_star
#else
               f_num_if (i, j, k, 2) =
     &              s_1_star*lambda_star_up
#endif /* MOMENTUM_EQUATION_MODE==0 */
               f_num_if (i, j, k, 3) =
     &              s_2_star*lambda_star_up
               f_num_if (i, j, k, 4) =
     &              s_3_star*lambda_star_up
               f_num_if (i, j, k, 5) =
     &              tau_star*lambda_star_up+
     &              p_star*lambda_star_up
#if MOMENTUM_EQUATION_MODE>0
            else
c     Recalculate HLLE momentum flux (without pressure contribution)
               f_num_r_if (2, 0) = s_1_if (0) *
     &              v_up_1_if (0) 
               f_num_r_if (2, 1) = s_1_if (1) *
     &              v_up_1_if (1) 
               f_num_if (i, j, k, 2) =
     &              (lambda_plus_max_if * f_num_r_if (2, 1) -
     &              lambda_minus_min_if * f_num_r_if (2, 0) +
     &              lambda_plus_max_if * lambda_minus_min_if *
     &              (s_1_if (0) - s_1_if (1))) / lambda_diff_if
c     Recalculate effective pressure at interface
               pav_r (i, j, k) =
     &              (lambda_plus_max_if * p_if (i    , j, k, 1) -
     &              lambda_minus_min_if * p_if (i + 1, j, k, 0)) / 
     &              lambda_diff_if
#endif /* MOMENTUM_EQUATION_MODE>0 */
            end if

            if (f_num_if(i,j,k,1) .gt. 0.0_rk) then 
               f_num_xnu_if (i, j, k, 1 : config%qn) = 
     &              f_num_if(i, j, k, 1) *
     &              xnu_if (i, j, k, 1, 1 : config%qn)
            else
               f_num_xnu_if (i, j, k, 1 : config%qn) = 
     &              f_num_if(i, j, k, 1) *
     &              xnu_if (i + 1, j, k, 0, 1 : config%qn)            
            end if

#ifdef CFC_MHD
            if (lscr1) then 
               v_1_upw = lambda_star_up
c               v_1_upw = s_1_star / d_cap_star
               v_2_upw = s_2_star / d_cap_star
               v_3_upw = s_3_star / d_cap_star
            else if (lambda_plus_max_if .le. 0.0_rk) then
               v_1_upw = v_up_1_if (0)
               v_2_upw = v_2_if (i + 1, j, k, 0)
               v_3_upw = v_3_if (i + 1, j, k, 0)
            else
               v_1_upw = v_up_1_if (1)
               v_2_upw = v_2_if (i, j, k, 1)
               v_3_upw = v_3_if (i, j, k, 1)
            end if
c     Compute v* and transform to orthonormal components
            v_1_upw = alpha_r (i,j,k) * v_1_upw
            v_2_upw = (alpha_r (i,j,k) * v_2_upw * g_up_22_r (i,j,k) -
     &           beta_up_2_r (i,j,k)) *
     &           r_if (i)
            v_3_upw = (alpha_r (i,j,k) * v_3_upw * g_up_33_r (i,j,k) -
     &           beta_up_3_r (i,j,k)) *
     &           r_if (i) * sinus_theta (j)

            if (v_1_upw .gt. 0.0_rk) then
               emf_face (i, j, k, 1, 1) =
     &              v_3_upw * b_1_if (i, j, k, 1) -
     &              v_1_upw * b_3_if (i, j, k, 1) 
               emf_face (i, j, k, 2, 1) =
     &              v_1_upw * b_2_if (i, j, k, 1) -
     &              v_2_upw * b_1_if (i, j, k, 1)
            else
               emf_face (i, j, k, 1, 1) =
     &              v_3_upw * b_1_if (i+1, j, k, 0) -
     &              v_1_upw * b_3_if (i+1, j, k, 0)
               emf_face (i, j, k, 2, 1) =
     &              v_1_upw * b_2_if (i+1, j, k, 0) -
     &              v_2_upw * b_1_if (i+1, j, k, 0)
            end if
            emf_face (i, j, k, 1, 1) = emf_face (i, j, k , 1, 1) *
     &           phi_r (i, j, k) ** 6
            emf_face (i, j, k, 2, 1) = emf_face (i, j, k , 2, 1) *
     &           phi_r (i, j, k) ** 6
#endif /* CFC_MHD */


#ifdef CONVECTION
c     ------------------------------------------------------------------
c     Convective energy and partial mass fluxes from mixing length
c     approach:

            rhotot=0.5_rk*(rho(i,j,k)+rho(i+1,j,k))
c            deltap=max(abs((p(i+1,j,k)-p(i,j,k))),
c     &           1.0e-8_rk*p(i,j,k))
            deltap = ABS (p(i+1,j,k)-p(i,j,k))

c     Local value of Ledoux criterion:
            cledoux=max(
     &           ((rho(i+1,j,k)-rho(i,j,k))-
     &           2.0_rk*(p(i+1,j,k)-p(i,j,k))/
     &           (c_sound_squared(i+1,j,k)+c_sound_squared(i,j,k)))/
     &           delta_r_if(i),
     &           0.0_rk)*filter(i)
            
c     Mixing length from pressure scale height:
c            lmix=0.5_rk*(p(i+1,j,k)+p(i,j,k))/deltap*delta_r_if(i)
            IF (deltap .NE. 0.0_rk) THEN
               lmix=0.5_rk*(p(i+1,j,k)+p(i,j,k))/
     &              deltap*delta_r_if(i)
            ELSE
               lmix = 0.0_rk
            ENDIF
            
c     Local gravitational acceleration:
            glocal=abs(phi_potential(i+1,j,k)-phi_potential(i,j,k))/
     &           delta_r_if(i)
            
c     Typical velocity of convective fluxes:
            vmix=lmix*sqrt(0.5_rk*glocal*cledoux/
     &           rhotot)
            vmix=min(vmix*CONV_PARAM_1,0.5_rk*
     &           (c_sound_if(0)+c_sound_if(1))*
     &           min(1.0_rk,delta_r_if(i)/lmix))
            
            f_num_xnu_if(i,j,k,:)=f_num_xnu_if(i,j,k,:)-
     &           0.25_rk*(rho(i,j,k)+rho(i+1,j,k))*
     &           vmix*(xnnu(i+1,j,k,:)-xnnu(i,j,k,:))*
     &           lmix/delta_r_if(i)
            f_num_if(i,j,k,5)=f_num_if(i,j,k,5)-
     &           0.25_rk*vmix*
     &           (rho(i,j,k)+rho(i+1,j,k))*
     &           ((eps(i+1,j,k)-eps(i,j,k))+
     &           0.5_rk*(p(i+1,j,k)+p(i,j,k))*
     &           (1.0_rk/rho(i+1,j,k)-1.0_rk/rho(i,j,k)))*
     &           lmix/delta_r_if(i)
c     ------------------------------------------------------------------
#endif /* CONVECTION */
            
         enddo
      enddo

      return

      END SUBROUTINE calculate_flux_r

c     ==================================================================

#else /* MHDFLX */
      
c     ==================================================================

      SUBROUTINE calculate_flux_r

c     ==================================================================

c     Computes HLLD/HLLE fluxes in r-direction (Newtonian MHD version)

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

#ifdef CONVECTION
      USE hydro_primitives_cfc, ONLY: rho,eps,p,gamm,
     &     c_sound_squared,xnnu,v_1
      USE hydro_areas_mod
#endif
      USE hydro_primitives_cfc

      use configure

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk, l, ll

      real (kind=rk) :: f_num_r_if (1 : 8, 0 : 1)
      real (kind=rk) :: f_num_xnu_r_if (1 : config%qn, 0 : 1)

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

      real (kind=rk) :: e_hll,s_1_hll,lambda_star,
     &     lambda_a_minus, lambda_a_plus,
     &     p_star,d_cap_star,v_1_star,v_2_star,v_3_star,
     &     b_1_star,b_2_star,b_3_star,psi_star,
     &     tau_star,scr1,vcln,vclni,vclnl,vclnr
      real (kind=rk) ::
     &     dcsl, dcsr, psl, psr, esl, esr,
     &     v2sl, v3sl, v2sr, v3sr,
     &     b2sl, b3sl, b2sr, b3sr,
     &     swl, swr, tmp1, tmp2
      real (kind=rk) :: slft, srgt, alft, argt
      real (kind=rk) :: v_1_upw, v_2_upw, v_3_upw
      logical :: lscr1, lscr2

#ifdef CONVECTION
      real (kind=rk) :: cledoux,lmix,vmix,glocal,rhotot,deltap
      integer :: i_sonic,filter(0:m),iarea
      integer (kind=ik) :: ixi,ixf,iox,iyi,iyf,ioy,izi,izf,ioz,isd
#endif


#ifdef CONVECTION
      i_sonic=m-1
#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,filter)
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         do i=2,m
            if (v_1(i,j,k)**2 .gt. c_sound_squared(i,j,k))
     &           i_sonic=min(i_sonic,i)
c            if (rho(i,j,k) .lt. 5e12_rk * pc_ggeo)
c     &           i_sonic=min(i_sonic,i)
         end do
      end do
      filter(0)=0
      filter(1:i_sonic-1)=1
      filter(i_sonic:m)=0
      
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
            
            if (ioy .ne.  config%qy) then
               filter (max(ixi-1,1):min(ixf+1,config%qx)) = 0
            endif
         enddo
      end if
#endif /* CONVECTION */


#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,i_n)
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         psi_if_glob (m    , j, k, 1) = 0.0_rk
         psi_if_glob (m + 1, j, k, 0) = 0.0_rk
 
         do i = 0, m

            lscr2 = .true.
                       
ccc   ==================================================================
ccc   left interface
               
            v_or_1_if (0) =
     &           sqrt (g_up_11_r (i, j, k)) * v_1_if (i + 1, j, k, 0)
            v_or_2_if (0) =
     &           sqrt (g_up_22_r (i, j, k)) * v_2_if (i + 1, j, k, 0)
            v_or_3_if (0) =
     &           sqrt (g_up_33_r (i, j, k)) * v_3_if (i + 1, j, k, 0)
ccc   calculate interface value of v^i = g^k v_i
            
            b_1_if (0) = b_1_if_glob (i + 1, j , k, 0)
            b_2_if (0) = b_2_if_glob (i + 1, j , k, 0)
            b_3_if (0) = b_3_if_glob (i + 1, j , k, 0)
            bsq_if (0) =
     &           b_1_if (0) ** 2 + b_2_if (0) ** 2 + b_3_if (0) ** 2
            psi_if (0) = psi_if_glob (i + 1, j, k, 0)
            
            d_cap_if (0) = rho_if (i + 1, j, k, 0)
            s_1_if (0) = d_cap_if (0) * v_or_1_if (0)
            s_2_if (0) = d_cap_if (0) * v_or_2_if (0)
            s_3_if (0) = d_cap_if (0) * v_or_3_if (0)
            tau_if (0) = d_cap_if (0) *
     &           h_if (i + 1, j, k, 0) - p_if (i + 1, j, k, 0)
     &           + 0.5_rk * bsq_if (0)

            d_cap_xnu_if (1 : config%qn, 0) = d_cap_if (0) *
     &           xnu_if (i + 1, j, k, 0, 1 : config%qn)

            c_sound_if (0) =
     &           sqrt (c_sound_squared_if (i + 1, j, k, 0))! +
c     &           v_or_3_if (0) **2)
ccc   calculate sound speed c_s
            cs2_if (0) = c_sound_squared_if (i + 1, j, k, 0) !+
c     &           v_or_3_if (0) **2

ccc   calculate eigenvalues --------------------------------------------
            tmp1 = cs2_if (0) * d_cap_if (0) + bsq_if (0)
            tmp2 = sqrt ((tmp1 + 
     &           sqrt (tmp1 ** 2 - 4.0_rk * cs2_if (0) * d_cap_if (0) *
     &           0*b_1_if (0) ** 2)) / (2.0_rk * d_cap_if (0)))
            lambda_minus_if (0) = v_or_1_if (0) - tmp2
ccc   lambda_-
            lambda_plus_if  (0) = v_or_1_if (0) + tmp2
ccc   lambda_+

ccc   ==================================================================
ccc   right interface
            
            v_or_1_if (1) =
     &           sqrt (g_up_11_r (i, j, k)) * v_1_if (i, j, k, 1)
            v_or_2_if (1) =
     &           sqrt (g_up_22_r (i, j, k)) * v_2_if (i, j, k, 1)
            v_or_3_if (1) =
     &           sqrt (g_up_33_r (i, j, k)) * v_3_if (i, j, k, 1)
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
            
            c_sound_if (1) = sqrt (c_sound_squared_if (i, j, k, 1))!+
c     &           v_or_3_if (1) **2)
ccc   calculate sound speed c_s
            cs2_if (1) = c_sound_squared_if (i, j, k, 1)! +
c     &           v_or_3_if (1) **2

ccc   calculate eigenvalues --------------------------------------------
            tmp1 = cs2_if (1) * d_cap_if (1) + bsq_if (1)
            tmp2 = sqrt ((tmp1 + 
     &           sqrt (tmp1 ** 2 - 4.0_rk * cs2_if (1) * d_cap_if (1) *
     &           0*b_1_if (1) ** 2)) / (2.0_rk * d_cap_if (1)))
            lambda_minus_if (1) = v_or_1_if (1) - tmp2
ccc   lambda_-
            lambda_plus_if  (1) = v_or_1_if (1) + tmp2
ccc   lambda_+

            
ccc   ==================================================================
ccc   calculate maximum of left and right interface eigenvalues
ccc   and their difference
            lambda_minus_min_if = min (lambda_minus_if (0),
     &           lambda_minus_if (1), 0.0_rk)
            lambda_plus_max_if = max (lambda_plus_if (0),
     &           lambda_plus_if (1), 0.0_rk)
            
            lambda_diff_if = lambda_plus_max_if - lambda_minus_min_if
            lambda_b_if = max (abs (v_or_1_if (0)), abs (v_or_1_if (1)))


ccc   ==================================================================
c     Treat 2 * 2 subsytem first

c advective flux of psi           
            f_num_if (i, j, k, 11) = 0.25_rk *
     &          ((psi_if (1) + psi_if (0)) * 
     &           (v_or_1_if (0) + v_or_1_if (1)) -
     &           (psi_if (1) - psi_if (0)) *
     &            abs (v_or_1_if (0) - v_or_1_if (1)))

            vcln = 0.25_rk * min (max (r_if (i), r (1)) * delta_theta,
     &           delta_r_if (i)) /
     &           delta_t
c            vcln = min (vcln, 0.5_rk *
c     &          (sqrt (cs2_if (0) + bsq_if (0) / d_cap_if (0)) +
c     &           sqrt (cs2_if (1) + bsq_if (1) / d_cap_if (1))))
c            if (r_if (i) .eq. 0.0_rk) vcln = 0.0_rk
c     vcln = 0.0_rk
            vclnr = min(vcln, sqrt(cs2_if(0) + bsq_if(0) / d_cap_if(0)))
            vclnl = min(vcln, sqrt(cs2_if(1) + bsq_if(1) / d_cap_if(1)))
c            vcln = 0.5_rk * (vclnl + vclnl)
            vcln = max (vclnl, vclnr)
            
c            srgt = 0.5_rk * (b_1_if (0) + psi_if (0))
            slft = 0.5_rk * (b_1_if (1) + psi_if (1))
            argt = 0.5_rk * (b_1_if (0) - psi_if (0))
c            alft = 0.5_rk * (b_1_if (1) - psi_if (1))

c            if (vcln + 0.5_rk * (v_or_1_if (0) + v_or_1_if (1))
c     &           .lt. 0.0_rk) then
c               b_1_star = b_1_if (0)
c               psi_star = psi_if (0)
c            else if (- vcln + 0.5_rk * (v_or_1_if (0) + v_or_1_if (1))
c     &              .gt. 0.0_rk) then
c               b_1_star = b_1_if (1)
cc               psi_star = psi_if (1)
c            else
            b_1_star = 0.5_rk *
     &           (b_1_if (0) + b_1_if (1) -
     &           (psi_if (0) - psi_if (1)))
            psi_star = 0.5_rk * 
     &           (psi_if (0) + psi_if (1) -
     &           (b_1_if (0) - b_1_if (1)))
c               psi_star =
c     &              (vclnl * slft - vclnr * argt) / vcln
c            end if
c            b_1_if (0) = b_1_star
c            b_1_if (1) = b_1_star

            
ccc   ==================================================================
ccc   calculate flux at left and right interface            
            do l = 0, 1
               f_num_r_if (1, l) = d_cap_if (l) * v_or_1_if (l)
               f_num_r_if (2, l) = s_1_if (l) * v_or_1_if (l) + 
     &              p_if (i + 1 - l, j, k, l) + 0.5_rk * bsq_if (l) -
     &              b_1_if (l) ** 2
               f_num_r_if (3, l) = s_2_if (l) * v_or_1_if (l) -
     &              b_1_if (l) * b_2_if (l)
               f_num_r_if (4, l) = s_3_if (l) * v_or_1_if (l) -
     &              b_1_if (l) * b_3_if (l)
               f_num_r_if (5, l) = (tau_if (l) +
     &              p_if (i + 1 - l, j, k, l) + 0.5_rk * bsq_if (l)) *
     &              v_or_1_if (l) -  b_1_if (l) *
     &              (v_or_1_if (l) * b_1_if (l) +
     &              v_or_2_if (l) * b_2_if (l) +
     &              v_or_3_if (l) * b_3_if (l))
               f_num_r_if (6, l) =   0.0_rk
               f_num_r_if (7, l) =   e_3_if (i + 1 - l, j, k, l)
               f_num_r_if (8, l) = - e_2_if (i + 1 - l, j, k, l)

               f_num_xnu_r_if (1 : config%qn , l) =
     &              d_cap_xnu_if (1 : config%qn, l) *
     &              v_or_1_if (l)
            enddo



ccc   ==================================================================
ccc   calculate numerical fluxes at interface
            f_num_if (i, j, k, 1) =
     &           (lambda_plus_max_if * f_num_r_if (1, 1) -
     &           lambda_minus_min_if * f_num_r_if (1, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (d_cap_if (0) - d_cap_if (1))) / lambda_diff_if
               
            f_num_if (i, j, k, 2) =
     &           (lambda_plus_max_if * f_num_r_if (2, 1) -
     &           lambda_minus_min_if * f_num_r_if (2, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_1_if (0) - s_1_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 3) =
     &           (lambda_plus_max_if * f_num_r_if (3, 1) -
     &           lambda_minus_min_if * f_num_r_if (3, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_2_if (0) - s_2_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 4) =
     &           (lambda_plus_max_if * f_num_r_if (4, 1) -
     &           lambda_minus_min_if * f_num_r_if (4, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (s_3_if (0) - s_3_if (1))) / lambda_diff_if
            
            f_num_if (i, j, k, 5) =
     &           (lambda_plus_max_if * f_num_r_if (5, 1) -
     &           lambda_minus_min_if * f_num_r_if (5, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (tau_if (0) - tau_if (1))) / lambda_diff_if

            f_num_if (i, j, k, 6) = psi_star * vcln
      
            f_num_if (i, j, k, 7) = !v_1 B_2 - B_2 v_1
     &           (lambda_plus_max_if * f_num_r_if (7, 1) -
     &           lambda_minus_min_if * f_num_r_if (7, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (b_2_if (0) - b_2_if (1))) / lambda_diff_if            

            f_num_if (i, j, k, 8) = !v_1 B_3 - B_3 v_1            
     &           (lambda_plus_max_if * f_num_r_if (8, 1) -
     &           lambda_minus_min_if * f_num_r_if (8, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (b_3_if (0) - b_3_if (1))) / lambda_diff_if

            f_num_if (i, j, k, 9) = b_1_star
            
            f_num_if (i, j, k, 10) =
     &           (lambda_plus_max_if * d_cap_if (1) * 
     &           eps_if (i, j, k, 1) * v_or_1_if (1) -
     &           lambda_minus_min_if * d_cap_if (0) *
     &           eps_if (i + 1, j, k, 0) * v_or_1_if (0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (d_cap_if (0) * eps_if (i + 1, j, k, 0) -
     &           d_cap_if (1) * eps_if (i, j, k, 1))) / lambda_diff_if

#ifdef CONSTRAINED            
            emf_face (i, j, k, 1, 1) = ! (v \times B)_2 = v_3 B_1 - v_1 B_3
     &           - (lambda_plus_max_if * f_num_r_if (7, 1) -
     &           lambda_minus_min_if * f_num_r_if (7, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (b_3_if (0) - b_3_if (1))) / lambda_diff_if

            emf_face (i, j, k, 2, 1) = ! (v \times B)_3 = v_1 B_2 - v_2 B_1
     &           (lambda_plus_max_if * f_num_r_if (6, 1) -
     &           lambda_minus_min_if * f_num_r_if (6, 0) +
     &           lambda_plus_max_if * lambda_minus_min_if *
     &           (b_2_if (0) - b_2_if (1))) / lambda_diff_if
#endif /* CONSTRAINED */
            
            if (f_num_if (i, j, k, 1) .gt. 0.0_rk) then
               v_1_star = f_num_if (i, j, k, 1) / d_cap_if (1)
            else
               v_1_star = f_num_if (i, j, k, 1) / d_cap_if (0)
            endif
            
c     calculate HLL state vector (only 2 components needed)
            e_hll = (lambda_plus_max_if * d_cap_if (0) -
     &           lambda_minus_min_if * d_cap_if (1) +
     &           f_num_r_if (1, 1) - f_num_r_if (1, 0)) /
     &           lambda_diff_if
            s_1_hll = (lambda_plus_max_if * s_1_if (0) -
     &           lambda_minus_min_if * s_1_if (1) +
     &           f_num_r_if (2, 1) - f_num_r_if (2, 0)) /
     &           lambda_diff_if
c     calculate physical root of (18)
            lambda_star=s_1_hll/e_hll

c     calculate state variables for intermediate states

c     right star state
            psr  =p_if(i+1,j,k,0)+0.5_rk*bsq_if(0)+
     &           d_cap_if(0)*
     &           (lambda_plus_max_if-v_or_1_if(0))*
     &           (lambda_star-v_or_1_if(0))
            dcsr = d_cap_if(0)*
     &           (lambda_plus_max_if-v_or_1_if(0))/
     &           (lambda_plus_max_if-lambda_star)
            scr1 = dcsr*
     &           (lambda_plus_max_if - lambda_star) ** 2 -
     &           b_1_if (0) ** 2
c            if (j.eq.64) then
c               print *,'flux',i,j,psr,dcsr,scr1,e_hll,s_1_hll
c               print *,'flux2',i,j,lambda_plus_max_if,
c     &              lambda_minus_min_if,lambda_star
c            end if
            if (scr1 .ge. 1e-5_rk * psr) then
               v2sr = v_or_2_if (0) -
     &              b_1_if (0) * b_2_if (0) *
     &              (lambda_star - v_or_1_if (0)) / scr1
               v3sr = v_or_3_if (0) -
     &              b_1_if (0) * b_3_if (0) *
     &              (lambda_star - v_or_1_if (0)) / scr1
               b2sr = b_2_if (0) * (d_cap_if (0) *
     &              (lambda_plus_max_if - v_or_1_if (0)) ** 2 -
     &              b_1_if (0) ** 2) / scr1
               b3sr = b_3_if (0) * (d_cap_if (0) *
     &              (lambda_plus_max_if - v_or_1_if (0)) ** 2 -
     &              b_1_if (0) ** 2) / scr1
            else
c               print *,'special case, flux_r (r)',i,j,k
               lscr2 = .false.
               v2sr = v_or_2_if (0)
               v3sr = v_or_3_if (0)
               b2sr = b_2_if (0)
               b3sr = b_3_if (0)
            end if
            esr  = (tau_if(0)*
     &           (lambda_plus_max_if-v_or_1_if(0))+
     &           psr*lambda_star-
     &           (p_if(i+1,j,k,0)+0.5_rk*bsq_if(0))*v_or_1_if(0)+
     &           b_1_if (0) * (
     &           v_or_1_if (0) * b_1_if (0) +
     &           v_or_2_if (0) * b_2_if (0) +
     &           v_or_3_if (0) * b_3_if (0) -
     &           lambda_star * b_1_if (0) - v2sr * b2sr - v3sr * b3sr))/
     &           (lambda_plus_max_if-lambda_star)
            lambda_a_plus = lambda_star + abs(b_1_if(0))/sqrt (dcsr)
            
c     left star state
            psl  =p_if(i,j,k,1)+0.5_rk*bsq_if(1)+
     &           d_cap_if(1)*
     &           (lambda_minus_min_if-v_or_1_if(1))*
     &           (lambda_star-v_or_1_if(1))
            dcsl = d_cap_if(1)*
     &           (lambda_minus_min_if-v_or_1_if(1))/
     &           (lambda_minus_min_if-lambda_star)
            scr1 = dcsl*
     &           (lambda_minus_min_if - lambda_star) ** 2 -
     &           b_1_if (1) ** 2

            if (scr1 .ge. 1e-5_rk * psl) then
               v2sl = v_or_2_if (1) -
     &              b_1_if (1) * b_2_if (1) *
     &              (lambda_star - v_or_1_if (1)) / scr1
               v3sl = v_or_3_if (1) -
     &              b_1_if (1) * b_3_if (1) *
     &              (lambda_star - v_or_1_if (1)) / scr1
               b2sl = b_2_if (1) * (d_cap_if (1) *
     &              (lambda_minus_min_if - v_or_1_if (1)) ** 2 -
     &              b_1_if (1) ** 2) / scr1
               b3sl = b_3_if (1) * (d_cap_if (1) *
     &              (lambda_minus_min_if - v_or_1_if (1)) ** 2 -
     &              b_1_if (1) ** 2) / scr1
            else
c               print *,'special case, flux_r (l)',i,j,k
               lscr2 = .false.
               v2sl = v_or_2_if (1)
               v3sl = v_or_3_if (1)
               b2sl = b_2_if (1)
               b3sl = b_3_if (1)               
            end if
            esl  = (tau_if(1)*
     &           (lambda_minus_min_if-v_or_1_if(1))+
     &           psl*lambda_star-
     &           (p_if(i,j,k,1)+0.5_rk*bsq_if(1))*v_or_1_if(1)+
     &           b_1_if (1) * (
     &           v_or_1_if (1) * b_1_if (1) +
     &           v_or_2_if (1) * b_2_if (1) +
     &           v_or_3_if (1) * b_3_if (1) -
     &           lambda_star * b_1_if (1) - v2sl * b2sl - v3sl * b3sl))/
     &           (lambda_minus_min_if-lambda_star)            
            lambda_a_minus = lambda_star - abs(b_1_if(1))/sqrt (dcsl)
            
            swl = 0.5_rk + sign (0.5_rk,  lambda_star)
            if (lambda_star .eq. 0.0_rk) swl = 0.5_rk
            swr = 1.0_rk - swl
            p_star     = psl * swl + psr * swr
            d_cap_star = dcsl * swl + dcsr * swr
            v_2_star   = v2sl * swl + v2sr * swr
            v_3_star   = v3sl * swl + v3sr * swr
c            b_1_star   = 0.5_rk * (b_1_if (0) + b_1_if (1))
            b_2_star   = b2sl * swl + b2sr * swr
            b_3_star   = b3sl * swl + b3sr * swr
            tau_star   = esl  * swl + esr * swr

#if MHD_SOLVER == HLLD
c     If  we're in one of the star states, this is already the state
c     that we need. If we're between the Alfven waves, we need
c     compute the double star state.
            if (lambda_a_minus * lambda_a_plus .le. 0.0_rk .and.
     &           lambda_a_minus .gt. lambda_minus_min_if .and.
     &           lambda_a_plus  .lt. lambda_plus_max_if) then
               dcsl = sqrt (dcsl)
               dcsr = sqrt (dcsr)
               scr1 = dcsl+dcsr
               v_2_star = (dcsl*v2sl+dcsr*v2sr+
     &              sign(b2sr-b2sl,b_1_star))/scr1
               v_3_star = (dcsl*v3sl+dcsr*v3sr+
     &              sign(b3sr-b3sl,b_1_star))/scr1
               b_2_star = (dcsl*b2sr+dcsr*b2sl+
     &              dcsl*dcsr*sign(v2sr-v2sl,b_1_star))/scr1
               b_3_star = (dcsl*b3sr+dcsr*b3sl+
     &              dcsl*dcsr*sign(v3sr-v3sl,b_1_star))/scr1
               tau_star = tau_star +
     &              swr  * dcsr *
     &              sign(v2sr*b2sr+v3sr*b3sr- !v_1 B_1 cancels
     &              v_2_star*b_2_star-v_3_star*b_3_star,b_1_star) -
     &              swl  * dcsl *
     &              sign(v2sl*b2sl+v3sl*b3sl- !v_1 B_1 cancels
     &              v_2_star*b_2_star-v_3_star*b_3_star,b_1_star)
            end if
#endif /* MHD_SOLVER == HLLD */
            
            pav_r (i, j, k) = p_star

            if (min (lambda_a_minus-lambda_minus_min_if,
     &           lambda_plus_max_if-lambda_a_plus)
     &           .lt. 5e-2_rk*lambda_diff_if) lscr2 = .false.
            
c     Use HLLC flux for subsonic flow
            lscr1=(lambda_plus_max_if*
     &           lambda_minus_min_if .lt. 0.0_rk).and.
     &           .not.(lshock(i,j,k).or.lshock(i+1,j,k))
#ifdef HLLE_AT_SONIC_POINT
            if (min(lambda_plus_max_if,-lambda_minus_min_if).lt.
     &           0.1_rk*(lambda_plus_max_if-lambda_minus_min_if))
     &           lscr1 = .false.
#endif /* HLLE_AT_SONIC_POINT */
#if MHD_SOLVER == HLLE
            lscr1 = .false.
#endif            
            if (lscr1 .and. lscr2) then
               f_num_if (i, j, k, 1) =
     &              d_cap_star*lambda_star
#if MOMENTUM_EQUATION_MODE==0
               f_num_if (i, j, k, 2) =
     &              d_cap_star*lambda_star**2+p_star-b_1_star**2
#else
               f_num_if (i, j, k, 2) =
     &              d_cap_star*lambda_star**2-b_1_star**2
#endif /* MOMENTUM_EQUATION_MODE==0 */
               f_num_if (i, j, k, 3) =
     &              d_cap_star*v_2_star*lambda_star-b_2_star*b_1_star
               f_num_if (i, j, k, 4) =
     &              d_cap_star*v_3_star*lambda_star-b_3_star*b_1_star
               f_num_if (i, j, k, 5) =
     &              tau_star*lambda_star+
     &              p_star*lambda_star - b_1_star *
     &              (lambda_star*b_1_star+
     &              v_2_star*b_2_star+v_3_star*b_3_star)

               v_1_upw = lambda_star
               v_2_upw = v_2_star
               v_3_upw = v_3_star

               f_num_if (i, j, k, 7) =
     &              v_1_upw * b_2_star - v_2_upw * b_1_star
               f_num_if (i, j, k, 8) =
     &              v_1_upw * b_3_star - v_3_upw * b_1_star

#ifdef CONSTRAINED               
               emf_face (i, j, k, 1, 1) =
     &              v_3_upw * b_1_star - v_1_upw * b_3_star
               emf_face (i, j, k, 2, 1) =
     &              v_1_upw * b_2_star - v_2_upw * b_1_star

               if (v_1_upw .gt. 0.0_rk) then
                  emf_face (i, j, k, 1, 1) =
     &                 e_2_if (i, j, k, 1) +
     &                 (emf_face (i, j, k, 1, 1) -
     &                 (v_or_3_if (1) * b_1_if (1) -
     &                 v_or_1_if (1) * b_3_if (1)))
                  emf_face (i, j, k, 2, 1) =
     &                 e_3_if (i, j, k, 1) +
     &                 (emf_face (i, j, k, 2, 1) -
     &                 (v_or_1_if (1) * b_2_if (1) -
     &                 v_or_2_if (1) * b_1_if (1)))
               else
                  emf_face (i, j, k, 1, 1) =
     &                 e_2_if (i + 1, j, k, 0) +
     &                 (emf_face (i, j, k, 1, 1) -
     &                 (v_or_3_if (0) * b_1_if (0) -
     &                  v_or_1_if (0) * b_3_if (0))
                  emf_face (i, j, k, 2, 1) =
     &                 e_3_if (i + 1, j, k, 0) +
     &                 (emf_face (i, j, k, 2, 1) -
     &                 (v_or_1_if (0) * b_2_if (0) -
     &                  v_or_2_if (0) * b_1_if (0))
               end if
#endif /* CONSTRAINED */
               
c     diffusive term for stabilising constraint terms
               f_num_if (i, j, k, 7) = f_num_if (i, j, k, 7) +
     &              0.5_rk * lambda_diff_if *
     &              (b_2_if (1) - b_2_if (0))
               f_num_if (i, j, k, 8) = f_num_if (i, j, k, 8) +
     &              0.5_rk * lambda_diff_if *
     &              (b_3_if (1) - b_3_if (0))
               
c               emf_face (i, j, k, 1, 1) = emf_face (i, j, k, 1, 1) -
c     &              2.0_rk *
c     &              max(abs (v_or_3_if (1)), abs(v_or_3_if (0))) *
c     &              (b_3_if (1) - b_3_if (0))

               v_1_star = lambda_star

#if MOMENTUM_EQUATION_MODE>0
               
            else
c     Recalculate HLLE momentum flux (without pressure contribution)
               f_num_r_if (2, 0) = s_1_if (0) * v_or_1_if (0) -
     &              b_1_if (0) ** 2
               f_num_r_if (2, 1) = s_1_if (1) * v_or_1_if (1) -
     &              b_1_if (1) ** 2
               f_num_if (i, j, k, 2) =
     &              (lambda_plus_max_if * f_num_r_if (2, 1) -
     &              lambda_minus_min_if * f_num_r_if (2, 0) +
     &              lambda_plus_max_if * lambda_minus_min_if *
     &              (s_1_if (0) - s_1_if (1))) / lambda_diff_if
c     Recalculate effective pressure at interface
               pav_r (i, j, k) =
     &              (lambda_plus_max_if * (p_if (i    , j, k, 1) +
     &              0.5_rk * bsq_if (1)) -
     &              lambda_minus_min_if * (p_if (i + 1, j, k, 0) +
     &              0.5_rk * bsq_if (0))) / 
     &              lambda_diff_if
#endif /* MOMENTUM_EQUATION_MODE>0 */
            end if

c     add energy flux in divergence cleaning mode
            tmp1 = 0.5_rk * (b_1_if (0) + b_1_if (1))
            tmp2 = 0.5_rk * (psi_if (0) + psi_if (1))
            f_num_if (i, j, k, 5) = f_num_if (i, j, k, 5) +
     &           vcln  * (tmp1 * tmp2 +
     &           0.5_rk * tmp1 * (b_1_if (0) - b_1_if (1)) +
     &           0.5_rk * tmp2 * (psi_if (0) - psi_if (1)))           
c     &           f_num_if (i, j, k, 6) * b_1_star +
c     &           f_num_if (i, j, k, 11) * psi_star
            
c     Transform to contravriant/mixed flux components in
c     spherical polar reference metric:
c            f_num_if (i, j, k, 1) = f_num_if (i, j, k, 1)
c            f_num_if (i, j, k, 2) = f_num_if (i, j, k, 2)
            f_num_if (i, j, k, 3) = f_num_if (i, j, k, 3) * r_if (i)
            f_num_if (i, j, k, 4) = f_num_if (i, j, k, 4) * r_if (i) *
     &           sinus_theta (j)
c            f_num_if (i, j, k, 5) = f_num_if (i, j, k, 5)
c            f_num_if (i, j, k, 6) = f_num_if (i, j, k, 6)
            f_num_if (i, j, k, 7) = f_num_if (i, j, k, 7) * r_if_inv (i)
            f_num_if (i, j, k, 8) = f_num_if (i, j, k, 8) *
     &           r_if_inv (i) * sinus_theta_inv (j)
c            f_num_if (i, j, k, 9) = f_num_if (i, j, k, 9)
c            f_num_if (i, j, k, 10) = f_num_if (i, j, k, 10)
c            f_num_if (i, j, k, 11) = f_num_if (i, j, k, 11)
            
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
     &              eps_if (i + 1, j, k, 0)
               f_num_xnu_if (i, j, k, 1 : config%qn) = 
     &              f_num_if(i, j, k, 1) *
     &              xnu_if (i + 1, j, k, 0, 1 : config%qn)            
            end if
c            f_num_if (i, j, k, 11) = v_1_star
            
#ifdef CONVECTION
c     ------------------------------------------------------------------
c     Convective energy and partial mass fluxes from mixing length
c     approach:

            rhotot=0.5_rk*(rho(i,j,k)+rho(i+1,j,k))
c            deltap=max(abs((p(i+1,j,k)-p(i,j,k))),
c     &           1.0e-8_rk*p(i,j,k))
            deltap = ABS (p(i+1,j,k)-p(i,j,k))

c     Local value of Ledoux criterion:
            cledoux=max(
     &           ((rho(i+1,j,k)-rho(i,j,k))-
     &           2.0_rk*(p(i+1,j,k)-p(i,j,k))/
     &           (c_sound_squared(i+1,j,k)+c_sound_squared(i,j,k)))/
     &           delta_r_if(i),
     &           0.0_rk)*filter(i)
            
c     Mixing length from pressure scale height:
c            lmix=0.5_rk*(p(i+1,j,k)+p(i,j,k))/deltap*delta_r_if(i)
            IF (deltap .NE. 0.0_rk) THEN
               lmix=0.5_rk*(p(i+1,j,k)+p(i,j,k))/
     &              deltap*delta_r_if(i)
            ELSE
               lmix = 1e-15_rk
            ENDIF
            
c     Local gravitational acceleration:
            glocal=abs(phi_potential(i+1,j,k)-phi_potential(i,j,k))/
     &           delta_r_if(i)
            
c     Typical velocity of convective fluxes:
            vmix=lmix*sqrt(0.5_rk*glocal*cledoux/
     &           rhotot)
            vmix=min(vmix*CONV_PARAM_1,0.5_rk*
     &           (c_sound_if(0)+c_sound_if(1))*
     &           min(1.0_rk,delta_r_if(i)/lmix))
            
            f_num_xnu_if(i,j,k,:)=f_num_xnu_if(i,j,k,:)-
     &           0.25_rk*(rho(i,j,k)+rho(i+1,j,k))*
     &           vmix*(xnnu(i+1,j,k,:)-xnnu(i,j,k,:))*
     &           lmix/delta_r_if(i)
            f_num_if(i,j,k,5)=f_num_if(i,j,k,5)-
     &           0.25_rk*vmix*
     &           (rho(i,j,k)+rho(i+1,j,k))*
     &           ((eps(i+1,j,k)-eps(i,j,k))+
     &           0.5_rk*(p(i+1,j,k)+p(i,j,k))*
     &           (1.0_rk/rho(i+1,j,k)-1.0_rk/rho(i,j,k)))*
     &           lmix/delta_r_if(i)
            f_num_if(i,j,k,10)=f_num_if(i,j,k,10)-
     &           0.25_rk*vmix*
     &           (rho(i,j,k)+rho(i+1,j,k))*
     &           ((eps(i+1,j,k)-eps(i,j,k))+
     &           0.5_rk*(p(i+1,j,k)+p(i,j,k))*
     &           (1.0_rk/rho(i+1,j,k)-1.0_rk/rho(i,j,k)))*
     &           lmix/delta_r_if(i)
c     ------------------------------------------------------------------
#endif /* CONVECTION */
            
         enddo
      enddo
c      stop
      return

      END SUBROUTINE calculate_flux_r

c     ==================================================================

#endif /* MHDFLX */

      
#endif /* CFC_TRANSPORT2 */




      END MODULE gr_flux_r

c     ==================================================================
