c     ==================================================================

      MODULE gr_initial_setup

c     ==================================================================

      IMPLICIT NONE

      CONTAINS




c     ==================================================================

      SUBROUTINE initialize_conserved_quantities

c     ==================================================================

      USE precision
      
      USE size_cfc
      USE parameters_cfc      
      USE conserved_cfc
      USE hydro_primitives_cfc
      USE metric_cfc
      USE grid_cfc
      USE perm_aux_cfc

      USE eos3d_routine

      USE gr_boundary
      USE gr_rescale

      use configure

      IMPLICIT NONE
      real(kind=rk)     :: eos3ds_self(2), eos3ds_children(2)
      integer (kind=ik) :: i, j, k, jk, i_nuc

      real (kind=rk)    :: rho_h_ww

      logical           :: ler

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

#ifdef CFC_TRANSPORT
         call eos3ds_CFC(3,j,k,ler,eos3ds_self,eos3ds_children)
#endif

         do i = 1, m

#ifdef CFC_TRANSPORT2
            rho_h_ww = rho (i, j, k) * h (i, j, k) *
     &           w (i, j, k) ** 2
            
            d_cap_hat (i, j, k) =  rho (i, j, k) * w (i, j, k) *
     &           sqrt_gamma (i, j, k)
            s_1_hat   (i, j, k) =  rho_h_ww * v_1 (i, j, k) *
     &           sqrt_gamma (i, j, k)
            s_2_hat   (i, j, k) =  rho_h_ww * v_2 (i, j, k) *
     &           sqrt_gamma (i, j, k)
            s_3_hat   (i, j, k) =  rho_h_ww * v_3 (i, j, k) *
     &           sqrt_gamma (i, j, k)
            tau_hat   (i, j, k) = (rho_h_ww - p (i, j, k)) *
     &           sqrt_gamma (i, j, k) - d_cap_hat (i, j, k)

#else /* CFC_TRANSPORT2*/
            v_squared (i, j, k) = 
     &           v_1 (i, j, k) ** 2 * g_up_11 (i, j, k)+ 
     &           v_2 (i, j, k) ** 2 * g_up_22 (i, j, k)+
     &           v_3 (i, j, k) ** 2 * g_up_33 (i, j, k)
            h         (i, j, k) = eps (i, j, k) + p (i, j, k) /
     &           rho (i, j, k) + 0.5_rk * v_squared (i, j, k)
            
            d_cap_hat (i, j, k) = rho (i, j, k)
            s_1_hat   (i, j, k) = rho (i, j, k) * v_1 (i, j, k)
            s_2_hat   (i, j, k) = rho (i, j, k) * v_2 (i, j, k)
            s_3_hat   (i, j, k) = rho (i, j, k) * v_3 (i, j, k)
            tau_hat   (i, j, k) = rho (i, j, k) * h (i, j, k) -
     &           p (i, j, k)
#if defined(CFC_MHD) &&  !defined(PASSIVE_MHD)
     &           + 0.5_rk * (b_1 (i, j, k) ** 2 + b_2 (i, j, k) ** 2 +
     &           b_3 (i, j, k) ** 2)
#ifndef CONSTRAINED
     &           + 0.5_rk * psi (i, j, k) ** 2
#endif
            
#endif /* CFC_MHD and not PASSIVE_MHD */            
#endif /* CFC_TRANSPORT2*/

#ifdef CFC_MHD
            b_cap_1_hat (i, j, k) = b_1 (i, j, k) * sqrt_gamma (i, j, k)
            b_cap_2_hat (i, j, k) = b_2 (i, j, k) * sqrt_gamma (i, j, k)
            b_cap_3_hat (i, j, k) = b_3 (i, j, k) * sqrt_gamma (i, j, k)
#endif /* CFC_MHD */     

#ifdef NEC_COMPILER               
!CDIR EXPAND=25
#endif
            do i_nuc = 1, config%qn
               d_cap_xnu_hat (i, j, k, i_nuc) = d_cap_hat (i, j, k) *
     &              xnnu (i, j, k, i_nuc)
            enddo

         enddo
      enddo

      call rescale_velocities (1)
      call hydro_boundary_conditions
      call rescale_velocities (2)


      return

      END SUBROUTINE initialize_conserved_quantities

c     ==================================================================




c     ==================================================================

      SUBROUTINE calculate_perm_aux_values

c     ==================================================================

      USE precision

      USE size_cfc
      USE grid_cfc
      USE parameters_cfc
      USE perm_aux_cfc
      USE totare_hy, ONLY: dvytot, dvztot

      USE gr_ppm

      use configure
      IMPLICIT NONE

      real (kind=rk) :: dvy
      integer (kind=ik) :: i, j, k, l

#if MOMENTUM_EQUATION_MODE==2

      if (config%nsdim .gt. 1) then
         cos_dtheta = cos (delta_theta / 2.0_rk)
         sin_dtheta = sin (delta_theta / 2.0_rk)
         cos_dvarphi = cos (delta_varphi / 2.0_rk)
         sin_dvarphi = sin (delta_varphi / 2.0_rk)
      else
         cos_dtheta = 1.0_rk
         sin_dtheta = 0.0_rk
         cos_dvarphi = 1.0_rk
         sin_dvarphi = 0.0_rk
      endif
#endif /* MOMENTUM_EQUATION_MODE==2 */

      delta_theta_inv  = 1.0_rk / delta_theta
      delta_ttheta_inv = 1.0_rk / delta_theta ** 2
      delta_theta_if_inv  = 1.0_rk / delta_theta_if

      delta_varphi_inv = 1.0_rk / delta_varphi
      delta_varphi_if_inv  = 1.0_rk / delta_varphi_if

      do i = - 3, m + 4
         r_inv (i) = 1.0_rk / r (i)
         rr (i) = r (i) ** 2
         rr_inv (i) = 1.0_rk / rr (i)
      enddo
      do i = 1, m + 1
         amp_inv (i) = 1.0_rk / amp (i)
         aamp (i) = amp (i) ** 2
         r_if_inv (i) = 1.0_rk / r_if (i)
         rr_if (i) = r_if (i) ** 2
         delta_r_if_inv (i) = 1.0_rk / delta_r_if (i)
         delta_r_amp_inv (i) =
     &     1.0_rk / (delta_r_if (i - 1) * amp (i) * (amp (i) + 1.0_rk))
         delta_rrr (i) = (r_if (i) ** 3 - r_if (i - 1) ** 3) / 3.0_rk
      enddo
      do i = 1, m
         rrr (i) = r (i) ** 3
         delta_r_inv (i) = 1.0_rk / delta_r (i)
         delta_r_mod_inv (i) = 3.0_rk /
     &        (r_if (i) ** 3 - r_if (i - 1) ** 3)
      enddo
      delta_r_if_inv (- 1) = 1.0_rk / delta_r_if (- 1)
      delta_r_if_inv (0) = 1.0_rk / delta_r_if (0)
      r_if_inv (0) = 0.0_rk
      rr_if (0) = r_if (0) ** 2
      r_inv (- 3) = 1.0_rk / r (- 3)
      r_inv (- 2) = 1.0_rk / r (- 2)
      r_inv (- 1) = 1.0_rk / r (- 1)
      r_inv (0) = 1.0_rk / r (0)
      r_inv (m + 2) = 1.0_rk / r (m + 2)
      r_inv (m + 3) = 1.0_rk / r (m + 3)
      r_inv (m + 4) = 1.0_rk / r (m + 4)

      do j = - 3, n + 4
         ssinus_theta_inv (j) = 1.0_rk / sin (theta (j)) ** 2
      enddo
      do j = 1, n
         cotan_theta (j) = 1.0_rk / tan (theta (j))
         ccotan_theta (j) = cotan_theta (j) ** 2
         sinus_theta (j) = sin (theta (j))
         cosin_theta (j) = cos (theta (j))
         ssinus_theta (j) = sinus_theta (j) ** 2
         ccosin_theta (j) = cosin_theta (j) ** 2
         sinus_theta_inv (j) = 1.0_rk / sinus_theta (j)
         sinus_theta_if (j) = sin (theta_if (j))
         ssinus_theta_if (j) = sinus_theta_if (j) ** 2
         cosin_theta_if (j) = cos (theta_if (j))
         sinuscosin_theta (j) = sinus_theta (j) * cosin_theta (j)
         delta_cosin_theta (j) = cos (theta_if (j - 1)) -
     &        cos (theta_if (j))
         
         if (config%nsdim .gt. 1) then
            delta_theta_mod_inv (j) = 1.0_rk /
     &           (- cos (theta_if (j)) + cos (theta_if (j - 1)))
         else
            delta_theta_mod_inv (j) = 2.0_rk
         endif
      enddo
      sinus_theta_if (- 1) = sin (theta_if (- 1))
      sinus_theta_if (0) = sin (theta_if (0))
      cosin_theta_if (0) = cos (theta_if (0))
      ssinus_theta_if (0) = sinus_theta_if (0) ** 2
      sinus_theta (- 3) = sin (theta (- 3))
      sinus_theta (- 2) = sin (theta (- 2))
      sinus_theta (- 1) = sin (theta (- 1))
      sinus_theta (0) = sin (theta (0))
      sinus_theta_inv (- 3) = 1.0_rk / sinus_theta (- 3)
      sinus_theta_inv (- 2) = 1.0_rk / sinus_theta (- 2)
      sinus_theta_inv (- 1) = 1.0_rk / sinus_theta (- 1)
      sinus_theta_inv (0) = 1.0_rk / sinus_theta (0)
      sinus_theta_if (n + 1) = sin (theta_if (n + 1))
      sinus_theta (n + 1) = sin (theta (n + 1))
      sinus_theta (n + 2) = sin (theta (n + 2))
      sinus_theta (n + 3) = sin (theta (n + 3))
      sinus_theta (n + 4) = sin (theta (n + 4))
      sinus_theta_inv (n + 1) = 1.0_rk / sinus_theta (n + 1)
      sinus_theta_inv (n + 2) = 1.0_rk / sinus_theta (n + 2)
      sinus_theta_inv (n + 3) = 1.0_rk / sinus_theta (n + 3)
      sinus_theta_inv (n + 4) = 1.0_rk / sinus_theta (n + 4)
      cosin_theta (- 3) = cos (theta (- 3))
      cosin_theta (- 2) = cos (theta (- 2))
      cosin_theta (- 1) = cos (theta (- 1))
      cosin_theta (0) = cos (theta (0))
      cosin_theta (n + 1) = cos (theta (n + 1))
      cosin_theta (n + 2) = cos (theta (n + 2))
      cosin_theta (n + 3) = cos (theta (n + 3))
      cosin_theta (n + 4) = cos (theta (n + 4))
      cotan_theta (0) = cotan_theta (1)
      cotan_theta (n + 1) = cotan_theta (n)

      do j = 1, n - 1
         sinus_theta_inv_if (j) = 1.0_rk / sinus_theta_if (j)
      end do
      sinus_theta_inv_if (0) = 0.0_rk 
      sinus_theta_inv_if (n) = 0.0_rk

#ifdef CFC_MHD
      sinus_theta (n/2+1:n+4)= sinus_theta (n/2:-3:-1)
      cosin_theta (n/2+1:n+4)=-cosin_theta (n/2:-3:-1)
      sinus_theta_inv (n/2+1:n+4)= sinus_theta_inv (n/2:-3:-1)
      sinus_theta_if (n/2+1:n+1)= sinus_theta_if (n/2-1:-1:-1)
      sinus_theta_if (n/2) = 1.0d0
      cosin_theta_if (n/2) = 0.0d0
#endif /* CFC_MHD */

      do k = 1, o
         sinus_varphi (k) = sin (varphi (k))
         ssinus_varphi (k) = sinus_varphi (k) ** 2
         cosin_varphi (k) = cos (varphi (k))
         ccosin_varphi (k) = cosin_varphi (k) ** 2
         sinuscosin_varphi (k) = sin (varphi (k)) * cos (varphi (k))
      enddo
      sinus_varphi (0) = sinus_varphi (o)
      sinus_varphi (o + 1) = sinus_varphi (1)
      cosin_varphi (0) = cosin_varphi (o)
      cosin_varphi (o + 1) = cosin_varphi (1)

#ifdef CFC_MHD
      sumdvy0(:) = 0.0_rk
      do k = 1, o
         do j = 1, n

            dvy=dvytot(j)*dvztot(k)
            sumdvy0(0) = sumdvy0(0) + dvy
            sumdvy0(1) = sumdvy0(1) + (ssinus_varphi (k) +
     &           ccosin_varphi (k) * ccosin_theta (j)) * dvy
            sumdvy0(2) = sumdvy0(2) + (ccosin_varphi (k) +
     &           ssinus_varphi (k) * ccosin_theta (j)) * dvy
            sumdvy0(3) = sumdvy0(3) + ssinus_theta (j) * dvy

         end do
      end do
#endif /* CFC_MHD */      

c     calculate PPM coefficients and set parameters
         
      call coefficients_r
      call coefficients_theta
      call coefficients_varphi
      
      eta_1 = 5.0_rk
      eta_2 = 0.05_rk
      epsilon_1 = 0.1_rk
      k_cap_0 = 1.0_rk
      epsilon_2_r = 0.33_rk
      epsilon_2_theta = 0.33_rk
      epsilon_2_varphi = 0.33_rk
      omega_1 = 0.75_rk
      omega_2 = 10.0_rk
         

      gw_factor_1 = 8.0_rk * sqrt (pi / 15.0_rk)
      if (equatorial_symmetry) gw_factor_1 = gw_factor_1 * 2.0_rk
      gw_factor_2 = c_light ** 2 / (32.0_rk * pi * g_grav) / m_solar
#ifdef CFC_TRANSPORT2
      gw_weight = 2.5_rk
#else
      gw_weight = 0.0_rk
#endif

      return 

      END SUBROUTINE calculate_perm_aux_values

c     ==================================================================





      END MODULE gr_initial_setup

c     ==================================================================
