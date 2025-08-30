c     ==================================================================

      MODULE metric_cfc

c     ==================================================================

      USE precision
      USE size_cfc, ONLY: m, n, o

      IMPLICIT NONE

      SAVE

      real (kind=rk), allocatable :: phi       (:,:,:)
      real (kind=rk), allocatable :: alpha     (:,:,:)
      real (kind=rk), allocatable :: beta_up_1 (:,:,:)
      real (kind=rk), allocatable :: beta_up_2 (:,:,:)
      real (kind=rk), allocatable :: beta_up_3 (:,:,:)
      real (kind=rk), allocatable :: g_11 (:,:,:)
      real (kind=rk), allocatable :: g_22 (:,:,:)
      real (kind=rk), allocatable :: g_33 (:,:,:)
      real (kind=rk), allocatable :: g_up_11 (:,:,:)
      real (kind=rk), allocatable :: g_up_22 (:,:,:)
      real (kind=rk), allocatable :: g_up_33 (:,:,:)

      real (kind=rk), allocatable :: phi_r       (:,:,:)
      real (kind=rk), allocatable :: alpha_r     (:,:,:)
      real (kind=rk), allocatable :: beta_up_1_r (:,:,:)
      real (kind=rk), allocatable :: beta_up_2_r (:,:,:)
      real (kind=rk), allocatable :: beta_up_3_r (:,:,:)
      real (kind=rk), allocatable :: g_11_r (:,:,:)
      real (kind=rk), allocatable :: g_22_r (:,:,:)
      real (kind=rk), allocatable :: g_33_r (:,:,:)
      real (kind=rk), allocatable :: g_up_11_r (:,:,:)
      real (kind=rk), allocatable :: g_up_22_r (:,:,:)
      real (kind=rk), allocatable :: g_up_33_r (:,:,:)

      real (kind=rk), allocatable :: phi_theta       (:,:,:)
      real (kind=rk), allocatable :: alpha_theta     (:,:,:)
      real (kind=rk), allocatable :: beta_up_1_theta (:,:,:)
      real (kind=rk), allocatable :: beta_up_2_theta (:,:,:)
      real (kind=rk), allocatable :: beta_up_3_theta (:,:,:)
      real (kind=rk), allocatable :: g_11_theta (:,:,:)
      real (kind=rk), allocatable :: g_22_theta (:,:,:)
      real (kind=rk), allocatable :: g_33_theta (:,:,:)
      real (kind=rk), allocatable :: g_up_11_theta (:,:,:)
      real (kind=rk), allocatable :: g_up_22_theta (:,:,:)
      real (kind=rk), allocatable :: g_up_33_theta (:,:,:)

      real (kind=rk), allocatable :: phi_varphi       (:,:,:)
      real (kind=rk), allocatable :: alpha_varphi     (:,:,:)
      real (kind=rk), allocatable :: beta_up_1_varphi (:,:,:)
      real (kind=rk), allocatable :: beta_up_2_varphi (:,:,:)
      real (kind=rk), allocatable :: beta_up_3_varphi (:,:,:)
      real (kind=rk), allocatable :: g_11_varphi (:,:,:)
      real (kind=rk), allocatable :: g_22_varphi (:,:,:)
      real (kind=rk), allocatable :: g_33_varphi (:,:,:)
      real (kind=rk), allocatable :: g_up_11_varphi (:,:,:)
      real (kind=rk), allocatable :: g_up_22_varphi (:,:,:)
      real (kind=rk), allocatable :: g_up_33_varphi (:,:,:)

      real (kind=rk), allocatable :: k_cap_hat (:,:,:,:,:)
      real (kind=rk), allocatable :: k_cap_hat_squared (:,:,:)

      real (kind=rk), allocatable :: sqrt_gamma (:,:,:)
      real (kind=rk), allocatable :: sqrt_gamma_inv (:,:,:)
      real (kind=rk), allocatable :: sqrt_gamma_r (:,:,:)
      real (kind=rk), allocatable :: sqrt_gamma_theta (:,:,:)
      real (kind=rk), allocatable :: sqrt_gamma_varphi (:,:,:)

      real (kind=rk), allocatable :: beta (:,:,:,:)
      real (kind=rk), allocatable :: w_tmp (:,:,:,:)
      real (kind=rk), allocatable :: dw_tmp_dr (:,:,:,:)
      real (kind=rk), allocatable :: dw_tmp_dtheta (:,:,:,:)

      real (kind=rk), allocatable :: s_cap_star      (:,:,:,:)
      real (kind=rk), allocatable :: s_cap           (:,:,:,:)
      real (kind=rk), allocatable :: scalar_tmp      (:,:,:)
      real (kind=rk), allocatable :: rho_hat_phi     (:,:,:)
      real (kind=rk), allocatable :: rho_hat_phi_1   (:,:,:)
      real (kind=rk), allocatable :: rho_hat_phi_2   (:,:,:)
      real (kind=rk), allocatable :: rho_hat_alpha   (:,:,:)
      real (kind=rk), allocatable :: rho_hat_alpha_1 (:,:,:)

      logical metric_exit

      integer metric_iteration
      integer initial_metric_iteration

!$OMP THREADPRIVATE (metric_iteration,initial_metric_iteration)

      real (kind=rk), allocatable :: phi_potential (:,:,:)

      real (kind=rk), allocatable :: dphi_potential_dr (:,:,:)
      real (kind=rk), allocatable :: dphi_potential_dtheta (:,:,:)
      real (kind=rk), allocatable :: dphi_potential_dvarphi (:,:,:)

      real (kind=rk), allocatable ::
     &     dphi_potential_alt_dr (:,:,:)
      real (kind=rk), allocatable ::
     &     dphi_potential_alt_dtheta (:,:,:)
      real (kind=rk), allocatable ::
     &     dphi_potential_alt_dvarphi (:,:,:)


      END MODULE metric_cfc

c     ==================================================================
