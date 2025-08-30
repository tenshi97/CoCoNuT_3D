c     ==================================================================
      
      MODULE parameters_cfc

c     ==================================================================

      USE precision
      
      IMPLICIT NONE

      SAVE

      real (kind=rk) :: pi, equatorial_symmetry_factor
      real (kind=rk), parameter :: m_solar = 1.98892e+33_rk

      integer (kind=ik) :: iteration, first_iteration, last_iteration

      real (kind=rk) :: conversion_tolerance, metric_tolerance
      real (kind=rk) :: relaxation_factor (1 : 5)
      
      real (kind=rk) :: delta_t, delta_t_old, delta_t_old_old,
     &     cfl_factor, t_total, delta_t_refine, initial_time_step,
     &     delta_t_r_inv, delta_t_theta_inv, delta_t_varphi_inv
      
      integer (kind=ik) :: metric_step_1, metric_step_2, metric_step_3,
     &        metric_step, metric_extrapolation, l_ex,
     &        newtonian_potential_form

      real (kind=rk) :: rho_nuc_12, rho_nuc_23, rho_nuc_34,
     &     rho_nuc_fraction

      integer (kind=ik) :: rk_order, rk_index

      real (kind=rk) :: eta_1, eta_2, epsilon_1, k_cap_0, epsilon_2_r,
     &       epsilon_2_theta, epsilon_2_varphi, omega_1, omega_2

      real (kind=rk) :: m_rest, m_rest_ini, 
     &       m_grav_1, m_grav_1_ini

      logical initial_step, rho_nuc_thr_passed,
     &       rho_nuc_thr_passed_old, metric_calculated,
     &       multiple_bounce,equatorial_symmetry

      real (kind=rk) :: g_grav, c_light,
     &       p_geom_factor, p_geom_factor_inv,
     &       rho_geom_factor, rho_geom_factor_inv,
     &       gw_factor_1, gw_factor_2, gw_weight

      real (kind=rk) :: a_quad_plus_equator_1
      real (kind=rk) :: a_quad_plus_equator_1_offset_1
      real (kind=rk) :: a_quad_plus_equator_1_offset_2
      real (kind=rk) :: a_quad_cross_equator_1
      real (kind=rk) :: a_quad_cross_equator_1_offset_1
      real (kind=rk) :: a_quad_cross_equator_1_offset_2
      real (kind=rk) :: a_quad_plus_pole_1
      real (kind=rk) :: a_quad_plus_pole_1_offset_1
      real (kind=rk) :: a_quad_plus_pole_1_offset_2
      real (kind=rk) :: a_quad_cross_pole_1
      real (kind=rk) :: a_quad_cross_pole_1_offset_1
      real (kind=rk) :: a_quad_cross_pole_1_offset_2
      real (kind=rk) :: a_quad_plus_equator_2
      real (kind=rk) :: a_quad_cross_equator_2
      real (kind=rk) :: a_quad_plus_pole_2
      real (kind=rk) :: a_quad_cross_pole_2
      real (kind=rk) :: ddi_cap_ddt_1
      real (kind=rk) :: ddi_cap_ddt_1_offset_1
      real (kind=rk) :: ddi_cap_ddt_1_offset_2
      real (kind=rk) :: ddi_cap_ddt_1_old
      real (kind=rk) :: ddi_cap_ddt_1_old_old
      real (kind=rk) :: ddi_cap_ddt_2
      real (kind=rk) :: ddi_cap_ddt_2_old
      real (kind=rk) :: ddi_cap_ddt_2_old_old
      real (kind=rk) :: ddi_cap_ddt_2_old_old_old
      real (kind=rk) :: ddi_cap_11_ddt_1
      real (kind=rk) :: ddi_cap_11_ddt_1_offset_1
      real (kind=rk) :: ddi_cap_11_ddt_1_offset_2
      real (kind=rk) :: ddi_cap_22_ddt_1
      real (kind=rk) :: ddi_cap_22_ddt_1_offset_1
      real (kind=rk) :: ddi_cap_22_ddt_1_offset_2
      real (kind=rk) :: ddi_cap_33_ddt_1
      real (kind=rk) :: ddi_cap_33_ddt_1_offset_1
      real (kind=rk) :: ddi_cap_33_ddt_1_offset_2
      real (kind=rk) :: ddi_cap_12_ddt_1
      real (kind=rk) :: ddi_cap_12_ddt_1_offset_1
      real (kind=rk) :: ddi_cap_12_ddt_1_offset_2
      real (kind=rk) :: ddi_cap_13_ddt_1
      real (kind=rk) :: ddi_cap_13_ddt_1_offset_1
      real (kind=rk) :: ddi_cap_13_ddt_1_offset_2
      real (kind=rk) :: ddi_cap_23_ddt_1
      real (kind=rk) :: ddi_cap_23_ddt_1_offset_1
      real (kind=rk) :: ddi_cap_23_ddt_1_offset_2
      real (kind=rk) :: ddi_cap_11_ddt_1_old
      real (kind=rk) :: ddi_cap_22_ddt_1_old
      real (kind=rk) :: ddi_cap_33_ddt_1_old
      real (kind=rk) :: ddi_cap_12_ddt_1_old
      real (kind=rk) :: ddi_cap_13_ddt_1_old
      real (kind=rk) :: ddi_cap_23_ddt_1_old
      real (kind=rk) :: ddi_cap_11_ddt_1_old_old
      real (kind=rk) :: ddi_cap_22_ddt_1_old_old
      real (kind=rk) :: ddi_cap_33_ddt_1_old_old
      real (kind=rk) :: ddi_cap_12_ddt_1_old_old
      real (kind=rk) :: ddi_cap_13_ddt_1_old_old
      real (kind=rk) :: ddi_cap_23_ddt_1_old_old
      real (kind=rk) :: ddi_cap_11_ddt_2
      real (kind=rk) :: ddi_cap_22_ddt_2
      real (kind=rk) :: ddi_cap_33_ddt_2
      real (kind=rk) :: ddi_cap_12_ddt_2
      real (kind=rk) :: ddi_cap_13_ddt_2
      real (kind=rk) :: ddi_cap_23_ddt_2
      real (kind=rk) :: ddi_cap_11_ddt_2_old
      real (kind=rk) :: ddi_cap_22_ddt_2_old
      real (kind=rk) :: ddi_cap_33_ddt_2_old
      real (kind=rk) :: ddi_cap_12_ddt_2_old
      real (kind=rk) :: ddi_cap_13_ddt_2_old
      real (kind=rk) :: ddi_cap_23_ddt_2_old
      real (kind=rk) :: ddi_cap_11_ddt_2_old_old
      real (kind=rk) :: ddi_cap_22_ddt_2_old_old
      real (kind=rk) :: ddi_cap_33_ddt_2_old_old
      real (kind=rk) :: ddi_cap_12_ddt_2_old_old
      real (kind=rk) :: ddi_cap_13_ddt_2_old_old
      real (kind=rk) :: ddi_cap_23_ddt_2_old_old
      real (kind=rk) :: ddi_cap_11_ddt_2_old_old_old
      real (kind=rk) :: ddi_cap_22_ddt_2_old_old_old
      real (kind=rk) :: ddi_cap_33_ddt_2_old_old_old
      real (kind=rk) :: ddi_cap_12_ddt_2_old_old_old
      real (kind=rk) :: ddi_cap_13_ddt_2_old_old_old
      real (kind=rk) :: ddi_cap_23_ddt_2_old_old_old
      real (kind=rk) :: di_cap_dt_2
      real (kind=rk) :: di_cap_dt_2_old
      real (kind=rk) :: di_cap_11_dt_2
      real (kind=rk) :: di_cap_22_dt_2
      real (kind=rk) :: di_cap_33_dt_2
      real (kind=rk) :: di_cap_12_dt_2
      real (kind=rk) :: di_cap_13_dt_2
      real (kind=rk) :: di_cap_23_dt_2
      real (kind=rk) :: di_cap_11_dt_2_old
      real (kind=rk) :: di_cap_22_dt_2_old
      real (kind=rk) :: di_cap_33_dt_2_old
      real (kind=rk) :: di_cap_12_dt_2_old
      real (kind=rk) :: di_cap_13_dt_2_old
      real (kind=rk) :: di_cap_23_dt_2_old
      real (kind=rk) :: dddi_cap_dddt_1
      real (kind=rk) :: dddi_cap_dddt_1_old
      real (kind=rk) :: dddi_cap_11_dddt_1
      real (kind=rk) :: dddi_cap_22_dddt_1
      real (kind=rk) :: dddi_cap_33_dddt_1
      real (kind=rk) :: dddi_cap_12_dddt_1
      real (kind=rk) :: dddi_cap_13_dddt_1
      real (kind=rk) :: dddi_cap_23_dddt_1
      real (kind=rk) :: dddi_cap_11_dddt_1_old
      real (kind=rk) :: dddi_cap_22_dddt_1_old
      real (kind=rk) :: dddi_cap_33_dddt_1_old
      real (kind=rk) :: dddi_cap_12_dddt_1_old
      real (kind=rk) :: dddi_cap_13_dddt_1_old
      real (kind=rk) :: dddi_cap_23_dddt_1_old
      real (kind=rk) :: dddi_cap_dddt_2
      real (kind=rk) :: dddi_cap_dddt_2_old
      real (kind=rk) :: dddi_cap_11_dddt_2
      real (kind=rk) :: dddi_cap_22_dddt_2
      real (kind=rk) :: dddi_cap_33_dddt_2
      real (kind=rk) :: dddi_cap_12_dddt_2
      real (kind=rk) :: dddi_cap_13_dddt_2
      real (kind=rk) :: dddi_cap_23_dddt_2
      real (kind=rk) :: dddi_cap_11_dddt_2_old
      real (kind=rk) :: dddi_cap_22_dddt_2_old
      real (kind=rk) :: dddi_cap_33_dddt_2_old
      real (kind=rk) :: dddi_cap_12_dddt_2_old
      real (kind=rk) :: dddi_cap_13_dddt_2_old
      real (kind=rk) :: dddi_cap_23_dddt_2_old
      real (kind=rk) :: e_grav_rad_1
      real (kind=rk) :: e_grav_rad_2

      END MODULE parameters_cfc

c     ==================================================================
