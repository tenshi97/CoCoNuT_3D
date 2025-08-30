c     ============================================================
      
      MODULE boundary_cfc

c     ============================================================

      USE precision
      USE size_cfc, ONLY: n

      IMPLICIT NONE
      
      SAVE

      integer (kind=ik) :: symm_v_r (1 : 3)
      integer (kind=ik) :: symm_v_theta_upper (1 : 3),
     &     symm_v_theta_lower (1 : 3)
      integer (kind=ik) :: symm_v_theta_lower_eq (1 : 3)
      integer (kind=ik) :: symm_v_theta_lower_noneq (1 : 3)

      integer (kind=ik) :: symm_m_r (1 : 5)
      integer (kind=ik) :: symm_m_theta_upper (1 : 5)
      integer (kind=ik) :: symm_m_theta_lower (1 : 5)

      data symm_v_r / - 1, - 1, - 1 /
      data symm_v_theta_upper / + 1, - 1, - 1 /
      data symm_v_theta_lower_eq / + 1, - 1, + 1 /
      data symm_v_theta_lower_noneq / + 1, - 1, - 1 /

      data symm_m_r / + 1, + 1, - 1, + 1, + 1 /
      data symm_m_theta_upper / + 1, + 1, + 1, - 1, + 1 /
      data symm_m_theta_lower / + 1, + 1, + 1, - 1, + 1/

      END MODULE boundary_cfc

c     ============================================================
