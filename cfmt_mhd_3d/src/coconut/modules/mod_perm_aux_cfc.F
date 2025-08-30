c     ==================================================================

      MODULE perm_aux_cfc

c     ==================================================================

      USE precision
      USE size_cfc, ONLY: m, n, o

      IMPLICIT NONE
      
      SAVE

      real (kind=rk) :: cos_dtheta
      real (kind=rk) :: sin_dtheta
      real (kind=rk) :: cos_dvarphi
      real (kind=rk) :: sin_dvarphi
      real (kind=rk) :: delta_theta_inv
      real (kind=rk) :: delta_ttheta_inv
      real (kind=rk) :: delta_theta_if_inv
      real (kind=rk), allocatable :: delta_theta_mod_inv (:)
      real (kind=rk) :: delta_varphi_inv
      real (kind=rk) :: delta_varphi_if_inv
      real (kind=rk), allocatable :: amp (:)
      real (kind=rk), allocatable :: amp_inv (:)
      real (kind=rk), allocatable :: aamp (:)
      real (kind=rk), allocatable :: rr (:)
      real (kind=rk), allocatable :: rrr (:)
      real (kind=rk), allocatable :: r_inv (:)
      real (kind=rk), allocatable :: r_if_inv (:)
      real (kind=rk), allocatable :: rr_inv (:)
      real (kind=rk), allocatable :: rr_if (:)
      real (kind=rk), allocatable :: delta_r_if_inv (:)
      real (kind=rk), allocatable :: delta_r_inv (:)
      real (kind=rk), allocatable :: delta_r_mod_inv (:)
      real (kind=rk), allocatable :: delta_r_amp_inv (:)
      real (kind=rk), allocatable :: cotan_theta (:)
      real (kind=rk), allocatable :: ccotan_theta (:)
      real (kind=rk), allocatable :: sinus_theta (:)
      real (kind=rk), allocatable :: cosin_theta (:)
      real (kind=rk), allocatable :: ssinus_theta (:)
      real (kind=rk), allocatable :: ccosin_theta (:)
      real (kind=rk), allocatable :: sinus_theta_inv (:)
      real (kind=rk), allocatable :: ssinus_theta_inv (:)
      real (kind=rk), allocatable :: sinus_theta_if (:)
      real (kind=rk), allocatable :: ssinus_theta_if (:)
      real (kind=rk), allocatable :: sinus_theta_inv_if (:)
      real (kind=rk), allocatable :: cosin_theta_if (:)
      real (kind=rk), allocatable :: sinuscosin_theta (:)
      real (kind=rk), allocatable :: sinus_varphi (:)
      real (kind=rk), allocatable :: ssinus_varphi (:)
      real (kind=rk), allocatable :: cosin_varphi (:)
      real (kind=rk), allocatable :: ccosin_varphi (:)
      real (kind=rk), allocatable :: sinuscosin_varphi (:)
      real (kind=rk), allocatable :: delta_rrr (:)
      real (kind=rk), allocatable :: delta_cosin_theta (:)
#ifdef CFC_MHD
      real (kind=rk) :: sumdvy0 (0:3)
#endif

      CONTAINS


c     ==================================================================

      SUBROUTINE allocate_perm_aux_cfc

c     ==================================================================

      IMPLICIT NONE

      allocate(delta_theta_mod_inv (1 : n),
     &         amp (1 : m + 3),
     &         amp_inv (1 : m + 1),
     &         aamp (1 : m + 1),
     &         rr (- 3 : m + 4),
     &         rrr (1 : m),
     &         r_inv (- 3 : m + 4),
     &         r_if_inv (0 : m + 1),
     &         rr_inv (- 3 : m + 4),
     &         rr_if (0 : m + 1),
     &         delta_r_if_inv (- 1 : m + 1),
     &         delta_r_inv (1 : m),
     &         delta_r_mod_inv (1 : m),
     &         delta_r_amp_inv (1 : m + 1),
     &         cotan_theta (0 : n + 1),
     &         ccotan_theta (1 : n),
     &         sinus_theta (- 3 : n + 4),
     &         cosin_theta (- 3 : n + 4),
     &         ssinus_theta (1 : n),
     &         ccosin_theta (1 : n),
     &         sinus_theta_inv (- 3 : n + 4),
     &         ssinus_theta_inv (- 3 : n + 4),
     &         sinus_theta_if (- 1 : n + 1),
     &         ssinus_theta_if (0 : n),
     &         sinus_theta_inv_if (0: n),
     &         cosin_theta_if (0 : n),
     &         sinuscosin_theta (1:n),
     &         sinus_varphi (0 : o + 1),
     &         ssinus_varphi (1 : o),
     &         cosin_varphi (0 : o + 1),
     &         ccosin_varphi (1 : o),
     &         sinuscosin_varphi (1 : o),
     &         delta_rrr (1 : m + 1),
     &         delta_cosin_theta (1 : n) )

      return

      END SUBROUTINE allocate_perm_aux_cfc

c     ==================================================================


      END MODULE perm_aux_cfc

c     ==================================================================
