c     ============================================================

      MODULE fluxes_cfc

c     ============================================================

      USE precision

      IMPLICIT NONE

      SAVE

      real (kind=rk), allocatable ::
     &     f_num_if (:,:,:,:), f_num_xnu_if (:,:,:,:)

      real (kind=rk), allocatable ::
     &     f_diff (:,:,:,:), f_diff_xnu (:,:,:,:)

      real (kind=rk), allocatable ::
     &     pav_r (:,:,:),
     &     pav_theta (:,:,:),
     &     pav_varphi (:,:,:)

#ifdef CFC_MHD
      real (kind=rk), allocatable ::
     &     emf_face (:,:,:,:,:)

      real (kind=rk), allocatable ::
     &     f_diff_b (:,:,:,:)
#endif /* CFC_MHD */

      END MODULE fluxes_cfc

c     ============================================================
