c     ==================================================================
      
      MODULE conserved_cfc

c     ==================================================================

      USE precision

      IMPLICIT NONE

      SAVE

      real (kind=rk), allocatable ::
     &     d_cap_hat (:,:,:),
     &     s_1_hat   (:,:,:),
     &     s_2_hat   (:,:,:),
     &     s_3_hat   (:,:,:),
     &     tau_hat   (:,:,:),
     &     d_cap_xnu_hat (:,:,:,:)

      real (kind=rk), allocatable ::
     &     d_cap_hat_old (:,:,:),
     &     s_1_hat_old   (:,:,:),
     &     s_2_hat_old   (:,:,:),
     &     s_3_hat_old   (:,:,:),
     &     tau_hat_old   (:,:,:),
     &     d_cap_xnu_hat_old (:,:,:,:)

#ifdef CFC_MHD
      real (kind=rk), allocatable :: b_cap_1_hat (:,:,:)
      real (kind=rk), allocatable :: b_cap_2_hat (:,:,:)
      real (kind=rk), allocatable :: b_cap_3_hat (:,:,:)

      real (kind=rk), allocatable :: b_cap_1_hat_old (:,:,:)
      real (kind=rk), allocatable :: b_cap_2_hat_old (:,:,:)
      real (kind=rk), allocatable :: b_cap_3_hat_old (:,:,:)
      real (kind=rk), allocatable :: psi_old (:,:,:)

      real (kind=rk), allocatable :: eps_aux (:,:,:)
      real (kind=rk), allocatable :: eps_old (:,:,:)
#endif /* CFC_MHD */


      END MODULE conserved_cfc

c     ==================================================================
