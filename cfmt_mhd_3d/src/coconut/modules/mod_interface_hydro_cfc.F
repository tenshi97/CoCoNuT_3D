c     ============================================================
      
      MODULE interface_hydro_cfc

c     ============================================================

      USE precision
      USE size_cfc, ONLY: m, n, o
      USE configure

      IMPLICIT NONE
      
      SAVE

      real (kind=rk), allocatable ::
     &    rho_if (:,:,:,:),
     &    v_1_if (:,:,:,:),
     &    v_2_if (:,:,:,:),
     &    v_3_if (:,:,:,:),
     &    eps_if (:,:,:,:),
     &    xnu_if (:,:,:,:,:)

      real (kind=rk), allocatable ::
     &    p_if (:,:,:,:),
     &    h_if (:,:,:,:),
     &    gamm_if (:,:,:,:),
     &    v_squared_if (:,:,:,:),
     &    c_sound_squared_if (:,:,:,:)

#ifdef CFC_MHD
      real (kind=rk), allocatable ::
     &    b_1_if (:,:,:,:),
     &    b_2_if (:,:,:,:),
     &    b_3_if (:,:,:,:),
     &    e_1_if (:,:,:,:),
     &    e_2_if (:,:,:,:),
     &    e_3_if (:,:,:,:),
     &    psi_if (:,:,:,:)
#endif /* CFC_MHD */

      END MODULE interface_hydro_cfc

c     ============================================================
