MODULE cfc_restart_vars

  IMPLICIT NONE

  SAVE

  real, allocatable ::  d_cap_hat (:,:,:),  s_1_hat (:,:,:), &
                        s_2_hat (:,:,:), s_3_hat (:,:,:),    &
                        tau_hat (:,:,:), xnnucfc (:,:,:,:), &
                        phi (:,:,:), alpha (:,:,:), beta_up_1 (:,:,:), &
                        beta_up_2 (:,:,:), beta_up_3 (:,:,:), &
                        enu(:,:,:),fnu(:,:,:)


  real, allocatable ::  &
           d_cap_hat_new (:,:,:), &
           s_1_hat_new (:,:,:), &
           s_2_hat_new (:,:,:), &
           s_3_hat_new (:,:,:), &
           tau_hat_new (:,:,:),     &
           xnnucfc_new (:,:,:,:), &
           phi_new (:,:,:), &
           alpha_new (:,:,:), &
           beta_up_1_new (:,:,:), &
           beta_up_2_new (:,:,:), &
           beta_up_3_new (:,:,:), &
           enu_new(:,:,:),fnu_new(:,:,:)

  real :: m_grav_1_ini,m_rest_ini

END MODULE cfc_restart_vars
