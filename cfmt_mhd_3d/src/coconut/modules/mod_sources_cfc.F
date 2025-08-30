c     ==================================================================
      
      MODULE sources_cfc

c     ==================================================================

      USE precision
      USE size_cfc, ONLY: m, n, o
      
      IMPLICIT NONE
      
      SAVE

      real (kind=rk), allocatable :: s_total (:,:,:,:)


      real (kind=rk), allocatable :: dln_sqrt_gamma_dr (:,:,:)
      real (kind=rk), allocatable :: dln_sqrt_gamma_dtheta (:,:,:)
      real (kind=rk), allocatable :: dln_sqrt_gamma_dvarphi (:,:,:)
      real (kind=rk), allocatable :: beta_up_k_k (:,:,:)

      real (kind=rk), allocatable :: dphi_dt_o_phi (:,:,:)

      real (kind=rk), allocatable :: alpha_dalpha_dr (:,:,:)
      real (kind=rk), allocatable :: alpha_dalpha_dtheta (:,:,:)
      real (kind=rk), allocatable :: alpha_dalpha_dvarphi (:,:,:)

      real (kind=rk), allocatable :: dalpha_dr (:,:,:)
      real (kind=rk), allocatable :: dalpha_dtheta (:,:,:)
      real (kind=rk), allocatable :: dalpha_dvarphi (:,:,:)

      real (kind=rk), allocatable :: dg_00_dr (:,:,:)
      real (kind=rk), allocatable :: dg_11_dr (:,:,:)
      real (kind=rk), allocatable :: dg_22_dr (:,:,:)
      real (kind=rk), allocatable :: dg_33_dr (:,:,:)
      real (kind=rk), allocatable :: dg_00_dtheta (:,:,:)
      real (kind=rk), allocatable :: dg_11_dtheta (:,:,:)
      real (kind=rk), allocatable :: dg_22_dtheta (:,:,:)
      real (kind=rk), allocatable :: dg_33_dtheta (:,:,:)
      real (kind=rk), allocatable :: dg_00_dvarphi (:,:,:)
      real (kind=rk), allocatable :: dg_11_dvarphi (:,:,:)
      real (kind=rk), allocatable :: dg_22_dvarphi (:,:,:)
      real (kind=rk), allocatable :: dg_33_dvarphi (:,:,:)

      real (kind=rk), allocatable :: dg_01_dr (:,:,:)
      real (kind=rk), allocatable :: dg_02_dr (:,:,:)
      real (kind=rk), allocatable :: dg_03_dr (:,:,:)
      real (kind=rk), allocatable :: dg_01_dtheta (:,:,:)
      real (kind=rk), allocatable :: dg_02_dtheta (:,:,:)
      real (kind=rk), allocatable :: dg_03_dtheta (:,:,:)
      real (kind=rk), allocatable :: dg_01_dvarphi (:,:,:)
      real (kind=rk), allocatable :: dg_02_dvarphi (:,:,:)
      real (kind=rk), allocatable :: dg_03_dvarphi (:,:,:)

      real (kind=rk), allocatable :: excrv (:,:,:,:,:)
      real (kind=rk), allocatable :: excrv2 (:,:,:)
      real (kind=rk), allocatable :: excrv_beta (:,:,:)

      END MODULE sources_cfc

c     ==================================================================
