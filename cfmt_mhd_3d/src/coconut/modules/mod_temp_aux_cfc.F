c     ==================================================================
      
      MODULE temp_aux_cfc

c     ==================================================================

      USE precision
      USE size_cfc, ONLY: m, n, o
      
      IMPLICIT NONE
      
      SAVE

      real (kind=rk), allocatable :: rho_mod (:,:,:)
      real (kind=rk), allocatable :: rho_gw (:,:,:)
      real (kind=rk), allocatable :: v_gw_r (:,:,:)
      real (kind=rk), allocatable :: v_gw_theta (:,:,:)
      real (kind=rk), allocatable :: v_gw_varphi (:,:,:)

      END MODULE temp_aux_cfc

c     ==================================================================
