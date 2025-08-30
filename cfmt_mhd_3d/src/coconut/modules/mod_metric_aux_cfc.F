c     ==================================================================

      MODULE metric_aux_cfc

c     ==================================================================
c     auxiliary variables that must be scoped as SHARED
c     ------------------------------------------------------------------

      USE precision

      IMPLICIT NONE

      SAVE

      real(kind=rk), allocatable :: alpha_phi6 (:,:,:)
      real(kind=rk), allocatable :: da_phi_dr (:,:,:)
      real(kind=rk), allocatable :: da_phi_dtheta (:,:,:)
      real(kind=rk), allocatable :: da_phi_dvarphi (:,:,:)
      real(kind=rk), allocatable :: drel_tmp(:,:)
      real(kind=rk) :: drel

      END MODULE metric_aux_cfc

c     ==================================================================
