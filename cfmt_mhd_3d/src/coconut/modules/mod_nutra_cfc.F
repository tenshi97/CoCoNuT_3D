c     ============================================================
      
      MODULE nutra_cfc

c     ============================================================

      USE precision
      
      IMPLICIT NONE
      
      SAVE

      real (kind=rk), allocatable ::
     &     enu (:,:,:), fnu (:,:,:), pnu (:,:,:)

      END MODULE nutra_cfc

c     ============================================================
