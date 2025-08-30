c     ==================================================================
      
      MODULE enecons_cfc

c     ==================================================================

      USE precision

      IMPLICIT NONE
      
      SAVE
      
      real (kind=rk), allocatable :: phi_potential_alt (:,:,:),
     &     alpha_alt (:,:,:)
      
      END MODULE enecons_cfc

c     ==================================================================
