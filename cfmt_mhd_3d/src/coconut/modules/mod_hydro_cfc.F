c     ==================================================================
      
      MODULE hydro_primitives_cfc

c     ==================================================================

      USE precision

      IMPLICIT NONE

      SAVE

      real (kind=rk), allocatable :: rho (:,:,:)
      real (kind=rk), allocatable :: v_1 (:,:,:)
      real (kind=rk), allocatable :: v_2 (:,:,:)
      real (kind=rk), allocatable :: v_3 (:,:,:)
      real (kind=rk), allocatable :: eps (:,:,:)

      real (kind=rk), allocatable :: p (:,:,:)
      real (kind=rk), allocatable :: v_squared (:,:,:)
      real (kind=rk), allocatable :: w (:,:,:)
      real (kind=rk), allocatable :: h (:,:,:)
      real (kind=rk), allocatable :: c_sound_squared (:,:,:)
      real (kind=rk), allocatable :: t (:,:,:)
      real (kind=rk), allocatable :: entropy (:,:,:)
      real (kind=rk), allocatable :: xnnu  (:,:,:,:)
      real (kind=rk), allocatable :: gamm  (:,:,:)
      real (kind=rk), allocatable :: cpot (:,:,:,:)

#ifdef CFC_MHD
      real (kind=rk), allocatable :: b_1 (:,:,:)
      real (kind=rk), allocatable :: b_2 (:,:,:)
      real (kind=rk), allocatable :: b_3 (:,:,:)
      real (kind=rk), allocatable :: e_1 (:,:,:)
      real (kind=rk), allocatable :: e_2 (:,:,:)
      real (kind=rk), allocatable :: e_3 (:,:,:)
      real (kind=rk), allocatable :: psi (:,:,:)
#endif /* CFC_MHD */

      END MODULE hydro_primitives_cfc

c     ==================================================================




c     ==================================================================
      
      MODULE shock_cfc

c     ==================================================================

      IMPLICIT NONE

      SAVE

      logical, allocatable :: lshock (:,:,:)

      END MODULE shock_cfc

c     ==================================================================
