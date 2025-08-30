c     ==================================================================
      
      MODULE ppm_cfc

c     ==================================================================

      USE precision
      USE size_cfc
      
      IMPLICIT NONE
      
      SAVE

      real (kind=rk), allocatable :: delta_q (:,:,:)
      real (kind=rk), allocatable :: f_tilde (:,:,:)
      real (kind=rk), allocatable :: delta_p (:,:,:)
      real (kind=rk), allocatable :: tmp_ppm (:,:,:)

      real (kind=rk), allocatable :: c_1_r (:), c_2_r (:), 
     &     c_3_r (:), c_4_r (:), c_5_r (:)
      real (kind=rk) :: c_1_theta, c_2_theta, c_3_theta,
     &     c_4_theta, c_5_theta
      real (kind=rk) :: c_1_varphi, c_2_varphi, c_3_varphi,
     &     c_4_varphi, c_5_varphi


      CONTAINS


c     ==================================================================

      SUBROUTINE allocate_ppm_cfc

c     ==================================================================

      IMPLICIT NONE

      integer :: istat

      allocate(
     &     delta_q (- 1 : m + 2, n_s - 2 : n_e + 2, o_s - 2 : o_e + 2),
     &     f_tilde (- 1 : m + 2, n_s - 2 : n_e + 2, o_s - 2 : o_e + 2),
     &     delta_p (- 1 : m + 2, n_s - 2 : n_e + 2, o_s - 2 : o_e + 2),
     &     tmp_ppm (- 2 : m + 2, n_s - 3 : n_e + 2, o_s - 3 : o_e + 2),
     &     stat=istat)

      allocate(c_1_r (- 1 : m + 2), c_2_r (- 1 : m + 2),
     &         c_3_r (- 1 : m + 1), c_4_r (- 1 : m + 1),
     &         c_5_r (- 1 : m + 1), stat=istat)
      
      return

      END SUBROUTINE allocate_ppm_cfc

c     ==================================================================


      END MODULE ppm_cfc

c     ==================================================================
