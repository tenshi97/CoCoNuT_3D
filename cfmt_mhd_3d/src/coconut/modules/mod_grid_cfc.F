c     ============================================================
      
      MODULE grid_cfc

c     ============================================================

      USE precision
      USE size_cfc, ONLY: m, n, o
      
      IMPLICIT NONE
      
      SAVE

      real (kind=rk), allocatable :: r (:), r_if (:),
     &       delta_r (:), delta_r_if (:)

      real (kind=rk), allocatable :: theta (:), theta_if (:)
      real(kind=rk) :: delta_theta, delta_theta_if

      real (kind=rk), allocatable :: varphi (:), varphi_if (:)
      real(kind=rk) ::       delta_varphi, delta_varphi_if

      CONTAINS 

c     ============================================================

      SUBROUTINE allocate_grid_cfc

c     ============================================================
      
      IMPLICIT NONE

      allocate(r (- 3 : m + 4), r_if (- 3 : m + 3),
     &         delta_r (- 2 : m + 3), delta_r_if (- 3 : m + 3),
     &         theta (- 3 : n + 4), theta_if (- 1 : n + 1),
     &         varphi (- 3 : o + 4), varphi_if (- 1 : o + 1))

      END SUBROUTINE allocate_grid_cfc

c     ============================================================


      END MODULE grid_cfc

c     ============================================================
