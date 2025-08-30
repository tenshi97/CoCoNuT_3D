c     ==================================================================

      MODULE modpot_cfc

c     ==================================================================

      USE precision

      IMPLICIT NONE

      SAVE

      real (kind=rk), allocatable :: rho_grav (:,:,:)

      real (kind=rk), allocatable :: rho_tangential_1_spherical (:)
      real (kind=rk), allocatable :: rho_tangential_2_spherical (:)
      real (kind=rk), allocatable :: rho_mod_spherical (:)
      real (kind=rk), allocatable :: eps_spherical (:)
      real (kind=rk), allocatable :: p_spherical (:)
      real (kind=rk), allocatable :: h_spherical (:)
      real (kind=rk), allocatable :: v_1_spherical (:)
      real (kind=rk), allocatable :: pnu_spherical (:)

      real (kind=rk), allocatable :: m_spherical_newton (:)
      real (kind=rk), allocatable :: m_spherical_newton_if (:)
      real (kind=rk), allocatable :: m_spherical_tov (:)
      real (kind=rk), allocatable :: m_spherical_tov_if (:)

      real (kind=rk), allocatable :: potential_spherical_newton (:)
      real (kind=rk), allocatable :: potential_spherical_tov (:)
      real (kind=rk), allocatable :: delta_potential_spherical (:)

      real (kind=rk) :: energy_factor
      real (kind=rk) :: velocity_factor

      real (kind=rk) :: gamma, gamma_mass, gamma_tangential
      real (kind=rk) :: v_tangential_1, v_tangential_2

      END MODULE modpot_cfc

c     ==================================================================
