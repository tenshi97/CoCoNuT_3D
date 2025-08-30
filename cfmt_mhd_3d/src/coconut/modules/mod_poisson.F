c     ==================================================================

      MODULE legendre

c     ==================================================================

      USE precision
      USE size_cfc

      IMPLICIT NONE

      SAVE

      real (kind=rk), allocatable :: p_leg(:,:)
      real (kind=rk), allocatable :: p_leg_int(:,:)
      real (kind=rk), allocatable :: p_leg_c(:,:)

      real (kind=rk), allocatable :: p1_leg(:,:)
      real (kind=rk), allocatable :: p1_leg_int(:,:)
      real (kind=rk), allocatable :: p1_leg_c(:,:)

#if defined(OPENMP_CFC) && (defined (OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP THREADPRIVATE(p_leg,p_leg_int,p_leg_c,p1_leg,p1_leg_int,p1_leg_c)
#endif

      CONTAINS


c     ==================================================================

      SUBROUTINE allocate_legendre

c     ==================================================================

      IMPLICIT NONE

      integer(kind=ik) :: istat
#if !(defined(PROGRAM_remap))
#if defined(OPENMP_CFC) && (defined (OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif
      allocate(p_leg( 0:n  ,0:nleg+2), p_leg_int( 1:n  ,0:nleg+1),
     &         p_leg_c(-1:n+2,0:nleg+2), 
     &         p1_leg( 0:n  ,1:nleg+2), p1_leg_int( 1:n  ,1:nleg+1),
     &         p1_leg_c(-1:n+2,1:nleg+2), stat=istat)

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_CFC) && (defined (OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp end parallel
#endif
#endif
      END SUBROUTINE allocate_legendre

c     ==================================================================



      END MODULE legendre

c     ==================================================================




c     ==================================================================

      MODULE mod_hilfsgroessen_r

c     ==================================================================

      USE size_cfc

      IMPLICIT NONE

      SAVE

      real (kind=rk), allocatable  :: r_hoch_l(:,:)
      real (kind=rk), allocatable:: r_hoch_l_r(:,:)

      real (kind=rk), allocatable :: gew_r_in(:,:), gew_r_ex(:,:)


      real, allocatable :: mom_rho(:,:), mom_rho2(:,:)
      real, allocatable :: phi_in(:,:), phi_ex(:,:)
      real, allocatable :: mom_s(:,:,:), mom_s2(:,:,:)
      real, allocatable :: phi_inv(:,:,:), phi_exv(:,:,:)

#if defined(OPENMP_CFC) && (defined (OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP THREADPRIVATE(r_hoch_l,r_hoch_l_r,gew_r_in,gew_r_ex)
#endif

      CONTAINS

c     ==================================================================

      SUBROUTINE allocate_hilfsgroessen_r

c     ==================================================================

      USE PRECISION

      IMPLICIT NONE

      integer(kind=ik) :: istat
#if !(defined(PROGRAM_remap))
#if defined(OPENMP_CFC) && (defined (OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif
      allocate(r_hoch_l(1:m, 0:nleg+2), r_hoch_l_r(1:m+1, 0:nleg+2),
     &         gew_r_in(1:m,0:nleg+1), gew_r_ex(1:m,0:nleg+1), 
     &         stat=istat)

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_CFC) && (defined (OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp end parallel
#endif
#endif
      

      allocate(
     &     mom_rho (1:m  ,0:nleg+1),
     &     mom_rho2(1:m  ,0:nleg+1),
     &     phi_in  (1:m+1,0:nleg+1),
     &     phi_ex  (1:m  ,0:nleg+1))

      allocate(
     &     mom_s  (1:m  ,0:nleg,4),
     &     mom_s2 (1:m  ,0:nleg,4),
     &     phi_inv(1:m+1,0:nleg,4),
     &     phi_exv(1:m  ,0:nleg,4))

      return

      END SUBROUTINE allocate_hilfsgroessen_r

c     ==================================================================


      END MODULE mod_hilfsgroessen_r

c     ==================================================================





c     ==================================================================

      MODULE mod_poisson_aux

c     ==================================================================

      USE precision

      IMPLICIT NONE

      SAVE

      real (kind=rk), allocatable :: s_x (:,:,:), s_y (:,:,:),
     &     s_z(:,:,:) ,s_tmp(:,:,:), beta_x (:,:,:), beta_y (:,:,:),
     &     beta_z(:,:,:), chi (:,:,:)


      END MODULE mod_poisson_aux

c     ==================================================================
