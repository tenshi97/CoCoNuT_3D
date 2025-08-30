c     ==================================================================

      MODULE induction

c     ==================================================================

      USE PRECISION

      IMPLICIT NONE

      CONTAINS

#ifdef CFC_MHD
c     ==================================================================

      SUBROUTINE emf_2d

c     ==================================================================

      USE precision

      USE size_cfc
      USE parameters_cfc
      USE conserved_cfc, ONLY: b_cap_1_hat, b_cap_2_hat, b_cap_3_hat
      USE fluxes_cfc
      USE grid_cfc
      USE hydro_primitives_cfc, ONLY: v_1, v_2, v_3
      USE metric_cfc
      USE perm_aux_cfc

      USE configure

      IMPLICIT NONE

      INTEGER (KIND=ik) :: i, j, k, jk

      REAL (KIND=rk) :: v_1_star, v_2_star, v_3_star

c     We need to set emf_face (:,:,: 1:2,3), i.e. the EMF on
c     varphi-interfaces by hand, since no fluxes in the
c     varphi-direction are computed.

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk,l)
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         k = 1

         do i = 1, m

            v_1_star = (alpha (i, j, k) * v_1 (i, j, k) *
     &           g_up_11 (i, j, k) - beta_up_1 (i, j, k))
            v_2_star = (alpha (i, j, k) * v_2 (i, j, k) *
     &           g_up_22 (i, j, k) - beta_up_2 (i, j, k)) *
     &           r (i)
            v_3_star = (alpha (i, j, k) * v_3 (i, j, k) *
     &           g_up_33 (i, j, k) - beta_up_3 (i, j, k)) *
     &           r (i) * sinus_theta (j)
            emf_face (i, j, k, 1, 3) = 
     &           0*v_2_star * b_cap_3_hat (i, j, k) -
     &           0*v_3_star * b_cap_2_hat (i, j, k)
            emf_face (i, j, k, 2, 3) = 
     &           0*v_3_star * b_cap_1_hat (i, j, k) -
     &           0*v_1_star * b_cap_3_hat (i, j, k)

         end do
      end do

      return
      
      END SUBROUTINE emf_2d

c     ==================================================================




c     ==================================================================

      SUBROUTINE calculate_db_dt (mode)

c     ==================================================================

      USE precision

      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
      USE fluxes_cfc
      USE metric_cfc
      USE perm_aux_cfc
      
      USE configure

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: mode

      INTEGER (KIND=ik) :: i, j, k, jk, ip1, im1, jp1, jm1, kp1, km1

      REAL (KIND=rk) :: emf_r (1:m, n_s-1:n_e, o_s-1:o_e)
      REAL (KIND=rk) :: emf_theta (0:m, n_s:n_e, o_s-1:o_e)
      REAL (KIND=rk) :: emf_varphi (0:m, n_s-1:n_e, o_s:o_e)

      REAL (KIND=rk) :: v_1_star, v_2_star, v_3_star

      if (config%nsdim .lt. 3) then
c     Replicate EMFs in the varphi-direction

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk,l)
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            
            emf_face (:,j,k+1,:,:) = emf_face (:,j,k,:,:) 
            emf_face (:,j,k-1,:,:) = emf_face (:,j,k,:,:) 
            
            if (j .eq. n_s) then
               emf_face (:,j-1,k+1,:,:) = emf_face (:,j-1,k,:,:) 
               emf_face (:,j-1,k-1,:,:) = emf_face (:,j-1,k,:,:) 
            end if
            if (j .eq. n_e) then
               emf_face (:,j+1,k+1,:,:) = emf_face (:,j+1,k,:,:) 
               emf_face (:,j+1,k-1,:,:) = emf_face (:,j+1,k,:,:) 
            end if
         
         end do

      end if

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk,l)
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         
c    E_r * dl
         do i = 1, m
            emf_r (i, j, k) = 0.5_rk *
     &           (emf_face (i, j  , k  , 1, 2) +
     &            emf_face (i, j  , k+1, 1, 2) + 
     &            emf_face (i, j  , k  , 1, 3) +
     &            emf_face (i, j+1, k  , 1, 3) ) *
     &           delta_r (i)
         end do

         if (j .eq. n_s) then
            do i = 1, m
               emf_r (i, j-1, k) = 0.5_rk *
     &              (emf_face (i, j-1, k  , 1, 2) +
     &               emf_face (i, j-1, k+1, 1, 2) + 
     &               emf_face (i, j-1, k  , 1, 3) +
     &               emf_face (i, j  , k  , 1, 3) ) *
     &              delta_r (i)
            end do
         end if

         if (k .eq. o_s) then
            do i = 1, m
               emf_r (i, j, k-1) = 0.5_rk *
     &              (emf_face (i, j  , k-1, 1, 2) +
     &               emf_face (i, j  , k  , 1, 2) + 
     &               emf_face (i, j  , k-1, 1, 3) +
     &               emf_face (i, j+1, k-1, 1, 3) ) *
     &              delta_r (i)
            end do
         end if

         if (j.eq. n_s .and. k .eq. o_s) then
            do i = 1, m
               emf_r (i, j-1, k-1) = 0.5_rk *
     &              (emf_face (i, j-1, k-1, 1, 2) +
     &               emf_face (i, j-1, k  , 1, 2) + 
     &               emf_face (i, j-1, k-1, 1, 3) +
     &               emf_face (i, j  , k-1, 1, 3)) *
     &              delta_r (i)
            end do
         end if

c     E_theta * dl
         do i = 0, m
            emf_theta (i, j, k) = 0.5_rk *
     &           (emf_face (i  , j, k  , 1, 1) +
     &            emf_face (i  , j, k+1, 1, 1) + 
     &            emf_face (i  , j, k  , 2, 3) +
     &            emf_face (i+1, j, k  , 2, 3)) *
     &           r_if (i) * !delta_theta
     &           2.0_rk * sin (0.5_rk * delta_theta)
         end do

         if (k .eq. o_s) then
            do i = 0, m
               emf_theta (i, j, k-1) = 0.5_rk *
     &              (emf_face (i  , j, k-1, 1, 1) +
     &               emf_face (i  , j, k  , 1, 1) + 
     &               emf_face (i  , j, k-1, 2, 3) +
     &               emf_face (i+1, j, k-1, 2, 3)) *
     &           r_if (i) * !delta_theta
     &           2.0_rk * sin (0.5_rk * delta_theta)
            end do
         end if

c     E_varphi * dl
         do i = 0, m
            emf_varphi (i, j, k) = 0.5_rk *
     &           (emf_face (i  , j  , k  , 2, 1) +
     &            emf_face (i  , j+1, k  , 2, 1) + 
     &            emf_face (i  , j  , k  , 2, 2) +
     &            emf_face (i+1, j  , k  , 2, 2)) *
     &           r_if (i) * sinus_theta_if (j) * delta_varphi
         end do

         if (j .eq. n_s) then
            do i = 0, m
               emf_varphi (i, j-1, k) = 0.5_rk *
     &              (emf_face (i  , j-1, k, 2, 1) +
     &               emf_face (i  , j  , k, 2, 1) + 
     &               emf_face (i  , j-1, k, 2, 2) +
     &               emf_face (i+1, j-1, k, 2, 2)) *
     &              r_if (i) * sinus_theta_if (j-1) * delta_varphi
            end do
         end if
      end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk,l)
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         jp1 = j+1
         jm1 = j-1
         if (o .ne. 1) then
            kp1 = k+1
            km1 = k-1
         else
            kp1 = 1
            km1 = 1
         endif

         do i = 1, m

            im1 = i-1
            ip1 = i+1
            
            f_diff_b (i, j, k, 1) = 0.5_rk * (
     &           (emf_theta  (i  , j  , km1) +
     &            emf_varphi (i  , j  , k  ) - 
     &            emf_theta  (i  , j  , k  ) -
     &            emf_varphi (i  , jm1, k  )) +
     &           (emf_theta  (im1, j  , km1) +
     &            emf_varphi (im1, j  , k  ) - 
     &            emf_theta  (im1, j  , k  ) -
     &            emf_varphi (im1, jm1, k  ))) *
     &           (rr_inv (i) * delta_theta_mod_inv (j) * 
     &           delta_varphi_inv)

c     Note: emf_theta terms cancel in 2D
            f_diff_b (i, j, k, 2) = 0.5_rk * (
     &           (emf_r      (i  , j  , k  ) -
     &            emf_varphi (i  , j  , k  ) - 
     &            emf_r      (i  , j  , km1) +
     &            emf_varphi (im1, j  , k  )) +
     &           (emf_r      (i  , jm1, k  ) -
     &            emf_varphi (i  , jm1, k  ) - 
     &            emf_r      (i  , jm1, km1) +
     &            emf_varphi (im1, jm1, k  ))) *
     &           2.0_rk / (rr_if (i) - rr_if (i-1)) *
     &           sinus_theta_inv (j) * delta_varphi_inv

c     Note: emf_r terms cancel in 2D
            f_diff_b (i, j, k, 3) = 0.5_rk * (
     &           (emf_r      (i  , jm1, k  ) +
     &            emf_theta  (i  , j  , k  ) - 
     &            emf_r      (i  , j  , k  ) -
     &            emf_theta  (im1, j  , k  )) +
     &           (emf_r      (i  , jm1, km1) +
     &            emf_theta  (i  , j  , km1) - 
     &            emf_r      (i  , j  , km1) -
     &            emf_theta  (im1, j  , km1))) *
     &           2.0_rk / (rr_if (i) - rr_if (i-1)) *
     &           delta_theta_inv            

         end do
      end do

      return

      END SUBROUTINE calculate_db_dt

c     ==================================================================



c     ==================================================================

      SUBROUTINE emf_centre

c     ==================================================================

      USE precision

      USE size_cfc

      USE hydro_primitives_cfc, ONLY: b_1, b_2, b_3, v_1, v_2, v_3,
     &     e_1, e_2, e_3

      IMPLICIT NONE

      INTEGER (KIND=ik) :: i, j, k, jk, ji, ki

c     We need to set emf_face (:,:,: 1:2,3), i.e. the EMF on
c     varphi-interfaces by hand, since no fluxes in the
c     varphi-direction are computed.

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk,l)
#endif
      DO jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         DO i = 1, m
            e_1 (i, j, k) =
     &           v_2 (i, j, k) * b_3 (i, j, k) - 
     &           v_3 (i, j, k) * b_2 (i, j, k)
            e_2 (i, j, k) =
     &           v_3 (i, j, k) * b_1 (i, j, k) - 
     &           v_1 (i, j, k) * b_3 (i, j, k)
            e_3 (i, j, k) =
     &           v_1 (i, j, k) * b_2 (i, j, k) - 
     &           v_2 (i, j, k) * b_1 (i, j, k)
         END DO

         IF (j .LE. n_s + 3) THEN
            ji = n_s - 1 - (j - n_s)
            DO i = 1, m
               e_1 (i, ji, k) =
     &              v_2 (i, ji, k) * b_3 (i, ji, k) - 
     &              v_3 (i, ji, k) * b_2 (i, ji, k)
               e_2 (i, ji, k) =
     &              v_3 (i, ji, k) * b_1 (i, ji, k) - 
     &              v_1 (i, ji, k) * b_3 (i, ji, k)
               e_3 (i, ji, k) =
     &              v_1 (i, ji, k) * b_2 (i, ji, k) - 
     &              v_2 (i, ji, k) * b_1 (i, ji, k)
            END DO
         END IF

         IF (j .GE. n_e - 3) THEN
            ji = n_e + 1 - (j - n_e)
            DO i = 1, m
               e_1 (i, ji, k) =
     &              v_2 (i, ji, k) * b_3 (i, ji, k) - 
     &              v_3 (i, ji, k) * b_2 (i, ji, k)
               e_2 (i, ji, k) =
     &              v_3 (i, ji, k) * b_1 (i, ji, k) - 
     &              v_1 (i, ji, k) * b_3 (i, ji, k)
               e_3 (i, ji, k) =
     &              v_1 (i, ji, k) * b_2 (i, ji, k) - 
     &              v_2 (i, ji, k) * b_1 (i, ji, k)
            END DO
         END IF

         IF (o .GT. 1) THEN
            IF (k .LE. o_s + 3) THEN
               ki = o_s - 1 - (k - o_s)
               DO i = 1, m
                  e_1 (i, j, ki) =
     &                 v_2 (i, j, ki) * b_3 (i, j, ki) - 
     &                 v_3 (i, j, ki) * b_2 (i, j, ki)
                  e_2 (i, j, ki) =
     &                 v_3 (i, j, ki) * b_1 (i, j, ki) - 
     &                 v_1 (i, j, ki) * b_3 (i, j, ki)
                  e_3 (i, j, ki) =
     &                 v_1 (i, j, ki) * b_2 (i, j, ki) - 
     &                 v_2 (i, j, ki) * b_1 (i, j, ki)
               END DO
            END IF
            
            IF (k .GE. o_e - 3) THEN
               ki = o_e + 1 - (k - o_e)
               DO i = 1, m
                  e_1 (i, j, ki) =
     &                 v_2 (i, j, ki) * b_3 (i, j, ki) - 
     &                 v_3 (i, j, ki) * b_2 (i, j, ki)
                  e_2 (i, j, ki) =
     &                 v_3 (i, j, ki) * b_1 (i, j, ki) - 
     &                 v_1 (i, j, ki) * b_3 (i, j, ki)
                  e_3 (i, j, ki) =
     &                 v_1 (i, j, ki) * b_2 (i, j, ki) - 
     &                 v_2 (i, j, ki) * b_1 (i, j, ki)
               END DO
            END IF
         END IF
         
      END DO

      RETURN

      END SUBROUTINE emf_centre

c     ==================================================================



c     ==================================================================

      SUBROUTINE setup_b_field (b_0_pol, b_0_tor, mode)

c     ==================================================================

      USE precision

      USE phycon, ONLY: pc_bgeo

      USE size_cfc
      USE parameters_cfc

      USE conserved_cfc, ONLY: b_cap_1_hat, b_cap_2_hat, b_cap_3_hat
      USE hydro_primitives_cfc, ONLY: b_1, b_2, b_3, psi
      USE grid_cfc
      USE fluxes_cfc
      USE metric_cfc
      USE perm_aux_cfc

      USE configure

      IMPLICIT NONE

      REAL (KIND=rk), INTENT (IN) :: b_0_pol, b_0_tor
      REAL (KIND=rk), PARAMETER :: r_0 = 1e8
      INTEGER (KIND=ik), INTENT (IN) :: mode
      INTEGER (KIND=ik) :: i, j, k, jk, ii1, ii2

      SELECT CASE (mode)
         
      CASE (1)                  !constant field
         DO k = o_s-1, o_e+1
            DO j = n_s-1, n_e+1
               DO i = 0, m + 1
                  
c     A_r
                  emf_face (i,j,k,1,2) = 0.0_rk
                  emf_face (i,j,k,1,3) = 0.0_rk
c     A_theta
                  emf_face (i,j,k,1,1) = 0.0_rk
                  emf_face (i,j,k,2,3) = 0.0_rk
c     A_varphi
                  emf_face (i,j,k,2,1) =
     &                 r_if (i) * sinus_theta (min(j,n+1-j))
                  emf_face (i,j,k,2,2) =
     &                 r (i) * sinus_theta_if (min(j,n-j))
                  
               END DO
            END DO
         END DO

      CASE (2) ! as in Obergaulinger et al. (2018), Eqs. (20,21)

         DO k = o_s-1, o_e+1
            DO j = n_s-1, n_e+1
               DO i = 0, m + 1
c     A_r
                  emf_face (i,j,k,1,2) = r (i) * !0.5_rk *
     &                 r_0 ** 3 / (r_0 ** 3 + r (i) ** 3) *
c     &                 sinus_theta_if (min(j,n-j))
     &                 cosin_theta_if (min(j,n-j))
                  emf_face (i,j,k,1,3) = r (i) * !0.5_rk *
     &                 r_0 ** 3 / (r_0 ** 3 + r (i) ** 3) *
c     &                 sinus_theta (min(j,n+1-j))
     &                 cosin_theta (min(j,n+1-j))
c     only if A_r~cos \theta
                  if (j .eq. n/2) then
                     emf_face (i,j,k,1,2) = 0.0_rk
                  end if
                  if (j. gt. n/2) then
                     emf_face (i,j,k,1,2) = - emf_face (i,j,k,1,2)
                     emf_face (i,j,k,1,3) = - emf_face (i,j,k,1,3)
                  end if
c     A_theta
                  emf_face (i,j,k,1,1) = 0.0_rk
                  emf_face (i,j,k,2,3) = 0.0_rk
c     A_varphi
                  emf_face (i,j,k,2,1) = 0.5_rk  * r_if (i) *
     &                 r_0 ** 3 / (r_0 ** 3 + r_if (i) ** 3) *
     &                 sinus_theta (min(j,n+1-j))
                  emf_face (i,j,k,2,2) = 0.5_rk  * r (i) *
     &                 r_0 ** 3 / (r_0 ** 3 + r (i) ** 3)  *
     &                 sinus_theta_if (min(j,n-j))

               END DO
            END DO
         END DO

      CASE (3) ! from progenitor

c         b_cap_1_hat (:,:,:) = b_cap_1_hat (:,:,:) * 1e1_rk
c         b_cap_3_hat (:,:,:) = b_cap_3_hat (:,:,:) * 1e1_rk
         
         DO k = o_s-1, o_e+1
            DO j = n_s-1, n_e+1
               DO i = 0, m + 1                  
                  ii1 = min(max(i  ,1),m)
                  ii2 = min(max(i+1,1),m)
c     A_r
                  emf_face (i,j,k,1,2) = r (i) *
     &                 b_cap_3_hat(ii1,n_s,o_s) *
     &                 cosin_theta_if (min(j,n-j))
                  emf_face (i,j,k,1,3) = r (i) * !0.5_rk *
     &                 b_cap_3_hat(ii1,n_s,o_s) *
     &                 cosin_theta (min(j,n+1-j))
c     only if A_r~cos \theta
                  if (j .eq. n/2) then
                     emf_face (i,j,k,1,2) = 0.0_rk
                  end if
                  if (j. gt. n/2) then
                     emf_face (i,j,k,1,2) = - emf_face (i,j,k,1,2)
                     emf_face (i,j,k,1,3) = - emf_face (i,j,k,1,3)
                  end if
                  
c     A_theta
                  emf_face (i,j,k,1,1) = 0.0_rk
                  emf_face (i,j,k,2,3) = 0.0_rk
c     A_varphi
                  emf_face (i,j,k,2,1) = r_if (i) * 0.25_rk *
     &                 (b_cap_1_hat(ii1,n_s,o_s)+
     &                 b_cap_1_hat(ii2,n_s,o_s)) *
     &                 sinus_theta (min(j,n+1-j))
                  emf_face (i,j,k,2,2) = r (i) * 0.5_rk  *
     &                 b_cap_1_hat(ii1,n_s,o_s) *
     &                 sinus_theta_if (min(j,n-j))

               END DO
            END DO
         END DO
         
      CASE DEFAULT
         PRINT *, 'Not imposing any magnetic field.'
         emf_face (:,:,:,:,:) = 0.0_rk
         
      END SELECT

c     Scale vector potential and Convert from cgs units to geometrized units
      emf_face (:,:,:,1,2) =
     &     b_0_tor * pc_bgeo * emf_face (:,:,:,1,2)
      emf_face (:,:,:,1,3) =
     &     b_0_tor * pc_bgeo * emf_face (:,:,:,1,3)
      
      emf_face (:,:,:,2,1) = 
     &     b_0_pol * pc_bgeo * emf_face (:,:,:,2,1)
      emf_face (:,:,:,2,2) = 
     &     b_0_pol * pc_bgeo * emf_face (:,:,:,2,2)

      emf_face (:,:,:,1,1) = 0.0_rk
      emf_face (:,:,:,2,3) = 0.0_rk
         

      CALL calculate_db_dt (0)

c     Need extra factor of 1/2 because of the way we assemble
c     the E field from the interface fluxes
      b_cap_1_hat (:,:,:) = 0.5_rk * f_diff_b (:,:,:,1)
      b_cap_2_hat (:,:,:) = 0.5_rk * f_diff_b (:,:,:,2)
      b_cap_3_hat (:,:,:) = 0.5_rk * f_diff_b (:,:,:,3)
      b_1 (1:m,n_s:n_e,o_s:o_e) =
     &     0.5_rk * f_diff_b (1:m,n_s:n_e,o_s:o_e,1)
      b_2 (1:m,n_s:n_e,o_s:o_e) =
     &     0.5_rk * f_diff_b (1:m,n_s:n_e,o_s:o_e,2)
      b_3 (1:m,n_s:n_e,o_s:o_e) =
     &     0.5_rk * f_diff_b (1:m,n_s:n_e,o_s:o_e,3)

      psi (1:m,n_s:n_e,o_s:o_e) = 0.0_rk

      RETURN
      
      END SUBROUTINE setup_b_field

c     ==================================================================


#endif /* CFC_MHD */
      
      END MODULE induction

c     ==================================================================
