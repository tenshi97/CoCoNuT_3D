c     ==================================================================

      MODULE gr_rescale

c     ==================================================================

      IMPLICIT NONE

      CONTAINS



c     ==================================================================

      SUBROUTINE rescale_velocities (flag)

c     ==================================================================

      USE precision

      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
      USE hydro_primitives_cfc
      USE conserved_cfc
      USE interface_hydro_cfc
      USE perm_aux_cfc

      IMPLICIT NONE

      integer (kind=ik), intent (in) :: flag
      integer (kind=ik) :: i, j, k, jk, n_zones, o_zones,
     &     n_start, o_start


ccc   1: change v_2, v_3 to v_theta, v_phi
ccc   2: change v_theta, v_phi to v_2, v_3
ccc   3: change v_theta, v_phi to v_2, v_3 at r interfaces
ccc   4: change v_theta, v_phi to v_2, v_3 at theta interfaces
ccc   5: change v_theta, v_phi to v_2, v_3 at phi interfaces
      
ccc   Caution: Also have to rescale S_i in 1, 2 for consistency with boundary conditions


      if (flag .eq. 1) then

         if (config%nsdim .ge. 2) then
            n_zones = n_loc! + 8
            n_start = n_s! - 4
         else
            n_zones = 1
            n_start = 1
         endif

         if (config%nsdim .eq. 3) then
            o_zones = o_loc! + 8
            o_start = o_s! - 4
         else
            o_zones = 1
            o_start = 1
         endif

#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_zones * o_zones
         
            k = int((jk + n_zones - 1) / n_zones )
            j = (n_start - 1) + (jk - (k - 1) * n_zones)
            k = k + o_start - 1
         
            do i = - 3, m + 4

               v_2 (i, j, k) = v_2 (i, j, k) *
     &              r_inv (i)
               v_3 (i, j, k) = v_3 (i, j, k) *
     &              r_inv (i) *  sinus_theta_inv (j)
c     &              (r_inv (i) *  sinus_theta_inv (j)) ** 2

            enddo
         enddo

      elseif (flag .eq. 2) then

         if (config%nsdim .ge. 2) then
            n_zones = n_loc  + 8
            n_start = n_s - 4
         else
            n_zones = 1
            n_start = 1
         endif

         if (config%nsdim .eq. 3) then
            o_zones = o_loc + 8
            o_start = o_s - 4
         else
            o_zones = 1
            o_start = 1
         endif

#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_zones * o_zones
         
            k = int((jk + n_zones - 1) / n_zones )
            j = (n_start - 1) + (jk - (k - 1) * n_zones)
            k = k + o_start - 1

            if ((j .ge. n_s .and. j .le. n_e) .and.
     &           (k .ge. o_s .and. k .le. o_e)) then

               do i = - 3, m + 4
                  v_2 (i, j, k) = v_2 (i, j, k) *
     &                 r (i)
                  v_3 (i, j, k) = v_3 (i, j, k) *
     &                 r (i) * sinus_theta (j)
c     &                 (r (i) * sinus_theta (j)) ** 2
               enddo

            else if ((j .ge. n_s .and. j .le. n_e) .or.
     &              (k .ge. o_s .and. k .le. o_e)) then
               
               do i = 1, m
                  v_2 (i, j, k) = v_2 (i, j, k) *
     &                 r (i)
                  v_3 (i, j, k) = v_3 (i, j, k) *
     &                 r (i) * sinus_theta (j)
c     &                 (r (i) * sinus_theta (j)) ** 2
               enddo
            end if
               
         enddo

      elseif (flag .eq. 3) then

         n_zones = n_loc
         n_start = n_s

         o_zones = o_loc
         o_start = o_s

#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_zones * o_zones
         
            k = int((jk + n_zones - 1) / n_zones )
            j = (n_start - 1) + (jk - (k - 1) * n_zones)
            k = k + o_start - 1

            do i = 0, m + 1

               v_2_if (i, j, k, 0) = v_2_if (i, j, k, 0) *
     &              r_if (i - 1)
               v_2_if (i, j, k, 1) = v_2_if (i, j, k, 1) *
     &              r_if (i)
               v_3_if (i, j, k, 0) = v_3_if (i, j, k, 0) *
     &              r_if (i - 1) * sinus_theta (j)
c     &              (r_if (i - 1) * sinus_theta (j)) ** 2
               v_3_if (i, j, k, 1) = v_3_if (i, j, k, 1) *
     &                 r_if (i) * sinus_theta (j)
c     &                 (r_if (i) * sinus_theta (j)) ** 2

            enddo
         enddo

      elseif (flag .eq. 4) then

         n_zones = n_loc + 2
         n_start = n_s - 1

         o_zones = o_loc
         o_start = o_s

#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_zones * o_zones
            
            k = int((jk + n_zones - 1) / n_zones )
            j = (n_start - 1) + (jk - (k - 1) * n_zones)
            k = k + o_start - 1

            do i = 1, m

               v_2_if (i, j, k, 0) = v_2_if (i, j, k, 0) *
     &              r (i)
               v_2_if (i, j, k, 1) = v_2_if (i, j, k, 1) *
     &              r (i)
               v_3_if (i, j, k, 0) = v_3_if (i, j, k, 0) *
     &              r (i) * sinus_theta_if (j - 1)
c     &              (r (i) * sinus_theta_if (j - 1)) ** 2
               v_3_if (i, j, k, 1) = v_3_if (i, j, k, 1) *
     &              r (i) * sinus_theta_if (j)
c     &              (r (i) * sinus_theta_if (j)) ** 2

            enddo
         enddo

      elseif (flag .eq. 5) then

         n_zones = n_loc
         n_start = n_s

         o_zones = o_loc + 2
         o_start = o_s - 1

#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
         do jk = 1, n_zones * o_zones
            
            k = int((jk + n_zones - 1) / n_zones )
            j = (n_start - 1) + (jk - (k - 1) * n_zones)
            k = k + o_start - 1

            do i = 1, m

               v_2_if (i, j, k, 0) = v_2_if (i, j, k, 0) *
     &              r (i)
               v_2_if (i, j, k, 1) = v_2_if (i, j, k, 1) *
     &              r (i)
               v_3_if (i, j, k, 0) = v_3_if (i, j, k, 0) *
     &              r (i) * sinus_theta (j)
c     &              (r (i) * sinus_theta (j)) ** 2
               v_3_if (i, j, k, 1) = v_3_if (i, j, k, 1) *
     &              r (i) * sinus_theta (j)
c     &              (r (i) * sinus_theta (j)) ** 2

            enddo
         enddo

      endif

      return

      END SUBROUTINE rescale_velocities

c     ==================================================================

      END MODULE gr_rescale

c     ==================================================================
