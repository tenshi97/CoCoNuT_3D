c     ==================================================================

      MODULE renormalize

c     ==================================================================

      IMPLICIT NONE

      CONTAINS

c     ==================================================================

      SUBROUTINE renormalize_mass_fractions

c     ==================================================================

      USE precision

      USE size_cfc
      USE conserved_cfc, ONLY: d_cap_hat, d_cap_xnu_hat
      USE hydro_primitives_cfc, ONLY: xnnu

      USE configure
      USE abort

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk, i_n
      real (kind=rk) :: xnu_sum (1 : m)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk,i_n,xnu_sum)
#endif
      do jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         
         xnu_sum (:) = 0.0_rk
           
         do i_n = 1, config%qn-1
            do i = 1, m
               d_cap_xnu_hat (i, j, k, i_n) = 
     &              MAX(d_cap_xnu_hat (i, j, k, i_n), 0.0_rk)
               xnu_sum(i) = xnu_sum (i) + d_cap_xnu_hat (i, j, k, i_n) 
            end do
         end do

         DO i = 1, m
            xnu_sum(i) = d_cap_hat (i, j, k) / xnu_sum(i)
         END DO

         do i_n = 1, config%qn-1
            do i = 1, m
               d_cap_xnu_hat (i, j, k, i_n) =
     &              d_cap_xnu_hat (i, j, k, i_n) * xnu_sum (i)
               xnnu (i, j, k, i_n) =
     &              d_cap_xnu_hat (i, j, k, i_n)/
     &              d_cap_hat     (i, j, k)
            end do
         enddo
         
      end do

      return

      END SUBROUTINE renormalize_mass_fractions

c     ==================================================================



      END MODULE renormalize

c     ==================================================================
