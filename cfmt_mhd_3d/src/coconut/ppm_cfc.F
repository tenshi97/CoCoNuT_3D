#define PPM_FIX_ORIGIN

c     ==================================================================

      MODULE gr_ppm

c     ==================================================================

      IMPLICIT NONE

      CONTAINS





c     ******************************************************************
c     ************************** r-direction ***************************
c     ******************************************************************

c     ==================================================================

      SUBROUTINE coefficients_r

c     ==================================================================
c
c     subroutine to calculate the coefficients for interpolation in
c     r direction
c
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE grid_cfc
      USE ppm_cfc

      IMPLICIT NONE
      
      real (kind=rk) :: tmp_1 (- 1 : m + 3), tmp_2 (- 1 : m + 3), 
     &     tmp_3 (- 1 : m + 3)
      real (kind=rk) :: tmp

      integer (kind=ik) :: i


#ifdef SEKORA_PPM
      c_1_r (:) = 0.5_rk
      c_2_r (:) = 0.5_rk
      c_3_r (:) = - 1.0_rk / 6.0_rk
      c_4_r (:) = 1.0_rk / 6.0_rk
      c_5_r (:) = 0.5_rk
#else /* SEKORA_PPM */
c     As in original Colella & Woodward paper
      do i = - 1, m + 3
ccc   calculate some delta_r for the coefficients for delta u (1.7)
         tmp_1 (i) = delta_r (i - 1) + delta_r (i)
         tmp_2 (i) = tmp_1 (i) + delta_r (i)
         tmp_3 (i) = tmp_1 (i) + delta_r (i - 1)
      enddo
      
      do i = - 1, m + 2
ccc   calculate the coefficients for delta u (1.7)
         tmp = delta_r (i) /
     &        (delta_r (i - 1) + delta_r (i) + delta_r (i + 1))
         c_1_r (i) = tmp * tmp_3 (i) / tmp_1 (i + 1)
         c_2_r (i) = tmp * tmp_2 (i + 1) / tmp_1 (i)
      enddo
      
      do i = - 1, m + 1
ccc   calculate the coefficients before delta u_i+1, delta u_i and u_i+1 - u_i (1.6)
         tmp = 1.0_rk / (tmp_1 (i) + tmp_1 (i + 2))
         c_3_r (i) = - tmp * delta_r (i) * tmp_1 (i) / tmp_3 (i + 1)
         c_4_r (i) = tmp * delta_r (i + 1) * tmp_1 (i + 2) / 
     &        tmp_2 (i + 1)
         c_5_r (i) = (delta_r (i) -
     &        2.0_rk * (delta_r (i + 1) * c_3_r (i) + 
     &        delta_r (i) * c_4_r (i))) / tmp_1 (i + 1)
      enddo
#endif /* SEKORA_PPM */      

      return
      
      END SUBROUTINE coefficients_r

c     ==================================================================




c     ==================================================================

      SUBROUTINE interpolate_r (q, q_if, is_mass_fraction, is_positive)

c     ==================================================================

c     subroutine to interpolate in r direction

c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE perm_aux_cfc
      USE ppm_cfc

      IMPLICIT NONE
      
      real (kind=rk), intent (in) ::
     &     q (- 3 : m + 4, n_s - 4 : n_e + 4, o_s - 4 : o_e + 4)
      real (kind=rk), intent (out) :: q_if (- 1 : m + 2,
     &     n_s - 2 : n_e + 2, o_s - 2 : o_e + 2, 0 : 1) 

      real (kind=rk) :: delta_q_sign
#ifdef SEKORA_PPM
      real (kind=rk), parameter :: c_ppm = 1.25_rk
      real (kind=rk) :: s2c, d2lim, d2p, d2m, d2c
#endif      
      logical, intent(in), optional :: is_mass_fraction, is_positive
      logical :: mass_fraction=.false.,positive=.false.
      
      integer (kind=ik) :: i, j, k, jk

      if (present(is_positive)) then
         positive=.true.
      else
         positive=.false.
      end if
      if (present(is_mass_fraction)) then
         positive=.true.
         mass_fraction=.true.
      else
         mass_fraction=.false.
      end if

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

c#ifndef SEKORA_PPM         
         do i = - 2, m + 2
ccc   a_i - a_i-1 (for delta a_i)
            tmp_ppm (i, j, k) = q (i + 1, j, k) - q (i, j, k)
         enddo

         do i = - 1, m + 2
ccc   delta a_i (C & W)
            delta_q (i, j, k) = c_1_r (i) * tmp_ppm (i, j, k) +
     &           c_2_r (i) * tmp_ppm (i - 1, j, k)
ccc   convert delta a_i to delta_m a_i
#ifdef PPM_FIX_ORIGIN
         enddo
         do i = 2, m + 2
#endif
            if (tmp_ppm (i, j, k) * tmp_ppm (i - 1, j, k) .gt.
     &           0.0_rk) then
               delta_q_sign =
     &              delta_q (i, j, k) / abs (delta_q (i, j, k))
               delta_q (i, j, k) =
     &              min (abs (delta_q (i, j, k)), 
     &              2.0_rk * min (abs (tmp_ppm (i - 1, j, k)), 
     &              abs (tmp_ppm (i, j, k)))) * delta_q_sign
            else
               delta_q (i, j, k) = 0.0_rk
            endif
         enddo

         do i = - 1, m + 1
ccc   a_i+1/2
            q_if (i, j, k, 1) = q (i, j, k) +
     &           c_5_r (i) * tmp_ppm (i, j, k) + 
     &              c_3_r (i) * delta_q (i + 1, j, k) +
     &           c_4_r (i) * delta_q (i, j, k)
            q_if (i + 1, j, k, 0) = q_if (i, j, k, 1)
         enddo
c#endif /* SEKORA_PPM */         
         
#ifdef SEKORA_PPM
         do i = - 1, m + 1
            q_if (i, j, k, 1) =
     &           (37.0_rk * (q (i, j, k) + q (i + 1, j, k)) -
     &           8.0_rk * (q (i - 1, j, k) + q (i + 2, j, k)) +
     &           (q (i - 2, j, k) + q (i + 3, j, k))) / 60.0_rk
     
            if ((q_if (i, j, k, 1) - q (i, j, k)) *
     &         (q (i + 1, j, k) - q_if (i, j, k, 1)) .lt. 0.0_rk) then
               d2c = 3.0_rk * (q (i, j, k) + q (i + 1, j, k)
     &              - 2.0_rk * q_if (i, j, k, 1))
               d2p = q (i    , j ,k) + q (i + 2, j, k)
     &              - 2.0_rk * q (i + 1, j, k)
               d2m = q (i - 1, j ,k) + q (i + 1, j, k)
     &              - 2.0_rk * q (i    , j, k)
               s2c = sign (1.0_rk, d2c)
               d2lim = max (min (c_ppm * s2c * d2m, c_ppm * s2c * d2p,
     &              s2c * d2c), 0.0_rk)
               q_if (i, j, k, 1) =
     &              0.5_rk * (q (i, j, k) + q (i + 1, j, k)) -
     &              s2c * d2lim / 6.0_rk !correct sign of correction term?
            end if
            q_if (i + 1, j, k, 0) = q_if (i, j, k, 1)
         end do
         if (positive) then
            do i = - 1, m + 1
               if (q_if (i, j, k, 1) .lt. 0.0_rk) then
                  q_if (i, j, k, 1) = q (i, j, k) +
     &                 c_5_r (i) * tmp_ppm (i, j, k) + 
     &                 c_3_r (i) * delta_q (i + 1, j, k) +
     &                 c_4_r (i) * delta_q (i, j, k)
                  q_if (i + 1, j, k, 0) = q_if (i, j, k, 1)
               end if
            end do
         end if
         if (mass_fraction) then
            do i = - 1, m + 1
               if (q_if (i, j, k, 1) .gt. 1.0_rk) then
                  q_if (i, j, k, 1) = q (i, j, k) +
     &                 c_5_r (i) * tmp_ppm (i, j, k) + 
     &                 c_3_r (i) * delta_q (i + 1, j, k) +
     &                 c_4_r (i) * delta_q (i, j, k)
                  q_if (i + 1, j, k, 0) = q_if (i, j, k, 1)
               end if
            end do
         end if
#endif /* SEKORA_PPM */      
      enddo
      
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END DO
#endif

      return

      END SUBROUTINE interpolate_r

c     ==================================================================




c     ==================================================================

      SUBROUTINE cd_steepening_r (xq, x_if)

c     ==================================================================
c
c     subroutine to detect contact discontinuity and sharpen
c     rho profile there
c
c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE grid_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE parameters_cfc
      USE ppm_cfc
      USE shock_cfc, ONLY: lshock

      IMPLICIT NONE

      real (kind=rk), intent (in) ::
     &     xq (- 3 : m + 4, n_s - 4 : n_e + 4, o_s - 4 : o_e + 4)
      real (kind=rk), intent (inout) :: x_if (- 1 : m + 2, 
     &     n_s - 2 : n_e + 2, o_s - 2 : o_e  + 2, 0 : 1)

      real (kind=rk) :: delta_2_rho_minus, delta_2_rho_plus
      real (kind=rk) :: eta, eta_tilde

      real (kind=rk) :: criterion_1, criterion_2,
     &     criterion_3, criterion_4, tmp
      
      real (kind=rk) :: gamma_eos_tmp
      
      integer (kind=ik) :: i, j, k, jk

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 0, m + 1

ccc   C & Q delta^2 a_i
            delta_2_rho_minus =
     &           ((xq (i, j, k) - xq (i - 1, j, k)) / 
     &           (delta_r (i - 1) + delta_r (i)) - 
     &           (xq (i - 1, j, k) - xq (i - 2, j, k)) / 
     &           (delta_r (i - 2) + delta_r (i - 1))) / 
     &           (delta_r (i - 2) + delta_r (i - 1) + delta_r (i))
            delta_2_rho_plus =
     &           ((xq (i + 2, j, k) - xq (i + 1, j, k)) / 
     &           (delta_r (i + 1) + delta_r (i + 2)) - 
     &           (xq (i + 1, j, k) - xq (i, j, k)) / 
     &           (delta_r (i) + delta_r (i + 1))) / 
     &           (delta_r (i) + delta_r (i + 1) + delta_r (i + 2))

            tmp = xq (i + 1, j, k) - xq (i - 1, j, k)
            if (tmp .eq. 0.0_rk) tmp = 1.0e-10_rk
ccc   check if denominator of second factor in (1.17 top) is 0
               
            criterion_1 = delta_2_rho_plus * delta_2_rho_minus
            criterion_2 = abs (tmp) - epsilon_1 * 
     &           min (abs (xq (i + 1, j, k)),
     &           abs (xq (i - 1, j, k)))
            
            if (criterion_1 .gt. 0.0_rk .or.
     &           criterion_2 .lt. 0.0_rk) then
ccc   second case in (1.17 top)
               eta_tilde = 0.0_rk
            else
ccc   first case in (1.17 top) (minus sign integrated in first term)
               eta_tilde = 
     &              (delta_2_rho_minus - delta_2_rho_plus) * 
     &              (delta_r_if (i - 1) ** 3 + delta_r_if (i) ** 3) / 
     &              ((delta_r_if (i) + delta_r_if (i - 1)) * tmp)
            endif
            
            eta = max (0.0_rk, min (eta_1 * (eta_tilde -
     &              eta_2), 1.0_rk))

ccc   (1.16)
#ifdef CFC_TRANSPORT2
            gamma_eos_tmp = c_sound_squared (i, j, k) *
     &           rho (i, j, k) * h (i, j, k) / p (i, j, k)
#else /* CFC_TRANSPORT2 */
            gamma_eos_tmp = c_sound_squared (i, j, k) *
     &           rho (i, j, k) / p (i, j, k)
#endif /* CFC_TRANSPORT2 */
            criterion_3 = gamma_eos_tmp * k_cap_0 * 
     &           abs (rho (i + 1, j, k) - rho (i - 1, j, k)) /
     &           min (rho (i + 1, j, k), rho (i - 1, j, k))
            criterion_4 =
     &           abs (p (i + 1, j, k) - p (i - 1, j, k)) / 
     &           min (p (i + 1, j, k), p (i - 1, j, k))
            
            if (criterion_3 .lt. criterion_4 .or.
     &           lshock (i, j, k)) eta = 0.0_rk
ccc   gamma * K_0 * density jump < pressure jump --> eta = 0
            x_if (i, j, k, 0) = x_if (i, j, k, 0) *
     &           (1.0_rk - eta) + (xq (i - 1, j, k) +
     &           0.5_rk * delta_q (i - 1, j, k)) * eta
            x_if (i, j, k, 1) = x_if (i, j, k, 1) *
     &           (1.0_rk - eta) + (xq (i + 1, j, k) -
     &           0.5_rk * delta_q (i + 1, j, k)) * eta
ccc   (1.15)
         enddo
      enddo
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END DO
#endif

      return
      
      END SUBROUTINE cd_steepening_r

c     ==================================================================



c     ==================================================================

      SUBROUTINE flatten_r

c     ==================================================================
c
c     subroutine to flatten profile near strong shocks
c
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE parameters_cfc
      USE ppm_cfc


      use configure

      IMPLICIT NONE

      real (kind=rk) :: q, tmp

      real (kind=rk) :: f, fx (0:m+1)

      real (kind=rk) :: delta_v, delta_p_2
      real (kind=rk) :: criterion

      integer (kind=ik) :: i, j, k, jk, i_nuc
      

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = - 1, m + 2
            delta_p (i, j, k) = p (i + 1, j, k) - p (i - 1, j, k)
            delta_v = v_1 (i + 1, j, k) - v_1 (i - 1, j, k)
ccc   can multiply with min (...) as it is positive ---> new criterion
            criterion = epsilon_2_r * min (p (i + 1, j, k), 
     &           p (i - 1, j, k)) - abs (delta_p (i, j, k))
            
            if (criterion .lt. 0.0_rk .and. delta_v .lt. 0.0_rk) then
ccc   set w_i to 1
               q = 1.0_rk
            else
               q = 0.0_rk
            endif
            
ccc   calculate f_i^tilde
            delta_p_2 = p (i + 2, j, k) - p (i - 2, j, k)
            if (delta_p_2 .eq. 0.0_rk) then
ccc   check if denominator is zero
               if (delta_p (i, j, k) .eq. 0.0_rk) then
ccc   this is somewhat arbitrary but obviously works
                  tmp = - omega_1
               else
                  tmp = 1.0_rk - omega_1
               endif
            else
ccc   regular
               tmp = delta_p (i, j, k) / delta_p_2 - omega_1
            endif
            f_tilde (i, j, k) =
     &           min (1.0_rk, q * max (0.0_rk, tmp * 
     &           omega_2))
         enddo

         do i = 0, m + 1
            if (delta_p (i, j, k) .lt. 0.0_rk) then
ccc   check for negative delta pressure (check in which direction shock travels)
               f = max (f_tilde (i, j, k), f_tilde (i + 1, j, k))
            else
               f = max (f_tilde (i, j, k), f_tilde (i - 1, j, k))
            endif
            fx (i) = f
               
            rho_if (i, j, k, 0) = rho (i, j, k) * f +
     &           rho_if (i, j, k, 0) * (1.0_rk - f)
            rho_if (i, j, k, 1) = rho (i, j, k) * f +
     &           rho_if (i, j, k, 1) * (1.0_rk - f)
            
            v_1_if (i, j, k, 0) = v_1 (i, j, k) * f +
     &           v_1_if (i, j, k, 0) * (1.0_rk - f)
            v_1_if (i, j, k, 1) = v_1 (i, j, k) * f +
     &           v_1_if (i, j, k, 1) * (1.0_rk - f)
            
            v_2_if (i, j, k, 0) = v_2 (i, j, k) * f +
     &           v_2_if (i, j, k, 0) * (1.0_rk - f)
            v_2_if (i, j, k, 1) = v_2 (i, j, k) * f +
     &           v_2_if (i, j, k, 1) * (1.0_rk - f)
            
            v_3_if (i, j, k, 0) = v_3 (i, j, k) * f +
     &           v_3_if (i, j, k, 0) * (1.0_rk - f)
            v_3_if (i, j, k, 1) = v_3 (i, j, k) * f +
     &           v_3_if (i, j, k, 1) * (1.0_rk - f)
            
            eps_if (i, j, k, 0) = eps (i, j, k) * f +
     &           eps_if (i, j, k, 0) * (1.0_rk - f)
            eps_if (i, j, k, 1) = eps (i, j, k) * f +
     &           eps_if (i, j, k, 1) * (1.0_rk - f)
            
            p_if (i, j, k, 0) = p (i, j, k) * f +
     &           p_if (i, j, k, 0) * (1.0_rk - f)
            p_if (i, j, k, 1) = p (i, j, k) * f +
     &           p_if (i, j, k, 1) * (1.0_rk - f)
            
            gamm_if (i, j, k, 0) =
     &           gamm (i, j, k) * f +
     &           gamm_if (i, j, k, 0) * (1.0_rk - f)
            gamm_if (i, j, k, 1) =
     &           gamm (i, j, k) * f +
     &           gamm_if (i, j, k, 1) * (1.0_rk - f)
            
         enddo
         do i_nuc = 1, config%qn
            do i = 0, m + 1
               xnu_if (i, j, k, 0, i_nuc) =
     &              xnnu (i, j, k, i_nuc) * fx (i) +
     &              xnu_if (i, j, k, 0, i_nuc) * (1.0_rk - fx (i))
               xnu_if (i, j, k, 1, i_nuc) =
     &              xnnu (i, j, k, i_nuc) * fx (i) +
     &              xnu_if (i, j, k, 1, i_nuc) * (1.0_rk - fx (i))
            enddo
         enddo
      enddo
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END DO
#endif

      return
      
      END SUBROUTINE flatten_r
      
c     ==================================================================




c     ==================================================================

      SUBROUTINE monotonize_r (q, q_if)

c     ==================================================================
c
c     subroutine to monotonize reconstructed quantities
c      
c     ------------------------------------------------------------------
      
      USE precision
      USE size_cfc

      IMPLICIT NONE

      real (kind=rk), intent (in)  ::
     &     q (- 3 : m + 4, n_s - 4 : n_e + 4, o_s - 4 : o_e + 4)
      real (kind=rk), intent (inout) :: q_if (- 1 : m + 2, 
     &     n_s - 2 : n_e + 2, o_s - 2 : o_e + 2, 0 : 1) 

#ifdef SEKORA_PPM
      real (kind=rk) :: qi_old (- 1 : m + 2,
     &     n_s - 2 : n_e + 2, o_s - 2 : o_e + 2, 0 : 1)
      real (kind=rk), parameter :: c_ppm = 1.25_rk
      real (kind=rk) :: s2ppm, d2lim, d2p, d2m, d2c, d2ppm, alpl, almn,
     &     almx, scri, qav
      real (kind=rk) :: dlft, drgt, dqlft, dqrgt
      logical :: bigp, bigm, extremum
#endif      
      
      integer (kind=ik) :: i, j, k, jk

#ifdef SEKORA_PPM
      qi_old (:,:,:,:) = q_if (:,:,:,:)
#endif

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         
#ifdef SEKORA_PPM
         do i = 0, m + 1

            alpl = q_if (i, j, k, 1) - q (i, j, k)
            almn = q_if (i, j, k, 0) - q (i, j, k)
            bigp = abs(alpl) .gt. 2.0_rk * abs(almn)
            bigm = abs(almn) .gt. 2.0_rk * abs(alpl)
            qav = 1e-10_rk * MAX(ABS(q (i, j, k)),
     &           ABS(q (i - 1, j, k)), ABS(q (i + 1, j, k)))
            extremum = .false.
            if (bigp .or. bigm) then
c     Condition for extrema as in CASTRO
               dlft = q (i, j, k) - qi_old (i, j, k, 0)
               drgt = qi_old (i, j, k, 1) - q (i, j, k)
               dqlft = q (i   , j, k) - q (i - 1, j, k)
               dqrgt = q (i +1, j, k) - q (i    , j, k)
               if (min (abs(dlft), abs (drgt)) .ge.
     &              min (abs(dqlft), abs(dqrgt))) then
                  extremum = dlft * drgt .le. 0.0_rk
               else
                  extremum = dqlft * dqrgt .le. 0.0_rk
               end if
            end if
            if (alpl * almn .ge. 0.0_rk .or. !extremum) then
     &           (q (i, j, k) - q (i - 1, j, k)) *
     &           (q (i, j, k) - q (i + 1, j, k)) .ge. 0.0_rk) then
!     local extremum
               d2ppm = 6.0_rk * (alpl + almn)
c     Use cell averages for estimating 2nd derivative
               d2c = q (i+1,j,k) + q (i-1,j,k) - 2.0_rk * q (i,j,k)
               d2m = q (i-2,j,k) + q (i,j,k) - 2.0_rk * q (i-1,j,k)
               d2p = q (i,j,k) + q (i+2,j,k) - 2.0_rk * q (i+1,j,k)
               s2ppm = sign (1.0_rk, d2ppm)
               d2lim = max (min (s2ppm * d2ppm, c_ppm * s2ppm * d2m,
     &              c_ppm * s2ppm * d2p, c_ppm * s2ppm * d2c), 0.0_rk)
               q_if (i, j, k, 1) = q (i, j, k) +
     &              alpl * abs (d2lim) / (abs(d2ppm) + 1e-99_rk)
               q_if (i, j, k, 0) = q (i, j, k) +
     &              almn * abs (d2lim) / (abs(d2ppm) + 1e-99_rk)

            else if (bigp) then
               almx = - 0.25_rk * alpl ** 2 / (alpl + almn)
               scri = q (i - 1, j, k) - q (i, j, k)
               if (almn*almx .ge. almn*scri) then
                  if (almn*(scri-almn) .ge. qav) then
                     alpl = -2.0_rk * scri - sign (2.0_rk, almn) *
     &                    sqrt (max(scri ** 2 - almn * scri, 0.0_rk))
                  else
                     alpl = - 2.0_rk * almn
                  end if
               end if
               q_if (i, j, k, 1) = q (i, j, k) + alpl
            else if (bigm) then
               almx = - 0.25_rk * almn ** 2 / (alpl + almn)
               scri = q (i + 1, j, k) - q (i, j, k)
               if (alpl*almx .ge. alpl*scri) then
                  if (alpl*(scri-alpl) .ge. qav) then
                     almn = - 2.0_rk * scri - sign (2.0_rk, alpl) *
     &                    sqrt (max(scri ** 2 - alpl * scri, 0.0_rk))
                  else
                     almn = - 2.0_rk * alpl
                  end if
               end if
               q_if (i, j, k, 0) = q (i, j, k) + almn
            end if
            
#else /* SEKORA_PPM */
         do i = 0, m + 1
c     Standard PPM  (Colella & Woodward)
            if ((q_if (i, j, k, 1) - q (i, j, k)) * 
     &           (q_if (i, j, k, 0) - q (i, j, k)) .ge. 0.0_rk) then
ccc   first case of (1.10)
               q_if (i, j, k, 0) = q (i, j, k)
               q_if (i, j, k, 1) = q (i, j, k)
            else
               if ((q_if (i, j, k, 1) - q_if (i, j, k, 0)) *
     &              (q_if (i, j, k, 0) - 3.0_rk * q (i, j, k) +
     &              2.0_rk * q_if (i, j, k, 1)) .lt. 0.0_rk) then
ccc   second case of (1.10) -- after a little algebra, bring rhs to the left
                  q_if (i, j, k, 0) = 3.0_rk * q (i, j, k) -
     &                 2.0_rk * q_if (i, j, k, 1)
               endif
               if ((q_if (i, j, k, 1) - q_if (i, j, k, 0)) * 
     &              (- q_if (i, j, k, 1) + 3.0_rk * q (i, j, k) -
     &              2.0_rk * q_if (i, j, k, 0)) .lt. 0.0_rk) then
ccc   third case of (1.10) -- after a little algebra, bring rhs to the left
                  q_if (i, j, k, 1) = 3.0_rk * q (i, j, k) -
     &                 2.0_rk * q_if (i, j, k, 0)
               endif
            endif
#endif /* SEKORA_PPM */         
         enddo
      enddo
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END DO
#endif

      return

      END SUBROUTINE monotonize_r

c     ==================================================================



c     ==================================================================

      SUBROUTINE flatten_xnu_r

c     ==================================================================
c
c     subroutine to additionally flatten mass fraction profiles
c     near strong shocks, and enforce \sum_i X_i=1.0_rk
c
c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE parameters_cfc

      IMPLICIT NONE

      real (kind=rk) :: f, f_xnu (- 1 : m + 2)
      real (kind=rk) :: xnu_sum (0 : m + 1, 0 : 1)

      integer (kind=ik) :: i, j, k, jk, i_nuc

#ifndef SEKORA_PPM
      do i_nuc = 1, config%qn

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc

            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            do i = - 1, m + 2
               
               if ((xnnu(i+1,j,k,i_nuc)-xnnu(i,j,k,i_nuc)) *
     &              (xnnu(i,j,k,i_nuc)-xnnu(i-1,j,k,i_nuc))
     &              .lt. 0.0_rk) then
                  f_xnu (i) = 1.0_rk
               else
                  f_xnu (i) = 0.0_rk
               endif
               
            enddo

            do i = 0, m + 1
            
               f = max (f_xnu (i - 1), f_xnu (i + 1))
               f = 0.5_rk * max (2.0_rk * f_xnu (i), f)
               
               xnu_if (i, j, k, 0, i_nuc) =
     &              xnnu (i, j, k, i_nuc) * f +
     &              xnu_if (i, j, k, 0, i_nuc) * (1.0_rk - f)
               xnu_if (i, j, k, 1, i_nuc) =
     &              xnnu (i, j, k, i_nuc) * f +
     &                 xnu_if (i, j, k, 1, i_nuc) * (1.0_rk - f)
               
            enddo
         enddo
      enddo
#endif /* SEKORA_PPM */

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         xnu_sum(:,:)=0.0_rk
         do i_nuc = 1, config%qn-1
            do i= 0, m + 1
               xnu_sum (i, 0) =
     &              xnu_sum (i, 0) + xnu_if (i, j, k, 0, i_nuc)
               xnu_sum (i, 1) =
     &              xnu_sum (i, 1) + xnu_if (i, j, k, 1, i_nuc)
            end do
         end do
         do i_nuc = 1, config%qn-1
            do i= 0, m + 1
               xnu_if (i, j, k, 0, i_nuc) =
     &              xnu_if (i, j, k, 0, i_nuc) / xnu_sum (i, 0)
               xnu_if (i, j, k, 1, i_nuc) =
     &              xnu_if (i, j, k, 1, i_nuc) / xnu_sum (i, 1)
            end do
         end do
      end do

      return

      END SUBROUTINE flatten_xnu_r

c     ==================================================================




c     ******************************************************************
c     ************************ theta-direction *************************
c     ******************************************************************

c     ==================================================================

      SUBROUTINE coefficients_theta

c     ==================================================================
c
c     subroutine to calculate the coefficients for interpolation
c     in theta direction
c
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE grid_cfc
      USE ppm_cfc

      IMPLICIT NONE
      
      c_1_theta = 0.5_rk
      c_2_theta = 0.5_rk
      c_3_theta = - 1.0_rk / 6.0_rk
      c_4_theta = 1.0_rk / 6.0_rk
      c_5_theta = 0.5_rk

      return

      END SUBROUTINE coefficients_theta

c     ==================================================================



c     ==================================================================

      SUBROUTINE interpolate_theta (q, q_if, is_mass_fraction,
     &     is_positive)

c     ==================================================================
c
c     subroutine to interpolate in theta direction
c
c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE perm_aux_cfc
      USE ppm_cfc

      IMPLICIT NONE
 
      real (kind=rk), intent (in) :: 
     &     q (- 3 : m + 4, n_s - 4 : n_e + 4, o_s - 4 : o_e + 4)
      real (kind=rk), intent (out) ::  q_if (- 1 : m + 2,
     &     n_s - 2 : n_e + 2, o_s - 2 : o_e + 2, 0 : 1) 

      real (kind=rk) :: delta_q_sign

#ifdef SEKORA_PPM
      real (kind=rk), parameter :: c_ppm = 1.25_rk
      real (kind=rk) :: s2c, d2lim, d2p, d2m, d2c
#endif      
      logical, intent(in), optional :: is_mass_fraction, is_positive
      logical :: mass_fraction=.false.,positive=.false.
      
      integer (kind=ik) :: i, j, k, jk

      if (present(is_positive)) then
         positive=.true.
      else
         positive=.false.
      end if
      if (present(is_mass_fraction)) then
         positive=.true.
         mass_fraction=.true.
      else
         mass_fraction=.false.
      end if

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, (n_loc + 5) * o_loc
         
         k = int((jk + n_loc + 5 - 1) / (n_loc + 5))
         j = (n_s - 4) + (jk - (k - 1) * (n_loc + 5))
         k = k + o_s - 1

         do i = 1, m
ccc   a_i - a_i-1 (for delta a_i)
            tmp_ppm (i, j, k) = q (i, j + 1, k) - q (i, j, k)
         enddo
      enddo


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, (n_loc + 4) * o_loc

         k = int((jk + n_loc + 4 - 1) / (n_loc + 4))
         j = (n_s - 3) + (jk - (k - 1) * (n_loc + 4))
         k = k + o_s - 1

         do i = 1, m
ccc   delta a_i (C & W)
            delta_q (i, j, k) = c_1_theta * tmp_ppm (i, j, k) +
     &           c_2_theta * tmp_ppm (i, j - 1, k)
ccc   convert delta a_i to delta_m a_i
            if (tmp_ppm (i, j, k) * tmp_ppm (i, j - 1, k) .gt.
     &           0.0_rk)
     &           then
               delta_q_sign =
     &              delta_q (i, j, k) / abs (delta_q (i, j, k))
               delta_q (i, j, k) =
     &              min (abs (delta_q (i, j, k)), 
     &              2.0_rk * min (abs (tmp_ppm (i, j - 1, k)), 
     &              abs (tmp_ppm (i, j, k)))) * delta_q_sign
            else
               delta_q (i, j, k) = 0.0_rk
            endif
         enddo
      enddo
      
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, (n_loc + 3) * o_loc

         k = int((jk + n_loc + 3 - 1) / (n_loc + 3))
         j = (n_s - 3) + (jk - (k - 1) * (n_loc + 3))
         k = k + o_s - 1

         do i = 1, m
#ifdef SEKORA_PPM
c     6th order interpolation
            q_if (i, j, k, 1) =
     &           (37.0_rk * (q (i, j,k) + q (i, j + 1, k)) -
     &           8.0_rk * (q (i, j - 1,k) + q (i, j + 2, k)) +
     &           (q (i, j - 2, k) + q (i, j + 3, k))) / 60.0_rk
            if ((q_if (i, j, k, 1) - q (i, j ,k)) *
     &         (q (i, j + 1, k) - q_if (i, j, k, 1)) .lt. 0.0_rk) then
               d2c = 3.0_rk * (q (i, j, k) + q (i, j + 1, k)
     &              - 2.0_rk * q_if (i, j, k, 1))
               d2p = q (i, j ,k) + q (i, j + 2, k)
     &              - 2.0_rk * q (i, j + 1, k)
               d2m = q (i, j - 1,k) + q (i, j + 1, k)
     &              - 2.0_rk * q (i, j, k)
               s2c = sign (1.0_rk, d2c)
               d2lim = max (min (c_ppm * s2c * d2m, c_ppm * s2c * d2p,
     &              s2c * d2c), 0.0_rk)
               q_if (i, j, k, 1) =
     &              0.5_rk * (q (i, j, k) + q (i, j + 1, k)) -
     &              s2c * d2lim / 6.0_rk !correct sign of correction term??
            end if
            q_if (i, j + 1, k, 0) = q_if (i, j, k, 1)
#else
c     Classical PPM
ccc   a_i+1/2
            q_if (i, j, k, 1) = q (i, j, k) +
     &           c_5_theta * tmp_ppm (i, j, k) + 
     &              c_3_theta * delta_q (i, j + 1, k) +
     &           c_4_theta * delta_q (i, j, k)
            q_if (i, j + 1, k, 0) = q_if (i, j, k, 1)
#endif /* SEKORA_PPM */
         enddo
#ifdef SEKORA_PPM
         if (positive) then
            do i = 1, m
               if (q_if (i, j, k, 1) .lt. 0.0_rk) then
                  q_if (i, j, k, 1) = q (i, j, k) +
     &                 c_5_theta * tmp_ppm (i, j, k) + 
     &                 c_3_theta * delta_q (i, j + 1, k) +
     &                 c_4_theta * delta_q (i, j, k)
                  q_if (i, j + 1, k, 0) = q_if (i, j, k, 1)
               end if
            end do
         end if
         if (mass_fraction) then
            do i = 1, m
               if (q_if (i, j, k, 1) .gt. 1.0_rk) then
                  q_if (i, j, k, 1) = q (i, j, k) +
     &                 c_5_theta * tmp_ppm (i, j, k) + 
     &                 c_3_theta * delta_q (i, j + 1, k) +
     &                 c_4_theta * delta_q (i, j, k)
                  q_if (i, j + 1, k, 0) = q_if (i, j, k, 1)
               end if
            end do
         end if
#endif /* SEKORA_PPM */

      enddo

      return

      END SUBROUTINE interpolate_theta

c     ==================================================================



c     ==================================================================

      SUBROUTINE cd_steepening_theta (xq, x_if)

c     ==================================================================
c
c     subroutine to detect contact discontinuity and sharpen rho profile there
c
c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE grid_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE parameters_cfc
      USE ppm_cfc
      USE shock_cfc, ONLY: lshock

      IMPLICIT NONE

      real (kind=rk), intent (in) ::
     &     xq (-3:m+4, n_s-4:n_e+4, o_s-4:o_e+4)
      real (kind=rk), intent (inout) ::
     &     x_if(-1:m+2, n_s-2:n_e+2, o_s-2:o_e+2, 0:1)

      real (kind=rk) :: delta_2_rho_minus, delta_2_rho_plus
      real (kind=rk) :: eta, eta_tilde

      real (kind=rk) :: criterion_1, criterion_2,
     &     criterion_3, criterion_4, tmp
      
      real (kind=rk) :: gamma_eos_tmp
      
      integer (kind=ik) :: i, j, k, jk


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, (n_loc + 2) * o_loc
         
         k = int((jk + n_loc + 2 - 1) / (n_loc + 2))
         j = (n_s - 2) + (jk - (k - 1) * (n_loc + 2))
         k = k + o_s - 1

         do i = 1, m

ccc   C & W delta^2 a_i
            delta_2_rho_minus =
     &           (xq (i, j, k) - 2.0_rk * xq (i, j - 1, k) + 
     &           xq (i, j - 2, k)) / 6.0_rk
            delta_2_rho_plus =
     &           (xq (i, j + 2, k) - 2.0_rk * xq (i, j + 1, k) + 
     &           xq (i, j, k)) / 6.0_rk
            
            tmp = xq (i, j + 1, k) - xq (i, j - 1, k)
            if (tmp .eq. 0.0_rk) tmp = 1.0e-10_rk
ccc   check if denominator of second factor in (1.17 top) is 0
               
            criterion_1 = delta_2_rho_plus * delta_2_rho_minus
            criterion_2 = abs (tmp) - epsilon_1 * 
     &           min (abs (xq (i, j + 1, k)),
     &           abs (xq (i, j - 1, k)))
            
            if (criterion_1 .gt. 0.0_rk .or.
     &           criterion_2 .lt. 0.0_rk) then
ccc   second case in (1.17 top)
               eta_tilde = 0.0_rk
            else
ccc   first case in (1.17 top) (minus sign integrated in first term)
               eta_tilde =
     &              (delta_2_rho_minus - delta_2_rho_plus) / tmp
            endif
            
            eta = max (0.0_rk, min (eta_1 * 
     &           (eta_tilde - eta_2), 1.0_rk))
ccc   (1.16)
#ifdef CFC_TRANSPORT2            
            gamma_eos_tmp = c_sound_squared (i, j, k) *
     &           rho (i, j, k) * h (i, j, k) / p (i, j, k)
#else /* CFC_TRANSPORT 2 */
            gamma_eos_tmp = c_sound_squared (i, j, k) *
     &           rho (i, j, k) / p (i, j, k)
#endif /* CFC_TRANSPORT2 */            
            criterion_3 = gamma_eos_tmp * k_cap_0 * 
     &           abs (rho (i, j + 1, k) - rho (i, j - 1, k)) /
     &           min (rho (i, j + 1, k), rho (i, j - 1, k))
            criterion_4 =
     &           abs (p (i, j + 1, k) - p (i, j - 1, k)) /
     &           min (p (i, j + 1, k), p (i, j - 1, k))
               
            if (criterion_3 .lt. criterion_4 .or.
     &           lshock (i, j, k)) eta = 0.0_rk
ccc   gamma * K_0 * density jump < pressure jump --> eta = 0
               
            x_if (i, j, k, 0) = x_if (i, j, k, 0) *
     &           (1.0_rk - eta) + (xq (i, j - 1, k) +
     &           0.5_rk * delta_q (i, j - 1, k)) * eta
            x_if (i, j, k, 1) = x_if (i, j, k, 1) *
     &           (1.0_rk - eta) + (xq (i, j + 1, k) -
     &           0.5_rk * delta_q (i, j + 1, k)) * eta
ccc   (1.15)
         enddo
      enddo
         
      return

      END SUBROUTINE cd_steepening_theta

c     ==================================================================



c     ==================================================================

      SUBROUTINE flatten_theta

c     ==================================================================
c
c     subroutine to flatten profile near strong shocks
c
c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE parameters_cfc
      USE ppm_cfc


      use configure
      IMPLICIT NONE

      real (kind=rk) :: q, tmp

      real (kind=rk) :: f, fx (1:m)

      real (kind=rk) :: delta_v, delta_p_2
      real (kind=rk) :: criterion

      integer (kind=ik):: i, j, k, jk, i_nuc


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, (n_loc + 4) * o_loc
         
         k = int((jk + n_loc + 4 - 1) / (n_loc + 4))
         j = (n_s - 3) + (jk - (k - 1) * (n_loc + 4))
         k = k + o_s - 1

         do i = 1, m
            delta_p (i, j, k) = p (i, j + 1, k) - p (i, j - 1, k)
            delta_v = v_2 (i, j + 1, k) - v_2 (i, j - 1, k)
ccc   can multiply with min (...) as it is positive ---> new criterion
            criterion = epsilon_2_theta * min (p (i, j + 1, k), 
     &           p (i, j - 1, k)) - abs (delta_p (i, j, k))
            
            if (criterion .lt. 0.0_rk.and.delta_v .lt. 0.0_rk) then
ccc   set w_i to 1
               q = 1.0_rk
            else
               q = 0.0_rk
            endif
            
ccc   calculate f_i^tilde
            delta_p_2 = p (i, j + 2, k) - p (i, j - 2, k)
            if (delta_p_2 .eq. 0.0_rk) then
ccc   check if denominator is zero
               if (delta_p (i, j, k) .eq. 0.0_rk) then
ccc   this is somewhat arbitrary but obviously works
                  tmp = - omega_1
               else
                  tmp = 1.0_rk - omega_1
               endif
            else
ccc   regular
               tmp = delta_p (i, j, k) / delta_p_2 - omega_1
            endif
            f_tilde (i, j, k) =
     &           min (1.0_rk, q * max (0.0_rk, tmp * 
     &           omega_2))
         enddo
      enddo

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, (n_loc + 2) * o_loc
         
         k = int((jk + n_loc + 2 - 1) / (n_loc + 2))
         j = (n_s - 2) + (jk - (k - 1) * (n_loc + 2))
         k = k + o_s - 1

         do i = 1, m
            if (delta_p (i, j, k) .lt. 0.0_rk) then
ccc   check for negative delta pressure (check in which direction shock travels)
               f = max (f_tilde (i, j, k), f_tilde (i, j + 1, k))
            else
               f = max (f_tilde (i, j, k), f_tilde (i, j - 1, k))
            endif
            fx (i) = f
            
            rho_if (i, j, k, 0) = rho (i, j, k) * f +
     &           rho_if (i, j, k, 0) * (1.0_rk - f)
            rho_if (i, j, k, 1) = rho (i, j, k) * f +
     &              rho_if (i, j, k, 1) * (1.0_rk - f)
            
            v_1_if (i, j, k, 0) = v_1 (i, j, k) * f +
     &           v_1_if (i, j, k, 0) * (1.0_rk - f)
            v_1_if (i, j, k, 1) = v_1 (i, j, k) * f +
     &           v_1_if (i, j, k, 1) * (1.0_rk - f)
            
            v_2_if (i, j, k, 0) = v_2 (i, j, k) * f +
     &           v_2_if (i, j, k, 0) * (1.0_rk - f)
            v_2_if (i, j, k, 1) = v_2 (i, j, k) * f +
     &           v_2_if (i, j, k, 1) * (1.0_rk - f)
            
            v_3_if (i, j, k, 0) = v_3 (i, j, k) * f +
     &           v_3_if (i, j, k, 0) * (1.0_rk - f)
            v_3_if (i, j, k, 1) = v_3 (i, j, k) * f +
     &              v_3_if (i, j, k, 1) * (1.0_rk - f) 
            
            eps_if (i, j, k, 0) = eps (i, j, k) * f +
     &           eps_if (i, j, k, 0) * (1.0_rk - f)
            eps_if (i, j, k, 1) = eps (i, j, k) * f +
     &           eps_if (i, j, k, 1) * (1.0_rk - f)
            
            p_if (i, j, k, 0) = p (i, j, k) * f +
     &           p_if (i, j, k, 0) * (1.0_rk - f)
            p_if (i, j, k, 1) = p (i, j, k) * f +
     &           p_if (i, j, k, 1) * (1.0_rk - f)
            
            gamm_if (i, j, k, 0) =
     &           gamm (i, j, k) * f +
     &           gamm_if (i, j, k, 0) * (1.0_rk - f)
            gamm_if (i, j, k, 1) =
     &           gamm (i, j, k) * f +
     &           gamm_if (i, j, k, 1) * (1.0_rk - f)
            
         enddo
         do i_nuc = 1, config%qn
            do i = 1, m
               xnu_if (i, j, k, 0, i_nuc) =
     &              xnnu (i, j, k, i_nuc) * fx (i) +
     &              xnu_if (i, j, k, 0, i_nuc) * (1.0_rk - fx (i))
               xnu_if (i, j, k, 1, i_nuc) =
     &              xnnu (i, j, k, i_nuc) * fx (i) +
     &              xnu_if (i, j, k, 1, i_nuc) * (1.0_rk - fx (i))
            enddo
         enddo
      enddo

      return

      END SUBROUTINE flatten_theta

c     ==================================================================
      


c     ==================================================================

      SUBROUTINE monotonize_theta (q, q_if)

c     ==================================================================
c
c     subroutine to monotonize reconstructed quantities
c
c     ------------------------------------------------------------------
      
      USE precision
      USE size_cfc

      IMPLICIT NONE

      real (kind=rk), intent (in) ::
     &     q (- 3 : m + 4, n_s - 4 : n_e + 4, o_s - 4 : o_e + 4)
      real (kind=rk), intent (inout) :: q_if (- 1 : m + 2,
     &     n_s - 2 : n_e + 2, o_s - 2 : o_e + 2, 0 : 1) 

#ifdef SEKORA_PPM
      real (kind=rk) :: qi_old (- 1 : m + 2,
     &     n_s - 2 : n_e + 2, o_s - 2 : o_e + 2, 0 : 1)
      real (kind=rk), parameter :: c_ppm = 1.25_rk
      real (kind=rk) :: s2ppm, d2lim, d2p, d2m, d2c, d2ppm, alpl, almn,
     &     almx, scri, qav
      real (kind=rk) :: dlft, drgt, dqlft, dqrgt
      logical :: bigp, bigm, extremum
#endif      
            
      integer (kind=ik) :: i, j, k, jk

#ifdef SEKORA_PPM
      qi_old (:,:,:,:) = q_if (:,:,:,:)
#endif
     
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, (n_loc + 2) * o_loc
         
         k = int((jk + n_loc + 2 - 1) / (n_loc + 2))
         j = (n_s - 2) + (jk - (k - 1) * (n_loc + 2))
         k = k + o_s - 1

         do i = 1, m

#ifdef SEKORA_PPM
            alpl = q_if (i, j, k, 1) - q (i, j, k)
            almn = q_if (i, j, k, 0) - q (i, j, k)
            bigp = abs(alpl) .gt. 2.0_rk * abs(almn)
            bigm = abs(almn) .gt. 2.0_rk * abs(alpl)
            qav = 1e-10_rk * MAX(ABS(q (i, j, k)),
     &           ABS(q (i, j - 1, k)), ABS(q (i, j + 1, k)))
            extremum = .false.
            if (bigp .or. bigm) then
c     Condition for extrema as in CASTRO
               dlft = q (i, j, k) - qi_old (i, j, k, 0)
               drgt = qi_old (i, j, k, 1) - q (i, j, k)
               dqlft = q (i, j    , k) - q (i, j - 1, k)
               dqrgt = q (i, j + 1, k) - q (i, j, k)
               if (min (abs(dlft), abs (drgt)) .ge.
     &              min (abs(dqlft), abs(dqrgt))) then
                  extremum = dlft * drgt .le. 0.0_rk
               else
                  extremum = dqlft * dqrgt .le. 0.0_rk
               end if
            end if
            if (alpl * almn .ge. 0.0_rk .or. !extremum) then
     &           (q (i, j, k) - q (i, j - 1, k)) *
     &           (q (i, j, k) - q (i, j + 1, k)) .ge. 0.0_rk) then
!     local extremum
               d2ppm = 6.0_rk * (alpl + almn)
c     Use cell averages for estimating 2nd derivative
               d2c = q (i,j+1,k) + q (i,j-1,k) - 2.0_rk * q (i,j  ,k)
               d2m = q (i,j-2,k) + q (i,j  ,k) - 2.0_rk * q (i,j-1,k)
               d2p = q (i,j  ,k) + q (i,j+2,k) - 2.0_rk * q (i,j+1,k)
               s2ppm = sign (1.0_rk, d2ppm)
               d2lim = max (min (s2ppm * d2ppm, c_ppm * s2ppm * d2m,
     &              c_ppm * s2ppm * d2p, c_ppm * s2ppm * d2c), 0.0_rk)
               q_if (i, j, k, 1) = q (i, j, k) +
     &              alpl * abs (d2lim) / (abs(d2ppm) + 1e-99_rk)
               q_if (i, j, k, 0) = q (i, j, k) +
     &              almn * abs (d2lim) / (abs(d2ppm) + 1e-99_rk)

            else if (bigp) then
               almx = - 0.25_rk * alpl ** 2 / (alpl + almn)
               scri = q (i, j - 1, k) - q (i, j, k)
               if (almn*almx .ge. almn*scri) then
                  if (almn*(scri-almn) .ge. qav) then
                     alpl = -2.0_rk * scri - sign (2.0_rk, almn) *
     &                    sqrt (max(scri ** 2 - almn * scri, 0.0_rk))
                  else
                     alpl = - 2.0_rk * almn
                  end if
               end if
               q_if (i, j, k, 1) = q (i, j, k) + alpl
            else if (bigm) then
               almx = - 0.25_rk * almn ** 2 / (alpl + almn)
               scri = q (i, j + 1, k) - q (i, j, k)
               if (alpl*almx .ge. alpl*scri) then
                  if (alpl*(scri-alpl) .ge. qav) then
                     almn = - 2.0_rk * scri - sign (2.0_rk, alpl) *
     &                    sqrt (max(scri ** 2 - alpl * scri, 0.0_rk))
                  else
                     almn = - 2.0_rk * alpl
                  end if
               end if
               q_if (i, j, k, 0) = q (i, j, k) + almn

            end if
            
#else /* SEKORA_PPM */            
            if ((q_if (i, j, k, 1) - q (i, j, k)) * 
     &           (q_if (i, j, k, 0) - q (i, j, k)) .ge. 0.0_rk) then
ccc   first case of (1.10)
               q_if (i, j, k, 0) = q (i, j, k)
               q_if (i, j, k, 1) = q (i, j, k)
            else
               if ((q_if (i, j, k, 1) - q_if (i, j, k, 0)) *
     &              (q_if (i, j, k, 0) - 3.0_rk * q (i, j, k) +
     &              2.0_rk * q_if (i, j, k, 1)) .lt. 0.0_rk) then
ccc   second case of (1.10) -- after a little algebra, bring rhs to the left
                  q_if (i, j, k, 0) = 3.0_rk * q (i, j, k) -
     &                 2.0_rk * q_if (i, j, k, 1)
               endif
               if ((q_if (i, j, k, 1) - q_if (i, j, k, 0)) * 
     &              (- q_if (i, j, k, 1) + 3.0_rk * q (i, j, k) -
     &              2.0_rk * q_if (i, j, k, 0)) .lt. 0.0_rk) then
ccc   third case of (1.10) -- after a little algebra, bring rhs to the left
                  q_if (i, j, k, 1) = 3.0_rk * q (i, j, k) -
     &                 2.0_rk * q_if (i, j, k, 0)
               endif
            endif
#endif /* SEKORA_PPM */            
         enddo
      enddo

      return

      END SUBROUTINE monotonize_theta

c     ==================================================================




c     ==================================================================

      SUBROUTINE flatten_xnu_theta

c     ==================================================================
c     subroutine to additionally flatten mass fraction profiles near
c     strong shocks, and enforce \sum_i X_i=1.0_rk. In contrast to
c     flatten_xnu_r, the discontinuity switch f_xnu is not stored in
c     a global array: This requires more operations, but consumes less
c     memory and allows more efficient cache utilization.
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE parameters_cfc

      use configure
      IMPLICIT NONE

      real (kind=rk) :: f_xnu, f_xnup, f_xnum, f
      real (kind=rk) :: xnu_sum (1 : m, 0 : 1)

      integer (kind=ik) :: i, j, k, jk, jp1, jm1, i_nuc


#ifndef SEKORA_PPM
      do i_nuc = 1, config%qn
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, (n_loc + 2) * o_loc
            
            k = int((jk + n_loc + 2 - 1) / (n_loc + 2))
            j = (n_s - 2) + (jk - (k - 1) * (n_loc + 2))
            k = k + o_s - 1
            
            do i = 1, m 

               if ((xnnu(i,j+1,k,i_nuc)-xnnu(i,j,k,i_nuc)) *
     &              (xnnu(i,j,k,i_nuc)-xnnu(i,j-1,k,i_nuc))
     &              .lt. 0.0_rk) then
                  f_xnu = 1.0_rk
               else
                  f_xnu = 0.0_rk
               endif
               if ((xnnu(i,j+2,k,i_nuc)-xnnu(i,j+1,k,i_nuc)) *
     &              (xnnu(i,j+1,k,i_nuc)-xnnu(i,j,k,i_nuc))
     &              .lt. 0.0_rk) then
                  f_xnup = 1.0_rk
               else
                  f_xnup = 0.0_rk
               endif
               if ((xnnu(i,j,k,i_nuc)-xnnu(i,j-1,k,i_nuc)) *
     &              (xnnu(i,j-1,k,i_nuc)-xnnu(i,j-2,k,i_nuc))
     &              .lt. 0.0_rk) then
                  f_xnum = 1.0_rk
               else
                  f_xnum = 0.0_rk
               endif

               f = max (f_xnum, f_xnup)
               f = 0.5_rk * max (2.0_rk * f_xnu, f)
               
               xnu_if (i, j, k, 0, i_nuc) =
     &                 xnnu (i, j, k, i_nuc) * f +
     &              xnu_if (i, j, k, 0, i_nuc) * (1.0_rk - f)
               xnu_if (i, j, k, 1, i_nuc) =
     &              xnnu (i, j, k, i_nuc) * f +
     &              xnu_if (i, j, k, 1, i_nuc) * (1.0_rk - f)
 
            end do
         end do
      end do
#endif /* SEKORA_PPM */

c     Enforce sum(xnu_if(i,j,k,1:config%qn-1,l))=1.0_rk:
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, (n_loc + 2) * o_loc
            
         k = int((jk + n_loc + 2 - 1) / (n_loc + 2))
         j = (n_s - 2) + (jk - (k - 1) * (n_loc + 2))
         k = k + o_s - 1

         xnu_sum(:,:)=0.0_rk
         do i_nuc = 1, config%qn-1
            do i= 1, m
               xnu_sum (i, 0) =
     &              xnu_sum (i, 0) + xnu_if (i, j, k, 0, i_nuc)
               xnu_sum (i, 1) =
     &              xnu_sum (i, 1) + xnu_if (i, j, k, 1, i_nuc)
            end do
         end do
         do i_nuc = 1, config%qn-1
            do i= 1, m
               xnu_if (i, j, k, 0, i_nuc) =
     &              xnu_if (i, j, k, 0, i_nuc) / xnu_sum (i, 0)
               xnu_if (i, j, k, 1, i_nuc) =
     &              xnu_if (i, j, k, 1, i_nuc) / xnu_sum (i, 1)
            end do
         end do
      end do
      
      return
      
      END SUBROUTINE flatten_xnu_theta

c     ==================================================================


c     ******************************************************************
c     ************************ varphi-direction ************************
c     ******************************************************************


c     ==================================================================

      SUBROUTINE coefficients_varphi

c     ==================================================================
c
c     subroutine to calculate the coefficients for interpolation
c     in varphi direction
c
c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE grid_cfc
      USE ppm_cfc
      
      IMPLICIT NONE
            
      c_1_varphi = 0.5_rk
      c_2_varphi = 0.5_rk
      c_3_varphi = - 1.0_rk / 6.0_rk
      c_4_varphi = 1.0_rk / 6.0_rk
      c_5_varphi = 0.5_rk

      return

      END SUBROUTINE coefficients_varphi

c     ==================================================================



c     ==================================================================

      SUBROUTINE interpolate_varphi (q, q_if, is_mass_fraction,
     &     is_positive)

c     ==================================================================
c
c     subroutine to interpolate in varphi direction
c
c     ------------------------------------------------------------------
      
      USE precision
      
      USE size_cfc
      USE perm_aux_cfc
      USE ppm_cfc
      
      IMPLICIT NONE
      
      real (kind=rk), intent (in) ::
     &     q (- 3 : m + 4, n_s - 4 : n_e + 4, o_s - 4 : o_e + 4)
      real (kind=rk), intent (inout) ::  q_if (- 1 : m + 2,
     &     n_s - 2 : n_e + 2, o_s -2 : o_e + 2, 0 : 1) 

      real (kind=rk) :: delta_q_sign

#ifdef SEKORA_PPM
      real (kind=rk), parameter :: c_ppm = 1.25_rk
      real (kind=rk) :: s2c, d2lim, d2p, d2m, d2c
#endif      
      logical, intent(in), optional :: is_mass_fraction, is_positive
      logical :: mass_fraction=.false.,positive=.false.
      
      integer (kind=ik) :: i, j, k, jk 

      if (present(is_positive)) then
         positive=.true.
      else
         positive=.false.
      end if
      if (present(is_mass_fraction)) then
         positive=.true.
         mass_fraction=.true.
      else
         mass_fraction=.false.
      end if

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * (o_loc + 5)
          
         k = int((jk + n_loc - 1) / n_loc)
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 4

         do i = 1, m

ccc   a_i - a_i-1 (for delta a_i)
            tmp_ppm (i, j, k) = q (i, j, k + 1) - q (i, j, k)
         enddo
      enddo

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * (o_loc + 4)
          
         k = int((jk + n_loc - 1) / n_loc)
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 3

         do i = 1, m
ccc   delta a_i (C & W)
            delta_q (i, j, k) = c_1_varphi * tmp_ppm (i, j, k) +
     &           c_2_varphi * tmp_ppm (i, j, k - 1)
ccc   convert delta a_i to delta_m a_i
            if (tmp_ppm (i, j, k) * tmp_ppm (i, j, k - 1) .gt.
     &           0.0_rk) then
               delta_q_sign =
     &              delta_q (i, j, k) / abs (delta_q (i, j, k))
               delta_q (i, j, k) =
     &              min (abs (delta_q (i, j, k)), 
     &              2.0_rk * min (abs (tmp_ppm (i, j, k - 1)), 
     &              abs (tmp_ppm (i, j, k)))) * delta_q_sign
            else
               delta_q (i, j, k) = 0.0_rk
            endif
         enddo
      enddo

      
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * (o_loc + 3)
          
         k = int((jk + n_loc - 1) / n_loc)
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 3

         do i = 1, m
#ifdef SEKORA_PPM
c     6th order interpolation
            q_if (i, j, k, 1) =
     &           (37.0_rk * (q (i, j, k) + q (i, j, k + 1)) -
     &           8.0_rk * (q (i, j, k - 1) + q (i, j, k + 2)) +
     &           (q (i, j, k - 2) + q (i, j, k + 3))) / 60.0_rk
            if ((q_if (i, j, k, 1) - q (i, j, k)) *
     &        (q (i, j, k + 1) - q_if (i, j, k, 1)) .lt. 0.0_rk) then
               d2c = 3.0_rk * (q (i, j, k) + q (i, j, k + 1)
     &              - 2.0_rk * q_if (i, j, k, 1))
               d2p = q (i, j ,k) + q (i, j, k + 2)
     &              - 2.0_rk * q (i, j, k + 1)
               d2m = q (i, j, k - 1) + q (i, j, k + 1)
     &              - 2.0_rk * q (i, j, k)
               s2c = sign (1.0_rk, d2c)
               d2lim = max (min (c_ppm * s2c * d2m, c_ppm * s2c * d2p,
     &              s2c * d2c), 0.0_rk)
               q_if (i, j, k, 1) =
     &              0.5_rk * (q (i, j, k) + q (i, j, k + 1)) -
     &              s2c * d2lim / 6.0_rk !correct sign of correction term??
            end if
            q_if (i, j, k + 1, 0) = q_if (i, j, k, 1)
#else /* SEKORA_PPM */
ccc   a_i+1/2
            q_if (i, j, k, 1) = q (i, j, k) +
     &           c_5_varphi * tmp_ppm (i, j, k) + 
     &           c_3_varphi * delta_q (i, j, k + 1) +
     &           c_4_varphi * delta_q (i, j, k)
            q_if (i, j, k + 1, 0) = q_if (i, j, k, 1)
#endif /* SEKORA_PPM */            
         enddo
#ifdef SEKORA_PPM
         if (positive) then
            do i = 1, m
               if (q_if (i, j, k, 1) .lt. 0.0_rk) then
                  q_if (i, j, k, 1) = q (i, j, k) +
     &                 c_5_varphi * tmp_ppm (i, j, k) + 
     &                 c_3_varphi * delta_q (i, j, k + 1) +
     &                 c_4_varphi * delta_q (i, j, k)
                  q_if (i, j, k + 1, 0) = q_if (i, j, k, 1)
               end if
            end do
         end if
         if (mass_fraction) then
            do i = 1, m
               if (q_if (i, j, k, 1) .gt. 1.0_rk) then
                  q_if (i, j, k, 1) = q (i, j, k) +
     &                 c_5_varphi * tmp_ppm (i, j, k) + 
     &                 c_3_varphi * delta_q (i, j, k + 1) +
     &                 c_4_varphi * delta_q (i, j, k)
                  q_if (i, j, k + 1, 0) = q_if (i, j, k, 1)
               end if
            end do
         end if
#endif /* SEKORA_PPM */            

      enddo

      return

      END SUBROUTINE interpolate_varphi

c     ==================================================================



c     ==================================================================

      SUBROUTINE cd_steepening_varphi (xq,x_if)

c     ==================================================================
c
c     subroutine to detect contact discontinuity and sharpen
c     rho profile there
c   
c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE grid_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE parameters_cfc
      USE ppm_cfc
      USE shock_cfc, ONLY: lshock
      
      IMPLICIT NONE

      real (kind=rk), intent (in) ::
     &     xq (-3:m+4, n_s-4:n_e+4, o_s-4:o_e+4)
      real (kind=rk), intent (inout) ::
     &     x_if(-1:m+2, n_s-2:n_e+2, o_s-2:o_e+2, 0:1)

      real (kind=rk) :: delta_2_rho_minus, delta_2_rho_plus
      real (kind=rk) :: eta, eta_tilde

      real (kind=rk) :: criterion_1, criterion_2,
     &     criterion_3, criterion_4, tmp
      
      real (kind=rk) :: gamma_eos_tmp
      
      integer (kind=ik) :: i, j, k, jk

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * (o_loc + 2)
         
         k = int((jk + n_loc - 1) / n_loc)
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 2

         do i = 1, m

ccc   C & W delta^2 a_i
            delta_2_rho_minus =
     &           (xq (i, j, k) - 2.0_rk * xq (i, j, k - 1) + 
     &           rho (i, j, k - 2)) / 6.0_rk
            delta_2_rho_plus =
     &           (xq (i, j, k + 2) - 2.0_rk * xq (i, j, k + 1) + 
     &           xq (i, j, k)) / 6.0_rk
            
            tmp = xq (i, j, k + 1) - xq (i, j, k - 1)
            if (tmp .eq. 0.0_rk) tmp = 1.0e-10_rk
ccc   check if denominator of second factor in (1.17 top) is 0
            
            criterion_1 = delta_2_rho_plus * delta_2_rho_minus
            criterion_2 = abs (tmp) - epsilon_1 * 
     &           min (abs (xq (i, j, k + 1)),
     &           abs (xq (i, j, k - 1)))
            
            if (criterion_1 .gt. 0.0_rk .or.
     &           criterion_2 .lt. 0.0_rk) then
ccc   second case in (1.17 top)
               eta_tilde = 0.0_rk
            else
ccc   first case in (1.17 top) (minus sign integrated in first term)
               eta_tilde =
     &              (delta_2_rho_minus - delta_2_rho_plus) / tmp
            endif
            
            eta = max (0.0_rk, min (eta_1 * 
     &           (eta_tilde - eta_2), 1.0_rk))
ccc   (1.16)
#ifdef CFC_TRANSPORT2            
            gamma_eos_tmp = c_sound_squared (i, j, k) *
     &           rho (i, j, k) * h (i, j, k) / p (i, j, k)
#else /* CFC_TRANSPORT 2 */
            gamma_eos_tmp = c_sound_squared (i, j, k) *
     &           rho (i, j, k) / p (i, j, k)
#endif /* CFC_TRANSPORT2 */            
            criterion_3 = gamma_eos_tmp * k_cap_0 * 
     &           abs (rho (i, j, k + 1) - rho (i, j, k - 1)) /
     &           min (rho (i, j, k + 1), rho (i, j, k - 1))
            criterion_4 =
     &           abs (p (i, j, k + 1) - p (i, j, k - 1)) /
     &           min (p (i, j, k + 1), p (i, j, k - 1))
            
            if (criterion_3 .lt. criterion_4 .or.
     &           lshock (i, j, k)) eta = 0.0_rk
ccc   gamma * K_0 * density jump < pressure jump --> eta = 0
            
            x_if (i, j, k, 0) = x_if (i, j, k, 0) *
     &           (1.0_rk - eta) + (xq (i, j, k - 1) +
     &           0.5_rk * delta_q (i, j, k - 1)) * eta
            x_if (i, j, k, 1) = x_if (i, j, k, 1) *
     &           (1.0_rk - eta) + (xq (i, j, k + 1) -
     &           0.5_rk * delta_q (i, j, k + 1)) * eta
ccc   (1.15)
         enddo
      enddo

      return

      END SUBROUTINE cd_steepening_varphi

c     ==================================================================



c     ==================================================================

      SUBROUTINE flatten_varphi

c     ==================================================================
c
c     subroutine to flatten profile near strong shocks
c
c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE parameters_cfc
      USE ppm_cfc


      use configure
      IMPLICIT NONE

      real (kind=rk) :: q, tmp

      real (kind=rk) :: f, fx (1:m)

      real (kind=rk) :: delta_v, delta_p_2
      real (kind=rk) :: criterion

      integer (kind=ik) :: i, j, k, jk, i_nuc

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * (o_loc + 4)
          
         k = int((jk + n_loc - 1) / n_loc)
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 3

         do i = 1, m
            delta_p (i, j, k) = p (i, j, k + 1) - p (i, j, k - 1)
            delta_v = v_3 (i, j, k + 1) - v_3 (i, j, k - 1)
ccc   can multiply with min (...) as it is positive ---> new criterion
            criterion = epsilon_2_varphi * min (p (i, j, k + 1), 
     &           p (i, j, k - 1)) - abs (delta_p (i, j, k))
            
            if (criterion .lt. 0.0_rk.and.delta_v .lt. 0.0_rk) then
ccc   set w_i to 1
               q = 1.0_rk
            else
               q = 0.0_rk
            endif
            
ccc   calculate f_i^tilde
            delta_p_2 = p (i, j, k + 2) - p (i, j, k - 2)
            if (delta_p_2 .eq. 0.0_rk) then
ccc   check if denominator is zero
               if (delta_p (i, j, k) .eq. 0.0_rk) then
ccc   this is somewhat arbitrary but obviously works
                  tmp = - omega_1
               else
                  tmp = 1.0_rk - omega_1
               endif
            else
ccc   regular
               tmp = delta_p (i, j, k) / delta_p_2 - omega_1
            endif
            f_tilde (i, j, k) =
     &           min (1.0_rk, q * max (0.0_rk, tmp * 
     &           omega_2))
         enddo
      enddo


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * (o_loc + 2)
         
         k = int((jk + n_loc - 1) / n_loc)
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 2

         do i = 1, m
            if (delta_p (i, j, k) .lt. 0.0_rk) then
ccc   check for negative delta pressure (check in which direction shock travels)
               f = max (f_tilde (i, j, k), f_tilde (i, j, k + 1))
            else
               f = max (f_tilde (i, j, k), f_tilde (i, j, k - 1))
            endif
            fx (i) = f
            
            rho_if (i, j, k, 0) = rho (i, j, k) * f +
     &           rho_if (i, j, k, 0) * (1.0_rk - f)
            rho_if (i, j, k, 1) = rho (i, j, k) * f +
     &           rho_if (i, j, k, 1) * (1.0_rk - f)
            
            v_1_if (i, j, k, 0) = v_1 (i, j, k) * f +
     &           v_1_if (i, j, k, 0) * (1.0_rk - f)
            v_1_if (i, j, k, 1) = v_1 (i, j, k) * f +
     &           v_1_if (i, j, k, 1) * (1.0_rk - f)
            
            v_2_if (i, j, k, 0) = v_2 (i, j, k) * f +
     &           v_2_if (i, j, k, 0) * (1.0_rk - f)
            v_2_if (i, j, k, 1) = v_2 (i, j, k) * f +
     &           v_2_if (i, j, k, 1) * (1.0_rk - f)
            
            v_3_if (i, j, k, 0) = v_3 (i, j, k) * f +
     &           v_3_if (i, j, k, 0) * (1.0_rk - f)
            v_3_if (i, j, k, 1) = v_3 (i, j, k) * f +
     &           v_3_if (i, j, k, 1) * (1.0_rk - f) 
            
            eps_if (i, j, k, 0) = eps (i, j, k) * f +
     &           eps_if (i, j, k, 0) * (1.0_rk - f)
            eps_if (i, j, k, 1) = eps (i, j, k) * f +
     &           eps_if (i, j, k, 1) * (1.0_rk - f)

            p_if (i, j, k, 0) = p (i, j, k) * f +
     &           p_if (i, j, k, 0) * (1.0_rk - f)
            p_if (i, j, k, 1) = p (i, j, k) * f +
     &           p_if (i, j, k, 1) * (1.0_rk - f)

            gamm_if (i, j, k, 0) = gamm (i, j, k) * f +
     &           gamm_if (i, j, k, 0) * (1.0_rk - f)
            gamm_if (i, j, k, 1) = gamm (i, j, k) * f +
     &           gamm_if (i, j, k, 1) * (1.0_rk - f)
         enddo
         do i_nuc = 1, config%qn
            do i = 1, m
               xnu_if (i, j, k, 0, i_nuc) =
     &              xnnu (i, j, k, i_nuc) * fx (i) +
     &              xnu_if (i, j, k, 0, i_nuc) * (1.0_rk - fx (i))
               xnu_if (i, j, k, 1, i_nuc) =
     &              xnnu (i, j, k, i_nuc) * fx (i) +
     &              xnu_if (i, j, k, 1, i_nuc) * (1.0_rk - fx (i))
            enddo
         enddo
      enddo

      return

      END SUBROUTINE flatten_varphi

c     ==================================================================
      
      
      
c     ==================================================================

      SUBROUTINE monotonize_varphi (q, q_if)

c     ==================================================================
c
c     subroutine to monotonize reconstructed quantities
c   
c     ------------------------------------------------------------------

      USE precision
      USE size_cfc

      IMPLICI TNONE

      real (kind=rk), intent (in) ::
     &     q (- 3 : m + 4, n_s - 4 : n_e + 4, o_s - 4 : o_e + 4)
      real (kind=rk), intent (inout) :: q_if (- 1 : m + 2,
     &     n_s - 2 : n_e + 2, o_s - 2 : o_e + 2, 0 : 1) 

#ifdef SEKORA_PPM
      real (kind=rk) :: qi_old (- 1 : m + 2,
     &     n_s - 2 : n_e + 2, o_s - 2 : o_e + 2, 0 : 1)
      real (kind=rk), parameter :: c_ppm = 1.25_rk
      real (kind=rk) :: s2ppm, d2lim, d2p, d2m, d2c, d2ppm, alpl, almn,
     &     almx, scri, qav
      real (kind=rk) :: dlft, drgt, dqlft, dqrgt
      logical :: bigp, bigm, extremum
#endif      

      integer (kind=ik) :: i, j, k, jk

#ifdef SEKORA_PPM
      qi_old (:,:,:,:) = q_if (:,:,:,:)
#endif
      
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * (o_loc + 2)
         
         k = int((jk + n_loc - 1) / n_loc)
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 2

         do i = 1, m

#ifdef SEKORA_PPM

            alpl = q_if (i, j, k, 1) - q (i, j, k)
            almn = q_if (i, j, k, 0) - q (i, j, k)
            bigp = abs(alpl) .gt. 2.0_rk * abs(almn)
            bigm = abs(almn) .gt. 2.0_rk * abs(alpl)
            qav = 1e-10_rk * MAX(ABS(q (i, j, k)),
     &           ABS(q (i, j, k - 1)), ABS(q (i, j, k + 1)))
            extremum = .false.
            if (bigp .or. bigm) then
c     Condition for extrema as in CASTRO
               dlft = q (i, j, k) - qi_old (i, j, k, 0)
               drgt = qi_old (i, j, k, 1) - q (i, j, k)
               dqlft = q (i, j, k    ) - q (i, j, k - 1)
               dqrgt = q (i, j, k + 1) - q (i, j, k)
               if (min (abs(dlft), abs (drgt)) .ge.
     &              min (abs(dqlft), abs(dqrgt))) then
                  extremum = dlft * drgt .le. 0.0_rk
               else
                  extremum = dqlft * dqrgt .le. 0.0_rk
               end if
            end if
            if (alpl * almn .ge. 0.0_rk .or. !extremum) then
     &           (q (i, j, k) - q (i, j, k - 1)) *
     &           (q (i, j, k) - q (i, j, k + 1)) .ge. 0.0_rk) then
!     local extremum
               d2ppm = 6.0_rk * (alpl + almn)
c     Use cell averages for estimating 2nd derivative
               d2c = q (i,j,k+1) + q (i,j,k-1) - 2.0_rk * q (i,j  ,k)
               d2m = q (i,j,k-2) + q (i,j  ,k) - 2.0_rk * q (i,j,k-1)
               d2p = q (i,j  ,k) + q (i,j,k+2) - 2.0_rk * q (i,j,k+1)
               s2ppm = sign (1.0_rk, d2ppm)
               d2lim = max (min (s2ppm * d2ppm, c_ppm * s2ppm * d2m,
     &              c_ppm * s2ppm * d2p, c_ppm * s2ppm * d2c), 0.0_rk)
               q_if (i, j, k, 1) = q (i, j, k) +
     &              alpl * abs (d2lim) / (abs(d2ppm) + 1e-99_rk)
               q_if (i, j, k, 0) = q (i, j, k) +
     &              almn * abs (d2lim) / (abs(d2ppm) + 1e-99_rk)

            else if (bigp) then
               
               almx = - 0.25_rk * alpl ** 2 / (alpl + almn)
               scri = q (i, j, k - 1) - q (i, j, k)
               if (almn*almx .ge. almn*scri) then
                  if (almn*(scri-almn) .ge. qav) then
                     alpl = -2.0_rk * scri - sign (2.0_rk, almn) *
     &                    sqrt (max(scri ** 2 - almn * scri, 0.0_rk))
                  else
                     alpl = - 2.0_rk * almn
                  end if
               end if
               q_if (i, j, k, 1) = q (i, j, k) + alpl
            else if (bigm) then
               
               almx = - 0.25_rk * almn ** 2 / (alpl + almn)
               scri = q (i, j, k + 1) - q (i, j, k)
               if (alpl*almx .ge. alpl*scri) then
                  if (alpl*(scri-alpl) .ge. qav) then
                     almn = - 2.0_rk * scri - sign (2.0_rk, alpl) *
     &                    sqrt (max(scri ** 2 - alpl * scri, 0.0_rk))
                  else
                     almn = - 2.0_rk * alpl
                  end if
               end if
               q_if (i, j, k, 0) = q (i, j, k) + almn

            end if
            
#else /* SEKORA_PPM */            

            if ((q_if (i, j, k, 1) - q (i, j, k)) * 
     &              (q_if (i, j, k, 0) - q (i, j, k)) .ge. 0.0_rk) then
ccc   first case of (1.10)
               q_if (i, j, k, 0) = q (i, j, k)
               q_if (i, j, k, 1) = q (i, j, k)
            else
               if ((q_if (i, j, k, 1) - q_if (i, j, k, 0)) *
     &              (q_if (i, j, k, 0) - 3.0_rk * q (i, j, k) +
     &              2.0_rk * q_if (i, j, k, 1)) .lt. 0.0_rk) then
ccc   second case of (1.10) -- after a little algebra, bring rhs to the left
                  q_if (i, j, k, 0) = 3.0_rk * q (i, j, k) -
     &                 2.0_rk * q_if (i, j, k, 1)
               endif
               if ((q_if (i, j, k, 1) - q_if (i, j, k, 0)) * 
     &              (- q_if (i, j, k, 1) + 3.0_rk * q (i, j, k) -
     &              2.0_rk * q_if (i, j, k, 0)) .lt. 0.0_rk) then
ccc   third case of (1.10) -- after a little algebra, bring rhs to the left
                  q_if (i, j, k, 1) = 3.0_rk * q (i, j, k) -
     &                 2.0_rk * q_if (i, j, k, 0)
               endif
            endif
#endif /* SEKORA_PPM */            
            
         enddo
      enddo

      return

      END SUBROUTINE monotonize_varphi

c     ==================================================================




c     ==================================================================

      SUBROUTINE flatten_xnu_varphi

c     ==================================================================
c     subroutine to additionally flatten mass fraction profiles near
c     strong shocks, and enforce \sum_i X_i=1.0_rk. In contrast to
c     flatten_xnu_r, the discontinuity switch f_xnu is not stored in
c     a global array: This requires more operations, but consumes less
c     memory and allows more efficient cache utilization.
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE parameters_cfc


      use configure

      IMPLICIT NONE

      real (kind=rk) :: f_xnu, f_xnup, f_xnum, f
      real (kind=rk) :: xnu_sum (1 : m, 0 : 1)

      integer (kind=ik) :: i, j, k, jk, jp1, jm1, i_nuc


#ifndef SEKORA_PPM
      do i_nuc = 1, config%qn
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * (o_loc + 2)
            
            k = int((jk + n_loc - 1) / n_loc)
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 2
            
            do i = 1, m 

               if ((xnnu(i,j,k+1,i_nuc)-xnnu(i,j,k,i_nuc)) *
     &              (xnnu(i,j,k,i_nuc)-xnnu(i,j,k-1,i_nuc))
     &              .lt. 0.0_rk) then
                  f_xnu = 1.0_rk
               else
                  f_xnu = 0.0_rk
               endif
               if ((xnnu(i,j,k+2,i_nuc)-xnnu(i,j,k+1,i_nuc)) *
     &              (xnnu(i,j,k+1,i_nuc)-xnnu(i,j,k,i_nuc))
     &              .lt. 0.0_rk) then
                  f_xnup = 1.0_rk
               else
                  f_xnup = 0.0_rk
               endif
               if ((xnnu(i,j,k,i_nuc)-xnnu(i,j,k-1,i_nuc)) *
     &              (xnnu(i,j,k-1,i_nuc)-xnnu(i,j,k-2,i_nuc))
     &              .lt. 0.0_rk) then
                  f_xnum = 1.0_rk
               else
                  f_xnum = 0.0_rk
               endif

               f = max (f_xnum, f_xnup)
               f = 0.5_rk * max (2.0_rk * f_xnu, f)
               
               xnu_if (i, j, k, 0, i_nuc) =
     &                 xnnu (i, j, k, i_nuc) * f +
     &              xnu_if (i, j, k, 0, i_nuc) * (1.0_rk - f)
               xnu_if (i, j, k, 1, i_nuc) =
     &              xnnu (i, j, k, i_nuc) * f +
     &              xnu_if (i, j, k, 1, i_nuc) * (1.0_rk - f)
 
            end do
         end do
      end do
#endif /* SEKORA_PPM */

c     Enforce sum(xnu_if(i,j,k,1:config%qn-1,l))=1.0_rk:
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * (o_loc + 2)
            
         k = int((jk + n_loc - 1) / n_loc)
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 2

         xnu_sum(:,:)=0.0_rk
         do i_nuc = 1, config%qn-1
            do i= 1, m
               xnu_sum (i, 0) =
     &              xnu_sum (i, 0) + xnu_if (i, j, k, 0, i_nuc)
               xnu_sum (i, 1) =
     &              xnu_sum (i, 1) + xnu_if (i, j, k, 1, i_nuc)
            end do
         end do
         do i_nuc = 1, config%qn-1
            do i= 1, m
               xnu_if (i, j, k, 0, i_nuc) =
     &              xnu_if (i, j, k, 0, i_nuc) / xnu_sum (i, 0)
               xnu_if (i, j, k, 1, i_nuc) =
     &              xnu_if (i, j, k, 1, i_nuc) / xnu_sum (i, 1)
            end do
         end do
      end do
      
      return
      
      END SUBROUTINE flatten_xnu_varphi

c     ==================================================================




      END MODULE gr_ppm

c     ==================================================================
