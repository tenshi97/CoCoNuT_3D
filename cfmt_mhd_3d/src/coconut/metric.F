c     ==================================================================

      MODULE metric

c     ==================================================================

      IMPLICIT NONE

      CONTAINS


c     ==================================================================

      SUBROUTINE calculate_metric

c     ==================================================================

      USE precision

      USE size_cfc
      USE metric_cfc, ONLY:s_cap_star,beta,w_tmp,rho_hat_phi,
     &     rho_hat_phi_1,rho_hat_phi_2,phi,
     &     rho_hat_alpha,rho_hat_alpha_1,alpha,s_cap,
     &     beta_up_1,beta_up_2,beta_up_3,scalar_tmp
      USE metric_aux_cfc, ONLY: drel_tmp,drel
      USE parameters_cfc, ONLY: metric_tolerance, initial_step
      USE hydro_areas_mod
      USE configure

      USE poisson
      USE cons_check
      USE recover_prim_vars

#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      IMPLICIT NONE

      integer (kind=ik) :: i,j,k,jk,iarea,ixi,ixf,ioy,ioz
      integer (kind=ik) :: iter

#ifdef MPI_HYDRO
      integer (kind=ik) :: ierr
#endif

      real (kind=rk) :: drel2

c     ------------------------------------------------------------------
c     Auxiliary vector field W (X in Cordero-Carrion et al., Eq. 30)

#ifdef CFC_TRANSPORT2
      call create_source_vector_w
      call loese_vektor_poisson (s_cap_star,w_tmp)
      call calculate_temp_aux_values_w
#else
      w_tmp = 0.0_rk
#endif /* CFC_TRANSPORT2 */
      call calculate_extrinsic_curvature

#ifdef CFC_TRANSPORT2
c     ------------------------------------------------------------------
c     Conformal factor (Eq. (31) in Cordero-Carrion et al.)

      call create_source_vector_phi

      drel=1.0_rk
      phi_iteration: do iter=1,100

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            scalar_tmp(1:m,j,k)=phi(1:m,j,k)
c            print *,j,k,jk,
c     &           s_cap_star(6,j,k,:),w_tmp(6,j,k,:),'debug: metric0'
         end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            do i=1,m
               rho_hat_phi(i,j,k)=
     &              rho_hat_phi_1(i,j,k)/phi(i,j,k)+
     &              rho_hat_phi_2(i,j,k)/phi(i,j,k)**7
            end do
         end do
         call loese_poisson (rho_hat_phi,phi)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            phi(:,j,k)=phi(:,j,k)+1.0_rk               
            if (j .eq. n_s) phi(:,j-1,k)=phi(:,j-1,k)+1.0_rk
            if (j .eq. n_e) phi(:,j+1,k)=phi(:,j+1,k)+1.0_rk
            if (k .eq. o_s) phi(:,j,k-1)=phi(:,j,k-1)+1.0_rk
            if (k .eq. o_e) phi(:,j,k+1)=phi(:,j,k+1)+1.0_rk
            drel_tmp(j,k)=maxval(abs(phi(1:m,j,k)-scalar_tmp(1:m,j,k))/
     &           phi(1:m,j,k))

         end do

         drel=0.0_rk
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO REDUCTION(max:drel)
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            drel=max(drel,drel_tmp(j,k))

         end do

#ifdef MPI_HYDRO
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
         call MPI_Allreduce (drel, drel2, 1_ik,
     &     MPI_DOUBLE_PRECISION, MPI_MAX, cart_comm, ierr)
         drel=drel2
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif
#endif /*MPI_HYDRO*/

         if (drel .lt. metric_tolerance) exit phi_iteration

      end do phi_iteration
#else
      phi = 1.0_rk
#endif /* CFC_TRANSPORT2 */
      
      call calculate_metric_coefficients (0)
      call calculate_metric_determinant (0)

#ifdef CFC_TRANSPORT2
c     ------------------------------------------------------------------
c     Transform to non-densitized quantities and recover primitives
      call recover_primitives

c     ------------------------------------------------------------------
c     Lapse function (Eq. (32) in Cordero-Carrion et al.)

      call create_source_vector_alpha
      alpha_iteration: do iter=1,100

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            scalar_tmp(1:m,j,k)=alpha(1:m,j,k)
         end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            do i=1,m
               rho_hat_alpha(i,j,k)=
     &              rho_hat_alpha_1(i,j,k)*alpha(i,j,k)
            end do
         end do

         call loese_poisson (rho_hat_alpha,alpha)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc

            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            alpha(:,j,k)=alpha(:,j,k)+1.0_rk
            if (j .eq. n_s) alpha(:,j-1,k)=alpha(:,j-1,k)+1.0_rk
            if (j .eq. n_e) alpha(:,j+1,k)=alpha(:,j+1,k)+1.0_rk
            if (k .eq. o_s) alpha(:,j,k-1)=alpha(:,j,k-1)+1.0_rk
            if (k .eq. o_e) alpha(:,j,k+1)=alpha(:,j,k+1)+1.0_rk
            drel_tmp(j,k)=
     &           maxval(abs(alpha(1:m,j,k)-scalar_tmp(1:m,j,k))/
     &           alpha(1:m,j,k))

         end do

         drel=0.0_rk

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO REDUCTION(max:drel)
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
         
            drel=max(drel,drel_tmp(j,k))

         end do


#ifdef MPI_HYDRO
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
         call MPI_Allreduce (drel, drel2, 1_ik,
     &     MPI_DOUBLE_PRECISION, MPI_MAX, cart_comm, ierr)
         drel=drel2
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif
#endif /*MPI_HYDRO*/

         if (drel .lt. metric_tolerance) exit alpha_iteration

      end do alpha_iteration

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         alpha(:,j,k)=alpha(:,j,k)/phi(:,j,k)               
         if (j .eq. n_s) alpha(:,j-1,k)=alpha(:,j-1,k)/phi(:,j-1,k)
         if (j .eq. n_e) alpha(:,j+1,k)=alpha(:,j+1,k)/phi(:,j+1,k)
         if (k .eq. o_s) alpha(:,j,k-1)=alpha(:,j,k-1)/phi(:,j,k-1)
         if (k .eq. o_e) alpha(:,j,k+1)=alpha(:,j,k+1)/phi(:,j,k+1)

      end do


c     ------------------------------------------------------------------
c     Shift vector (Eq. (33) in Cordero-Carrion et al.)

      call create_source_vector_beta
      call loese_vektor_poisson (s_cap,beta)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         beta_up_1(:,j,k)=beta(:,j,k,1)
         beta_up_2(:,j,k)=beta(:,j,k,2)
         beta_up_3(:,j,k)=beta(:,j,k,3)
         if (j .eq. n_s) then
            beta_up_1(:,j-1,k)=beta(:,j-1,k,1)
            beta_up_2(:,j-1,k)=beta(:,j-1,k,2)
            beta_up_3(:,j-1,k)=beta(:,j-1,k,3)
         end if
         if (j .eq. n_e) then
            beta_up_1(:,j+1,k)=beta(:,j+1,k,1)
            beta_up_2(:,j+1,k)=beta(:,j+1,k,2)
            beta_up_3(:,j+1,k)=beta(:,j+1,k,3)
         end if
         if (k .eq. o_s) then
            beta_up_1(:,j,k-1)=beta(:,j,k-1,1)
            beta_up_2(:,j,k-1)=beta(:,j,k-1,2)
            beta_up_3(:,j,k-1)=beta(:,j,k-1,3)
         end if
         if (k .eq. o_e) then
            beta_up_1(:,j,k+1)=beta(:,j,k+1,1)
            beta_up_2(:,j,k+1)=beta(:,j,k+1,2)
            beta_up_3(:,j,k+1)=beta(:,j,k+1,3)
         end if

      end do

c     ------------------------------------------------------------------
c     Compute everything that is required for the hydrodynamics...

      do iarea=1,areas%are_nu
         ixi    = areas%ix_are(iarea, 1)
         ixf    = areas%ix_are(iarea, 2)
         ioy    = areas%ix_are(iarea, 6)
         ioz    = areas%ix_are(iarea, 9)
         if (ioy .ne. 1 .and. ioy .eq.  config%qy) then
            if (ixi .eq. 1) ixi = ixi-1
            if (ixf .eq. m) ixi = ixi+1
            beta_up_2(ixi:ixf,:,:) = 0.0_rk
            beta_up_3(ixi:ixf,:,:) = 0.0_rk
         end if
      end do
#else
      alpha = 1.0_rk
      beta_up_1 = 0.0_rk
      beta_up_2 = 0.0_rk
      beta_up_3 = 0.0_rk
#endif /* CFC_TRANSPORT2 */

      call calculate_metric_coefficients (1)
      call calculate_metric_determinant (1)
      call calculate_metric_derivatives
      
      return

      END SUBROUTINE calculate_metric

c     ==================================================================




c     ==================================================================

      SUBROUTINE calculate_temp_aux_values_w

c     ==================================================================
c
c     Comment: The varphi-derivative of W will need to be computed
c     once the metric is no longer assumed to be axisymmetric
c
c     ------------------------------------------------------------------
      
      USE precision

      USE size_cfc
      USE metric_cfc, ONLY: w_tmp,dw_tmp_dr,dw_tmp_dtheta
      USE perm_aux_cfc, ONLY: aamp, delta_r_amp_inv, delta_theta_inv

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk, l

      real (kind=rk) :: w_tmp_1,w_tmp_2,w_tmp_3,w_tmp_4,w_tmp_5,
     &     w_tmp_6,w_tmp_7,w_tmp_8,w_tmp_9

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do l = 1, 3
            do i = 1, m 
               
c               w_tmp_1 = w_tmp (i - 1, j - 1, k, l)
               w_tmp_2 = w_tmp (i - 1, j  , k, l)
c               w_tmp_3 = w_tmp (i - 1, j + 1, k, l)
               w_tmp_4 = w_tmp (i    , j - 1, k, l)
               w_tmp_5 = w_tmp (i    , j  , k, l)
               w_tmp_6 = w_tmp (i    , j + 1, k, l)
c               w_tmp_7 = w_tmp (i + 1, j - 1, k, l)
               w_tmp_8 = w_tmp (i + 1, j  , k, l)
c               w_tmp_9 = w_tmp (i + 1, j + 1, k, l)
               
               dw_tmp_dr (i, j, k, l) = delta_r_amp_inv (i) *
     &              (w_tmp_8 + (aamp (i) - 1.0_rk) * w_tmp_5 -
     &              aamp (i) * w_tmp_2)
               dw_tmp_dtheta (i, j, k, l) = (w_tmp_6 - w_tmp_4) *
     &              0.5_rk * delta_theta_inv
                  
            enddo
         enddo
      enddo

      return

      END SUBROUTINE calculate_temp_aux_values_w

c     ==================================================================

      


c     ==================================================================

      SUBROUTINE calculate_extrinsic_curvature

c     ==================================================================

      USE precision

      USE size_cfc
      USE metric_cfc, ONLY: dw_tmp_dr, dw_tmp_dtheta, k_cap_hat,
     &     k_cap_hat_squared, w_tmp
      USE perm_aux_cfc, ONLY: r_inv, cotan_theta, rr, rr_inv,
     &     ssinus_theta, ssinus_theta_inv

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk
      

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m 
            
            k_cap_hat (i, j, k, 1, 1) =
     &           2.0_rk / 3.0_rk *
     &           (2.0_rk * dw_tmp_dr (i, j, k, 1) -
     &           dw_tmp_dtheta (i, j, k, 2) -
     &           2.0_rk * r_inv (i) * w_tmp (i, j, k, 1) -
     &           cotan_theta (j) * w_tmp (i, j, k, 2))

            k_cap_hat (i, j, k, 2, 2) =
     &           2.0_rk / 3.0_rk * rr_inv (i) *
     &           (- dw_tmp_dr (i, j, k, 1) +
     &           2.0_rk * dw_tmp_dtheta (i, j, k, 2) +
     &           r_inv (i) * w_tmp (i, j, k, 1) -
     &           cotan_theta (j) * w_tmp (i, j, k, 2))
            
            k_cap_hat (i, j, k, 3, 3) =
     &           2.0_rk / 3.0_rk * rr_inv (i) * ssinus_theta_inv (j) *
     &           (- dw_tmp_dr (i, j, k, 1) -
     &           dw_tmp_dtheta (i, j, k, 2) +
     &           r_inv (i) * w_tmp (i, j, k, 1) +
     &           2.0_rk * cotan_theta (j) * w_tmp (i, j, k, 2))
            
            k_cap_hat (i, j, k, 1, 2) =
     &           (rr_inv (i) * dw_tmp_dtheta (i, j, k, 1) +
     &           dw_tmp_dr (i, j, k, 2))
            
            k_cap_hat (i, j, k, 1, 3) =
     &           dw_tmp_dr (i, j, k, 3)
            
            k_cap_hat (i, j, k, 2, 3) =
     &           rr_inv (i) * dw_tmp_dtheta (i, j, k, 3)
            
         enddo
      enddo

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif      
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m
            
            k_cap_hat_squared (i, j, k) =
     &           k_cap_hat (i, j, k, 1, 1) ** 2 +
     &           (rr (i) * k_cap_hat (i, j, k, 2, 2)) ** 2 +
     &           (rr (i) * ssinus_theta (j) *
     &           k_cap_hat (i, j, k, 3, 3)) ** 2 +
     &           2.0_rk * (rr (i) * k_cap_hat (i, j, k, 1, 2) ** 2 +
     &           rr (i) * ssinus_theta (j) *
     &           k_cap_hat (i, j, k, 1, 3) ** 2 +
     &           rr (i) * rr (i) * ssinus_theta (j) *
     &           k_cap_hat (i, j, k, 2, 3) ** 2)
            
         enddo
      enddo


      return

      END SUBROUTINE calculate_extrinsic_curvature

c     ==================================================================




c     ==================================================================

      SUBROUTINE create_source_vector_w

c     ==================================================================

      USE precision

      USE size_cfc
      USE conserved_cfc, ONLY: s_1_hat, s_2_hat, s_3_hat
      USE metric_cfc, ONLY: s_cap_star
      USE nutra_cfc, ONLY: fnu
      USE parameters_cfc, ONLY: pi
      USE perm_aux_cfc, ONLY: rr_inv, ssinus_theta_inv

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m

            s_cap_star (i, j, k, 1) = 8.0_rk * pi *
     &           (s_1_hat (i, j, k) + fnu (i, j, k))
            
            s_cap_star (i, j, k, 2) = 8.0_rk * pi *
     &           rr_inv (i) *
     &           s_2_hat (i, j, k)

            s_cap_star (i, j, k, 3) = 8.0_rk * pi *
     &           rr_inv (i) * ssinus_theta_inv (j) *
     &           s_3_hat (i, j, k)

         end do
      end do

      return

      END SUBROUTINE create_source_vector_w

c     ==================================================================




c     ==================================================================

      SUBROUTINE create_source_vector_phi

c     ==================================================================

      USE precision

      USE size_cfc
      USE conserved_cfc, ONLY: d_cap_hat, tau_hat
      USE metric_cfc, ONLY: k_cap_hat_squared, rho_hat_phi_1,
     &     rho_hat_phi_2
      USE nutra_cfc, ONLY: enu
      USE parameters_cfc, ONLY: pi

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m

            rho_hat_phi_1 (i, j, k) = - 2.0_rk * pi *
     &           (tau_hat (i, j, k) + d_cap_hat (i, j, k) + !matter
     &           enu (i, j, k)) !neutrinos
            rho_hat_phi_2 (i, j, k) = - 0.125_rk *
     &           k_cap_hat_squared (i, j, k)

         end do
      end do

      return

      END SUBROUTINE create_source_vector_phi

c     ==================================================================




c     ==================================================================

      SUBROUTINE create_source_vector_alpha

c     ==================================================================

      USE precision

      USE size_cfc
      USE hydro_primitives_cfc, ONLY: rho, h, w, p
      USE metric_cfc, ONLY: k_cap_hat_squared, phi, rho_hat_alpha_1
      USE nutra_cfc, ONLY: enu
      USE parameters_cfc, ONLY: pi

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m

            rho_hat_alpha_1 (i, j, k) = 2.0_rk * pi *
     &           ((rho (i, j, k) * h (i, j, k) *
     &           (3.0_rk * w (i, j, k) ** 2 - 2.0_rk) +
     &           5.0_rk * p (i, j, k)) *
     &           phi (i, j, k) ** 4 + 
     &           3.0_rk *  enu (i, j, k) / phi (i, j, k) ** 2) +
     &           7.0_rk * k_cap_hat_squared (i, j, k) /
     &           (8.0_rk * phi (i, j, k) ** 8)

         end do
      end do

      return

      END SUBROUTINE create_source_vector_alpha

c     ==================================================================



      
c     ==================================================================

      SUBROUTINE create_source_vector_beta

c     ==================================================================

      USE precision

      USE size_cfc
      USE conserved_cfc, ONLY: s_1_hat, s_2_hat, s_3_hat
      USE metric_cfc, ONLY: alpha, phi, sqrt_gamma_inv, s_cap, k_cap_hat
      USE nutra_cfc, ONLY: fnu
      USE parameters_cfc, ONLY: pi
      USE perm_aux_cfc, ONLY: aamp, delta_r_amp_inv, delta_theta_inv,
     &     rr_inv, ssinus_theta_inv, delta_varphi_inv

      USE metric_aux_cfc, ONLY: alpha_phi6, da_phi_dr, da_phi_dtheta,
     &     da_phi_dvarphi

      use configure
      IMPLICIT NONE
           
      integer (kind=ik) :: i, j, k, jk

c     ------------------------------------------------------------------
c     The derivatives of alpha/phi^6 is needed to construct the
c     source density:


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         alpha_phi6 (0 : m + 1, j, k) = alpha (0 : m + 1, j, k) *
     &        sqrt_gamma_inv (0 : m + 1, j, k)
         if (j .eq. n_s) alpha_phi6 (0 : m + 1, j - 1, k) =
     &        alpha (0 : m + 1, j - 1, k) *
     &        sqrt_gamma_inv (0 : m + 1, j - 1, k)
         if (j .eq. n_e) alpha_phi6 (0 : m + 1, j + 1, k) =
     &        alpha (0 : m + 1, j + 1, k) *
     &        sqrt_gamma_inv (0 : m + 1, j + 1, k)
         if (k .eq. o_s) alpha_phi6 (0 : m + 1, j, k - 1) =
     &        alpha (0 : m + 1, j, k - 1) *
     &        sqrt_gamma_inv (0 : m + 1, j, k - 1)
         if (k .eq. o_e) alpha_phi6 (0 : m + 1, j, k + 1) =
     &        alpha (0 : m + 1, j, k + 1) *
     &        sqrt_gamma_inv (0 : m + 1, j, k + 1)

      end do


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m
            da_phi_dr (i, j, k) = delta_r_amp_inv (i) *
     &           (alpha_phi6 (i + 1, j, k) +
     &           (aamp (i) - 1.0_rk) * alpha_phi6 (i, j ,k) -
     &           aamp (i) * alpha_phi6 (i - 1, j, k))
            da_phi_dtheta (i, j, k) =
     &           (alpha_phi6 (i, j + 1, k) -
     &           alpha_phi6 (i, j - 1, k)) *
     &           0.5_rk * delta_theta_inv

            if (config%nsdim .lt. 3) then
               da_phi_dvarphi (i, j, k) =  0.0_rk
            else
            da_phi_dvarphi (i, j, k) =  
     &           (alpha_phi6 (i, j, k + 1) -
     &           alpha_phi6 (i, j, k - 1)) *
     &           0.5_rk * delta_varphi_inv
         endif
         end do
      end do


c     ------------------------------------------------------------------
c     Now construct the source density according to Eq. (33) in 
c     Cordero-Carrion et al. Note: no second derivative of the
c     auxiliary vector field X appears!


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m

            s_cap (i, j, k, 1) = 16.0_rk * pi * alpha (i, j, k) *
     &           ((s_1_hat (i, j, k) + fnu (i, j, k)) /
     &           phi (i, j,k) ** 6) +
     &           2.0_rk * 
     &           (da_phi_dr      (i, j, k) * k_cap_hat (i, j, k, 1, 1) +
     &            da_phi_dtheta  (i, j, k) * k_cap_hat (i, j, k, 1, 2) +
     &            da_phi_dvarphi (i, j, k) * k_cap_hat (i, j, k, 1, 3))

            s_cap (i, j, k, 2) = 16.0_rk * pi * alpha (i, j, k) *
     &           rr_inv (i) *
     &           s_2_hat (i, j, k) / phi (i, j, k) ** 6 +
     &           2.0_rk *
     &           (da_phi_dr      (i, j, k) * k_cap_hat (i, j, k, 1, 2) +
     &            da_phi_dtheta  (i, j, k) * k_cap_hat (i, j, k, 2, 2) +
     &            da_phi_dvarphi (i, j, k) * k_cap_hat (i, j, k, 2, 3))

            s_cap (i, j, k, 3) = 16.0_rk * pi * alpha (i, j, k) *
     &           rr_inv (i) * ssinus_theta_inv (j) *
     &           s_3_hat (i, j, k) / phi (i, j, k) ** 6 + 
     &           2.0_rk *
     &           (da_phi_dr      (i, j, k) * k_cap_hat (i, j, k, 1, 3) +
     &            da_phi_dtheta  (i, j, k) * k_cap_hat (i, j, k, 2, 3) +
     &            da_phi_dvarphi (i, j, k) * k_cap_hat (i, j, k, 3, 3))


         end do
      end do

      return

      END SUBROUTINE create_source_vector_beta

c     ==================================================================




c     ==================================================================

      SUBROUTINE calculate_metric_determinant (mode)

c     ==================================================================
c
c     mode = 0: cell centres
c     mode = 1: cell interfaces
c
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE metric_cfc
      USE grid_cfc, ONLY: r, delta_theta
      USE perm_aux_cfc

      use configure
      IMPLICIT NONE

      integer (kind=ik), intent (in) :: mode

      integer (kind=ik) :: i, j, k, jk, l

      if (mode .eq. 0) then
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            sqrt_gamma(0 : m + 1, j, k) = phi (0 : m + 1, j, k) ** 6
            sqrt_gamma_inv (0 : m + 1, j, k) = 1.0_rk /
     &           sqrt_gamma (0 : m + 1, j, k)
            
            if (j .eq. n_s) then
               sqrt_gamma    (0 : m + 1, j - 1, k) =
     &              phi (0 : m  + 1, j - 1, k) ** 6
               sqrt_gamma_inv(0 : m + 1, j - 1, k) = 1.0_rk /
     &              sqrt_gamma    (0 : m + 1, j - 1, k)
            end if
            if (j .eq. n_e) then
               sqrt_gamma    (0 : m + 1, j + 1, k) =
     &              phi (0 : m  + 1, j + 1, k) ** 6
               sqrt_gamma_inv(0 : m + 1, j + 1, k) = 1.0_rk /
     &              sqrt_gamma    (0 : m + 1, j + 1, k)
            end if

            if (k .eq. o_s) then
               sqrt_gamma    (0 : m + 1, j, k - 1) =
     &              phi (0 : m  + 1, j, k - 1) ** 6
               sqrt_gamma_inv(0 : m + 1, j, k - 1) = 1.0_rk /
     &              sqrt_gamma    (0 : m + 1, j, k - 1)
            end if
            if (k .eq. o_e) then
               sqrt_gamma    (0 : m + 1, j, k + 1) =
     &              phi (0 : m  + 1, j, k + 1) ** 6
               sqrt_gamma_inv(0 : m + 1, j, k + 1) = 1.0_rk /
     &              sqrt_gamma    (0 : m + 1, j, k + 1)
            end if

         enddo
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END DO NOWAIT
#endif
      endif

      if (mode .eq. 1) then
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            do i = 0, m
c     phi^6*dA_r/(dvytot/dvztot)
               sqrt_gamma_r (i, j, k) =
#if MOMENTUM_EQUATION_MODE==2
     &              cos_dtheta *
#endif
     &              phi_r (i, j, k) ** 6 *
     &              rr_if (i)
            enddo
         enddo
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END DO NOWAIT
#endif


         if (config%nsdim .gt. 1) then
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
            do jk = 1, n_loc * o_loc
               
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               
               do i = 1, m
c     phi^6*dA_theta/(dvxtot/dvztot)*r
                  sqrt_gamma_theta (i, j, k) =
     &                 phi_theta (i, j, k) ** 6 *
     &                 r (i) * sinus_theta_if (j) *
     &                 0.5_rk * (rr_if (i) - rr_if (i - 1)) *
     &                 delta_r_mod_inv (i)
               enddo
               
               if (j .eq. n_s) then
                  do i = 1, m
c     phi^6*dA_theta/(dvxtot/dvztot)*r
                     sqrt_gamma_theta (i, j - 1, k) =
     &                    phi_theta (i, j - 1, k) ** 6 *
     &                    r (i) * sinus_theta_if (j - 1) *
     &                    0.5_rk * (rr_if (i) - rr_if (i - 1)) *
     &                    delta_r_mod_inv (i)
                  enddo
               end if
               
            enddo
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END DO NOWAIT
#endif

         else ! config%nsdim .gt 1
            j = 1 ; k = 1
            do i = 1, m
               sqrt_gamma_theta (i, j, k) =
     &              sqrt_gamma (i, j, k) * r (i) *
     &              0.5_rk * (rr_if (i) - rr_if (i - 1)) *
     &              delta_r_mod_inv (i)
            end do
         endif ! config%nsdim .gt 1



         if (config%nsdim .eq. 3) then
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
            do jk = 1, n_loc * o_loc
               
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               
               do i = 1, m

c     phi^6*dA_varphi/(dvxtot/dvytot)*(r sin \theta)
                  sqrt_gamma_varphi (i, j, k) =
     &                 phi_varphi (i, j, k) ** 6 *
     &                 0.5_rk * (rr_if (i) - rr_if (i - 1)) *
     &                 delta_theta *
     &                 r (i) * sinus_theta (j) * 
c#if MOMENTUM_EQUATION_MODE==2
c     &                 cos_dtheta * sin_dtheta *
c#endif
     &                 delta_theta_mod_inv (j) *
     &                 delta_r_mod_inv (i)
               enddo
               
               if (k .eq. o_s) then
                  do i = 1, m
                     sqrt_gamma_varphi (i, j, k - 1) =
     &                    phi_varphi (i, j, k - 1) ** 6 *
     &                    0.5_rk * (rr_if (i) - rr_if (i - 1)) *
     &                    delta_theta *
     &                    r (i) * sinus_theta (j) * 
c#if MOMENTUM_EQUATION_MODE==2
c     &                    cos_dtheta * sin_dtheta *
c#endif
     &                    delta_theta_mod_inv (j) *
     &                    delta_r_mod_inv (i)
                  enddo
               end if
               
            enddo
         else ! config%nsdim .eq. 3
!$OMP DO
            do j = n_s, n_e
               do i = 1, m
                  sqrt_gamma_varphi (i, j, 0:1) =
     &                 phi (i, j, k) ** 6 *
     &                 0.5_rk * (rr_if (i) - rr_if (i - 1)) *
     &                 delta_theta *
     &                 r (i) * sinus_theta (j) *
#if MOMENTUM_EQUATION_MODE==2
c     &                 cos_dtheta * sin_dtheta *
#endif
     &                 delta_theta_mod_inv (j) *
     &                 delta_r_mod_inv (i)
               end do
            end do
         endif ! config%nsdim .eq. 3

      endif
         
      return
 
      END SUBROUTINE calculate_metric_determinant

c     ==================================================================




c     ==================================================================

      SUBROUTINE calculate_metric_coefficients (mode)

c     ==================================================================
c
c     mode = 0: only three-metric at cell centres
c     mode = 1: only metric at cell interfaces
c
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE grid_cfc
      USE metric_cfc
      USE perm_aux_cfc
      USE parameters_cfc, ONLY: equatorial_symmetry

      use configure
      IMPLICIT NONE

      integer (kind=ik), intent (in) :: mode

      real (kind=rk) :: delta_r_tmp (0 : m)

      integer (kind=ik) :: i, j, k, jk

      if (mode .eq. 0) then

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1            

            do i = 1, m
               
               g_11 (i, j, k) = phi (i, j, k) ** 4
               g_22 (i, j, k) = g_11 (i, j, k) * rr (i)
               g_33 (i, j, k) = g_22 (i, j, k) * ssinus_theta (j)
               
               g_up_11 (i, j, k) = 1.0_rk / g_11 (i, j, k)
               g_up_22 (i, j, k) = 1.0_rk / g_22 (i, j, k)
               g_up_33 (i, j, k) = 1.0_rk / g_33 (i, j, k)
               
            enddo
         enddo
      endif

      if (mode .eq. 1) then
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

c     ------------------------------------------------------------------
c     r interface zones
            do i = 0, m
                  
               delta_r_tmp (i) = r_if (i) - r (i)
               
               phi_r (i, j, k) = phi (i, j, k) +
     &              delta_r_tmp (i) * delta_r_if_inv (i) * 
     &              (phi (i + 1, j, k) - phi (i, j, k))
               alpha_r (i, j, k) = alpha (i, j, k) +
     &              delta_r_tmp (i) * delta_r_if_inv (i) * 
     &              (alpha (i + 1, j, k) - alpha (i, j, k)) 
               
               beta_up_1_r (i, j, k) = beta_up_1 (i, j, k) + 
     &              delta_r_tmp (i) * delta_r_if_inv (i) *
     &              (beta_up_1 (i + 1, j, k) - beta_up_1 (i, j, k))
               beta_up_2_r (i, j, k) = beta_up_2 (i, j, k) + 
     &              delta_r_tmp (i) * delta_r_if_inv (i) *
     &              (beta_up_2 (i + 1, j, k) - beta_up_2 (i, j, k))
               beta_up_3_r (i, j, k) = beta_up_3 (i, j, k) + 
     &              delta_r_tmp (i) * delta_r_if_inv (i) *
     &              (beta_up_3 (i + 1, j, k) - beta_up_3 (i, j, k))
               
               g_11_r (i, j, k) = phi_r (i, j, k) ** 4
               g_22_r (i, j, k) = g_11_r (i, j, k) * rr_if (i)
               g_33_r (i, j, k) = g_22_r (i, j, k) * ssinus_theta (j)
               
            enddo

            do i = 1, m
               
               g_up_11_r (i, j, k) = 1.0_rk / g_11_r (i, j, k)
               g_up_22_r (i, j, k) = 1.0_rk / g_22_r (i, j, k)
               g_up_33_r (i, j, k) = 1.0_rk / g_33_r (i, j, k)
               
            enddo
            
            g_up_11_r (0, j, k) = 1.0_rk / g_11_r (0, j, k)
            g_up_22_r (0, j, k) = 0.0_rk
            g_up_33_r (0, j, k) = 0.0_rk
            

            if (config%nsdim .gt. 1) then
c     ------------------------------------------------------------------
c     theta interface zones

               do i = 1, m
                  
                  phi_theta (i, j, k) =
     &                 0.5_rk * (phi (i, j, k) + phi (i, j + 1, k))
                  alpha_theta (i, j, k) =
     &                 0.5_rk * (alpha (i, j, k) + alpha (i, j + 1, k))
                  
                  beta_up_1_theta (i, j, k) = 0.5_rk *
     &                 (beta_up_1 (i, j, k) + beta_up_1 (i, j + 1, k))
                  beta_up_2_theta (i, j, k) = 0.5_rk *
     &                 (beta_up_2 (i, j, k) + beta_up_2 (i, j + 1, k))
                  beta_up_3_theta (i, j, k) = 0.5_rk *
     &                 (beta_up_3 (i, j, k) + beta_up_3 (i, j + 1, k))
                  
                  g_11_theta (i, j, k) = phi_theta (i, j, k) ** 4
                  g_22_theta (i, j, k) = g_11_theta (i, j, k) * rr (i)
                  g_33_theta (i, j, k) = g_22_theta (i, j, k) *
     &                 ssinus_theta_if (j)
                  
               enddo
               
               if (j .eq. n_s) then
                  do i = 1, m
                     
                     phi_theta (i, j - 1, k) =
     &                    0.5_rk * (phi (i, j - 1, k) + phi (i, j, k))
                     alpha_theta (i, j - 1, k) =
     &                  0.5_rk * (alpha (i, j - 1, k) + alpha (i, j, k))
                     
                     beta_up_1_theta (i, j - 1, k) = 0.5_rk *
     &                   (beta_up_1 (i, j - 1, k) + beta_up_1 (i, j, k))
                     beta_up_2_theta (i, j - 1, k) = 0.5_rk *
     &                   (beta_up_2 (i, j - 1, k) + beta_up_2 (i, j, k))
                     beta_up_3_theta (i, j - 1, k) = 0.5_rk *
     &                 (beta_up_3 (i, j - 1, k) + beta_up_3 (i, j, k))
               
                     g_11_theta (i, j - 1, k) =
     &                    phi_theta (i, j - 1, k) ** 4
                     g_22_theta (i, j - 1, k) =
     &                    g_11_theta (i, j - 1, k) * rr (i)
                     g_33_theta (i, j - 1, k) =
     &                    g_22_theta (i, j - 1, k) *
     &                    ssinus_theta_if (j - 1)
                     
                  enddo
               endif
               
               if (j .eq. n_s .and. j .ne. 1) then
                  do i = 1, m
                     
                     g_up_11_theta (i, j - 1, k) =
     &                    1.0_rk /g_11_theta (i, j - 1, k)
                     g_up_22_theta (i, j - 1, k) =
     &                    1.0_rk / g_22_theta (i, j - 1, k)
                     g_up_33_theta (i, j - 1, k) =
     &                    1.0_rk / g_33_theta (i, j - 1, k)
                     
                  enddo
               endif
               
               if (j. ne. n) then  
                  do i = 1, m
                     
                     g_up_11_theta (i, j, k) =
     &                    1.0_rk /g_11_theta (i, j, k)
                     g_up_22_theta (i, j, k) =
     &                    1.0_rk / g_22_theta (i, j, k)
                     g_up_33_theta (i, j, k) =
     &                    1.0_rk / g_33_theta (i, j, k)
                     
                  enddo
               endif
               
               if (j .eq. 1) then
                  do i = 1, m
                     
                     g_up_11_theta (i, 0, k) =
     &                    1.0_rk / g_11_theta (i, 0, k)
                     g_up_22_theta (i, 0, k) =
     &                    1.0_rk / g_22_theta (i, 0, k)
                     g_up_33_theta (i, 0, k) = 0.0_rk
                     
                  enddo
               endif
               
               if (j .eq. n .and. equatorial_symmetry) then
                  do i = 1, m
                     g_up_11_theta (i, n, k) =
     &                    1.0_rk / g_11_theta (i, n, k)
                     g_up_22_theta (i, n, k) =
     &                    1.0_rk / g_22_theta (i, n, k)
                     g_up_33_theta (i, n, k) =
     &                    1.0_rk / g_33_theta (i, n, k)
                  enddo
               else if (j. eq. n .and. .not. equatorial_symmetry) then
                  do i = 1, m
                     g_up_11_theta (i, n, k) = 1.0_rk /
     &                    g_11_theta (i, n, k)
                     g_up_22_theta (i, n, k) = 1.0_rk /
     &                    g_22_theta (i, n, k)
                     g_up_33_theta (i, n, k) = 0.0_rk
                  enddo
               endif
               
            endif ! config%nsdim .gt. 1



            if (config%nsdim .eq. 3) then
c     varphi interface zones

            do i = 1, m
                  
               phi_varphi (i, j, k) =
     &              0.5_rk * (phi (i, j, k) + phi (i, j, k + 1))
               alpha_varphi (i, j, k) =
     &              0.5_rk * (alpha (i, j, k) + alpha (i, j, k + 1))
               
               beta_up_1_varphi (i, j, k) = 0.5_rk *
     &              (beta_up_1 (i, j, k) + beta_up_1 (i, j, k + 1))
               beta_up_2_varphi (i, j, k) = 0.5_rk *
     &              (beta_up_2 (i, j, k) + beta_up_2 (i, j, k + 1))
               beta_up_3_varphi (i, j, k) = 0.5_rk *
     &              (beta_up_3 (i, j, k) + beta_up_3 (i, j, k + 1))
               
               g_11_varphi (i, j, k) = phi_varphi (i, j, k) ** 4
               g_22_varphi (i, j, k) = g_11_varphi (i, j, k) * rr (i)
               g_33_varphi (i, j, k) = g_22_varphi (i, j, k) *
     &              ssinus_theta (j)
               
               g_up_11_varphi (i, j, k) =
     &              1.0_rk / g_11_varphi (i, j, k)
               g_up_22_varphi (i, j, k) =
     &              1.0_rk / g_22_varphi (i, j, k)
               g_up_33_varphi (i, j, k) =
     &              1.0_rk / g_33_varphi (i, j, k)

            enddo

            if (k .eq. o_s) then
               do i = 1, m
                  
                  phi_varphi (i, j, k - 1) =
     &                 0.5_rk * (phi (i, j, k - 1) + phi (i, j, k))
                  alpha_varphi (i, j, k - 1) =
     &                 0.5_rk * (alpha (i, j, k - 1) + alpha (i, j, k))
                  
                  beta_up_1_varphi (i, j, k - 1) = 0.5_rk *
     &                 (beta_up_1 (i, j, k - 1) + beta_up_1 (i, j, k))
                  beta_up_2_varphi (i, j, k - 1) = 0.5_rk *
     &                 (beta_up_2 (i, j, k - 1) + beta_up_2 (i, j, k))
                  beta_up_3_varphi (i, j, k - 1) = 0.5_rk *
     &                 (beta_up_3 (i, j, k - 1) + beta_up_3 (i, j, k))
                  
                  g_11_varphi (i, j, k - 1) =
     &                 phi_varphi (i, j, k - 1) ** 4
                  g_22_varphi (i, j, k - 1) = 
     &                 g_11_varphi (i, j, k - 1) * rr (i)
                  g_33_varphi (i, j, k - 1) =
     &                 g_22_varphi (i, j, k - 1) * ssinus_theta (j)
                  
                  g_up_11_varphi (i, j, k - 1) =
     &                 1.0_rk / g_11_varphi (i, j, k - 1)
                  g_up_22_varphi (i, j, k - 1) =
     &                 1.0_rk / g_22_varphi (i, j, k - 1)
                  g_up_33_varphi (i, j, k - 1) =
     &                 1.0_rk / g_33_varphi (i, j, k - 1)

               enddo
            endif
 
         endif ! config%nsdim .eq. 3

         enddo
      endif

      return

      END SUBROUTINE calculate_metric_coefficients

c     ==================================================================




c     ==================================================================

      SUBROUTINE calculate_metric_derivatives

c     ==================================================================

      USE precision

      USE size_cfc
      USE grid_cfc
      USE metric_cfc
      USE perm_aux_cfc
      USE sources_cfc

      use configure

      IMPLICIT NONE

      real (kind=rk) :: dphi_dr, dphi_dtheta, dphi_dvarphi
      real (kind=rk) :: dphi_dr_o_phi, dphi_dtheta_o_phi,
     &     dphi_dvarphi_o_phi

      real (kind=rk) :: dbeta_up_1_dr, dbeta_up_1_dtheta,
     &     dbeta_up_1_dvarphi
      real (kind=rk) :: dbeta_up_2_dr, dbeta_up_2_dtheta,
     &     dbeta_up_2_dvarphi
      real (kind=rk) :: dbeta_up_3_dr, dbeta_up_3_dtheta,
     &     dbeta_up_3_dvarphi

      integer (kind=ik) :: i, j, k, jk

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m
               
c     calculate metric derivatives
            
            dphi_dr = delta_r_inv (i) *
     &           (phi_r (i, j, k) -
     &           phi_r (i - 1, j, k))
            dphi_dtheta = delta_theta_inv *
     &           (phi_theta (i, j, k) -
     &           phi_theta (i, j - 1, k))

            if (config%nsdim .eq. 3) then

            dphi_dvarphi = delta_varphi_inv *
     &           (phi_varphi (i, j, k) -
     &           phi_varphi (i, j, k - 1))
            else
               dphi_dvarphi = 0.0_rk
            endif
            
            dalpha_dr (i, j, k) = delta_r_inv (i) *
     &           (alpha_r (i, j, k) -
     &           alpha_r (i - 1, j, k))
            dalpha_dtheta (i, j, k) = delta_theta_inv *
     &           (alpha_theta (i, j, k) -
     &           alpha_theta (i, j - 1, k))
            if (config%nsdim .eq. 3) then
            dalpha_dvarphi (i, j, k) = delta_varphi_inv *
     &           (alpha_varphi (i, j, k) -
     &           alpha_varphi (i, j, k - 1))
            else
            dalpha_dvarphi = 0.0_rk
            endif
            
            dbeta_up_1_dr = delta_r_inv (i) *
     &           (beta_up_1_r (i, j, k) -
     &           beta_up_1_r (i - 1, j, k))
            dbeta_up_1_dtheta = delta_theta_inv *
     &           (beta_up_1_theta (i, j, k) -
     &           beta_up_1_theta (i, j - 1, k))
            if (config%nsdim .eq. 3) then
            dbeta_up_1_dvarphi = delta_varphi_inv *
     &           (beta_up_1_varphi (i, j, k) -
     &           beta_up_1_varphi (i, j, k - 1))

            else
            dbeta_up_1_dvarphi = 0.0_rk
         endif
            
            dbeta_up_2_dr = delta_r_inv (i) *
     &           (beta_up_2_r (i, j, k) -
     &           beta_up_2_r (i - 1, j, k))
            dbeta_up_2_dtheta = delta_theta_inv *
     &           (beta_up_2_theta (i, j, k) -
     &           beta_up_2_theta (i, j - 1, k))
            if (config%nsdim .eq. 3) then
            dbeta_up_2_dvarphi = delta_varphi_inv *
     &           (beta_up_2_varphi (i, j, k) -
     &           beta_up_2_varphi (i, j, k - 1))
            else
            dbeta_up_2_dvarphi = 0.0_rk
            endif
            
            dbeta_up_3_dr = delta_r_inv (i) *
     &           (beta_up_3_r (i, j, k) -
     &           beta_up_3_r (i - 1, j, k))
            dbeta_up_3_dtheta = delta_theta_inv *
     &           (beta_up_3_theta (i, j, k) -
     &           beta_up_3_theta (i, j - 1, k))
            if (config%nsdim .eq. 3) then
            dbeta_up_3_dvarphi = delta_varphi_inv *
     &           (beta_up_3_varphi (i, j, k) -
     &           beta_up_3_varphi (i, j, k - 1))
            else
            dbeta_up_3_dvarphi = 0.0_rk
         endif
            
            dphi_dr_o_phi = dphi_dr / phi (i, j, k)
            dphi_dtheta_o_phi = dphi_dtheta / phi (i, j, k)
            dphi_dvarphi_o_phi = dphi_dvarphi / phi (i, j, k)
            
            alpha_dalpha_dr (i, j, k) =
     &           alpha (i, j, k) * dalpha_dr (i, j, k)
            alpha_dalpha_dtheta (i, j, k) =
     &           alpha (i, j, k) * dalpha_dtheta (i, j, k)
            alpha_dalpha_dvarphi (i, j, k) =
     &           alpha (i, j, k) * dalpha_dvarphi (i, j, k)
            
            dln_sqrt_gamma_dr (i, j, k) =
     &           2.0_rk * r_inv (i) + 6.0_rk * dphi_dr_o_phi
            dln_sqrt_gamma_dtheta (i, j, k) =
     &           cotan_theta (j) + 6.0_rk * dphi_dtheta_o_phi
            dln_sqrt_gamma_dvarphi (i, j, k) =
     &           6.0_rk * dphi_dvarphi_o_phi
            
            beta_up_k_k (i, j, k) =
     &           dbeta_up_1_dr  +
     &           dbeta_up_2_dtheta +
     &           dbeta_up_3_dvarphi +
     &           dln_sqrt_gamma_dr (i, j, k) *
     &           beta_up_1 (i, j, k) +
     &           dln_sqrt_gamma_dtheta (i, j, k) *
     &           beta_up_2 (i, j, k) +
     &           dln_sqrt_gamma_dvarphi (i, j, k) *
     &           beta_up_3 (i, j, k)
            
            
            
            excrv (i, j, k, 1, 1) =
     &           g_11 (i, j, k) / (3.0_rk * alpha (i, j, k)) *
     &           (2.0_rk * dbeta_up_1_dr -
     &           dbeta_up_2_dtheta -
     &           dbeta_up_3_dvarphi -
     &           2.0_rk * r_inv (i) * beta_up_1 (i, j, k) -
     &           cotan_theta (j) * beta_up_2 (i, j, k))
            excrv (i, j, k, 2, 2) =
     &           g_22 (i, j, k) / (3.0_rk * alpha (i, j, k)) *
     &           (- dbeta_up_1_dr +
     &           2.0_rk * dbeta_up_2_dtheta -
     &           dbeta_up_3_dvarphi +
     &           r_inv (i) * beta_up_1 (i, j, k) -
     &           cotan_theta (j) * beta_up_2 (i, j, k))
            excrv (i, j, k, 3, 3) =
     &           g_33 (i, j, k) / (3.0_rk * alpha (i, j, k)) *
     &           (- dbeta_up_1_dr -
     &           dbeta_up_2_dtheta +
     &           2.0_rk * dbeta_up_3_dvarphi +
     &           r_inv (i) * beta_up_1 (i, j, k) +
     &           2.0_rk * cotan_theta (j) * beta_up_2 (i, j, k))
            excrv (i, j, k, 1, 2) =
     &           g_22 (i, j, k) / (2.0_rk * alpha (i, j, k)) *
     &           (rr_inv (i) * dbeta_up_1_dtheta +
     &           dbeta_up_2_dr)
            excrv (i, j, k, 1, 3) =
     &           g_33 (i, j, k) / (2.0_rk * alpha (i, j, k)) *
     &           (dbeta_up_3_dr +
     &           rr_inv (i) * ssinus_theta_inv (j) *
     &           dbeta_up_1_dvarphi)
            excrv (i, j, k, 2, 3) =
     &           g_33 (i, j, k) / (2.0_rk * alpha (i, j, k)) *
     &           (dbeta_up_3_dtheta +
     &           ssinus_theta_inv (j) * dbeta_up_2_dvarphi)
            excrv (i, j, k, 2, 1) =
     &           excrv (i, j, k, 1, 2)
            excrv (i, j, k, 3, 1) =
     &           excrv (i, j, k, 1, 3)
            excrv (i, j, k, 3, 2) =
     &           excrv (i, j, k, 2, 3)
            
               
               
            excrv_beta (i, j, k) = 
     &           excrv (i, j, k, 1, 1) *
     &           beta_up_1 (i, j, k) ** 2 +
     &           excrv (i, j, k, 2, 2) *
     &           beta_up_2 (i, j, k) ** 2 +
     &           excrv (i, j, k, 3, 3) *
     &           beta_up_3 (i, j, k) ** 2 +
     &           excrv (i, j, k, 1, 2) *
     &           beta_up_1 (i, j, k) * beta_up_2 (i, j, k) +
     &           excrv (i, j, k, 1, 3) *
     &           beta_up_1 (i, j, k) * beta_up_3 (i, j, k) +
     &           excrv (i, j, k, 2, 3) *
     &           beta_up_2 (i, j, k) * beta_up_3 (i, j, k) +
     &           excrv (i, j, k, 2, 1) *
     &           beta_up_2 (i, j, k) * beta_up_1 (i, j, k) +
     &           excrv (i, j, k, 3, 1) *
     &           beta_up_3 (i, j, k) * beta_up_1 (i, j, k) +
     &           excrv (i, j, k, 3, 2) *
     &           beta_up_3 (i, j, k) * beta_up_2 (i, j, k)
            
            excrv2 (i, j, k) =
     &           g_up_11 (i, j, k) ** 2 *
     &           excrv (i, j, k, 1, 1) ** 2 +
     &           g_up_22 (i, j, k) ** 2 *
     &           excrv (i, j, k, 2, 2) ** 2 +
     &           g_up_33 (i, j, k) ** 2 *
     &           excrv (i, j, k, 3, 3) ** 2 +
     &           2.0_rk * g_up_11 (i, j, k) * g_up_22 (i, j, k) *
     &           excrv (i, j, k, 1, 2) ** 2 +
     &           2.0_rk * g_up_11 (i, j, k) * g_up_33 (i, j, k) *
     &           excrv (i, j, k, 1, 3) ** 2 +
     &           2.0_rk * g_up_22 (i, j, k) * g_up_33 (i, j, k) *
     &           excrv (i, j, k, 2, 3) ** 2
            
            dg_00_dr (i, j, k) =
     &           - 2.0_rk * alpha_dalpha_dr (i, j, k) +
     &           2.0_rk * g_11 (i, j, k) *
     &           (beta_up_1 (i, j, k) *
     &           dbeta_up_1_dr +
     &           rr (i) * beta_up_2 (i, j, k) *
     &           dbeta_up_2_dr +
     &           rr (i) * ssinus_theta (j) * beta_up_3 (i, j, k) *
     &           dbeta_up_3_dr +
     &           r (i) * beta_up_2 (i, j, k) ** 2 +
     &           r (i) * ssinus_theta (j) *
     &           beta_up_3 (i, j, k) ** 2 +
     &           2.0_rk * dphi_dr_o_phi *
     &           (beta_up_1 (i, j, k) ** 2 +
     &           rr (i) * beta_up_2 (i, j, k) ** 2 +
     &           rr (i) * ssinus_theta (j) *
     &           beta_up_3 (i, j, k) ** 2))
            dg_11_dr (i, j, k) =
     &           4.0_rk * g_11 (i, j, k) * dphi_dr_o_phi
            dg_22_dr (i, j, k) =
     &           4.0_rk * g_22 (i, j, k) *
     &           (dphi_dr_o_phi + 0.5_rk * r_inv (i))
            dg_33_dr (i, j, k) =
     &           4.0_rk * g_33 (i, j, k) *
     &           (dphi_dr_o_phi + 0.5_rk * r_inv (i))
            
            dg_00_dtheta (i, j, k) =
     &           - 2.0_rk * alpha_dalpha_dtheta (i, j, k) +
     &           2.0_rk * g_11 (i, j, k) *
     &           (beta_up_1 (i, j, k) *
     &           dbeta_up_1_dtheta +
     &           rr (i) * beta_up_2 (i, j, k) *
     &           dbeta_up_2_dtheta +
     &           rr (i) * ssinus_theta (j) * beta_up_3 (i, j, k) *
     &           dbeta_up_3_dtheta +
     &           rr (i) * sinus_theta (j) * cosin_theta (j) *
     &           beta_up_3 (i, j, k) ** 2 +
     &           2.0_rk * dphi_dtheta_o_phi *
     &           (beta_up_1 (i, j, k) ** 2 +
     &           rr (i) * beta_up_2 (i, j, k) ** 2 +
     &           rr (i) * ssinus_theta (j) *
     &           beta_up_3 (i, j, k) ** 2))
            dg_11_dtheta (i, j, k) =
     &           4.0_rk * g_11 (i, j, k) *
     &           dphi_dtheta_o_phi
            dg_22_dtheta (i, j, k) =
     &           4.0_rk * g_22 (i, j, k) *
     &           dphi_dtheta_o_phi
            dg_33_dtheta (i, j, k) =
     &           4.0_rk * g_33 (i, j, k) *
#if MOMENTUM_EQUATION_MODE==2
     &           (dphi_dtheta_o_phi + 0.5_rk * cotan_theta (j))
#else
     &           (dphi_dtheta_o_phi + 0.5_rk *
     &           (sinus_theta_if (j) - sinus_theta_if (j - 1)) *
     &           delta_theta_inv * sinus_theta_inv (j))
#endif /* MOMENTUM_EQUATION_MODE==2 */
            
            dg_00_dvarphi (i, j, k) =
     &           - 2.0_rk * alpha_dalpha_dvarphi (i, j, k) +
     &           2.0_rk * g_11 (i, j, k) *
     &           (beta_up_1 (i, j, k) *
     &           dbeta_up_1_dvarphi +
     &           rr (i) * beta_up_2 (i, j, k) *
     &           dbeta_up_2_dvarphi +
     &           rr (i) * ssinus_theta (j) * beta_up_3 (i, j, k) *
     &           dbeta_up_3_dvarphi +
     &           2.0_rk * dphi_dvarphi_o_phi *
     &           (beta_up_1 (i, j, k) ** 2 +
     &           rr (i) * beta_up_2 (i, j, k) ** 2 +
     &           rr (i) * ssinus_theta (j) *
     &           beta_up_3 (i, j, k) ** 2))
            dg_11_dvarphi (i, j, k) =
     &           4.0_rk * g_11 (i, j, k) *
     &           dphi_dvarphi_o_phi
            dg_22_dvarphi (i, j, k) =
     &           4.0_rk * g_22 (i, j, k) *
     &           dphi_dvarphi_o_phi
            dg_33_dvarphi (i, j, k) =
     &           4.0_rk * g_33 (i, j, k) *
     &           dphi_dvarphi_o_phi
            
            dg_01_dr (i, j, k) =
     &           g_11 (i, j, k) * (dbeta_up_1_dr +
     &           4.0_rk * dphi_dr_o_phi * beta_up_1 (i, j, k))
            dg_02_dr (i, j, k) =
     &           g_22 (i, j, k) * (dbeta_up_2_dr +
     &           4.0_rk * dphi_dr_o_phi * beta_up_2 (i, j, k) +
     &           2.0_rk * r_inv (i) * beta_up_2 (i, j, k))
            dg_03_dr (i, j, k) =
     &           g_33 (i, j, k) * (dbeta_up_3_dr +
     &           4.0_rk * dphi_dr_o_phi * beta_up_3 (i, j, k) +
     &           2.0_rk * r_inv (i) * beta_up_3 (i, j, k))
            
            dg_01_dtheta (i, j, k) =
     &           g_11 (i, j, k) * (dbeta_up_1_dtheta +
     &           4.0_rk * dphi_dtheta_o_phi * beta_up_1 (i, j, k))
            dg_02_dtheta (i, j, k) =
     &           g_22 (i, j, k) * (dbeta_up_2_dtheta +
     &           4.0_rk * dphi_dtheta_o_phi * beta_up_2 (i, j, k))
            dg_03_dtheta (i, j, k) =
     &           g_33 (i, j, k) * (dbeta_up_3_dtheta +
     &           4.0_rk * dphi_dtheta_o_phi * beta_up_3 (i, j, k) +
     &           2.0_rk * cotan_theta (j) * beta_up_3 (i, j, k))
            
            if (config%nsdim .eq. 3) then
            dg_01_dvarphi (i, j, k) =
     &           g_11 (i, j, k) * (dbeta_up_1_dvarphi +
     &           4.0_rk * dphi_dvarphi_o_phi * beta_up_1 (i, j, k))
            dg_02_dvarphi (i, j, k) =
     &           g_22 (i, j, k) * (dbeta_up_2_dvarphi +
     &           4.0_rk * dphi_dvarphi_o_phi * beta_up_2 (i, j, k))
            dg_03_dvarphi (i, j, k) =
     &           g_33 (i, j, k) * (dbeta_up_3_dvarphi +
     &           4.0_rk * dphi_dvarphi_o_phi * beta_up_3 (i, j, k))
            endif

         enddo
      enddo

      return
      
      END SUBROUTINE calculate_metric_derivatives

c     ==================================================================




      END MODULE metric

c     ===============================================================
