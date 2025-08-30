c     ==================================================================

      MODULE gr_reconstruction

c     ==================================================================

      IMPLICIT NONE

      CONTAINS



c     ==================================================================

      SUBROUTINE reconstruct_interface_r

c     ==================================================================
c     reconstructs interface values in r direction
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE grid_cfc
      USE parameters_cfc
      USE hydro_primitives_cfc      
      USE metric_cfc      
      USE interface_hydro_cfc

      USE gr_ppm
      USE gr_rescale

      use configure
      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk, l, i_nuc, m0, m1

      call interpolate_r (rho, rho_if, is_positive=.true.)
#ifdef CD_STEEPENING
      call cd_steepening_r (rho, rho_if)
#endif
      
      call interpolate_r (v_1, v_1_if)
      call interpolate_r (v_2, v_2_if)
      call interpolate_r (v_3, v_3_if)
      call interpolate_r (eps, eps_if)!, is_positive=.true.)
      do i_nuc=1,config%qn
c         call interpolate_r(xnnu(:,:,:,i_nuc),
c     &        xnu_if(:,:,:,:,i_nuc))
         call interpolate_r(xnnu(-3,n_s-4,o_s-4,i_nuc),
     &        xnu_if(-1,n_s-2,o_s-2,0,i_nuc),is_mass_fraction=.true.)
#ifdef CD_STEEPENING
c         call cd_steepening_r (xnnu(:,:,:,i_nuc),
c     &        xnu_if(:,:,:,:,i_nuc))
         call cd_steepening_r (xnnu(-3,n_s-4,o_s-4,i_nuc),
     &        xnu_if(-1,n_s-2,o_s-2,0,i_nuc))
#endif
      end do

      call interpolate_r (p,p_if, is_positive=.true.)
      call interpolate_r (gamm,gamm_if, is_positive=.true.)

      call flatten_r
         
      call monotonize_r (rho, rho_if)
      call monotonize_r (v_1, v_1_if)
      call monotonize_r (v_2, v_2_if)
      call monotonize_r (v_3, v_3_if)
      call monotonize_r (eps, eps_if)
      do i_nuc=1,config%qn
c         call monotonize_r (xnnu(:,:,:,i_nuc),
c     &        xnu_if(:,:,:,:,i_nuc))
         call monotonize_r (xnnu(-3,n_s-4,o_s-4,i_nuc),
     &        xnu_if(-1,n_s-2,o_s-2,0,i_nuc))
      end do
      call flatten_xnu_r

      call monotonize_r (p,p_if)
      call monotonize_r (gamm,gamm_if)

#ifdef CFC_MHD
      call interpolate_r (b_1, b_1_if)
      call interpolate_r (b_2, b_2_if)
      call interpolate_r (b_3, b_3_if)

#ifndef CONSTRAINED
      call monotonize_r (b_1, b_1_if)
#endif      
      call monotonize_r (b_2, b_2_if)
      call monotonize_r (b_3, b_3_if)

      call interpolate_r (e_1, e_1_if)
      call interpolate_r (e_2, e_2_if)
      call interpolate_r (e_3, e_3_if)

      call monotonize_r (e_1, e_1_if)
      call monotonize_r (e_2, e_2_if)
      call monotonize_r (e_3, e_3_if)

      call interpolate_r (psi, psi_if)
      call monotonize_r (psi, psi_if)
#endif /* CFC_MHD */

ccc   set v_r, v_theta, and v_phi explicitely to zero at center
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         
         v_1_if (1, j, k, 0) = 0.0_rk
         v_1_if (0, j, k, 1) = 0.0_rk
         v_2_if (1, j, k, 0) = 0.0_rk
         v_2_if (0, j, k, 1) = 0.0_rk
         v_3_if (1, j, k, 0) = 0.0_rk
         v_3_if (0, j, k, 1) = 0.0_rk
#ifdef ORIGIN
         v_1_if (1, j, k, 1) =
     &        (r_if (1) / r (2)) * v_1 (2, j, k)
         v_1_if (2, j, k, 0) =
     &        (r_if (1) / r (2)) * v_1 (2, j, k)
#endif

#ifdef CFC_MHD
         b_1_if (1, j, k, 0) = 0.0_rk
         b_1_if (0, j, k, 1) = 0.0_rk
         b_2_if (1, j, k, 0) = 0.0_rk
         b_2_if (0, j, k, 1) = 0.0_rk
         b_3_if (1, j, k, 0) = 0.0_rk
         b_3_if (0, j, k, 1) = 0.0_rk
#endif /* CFC_MHD*/

      enddo

      call rescale_velocities (3)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         
         do l = 0, 1
            do i = 0, m

ccc   calculate fluid velocity v^2
               v_squared_if (i + 1 - l, j, k, l) = 
     &              g_up_11_r (i, j, k) *
     &              v_1_if (i + 1 - l, j, k, l) ** 2 +
     &              g_up_22_r (i, j, k) *
     &                 v_2_if (i + 1 - l, j, k, l) ** 2 +
     &              g_up_33_r (i, j, k) *
     &              v_3_if (i + 1 - l, j, k, l) ** 2
            end do
               
            m0=1-l
            m1=m+m0

            do i = m0,m1
#ifdef CFC_TRANSPORT2
               h_if(i,j,k,l) = 1.0_rk+
     &              eps_if(i,j,k,l)+
     &              p_if(i,j,k,l)/rho_if(i,j,k,l)
               c_sound_squared_if(i,j,k,l) =
     &              gamm_if(i,j,k,l)*p_if(i,j,k,l)/
     &              (rho_if(i,j,k,l)*h_if(i,j,k,l))
#else
               h_if(i,j,k,l) =
     &              eps_if(i,j,k,l)+
     &              p_if(i,j,k,l)/rho_if(i,j,k,l)+
     &              0.5_rk*v_squared_if(i,j,k,l)
               c_sound_squared_if(i,j,k,l)=
     &              gamm_if(i,j,k,l)*p_if(i,j,k,l)/
     &              rho_if(i,j,k,l)
#endif
            enddo
         enddo
      enddo

      return

      END SUBROUTINE reconstruct_interface_r

c     ==================================================================




c     ==================================================================

      SUBROUTINE reconstruct_interface_theta

c     ==================================================================
c     reconstructs interface values in theta direction
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE grid_cfc
      USE parameters_cfc
      USE hydro_primitives_cfc      
      USE metric_cfc      
      USE interface_hydro_cfc
      USE nutra_cfc

      USE gr_ppm
      USE gr_rescale

      IMPLICIT NONE

      integer (kind=ik) :: i, j, jj, k, jk, l, i_nuc, k1, k2

#ifdef AXIS
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         k1 = k
         k2 = k
         if (k. eq. o_s) k1 = o_s - 4
         if (k. eq. o_e) k2 = o_e + 4
         if (j. eq. 1) v_2 (:, 1, k1:k2) =
     &        0.333333_rk * v_2 (:,     2, k1:k2)
         if (j. eq. n) v_2 (:, n, k1:k2) =
     &        0.333333_rk * v_2 (:, n - 1, k1:k2)
         if (j. eq. 1) v_3 (:, 1, k1:k2) =
     &        0.333333_rk * v_3 (:,     2, k1:k2)
         if (j. eq. n) v_3 (:, n, k1:k2) =
     &        0.333333_rk * v_3 (:, n - 1, k1:k2)
      end do
#endif /* AXIS */

      call interpolate_theta (rho, rho_if, is_positive=.true.)
#ifdef CD_STEEPENING
      call cd_steepening_theta (rho, rho_if)
#endif
      
      call interpolate_theta (v_1, v_1_if)
      call interpolate_theta (v_2, v_2_if)
      call interpolate_theta (v_3, v_3_if)
      call interpolate_theta (eps, eps_if)!, is_positive=.true.)

      do i_nuc=1,config%qn
c         call interpolate_theta (xnnu(:,:,:,i_nuc),
c     &        xnu_if(:,:,:,:,i_nuc))
         call interpolate_theta (xnnu(-3,n_s-4,o_s-4,i_nuc),
     &        xnu_if(-1,n_s-2,o_s-2,0,i_nuc),is_mass_fraction=.true.)
#ifdef CD_STEEPENING
c         call cd_steepening_theta (xnnu(:,:,:,i_nuc),
c     &        xnu_if(:,:,:,:,i_nuc))
         call cd_steepening_theta (xnnu(-3,n_s-4,o_s-4,i_nuc),
     &        xnu_if(-1,n_s-2,o_s-2,0,i_nuc))
#endif
      end do

      call interpolate_theta (p,p_if, is_positive=.true.)
      call interpolate_theta (gamm,gamm_if, is_positive=.true.)

      call flatten_theta
      
      call monotonize_theta (rho, rho_if)
      call monotonize_theta (v_1, v_1_if)
      call monotonize_theta (v_2, v_2_if)
      call monotonize_theta (v_3, v_3_if)
      call monotonize_theta (eps, eps_if)
      do i_nuc=1,config%qn
c         call monotonize_theta (xnnu(:,:,:,i_nuc),
c     &        xnu_if(:,:,:,:,i_nuc))
         call monotonize_theta (xnnu(-3,n_s-4,o_s-4,i_nuc),
     &        xnu_if(-1,n_s-2,o_s-2,0,i_nuc))
      end do

      call flatten_xnu_theta

      call monotonize_theta (p,p_if)
      call monotonize_theta (gamm,gamm_if)

#ifdef CFC_MHD
      call interpolate_theta (b_1, b_1_if)
      call interpolate_theta (b_2, b_2_if)
      call interpolate_theta (b_3, b_3_if)

      call monotonize_theta (b_1, b_1_if)
#ifndef CONSTRAINED
      call monotonize_theta (b_2, b_2_if)
#endif
      call monotonize_theta (b_3, b_3_if)

      call interpolate_theta (e_1, e_1_if)
      call interpolate_theta (e_2, e_2_if)
      call interpolate_theta (e_3, e_3_if)

      call monotonize_theta (e_1, e_1_if)
      call monotonize_theta (e_2, e_2_if)
      call monotonize_theta (e_3, e_3_if)

      call interpolate_theta (psi, psi_if)
      call monotonize_theta (psi, psi_if)
#endif /* CFC_MHD */

ccc   set v_theta explicitely to zero at pole and equator
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         if (j.eq.1) then
            do i = 1, m
               v_2_if (i, 1, k, 0) = 0.0_rk
               v_2_if (i, 0, k, 1) = 0.0_rk
               v_3_if (i, 1, k, 0) = 0.0_rk
               v_3_if (i, 0, k, 1) = 0.0_rk
#ifdef CFC_MHD
               b_2_if (i, 1, k, 0) = 0.0_rk
               b_2_if (i, 0, k, 1) = 0.0_rk
               b_3_if (i, 1, k, 0) = 0.0_rk
               b_3_if (i, 0, k, 1) = 0.0_rk
#endif /* CFC_MHD */
            enddo
         endif
         if (j.eq.n) then
            do i = 1, m
               v_2_if (i, n + 1, k, 0) = 0.0_rk
               v_2_if (i, n, k, 1) = 0.0_rk
            enddo
         endif
         if (j.eq.n .and. .not. equatorial_symmetry) then
            do i = 1, m
               v_3_if (i, n + 1, k, 0) = 0.0_rk
               v_3_if (i, n, k, 1) = 0.0_rk
#ifdef CFC_MHD
c B-field does not need to be zero at equator since it is
c a pseudovector
               b_2_if (i, n + 1, k, 0) = 0.0_rk
               b_2_if (i, n, k, 1) = 0.0_rk
               b_3_if (i, n + 1, k, 0) = 0.0_rk
               b_3_if (i, n, k, 1) = 0.0_rk
#endif /* CFC_MHD */
            end do
         end if
      enddo

      call rescale_velocities (4)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO      
#endif
      do jk = 1, (n_loc + 1) * o_loc
         
         k = int((jk + (n_loc + 1) - 1) / (n_loc + 1))
         j = (n_s - 2) + (jk - (k - 1) * (n_loc + 1))
         k = k + o_s - 1

         do l = 0, 1
            jj=j+1-l
            do i = 1, m
ccc   calculate fluid velocity v^2
               v_squared_if (i, jj, k, l) = 
     &              g_up_11_theta (i, j, k) *
     &              v_1_if (i, jj, k, l) ** 2 +
     &              g_up_22_theta (i, j, k) *
     &              v_2_if (i, jj, k, l) ** 2 +
     &              g_up_33_theta (i, j, k) *
     &              v_3_if (i, jj, k, l) ** 2
            end do
            
            do i=1,m
#ifdef CFC_TRANSPORT2
               h_if(i,jj,k,l) = 1.0_rk+
     &              eps_if(i,jj,k,l)+
     &              p_if(i,jj,k,l)/rho_if(i,jj,k,l)
               c_sound_squared_if(i,jj,k,l) =
     &              gamm_if(i,jj,k,l)*p_if(i,jj,k,l)/
     &              (rho_if(i,jj,k,l)*h_if(i,jj,k,l))
#else 
               h_if(i,jj,k,l) =
     &              eps_if(i,jj,k,l)+
     &              p_if(i,jj,k,l)/rho_if(i,jj,k,l)+
     &              0.5_rk*v_squared_if(i,jj,k,l)
               c_sound_squared_if(i,jj,k,l)=
     &              gamm_if(i,jj,k,l)*p_if(i,jj,k,l)/
     &              rho_if(i,jj,k,l)
#endif
            enddo
         enddo
      enddo

      return

      END SUBROUTINE reconstruct_interface_theta

c     ==================================================================


c     ==================================================================

      SUBROUTINE reconstruct_interface_varphi

c     ==================================================================
c     reconstructs interface values in varpi direction
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE grid_cfc
      USE parameters_cfc
      USE hydro_primitives_cfc      
      USE metric_cfc      
      USE interface_hydro_cfc

      USE gr_ppm
      USE gr_rescale

      use configure

      real (kind=rk) :: den_tmp(m),ei_tmp(m),tem_tmp(m),
     &     pre_tmp(m),s_tmp(m),gam_tmp(m),xnu_tmp(m,config%qn),
     &     kap_tmp(m)
      real (kind=rk) ::  ccu(m), cce(m), ccn(m), ccp(m),
     &     dmy1(m), dmy2(m,2) !dummys

      integer (kind=ik) :: i, j, k, kk, jk, l, i_nuc

      call interpolate_varphi (rho, rho_if, is_positive=.true.)
#ifdef CD_STEEPENING
      call cd_steepening_varphi (rho, rho_if)
#endif
      
      call interpolate_varphi (v_1, v_1_if)
      call interpolate_varphi (v_2, v_2_if)
      call interpolate_varphi (v_3, v_3_if)
      call interpolate_varphi (eps, eps_if)!, is_positive=.true.)

      do i_nuc=1,config%qn
c         call interpolate_varphi (xnnu(:,:,:,i_nuc),
c     &        xnu_if(:,:,:,:,i_nuc))
         call interpolate_varphi (xnnu(-3,n_s-4,o_s-4,i_nuc),
     &        xnu_if(-1,n_s-2,o_s-2,0,i_nuc),is_mass_fraction=.true.)
#ifdef CD_STEEPENING
c         call cd_steepening_varphi (xnnu(:,:,:,i_nuc),
c     &        xnu_if(:,:,:,:,i_nuc))
         call cd_steepening_varphi (xnnu(-3,n_s-4,o_s-4,i_nuc),
     &        xnu_if(-1,n_s-2,o_s-2,0,i_nuc))
#endif         
      end do

      call interpolate_varphi (p,p_if, is_positive=.true.)
      call interpolate_varphi (gamm,gamm_if, is_positive=.true.)


      call flatten_varphi
      
      call monotonize_varphi (rho, rho_if)
      call monotonize_varphi (v_1, v_1_if)
      call monotonize_varphi (v_2, v_2_if)
      call monotonize_varphi (v_3, v_3_if)
      call monotonize_varphi (eps, eps_if)
      do i_nuc=1,config%qn
c         call monotonize_varphi (xnnu(:,:,:,i_nuc),
c     &        xnu_if(:,:,:,:,i_nuc))
         call monotonize_varphi (xnnu(-3,n_s-4,o_s-4,i_nuc),
     &        xnu_if(-1,n_s-2,o_s-2,0,i_nuc))
      end do

      call flatten_xnu_varphi
      call monotonize_varphi (p,p_if)
      call monotonize_varphi (gamm,gamm_if)

#ifdef CFC_MHD
      call interpolate_varphi (b_1, b_1_if)
      call interpolate_varphi (b_2, b_2_if)
      call interpolate_varphi (b_3, b_3_if)

      call monotonize_varphi (b_1, b_1_if)
      call monotonize_varphi (b_2, b_2_if)
#ifndef CONSTRAINED
      call monotonize_varphi (b_3, b_3_if)
#endif
      
      call interpolate_varphi (e_1, e_1_if)
      call interpolate_varphi (e_2, e_2_if)
      call interpolate_varphi (e_3, e_3_if)

      call monotonize_varphi (e_1, e_1_if)
      call monotonize_varphi (e_2, e_2_if)
      call monotonize_varphi (e_3, e_3_if)

      call interpolate_varphi (psi, psi_if)
      call monotonize_varphi (psi, psi_if)
#endif /* CFC_MHD */
      
      call rescale_velocities (5)


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO      
#endif
      do jk = 1, n_loc * (o_loc + 1)

         k = int((jk + n_loc - 1) / n_loc)
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 2

         do l = 0, 1
            
            kk=k+1-l

            do i = 1, m
ccc   calculate fluid velocity v^2
               v_squared_if (i, j, kk, l) = 
     &              g_up_11_varphi (i, j, k) *
     &              v_1_if (i, j, kk, l) ** 2 +
     &              g_up_22_varphi (i, j, k) *
     &              v_2_if (i, j, kk, l) ** 2 +
     &              g_up_33_varphi (i, j, k) *
     &                 v_3_if (i, j, kk, l) ** 2
            end do
            
            do i=1,m
#ifdef CFC_TRANSPORT2
               h_if(i,j,kk,l) = 1.0_rk+
     &              eps_if(i,j,kk,l)+
     &              p_if(i,j,kk,l)/rho_if(i,j,kk,l)
               c_sound_squared_if(i,j,kk,l) =
     &              gamm_if(i,j,kk,l)*p_if(i,j,kk,l)/
     &              (rho_if(i,j,kk,l)*h_if(i,j,kk,l))
#else 
               h_if(i,j,kk,l) =
     &              eps_if(i,j,kk,l)+
     &              p_if(i,j,kk,l)/rho_if(i,j,kk,l)+
     &              0.5_rk*v_squared_if(i,j,kk,l)
               c_sound_squared_if(i,j,kk,l)=
     &              gamm_if(i,j,kk,l)*p_if(i,j,kk,l)/
     &              rho_if(i,j,kk,l)
#endif
            enddo
         enddo
      enddo

      return

      END SUBROUTINE reconstruct_interface_varphi

c     ==================================================================


      END MODULE gr_reconstruction

c     ==================================================================
