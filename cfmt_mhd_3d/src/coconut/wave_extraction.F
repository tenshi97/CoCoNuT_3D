c     ============================================================

      MODULE gw_extraction

c     ============================================================

      IMPLICIT NONE

      CONTAINS



c     ============================================================

      SUBROUTINE extract_gravitational_waves

c     ============================================================
c
c     WARNING: note that only the amplitude of the plus-polarized
c     waves is computed in 3D!
c
c     ------------------------------------------------------------

      USE precision

      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
      USE hydro_primitives_cfc
      USE conserved_cfc
      USE metric_cfc
      USE perm_aux_cfc
      USE temp_aux_cfc
      USE boundary_cfc
      USE grav_wave_cfc

#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      IMPLICIT NONE

      logical, parameter :: use_rho_star = .true.,
     &     use_stress_formula = .true., cowling = .false.

      logical, parameter :: newtonian = .false.

      real (kind=rk) :: delta_t_tmp, delta_t_tmp_old, dummy

      integer (kind=ik) :: i, j, k, jk, l
#ifdef MPI_HYDRO
      integer (kind=ik) :: ierr
#endif

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif

      if (iteration.eq.1) then
         
         di_cap_dt_2 = 0.0_rk

         di_cap_11_dt_2 = 0.0_rk
         di_cap_22_dt_2 = 0.0_rk
         di_cap_33_dt_2 = 0.0_rk
         di_cap_12_dt_2 = 0.0_rk
         di_cap_13_dt_2 = 0.0_rk
         di_cap_23_dt_2 = 0.0_rk

         ddi_cap_ddt_1 = 0.0_rk

         ddi_cap_11_ddt_1 = 0.0_rk
         ddi_cap_22_ddt_1 = 0.0_rk
         ddi_cap_33_ddt_1 = 0.0_rk
         ddi_cap_12_ddt_1 = 0.0_rk
         ddi_cap_13_ddt_1 = 0.0_rk
         ddi_cap_23_ddt_1 = 0.0_rk

         ddi_cap_ddt_2 = 0.0_rk

         ddi_cap_11_ddt_2 = 0.0_rk
         ddi_cap_22_ddt_2 = 0.0_rk
         ddi_cap_33_ddt_2 = 0.0_rk
         ddi_cap_12_ddt_2 = 0.0_rk
         ddi_cap_13_ddt_2 = 0.0_rk
         ddi_cap_23_ddt_2 = 0.0_rk

         ddi_cap_ddt_1_old = 0.0_rk

         ddi_cap_11_ddt_1_old = 0.0_rk
         ddi_cap_22_ddt_1_old = 0.0_rk
         ddi_cap_33_ddt_1_old = 0.0_rk
         ddi_cap_12_ddt_1_old = 0.0_rk
         ddi_cap_13_ddt_1_old = 0.0_rk
         ddi_cap_23_ddt_1_old = 0.0_rk

         ddi_cap_ddt_2_old_old = 0.0_rk

         ddi_cap_11_ddt_2_old_old = 0.0_rk
         ddi_cap_22_ddt_2_old_old = 0.0_rk
         ddi_cap_33_ddt_2_old_old = 0.0_rk
         ddi_cap_12_ddt_2_old_old = 0.0_rk
         ddi_cap_13_ddt_2_old_old = 0.0_rk
         ddi_cap_23_ddt_2_old_old = 0.0_rk

         ddi_cap_ddt_2_old = 0.0_rk

         ddi_cap_11_ddt_2_old = 0.0_rk
         ddi_cap_22_ddt_2_old = 0.0_rk
         ddi_cap_33_ddt_2_old = 0.0_rk
         ddi_cap_12_ddt_2_old = 0.0_rk
         ddi_cap_13_ddt_2_old = 0.0_rk
         ddi_cap_23_ddt_2_old = 0.0_rk

         dddi_cap_dddt_1 = 0.0_rk

         dddi_cap_11_dddt_1 = 0.0_rk
         dddi_cap_22_dddt_1 = 0.0_rk
         dddi_cap_33_dddt_1 = 0.0_rk
         dddi_cap_12_dddt_1 = 0.0_rk
         dddi_cap_13_dddt_1 = 0.0_rk
         dddi_cap_23_dddt_1 = 0.0_rk

         dddi_cap_dddt_2 = 0.0_rk

         dddi_cap_11_dddt_2 = 0.0_rk
         dddi_cap_22_dddt_2 = 0.0_rk
         dddi_cap_33_dddt_2 = 0.0_rk
         dddi_cap_12_dddt_2 = 0.0_rk
         dddi_cap_13_dddt_2 = 0.0_rk
         dddi_cap_23_dddt_2 = 0.0_rk

      endif

      di_cap_dt_2_old = di_cap_dt_2

      di_cap_11_dt_2_old = di_cap_11_dt_2
      di_cap_22_dt_2_old = di_cap_22_dt_2
      di_cap_33_dt_2_old = di_cap_33_dt_2
      di_cap_12_dt_2_old = di_cap_12_dt_2
      di_cap_13_dt_2_old = di_cap_13_dt_2
      di_cap_23_dt_2_old = di_cap_23_dt_2

      ddi_cap_ddt_1_old_old = ddi_cap_ddt_1_old

      ddi_cap_11_ddt_1_old_old = ddi_cap_11_ddt_1_old
      ddi_cap_22_ddt_1_old_old = ddi_cap_22_ddt_1_old
      ddi_cap_33_ddt_1_old_old = ddi_cap_33_ddt_1_old
      ddi_cap_12_ddt_1_old_old = ddi_cap_12_ddt_1_old
      ddi_cap_13_ddt_1_old_old = ddi_cap_13_ddt_1_old
      ddi_cap_23_ddt_1_old_old = ddi_cap_23_ddt_1_old

      ddi_cap_ddt_2_old_old_old = ddi_cap_ddt_2_old_old

      ddi_cap_11_ddt_2_old_old_old = ddi_cap_11_ddt_2_old_old
      ddi_cap_22_ddt_2_old_old_old = ddi_cap_22_ddt_2_old_old
      ddi_cap_33_ddt_2_old_old_old = ddi_cap_33_ddt_2_old_old
      ddi_cap_12_ddt_2_old_old_old = ddi_cap_12_ddt_2_old_old
      ddi_cap_13_ddt_2_old_old_old = ddi_cap_13_ddt_2_old_old
      ddi_cap_23_ddt_2_old_old_old = ddi_cap_23_ddt_2_old_old

      ddi_cap_ddt_2_old_old = ddi_cap_ddt_2_old

      ddi_cap_11_ddt_2_old_old = ddi_cap_11_ddt_2_old
      ddi_cap_22_ddt_2_old_old = ddi_cap_22_ddt_2_old
      ddi_cap_33_ddt_2_old_old = ddi_cap_33_ddt_2_old
      ddi_cap_12_ddt_2_old_old = ddi_cap_12_ddt_2_old
      ddi_cap_13_ddt_2_old_old = ddi_cap_13_ddt_2_old
      ddi_cap_23_ddt_2_old_old = ddi_cap_23_ddt_2_old

      ddi_cap_ddt_1_old = ddi_cap_ddt_1

      ddi_cap_11_ddt_1_old = ddi_cap_11_ddt_1
      ddi_cap_22_ddt_1_old = ddi_cap_22_ddt_1
      ddi_cap_33_ddt_1_old = ddi_cap_33_ddt_1
      ddi_cap_12_ddt_1_old = ddi_cap_12_ddt_1
      ddi_cap_23_ddt_1_old = ddi_cap_23_ddt_1
      ddi_cap_33_ddt_1_old = ddi_cap_33_ddt_1

      ddi_cap_ddt_2_old = ddi_cap_ddt_2
      
      ddi_cap_11_ddt_2_old = ddi_cap_11_ddt_2
      ddi_cap_22_ddt_2_old = ddi_cap_22_ddt_2
      ddi_cap_33_ddt_2_old = ddi_cap_33_ddt_2
      ddi_cap_12_ddt_2_old = ddi_cap_12_ddt_2
      ddi_cap_13_ddt_2_old = ddi_cap_13_ddt_2
      ddi_cap_23_ddt_2_old = ddi_cap_23_ddt_2

      dddi_cap_dddt_1_old = dddi_cap_dddt_1

      dddi_cap_11_dddt_1_old = dddi_cap_11_dddt_1
      dddi_cap_22_dddt_1_old = dddi_cap_22_dddt_1
      dddi_cap_33_dddt_1_old = dddi_cap_33_dddt_1
      dddi_cap_12_dddt_1_old = dddi_cap_12_dddt_1
      dddi_cap_13_dddt_1_old = dddi_cap_13_dddt_1
      dddi_cap_23_dddt_1_old = dddi_cap_23_dddt_1

      dddi_cap_dddt_2_old = dddi_cap_dddt_2

      dddi_cap_11_dddt_2_old = dddi_cap_11_dddt_2
      dddi_cap_22_dddt_2_old = dddi_cap_22_dddt_2
      dddi_cap_33_dddt_2_old = dddi_cap_33_dddt_2
      dddi_cap_12_dddt_2_old = dddi_cap_12_dddt_2
      dddi_cap_13_dddt_2_old = dddi_cap_13_dddt_2
      dddi_cap_23_dddt_2_old = dddi_cap_23_dddt_2

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE NOWAIT
#endif




c     calculate Newtonian potential and its spatial derivatives
      if (use_stress_formula) then
         
c     calculate the potential by solving the Poisson equation
c     with expansion into Legendre polynoms
         
c     calculate the potential Phi = 1/2 (alpha phi^6 r^2 sin theta) from metric functions
#if defined(OPENMP_CFC) && ( defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc

           k = int((jk + n_loc - 1) / n_loc )
           j = (n_s - 1) + (jk - (k - 1) * n_loc)
           k = k + o_s - 1

           do i = 1, m
c                  dphi_potential_alt_dr (i, j, k) =
c     &                 - 0.5_rk * alpha (i, j, k) * phi (i, j, k) ** 6 *
c     &                 ((alpha_r (i, j, k) -
c     &                 alpha_r (i - 1, j, k)) /
c     &                 alpha (i, j, k) +
c     &                 6.0_rk * (phi_r (i, j, k) -
c     &                 phi_r (i - 1, j, k)) /
c     &                 phi (i, j, k)) *
c     &                 delta_r_inv (i)
c                  dphi_potential_alt_dtheta (i, j, k) =
c     &                 - 0.5_rk * alpha (i, j, k) * phi (i, j, k) ** 6 *
c     &                 ((alpha_theta (i, j, k) -
c     &                 alpha_theta (i, j - 1, k)) /
c     &                 alpha (i, j, k) +
c     &                 6.0_rk * (phi_theta (i, j, k) -
c     &                 phi_theta (i, j - 1, k)) /
c     &                 phi (i, j, k)) *
c     &                 delta_theta_inv
c                  dphi_potential_3_dvarphi (i, j, k) =
c     &                 - 0.5_rk * alpha (i, j, k) * phi (i, j, k) ** 6 *
c     &                 ((alpha_varphi (i, j, k) -
c     &                 alpha_varphi (i, j, k - 1)) /
c     &                 alpha (i, j, k) +
c     &                 6.0_rk * (phi_varphi (i, j, k) -
c     &                 phi_varphi (i, j, k - 1)) /
c     &                 phi (i, j, k)) *
c     &                 delta_varphi_inv
                  
              dphi_potential_alt_dr (i, j, k) = 
     &             - 2.0_rk * phi (i, j, k) ** 3 *
     &             (phi_r (i, j, k) - phi_r (i - 1, j, k)) *
     &             delta_r_inv (i)
              dphi_potential_alt_dtheta (i, j, k) =
     &             - 2.0_rk * phi (i, j, k) ** 3 *
     &             (phi_theta (i, j, k) - phi_theta (i, j - 1, k)) *
     &          delta_theta_inv
              dphi_potential_alt_dvarphi (i, j, k) =
     &             - 2.0_rk * phi (i, j, k) ** 3 *
     &             (phi_varphi (i, j, k) -
     &             phi_varphi (i, j, k - 1)) *
     &             delta_varphi_inv
              
           enddo
        enddo
         
      endif

c     calculate physical velocities
#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m
            
            v_gw_r (i, j, k) =
     &           v_1 (i, j, k) * sqrt (g_up_11 (i, j, k))
            v_gw_theta (i, j, k) =
     &           v_2 (i, j, k) * sqrt (g_up_22 (i, j, k))
            v_gw_varphi (i, j, k) =
     &           v_3 (i, j, k) * sqrt (g_up_33 (i, j, k))
            
         enddo
      enddo



c     calculate rho_gw = rho or rho*, depending on flag
      if (use_rho_star) then
#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            do i = 1, m
               rho_gw (i, j, k) =
     &              rho (i, j, k) * phi (i, j, k) ** 6 * w (i, j, k)
            enddo
         enddo
      else
#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            do i = 1, m
               rho_gw (i, j, k) = rho (i, j, k)
            enddo
         enddo
      endif



c     sum up quadrupole moments

#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D)    )
!$OMP SINGLE
#endif
         
      di_cap_dt_2 = 0.0_rk
      
#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D)   )
!$OMP END SINGLE
#endif
         
#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO reduction (+: di_cap_dt_2)
#endif
      do jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m
            
            di_cap_dt_2 = di_cap_dt_2 +
     &           rrr (i) * sinus_theta (j) *
     &           delta_r (i) * delta_theta * delta_varphi *
     &           rho_gw (i, j, k) * (v_gw_r (i, j, k) *
     &           (3.0_rk * cosin_theta (j) ** 2 - 1.0_rk) -
     &           3.0_rk * v_gw_theta (i, j, k) *
     &           sinuscosin_theta (j))
                  
         enddo
      enddo
      
#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D)      )
!$OMP SINGLE
#endif
      
#ifdef MPI_HYDRO
c     Only process 0 needs di_cap_dt_2 -> MPI_Reduce is enough
      call MPI_Reduce(di_cap_dt_2,dummy,1_ik,
     $     MPI_DOUBLE_PRECISION,MPI_SUM,0_ik,cart_comm,ierr)
      di_cap_dt_2=dummy
#endif      

      ddi_cap_ddt_2 = (di_cap_dt_2 - di_cap_dt_2_old) / delta_t
      
      a_quad_plus_equator_2 = gw_factor_1 * ddi_cap_ddt_2
      
      delta_t_tmp = 0.5_rk * (delta_t + delta_t_old) + 1e-80_rk
      delta_t_tmp_old = 0.5_rk * (delta_t_old + delta_t_old_old) + 
     &     1e-80_rk
      
      dddi_cap_dddt_2 =
     &     ((ddi_cap_ddt_2 - ddi_cap_ddt_2_old) *
     &     (delta_t_tmp_old / delta_t_tmp) +
     &     (ddi_cap_ddt_2_old - ddi_cap_ddt_2_old_old) *
     &     (delta_t_tmp / delta_t_tmp_old)) /
     &     (delta_t_tmp + delta_t_tmp_old)
c     simpler form
c         dddi_cap_dddt_2 =
c     &        (ddi_cap_ddt_2 - ddi_cap_ddt_2_old) / delta_t
         
      if (iteration.lt.3) then
         e_grav_rad_2 = 0.0_rk
      else
         e_grav_rad_2 = e_grav_rad_2 +
     &        0.5_rk * gw_factor_1 ** 2 * gw_factor_2 *
     &        delta_t_tmp_old *
     &        (dddi_cap_dddt_2 ** 2 + dddi_cap_dddt_2_old ** 2)
c     simpler form
c            e_grav_rad_2 = e_grav_rad_2 +
c     &           gw_factor_1 ** 2 * gw_factor_2 * delta_t *
c     &           dddi_cap_dddt_2 ** 2
      endif

      a_quad_plus_equator_out = a_quad_plus_equator_2
      a_quad_plus_pole_out    = a_quad_plus_pole_2
      a_quad_cross_pole_out   = a_quad_cross_pole_2

#if  defined(OPENMP_CFC) &&  (defined(OPEN_MP_2D) || defined(OPEN_MP_3D)     )
!$OMP END SINGLE
#endif         
      

      return

      END SUBROUTINE extract_gravitational_waves

c     ============================================================



      END MODULE gw_extraction

c     ============================================================
