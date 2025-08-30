c     ==================================================================

      MODULE gr_tstep

c     ==================================================================

      use precision

      IMPLICIT NONE
      
      real (kind=rk), allocatable, save :: lambda_p_r (:,:,:)
      real (kind=rk), allocatable, save :: lambda_m_r (:,:,:)
      real (kind=rk), allocatable, save :: lambda_p_theta (:,:,:)
      real (kind=rk), allocatable, save :: lambda_m_theta (:,:,:)
      real (kind=rk), allocatable, save :: lambda_p_varphi (:,:,:)
      real (kind=rk), allocatable, save :: lambda_m_varphi (:,:,:)


      CONTAINS


c     ==================================================================
      
      SUBROUTINE allocate_gr_tstep

c     ==================================================================

      USE size_cfc

      IMPLICIT NONE

      allocate(lambda_p_r (1 : m, n_s : n_e, o_s : o_e), 
     &         lambda_m_r (1 : m, n_s : n_e, o_s : o_e), 
     &         lambda_p_theta (1 : m, n_s : n_e, o_s : o_e), 
     &         lambda_m_theta (1 : m, n_s : n_e, o_s : o_e),
     &         lambda_p_varphi (1 : m, n_s : n_e, o_s : o_e),
     &         lambda_m_varphi (1 : m, n_s : n_e, o_s : o_e) )

      END SUBROUTINE allocate_gr_tstep

c     ==================================================================



c     ==================================================================

      SUBROUTINE compute_time_step

c     ==================================================================
c      
c     This subroutine computes the time step from the CFL criterion
c
c     ------------------------------------------------------------------

      USE precision
      
      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
      USE hydro_primitives_cfc
      USE conserved_cfc
      USE metric_cfc
      USE perm_aux_cfc


#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      use hydro_areas_mod
      use configure

      IMPLICIT NONE

      real (kind=rk) :: v_up_1, v_up_2, v_up_3

      real (kind=rk) :: c_sound
      
      real (kind=rk) :: delta_t_temp_inv

      real (kind=rk) :: tmp_1, tmp_2, tmp_3, dummy

      real (kind=rk) :: cfac, rspher
      integer (kind=ik) :: clevel

      integer (kind=ik) :: i, j, k, jk, iarea,ixi,ixf,ioy


#ifdef MPI_HYDRO
      integer (kind=ik) :: ierr
#endif

      ixf    = areas%ix_are(1, 2)
      ioy    = areas%ix_are(1, 6)
      if (ioy.eq.n) then
         rspher = r (ixf + 1)
      else
         rspher = 1e6
      end if


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif

ccc   store old time step
      delta_t_old_old = delta_t_old
      delta_t_old = delta_t

ccc   initialization of some quantities
      delta_t_r_inv = 0.0_rk 
      delta_t_theta_inv = 0.0_rk
      delta_t_varphi_inv = 0.0_rk
      
      delta_t = 1.0e+100_rk
      
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE NOWAIT
#endif
         

ccc   calculate eigenvalues lambda_+ and lambda_-
ccc   for the characteristic speeds:
#ifdef CFC_TRANSPORT2
c     Relativistic case:
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         if (o. gt. 1) then
            clevel = INT (- LOG (sinus_theta (j) * delta_varphi /
     &           delta_theta) / LOG (2.0_rk))
         else
            clevel = 0
         end if
         cfac = 1.0_rk / 2.0_rk ** clevel

         do i = 1, m
            
            v_up_1 = g_up_11 (i, j, k) * v_1 (i, j, k)
            v_up_2 = g_up_22 (i, j, k) * v_2 (i, j, k)
            v_up_3 = g_up_33 (i, j, k) * v_3 (i, j, k)

            c_sound = sqrt (c_sound_squared (i, j, k))
            
            tmp_1 =
     &           v_up_1 * (1.0_rk - c_sound_squared (i, j, k))
            tmp_2 = (1.0_rk - v_squared (i, j, k)) *
     &           (g_up_11 (i, j, k) * (1.0_rk -
     &           c_sound_squared (i, j, k) *
     &           v_squared (i, j, k)) - v_up_1 ** 2 * (1.0_rk -
     &           c_sound_squared (i, j, k)))
            tmp_3 = alpha (i, j, k) / (1.0_rk -
     &           c_sound_squared (i, j, k) *
     &           v_squared (i, j, k))
            
            lambda_p_r (i, j, k) =
     &           tmp_3 * (tmp_1 + c_sound * 
     &           sqrt (tmp_2)) - beta_up_1 (i, j, k)
            lambda_m_r (i, j, k) =
     &           tmp_3 * (tmp_1 - c_sound * 
     &           sqrt (tmp_2)) - beta_up_1 (i, j, k)
            
            tmp_1 =
     &           v_up_2 * (1.0_rk - c_sound_squared (i, j, k))
            tmp_2 = (1.0_rk - v_squared (i, j, k)) *
     &           (g_up_22 (i, j, k) * (1.0_rk -
     &           c_sound_squared (i, j, k) *
     &           v_squared (i, j, k)) - v_up_2 ** 2 * (1.0_rk -
     &           c_sound_squared (i, j, k)))
            tmp_3 = alpha (i, j, k) / (1.0_rk -
     &           c_sound_squared (i, j, k) *
     &           v_squared (i, j, k))
            
            lambda_p_theta (i, j, k) =
     &           tmp_3 * (tmp_1 + c_sound * 
     &           sqrt (tmp_2)) - beta_up_2 (i, j, k)
            lambda_m_theta (i, j, k) =
     &           tmp_3 * (tmp_1 - c_sound * 
     &           sqrt (tmp_2)) - beta_up_2 (i, j, k)
            
            tmp_1 =
     &           v_up_3 * (1.0_rk - c_sound_squared (i, j, k))
            tmp_2 = (1.0_rk - v_squared (i, j, k)) *
     &           (g_up_33 (i, j, k) * (1.0_rk -
     &           c_sound_squared (i, j, k) *
     &           v_squared (i, j, k)) - v_up_3 ** 2 * (1.0_rk -
     &           c_sound_squared (i, j, k)))
            tmp_3 = alpha (i, j, k) / (1.0_rk -
     &           c_sound_squared (i, j, k) *
     &           v_squared (i, j, k))
 
            lambda_p_varphi (i, j, k) =
     &           tmp_3 * (tmp_1 + c_sound * 
     &           sqrt (tmp_2)) - beta_up_3 (i, j, k)
            lambda_m_varphi (i, j, k) =
     &           tmp_3 * (tmp_1 - c_sound * 
     &           sqrt (tmp_2)) - beta_up_3 (i, j, k)
#ifdef COARSENING
c     mesh coarsening
            lambda_p_varphi (i, j, k) = lambda_p_varphi (i, j, k) * 
     &           cfac
            lambda_m_varphi (i, j, k) = lambda_m_varphi (i, j, k) *
     &           cfac
#else /* COARSENING */
c     filtering
            lambda_p_varphi (i, j, k) = lambda_p_varphi (i, j, k) * 
     &           MIN(1.0_rk, 2.0_rk*sinus_theta (j) * r(i) / rspher)
            lambda_m_varphi (i, j, k) = lambda_m_varphi (i, j, k) *
     &           MIN(1.0_rk, 2.0_rk*sinus_theta (j) * r(i) / rspher)
#endif
           
         enddo
      enddo
#else /* CFC_TRANSPORT2 */
c     Newtonian case:
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m
            
            v_up_1 = g_up_11 (i, j, k) * v_1 (i, j, k)
            v_up_2 = g_up_22 (i, j, k) * v_2 (i, j, k)
            v_up_3 = g_up_33 (i, j, k) * v_3 (i, j, k)

#if defined(CFC_MHD) && !defined(PASSIVE_MHD)
            tmp_1 = b_1 (i, j,k) ** 2 + b_2 (i, j, k) ** 2 +
     &           b_3 (i, j, k) ** 2
            c_sound = sqrt (tmp_1 / rho (i, j, k) +
     &           c_sound_squared (i, j, k))
#if DIMN<=3
            c_sound = c_sound + sqrt (v_up_3 * v_3 (i, j, k))
#endif            
#else
            c_sound = sqrt (c_sound_squared (i, j, k))
#endif
            
            lambda_p_r (i, j, k) =
     &           v_up_1 + sqrt (g_up_11 (i, j, k)) * c_sound
            lambda_m_r (i, j, k) =
     &           v_up_1 - sqrt (g_up_11 (i, j, k)) * c_sound
            
            lambda_p_theta (i, j, k) =
     &           v_up_2 + sqrt (g_up_22 (i, j, k)) * c_sound
            lambda_m_theta (i, j, k) =
     &              v_up_2 - sqrt (g_up_22 (i, j, k)) * c_sound
            
            lambda_p_varphi (i, j, k) =
     &           v_up_3 + sqrt (g_up_33 (i, j, k)) * c_sound
            lambda_m_varphi (i, j, k) =
     &           v_up_3 - sqrt (g_up_33 (i, j, k)) * c_sound

#if 1
c     mesh coarsening
            lambda_p_varphi (i, j, k) = lambda_p_varphi (i, j, k) * 
c     &           cfac
     &           MIN(1.0_rk, 2.0_rk*sinus_theta (j) * r(i) / rspher)
            lambda_m_varphi (i, j, k) = lambda_m_varphi (i, j, k) *
c     &           cfac
     &           MIN(1.0_rk, 2.0_rk*sinus_theta (j) * r(i) / rspher)
#endif
            
         enddo
      enddo
#endif /*CFC_TRANSPORT2*/


ccc   calculate dt in r direction
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO REDUCTION (max: delta_t_r_inv)
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m - 1
ccc   do for all gridpoints -- search maximum of
ccc   |lambda_p_r - lambda_m_r| / delta_r
            delta_t_temp_inv = 
     &           (abs (lambda_p_r (i, j, k)) + 
     &           abs (lambda_m_r (i + 1, j, k))) *
     &           delta_r_inv (i)
            delta_t_r_inv =
     &           max (delta_t_temp_inv, delta_t_r_inv)
         enddo
      enddo

c     Regions that are computed in spherical symmetry:
      do iarea=1,areas%are_nu
         ixi    = areas%ix_are(iarea, 1)
         ixf    = areas%ix_are(iarea, 2)
         ioy    = areas%ix_are(iarea, 6)
         if (ioy.eq.n) then
            lambda_p_theta (ixi:ixf,:,:)=0.0_rk
            lambda_m_theta (ixi:ixf,:,:)=0.0_rk
            lambda_p_varphi (ixi:ixf,:,:)=0.0_rk
            lambda_m_varphi (ixi:ixf,:,:)=0.0_rk
         end if
      end do


ccc   calculate dt in theta direction      
      if (config%nsdim .ge. 2) then
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO REDUCTION (max: delta_t_theta_inv)
#endif
      do jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

c      do jk = 1, (n_loc - 1) * o_loc
c         k = int((jk + (n_loc - 1) - 1) / (n_loc - 1))
c         j = (n_s - 1) + (jk - (k - 1) * (n_loc - 1))
c         k = k + o_s - 1

         do i = 1, m
ccc   do for all gridpoints -- search maximum of
ccc   |lambda_p_theta - lambda_m_theta| / delta_theta
            delta_t_temp_inv =
     &           (abs (lambda_p_theta (i, j, k)) +
     &           abs (lambda_m_theta (i, j, k))) *
     &           delta_theta_inv
            delta_t_theta_inv =
     &           max (delta_t_temp_inv, delta_t_theta_inv)
         enddo
      enddo

      if (config%nsdim .eq. 3) then
ccc   calculate dt in varphi direction
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO REDUCTION (max: delta_t_varphi_inv)
#endif
      do jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

c      do jk = 1, n_loc * (o_loc - 1)
c         k = int((jk + n_loc - 1) / n_loc)
c         j = (n_s - 1) + (jk - (k - 1) * n_loc)
c         k = k + o_s - 1

         do i = 1, m
ccc   do for all gridpoints -- search maximum of
ccc   |lambda_p_varphi - lambda_m_varphi| / delta_varphi
            delta_t_temp_inv =
     &           (abs (lambda_p_varphi (i, j, k)) +
     &           abs (lambda_m_varphi (i, j, k))) *
     &           delta_varphi_inv
            delta_t_varphi_inv =
     &           max (delta_t_temp_inv, delta_t_varphi_inv)
         enddo
      enddo
      endif ! config%nsdim = 3
      endif ! config%nsdim .ge. 2


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif

      delta_t = min (delta_t, 1.0_rk / delta_t_r_inv)

      if (config%nsdim .ge. 2) then
         delta_t = min (delta_t, 1.0_rk / delta_t_theta_inv)
      endif
      if (config%nsdim .eq. 3) then
         delta_t = min (delta_t, 1.0_rk / delta_t_varphi_inv)
      endif

c      print *, 'timestep',1.0/delta_t_r_inv,
c     &     1.0/delta_t_theta_inv, 1.0/delta_t_varphi_inv


ccc   new dt = CFL factor * dt
      delta_t = cfl_factor * delta_t

#ifdef MPI_HYDRO
      call MPI_Allreduce(delta_t, dummy, 1_ik,
     &     MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      delta_t=dummy
#endif      
      
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif

      return

      END SUBROUTINE compute_time_step

c     ==================================================================



      END MODULE gr_tstep

c     ==================================================================
