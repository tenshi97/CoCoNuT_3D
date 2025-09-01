#ifdef NOMOTO_MODEL
#define RESCALE_CONSERVED_QUANTITIES
#define BURN_2ND_ORDER
#else
#undef RESCALE_CONSERVED_QUANTITIES
#undef BURN_2ND_ORDER
#endif /* NOMOTO_MODEL */

c     ==================================================================

      MODULE coconut

c     ==================================================================

      USE PRECISION

      IMPLICIT NONE

c     real (kind=rk) :: qtot(1:m,5),xtot(1:m,config%qn)
      real(kind=rk), allocatable, save :: qtot(:,:), xtot(:,:)
c     real (kind=rk) ::  sumdvy(1:m), sumdvy2(1:m)

      real(kind=rk), allocatable :: sumdvy(:), sumdvy2(:)
      real(kind=rk), save :: l_z, om_z, p_z, v_z

      CONTAINS

c     ==================================================================

      SUBROUTINE init_coconut(restrt)

c     ==================================================================

      USE precision

      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
      USE conserved_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE fluxes_cfc
      USE sources_cfc
      USE metric_cfc
      USE boundary_cfc
      USE perm_aux_cfc
      USE temp_aux_cfc

      USE coconut_main
      USE enecons_cfc

      USE param_rt

      USE filtering 
      USE mesh_coarsening
      USE metric
      USE gr_initial_setup
#ifndef CFC_TRANSPORT2
      USE gr_sources
#endif /** CFC_TRANSPORT2 */
      USE delep_scheme
      USE gr_tstep
      USE cons_check
      USE renormalize

      USE tracer
      USE tracer_cfc
#if 0
      USE timmes_network
#endif
#ifdef EXTRACT_GW
      USE gw3d
#endif /* EXTRACT_GW */

      USE configure

      IMPLICIT NONE

      integer (kind=ik), intent(in) :: restrt
      
#ifdef RESCALE_CONSERVED_QUANTITIES
      real (kind=rk) :: d_cap (1:m,n_s:n_e,o_s:o_e)
      real (kind=rk) :: s_1   (1:m,n_s:n_e,o_s:o_e)
      real (kind=rk) :: s_2   (1:m,n_s:n_e,o_s:o_e)
      real (kind=rk) :: s_3   (1:m,n_s:n_e,o_s:o_e)
      real (kind=rk) :: tau   (1:m,n_s:n_e,o_s:o_e)
#endif
      
      integer (kind=ik) :: i, j, k, jk, l, i_n

      logical :: restart

      restart=restrt .ne. 0

      if (config%nsdim .eq. 3) then
#ifdef COARSENING         
         print *,'Setup: mesh coarsening...'
         call setup_coarse_grid
         print *,'Done!'
         call coarse_grid_average (0)
         call prolongation
         call initialize_conserved_quantities
#endif         
c         call polar_filter (2)
      end if
      if (config%p_ntr.eq.0) call read_in_y_e_profile

      if (.not.restart) then
         initial_step = .true.
         iteration=1
      end if

ccc   input control parameters
ccc   set Runge-Kutta parameters
      if (rk_order .eq. 2) then
ccc   2nd order Runge-Kutta
         alpha_rk = 2.0_rk
         beta_rk = 1.0_rk
      else
ccc   3nd order Runge-Kutta
         alpha_rk = 4.0_rk
         beta_rk = 3.0_rk
      endif

      metric_step = 1

ccc   initial metric calculation 

      print *,'First metric calculation...'
      call renormalize_mass_fractions
      call calculate_metric
      print *,'Done.'

#ifdef RESCALE_CONSERVED_QUANTITIES
      if (.not.restart) then
         d_cap (1:m,n_s:n_e,o_s:o_e) = d_cap_hat (1:m,n_s:n_e,o_s:o_e)
         s_1   (1:m,n_s:n_e,o_s:o_e) = s_1_hat   (1:m,n_s:n_e,o_s:o_e)
         s_2   (1:m,n_s:n_e,o_s:o_e) = s_2_hat   (1:m,n_s:n_e,o_s:o_e)
         s_3   (1:m,n_s:n_e,o_s:o_e) = s_3_hat   (1:m,n_s:n_e,o_s:o_e)
         tau   (1:m,n_s:n_e,o_s:o_e) = tau_hat   (1:m,n_s:n_e,o_s:o_e)

         alpha_alt(1:m,n_s:n_e,o_s:o_e)=1.0_rk
         do l = 1, 8
            d_cap_hat(1:m,n_s:n_e,o_s:o_e)=
     &           d_cap_hat(1:m,n_s:n_e,o_s:o_e)*
     &           (sqrt_gamma(1:m,n_s:n_e,o_s:o_e)/
     &           alpha_alt(1:m,n_s:n_e,o_s:o_e))
            s_1_hat(1:m,n_s:n_e,o_s:o_e)=
     &           s_1_hat(1:m,n_s:n_e,o_s:o_e)*
     &           (sqrt_gamma(1:m,n_s:n_e,o_s:o_e)/
     &           alpha_alt(1:m,n_s:n_e,o_s:o_e))
            s_2_hat(1:m,n_s:n_e,o_s:o_e)=
     &           s_2_hat(1:m,n_s:n_e,o_s:o_e)*
     &           (sqrt_gamma(1:m,n_s:n_e,o_s:o_e)/
     &           alpha_alt(1:m,n_s:n_e,o_s:o_e))
            s_3_hat(1:m,n_s:n_e,o_s:o_e)=
     &           s_3_hat(1:m,n_s:n_e,o_s:o_e)*
     &           (sqrt_gamma(1:m,n_s:n_e,o_s:o_e)/
     &           alpha_alt(1:m,n_s:n_e,o_s:o_e))
            tau_hat(1:m,n_s:n_e,o_s:o_e)=
     &           tau_hat(1:m,n_s:n_e,o_s:o_e)*
     &           (sqrt_gamma(1:m,n_s:n_e,o_s:o_e)/
     &           alpha_alt(1:m,n_s:n_e,o_s:o_e))
            do i_n=1,config%qn
               d_cap_xnu_hat(1:m,n_s:n_e,o_s:o_e,i_n)=
     &              d_cap_xnu_hat(1:m,n_s:n_e,o_s:o_e,i_n)*
     &           (sqrt_gamma(1:m,n_s:n_e,o_s:o_e)/
     &           alpha_alt(1:m,n_s:n_e,o_s:o_e))
            end do
            alpha_alt(1:m,n_s:n_e,o_s:o_e)=
     &           sqrt_gamma(1:m,n_s:n_e,o_s:o_e)
         end do
      end if
      call calculate_metric
#endif /* RESCALE_CONSERVED_QUANTITIES */

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP PARALLEL
#endif
#ifndef CFC_TRANSPORT2
      call calculate_newtonian_potential
      call add_newtonian_sources
#endif

      if (.not.restart) then

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif

         delta_t_old = 0.0_rk
         delta_t = c_light * 1.0e-5_rk

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D)         )
!$OMP END SINGLE
!$OMP BARRIER
#endif

ccc   calculate integral quantities
         print *,'Conservation check...'
         call do_conservation_check
         

ccc   compute initial time step
         
         print *,'Compute intitial time step...'
         call compute_time_step

#ifdef EXTRACT_GW
      if (config%nsdim .eq. 3) then
         call extract_gw (delta_t)
      end if
#endif /* EXTRACT_GW */
         
ccc   set old time steps as dummy on the first time slice
         delta_t_old_old = delta_t_old
         delta_t_old = delta_t

      endif

#ifdef CFC_TRANSPORT2
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         do i = 1, m
            alpha_alt (i, j, k) = alpha(i, j, k)
         end do
      end do

#else /* CFC_TRANSPORT2 */

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         do i = 1, m
            phi_potential_alt (i, j, k)=
     &           phi_potential (i, j, k)
         end do
      end do
#endif /* CFC_TRANSPORT2 */


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END PARALLEL
#endif

      initial_step = .false.


#if 0
      CALL init_network
#endif

      return

      END SUBROUTINE init_coconut

c     ==================================================================




c     ==================================================================

      SUBROUTINE coconut_hydro(dt_act,dt_new)

c     ==================================================================

      USE precision


      USE param_rt

      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
      USE conserved_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE fluxes_cfc
      USE sources_cfc
      USE metric_cfc
      USE boundary_cfc
      USE perm_aux_cfc
      USE temp_aux_cfc
      USE enecons_cfc
      USE coconut_main

      USE nutrio_hy, ONLY: qentot,qmotot,qmytot,qyetot 
      USE phycon

      USE metric
      USE grids
      USE cons_check
      USE gr_boundary
      USE gr_ppm
      USE gr_flux_r
      USE gr_flux_theta
      USE gr_flux_varphi
      USE gr_reconstruction
      USE gr_rescale
      USE gr_sources
      USE gr_tstep
#ifdef CFC_MHD
      USE induction
#endif /* CFC_MHD */
      USE recover_prim_vars
      USE renormalize
      USE delep_scheme
      USE gw_extraction
      USE mesh_coarsening
      USE filtering
      USE configure
      USE hydro_areas_mod

      USE mo_mpi

      USE hydro_areas_mod
#ifdef EXTRACT_GW
      USE gw3d
#endif /* EXTRACT_GW */

      IMPLICIT NONE

      real (kind=rk), intent(in) :: dt_act
      real (kind=rk), intent(out) :: dt_new
      real (kind=rk) :: dt1,dt2,dt_tmp
      real :: wm1_v2w,v2,winv
      
      logical :: bounce, bounce2

      real (kind=rk), parameter :: rho_trapping = 2.0e12_rk

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
!$    integer OMP_GET_NUM_THREADS
#endif

      integer i, j, k, jk, l, i_n, ierr
      
      delta_t = dt_act * c_light

#ifndef NOTRA
      call map_nutra2cfc
#endif

c      bounce = any(rho(1:m,n_s:n_e,o_s:o_e) .gt. rho_trapping*pc_ggeo
c     &     .and. entropy(1:m,n_s:n_e,o_s:o_e) .gt. 3.0_rk)
c      bounce = bounce .and. any(rho(1:m,n_s:n_e,o_s:o_e) .gt. 2.7e14_rk*
c     &     pc_ggeo)
      bounce = .false.
      if (n_s.le.n/2+1.and.n_e.ge.n/2+1) then
         bounce = any(rho(1:m,n/2+1,o_s:o_e) .gt. rho_trapping*pc_ggeo
     &        .and. entropy(1:m,n/2+1,o_s:o_e) .gt. 3.0_rk)
         bounce = bounce .and. any(rho(1:m,n/2+1,o_s:o_e).gt.2.7e14_rk*
     &        pc_ggeo)
      end if

#ifdef MPI_HYDRO
      call MPI_Allreduce (bounce, bounce2, 1_ik,
     &     MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      bounce = bounce2
#endif

#ifdef CONVECTION_1D
c     Adjust spherical core
      if (areas%ix_are(1, 6) .eq. config%qy) then
         core_loop: do i = 1, config%qx
            if (MINVAL(rho(i,n_s:n_e,o_s:o_e)) .lt. 
     &        5e11_rk*pc_ggeo) then
               areas%ix_are(1, 2) = MIN(i,72)
               exit core_loop
            end if
         end do core_loop
         if (MAXVAL(rho(1:m,n_s:n_e,o_s:o_e)) .lt. 1e14_rk*pc_ggeo) then
            areas%ix_are(1,2)=120
         endif
         i = areas%ix_are(1, 2)
         call MPI_Allreduce (i, areas%ix_are(1,2), 1_ik,
     &        MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
c         print *,'Core radius [km]:',r(areas%ix_are(1,2))/1e5_rk
c         areas%ix_are(1,2)=MAX(areas%ix_are(1,2),32)
c         areas%ix_are(1,2)=MAX(32,
         areas%ix_are(1,2)=MAX(16,areas%ix_are(1,2))
         areas%ix_are(2,1)=MIN(areas%ix_are(1,2)+1,config%qx)
      end if
#endif /* CONVECTION_1D */

#if defined(OPENMP_CFC) && (defined (OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP PARALLEL 
!$OMP& PRIVATE (i,j,jk,i_n,k,rk_index)
!$OMP& SHARED(delta_t,dt1,dt2,dt_act,dt_new)

!$OMP SINGLE
#endif
      t_total = t_total + delta_t

#if defined(OPENMP_CFC) && (defined (OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
!$OMP BARRIER
#endif

ccc   store hydro and metric on old time slice

      call smooth_entropy_holes
      call renormalize_mass_fractions
      
      if (config%p_ntr.eq.0) call fast_neutrino (bounce)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m

            d_cap_hat_old (i, j, k) = d_cap_hat (i, j, k)
            s_1_hat_old (i, j, k) = s_1_hat (i, j, k)
            s_2_hat_old (i, j, k) = s_2_hat (i, j, k)
            s_3_hat_old (i, j, k) = s_3_hat (i, j, k)
            tau_hat_old (i, j, k) = tau_hat (i, j, k)
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
            d_cap_xnu_hat_old (i, j, k, 1 : config%qn) =
     &           d_cap_xnu_hat (i, j, k, 1 : config%qn)
#ifdef CFC_MHD
            b_cap_1_hat_old (i, j, k) = b_cap_1_hat (i, j, k)
            b_cap_2_hat_old (i, j, k) = b_cap_2_hat (i, j, k)
            b_cap_3_hat_old (i, j, k) = b_cap_3_hat (i, j, k)
            psi_old (i, j, k) = psi (i, j, k)
#ifdef CFC_TRANSPORT2
#error
#endif            
            eps_old (i, j, k) = eps (i, j, k) * d_cap_hat (i, j, k)
#endif /* CFC_MHD */
         enddo
      enddo

ccc   start Runge-Kutta step -------------------------------------------
      do rk_index = 1, rk_order

c         print *,'Conservation check 1...'
c         call do_conservation_check
         
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
#ifdef COARSENING
          if (config%nsdim .eq. 3)   call prolongation
#endif
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif
         
         if (rk_index .eq. 1) then 
            dt_tmp=dt_act
         else
            dt_tmp=0.5_rk*dt_act
         endif

         call rescale_velocities (1)
         
#ifdef MPI_HYDRO
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
         call set_ghost_zones_mpi
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif
#endif /*MPI_HYDRO*/

ccc   enforce boundary conditions
         call hydro_boundary_conditions

         call detect_sh_cfc
#ifdef CFC_MHD
         call emf_centre
#endif
         
ccc   reconstruction in r direction
         call reconstruct_interface_r

ccc   numerical fluxes in r direction
         call calculate_flux_r
         call calculate_flux_differences (1)

         if (config%nsdim .ge. 2) then
ccc   reconstruction in theta direction
            call reconstruct_interface_theta


ccc   numerical fluxes in theta direction
            call calculate_flux_theta
            call calculate_flux_differences (2)

            if (config%nsdim .eq. 3) then
c#ifdef COARSENING
c               call set_ghost_zones_mpi
c               call prolongation
c               call emf_centre
c#endif /* COARSENING */               
ccc   reconstruction in varphi direction
               call reconstruct_interface_varphi

ccc   numerical fluxes in varphi direction
               call calculate_flux_varphi
c               call polar_filter (0)
               call calculate_flux_differences (3)
            endif !  config%nsdim .eq. 3
         endif ! config%nsdim .ge. 2

         call rescale_velocities (2)

#ifdef CONSTRAINED
#ifdef CFC_MHD
         if (config%nsdim .eq. 2)
     &        call emf_2d
c#endif /* CFC_MHD */
c#ifdef CFC_MHD
         call induction_boundary_conditions
         call set_ghost_zones_induction_mpi
         call rescale_velocities (1)
         call calculate_db_dt (0)
         call rescale_velocities (2)
#endif /* CFC_MHD */
#endif /* CONSTRAINED */         

         call calculate_flux_differences (4)
c         print *,'Conservation check 1a...'
c         call do_conservation_check
         call avare_hy (1)
c         print *,'Conservation check 1b...'
c         call do_conservation_check
ccc   compute source terms
         call calculate_sources


#ifdef AXIS
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            if (j .eq. 1 .or. j .eq. n) then
               do i = 1, m
                  f_diff  (i, j, k, 3) = 0.0_rk
                  s_total (i, j, k, 3) = 0.0_rk
                  f_diff  (i, j, k, 4) = 0.0_rk
                  s_total (i, j, k, 4) = 0.0_rk
               end do
            end if
         end do
#endif /* AXIS */

         if (rk_index .eq. 1) then
ccc   first Runge-Kutta step

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
            do jk = 1, n_loc * o_loc

               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               k = k + o_s - 1

               do i = 1, m
                     
                  d_cap_hat (i, j, k) =
     &                 d_cap_hat_old (i, j, k) -
     &                 delta_t * f_diff (i, j, k, 1)
                  
                  s_1_hat (i, j, k) = 
     &                 s_1_hat_old (i, j, k) +
     &                 delta_t * (- f_diff (i, j, k, 2) +
     &                 alpha (i, j, k) * s_total (i, j, k, 2))

                  s_2_hat (i, j, k) = 
     &                 s_2_hat_old (i, j, k) +
     &                 delta_t * (- f_diff (i, j, k, 3) +
     &                 alpha (i, j, k) * s_total (i, j, k, 3))

                  if (config%nsdim .lt. 3) then
                     s_3_hat (i, j, k) =
     &                   s_3_hat_old (i, j, k) -
     &                   delta_t * f_diff (i, j, k, 4) 
                  else
                      s_3_hat (i, j, k) =
     &                    s_3_hat_old (i, j, k) -
     &                    delta_t * (f_diff (i, j, k, 4) -
     &                    alpha (i, j, k) * s_total (i, j, k, 4))
                  endif

                  tau_hat (i, j, k) =
     &                 tau_hat_old (i, j, k) +
     &                 delta_t * (- f_diff (i, j, k, 5)+
     &                 alpha (i, j, k) * s_total (i, j, k, 5))
#ifdef ENECONS
     &                 / alpha (i, j, k)
#endif
                  
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
                  do i_n = 1, config%qn 
                     d_cap_xnu_hat (i, j, k, i_n) =
     &                    d_cap_xnu_hat_old (i, j, k, i_n) - 
     &                    delta_t * f_diff_xnu (i, j, k, i_n)
                  enddo

#ifdef CFC_MHD
#ifndef CONSTRAINED
                  b_cap_1_hat (i, j, k) = 
     &                 b_cap_1_hat_old (i, j, k) +
     &                 delta_t * (- f_diff (i, j, k, 6) +
     &                 alpha (i, j, k) * s_total (i, j, k, 6))
                  
                  b_cap_2_hat (i, j, k) = 
     &                 b_cap_2_hat_old (i, j, k) +
     &                 delta_t * (- f_diff (i, j, k, 7) +
     &                 alpha (i, j, k) * s_total (i, j, k, 7)) *
     &                 r (i)

                  b_cap_3_hat (i, j, k) = 
     &                 b_cap_3_hat_old (i, j, k) +
     &                 delta_t * (- f_diff (i, j, k, 8) +
     &                 alpha (i, j, k) * s_total (i, j, k, 8)) *
     &                 r (i) * sinus_theta (j)

                  psi (i, j, k) = 
     &                 psi_old (i, j, k) +
     &                 delta_t * (- f_diff (i, j, k, 9) +
     &                 alpha (i, j, k) * s_total (i, j, k, 9))

#else /* CONSTRAINED */
                  b_cap_1_hat (i, j, k) = 
     &                 b_cap_1_hat_old (i, j, k) +
     &                 delta_t * f_diff_b (i, j, k, 1)

                  b_cap_2_hat (i, j, k) = 
     &                 b_cap_2_hat_old (i, j, k) +
     &                 delta_t * f_diff_b (i, j, k, 2)

                  b_cap_3_hat (i, j, k) =
     &                 b_cap_3_hat_old (i, j, k) +
     &                 delta_t * f_diff_b (i, j, k, 3) 
#endif /* CONSTRAINED */
#ifdef CFC_TRANSPORT2
#error
#endif
                  eps_aux (i, j, k) =
     &                 eps_old (i, j, k) -
     &                 delta_t * (f_diff (i, j, k, 10) -
     &                 p (i, j, k) * f_diff (i, j, k, 1) /
     &                 d_cap_hat (i , j, k))
#endif /* CFC_MHD */

               enddo            !i-loop
            enddo               !j,k-loop
            
         elseif (rk_index .eq. 2) then
ccc   second Runge-Kutta step

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
            do jk = 1, n_loc * o_loc

               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               k = k + o_s - 1

               do i = 1, m
              
                  d_cap_hat (i, j, k) =
     &                 (beta_rk * d_cap_hat_old (i, j, k) +
     &                 d_cap_hat (i, j, k) - delta_t * 
     &                 f_diff (i, j, k, 1)) / alpha_rk
                  
                  s_1_hat (i, j, k) =
     &                 (beta_rk * s_1_hat_old (i, j, k) +
     &                 s_1_hat (i, j, k) + delta_t *
     &                 (- f_diff (i, j, k, 2) +
     &                 alpha (i, j, k) * s_total (i, j, k, 2))) / 
     &                 alpha_rk

                  s_2_hat (i, j, k) = 
     &                 (beta_rk * s_2_hat_old (i, j, k) +
     &                 s_2_hat (i, j, k) + delta_t *
     &                 (- f_diff (i, j, k, 3) +
     &                 alpha (i, j, k) * s_total (i, j, k, 3))) / 
     &                 alpha_rk

                  if (config%nsdim .lt. 3) then
                     s_3_hat (i, j, k) =
     &                   (beta_rk * s_3_hat_old (i, j, k) +
     &                   s_3_hat (i, j, k) - delta_t *
     &                   f_diff (i, j, k, 4)) /
     &                   alpha_rk
                  else
                     s_3_hat (i, j, k) =
     &                   (beta_rk * s_3_hat_old (i, j, k) +
     &                   s_3_hat (i, j, k) - delta_t *
     &                   (f_diff (i, j, k, 4) -
     &                   alpha (i, j, k) * s_total (i, j, k, 4))) / 
     &                   alpha_rk
                  endif


                  tau_hat (i, j, k) =
     &                 (beta_rk * tau_hat_old (i, j, k) +
     &                 tau_hat (i, j, k) + delta_t * 
     &                 (- f_diff (i, j, k, 5) +
     &                 alpha (i, j, k) *
     &                 s_total (i, j, k, 5))
#ifdef ENECONS
     &                 / alpha (i, j, k)
#endif
     &                 ) / alpha_rk

#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
                  do i_n = 1, config%qn
                     d_cap_xnu_hat (i, j, k, i_n) =
     &                    (beta_rk *
     &                    d_cap_xnu_hat_old (i, j, k, i_n) +
     &                    d_cap_xnu_hat (i, j, k, i_n) - delta_t * 
     &                    f_diff_xnu (i, j, k, i_n)) / alpha_rk                     
                  enddo        

#ifdef CFC_MHD
#ifndef CONSTRAINED
                  b_cap_1_hat (i, j, k) =
     &                 (beta_rk * b_cap_1_hat_old (i, j, k) +
     &                 b_cap_1_hat (i, j, k) + delta_t *
     &                 (- f_diff (i, j, k, 6) +
     &                 alpha (i, j, k) * s_total (i, j, k, 6))) / 
     &                 alpha_rk

                  b_cap_2_hat (i, j, k) =
     &                 (beta_rk * b_cap_2_hat_old (i, j, k) +
     &                 b_cap_2_hat (i, j, k) + delta_t *
     &                 (- f_diff (i, j, k, 7) +
     &                 alpha (i, j, k) * s_total (i, j, k, 7)) *
     &                 r (i)) / 
     &                 alpha_rk

                  b_cap_3_hat (i, j, k) =
     &                 (beta_rk * b_cap_3_hat_old (i, j, k) +
     &                 b_cap_3_hat (i, j, k) + delta_t *
     &                 (- f_diff (i, j, k, 8) +
     &                 alpha (i, j, k) * s_total (i, j, k, 8)) *
     &                 r (i) * sinus_theta (j)) / 
     &                 alpha_rk                  

                  psi (i, j, k) =
     &                 (beta_rk * psi_old (i, j, k) +
     &                 psi (i, j, k) + delta_t *
     &                 (- f_diff (i, j, k, 9) +
     &                 alpha (i, j, k) * s_total (i, j, k, 9))) / 
     &                 alpha_rk                  

#else /* CONSTRAINED */
                  b_cap_1_hat (i, j, k) =
     &                 (beta_rk * b_cap_1_hat_old (i, j, k) +
     &                 b_cap_1_hat (i, j, k) + delta_t *
     &                 f_diff_b (i, j, k, 1)) /
     &                 alpha_rk

                  b_cap_2_hat (i, j, k) =
     &                 (beta_rk * b_cap_2_hat_old (i, j, k) +
     &                 b_cap_2_hat (i, j, k) + delta_t *
     &                 f_diff_b (i, j, k, 2)) /
     &                 alpha_rk

                  b_cap_3_hat (i, j, k) =
     &                 (beta_rk * b_cap_3_hat_old (i, j, k) +
     &                 b_cap_3_hat (i, j, k) + delta_t *
     &                 f_diff_b (i, j, k, 3)) /
     &                 alpha_rk
#endif /* CONSTRAINED */
                  
                  eps_aux (i, j, k) =
     &                 (beta_rk * eps_old (i, j, k) +
     &                 eps_aux (i, j, k) - delta_t *
     &                 (f_diff (i, j, k, 10) -
     &                 p (i, j, k) * f_diff (i, j, k, 1) /
     &                 d_cap_hat (i , j, k))) /
     &                 alpha_rk
#endif /* CFC_MHD */                  
                  
               enddo            !i-loop
            enddo               !j,k-loop

         elseif (rk_index .eq. 3) then
ccc   third Runge-Kutta step

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
            do jk = 1, n_loc * o_loc

               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               k = k + o_s - 1

               do i = 1, m
                     
                  d_cap_hat (i, j, k) =
     &                 (d_cap_hat_old (i, j, k) +
     &                 2.0_rk *
     &                 (d_cap_hat (i, j, k) - delta_t *
     &                 f_diff (i, j, k, 1))) / 3.0_rk

                  s_1_hat (i, j, k) =
     &                    (s_1_hat_old (i, j, k) + 2.0_rk *
     &                    (s_1_hat (i, j, k) + delta_t *
     &                    (- f_diff (i, j, k, 2) +
     &                    alpha (i, j, k) * s_total (i, j, k, 2)))) /
     &                    3.0_rk

                  s_2_hat (i, j, k) =
     &                 (s_2_hat_old (i, j, k) + 2.0_rk *
     &                 (s_2_hat (i, j, k) + delta_t *
     &                 (- f_diff (i, j, k, 3) +
     &                 alpha (i, j, k) * s_total (i, j, k, 3)))) /
     &                 3.0_rk

                  if (config%nsdim .lt. 3) then
                     s_3_hat (i, j, k) =
     &                  (s_3_hat_old (i, j, k) + 2.0_rk *
     &                  (s_3_hat (i, j, k) - delta_t * 
     &                  f_diff (i, j, k, 4))) /
     &                  3.0_rk
                  else
                     s_3_hat (i, j, k) =
     &                   (s_3_hat_old (i, j, k) + 2.0_rk *
     &                   (s_3_hat (i, j, k) - delta_t * 
     &                   (f_diff (i, j, k, 4) -
     &                   alpha (i, j, k) * s_total (i, j, k, 4)))) /
     &                  3.0_rk
                  endif


                  tau_hat (i, j, k) =
     &                 (tau_hat_old (i, j, k) + 2.0_rk *
     &                 (tau_hat (i, j, k) + delta_t * 
     &                 (- f_diff (i, j, k, 5) +
     &                 alpha (i, j, k) *
     &                 s_total (i, j, k, 5))
#ifdef ENECONS
     &                 / alpha (i, j, k)
#endif
     &                 )) / 3.0_rk

#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
                  do i_n = 1, config%qn
                     d_cap_xnu_hat (i, j, k, i_n) =
     &                    (d_cap_xnu_hat_old (i, j, k, i_n) +
     &                    2.0_rk *
     &                    (d_cap_xnu_hat (i, j, k, i_n) -
     &                    delta_t * f_diff_xnu (i, j, k, i_n))) / 
     &                    3.0_rk
                  enddo

#ifdef CFC_MHD
#ifndef CONSTRAINED
                  b_cap_1_hat (i, j, k) =
     &                 (b_cap_1_hat_old (i, j, k) + 2.0_rk *
     &                 (b_cap_1_hat (i, j, k) + delta_t *
     &                 (- f_diff (i, j, k, 6) +
     &                 alpha (i, j, k) * s_total (i, j, k, 6)))) /
     &                 3.0_rk

                  b_cap_2_hat (i, j, k) =
     &                 (b_cap_2_hat_old (i, j, k) + 2.0_rk *
     &                 (b_cap_2_hat (i, j, k) + delta_t *
     &                 (- f_diff (i, j, k, 7) +
     &                 alpha (i, j, k) * s_total (i, j, k, 7)) *
     &                 r (i))) /
     &                    3.0_rk

                  b_cap_3_hat (i, j, k) =
     &                 (b_cap_3_hat_old (i, j, k) + 2.0_rk *
     &                 (b_cap_3_hat (i, j, k) + delta_t *
     &                 (- f_diff (i, j, k, 8) +
     &                 alpha (i, j, k) * s_total (i, j, k, 8)) *
     &                 r (i) * sinus_theta (j))) /
     &                 3.0_rk                  

                  psi (i, j, k) =
     &                 (psi_old (i, j, k) + 2.0_rk *
     &                 (psi (i, j, k) + delta_t *
     &                 (- f_diff (i, j, k, 9) +
     &                 alpha (i, j, k) * s_total (i, j, k, 9)))) /
     &                 3.0_rk
#else /* CONSTRAINED */                  
                  b_cap_1_hat (i, j, k) =
     &                 (b_cap_1_hat_old (i, j, k) + 2.0_rk *
     &                 (b_cap_1_hat (i, j, k) + delta_t * 
     &                 f_diff_b (i, j, k, 1))) /
     &                 3.0_rk

                  b_cap_2_hat (i, j, k) =
     &                 (b_cap_2_hat_old (i, j, k) + 2.0_rk *
     &                 (b_cap_2_hat (i, j, k) + delta_t * 
     &                 f_diff_b (i, j, k, 2))) /
     &                 3.0_rk

                  b_cap_3_hat (i, j, k) =
     &                 (b_cap_3_hat_old (i, j, k) + 2.0_rk *
     &                 (b_cap_3_hat (i, j, k) + delta_t * 
     &                 f_diff_b (i, j, k, 3))) /
     &                 3.0_rk
#endif /* CONSTRAINED */
                  
                  eps_aux (i, j, k) =
     &                 (eps_old (i, j, k) + 2.0_rk *
     &                 (eps_aux (i, j, k) - delta_t * 
     &                 (f_diff (i, j, k, 10) -
     &                 p (i, j, k) * f_diff (i, j, k, 1) /
     &                 d_cap_hat (i , j, k)))) /
     &                 3.0_rk
#endif /* CFC_MHD */

               enddo            !i-loop
            enddo               !j,k-loop

         endif

c     ------------------------------------------------------------------
c     Apply neutrino source terms

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            do i = 1, m

#ifdef CFC_TRANSPORT2               
c     Relativistic version (Newtonian version still missing!!!):
               if (config%nsdim .ge. 2) then
               
                  v2=v_squared (i, j, k)
                  wm1_v2w=0.5_rk*(1.0_rk+v2*(0.25_rk+v2*(0.125_rk+
     &                 v2*(0.078125_rk+v2*(0.0546875_rk+
     &                 v2*(0.041015625_rk+
     &                 v2*(0.0322265625_rk+v2*0.02618408203125_rk)))))))
                  if (v2 .gt. 0.01_rk) wm1_v2w=
     &                 (w(i,j,k)-1.0_rk)/(v2*w(i,j,k))
                  winv=1.0_rk/w(i,j,k)

               s_2_hat(i,j,k)= s_2_hat(i,j,k)
     &              + dt_tmp*pc_egeo*
     &              (pc_cl*qmytot(i,j,k)*phi(i,j,k)**2*r(i)
     &              +pc_cl*qmotot(i,j,k)*
     &              v_1(i,j,k)/phi(i,j,k)**2*
     &              v_2(i,j,k)*wm1_v2w
     &              +v_2(i,j,k)*qentot(i,j,k))

               s_3_hat(i,j,k)= s_3_hat(i,j,k)+
     &              dt_tmp*pc_egeo*(
     &              +pc_cl*qmotot(i,j,k)*
     &              v_1(i,j,k)/phi(i,j,k)**2*
     &              v_3(i,j,k)*wm1_v2w
     &              +v_3(i,j,k)*qentot(i,j,k))
            endif ! config%nsdim .ge. 2

               if (config%nsdim .ge. 2) then
                  s_1_hat(i,j,k)= s_1_hat(i,j,k)
     &                + dt_tmp*pc_egeo*
     &                (pc_cl*qmotot(i,j,k)*
     &                (winv+v_1(i,j,k)**2*wm1_v2w)
     &                +v_1(i,j,k)*qentot(i,j,k))

                  tau_hat(i,j,k)= tau_hat(i,j,k)
     &                + dt_tmp*pc_egeo*
     &                (             qentot(i,j,k)+
     &                v_2(i,j,k)/(phi(i,j,k)**2*r(i))*
     &                qmytot(i,j,k)*pc_cl+
     &                v_1(i,j,k)/(phi(i,j,k)**2)*
     &                qmotot(i,j,k)*pc_cl)
               else
                  s_1_hat(i,j,k)= s_1_hat(i,j,k)
     &                + dt_tmp*pc_egeo*
     &               (pc_cl*qmotot(i,j,k)*phi(i,j,k)**2
     &               +v_1(i,j,k)*qentot(i,j,k))

                  tau_hat(i,j,k)= tau_hat(i,j,k)
     &                + dt_tmp*pc_egeo*
     &                (             qentot(i,j,k)+
     &                v_1(i,j,k)/(phi(i,j,k)**2)*
     &                qmotot(i,j,k)*pc_cl)
               endif ! config%nsdim .ge. 2
#else /* CFC_TRANSPORT2 */
c     Newtonian version
               s_1_hat(i,j,k)= s_1_hat(i,j,k)
     &              + dt_tmp*pc_egeo*
     &              pc_cl*qmotot(i,j,k)
               if (config%nsdim .ge. 2) then
                  s_2_hat(i,j,k)= s_2_hat(i,j,k)
     &                 + dt_tmp*pc_egeo*
     &                 pc_cl*qmytot(i,j,k)*r(i)
c                  s_3_hat(i,j,k)= s_3_hat(i,j,k)
                  tau_hat(i,j,k)= tau_hat(i,j,k)
     &                + dt_tmp*pc_egeo*
     &                 (qentot(i,j,k)+
     &                 v_2(i,j,k)*r_inv(i)*
     &                 qmytot(i,j,k)*pc_cl+
     &                 v_1(i,j,k)*
     &                 qmotot(i,j,k)*pc_cl)
               else
                  tau_hat(i,j,k)= tau_hat(i,j,k)
     &                + dt_tmp*pc_egeo*
     &                (qentot(i,j,k)+
     &                v_1(i,j,k)*qmotot(i,j,k)*pc_cl)
               endif ! config%nsdim .ge. 2               
#endif /* CFC_TRANSPORT2 */               

               
               d_cap_xnu_hat(i,j,k,config%qn)=
     &              d_cap_xnu_hat(i,j,k,config%qn)
     &              + dt_tmp*pc_ggeo*qyetot(i,j,k,1)

            enddo
#ifdef ORIGIN
            s_1_hat(1,j,k)=0.0_rk
#endif
#ifdef AXIS
            if (j.eq. 1 .or. j .eq. n) then
               s_2_hat(:,j,k)=0.0_rk
               s_3_hat(:,j,k)=0.0_rk
            end if
#endif /* AXIS */            
         enddo

c#if 0 /* geht noch nicht */
         if (config%nsdim .eq. 3) then
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
#ifdef COARSENING
c            print *,'Conservation check 2...'
c            call do_conservation_check
            call coarse_grid_average (0)
c            print *,'Conservation check 3...'
c            call do_conservation_check
#else            
            call polar_filter (1)
#endif
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif
         endif ! config%nsdim .eq. 3
c#endif 
         
         call avare_hy(0)

#ifdef BURN_2ND_ORDER
         call burn3d(0, rk_index, dt_act)
#endif

         call renormalize_mass_fractions
         call bndry_flx(rk_index,rk_order)

#if !defined(CFC_TRANSPORT2) && defined (ENECONS)
c     ------------------------------------------------------------------
c     update energy density after the update of the grav. potential

         if (rk_index .eq. rk_order) then
            call calculate_newtonian_potential
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
            do jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               do i=1,m
                  tau_hat(i,j,k)=tau_hat(i,j,k)+
     &                 0.5_rk*
     &                 (d_cap_hat(i,j,k)-
     &                 d_cap_hat_old(i,j,k))*
     &                 (phi_potential_alt(i,j,k)-
     &                 phi_potential(i,j,k))
                  phi_potential_alt(i,j,k)=
     &                 phi_potential(i,j,k)
               end do
            end do

         end if
#endif /* ! CFC_TRANSPORT2 && ENECONS */
         
ccc   recover the primitive (physical) variables rho, v^i, W, epsilon
         call recover_primitives

      enddo
ccc   end of Runge-Kutta loop


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END PARALLEL
#endif
      if (config%p_ntr .eq. 0 .and. .not. bounce) call
     &     deleptonization (bounce)

#ifndef BURN_2ND_ORDER
      call burn3d(1, 1, dt_act)
#endif

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP PARALLEL
!$OMP& PRIVATE (i,j,jk,i_n,k,rk_index)
!$OMP& SHARED (delta_t,dt1,dt2,dt_act,dt_new)
ccc   calculate new metric
#endif
      metric_calculated = .false.
#ifdef CFC_TRANSPORT2
      call calculate_metric

#ifdef ENECONS
c     ------------------------------------------------------------------
c     update energy density after calculating the new lapse function

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc
         
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m

c            tau_hat (i, j, k) =
c     &           (2 * alpha_alt (i, j, k) * tau_hat (i, j, k) +
c     &           (alpha (i, j, k) - alpha_alt (i, j, k)) *
c     &           (d_cap_hat_old (i, j, k) - d_cap_hat (i, j,k) +
c     &           tau_hat_old (i, j, k))) /
c     &           (alpha (i, j, k) + alpha_alt (i, j, k))
            tau_hat (i, j, k) = tau_hat (i, j, k) +
     &           0.5_rk *
     &           ((d_cap_hat (i, j, k) +
     &           tau_hat (i, j, k)) -
     &           (d_cap_hat_old (i, j, k) +
     &           tau_hat_old (i, j, k)))*
     &           (alpha_alt (i, j, k) -
     &           alpha     (i, j, k))/
     &           alpha_alt (i, j, k)
            alpha_alt(i, j, k)=alpha(i, j, k)

         end do
      end do
#endif /* ENECONS */
         
#endif /*CFC_TRANSPORT2*/

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
#ifdef COARSENING
c          if (config%nsdim .eq. 3)   call prolongation
#endif
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif

      
ccc   extract gravitational waves

c      if (config%nsdim .ge. 2) then
c         call extract_gravitational_waves
c      endif
      
#ifdef EXTRACT_GW
c      if (config%nsdim .eq. 3) then
         call extract_gw (delta_t)
c      end if
#endif /* EXTRACT_GW */
      
ccc   calculate time step
      call compute_time_step
      
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END PARALLEL
#endif
      
      dt_new=delta_t/c_light
      iteration=iteration+1

      return
      
      END SUBROUTINE coconut_hydro

c     ==================================================================




c     ==================================================================

      SUBROUTINE calculate_flux_differences (direction)

c     ==================================================================

      USE precision

      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
#ifdef CFC_MHD      
      USE hydro_primitives_cfc, ONLY: p, v_1, v_2, v_3, rho, v_squared,
     &     c_sound_squared, b_1, b_2, b_3
      USE conserved_cfc, ONLY:  b_cap_1_hat, b_cap_2_hat, b_cap_3_hat
#endif      
      USE fluxes_cfc
      USE metric_cfc
      USE perm_aux_cfc
      USE sources_cfc, ONLY: s_total
      USE shock_cfc
      
      USE configure

      IMPLICIT NONE

      integer (kind=ik), intent (in) :: direction

      real (kind=rk) :: fxr, fyr, fzr, fxl, fyl, fzl
      real (kind=rk) :: scr1, scr2, scr3
      real (kind=rk) :: delta_t_inv, vcln

      integer (kind=ik) :: i, j, k, jk, l


      select case (direction)

      case(1) !r-direction

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk,l)
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         delta_t_inv = 1.0_rk / delta_t

         do i = 1, m

               f_diff (i, j, k, 1) =
     &           (alpha_r (i, j, k) *
     &           sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 1) -
     &           alpha_r (i - 1, j, k) *
     &           sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 1)) *
     &           delta_r_mod_inv (i)
               
               f_diff (i, j, k, 2) =
     &           (alpha_r (i, j, k) *
     &           sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 2) -
     &           alpha_r (i - 1, j, k) *
     &           sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 2)) *
     &           delta_r_mod_inv (i) 

#if MOMENTUM_EQUATION_MODE>0
               f_diff (i, j, k, 2) = f_diff (i, j, k, 2) +
     &           (pav_r (i    , j, k) -
     &           pav_r (i - 1, j, k)) *
     &           alpha (i, j, k) *
     &           phi (i, j, k) ** 6 *
     &           delta_r_inv (i)
#endif /* MOMENTUM_EQUATION_MODE>0 */
        
#if MOMENTUM_EQUATION_MODE==2
               f_diff (i, j, k, 3) =
     &           (alpha_r (i, j, k) *
     &           sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 3) * r_if_inv (i) -
     &           alpha_r (i - 1, j, k) *
     &           sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 3) * r_if_inv (i - 1)) *
     &           delta_r_mod_inv (i) *
     &           r (i)
#else /* MOMENTUM_EQUATION_MODE==2 */
               f_diff (i, j, k, 3) =
     &           (alpha_r (i, j, k) *
     &           sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 3) -
     &           alpha_r (i - 1, j, k) *
     &           sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 3)) *
     &           delta_r_mod_inv (i)
#endif /* MOMENTUM_EQUATION_MODE==2 */

#if 0
c     Advect orthonormal components of S_3 (currently not recommended)
               f_diff (i, j, k, 4) =
     &              (alpha_r (i, j, k) *
     &              sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 4) * r_if_inv (i) -
     &              alpha_r (i - 1, j, k) *
     &              sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 4)  * r_if_inv (i - 1)) *
     &              delta_r_mod_inv (i) * r (i)
#else
               f_diff (i, j, k, 4) =
     &              (alpha_r (i, j, k) *
     &              sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 4) -
     &              alpha_r (i - 1, j, k) *
     &              sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 4)) *
     &           delta_r_mod_inv (i)
#endif

#ifdef ENECONS
#ifdef CFC_TRANSPORT2
c     improved formulation of GR energy equation
               f_diff (i, j, k, 5) =
     &           ((alpha_r (i, j, k) ** 2 *
     &           sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 5) -
     &           alpha_r (i - 1, j, k) ** 2 *
     &           sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 5)) +
     &           (alpha_r (i, j, k) * 
     &           sqrt_gamma_r (i, j, k) *
     &           (alpha_r (i, j, k) - alpha (i, j, k)) *
     &              f_num_if (i, j, k, 1) -
     &           alpha_r (i - 1, j, k) *
     &           sqrt_gamma_r (i - 1, j, k) *
     &           (alpha_r (i - 1, j, k) - alpha (i, j, k)) *
     &              f_num_if (i - 1, j, k, 1))) *
     &           delta_r_mod_inv (i)
#else
c     improved formulation of Newtonian energy equation
            f_diff (i, j, k, 5) =
     &           (sqrt_gamma_r (i, j, k) *
     &           (f_num_if (i, j, k, 5) + 0.5_rk *
     &           (phi_potential (i + 1, j, k) -
     &           phi_potential (i, j, k)) *
     &           f_num_if (i, j, k, 1)) -
     &           sqrt_gamma_r (i - 1, j, k) *
     &           (f_num_if (i - 1, j, k, 5) + 0.5_rk *
     &           (phi_potential (i - 1, j, k) -
     &           phi_potential (i, j, k)) *
     &           f_num_if (i - 1, j, k, 1))) *
     &              delta_r_mod_inv (i)
#endif /* CFC_TRANSPORT2 */

#else

c Standard formulation of energy equation:

               f_diff (i, j, k, 5) =
     &           (alpha_r (i, j, k) *
     &           sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 5) -
     &           alpha_r (i - 1, j, k) *
     &           sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 5)) *
     &           delta_r_mod_inv (i)

#endif /* ENECONS */
               
#ifdef CFC_MHD
               f_diff (i, j, k, 6) = 0.0_rk

c     It is more convenient to put the grad psi term into
c     s_total, but compute it here where f_num_if is available
               s_total (i, j, k, 6) =
     &              - sqrt_gamma (i, j, k) * g_up_11 (i, j, k) *
     &              (f_num_if (i, j, k, 6) -
     &              f_num_if (i - 1, j, k, 6)) *
     &              delta_r_inv (i) !* 0.0_rk
               
               f_diff (i, j, k, 7) =
     &              (alpha_r (i, j, k) *
     &              sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 7) -
     &              alpha_r (i - 1, j, k) *
     &              sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 7)) *
     &              delta_r_mod_inv (i) !*
c     &              g_22 (i, j, k)

               f_diff (i, j, k, 8) =
     &              (alpha_r (i, j, k) *
     &              sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 8) -
     &              alpha_r (i - 1, j, k) *
     &              sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 8)) *
     &              delta_r_mod_inv (i) 
c     &              g_33 (i, j, k)
               
               f_diff (i, j, k, 9) =
     &              (alpha_r (i, j, k) *
     &              sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 9) -
     &              alpha_r (i - 1, j, k) *
     &              sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 9)) * 
     &              delta_r_mod_inv (i)
               
               f_diff (i, j, k, 11) =
     &              (alpha_r (i, j, k) *
     &              sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 11) -
     &              alpha_r (i - 1, j, k) *
     &              sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 11)) * 
     &              delta_r_mod_inv (i)
c               if (j.eq.n_s.and.k.eq.o_e) print *,'fdiff1',
c     &              i,j,k,f_diff (i, j, k, 9), f_diff (i, j, k, 11)
               
#ifdef CFC_TRANSPORT2
#error
#endif               
c     internal energy equation
               f_diff (i, j, k, 10) =
     &              (alpha_r (i, j, k) *
     &              sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 10) -
     &              alpha_r (i - 1, j, k) *
     &              sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 10)) *
     &              delta_r_mod_inv (i) + 0.25_rk *
     &              (alpha_r (i, j, k) *
     &              sqrt_gamma_r (i, j, k) *
     &              f_num_if (i, j, k, 1) *
     &              (2*p (i, j, k) + 0*p (i + 1, j, k)) *
     &              (1.0_rk / rho (i, j, k) -
     &              1.0_rk / rho (i + 1, j, k)) -
     &              alpha_r (i - 1, j, k) *
     &              sqrt_gamma_r (i - 1, j, k) *
     &              f_num_if (i - 1, j, k, 1) *
     &              (2*p (i, j, k) + 0*p (i - 1, j, k)) *
     &              (1.0_rk / rho (i, j, k) -
     &              1.0_rk / rho (i - 1, j, k))) *
     &              delta_r_mod_inv (i) * rho (i, j, k)              
#endif /* CFC_MHD */


#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
            do l = 1, config%qn
                  f_diff_xnu (i, j, k, l) =
     &              (alpha_r (i, j, k) *
     &              sqrt_gamma_r (i, j, k) *
     &                 f_num_xnu_if (i, j, k, l) -
     &              alpha_r (i - 1, j, k) *
     &              sqrt_gamma_r (i - 1, j, k) *
     &                 f_num_xnu_if (i - 1, j, k, l)) *
     &              delta_r_mod_inv (i)
            enddo
            
         enddo                  !i-loop

#ifdef ORIGIN
            f_diff (1, j, k, 2) = 0.0_rk
#endif
      enddo                     !j,k-loop


      case(2) !theta-direction
      
         call avare_hy(2)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk,l)
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            do i = 1, m

               f_diff (i, j, k, 1) = f_diff (i, j, k, 1) +
     &             (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 1) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 1)) *
     &              delta_theta_mod_inv (j)
               
#if MOMENTUM_EQUATION_MODE==2
               f_diff (i, j, k, 2) = f_diff (i, j, k, 2) +
     &              (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              (f_num_if (i, j, k, 2) * cos_dtheta -
     &              f_num_if (i, j, k, 3) * sin_dtheta *
     &              r_inv (i)) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              (f_num_if (i, j - 1, k, 2) * cos_dtheta +
     &              f_num_if (i, j - 1, k, 3) * sin_dtheta *
     &              r_inv (i))) *
     &              delta_theta_mod_inv (j)
               
               f_diff (i, j, k, 3) = f_diff (i, j, k, 3) +
     &              (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              (f_num_if (i, j, k, 3) * cos_dtheta +
     &              f_num_if (i, j, k, 2) * sin_dtheta * r (i)) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              (f_num_if (i, j - 1, k, 3) * cos_dtheta -
     &              f_num_if (i, j - 1, k, 2) * sin_dtheta * r (i)))
     &              * delta_theta_mod_inv (j)
#else /* MOMENTUM_EQUATION_MODE==2 */              
               f_diff (i, j, k, 2) = f_diff (i, j, k, 2) +
     &              (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 2) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 2)) *
     &              delta_theta_mod_inv (j)
               
               f_diff (i, j, k, 3) = f_diff (i, j, k, 3) +
     &              (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 3) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 3)) *
     &              delta_theta_mod_inv (j) 
#endif /* MOMENTUM_EQUATION_MODE==2 */


#if MOMENTUM_EQUATION_MODE>0
               f_diff (i, j, k, 3) = f_diff (i, j, k, 3) +
     &              (pav_theta (i, j    , k) - 
     &              pav_theta (i, j - 1, k)) *
     &              alpha (i, j, k) *
     &              phi (i, j, k) ** 6 *
     &              delta_theta_inv
#endif /* MOMENTUM_EQUATION_MODE>0 */

#if 0
               f_diff (i, j, k, 4) = f_diff (i, j, k, 4) +
     &              (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 4) *
     &              sinus_theta_inv_if (j) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 4) *
     &              sinus_theta_inv_if (j - 1)) *
     &              delta_theta_mod_inv (j) * sinus_theta (j)
#else
               f_diff (i, j, k, 4) = f_diff (i, j, k, 4) +
     &              (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 4) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 4)) *
     &              delta_theta_mod_inv (j)
#endif

#ifdef ENECONS
#ifdef CFC_TRANSPORT2
c     improved formulation of GR energy equation
               f_diff (i, j, k, 5) = f_diff (i, j, k, 5) +
     &              ((alpha_theta (i, j, k) ** 2 *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 5) -
     &              alpha_theta (i, j - 1, k) ** 2 *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 5)) +
     &              (alpha_theta (i, j, k) * 
     &              sqrt_gamma_theta (i, j, k) *
     &              (alpha_theta (i, j, k) - alpha (i, j, k)) *
     &              f_num_if (i, j, k, 1) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              (alpha_theta (i, j - 1, k) - alpha (i, j, k)) *
     &              f_num_if (i, j - 1, k, 1))) *
     &              delta_theta_mod_inv (j)
#else /* CFC_TRANSPORT2 */
c     improved formulation of Newtonian energy equation
               f_diff (i, j, k, 5) = f_diff (i, j, k, 5) +
     &              (sqrt_gamma_theta (i, j, k) *
     &              (f_num_if (i, j, k, 5) +
     &              (phi_potential (i, j + 1, k) -
     &              phi_potential (i, j, k)) *
     &              f_num_if (i, j, k, 1)) -
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              (f_num_if (i, j - 1, k, 5) +
     &              (phi_potential (i, j - 1, k) -
     &              phi_potential (i, j, k)) *
     &              f_num_if (i, j - 1, k, 1))) *
     &              delta_theta_mod_inv (j)               
#endif /* CFC_TRANSPORT2 */

#else /* ENECONS */
c     standard formulation of energy equation
               f_diff (i, j, k, 5) = f_diff (i, j, k, 5) +
     &             (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 5) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 5)) *
     &              delta_theta_mod_inv (j)
               
#endif /* ENECONS*/
               
#ifdef CFC_MHD
               f_diff (i, j, k, 6) = f_diff (i, j, k, 6) +
     &             (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 6) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 6)) *
     &              delta_theta_mod_inv (j)

               s_total (i, j, k, 7) =
     &              - sqrt_gamma (i, j, k) * g_up_22 (i, j, k) *
     &              (f_num_if (i, j, k, 7) -
     &              f_num_if (i, j - 1, k, 7)) *
     &              delta_theta_inv
               
               f_diff (i, j, k, 8) = f_diff (i, j, k, 8) +
     &             (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 8) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 8)) *
     &              delta_theta_mod_inv (j)
c              if (i.eq.100.and.j.eq.n_e)
c     &              print *,'fdiffth',k,f_num_if (i,j,k,6:8),
c     &              f_num_if (i,j-1,k,6:8)

               f_diff (i, j, k, 9) = f_diff (i, j, k, 9) +
     &             (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 9) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 9)) *
     &              delta_theta_mod_inv (j)

               f_diff (i, j, k, 11) = f_diff (i, j, k, 11) +
     &             (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 11) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 11)) *
     &              delta_theta_mod_inv (j)
c               if (j.eq.n_s.and.k.eq.o_e) print *,'fdiff2',
c     &              i,j,k,f_diff (i, j, k, 9),f_diff (i, j, k, 11)
               
#ifdef CFC_TRANSPORT2
#error
#endif               
c     internal energy equation
               f_diff (i, j, k, 10) = f_diff (i, j, k, 10) +
     &              (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 10) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 10)) *
     &              delta_theta_mod_inv (j) + 0.25_rk *
     &              (alpha_theta (i, j, k) *
     &              sqrt_gamma_theta (i, j, k) *
     &              f_num_if (i, j, k, 1) *
     &              (2*p (i, j, k) + 0*p (i, j + 1, k)) *
     &              (1.0_rk / rho (i, j, k) -
     &              1.0_rk / rho (i, j + 1, k)) -
     &              alpha_theta (i, j - 1, k) *
     &              sqrt_gamma_theta (i, j - 1, k) *
     &              f_num_if (i, j - 1, k, 1) *
     &              (2*p (i, j, k) + 0*p (i, j - 1, k)) *
     &              (1.0_rk / rho (i, j, k) -
     &              1.0_rk / rho (i, j - 1, k))) *
     &              delta_theta_mod_inv (j) * rho (i, j, k)
#endif /* CFC_MHD */               

               
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
               do l = 1, config%qn
                  f_diff_xnu (i, j, k, l) = f_diff_xnu (i, j, k, l) +
     &                (alpha_theta (i, j, k) *
     &                 sqrt_gamma_theta (i, j, k) *
     &                 f_num_xnu_if (i, j, k, l) -
     &                 alpha_theta (i, j - 1, k) *
     &                 sqrt_gamma_theta (i, j - 1, k) *
     &                 f_num_xnu_if (i, j - 1, k, l)) *
     &                 delta_theta_mod_inv (j)                  
               enddo

            enddo
         enddo


      case(3) !varphi-direction

         call avare_hy(3)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk,l)
#endif
            do jk = 1, n_loc * o_loc
            
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               k = k + o_s - 1
               
               do i = 1, m
            
                  f_diff (i, j, k, 1) = f_diff (i, j, k, 1) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 f_num_if (i, j, k, 1) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_if (i, j, k - 1, 1)) *
     &                 delta_varphi_inv

#if 0
c     Eliminate geometric source terms arising from varphi component
c     of velocity (currently not recommended)
                  fxr = f_diff (i, j, k, 2) * sinus_theta (j)
     &              + f_diff (i, j, k, 3) * cosin_theta (j) *
     &                 r_inv (i)
                  fyr = f_diff (i, j, k, 4) * sinus_theta_inv (j) *
     &                 r_inv (i)
                  fzr = f_diff (i, j, k, 2) * cosin_theta (j)
     &                 - f_diff (i, j, k, 3) * sinus_theta (j) *
     &                 r_inv (i)
                  
                  fxl = f_diff (i, j, k - 1, 2) * sinus_theta (j)
     &                 + f_diff (i, j, k - 1, 3) * cosin_theta (j) *
     &                 r_inv (i)
                  fyl = f_diff (i, j, k - 1, 4) * sinus_theta_inv (j) *
     &                 r_inv (i)
                  fzl = f_diff (i, j, k - 1, 2) * cosin_theta (j)
     &              - f_diff (i, j, k - 1, 3) * sinus_theta (j) *
     &                 r_inv (i)
                  
                  f_diff (i, j, k, 2) = f_diff (i, j, k, 2) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 ((fxr * cos_dvarphi - fyr * sin_dvarphi) *
     &                 sinus_theta (j) + fzr * cosin_theta (j)) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 ((fxl * cos_dvarphi + fyl * sin_dvarphi) *
     &              sinus_theta (j) + fzl * cosin_theta (j))) *
     &                 delta_varphi_inv
                  
                  f_diff (i, j, k, 3) = f_diff (i, j, k, 3) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 ((fxr * cos_dvarphi - fyr * sin_dvarphi) *
     &                 cosin_theta (j) - fzr * sinus_theta (j)) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 ((fxl * cos_dvarphi + fyl * sin_dvarphi) *
     &                 cosin_theta (j) - fzl * sinus_theta (j))) *
     &                 delta_varphi_inv * r (i)
                  
                  f_diff (i, j, k, 4) = f_diff (i, j, k, 4) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 (fyr * cos_dvarphi + fxr * sin_dvarphi) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 (fyl * cos_dvarphi - fxl * sin_dvarphi)) *
     &                 delta_varphi_inv * r (i) * sinus_theta (j)
#else               
                  f_diff (i, j, k, 2) = f_diff (i, j, k, 2) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 f_num_if (i, j, k, 2) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_if (i, j, k - 1, 2)) *
     &                 delta_varphi_inv
                  
                  f_diff (i, j, k, 3) = f_diff (i, j, k, 3) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 f_num_if (i, j, k, 3) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_if (i, j, k - 1, 3)) *
     &                 delta_varphi_inv
                  
                  f_diff (i, j, k, 4) = f_diff (i, j, k, 4) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 f_num_if (i, j, k, 4) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_if (i, j, k - 1, 4)) *
     &                 delta_varphi_inv
#endif
         
#if MOMENTUM_EQUATION_MODE>0
                  f_diff (i, j, k, 4) = f_diff (i, j, k, 4) +
     &                 (pav_varphi (i, j, k) - 
     &                 pav_varphi (i, j, k -1)) *
     &                 alpha (i, j, k) * 
     &                 phi (i, j, k) ** 6 *
     &                 delta_varphi_inv
#endif /* MOMENTUM_EQUATION_MODE>0 */

#ifdef ENECONS
#ifdef CFC_TRANSPORT2
                  f_diff (i, j, k, 5) = f_diff (i, j, k, 5) +
     &                 ((alpha_varphi (i, j, k) ** 2 *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 f_num_if (i, j, k, 5) -
     &                 alpha_varphi (i, j, k - 1) ** 2 *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_if (i, j, k - 1, 5)) +
     &                 (alpha_varphi (i, j, k) * 
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 (alpha_varphi (i, j, k) - alpha (i, j, k)) *
     &                 f_num_if (i, j, k, 1) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 (alpha_varphi (i, j, k - 1) - alpha (i, j, k)) *
     &                 f_num_if (i, j, k - 1, 1))) *
     &                 delta_varphi_inv
#else /*CFC_TRANSPORT2 */
                  f_diff (i, j, k, 5) = f_diff (i, j, k, 5) +
     &                 (sqrt_gamma_varphi (i, j, k) *
     &                 (f_num_if (i, j, k, 5) +
     &                 (phi_potential (i, j, k + 1) -
     &                 phi_potential (i, j, k)) *
     &                 f_num_if (i, j, k, 1)) -
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 (f_num_if (i, j, k - 1, 5) +
     &                 (phi_potential (i, j, k - 1) -
     &                 phi_potential (i, j, k)) *
     &                 f_num_if (i, j, k - 1, 1))) *
     &                 delta_varphi_inv
#endif /* CFC_TRANSPORT2 */

#endif /* ENECONS */
                  
#ifdef CFC_MHD
                  f_diff (i, j, k, 6) = f_diff (i, j, k, 6) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 f_num_if (i, j, k, 6) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_if (i, j, k - 1, 6)) *
     &                 delta_varphi_inv
                  
                  f_diff (i, j, k, 7) = f_diff (i, j, k, 7) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 f_num_if (i, j, k, 7) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_if (i, j, k - 1, 7)) *
     &                 delta_varphi_inv
                  
                  s_total (i, j, k, 8) = 
     &                 - sqrt_gamma (i, j, k) * g_up_33 (i, j, k) *
     &                 (f_num_if (i, j, k, 8) -
     &                 f_num_if (i, j, k - 1, 8)) *
     &                 delta_varphi_inv
                  
                  f_diff (i, j, k, 9) = f_diff (i, j, k, 9) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 f_num_if (i, j, k, 9) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_if (i, j, k - 1, 9)) *
     &                 delta_varphi_inv


                  f_diff (i, j, k, 11) = f_diff (i, j, k, 11) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 f_num_if (i, j, k, 11) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_if (i, j, k - 1, 11)) *
     &                 delta_varphi_inv
                  
#ifdef CFC_TRANSPORT2
#error
#endif               
c     internal energy equation                  
                  f_diff (i, j, k, 10) = f_diff (i, j, k, 10) +
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 f_num_if (i, j, k, 10) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_if (i, j, k - 1, 10)) *
     &                 delta_varphi_inv + 0.25_rk * 
     &                 (alpha_varphi (i, j, k) *
     &                 sqrt_gamma_varphi (i, j, k) *
     &                 f_num_if (i, j, k, 1) *
     &                 (2*p (i, j, k) + 0*p (i, j, k + 1)) *
     &                 (1.0_rk / rho (i, j, k) -
     &                 1.0_rk / rho (i, j, k + 1)) -
     &                 alpha_varphi (i, j, k - 1) *
     &                 sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_if (i, j, k - 1, 1) *
     &                 (2*p (i, j, k) + 0*p (i, j, k - 1)) *
     &                 (1.0_rk / rho (i, j, k) -
     &                 1.0_rk / rho (i, j, k - 1))) *
     &                 delta_varphi_inv * rho (i, j, k)                  
#endif /* CFC_MHD */               
               enddo

               do l = 1, config%qn
                  do i = 1, m
                     
                  f_diff_xnu (i, j, k, l) = f_diff_xnu (i, j, k, l) +
     &                    (alpha_varphi (i, j, k) *
     &                    sqrt_gamma_varphi (i, j, k) *
     &                 f_num_xnu_if (i, j, k, l) -
     &                    alpha_varphi (i, j, k - 1) *
     &                    sqrt_gamma_varphi (i, j, k - 1) *
     &                 f_num_xnu_if (i, j, k - 1, l)) *
     &                    delta_varphi_inv
                     
                  enddo
               enddo
               
            enddo               !j,k-loop

#ifdef CFC_MHD
c     rescale divergence in psi-equation and add divergence source in
c     momentum & energy equation
      case(4)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk,l)
#endif
         do jk = 1, n_loc * o_loc
               
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            
            delta_t_inv = 1.0_rk / delta_t
            
            do i = 1, m

               scr1 = 0*f_diff (i, j, k, 9) *
     &              b_cap_1_hat (i, j, k) * phi (i, j, k) ** 2
               scr2 = 0*f_diff (i, j, k, 9) *
     &              b_cap_2_hat (i, j, k) * phi (i, j, k) ** 2 * r (i)
               scr3 = 0*f_diff (i, j, k, 9) *
     &              b_cap_3_hat (i, j, k) * phi (i, j, k) ** 2 * r (i) *
     &              sinus_theta (j)
               f_diff (i, j, k, 2) = f_diff (i, j, k, 2) - scr1
               f_diff (i, j, k, 3) = f_diff (i, j, k, 3) - scr2
               f_diff (i, j, k, 4) = f_diff (i, j, k, 4) - scr3
c               f_diff (i, j, k, 5) = f_diff (i, j, k, 5) - 
c     &             (scr1 * v_1 (i, j, k) * g_up_11 (i, j, k) +
c     &              scr2 * v_2 (i, j, k) * g_up_22 (i, j, k) +
c     &              scr3 * v_3 (i, j, k) * g_up_33 (i, j, k))
c               s_total (i, j, k, 6) = 0.0_rk
c               s_total (i, j, k, 7) = 0.0_rk
c               s_total (i, j, k, 8) = 0.0_rk
               vcln = 0.25_rk * min (r (i) * delta_theta, delta_r (i)) *
     &              delta_t_inv !cleaning speed
               vcln = min (vcln, sqrt ((b_1 (i, j, k) ** 2 +
     &              b_2 (i, j, k) ** 2  + b_3 (i, j, k) ** 2) /
     &              rho (i, j, k) + c_sound_squared (i, j, k)))
c               vcln = 0.0_rk
               f_diff (i, j, k, 9) = f_diff (i, j, k, 9) * vcln
c     &              + f_diff (i, j, k, 11)
               f_diff (i, j, k, 10) = f_diff (i, j, k, 10) -
     &              f_diff (i, j, k, 1) * p (i, j, k)

#if 0
c     add corrections to satisfy internal energy equation
               if (p(i,j,k) .lt. 0.5_rk *
     &              (b_cap_1_hat(i,j,k)**2+b_cap_2_hat(i,j,k)**2+
     &              b_cap_3_hat(i,j,k))**2)
     &              f_diff (i, j, k, 5) = f_diff (i, j, k, 5) -
     &              0.5_rk * delta_t *
     &             ((- f_diff (i, j, k, 6) +
     &              alpha (i, j, k) * s_total (i, j, k, 6)) ** 2 *
     &              g_11 (i, j, k) +               
     &              (- f_diff (i, j, k, 7) +
     &              alpha (i, j, k) * s_total (i, j, k, 7)) ** 2 *
     &              g_22 (i, j, k) +
     &              (- f_diff (i, j, k, 8) +
     &              alpha (i, j, k) * s_total (i, j, k, 8)) ** 2 *
     &              g_33 (i, j, k))
#endif               
#if 0               
               if (p (i, j, k) .lt. 0.15_rk * rho(i, j, k) *
     &              v_squared (i, j, k))
     &              f_diff (i, j, k, 5) = f_diff (i, j, k, 5) -
     &              0.5_rk * delta_t *
     &             ((- f_diff (i, j, k, 2) +
     &              alpha (i, j, k) * s_total (i, j, k, 2)) ** 2 *
     &              g_up_11 (i, j, k) +               
     &              (- f_diff (i, j, k, 3) +
     &              alpha (i, j, k) * s_total (i, j, k, 3)) ** 2 *
     &              g_up_22 (i, j, k) +
     &              (- f_diff (i, j, k, 4) +
     &              alpha (i, j, k) * s_total (i, j, k, 4)) ** 2 *
     &              g_up_33 (i, j, k))
#endif               

            enddo
         enddo
#endif /* CFC_MHD */

      end select

      return

      END SUBROUTINE calculate_flux_differences

c     ==================================================================




c     ==================================================================

      SUBROUTINE bndry_flx(rk_step,rk_order)

c     ==================================================================
c     Integrate fluxes through outer boundary.
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE grid_cfc
      USE fluxes_cfc
      USE metric_cfc
      USE perm_aux_cfc


      USE massio_hy
      USE phycon

      USE configure
      USE hydro_areas_mod
      USE abort

      IMPLICIT NONE

      integer (kind=ik), intent(in) :: rk_step, rk_order

      real (kind=rk) :: alpha_rk,beta_rk
      integer (kind=ik) :: i, j, k, jk

c     Zero arrays:
      if (rk_step .eq. 1) then
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk)
#endif
         do jk = 1, n_loc * o_loc
            
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            dflxtot (2,j,k,areas%are_nu)           = 0.0_rk
            eflxtot (2,j,k,areas%are_nu)           = 0.0_rk
            xnflxtot(2,j,k,config%qn,areas%are_nu) = 0.0_rk

         end do
      end if

c     Set weights for Runge-Kutta stepping:
      if (rk_order .eq. 2) then
         if (rk_step .eq. 1) then
            alpha_rk=1.0_rk
            beta_rk =1.0_rk
         else if (rk_step .eq. 2) then
            alpha_rk=2.0_rk
            beta_rk =1.0_rk
         else
            raise_abort('ERROR in bndry_flx')
         end if
      else
         if (rk_step .eq. 1) then
            alpha_rk=1.0_rk
            beta_rk =1.0_rk
         else if (rk_step .eq. 2) then
            alpha_rk=4.0_rk
            beta_rk =3.0_rk
         else if (rk_step .eq. 3) then
            alpha_rk=3.0_rk
            beta_rk =2.0_rk
         else
            raise_abort('ERROR in bndry_flx')
         end if
      end if

c     Integrate surface fluxes:

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk)
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         dflxtot (2,j,k,areas%are_nu) =(dflxtot(2,j,k,areas%are_nu)    +
     &        beta_rk*
     &        alpha_r(m,j,k)*sqrt_gamma_r(m,j,k)*
     &        f_num_if    (m,j,k,1))/alpha_rk
         eflxtot (2,j,k,areas%are_nu) =(eflxtot(2,j,k,areas%are_nu)    +
     &        beta_rk*
     &        alpha_r(m,j,k)*sqrt_gamma_r(m,j,k)*
     &        f_num_if    (m,j,k,5))/alpha_rk
         xnflxtot(2,j,k,config%qn,areas%are_nu)= 
     &       (xnflxtot(2,j,k,config%qn,areas%are_nu)+
     &        beta_rk*
     &        alpha_r(m,j,k)*sqrt_gamma_r(m,j,k)*
     &        f_num_xnu_if(m,j,k,config%qn))/alpha_rk

      end do


c     Convert to cgs units:
      if (rk_step .eq. rk_order) then

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk)
#endif
         do jk = 1, n_loc * o_loc
         
            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1

            dflxtot (2,j,k,areas%are_nu)   =
     &           dflxtot (2,j,k,areas%are_nu)   *pc_geog*pc_cl
            eflxtot (2,j,k,areas%are_nu)   =
     &           (eflxtot (2,j,k,areas%are_nu)  *pc_geoe*pc_cl
     &           +dflxtot (2,j,k,areas%are_nu)  *        pc_cl**2)*
     &           alpha_r(m,j,k)
            xnflxtot(2,j,k,config%qn,areas%are_nu)=
     &           xnflxtot(2,j,k,config%qn,areas%are_nu)*pc_geog*pc_cl

         end do
      end if

      return

      END SUBROUTINE bndry_flx

c     ==================================================================




c     ==================================================================

      SUBROUTINE avare_hy(mode)

c     ==================================================================

      USE precision


      USE totare_hy

      USE size_cfc
      USE parameters_cfc
      USE grid_cfc
      USE hydro_primitives_cfc
      USE conserved_cfc
      USE fluxes_cfc
      USE sources_cfc
      USE metric_cfc
      USE perm_aux_cfc

#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      USE configure

      USE hydro_areas_mod
      USE abort

      IMPLICIT NONE

      integer (kind=ik), intent(in) :: mode


      real(kind=rk) :: dvy

      integer (kind=ik) :: ixi,ixf,iox,iyi,iyf,ioy,izi,izf,ioz,isd
      integer (kind=ik) :: iarea,i,j,k,jk,nvars

#ifdef MPI_HYDRO
      real (kind=rk), allocatable :: mpibuf(:,:),mpibuf2(:,:)
      real (kind=rk) :: scr1
      integer (kind=ik) :: ierr,istat
#endif

      if (config%nsdim .eq. 1) return

#ifdef CFC_MHD
      nvars = 13+config%qn
#else            
      nvars = 7+config%qn
#endif /* CFC_MHD */               

      
      do iarea=1,areas%are_nu

         ixi    = areas%ix_are(iarea, 1)
         ixf    = areas%ix_are(iarea, 2)
         iox    = areas%ix_are(iarea, 3)
         iyi    = areas%ix_are(iarea, 4)
         iyf    = areas%ix_are(iarea, 5)
         ioy    = areas%ix_are(iarea, 6)
         izi    = areas%ix_are(iarea, 7)
         izf    = areas%ix_are(iarea, 8)
         ioz    = areas%ix_are(iarea, 9)
         isd    = areas%ix_are(iarea,10)

         if (ioy .ne. 1 .and. ioy .eq.  config%qy) then
            select case(mode)

            case(2) 
c     Zero fluxes in theta-direction

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk)
#endif
               do jk = 1, n_loc * o_loc
                  
                  k = int((jk + n_loc - 1) / n_loc )
                  j = (n_s - 1) + (jk - (k - 1) * n_loc)
                  k = k + o_s - 1
         
                  f_num_if      (ixi:ixf,j,k,1:5)=0.0_rk
#ifdef CFC_MHD
                  f_num_if      (ixi:ixf,j,k,10 )=0.0_rk
#endif                  
                  f_num_xnu_if  (ixi:ixf,j,k,:)=0.0_rk
                  s_total          (ixi:ixf,j,k,3)=0.0_rk

                  if (j .eq. n_s) then
                      f_num_if      (ixi:ixf,j-1,k,1:5)=0.0_rk
                      f_num_if      (ixi:ixf,j-1,k,10 )=0.0_rk
                      f_num_xnu_if  (ixi:ixf,j-1,k,:)=0.0_rk
                   end if
               end do

            case(3) 
c     Zero flux in varphi-direction

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk)
#endif
               do jk = 1, n_loc * o_loc
                  
                  k = int((jk + n_loc - 1) / n_loc )
                  j = (n_s - 1) + (jk - (k - 1) * n_loc)
                  k = k + o_s - 1
         
                  f_num_if     (ixi:ixf,j,k,1:5)=0.0_rk
#ifdef CFC_MHD
                  f_num_if     (ixi:ixf,j,k,10 )=0.0_rk
#endif                  
                  f_num_xnu_if (ixi:ixf,j,k,:)=0.0_rk
                  s_total          (ixi:ixf,j,k,4)=0.0_rk

                  if (k .eq. o_s) then
                     f_num_if     (ixi:ixf,j,k-1,1:5)=0.0_rk
#ifdef CFC_MHD
                     f_num_if     (ixi:ixf,j-1,k,10 )=0.0_rk
#endif                  
                     f_num_xnu_if (ixi:ixf,j,k-1,:)=0.0_rk
                  end if
               end do

            case(1) 
c     Zero sources for S_theta and S_varphi

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk)
#endif
               do jk = 1, n_loc * o_loc
                  
                  k = int((jk + n_loc - 1) / n_loc )
                  j = (n_s - 1) + (jk - (k - 1) * n_loc)
                  k = k + o_s - 1
         
                  s_total          (ixi:ixf,j,k,3)=0.0_rk
                  s_total          (ixi:ixf,j,k,4)=0.0_rk

               end do

            case(0)
c     Average conserved variables:
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif

               qtot   (ixi:ixf,:)=0.0_rk
               xtot   (ixi:ixf,:)=0.0_rk
               sumdvy (ixi:ixf)  =0.0_rk
               sumdvy2(ixi:ixf)  =0.0_rk

               do jk = 1, n_loc * o_loc
            
                  k = int((jk + n_loc - 1) / n_loc )
                  j = (n_s - 1) + (jk - (k - 1) * n_loc)
                  k = k + o_s - 1
            
                  dvy=dvytot(j)*dvztot(k)

                  do i=ixi,ixf

#ifdef CFC_MHD
c     Subtract magnetic energy before averaging
c     (STILL NEED TO ADAPT THIS TO GR CASE):                     
                     tau_hat (i, j, k) = tau_hat (i, j, k) - 0.5_rk *
     &                    (b_1 (i, j, k) ** 2 +
     &                    b_2 (i, j, k) ** 2 + 
     &                    b_3 (i, j, k) ** 2) 
#endif /* CFC_MHD */
            
                     qtot(i,1)   =qtot(i,1)   +
     &                    d_cap_hat(i,j,k)   *dvy
                     qtot(i,2)   =qtot(i,2)   +
     &                    s_1_hat(i,j,k)     *dvy
c                     qtot(i,3)   =qtot(i,3)   +
c     &                    s_2_hat(i,j,k)     *dvy
                     qtot(i,4)   =qtot(i,4)   +
     &                    s_3_hat(i,j,k)     *dvy
                     qtot(i,5)   =qtot(i,5)   +
     &                    tau_hat(i,j,k)     *dvy
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
                     xtot(i,1:config%qn)=xtot(i,1:config%qn)+
     &                    d_cap_xnu_hat(i,j,k,1:config%qn)*dvy

#ifdef CFC_MHD
#if DIMN==3
                     scr1 =
     &                    (b_cap_1_hat(i,j,k) * sinus_theta (j) +
     &                    b_cap_2_hat(i,j,k) * cosin_theta (j))
                     qtot(i,6)   =qtot(i,6)   + !B_x
     &                    (scr1 * cosin_varphi (k) -
     &                    b_cap_3_hat(i,j,k) * sinus_varphi (k))*dvy
                     qtot(i,7)   =qtot(i,7)   + !B_x
     &                    (scr1 * sinus_varphi (k) +
     &                    b_cap_3_hat(i,j,k) * cosin_varphi (k))*dvy
#endif /* DIMN==3 */
                     qtot(i,8)   =qtot(i,8)   + !B_z
     &                    (b_cap_1_hat(i,j,k) * cosin_theta (j) -
     &                    b_cap_2_hat(i,j,k) * sinus_theta (j))*dvy
#if DIMN==3                     
                     qtot(i,9)   =qtot(i,9)   + !rot B ~ e_x
     &                    (b_cap_2_hat(i,j,k) * sinus_varphi (k) +
     &                    b_cap_3_hat(i,j,k) * cosin_varphi (k) *
     &                    cosin_theta (j)) * dvy
                     qtot(i,10)   =qtot(i,10)   + !rot B ~ e_y
     &                    (b_cap_2_hat(i,j,k) * cosin_varphi (k) -
     &                    b_cap_3_hat(i,j,k) * sinus_varphi (k) *
     &                    cosin_theta (j)) * dvy
#endif /* DIMN==3 */
                     qtot(i,11)   =qtot(i,11)   + ! rot_B ~ e_z
     &                    b_cap_3_hat(i,j,k ) * sinus_theta (j) * dvy
#endif /* CFC_MHD */

                     sumdvy (i)=sumdvy (i)+dvy*
     &                    sqrt_gamma(i,j,k)
#ifdef CFC_TRANSPORT2
                     sumdvy2(i)=sumdvy2(i)+dvy*
     &                    d_cap_hat(i,j,k)*h(i,j,k)*w(i,j,k)*
     &                    rr(i)*ssinus_theta(j)*phi(i,j,k)**4
#else /* CFC_TRANSPORT2 */
                     sumdvy2(i)=sumdvy2(i)+dvy*
     &                    d_cap_hat(i,j,k)*
     &                    rr(i)*ssinus_theta(j)*phi(i,j,k)**4
#endif /* CFC_TRANSPORT2 */
                     
                  end do
               end do

#ifdef MPI_HYDRO
               allocate(mpibuf (1:ixf+1-ixi,nvars),stat=istat)
               allocate(mpibuf2(1:ixf+1-ixi,nvars),stat=istat)
               mpibuf(1:ixf+1-ixi,1:5   )=qtot  (ixi:ixf,1:5)
#ifdef CFC_MHD
               mpibuf(1:ixf+1-ixi,8+config%qn:13+config%qn)=
     &              qtot  (ixi:ixf,6:11)
#endif /* CFC_MHD */
               mpibuf(1:ixf+1-ixi,  6   )=sumdvy(ixi:ixf)
               mpibuf(1:ixf+1-ixi,7:6+config%qn)=xtot  (ixi:ixf, 
     &                                                  1:config%qn)
               mpibuf(1:ixf+1-ixi,  7+config%qn)=sumdvy2(ixi:ixf)

               call MPI_Allreduce(mpibuf,mpibuf2,(ixf+1_ik-ixi)*
     &              nvars,MPI_DOUBLE_PRECISION,MPI_SUM,
     &              MPI_COMM_WORLD,ierr)
               mpibuf=mpibuf2

               qtot   (ixi:ixf,1:5 )=mpibuf(1:ixf+1-ixi,1:5   )
#ifdef CFC_MHD
               qtot  (ixi:ixf,6:11)=
     &              mpibuf(1:ixf+1-ixi,8+config%qn:13+config%qn)
#endif /* CFC_MHD */
               sumdvy (ixi:ixf     )=mpibuf(1:ixf+1-ixi,  6   )
               xtot   (ixi:ixf,1:config%qn)=
     &                        mpibuf(1:ixf+1-ixi,7:6+config%qn)
               sumdvy2(ixi:ixf     )=mpibuf(1:ixf+1-ixi,7+config%qn)
               deallocate(mpibuf,mpibuf2,stat=istat)
#endif

               do i=ixi,ixf

                  qtot(i,1)=qtot(i,1)/sumdvy(i)
                  qtot(i,2)=qtot(i,2)/sumdvy(i)
                  qtot(i,3)=qtot(i,3)/sumdvy(i)
                  qtot(i,5)=qtot(i,5)/sumdvy(i)
!                  qtot(i,4)=qtot(i,4)/sumdvy(i)
#ifdef CFC_MHD
                  qtot(i, 6)=qtot(i, 6)/sumdvy0(0)
                  qtot(i, 7)=qtot(i, 7)/sumdvy0(0)
                  qtot(i, 8)=qtot(i, 8)/sumdvy0(0)
                  qtot(i, 9)=qtot(i, 9)/sumdvy0(1)
                  qtot(i,10)=qtot(i,10)/sumdvy0(2)
                  qtot(i,11)=qtot(i,11)/sumdvy0(3)
#endif /* CFC_MHD */

#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
                  xtot(i,1:config%qn)=xtot(i,1:config%qn)/sumdvy(i)

               end do

               l_z = SUM (qtot(ixi:ixf,4)*dvxtot(ixi:ixf))
               om_z = l_z / SUM (sumdvy2(ixi:ixf)*dvxtot(ixi:ixf))

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE

!$OMP DO
!$OMP& PRIVATE (i,j,k,jk)
#endif
               do jk = 1, n_loc * o_loc
            
                  k = int((jk + n_loc - 1) / n_loc )
                  j = (n_s - 1) + (jk - (k - 1) * n_loc)
                  k = k + o_s - 1

                  do i=ixi,ixf
                     
                     d_cap_hat    (i,j,k)     =qtot(i,1)*
     &                    sqrt_gamma(i,j,k)
                     s_1_hat      (i,j,k)     =qtot(i,2)*
     &                    sqrt_gamma(i,j,k)
                     s_2_hat      (i,j,k)     =0.0_rk !qtot(i,3)*
c     &                    sqrt_gamma(i,j,k)
c                     s_3_hat      (i,j,k)     =qtot(i,4)!*
c     &                    (sinus_theta(j)**2/sumdvy2(i))*
c     &                    sqrt_gamma(i,j,k)
c                     s_3_hat      (i,j,k)     = qtot(i,4)*
c     &                    (sinus_theta(j)**2/sumdvy2(i))*
c     &                    sqrt_gamma(i,j,k)
                     tau_hat      (i,j,k)     =qtot(i,5)*
     &                    sqrt_gamma(i,j,k)
#ifdef CFC_MHD
                     b_cap_1_hat (i, j, k) =
     &                    (qtot(i,6) * cosin_varphi (k) +
     &                     qtot(i,7) * sinus_varphi (k)) *
     &                    sinus_theta (j) +
     &                    qtot(i, 8) * cosin_theta (j)
                     b_cap_2_hat (i, j, k) =
     &                    (qtot(i,6) * cosin_varphi (k) +
     &                     qtot(i,7) * sinus_varphi (k)) *
     &                    cosin_theta (j) -
     &                    qtot(i, 8) * sinus_theta (j) +
     &                    qtot(i, 9) * sinus_varphi (k) +
     &                    qtot(i,10) * cosin_varphi (k)
                     b_cap_3_hat (i, j, k) =
     &                    -qtot(i,6) * sinus_varphi (k) +
     &                     qtot(i,7) * cosin_varphi (k) +
     &                    (qtot(i,9) * cosin_varphi (k) *
     &                    qtot(i,10) * sinus_varphi (k)) *
     &                    cosin_theta (j) +
     &                    qtot(i,11) * sinus_theta (j)

                     b_1 (i, j, k) = b_cap_1_hat (i, j, k) *
     &                    sqrt_gamma_inv (i, j, k)
                     b_2 (i, j, k) = b_cap_2_hat (i, j, k) *
     &                    sqrt_gamma_inv (i, j, k)
                     b_3 (i, j, k) = b_cap_3_hat (i, j, k) *
     &                    sqrt_gamma_inv (i, j, k)
                     
c     Add magnetic energy again:
                     tau_hat (i, j, k) = tau_hat (i, j, k) + 0.5_rk *
     &                    (b_1 (i, j, k) ** 2 +
     &                    b_2 (i, j, k) ** 2 + 
     &                    b_3 (i, j, k) ** 2) 
#endif /* CFC_MHD */
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
                     d_cap_xnu_hat(i,j,k,1:config%qn)= 
     &                    xtot(i,1:config%qn)*
     &                    sqrt_gamma(i,j,k)

#ifdef CFC_TRANSPORT2
                     s_3_hat      (i,j,k)     = om_z *
     &                    d_cap_hat(i,j,k)*h(i,j,k)*w(i,j,k)*
     &                    rr(i)*ssinus_theta(j)*phi(i,j,k)**4
#else /* CFC_TRANSPORT2 */
                     s_3_hat      (i,j,k)     = om_z *
     &                    d_cap_hat(i,j,k)*
     &                    rr(i)*ssinus_theta(j)*phi(i,j,k)**4
#endif /* CFC_TRANSPORT2 */

                  end do        !i-loop
               end do           !j,k-loop

            end select
         else if (ioy .ne. 1) then
            raise_abort('Error: Y-resolution not supported by CoCoNuT')
         end if

      end do

      return

      END SUBROUTINE avare_hy

c     ==================================================================




c     ==================================================================

      SUBROUTINE burn3d(mode, rk_index, dt_hyd)

c     ==================================================================
c     Driver for nuclear burning routine in sweeps:
c     - Copy hydro variable from CoCoNuT to densty, xnnu, etc. as used
c       in PROMETHEUS
c     - Parallel loop over 1D stencils
c     - Copy updated variables back to CoCoNuT
c     - optionally recover primitive quantities (for mode=1)
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE conserved_cfc, ONLY: d_cap_hat,d_cap_xnu_hat,tau_hat,s_1_hat,
     &     s_2_hat,s_3_hat
      USE hydro_primitives_cfc
      USE metric_cfc
c      USE parameters_cfc

!      USE gfloat_hy, ONLY: dt
      USE intgrs_hy
      USE vnew_hy
      USE phycon

      USE gr_initial_setup
      USE recover_prim_vars
      USE cons_check
      USE burn_mod

      USE configure
      USE state
      USE shock_cfc

      USE eos_sn2, ONLY: lsrolo
#if 0
      USE timmes_network
#endif

      IMPLICIT NONE

      integer (kind=ik), intent(in) :: mode, rk_index
      real (kind=rk), intent(in) :: dt_hyd
      real(kind=rk)              :: burn_self(2), burn_children(2)
      real (kind=rk) :: alpha_rk
      real (kind=rk) :: delta_rho_w2_eps_hat
      real (kind=rk) :: xin(19), xout(19)

      integer (kind=ik) :: i,j,k,jk,i_n
      logical           :: ler


      nzn  = m
      nzn1 = nzn + 1
      nzn2 = nzn + 2
      nzn3 = nzn + 3
      nzn4 = nzn + 4
      nzn5 = nzn + 5
      nzn6 = nzn + 6
      nzn7 = nzn + 7
      nzn8 = nzn + 8
      
      hydro%dt = dt_hyd

!      alpha_rk = 1.0_rk / real(rk_index,kind=rk)

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
#ifdef BURN_2ND_ORDER
c     In this case, burn_3d is called from within a parallel region
c     -> use orphaned OpenMP construct, or some modern compiler
c     might think we want a nested parallel loop.
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#else /* BURN_2ND_ORDER */
c     In this case, we have to open a new parallel region.
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,jk)
!$OMP& SHARED(n_s,o_s,n_loc,o_loc,alpha_rk,dt_hyd,
!$OMP& gammac,gamm,press,p,velx,v_1,v_2,v_3,densty,rho,phi,w,
!$OMP& delta_rho_w2_eps_hat,s_1_hat,s_2_hat,s_3_hat,tau_hat,sqrt_gamma,
!$OMP& temp,t,energy,eps,xnuc,xnnu,d_cap_hat,d_cap_xnu_hat,config,
!$OMP& burn_self, burn_children,m, lshock, ishck)
!$OMP& DEFAULT(none)
!$OMP& SCHEDULE(static)
#endif /* BURN_2ND_ORDER */
#endif /* OPENMP_CFC etc. */
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

c     get arrays:
         do i=1,m

            if (lshock(i,j,k)) then
              ishck(i,j,k) = 1
            else
              ishck(i,j,k) = 0
            end if

#if defined(BURN_SS) || defined(BURN_NETW_NOA) || defined(TAK_RATES) || defined(BURN_NETWORK) || defined(NO_FLASHING_AT_SONIC_POINT)
            gammac(i,j,k) = gamm  (i,j,k)
            press (i,j,k) = p     (i,j,k)*pc_geoe 
            energy(i,j,k) = eps   (i,j,k)*rho(i,j,k)*pc_geoe
            velx  (i,j,k) = v_1   (i,j,k)/phi(i,j,k)**2*pc_cl
#endif
            densty(i,j,k) = rho   (i,j,k)*pc_geog
            temp  (i,j,k) = t     (i,j,k)
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
            xnuc(i,j,k,1:config%qn) = xnnu(i,j,k,1:config%qn)

         end do

#if defined(BURN_NETW_NOA) || defined(BURN_NETWORK)
         do i=2,m-1
            temp(i,j,k) = min (temp(i,j,k),
     &           0.5_rk*(t(i-1,j,k)+t(i+1,j,k)))
         end do
#endif /* BURN_NETW_NOA || BURN_NETWORK */


c    do the burning:
#if 1
         call burn(j,k, burn_self, burn_children)
#else
         do i=1,m
            
c               IF (
            IF (temp(i,j,k) .lt. config%tkok*pc_mevk .and.
     &           densty(i,j,k) .lt. lsrolo .and.
     &           .not. lshock(i,j,k)) then
               xin (18) = xnnu (i,j,k,1) !neutrons
               xin (19) = xnnu (i,j,k,2) !protons
               xin ( 1) = 1d-30 !no hydrogen
               xin ( 2) = xnnu (i,j,k,5)+ 
     &              xnnu (i,j,k,3)+xnnu (i,j,k,4)
                                ! lump deuterons, tritons,
                                ! and 3He in 3He
               xin ( 3) = xnnu (i,j,k,6) !4He
               xin ( 4) = xnnu (i,j,k,7) !12C
               xin ( 5) = xnnu (i,j,k,8) !14N
               xin ( 6) = xnnu (i,j,k,9) !16O
               xin ( 7) = xnnu (i,j,k,10) !20Ne
               xin ( 8) = xnnu (i,j,k,11) !24Mg
               xin ( 9) = xnnu (i,j,k,12) !28Si
               xin (10) = xnnu (i,j,k,13) !32S
               xin (11) = xnnu (i,j,k,14) !36Ar
               xin (12) = xnnu (i,j,k,15) !40Ca
               xin (13) = xnnu (i,j,k,16) !44Ti
               xin (14) = xnnu (i,j,k,17) !48Cr
               xin (15) = xnnu (i,j,k,19) !52Fe
               xin (16) = 0.0_rk !We don't have 54Fe in VERTEX
               xin (17) = xnnu (i,j,k,22) !56Ni
c     &                 +xnnu (i,j,k,18)+xnnu (i,j,k,21)
c     &                 +xnnu (i,j,k,23)+
c     &                 xnnu (i,j,k,24)+xnnu (i,j,k,25)
                  
c                  PRINT *,'burning',i,densty(i,j,k),
c     &                 temp(i,j,k),xnuc(i,j,k,config%qn)
c                  call second_v(tim1)
               CALL burn_driver (dt_hyd,temp(i,j,k),densty(i,j,k),
     &              energy(i,j,k),press(i,j,k),xin(:),
     &              xnuc(i,j,k,config%qn),
     &              xout(:))
c                  PRINT *,xout-xin
c                  call second_v(tim2)
c                  print *,tim2-tim1
c                  if (any(xout.ne.xout)) then
c                     print *,'network solver: NaN'
c                     print *,i,log10(densty(i,j,k)),temp(i,j,k)/1d9,
c     &                    xnuc(i,j,k,config%qn)
c                  end if
               xout(:)=xout(:)*sum(xin)/sum(xout)
               xnuc(i,j,k,1) = xout (18)
               xnuc(i,j,k,2) = xout (19)+xout(1)
               xnuc(i,j,k,3) = 0.0_rk
               xnuc(i,j,k,4) = 0.0_rk
               xnuc(i,j,k,5) = xout (2)
               xnuc(i,j,k,6) = xout (3)
               xnuc(i,j,k,7) = xout (4)
               xnuc(i,j,k,8) = xout (5)
               xnuc(i,j,k,9) = xout (6)
               xnuc(i,j,k,10) = xout (7)
[6~               xnuc(i,j,k,11) = xout (8)
               xnuc(i,j,k,12) = xout (9)
               xnuc(i,j,k,13) = xout (10)
               xnuc(i,j,k,14) = xout (11)
               xnuc(i,j,k,15) = xout (12)
               xnuc(i,j,k,16) = xout (13)
               xnuc(i,j,k,17) = xout (14)
               xnuc(i,j,k,18) = xout (16)
               xnuc(i,j,k,19) = xout (15)
!                  xnuc(i,j,k,20) = 0.0d0
!                  xnuc(i,j,k,21) = 0.0d0
               xnuc(i,j,k,22) = xout (17)
!                  xnuc(i,j,k,23) = 0.0d0
!                  xnuc(i,j,k,24) = 0.0d0
!                  xnuc(i,j,k,25) = 0.0d0
c$$$                  print '("vorher ",i4,25(x,e19.7))',i,
c$$$     &                 xnnu(i,j,k,1:config%qn-1)
c$$$                  print '("nachher",i4,25(x,e19.7))',i,
c$$$     &                 xnuc(i,j,k,1:config%qn-1)
c$$$                  print *,'delta_energy',sum((xnnu(i,j,k,1:config%qn-1)-
c$$$     &                 xnuc(i,j,k,1:config%qn-1))*
c$$$     &                 pc_nuc(1:config%qn-1,3)/pc_nuc(1:config%qn-1,2))/
c$$$     &                 pc_mb*pc_meverg
c$$$                  print *,'delta_energy2',(xnnu(i,j,k,1:config%qn-1)-
c$$$     &                 xnuc(i,j,k,1:config%qn-1))*
c$$$     &                 pc_nuc(1:config%qn-1,3)/pc_nuc(1:config%qn-1,2)/
c$$$     &                 pc_mb*pc_meverg
            end if

         end do
#endif


         do i=1,m
c     update hydro arrays:
            t     (i,j,k) = temp(i,j,k)
#if defined(BURN_SS) || defined(BURN_NETW_NOA) || defined(TAK_RATES) || defined(BURN_NETWORK)
            delta_rho_w2_eps_hat =
     &           w(i,j,k)**2*sqrt_gamma (i,j,k)*
     &           (energy(i,j,k)/pc_geoe-eps(i,j,k)*rho(i,j,k))
            tau_hat (i,j,k) = tau_hat (i,j,k) + delta_rho_w2_eps_hat
            s_1_hat (i,j,k) = s_1_hat (i,j,k) +
     &           v_1(i,j,k)*delta_rho_w2_eps_hat
            s_2_hat (i,j,k) = s_2_hat (i,j,k) +
     &           v_2(i,j,k)*delta_rho_w2_eps_hat
            s_3_hat (i,j,k) = s_3_hat (i,j,k) +
     &           v_3(i,j,k)*delta_rho_w2_eps_hat
            eps   (i,j,k) = energy(i,j,k)/pc_geoe/rho(i,j,k)
#endif
#ifdef NEC_COMPILER
!CDIR EXPAND=25
#endif
c            xnnu(i,j,k,1:config%qn) = xnuc(i,j,k,1:config%qn)

         end do

         do i_n = 1, config%qn
            do i = 1, m
               d_cap_xnu_hat(i,j,k,i_n)=d_cap_xnu_hat(i,j,k,i_n)+
     &              d_cap_hat(i,j,k)*
     &              (xnuc(i,j,k,i_n)-xnnu(i,j,k,i_n))
            end do
         end do

      end do

      if (mode .eq. 1) then
         call recover_primitives
      else if (mode .ne. 0) then
         call initialize_conserved_quantities 
      end if

      return

      END SUBROUTINE burn3d

c     ==================================================================




c     ==================================================================

      SUBROUTINE detect_sh_cfc

c     ==================================================================
c     Provides shock detection algorithm needed for adaptive switching
c     from HLLC to HLL.
c     ------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE hydro_primitives_cfc
      USE shock_cfc

      USE hydro_areas_mod
      USE configure

      IMPLICIT NONE

      real (kind=rk) :: pr, pl
      real (kind=rk), parameter :: epsiln=0.6_rk

      integer (kind=ik) :: i,j,k,jk,jp1,jm1,kp1,km1
#ifdef CFC_MHD
      integer (kind=ik) :: iarea
      integer (kind=ik) :: ixi,ixf,iox,iyi,iyf,ioy,izi,izf,ioz,isd
#endif

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk,jp1,jm1,kp1,km1)
#endif
#if DIMN==1
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1         
#endif
#if DIMN==2
      do jk = 1, (n_loc + 2) * o_loc
         k = int((jk + (n_loc + 2) - 1) / (n_loc + 2) )
         j = (n_s - 2) + (jk - (k - 1) * (n_loc + 2))
         k = k + o_s - 1         
#endif
#if DIMN==3
      do jk = 1, (n_loc + 2) * (o_loc + 2)
         k = int((jk + (n_loc + 2) - 1) / (n_loc + 2) )
         j = (n_s - 2) + (jk - (k - 1) * (n_loc + 2))
         k = k + o_s - 2         
#endif
         
         if (config%nsdim .ge. 2) then
            jp1 = j + 1
            jm1 = j - 1
         endif

         if (config%nsdim .eq. 3) then
            kp1 = k + 1
            km1 = k - 1
         endif
         
         do i=2,m-1
            
#ifdef CFC_MHD
c     less restrictive criterion in MHD case
            ! radial pressure gradient
#if 0
            pr = p(i+1,j,k) + 0.5_rk *
     &           (b_1(i+1,j,k)**2+b_2(i+1,j,k)**2+b_3(i+1,j,k)**2)
            pl = p(i-1,j,k) + 0.5_rk *
     &           (b_1(i-1,j,k)**2+b_2(i-1,j,k)**2+b_3(i-1,j,k)**2)
            lshock(i,j,k)= abs(pr-pl) .gt. epsiln*min(pl, pr)

            if (config%nsdim .ge. 2) then
               ! theta pressure gradient
               pr = p(i,j+1,k) +  0.5_rk *
     &              (b_1(i,j+1,k)**2+b_2(i,j+1,k)**2+b_3(i,j+1,k)**2)
               pl = p(i,j-1,k) +  0.5_rk *
     &              (b_1(i,j-1,k)**2+b_2(i,j-1,k)**2+b_3(i,j-1,k)**2)
               lshock(i,j,k)=lshock(i,j,k).or.
     &              (abs(pr-pl) .gt. epsiln*min(pr,pl))
            endif

            if (config%nsdim .eq. 3) then
               ! phi pressure gradient
               pr = p(i,j,k+1) + 0.5_rk *
     &              (b_1(i,j,k+1)**2+b_2(i,j,k+1)**2+b_3(i,j,k+1)**2)
               pl = p(i,j,k-1) + 0.5_rk *
     &              (b_1(i,j,k-1)**2+b_2(i,j,k-1)**2+b_3(i,j,k-1)**2)
               lshock(i,j,k)=lshock(i,j,k).or.
     &              (abs(pr-pl) .gt. epsiln*min(pr,pl))
            endif
#else
            lshock(i,j,k)=
     &           abs(p(i+1,j,k)-p(i-1,j,k))+
     &           abs(b_1(i+1,j,k)**2-b_1(i-1,j,k)**2)+
     &           abs(b_2(i+1,j,k)**2-b_2(i-1,j,k)**2)+
     &           abs(b_3(i+1,j,k)**2-b_3(i-1,j,k)**2)+
     &           abs(rho(i+1,j,k)*v_1(i+1,j,k)**2 -
     &               rho(i-1,j,k)*v_1(i-1,j,k)**2)+
     &           abs(rho(i+1,j,k)*v_2(i+1,j,k)**2 -
     &               rho(i-1,j,k)*v_2(i-1,j,k)**2)+
     &           abs(rho(i+1,j,k)*v_3(i+1,j,k)**2 -
     &               rho(i-1,j,k)*v_3(i-1,j,k)**2)
     &           .gt.
     &           epsiln*min(p(i+1,j,k),p(i,j,k),p(i-1,j,k))

            if (config%nsdim .ge. 2) then
            ! theta pressure gradient
               lshock(i,j,k)=lshock(i,j,k).or.
     &              (abs(p(i,jp1,k)-p(i,jm1,k))+
     &              abs(b_1(i,jp1,k)**2-b_1(i,jm1,k)**2)+
     &              abs(b_2(i,jp1,k)**2-b_2(i,jm1,k)**2)+
     &              abs(b_3(i,jp1,k)**2-b_3(i,jm1,k)**2)+
     &              abs(rho(i,jp1,k)*v_1(i,jp1,k)**2 -
     &                  rho(i,jm1,k)*v_1(i,jm1,k)**2)+
     &              abs(rho(i,jp1,k)*v_2(i,jp1,k)**2 -
     &                  rho(i,jm1,k)*v_2(i,jm1,k)**2)+
     &              abs(rho(i,jp1,k)*v_3(i,jp1,k)**2 -
     &                  rho(i,jm1,k)*v_3(i,jm1,k)**2)
     &              .gt.
     &              epsiln*min(p(i,jp1,k),p(i,j,k),p(i,jm1,k)))
            endif
            
            if (config%nsdim .eq. 3) then
            ! phi pressure gradient
               lshock(i,j,k)=lshock(i,j,k).or.
     &              (abs(p(i,j,kp1)-p(i,j,k))+
     &              abs(b_1(i,j,kp1)**2-b_1(i,j,k)**2)+
     &              abs(b_2(i,j,kp1)**2-b_2(i,j,k)**2)+
     &              abs(b_3(i,j,kp1)**2-b_3(i,j,k)**2)+
     &              abs(rho(i,j,kp1)*v_1(i,j,kp1)**2 -
     &                  rho(i,j,km1)*v_1(i,j,km1)**2)+
     &              abs(rho(i,j,kp1)*v_2(i,j,kp1)**2 -
     &                  rho(i,j,km1)*v_2(i,j,km1)**2)+
     &              abs(rho(i,j,kp1)*v_3(i,j,kp1)**2 -
     &                  rho(i,j,km1)*v_3(i,j,km1)**2)
     &              .gt.
     &              epsiln*min(p(i,j,kp1),p(i,j,k),p(i,j,km1)))
            endif
#endif

#else /* CFC_MHD*/
c     hydro only
            ! radial pressure gradient
            lshock(i,j,k)=
     &           abs(p(i+1,j,k)-p(i-1,j,k)) .gt.
     &           epsiln*min(p(i+1,j,k),p(i-1,j,k))
            if (config%nsdim .ge. 2) then
               ! theta pressure gradient
               lshock(i,j,k)=lshock(i,j,k).or.
     &              (abs(p(i,jp1,k)-p(i,jm1,k)) .gt.
     &              epsiln*min(p(i,jp1,k),p(i,jm1,k)))
            endif

            if (config%nsdim .eq. 3) then
               ! phi pressure gradient
               lshock(i,j,k)=lshock(i,j,k).or.
     &              (abs(p(i,j,kp1)-p(i,j,km1)) .gt.
     &              epsiln*min(p(i,j,kp1),p(i,j,km1)))
            endif
#endif /* CFC_MHD */
         end do
      end do

#if 0      
#ifdef CFC_MHD
c     We need a more dissipative solver at the interface
c     between areas of different resolution
      do iarea=1,areas%are_nu-1

         ixi    = areas%ix_are(iarea, 1)
         ixf    = areas%ix_are(iarea, 2)
         ioy    = areas%ix_are(iarea, 6)

         if (ioy .ne. 1 .and. ioy .eq.  config%qy) then
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE (i,j,k,jk,jp1,jm1,kp1,km1)
#endif
            do jk = 1, n_loc * o_loc
               k = int((jk + n_loc - 1) / n_loc )
               j = (n_s - 1) + (jk - (k - 1) * n_loc)
               k = k + o_s - 1         
               lshock(ixf:ixf+1,j,k) = .true.
            end do
         end if
      end do
#endif /* CFC_MHD */
#endif
      
      return

      END SUBROUTINE detect_sh_cfc

c     ==================================================================




c     ======================================================================

      SUBROUTINE smooth_entropy_holes

c     ======================================================================
c     Detects and eliminates cells with anomalous entropy values (s<0).
c     ----------------------------------------------------------------------

      USE precision

      USE size_cfc
      USE hydro_primitives_cfc
      USE conserved_cfc
      USE metric_cfc, ONLY: sqrt_gamma, sqrt_gamma_inv
      USE phycon, ONLY: pc_geog, pc_geoe, pc_egeo


      USE eos_sn2
      USE configure

      IMPLICIT NONE

      real :: eos_self(2), eos_children(2)
      real :: den(1:m), tmp(1:m), ek(1:m), ei(1:m), pre(1:m),
     &     e(1:m), s(1:m), gamc(1:m), xn(1:m,config%qn), kap(1:m)

! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real :: ccu(1:m),  cce(1:m),  ccn(1:m),  ccp(1:m)
      real :: dmy1(1:m),  dmy2(1:m,2) !dummys
! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real :: rho_h_ww

      integer :: i,j,k,jk,i_nuc,cntr,index_list(1:m)

      logical :: use_ls_eos(1:m)
      logical :: ler 

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
#endif
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         cntr=0
         do i=2,m-1
c     -------- Gather -------------------------------------------------
c               if (entropy(i,j,k).lt.0.0_rk) then
            if (entropy(i,j,k) .lt. 0.3_rk*
     &           min(entropy(i+1,j,k),entropy(i-1,j,k)) .and.
     &           entropy(i,j,k) .lt. 5.0_rk) then
c     &           .or.
c     &           t(i,j,k) .lt. 0.4_rk *
c     &           max(t(i+1,j,k),t(i-1,j,k))) then
               cntr=cntr+1
               index_list(cntr)=i
               rho(i,j,k) = 0.5_rk*(rho(i+1,j,k)+rho(i-1,j,k))
               den(cntr)   = rho(i,j,k)*pc_geog
               ei (cntr)   = 0.5_rk*rho(i,j,k)*pc_geoe*
     &              (eps(i+1,j,k)+eps(i-1,j,k))
               tmp(cntr)   =  max(t(i+1,j,k),t(i-1,j,k))
               xn (cntr,:) = 0.5*(xnnu(i+1,j,k,:)+xnnu(i-1,j,k,:))
c               tmp(cntr)   = 0.5_rk*
c     &              (t(i+1,j,k)+t(i-1,j,k))
c               xn (cntr,:) = xnnu(i,j,k,:)
c               s  (cntr)   = max(
c     &              0.5_rk*(entropy(i+1,j,k)+entropy(i-1,j,k)),
c     &              0.5_rk*(entropy(i,j-1,k)+entropy(i,j+1,k)),
c     &              0.5_rk*(entropy(i,j,k-1)+entropy(i,j,k+1)))
               s  (cntr)   = max(entropy(i+1,j,k),
     &              entropy(i-1,j,k))
c
               print *,'Smoothing cell',i,j,k,s(cntr),entropy(i,j,k),
     &              den(cntr),tmp(cntr),xn(cntr,config%qn)
            end if
         end do

c     -------- Call EoS with entropy as argument ----------------------
c            call eos(den(1:cntr),tmp(1:cntr),xn(1:cntr,:),
c     &           dmy1(1:cntr),dmy2(1:cntr,:),ei(1:cntr),
c     &           pre(1:cntr),gamc(1:cntr),s(1:cntr),
c     &           ccu(1:cntr),cce(1:cntr),ccn(1:cntr),ccp(1:cntr),
c     &           mode=1,nsemode=0,ler=ler)
         if (cntr .GT. 0) then
            use_ls_eos(1:cntr)=.FALSE.
c            call eos(den(1:cntr),tmp(1:cntr),xn(1:cntr,:),
c     &           dmy1(1:cntr),dmy2(1:cntr,:),ei(1:cntr),
c     &           pre(1:cntr),gamc(1:cntr),s(1:cntr),
c     &           ccu(1:cntr),cce(1:cntr),ccn(1:cntr),ccp(1:cntr),
c     &           eos_self,eos_children,mode=1,nsemode=0,ler=ler)
            call eos(den(1:cntr),tmp(1:cntr),xn(1:cntr,:),
     &           dmy1(1:cntr),dmy2(1:cntr,:),ei(1:cntr),
     &           pre(1:cntr),gamc(1:cntr),s(1:cntr),
     &           ccu(1:cntr),cce(1:cntr),ccn(1:cntr),ccp(1:cntr),
     &              eos_self, eos_children, mode=4,nsemode=0,ler=ler)
         end if
         

c     -------- Scatter ------------------------------------------------
         do i=1,cntr
            print *,'Zone/entropy:',index_list(i),s(i)
            eps(index_list(i),j,k) =
     &           ei(i)*pc_egeo/rho(index_list(i),j,k)         
            p      (index_list(i),j,k) = pre(i)*pc_egeo
            t      (index_list(i),j,k) = tmp(i)
            entropy(index_list(i),j,k) = s  (i)
            gamm   (index_list(i),j,k) = gamc(i)
c     ! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            cpot (index_list(i),j,k,1)  = ccu(i)
            cpot (index_list(i),j,k,2)  = cce(i)
            cpot (index_list(i),j,k,3)  = ccn(i)
            cpot (index_list(i),j,k,4)  = ccp(i)
c     ! Quick fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CFC_TRANSPORT2            
c     Relativistic version
            h(index_list(i),j,k) =
     &           1.0_rk+eps(index_list(i),j,k)+
     &           p(index_list(i),j,k)/rho(index_list(i),j,k)
            c_sound_squared(index_list(i),j,k) =
     &           gamc(i)*p(index_list(i),j,k)/
     &           (rho(index_list(i),j,k)*h(index_list(i),j,k))

            rho_h_ww = rho (index_list(i), j, k) *
     &           h (index_list(i), j, k) *
     &           w (index_list(i), j, k) ** 2
            
            s_1_hat (index_list(i), j, k) =
     &           rho_h_ww * v_1 (index_list (i), j, k) *
     &           sqrt_gamma (index_list(i), j, k)

            s_2_hat (index_list(i), j, k) =
     &           rho_h_ww * v_2 (index_list (i), j, k) *
     &           sqrt_gamma (index_list(i), j, k)

            s_3_hat (index_list(i), j, k) =
     &           rho_h_ww * v_3 (index_list (i), j, k) *
     &           sqrt_gamma (index_list(i), j, k)
            
            tau_hat (index_list(i), j, k) =
     &           (rho_h_ww -
     &           d_cap_hat (index_list (i), j, k) * 
     &           sqrt_gamma_inv (index_list(i), j, k) -
     &           p (index_list(i), j, k))  *
     &           sqrt_gamma (index_list(i), j, k)

            tau_hat (index_list(i), j, k) =
     &           (rho_h_ww -
     &           d_cap_hat (index_list (i), j, k) * 
     &           sqrt_gamma_inv (index_list(i), j, k) -
     &           p (index_list(i), j, k))  *
     &           sqrt_gamma (index_list(i), j, k)
                        
            tau_hat (index_list(i), j, k) =
     &           (rho_h_ww -
     &           d_cap_hat (index_list (i), j, k) * 
     &           sqrt_gamma_inv (index_list(i), j, k) -
     &           p (index_list(i), j, k))  *
     &           sqrt_gamma (index_list(i), j, k)

#if defined(CFC_MHD) &&         !defined(PASSIVE_MHD)
            STOP 'smooth_entropy_holes() does not work yet in GRMHD.'
#endif /* CFC_MHD and not PASSIVE_MHD */
         
#else /* CFC_TRANSPORT2 */
c     Newtonian version
            h(index_list(i),j,k) =
     &           eps(index_list(i),j,k)+
     &           p(index_list(i),j,k)/rho(index_list(i),j,k)+
     &           0.5_rk*v_squared(index_list(i),j,k)
            c_sound_squared(index_list(i),j,k) =
     &           gamc(i)*p(index_list(i),j,k)/
     &           (rho(index_list(i),j,k))
            
c            s_1_hat (index_list(i), j, k) =
c     &           rho (index_list(i), j, k) * v_1 (index_list(i), j, k)
c            s_2_hat (index_list(i), j, k) =
c     &           rho (index_list(i), j, k) * v_2 (index_list(i), j, k)
c            s_3_hat (index_list(i), j, k) =
c     &           rho (index_list(i), j, k) * v_3 (index_list(i), j, k)
            tau_hat (index_list(i), j, k) =
     &           rho (index_list(i), j, k) * h (index_list(i), j, k) -
     &           p (index_list(i), j, k)
#if defined(CFC_MHD) && !defined(PASSIVE_MHD)
     &           + 0.5_rk *
     &           (b_1 (index_list(i), j, k) ** 2 +
     &           b_2 (index_list(i), j, k) ** 2 +
     &           b_3 (index_list(i), j, k) ** 2)
#ifndef CONSTRAINED
     &           + 0.5_rk * psi (index_list (i), j, k) ** 2
#endif
#endif /* CFC_MHD and not PASSIVE_MHD */
            
#endif /* CFC_TRANSPORT2 */

        
         end do
         
         do i_nuc = 1, config%qn
            do i=1,cntr
               xnnu (index_list(i), j, k, i_nuc)= xn (i,i_nuc)
               d_cap_xnu_hat (index_list(i), j, k, i_nuc) =
     &              d_cap_hat (index_list(i), j, k) * xn (i, i_nuc)
            end do
         end do
         
      end do
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END DO
#endif

      return

      END SUBROUTINE smooth_entropy_holes

c     ======================================================================




c     ==================================================================

      SUBROUTINE first_touch_coconut

c     ==================================================================

c      USE precision

      USE size_cfc

      USE conserved_cfc
      USE enecons_cfc
      USE fluxes_cfc
      USE hydro_primitives_cfc
      USE interface_hydro_cfc
      USE shock_cfc
      USE metric_cfc
      USE nutra_cfc
      USE sources_cfc
      USE temp_aux_cfc
#ifndef CFC_TRANSPORT2
      USE modpot_cfc
#endif

      USE configure

      IMPLICIT NONE

      integer :: i,j,k,jk,i_n 

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,jk,i_n)
#endif
      do jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
         
         d_cap_hat(:,j,k)=0.0_rk
         s_1_hat  (:,j,k)=0.0_rk
         s_2_hat  (:,j,k)=0.0_rk
         s_3_hat  (:,j,k)=0.0_rk
#ifdef CFC_MHD
         b_cap_1_hat  (:,j,k)=0.0_rk
         b_cap_2_hat  (:,j,k)=0.0_rk
         b_cap_3_hat  (:,j,k)=0.0_rk
#endif /* CFC_MHD */
         tau_hat  (:,j,k)=0.0_rk
         d_cap_xnu_hat(:,j,k,:)=0.0_rk
         
         d_cap_hat_old(:,j,k)=0.0_rk
         s_1_hat_old  (:,j,k)=0.0_rk
         s_2_hat_old  (:,j,k)=0.0_rk
         s_3_hat_old  (:,j,k)=0.0_rk
#ifdef CFC_MHD
         b_cap_1_hat_old  (:,j,k)=0.0_rk
         b_cap_2_hat_old  (:,j,k)=0.0_rk
         b_cap_3_hat_old  (:,j,k)=0.0_rk
#endif /* CFC_MHD */
         tau_hat_old  (:,j,k)=0.0_rk
         d_cap_xnu_hat_old(:,j,k,:)=0.0_rk
         

         phi_potential_alt    (:,j,k)=0.0_rk
         alpha_alt(:,j,k)=0.0_rk

         pav_r     (:,j,k)=0.0_rk
         
         f_diff       (:,j,k,:)=0.0_rk
         
         rho       (:,j,k)=0.0_rk
         v_1       (:,j,k)=0.0_rk
         v_2       (:,j,k)=0.0_rk
         v_3       (:,j,k)=0.0_rk
#ifdef CFC_MHD
         b_1       (:,j,k)=0.0_rk
         b_2       (:,j,k)=0.0_rk
         b_3       (:,j,k)=0.0_rk
#endif /*CFC_MHD */
         eps       (:,j,k)=0.0_rk
         p         (:,j,k)=0.0_rk
         v_squared (:,j,k)=0.0_rk
         w         (:,j,k)=0.0_rk
         h         (:,j,k)=0.0_rk
         t         (:,j,k)=0.0_rk
         entropy   (:,j,k)=0.0_rk
         gamm      (:,j,k)=0.0_rk
         cpot      (:,j,k,:)=0.0_rk
         xnnu      (:,j,k,:)=0.0_rk
         c_sound_squared (:,j,k)=0.0_rk
         enu (:,j,k)=0.0_rk
         fnu (:,j,k)=0.0_rk
         pnu (:,j,k)=0.0_rk
#ifndef CFC_TRANSPORT2
         rho_grav (:,j,k)=0.0_rk
#endif

         if (j. le. n_s+3) then
            rho       (:,2*n_s-1-j,k)=0.0_rk
            v_1       (:,2*n_s-1-j,k)=0.0_rk
            v_2       (:,2*n_s-1-j,k)=0.0_rk
            v_3       (:,2*n_s-1-j,k)=0.0_rk
#ifdef CFC_MHD
            b_1       (:,2*n_s-1-j,k)=0.0_rk
            b_2       (:,2*n_s-1-j,k)=0.0_rk
            b_3       (:,2*n_s-1-j,k)=0.0_rk
#endif /* CFC_MHD */
            eps       (:,2*n_s-1-j,k)=0.0_rk
            p         (:,2*n_s-1-j,k)=0.0_rk
            v_squared (:,2*n_s-1-j,k)=0.0_rk
            w         (:,2*n_s-1-j,k)=0.0_rk
            h         (:,2*n_s-1-j,k)=0.0_rk
            t         (:,2*n_s-1-j,k)=0.0_rk
            gamm      (:,2*n_s-1-j,k)=0.0_rk
            xnnu      (:,2*n_s-1-j,k,:)=0.0_rk
            c_sound_squared (:,2*n_s-1-j,k)=0.0_rk
         end if
         if (j. ge. n-3) then
            rho       (:,2*n_e+1-j,k)=0.0_rk
            v_1       (:,2*n_e+1-j,k)=0.0_rk
            v_2       (:,2*n_e+1-j,k)=0.0_rk
            v_3       (:,2*n_e+1-j,k)=0.0_rk
#ifdef CFC_MHD
            b_1       (:,2*n_e+1-j,k)=0.0_rk
            b_2       (:,2*n_e+1-j,k)=0.0_rk
            b_3       (:,2*n_e+1-j,k)=0.0_rk
#endif /* CFC_MHD */
            eps       (:,2*n_e+1-j,k)=0.0_rk
            p         (:,2*n_e+1-j,k)=0.0_rk
            v_squared (:,2*n_e+1-j,k)=0.0_rk
            w         (:,2*n_e+1-j,k)=0.0_rk
            h         (:,2*n_e+1-j,k)=0.0_rk
            t         (:,2*n_e+1-j,k)=0.0_rk
            gamm      (:,2*n_e+1-j,k)=0.0_rk
            xnnu      (:,2*n_e+1-j,k,:)=0.0_rk
            c_sound_squared (:,2*n_e+1-j,k)=0.0_rk
         end if
         if (k. le. o_s+3) then
            rho       (:,j,2*o_s-1-k)=0.0_rk
            v_1       (:,j,2*o_s-1-k)=0.0_rk
            v_2       (:,j,2*o_s-1-k)=0.0_rk
            v_3       (:,j,2*o_s-1-k)=0.0_rk
#ifdef CFC_MHD
            b_1       (:,j,2*o_s-1-k)=0.0_rk
            b_2       (:,j,2*o_s-1-k)=0.0_rk
            b_3       (:,j,2*o_s-1-k)=0.0_rk
#endif /* CFC_MHD */
            eps       (:,j,2*o_s-1-k)=0.0_rk
            p         (:,j,2*o_s-1-k)=0.0_rk
            v_squared (:,j,2*o_s-1-k)=0.0_rk
            w         (:,j,2*o_s-1-k)=0.0_rk
            h         (:,j,2*o_s-1-k)=0.0_rk
            t         (:,j,2*o_s-1-k)=0.0_rk
            gamm      (:,j,2*o_s-1-k)=0.0_rk
            xnnu      (:,j,2*o_s-1-k,:)=0.0_rk
            c_sound_squared (:,j,2*o_s-1-k)=0.0_rk
         end if
         if (k. ge. o-3) then
            rho       (:,j,2*o_e+1-k)=0.0_rk
            v_1       (:,j,2*o_e+1-k)=0.0_rk
            v_2       (:,j,2*o_e+1-k)=0.0_rk
            v_3       (:,j,2*o_e+1-k)=0.0_rk
#ifdef CFC_MHD
            b_1       (:,j,2*o_e+1-k)=0.0_rk
            b_2       (:,j,2*o_e+1-k)=0.0_rk
            b_3       (:,j,2*o_e+1-k)=0.0_rk
#endif /* CFC_MHD */
            eps       (:,j,2*o_e+1-k)=0.0_rk
            p         (:,j,2*o_e+1-k)=0.0_rk
            v_squared (:,j,2*o_e+1-k)=0.0_rk
            w         (:,j,2*o_e+1-k)=0.0_rk
            h         (:,j,2*o_e+1-k)=0.0_rk
            t         (:,j,2*o_e+1-k)=0.0_rk
            gamm      (:,j,2*o_e+1-k)=0.0_rk
            xnnu      (:,j,2*o_e+1-k,:)=0.0_rk
            c_sound_squared (:,j,2*o_e+1-k)=0.0_rk
         end if


         lshock (:,j,k)=.false.
         if (j.eq.n_s) lshock (:,j-1,:)=.false.
         if (j.eq.n_e) lshock (:,j+1,:)=.false.
         if (k.eq.o_s) lshock (:,:,k-1)=.false.
         if (k.eq.o_e) lshock (:,:,k+1)=.false.
         
         
         rho_if (:,j,k,:)=0.0_rk
         v_1_if (:,j,k,:)=0.0_rk
         v_2_if (:,j,k,:)=0.0_rk
         v_3_if (:,j,k,:)=0.0_rk
#ifdef CFC_MHD
         b_1_if (:,j,k,:)=0.0_rk
         b_2_if (:,j,k,:)=0.0_rk
         b_3_if (:,j,k,:)=0.0_rk
#endif /*CFC_MHD */
         eps_if (:,j,k,:)=0.0_rk
         xnu_if (:,j,k,:,:)=0.0_rk
         p_if   (:,j,k,:)=0.0_rk
         v_squared_if(:,j,k,:)=0.0_rk
         h_if   (:,j,k,:)=0.0_rk
         gamm_if(:,j,k,:)=0.0_rk
         c_sound_squared_if (:,j,k,:)=0.0_rk
         
         if (j. le. n_s+1) then
            rho_if (:,2*n_s-1-j,k,:)=0.0_rk
            v_1_if (:,2*n_s-1-j,k,:)=0.0_rk
            v_2_if (:,2*n_s-1-j,k,:)=0.0_rk
            v_3_if (:,2*n_s-1-j,k,:)=0.0_rk
#ifdef CFC_MHD
            b_1_if (:,2*n_s-1-j,k,:)=0.0_rk
            b_2_if (:,2*n_s-1-j,k,:)=0.0_rk
            b_3_if (:,2*n_s-1-j,k,:)=0.0_rk
#endif /* CFC_MHD */
            eps_if (:,2*n_s-1-j,k,:)=0.0_rk
            xnu_if (:,2*n_s-1-j,k,:,:)=0.0_rk
            p_if   (:,2*n_s-1-j,k,:)=0.0_rk
            v_squared_if(:,2*n_s-1-j,k,:)=0.0_rk
            h_if   (:,2*n_s-1-j,k,:)=0.0_rk
            gamm_if(:,2*n_s-1-j,k,:)=0.0_rk
            c_sound_squared_if (:,2*n_s-1-j,k,:)=0.0_rk
         end if
         if (j. ge. n_e-1) then
            rho_if (:,2*n_e+1-j,k,:)=0.0_rk
            v_1_if (:,2*n_e+1-j,k,:)=0.0_rk
            v_2_if (:,2*n_e+1-j,k,:)=0.0_rk
            v_3_if (:,2*n_e+1-j,k,:)=0.0_rk
#ifdef CFC_MHD
            b_1_if (:,2*n_e+1-j,k,:)=0.0_rk
            b_2_if (:,2*n_e+1-j,k,:)=0.0_rk
            b_3_if (:,2*n_e+1-j,k,:)=0.0_rk
#endif /* CFC_MHD */
            eps_if (:,2*n_e+1-j,k,:)=0.0_rk
            xnu_if (:,2*n_e+1-j,k,:,:)=0.0_rk
            p_if   (:,2*n_e+1-j,k,:)=0.0_rk
            v_squared_if(:,2*n_e+1-j,k,:)=0.0_rk
            h_if   (:,2*n_e+1-j,k,:)=0.0_rk
            gamm_if(:,2*n_e+1-j,k,:)=0.0_rk
            c_sound_squared_if (:,2*n_e+1-j,k,:)=0.0_rk
         end if
         if (k. le. o_s+1) then
            rho_if (:,j,2*o_s-1-j,:)=0.0_rk
            v_1_if (:,j,2*o_s-1-j,:)=0.0_rk
            v_2_if (:,j,2*o_s-1-j,:)=0.0_rk
            v_3_if (:,j,2*o_s-1-j,:)=0.0_rk
#ifdef CFC_MHD
            b_1_if (:,j,2*o_s-1-j,:)=0.0_rk
            b_2_if (:,j,2*o_s-1-j,:)=0.0_rk
            b_3_if (:,j,2*o_s-1-j,:)=0.0_rk
#endif /* CFC_MHD */
            eps_if (:,j,2*o_s-1-j,:)=0.0_rk
            xnu_if (:,j,2*o_s-1-j,:,:)=0.0_rk
            p_if   (:,j,2*o_s-1-j,:)=0.0_rk
            v_squared_if(:,j,2*o_s-1-j,:)=0.0_rk
            h_if   (:,j,2*o_s-1-j,:)=0.0_rk
            gamm_if(:,j,2*o_s-1-j,:)=0.0_rk
            c_sound_squared_if (:,j,2*o_s-1-j,:)=0.0_rk
         end if
         if (k. ge. o_e-1) then
            rho_if (:,j,2*o_e+1-j,:)=0.0_rk
            v_1_if (:,j,2*o_e+1-j,:)=0.0_rk
            v_2_if (:,j,2*o_e+1-j,:)=0.0_rk
            v_3_if (:,j,2*o_e+1-j,:)=0.0_rk
#ifdef CFC_MHD
            b_1_if (:,j,2*o_e+1-j,:)=0.0_rk
            b_2_if (:,j,2*o_e+1-j,:)=0.0_rk
            b_3_if (:,j,2*o_e+1-j,:)=0.0_rk
#endif /* CFC_MHD */
            eps_if (:,j,2*o_e+1-j,:)=0.0_rk
            xnu_if (:,j,2*o_e+1-j,:,:)=0.0_rk
            p_if   (:,j,2*o_e+1-j,:)=0.0_rk
            v_squared_if(:,j,2*o_e+1-j,:)=0.0_rk
            h_if   (:,j,2*o_e+1-j,:)=0.0_rk
            gamm_if(:,j,2*o_e+1-j,:)=0.0_rk
            c_sound_squared_if (:,j,2*o_e+1-j,:)=0.0_rk
         end if

         
         rho_gw      (:,j,k)=0.0_rk
         v_gw_r      (:,j,k)=0.0_rk
         v_gw_theta  (:,j,k)=0.0_rk
         v_gw_varphi (:,j,k)=0.0_rk
         

         s_total(:,j,k,:)=0.0_rk
         dln_sqrt_gamma_dr(:,j,k)=0.0_rk
         dln_sqrt_gamma_dtheta(:,j,k)=0.0_rk
         dln_sqrt_gamma_dvarphi(:,j,k)=0.0_rk
         beta_up_k_k(:,j,k)=0.0_rk
         dphi_dt_o_phi(:,j,k)=0.0_rk
         alpha_dalpha_dr(:,j,k)=0.0_rk
         alpha_dalpha_dtheta(:,j,k)=0.0_rk
         alpha_dalpha_dvarphi(:,j,k)=0.0_rk
         dalpha_dr(:,j,k)=0.0_rk
         dalpha_dtheta(:,j,k)=0.0_rk
         dalpha_dvarphi(:,j,k)=0.0_rk
         dg_00_dr(:,j,k)=0.0_rk
         dg_11_dr(:,j,k)=0.0_rk
         dg_22_dr(:,j,k)=0.0_rk
         dg_33_dr(:,j,k)=0.0_rk
         dg_00_dtheta(:,j,k)=0.0_rk
         dg_11_dtheta(:,j,k)=0.0_rk
         dg_22_dtheta(:,j,k)=0.0_rk
         dg_33_dtheta(:,j,k)=0.0_rk
         dg_00_dvarphi(:,j,k)=0.0_rk
         dg_11_dvarphi(:,j,k)=0.0_rk
         dg_22_dvarphi(:,j,k)=0.0_rk
         dg_33_dvarphi(:,j,k)=0.0_rk
         dg_01_dr(:,j,k)=0.0_rk
         dg_02_dr(:,j,k)=0.0_rk
         dg_03_dr(:,j,k)=0.0_rk
         dg_01_dtheta(:,j,k)=0.0_rk
         dg_02_dtheta(:,j,k)=0.0_rk
         dg_03_dtheta(:,j,k)=0.0_rk
         dg_01_dvarphi(:,j,k)=0.0_rk
         dg_02_dvarphi(:,j,k)=0.0_rk
         dg_03_dvarphi(:,j,k)=0.0_rk
         excrv(:,j,k,:,:)=0.0_rk
         excrv2(:,j,k)=0.0_rk
         excrv_beta(:,j,k)=0.0_rk
         

         phi       (:,j,k)=1.0_rk
         alpha     (:,j,k)=1.0_rk
         beta_up_1 (:,j,k)=0.0_rk
         beta_up_2 (:,j,k)=0.0_rk
         beta_up_3 (:,j,k)=0.0_rk
         g_11 (:,j,k)=0.0_rk
         g_22 (:,j,k)=0.0_rk
         g_33 (:,j,k)=0.0_rk
         g_up_11 (:,j,k)=0.0_rk
         g_up_22 (:,j,k)=0.0_rk
         g_up_33 (:,j,k)=0.0_rk
         phi_r       (:,j,k)=0.0_rk
         alpha_r     (:,j,k)=0.0_rk
         beta_up_1_r (:,j,k)=0.0_rk
         beta_up_2_r (:,j,k)=0.0_rk
         beta_up_3_r (:,j,k)=0.0_rk
         g_11_r (:,j,k)=0.0_rk
         g_22_r (:,j,k)=0.0_rk
         g_33_r (:,j,k)=0.0_rk
         g_up_11_r (:,j,k)=0.0_rk
         g_up_22_r (:,j,k)=0.0_rk
         g_up_33_r (:,j,k)=0.0_rk
         phi_theta       (:,j,k)=0.0_rk
         alpha_theta     (:,j,k)=0.0_rk
         beta_up_1_theta (:,j,k)=0.0_rk
         beta_up_2_theta (:,j,k)=0.0_rk
         beta_up_3_theta (:,j,k)=0.0_rk
         g_11_theta (:,j,k)=0.0_rk
         g_22_theta (:,j,k)=0.0_rk
         g_33_theta (:,j,k)=0.0_rk
         g_up_11_theta (:,j,k)=0.0_rk
         g_up_22_theta (:,j,k)=0.0_rk
         g_up_33_theta (:,j,k)=0.0_rk
         phi_varphi       (:,j,k)=0.0_rk
         alpha_varphi     (:,j,k)=0.0_rk
         beta_up_1_varphi (:,j,k)=0.0_rk
         beta_up_2_varphi (:,j,k)=0.0_rk
         beta_up_3_varphi (:,j,k)=0.0_rk
         g_11_varphi (:,j,k)=0.0_rk
         g_22_varphi (:,j,k)=0.0_rk
         g_33_varphi (:,j,k)=0.0_rk
         g_up_11_varphi (:,j,k)=0.0_rk
         g_up_22_varphi (:,j,k)=0.0_rk
         g_up_33_varphi (:,j,k)=0.0_rk
         k_cap_hat (:,j,k,:,:)=0.0_rk
         k_cap_hat_squared (:,j,k)=0.0_rk
         sqrt_gamma (:,j,k)=1.0_rk
         sqrt_gamma_inv (:,j,k)=1.0_rk
         sqrt_gamma_r (:,j,k)=0.0_rk
         sqrt_gamma_theta (:,j,k)=0.0_rk
         sqrt_gamma_varphi (:,j,k)=0.0_rk
         w_tmp (:,j,k,:)=0.0_rk
         dw_tmp_dr (:,j,k,:)=0.0_rk
         dw_tmp_dtheta (:,j,k,:)=0.0_rk
         s_cap_star      (:,j,k,:)=0.0_rk
         s_cap           (:,j,k,:)=0.0_rk
         rho_hat_phi     (:,j,k)=0.0_rk
         rho_hat_phi_1   (:,j,k)=0.0_rk
         rho_hat_phi_2   (:,j,k)=0.0_rk
         rho_hat_alpha   (:,j,k)=0.0_rk
         rho_hat_alpha_1 (:,j,k)=0.0_rk
#ifndef CFC_TRANSPORT2
         phi_potential (:,j,k)=0.0_rk
         dphi_potential_dr (:,j,k)=0.0_rk
         dphi_potential_dtheta (:,j,k)=0.0_rk
         dphi_potential_dvarphi (:,j,k)=0.0_rk
         dphi_potential_alt_dr (:,j,k)=0.0_rk
         dphi_potential_alt_dtheta (:,j,k)=0.0_rk
         dphi_potential_alt_dvarphi (:,j,k)=0.0_rk
#endif /* CFC_TRANSPORT2 */


         if (j.eq.n_s) then
            phi       (:,j-1,k)=0.0_rk
            alpha     (:,j-1,k)=0.0_rk
            beta_up_1 (:,j-1,k)=0.0_rk
            beta_up_2 (:,j-1,k)=0.0_rk
            beta_up_3 (:,j-1,k)=0.0_rk
            sqrt_gamma (:,j-1,k)=0.0_rk
            sqrt_gamma_inv (:,j-1,k)=0.0_rk
            w_tmp (:,j-1,k,:)=0.0_rk
         end if
         if (j.eq.n_e) then
            phi       (:,j+1,k)=0.0_rk
            alpha     (:,j+1,k)=0.0_rk
            beta_up_1 (:,j+1,k)=0.0_rk
            beta_up_2 (:,j+1,k)=0.0_rk
            beta_up_3 (:,j+1,k)=0.0_rk
            sqrt_gamma (:,j+1,k)=0.0_rk
            sqrt_gamma_inv (:,j+1,k)=0.0_rk
            w_tmp (:,j+1,k,:)=0.0_rk
         endif
         if (k.eq.o_s) then
            phi       (:,j,k-1)=0.0_rk
            alpha     (:,j,k-1)=0.0_rk
            beta_up_1 (:,j,k-1)=0.0_rk
            beta_up_2 (:,j,k-1)=0.0_rk
            beta_up_3 (:,j,k-1)=0.0_rk
            sqrt_gamma (:,j,k-1)=0.0_rk
            sqrt_gamma_inv (:,j,k-1)=0.0_rk
            w_tmp (:,j,k-1,:)=0.0_rk
         end if
         if (k.eq.o_e) then
            phi       (:,j,k+1)=0.0_rk
            alpha     (:,j,k+1)=0.0_rk
            beta_up_1 (:,j,k+1)=0.0_rk
            beta_up_2 (:,j,k+1)=0.0_rk
            beta_up_3 (:,j,k+1)=0.0_rk
            sqrt_gamma (:,j,k+1)=0.0_rk
            sqrt_gamma_inv (:,j,k+1)=0.0_rk
            w_tmp (:,j,k+1,:)=0.0_rk
         endif         

         if (j.eq.n_s) then
            phi_theta       (:,j-1,k)=0.0_rk
            alpha_theta     (:,j-1,k)=0.0_rk
            beta_up_1_theta (:,j-1,k)=0.0_rk
            beta_up_2_theta (:,j-1,k)=0.0_rk
            beta_up_3_theta (:,j-1,k)=0.0_rk
            g_11_theta (:,j-1,k)=0.0_rk
            g_22_theta (:,j-1,k)=0.0_rk
            g_33_theta (:,j-1,k)=0.0_rk
            g_up_11_theta (:,j-1,k)=0.0_rk
            g_up_22_theta (:,j-1,k)=0.0_rk
            g_up_33_theta (:,j-1,k)=0.0_rk
            sqrt_gamma_theta (:,j-1,k)=0.0_rk
         endif
         if (k.eq.o_s) then
            phi_varphi       (:,j,k-1)=0.0_rk
            alpha_varphi     (:,j,k-1)=0.0_rk
            beta_up_1_varphi (:,j,k-1)=0.0_rk
            beta_up_2_varphi (:,j,k-1)=0.0_rk
            beta_up_3_varphi (:,j,k-1)=0.0_rk
            g_11_varphi (:,j,k-1)=0.0_rk
            g_22_varphi (:,j,k-1)=0.0_rk
            g_33_varphi (:,j,k-1)=0.0_rk
            g_up_11_varphi (:,j,k-1)=0.0_rk
            g_up_22_varphi (:,j,k-1)=0.0_rk
            g_up_33_varphi (:,j,k-1)=0.0_rk
            sqrt_gamma_varphi (:,j,k-1)=0.0_rk
         endif

      end do
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END PARALLEL DO
#endif


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,jk,i_n)
#endif
      do jk = 1, (n_loc + 1) * o_loc
         k = int((jk + (n_loc + 1) - 1) / (n_loc + 1))
         j = (n_s - 2) + (jk - (k - 1) * (n_loc + 1))
         k = k + o_s - 1

         pav_theta (:,j,k)=0.0_rk

      end do

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,jk,i_n)
#endif
      do jk = 1, n_loc * (o_loc + 1)
         k = int((jk + n_loc - 1) / n_loc)
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 2

         pav_varphi (:,j,k)=0.0_rk

      end do

      END SUBROUTINE first_touch_coconut

c     ==================================================================




c     ==================================================================

      SUBROUTINE fast_neutrino (bounce)

c     ==================================================================

#define BLOCKING
#define STRANGE
#define HOROWITZ
#define WEAKMAG
      
      USE precision
      USE machcons
      USE grid_cfc
#ifdef CFC_MHD
      USE hydro_primitives_cfc, psi_mag => psi
#else
      USE hydro_primitives_cfc
#endif /* CFC_MHD */
      USE metric_cfc
      USE perm_aux_cfc, ONLY: r_inv, delta_r_if_inv
      USE size_cfc
      USE sources_cfc, ONLY: dalpha_dr
      USE nutrio_hy, ONLY: qentot,qmotot,qmytot,qyetot,fnutot,enutot,
     &     dnutot  
      USE phycon
      USE nucparam
      USE configure
#ifndef CFC_TRANSPORT2
      USE totgrq_hy, ONLY: ephtot
#endif /* CFC_TRANSPORT2 */
      
      USE eos_sn2
      USE mod_degeneracy_parameter, ONLY: etapike

      IMPLICIT NONE

      LOGICAL, INTENT (IN) :: bounce

      INTEGER (KIND=ik), PARAMETER :: nemid = 21!12 !20 !12
      REAL (KIND=rk), PARAMETER :: eran_1 = 4.0_rk,
     &     eran_max = 240.0_rk
      
#ifdef STRANGE      
      REAL (KIND=rk), PARAMETER :: ga_s = -0.1_rk ! nucleon strangeness
#else
      REAL (KIND=rk), PARAMETER :: ga_s = 0.0_rk ! nucleon strangeness
#endif /* STRANGE */

      REAL (KIND=rk) :: emid0 (1:nemid), eran (0:nemid),  emid (1:m),
     &     epl (1:m), emn (1:m) 
      REAL (KIND=rk) :: kapa (0:m+1,3), kaps (0:m+1,3), kapeff (0:m+1,3)
      REAL (KIND=rk) :: dare (0:m), dvol (1:m)
      REAL (KIND=rk) :: f_eq (0:m+1,3), f_nu (0:m+1,0:1,3)
      REAL (KIND=rk) :: f_k (1:m,3), f_h (0:m,3), f_h_in (0:m,3), 
     &     tau (1:m,3)
      REAL (KIND=rk) :: psi (0:m), chi (0:m)
      REAL (KIND=rk) :: scr1 (3) 
      REAL (KIND=rk) :: lum (0:m,3), lumtot (0:m,3)

      REAL (kind=rk) :: tns (3) ! temperature of the neutrinosphere
      REAL (kind=rk) :: rns (3) ! neutrinosphere radius
      REAL (kind=rk) :: nden (3) ! power-law index

      REAL(KIND=rk), PARAMETER :: wc_mnc=0.5_rk*(wc_mn+wc_mp)
      REAL(KIND=rk), PARAMETER :: fac=3._rk*pc_pi**2*wc_hc**3
      REAL(kind=rk), PARAMETER :: pcon=0.5_rk*wc_hc*wc_hc/wc_mb
      REAL(kind=rk), PARAMETER :: pi32=3.0_rk*pc_pi*pc_pi
      REAL(KIND=rk), PARAMETER :: bnorm=pc_cl/(2.0_rk*pc_pi*wc_hc)**3*
     &     4.0_rk*pc_pi
      REAL(KIND=rk), PARAMETER :: delta=3.0_rk

      REAL(KIND=rk) :: eargn, cel_m, cnu_m, con_en, con_ean, fn, fe
      REAL(KIND=rk) :: earge, scre, epn1, enp1, escr, etnucp, dp0,
     &     xi, etnucn, detnucnp, dn0, fnh, fnp, scr2, scr3
      REAL(KIND=rk), DIMENSION(1:m) :: etapp, etann, etanp, etapn, dup,
     &     den, nby, tem, tem1, tmev, cuq, ceq, cnq, cpq, pre, gamc, s,
     &     xhrepq, ein, qprim, dnuc, s_a
      REAL (KIND=rk), DIMENSION(1:m,2) :: zah
      REAL (KIND=rk), DIMENSION(1:m, 1:config%qn) :: xiq
      REAL (KIND=rk), DIMENSION(1:m,1:nemid) :: cbrems
#ifdef WEAKMAG
      REAL (KIND=rk), PARAMETER :: f2c = 3.706_rk
      REAL (KIND=rk), PARAMETER :: f2np = 1.019_rk
      REAL (KIND=rk), PARAMETER :: f2nn = -0.963_rk
      REAL (KIND=rk), PARAMETER :: cvpf2ca = (1.0_rk+f2c) * SQRT(wc_ga2)
      REAL (KIND=rk), PARAMETER :: cvp = 1.0_rk - wc_cv
      REAL (KIND=rk), PARAMETER :: cap = 0.5_rk * (SQRT(wc_ga2)-ga_s)
      REAL (KIND=rk), PARAMETER :: cvn = -0.5_rk
      REAL (KIND=rk), PARAMETER :: can = 0.5_rk * (-SQRT(wc_ga2)-ga_s)

      REAL (KIND=rk) :: ee, ee2, tep1, ftep1, ln2ep1
      REAL (KIND=rk) :: wkc (2), wknp (3), wknn (3)
#endif /* WEAKMAG */

      REAL(KIND=rk), PARAMETER :: g_nu (1:3) = (/1,1,4/)
      REAL(KIND=rk), PARAMETER :: siglep (1:3) = (/1,-1,0/)

      REAL (KIND=RK) :: eos_self(2), eos_children(2)
      LOGICAL :: ler 

      INTEGER (KIND=ik) :: i,im1,ip1,ie,ie1,is,i_n,j,k,jk,mbrems

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO
!$OMP& PRIVATE(i,j,k,jk)
#endif
      DO jk = 1, n_loc * o_loc

         k = INT ((jk + n_loc - 1) / n_loc)
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

#ifndef CFC_TRANSPORT2
c     Newtonian hydro: decide whether or not to include redshift terms.
c     If so, we temporarily set alpha to the lapse function as computed
c     in sources.F.
      IF (config%i_grtr .NE. 0) THEN
         DO i = 1, m
            alpha (i, j, k) = ephtot (i)
         END DO
         alpha (0  , j, k) = alpha (1  , j, k)
         alpha (m+1, j, k) = 1.0_rk - (1.0_rk - alpha (m, j, k)) *
     &        (r (m) / r (m + 1))
         DO i = 0, m
            scre = r_if (i) - r (i)
            alpha_r (i, j, k) = alpha (i, j, k) +
     &           scre * delta_r_if_inv (i) * 
     &           (alpha (i + 1, j, k) - alpha (i, j, k)) 
         END DO
      ENDIF
#endif /* CFC_TRANSPORT2 */         
         
      DO i = 0, m
         dare (i) = alpha_r (i,j,k) * phi_r (i,j,k) ** 4 * 
     &        r_if (i) ** 2 / alpha_r (i,j,k) ** 3
c        dare (i) = 1.0_rk
      END DO
      DO i = 1, m
         dvol (i) = alpha (i,j,k) * phi (i,j,k) ** 6 * 
     &        (r_if (i) ** 3 - r_if (i-1) ** 3) / 3.0_rk /
     &        alpha (i,j,k) ** 3
c         dvol (i) = phi (i,j,k) ** 2 * delta_r (i)
      END DO

      eran (0) = 0.0_rk
      DO ie = 1, nemid
         eran (ie) = eran_1 * (eran_max / eran_1) ** 
     &        (REAL (ie-1, KIND=rk)  / REAL (nemid -1, KIND=rk))
c         emid0 (ie) = 0.5_rk * (eran (ie-1) + eran (ie))
         emid0 (ie) = (eran (ie) ** 4 - eran (ie-1) ** 4) /
     &        (eran (ie) ** 3 - eran (ie-1) ** 3) * 3.0_rk / 4.0_rk
      END DO

      qentot(:,j,k) = 0.0_rk 
      qyetot(:,j,k,1) = 0.0_rk
      qmotot(:,j,k) = 0.0_rk
      dnutot(:,j,k,:) = 0.0_rk
      enutot(:,j,k,:) = 0.0_rk
      fnutot(:,j,k,:) = 0.0_rk

c      lumtot (:,:) = 0.0_rk

c     precompute auxiliary variables for opacities
      DO i = 1, m

         den (i) = rho (i,j,k) * pc_geog
         nby (i) = den (i) / pc_mb
         tem (i) = t (i,j,k)
         tmev (i) = tem (i) * pc_kmev
         tem1 (i) = 1.0_rk / tmev (i)
                     
         etanp(i) = nby (i) * xnnu (i,j,k,1)
         etapn(i) = nby (i) * xnnu (i,j,k,2)

#ifdef BLOCKING
         scr2 = fac / SQRT (2.0_rk * wc_mnc * tmev(i)) ** 3 *
     &        MAX (etanp(i), 1.0e-16_rk)
         scr3 = fac / SQRT (2.0_rk * wc_mnc * tmev(i)) ** 3 *
     &        MAX (etapn(i), 1.0e-16_rk)
         etnucn = MIN(etapike(scr2), bigeexp)
         etnucp = MIN(etapike(scr3), bigeexp)
         detnucnp=etnucn-etnucp

         dup  (i) = (cpot(i,j,k,3) - cpot(i,j,k,4)) -
     &        tmev (i) * detnucnp

         if (abs(detnucnp) .gt. 1.e-7_rk) then
            escr=exp(detnucnp)
            scr2=escr       -1.0_rk
            scr3=1.0_rk/escr-1.0_rk
         else
            scr2= detnucnp+0.5_rk*detnucnp**2
            scr3=-detnucnp+0.5_rk*detnucnp**2
         endif

         scr2 = sign( MAX(abs(scr2),1.0e-15_rk) , scr2 )
         scr3 = sign( MAX(abs(scr3),1.0e-15_rk) , scr3 )

         scr2 = (etanp (i) - etapn (i)) / scr2
         scr3 = (etapn (i) - etanp (i)) / scr3

         etapn(i) = MAX (MIN (etapn(i), scr2), 0.0_rk)
         etanp(i) = MAX (MIN (etanp(i), scr3), 0.0_rk)
#else
         dup (i) = 0.0_rk
#endif /* BLOCKING */

         etapp(i) = nby (i) * xnnu (i,j,k,2)
         etann(i) = nby (i) * xnnu (i,j,k,1)

#ifdef BLOCKING
         scr2 = MAX (etapp (i),1.0e-15_rk) ! 1e-15 protons/ccm is a small number
         scr2 = 1.5_rk * tmev(i) /
     &        (pcon * (scr2 * pi32) ** (2.0_rk / 3.0_rk))
         etapp(i) = etapp (i) * scr2 / SQRT (1.0_rk+scr2**2) 

         scr3 = MAX (etann (i),1.0e-15_rk) ! 1e-15 neutrons/ccm is a small number
         scr3 = 1.5_rk * tmev(i) /
     &        (pcon * (scr3 * pi32) ** (2.0_rk / 3.0_rk))
         etann(i) = etann (i) * scr3 / SQRT (1.0_rk+scr3**2) 
#endif     
  
      END DO
      DO i_n = 1, config%qn
            xiq (1:m,i_n) = xnnu (1:m,j,k,i_n)
      END DO
      
c     Determine where we need to calculate Bremsstrahlung
      mbrems = 1
      DO i = 2, m
         IF (den (i) .GT. 1e9_rk) THEN
            mbrems = i
         END IF
      END DO

      ceq(1:m)=cpot(1:m,j,k,2)
      CALL eos (den,tem,xiq,xhrepq,zah,ein,pre,gamc,s,cuq,ceq,cnq,cpq,
     &     eos_self,eos_children,mode=1,nsemode=1,ler=ler)

c      cuq(:) = MAX (ceq (:) + cpq (:) -cnq (:) - wc_mq, 0.0_rk)
      cuq(:) = ceq (:) + cpq (:) -cnq (:) - wc_mq
      qprim(:) = (cnq (:) - cpq(:) + wc_mq) + delta

      WHERE (zah(:,2) .NE. 0.0_rk)
         dnuc(:)=den(:)*xhrepq(:)/(pc_mb*zah(:,2))
      ELSE WHERE
         dnuc(:)=0.0_rk
      END WHERE

      DO i = 1, m
#ifdef HOROWITZ
         scr2 = nby (i) * 1e-39
         s_a (i) = 1.0_rk / (1.0_rk +
     &        920.0_rk * scr2 * !A
     &        (1.0 - xiq (i, config%qn) + xiq (i, config%qn) ** 2) /
     &        tmev (i) ** 1.22_rk *
     &        (1.0_rk + 3.05_rk / SQRT(SQRT(tmev (i))) ** 3 / !B
     &        EXP (MIN(6140.0_rk * scr2 * xiq (i, config%qn) * !C
     &        (1.0_rk - xiq (i, config%qn)) / SQRT (tmev (i)) +
     &        1.5e13_rk * scr2 ** 4 / tmev (i) ** 6,bigeexp))))
#else
         s_a (i) = 1.0_rk
#endif /* HOROWITZ */
      END DO

      DO ie = 1, nemid

         cbrems (:,ie) = 0.0d0
         
         DO ie1 = 1, nemid
            DO i = 1, mbrems

               scr2 = alpha (i,j,k) * (emid0 (ie) + emid0 (ie1)) /
     &              tmev (i)
               scr3 = 1.0d0 / 
     &           (1.0_rk + EXP (MIN (emid0(ie1) *
     &              alpha (i,j,k) / tmev (i), bigeexp)))

               cbrems (i,ie) = cbrems (i,ie) +
     &              (eran (ie1) ** 3 - eran (ie1-1) ** 3) / 3.0d0 /
     &              (emid0 (ie) + emid0 (ie1)) * alpha (i,j,k) ** 2 *
     &              SQRT (0.5d0 * pc_pi / scr2) * 
     &              ((1.0d0 - scr3)  * EXP (-scr2) + scr3)
               
            END DO
         END DO

         DO i = 1, mbrems
            cbrems (i,ie) = cbrems (i,ie) * 0.25d0 *
     &           0.5d0  * wc_s0 * wc_ga2 / SQRT(pc_pi) ** 5 *
     &           (1.0d0 / 135.0d0) ** 4 * (wc_mn / wc_me) ** 2 /
     &           SQRT (wc_mn * tmev (i)) * wc_hc ** 3 *  
     &           (etapp (i) ** 2 + etann (i) ** 2 + 28.0d0 / 3.0d0 *
     &        etapp (i) * etann (i))
         END DO

      END DO
 

      DO ie = 1, nemid

         f_nu (:,:,:) = 0.0_rk
         f_nu (:,:,:) = 0.0_rk

         kapa(:,:) = 0.0_rk
         kaps(:,:) = 0.0_rk
         f_h (0,:) = 0.0_rk
         f_h (m,:) = 1.0_rk
         f_h_in (0,:) = 0.0_rk
         f_h_in (m,:) = 0.0_rk

         DO i = 1, m

c     set up energy grid & calculate equilibrium distribution function
            emid (i) = emid0 (ie)  / alpha (i,j,k)
            epl  (i)  = eran (ie  ) / alpha_r (i,j,k)
            emn  (i)  = eran (ie-1) / alpha_r (i,j,k)
            f_eq (i,1) = bnorm /
c     &           (1.0_rk + EXP(MIN((emid(i) - cpot(i,j,k,1)) * tem1 (i),
     &           (1.0_rk + EXP(MIN((emid(i) - cuq(i)) * tem1 (i),
     &           bigeexp)))
            f_eq (i,2) = bnorm /
c     &           (1.0_rk + EXP(MIN((emid(i) + cpot(i,j,k,1)) * tem1 (i),
     &           (1.0_rk + EXP(MIN((emid(i) + cuq (i)) * tem1 (i),
     &           bigeexp)))
            f_eq (i,3) = bnorm /
     &           (1.0_rk + EXP(MIN(emid(i) * tem1 (i), bigeexp)))

c     -----------------------------------------------------------------
c     opacities
c     -----------------------------------------------------------------
#ifdef WEAKMAG
c     weak magnetism correction
            ee = emid (i) / wc_mb ! assumes nucleon mass in vacuum
            tep1 = 2.0_rk*ee+1.0_rk
            ln2ep1 = log(tep1)
            ftep1 = 1.333333333_rk*ee+1.0_rk
            ee2 = ee**2
            scr2=(1.0_rk+4.0_rk*ee+5.333333333_rk*ee2)+
     &           3.0_rk*wc_ga2*ftep1**2+
     &           2.666666666_rk*f2c*ee2+
     &           1.666666666_rk*ee2*(1.0_rk+0.4_rk*ee)*f2c**2
            scr3=4.0_rk*cvpf2ca*ee*ftep1
            scre = ((1.0_rk+3.0_rk*wc_ga2)*tep1**3)
            wkc (1) = (scr2+scr3)/scre
            wkc (2) = (scr2-scr3)/scre
            
            scr2 = cvp**2*
     &           ((ee-1.0_rk)/(2.0_rk*ee2*ee)*ln2ep1+
     &           (3.0_rk+12.0_rk*ee+9.0_rk*ee2-
     &           10.0_rk*ee2*ee)/(3.0_rk*ee2*tep1**3))+
     &           cap**2*
     &           ((1.0_rk+ee)/(2.0_rk*ee2*ee)*ln2ep1-
     &           (10.0_rk*ee2*ee+27.0_rk*ee2+
     &           18.0_rk*ee+3.0_rk)/(3.0_rk*ee2*tep1**3))+
     &           cvp*f2np*(ln2ep1/ee2-0.666666666_rk*
     &           (3.0_rk+15.0_rk*ee+22.0_rk*ee2)/
     &           (ee*tep1**3))+
     &           f2np**2*
     &           (0.25_rk/ee2*ln2ep1+
     &           (8.0_rk*ee*ee2-22.0_rk*ee2-15.0_rk*ee-3.0_rk)/
     &           (6.0_rk*ee*tep1**3))
            scr3 = (cvp+f2np)*cap*
     &           (ln2ep1/ee2-(2.0_rk+10.0_rk*ee+9.33333333_rk*ee2)/
     &           (ee*tep1**3))
            scre = 1.5_rk/(cvp**2+5.0_rk*cap**2)
            wknp (1) = (scr2+scr3)*scre
            wknp (2) = (scr2-scr3)*scre
            wknp (3) = 2.0_rk*wknp(1)*wknp(2)/(wknp(1)+wknp(2))
            
            scr2 = cvn**2*
     &           ((ee-1.0_rk)/(2.0_rk*ee2*ee)*ln2ep1+
     &           (3.0_rk+12.0_rk*ee+9.0_rk*ee2-
     &           10.0_rk*ee2*ee)/(3.0_rk*ee2*tep1**3))+
     &           can**2*
     &           ((1.0_rk+ee)/(2.0_rk*ee2*ee)*ln2ep1-
     &           (10.0_rk*ee2*ee+27.0_rk*ee2+
     &           18.0_rk*ee+3.0_rk)/(3.0_rk*ee2*tep1**3))+
     &        cvn*f2nn*(ln2ep1/ee2-0.666666666_rk*
     &           (3.0_rk+15.0_rk*ee+22.0_rk*ee2)/
     &           (ee*tep1**3))+
     &           f2nn**2*
     &           (0.25_rk/ee2*ln2ep1+
     &           (8.0_rk*ee*ee2-22.0_rk*ee2-15.0_rk*ee-3.0_rk)/
     &           (6.0_rk*ee*tep1**3))
            scr3 = (cvn+f2nn)*can*
     &           (ln2ep1/ee2-(2.0_rk+10.0_rk*ee+9.33333333_rk*ee2)/
     &           (ee*tep1**3))
            scre = 1.5_rk/(cvn**2+5.0_rk*can**2)
            wknn (1) = (scr2+scr3)*scre
            wknn (2) = (scr2-scr3)*scre
            wknn (3) = 2.0_rk*wknn(1)*wknn(2)/(wknn(1)+wknn(2))
#endif /* WEAKMAG */


c     absorption on neutrons
            scre = 0.25_rk * wc_s0 *
     &           MAX(emid(i) + wc_mq + dup (i), 0.0_rk) / wc_me ** 2 * 
     &           SQRT(MAX(MAX(emid(i) + wc_mq + dup(i), 0.0_rk) ** 2 -
     &           wc_me ** 2, 0.0_rk))
c            scre = 0.25_rk * wc_s0 * (emid(i) + wc_mq) / wc_me ** 2 * 
c     &           SQRT((emid(i) + wc_mq) ** 2 - wc_me ** 2)
            con_en = (1._rk+3._rk*wc_ga2) * scre

c            etanp(i) = nby (i) * xnnu (i,j,k,1)

c            cnu_m = cpot(i,j,k,2)+cpot(i,j,k,4)-cpot(i,j,k,3) - wc_mq  
c            cel_m = cpot(i,j,k,2)          
            cnu_m = cuq (i)
            cel_m = ceq (i)

            eargn = (emid(i)-cnu_m      )*tem1(i) 
            earge = (emid(i)-cel_m+wc_mq+dup(i))*tem1(i) 
c            earge = (emid(i)-cel_m+wc_mq)*tem1(i) 
            
            IF (earge .GT. 35._rk) THEN
               fe=1._rk
            ELSE
               fe=EXP(earge)/( 1._rk + EXP(earge) )
            END IF
                  
            fn=1._rk + EXP(MIN(-eargn,bigeexp))
            kapa (i,1) = con_en*fe*fn*etanp(i)
#ifdef WEAKMAG
     &           *wkc(1)
#endif /* WEAKMAG */

c     absorption on protons
            scre = 0.25_rk * wc_s0 * 
     &           MAX(emid(i) - wc_mq - dup (i), 0.0_rk) / wc_me ** 2 *
     &           SQRT(MAX(MAX(emid(i) - wc_mq - dup (i), 0.0_rk) ** 2 -
     &           wc_me ** 2, 0.0_rk))
c            scre = 0.25_rk * wc_s0 * 
c     &           MAX(emid(i) - wc_mq, 0.0_rk) / wc_me ** 2 
c     &           * SQRT (MAX ((emid(i) - wc_mq) ** 2 - wc_me ** 2,
c     &           0.0_rk))

c            scre = 0.25_rk*wc_s0* (emid(i)-wc_mq)/wc_me**2 *
c     &           SQRT( (emid(i)-wc_mq)**2-wc_me**2 )
            con_ean = (1._rk+3._rk*wc_ga2)*scre
            
c            etapn(i) = nby (i) * xnnu (i,j,k,2)
            
C            cnu_m = -(cpot(i,j,k,2)+cpot(i,j,k,4)-cpot(i,j,k,3) - wc_mq)
C            cel_m = -cpot(i,j,k,2) 
            cnu_m = -cuq(i)
            cel_m = -ceq(i)
            
            eargn = (emid(i)-cnu_m      )*tem1(i) 
            earge = (emid(i)-cel_m-wc_mq-dup(i))*tem1(i) 
c            earge = (emid(i)-cel_m-wc_mq)*tem1(i) 
            
            IF (earge .GT. 35._rk) then
               fe=1._rk
            ELSE
               fe=EXP(earge)/( 1._rk + EXP(earge) )
            END IF

            fn=1._rk + EXP(MIN(-eargn,bigeexp))
            kapa(i,2) = con_ean*fe*fn*etapn(i)
#ifdef WEAKMAG
     &           *wkc(2)
#endif /* WEAKMAG */
           
c     scattering on nucleons
            scre = 0.25_rk*wc_s0*(emid(i)/wc_me)**2
            
c     0th Legendre moment of scattering kernel
#ifdef WEAKMAG
#ifdef STRANGE
            kaps (i,:) = scre * (                     
     &           etapp(i) * ((1.0_rk - wc_cv)**2 + 0.75_rk * 
     &           (SQRT(wc_ga2)-ga_s) ** 2 * s_a (i)) *
     &           wknp(:) +
     &           etann(i) * (0.25_rk * (1._rk + 3._rk * 
     &           (-SQRT(wc_ga2)-ga_s) ** 2) * s_a (i))) *
     &           wknn(:)
#else /* STRANGE */
            kaps (i,:) = scre * (                     
     &           etapp(i) * ((1.0_rk - wc_cv)**2 + 0.75_rk * wc_ga2 *
     &           s_a (i)) * wknp (:) +
     &           etann(i) * (0.25_rk * (1._rk + 3._rk * wc_ga2 *
     &           s_a (i)) * wknn(:)))
#endif /* STRANGE */
#else /* WEAKMAG */
c no weak magnetism
#ifdef STRANGE
            kaps (i,:) = scre * (                     
     &           etapp(i) * ((1.0_rk - wc_cv)**2 + 0.75_rk * 
     &           (SQRT(wc_ga2)-ga_s) ** 2 * s_a (i)) +
     &           etann(i) * (0.25_rk * (1._rk + 3._rk * 
     &           (-SQRT(wc_ga2)-ga_s) ** 2) * s_a (i)))
#else /* STRANGE */
            kaps (i,:) = scre * (                     
     &           etapp(i) * ((1.0_rk - wc_cv)**2 + 0.75_rk * wc_ga2 *
     &           s_a (i)) +
     &           etann(i) * (0.25_rk * (1._rk + 3._rk * wc_ga2 *
     &           s_a (i))))
#endif /* STRANGE */
#endif /* WEAKMAG */
            
c     0th and 1st Legendre moment of scattering kernel
c            kaps (i,:) = scre * (                     
c     &           etapp(i) * ((1.0_rk - wc_cv)**2 + 0.75_rk * wc_ga2 -
c     &           ((1.0_rk-wc_cv)**2 - 0.25_rk * wc_ga2) / 3.0_rk) +
c     &           etann(i) * (0.25_rk * (1._rk + 3._rk * wc_ga2) -
c     &           0.25_rk*(1.0_rk-wc_ga2) / 3.0_rk))

         END DO

#if 0
c     scattering on nuclei
         DO i_n = 3, config%qn - 1

c            fnp = MAX (pc_nuc (i_n,1)-20.0_rk, 0.0_rk)
c            fnp = MIN (fnp , 8.0_rk)
c            fnh = MAX (40.0_rk - 
c     &           (pc_nuc (i_n,2) - pc_nuc (i_n,1)) , 0.0_rk)
c            fnh = MIN (fnh , 6.0_rk)
c            scr3 = fnp * fnh

            DO i = 1, m

               scre = 0.125_rk*wc_s0*(emid(i)/wc_me)**2
               scr2 = nby (i) * xnnu (i,j,k,i_n) / pc_nuc (i_n,2)
               scre = scre * 0.25_rk * pc_nuc (i_n,2) ** 2 *
     &              (wc_ca - wc_cv + (2.0_rk - wc_ca - wc_cv) *
     &              (2.0_rk * pc_nuc (i_n,1) - pc_nuc (i_n,2)) / 
     &              pc_nuc (i_n,2)) ** 2
               IF (zah (i,2) .EQ. 0.0_rk) 
     &              kaps(i,:) = kaps(i,:) + scr2 * scre 

            END DO
         END DO

         DO i = 1, m

            scre = 0.125_rk*wc_s0*(emid(i)/wc_me)**2
            scr2 = dnuc (i)
            IF (zah (i,2) .NE. 0.0_rk) 
     &           kaps(i,:) = kaps(i,:) + scr2 *
     &           scre * 0.25_rk * zah (i,2) ** 2 *
     &           (wc_ca - wc_cv + (2.0_rk - wc_ca - wc_cv) *
     &           (2.0_rk * zah (i,1) - zah (i,2)) / 
     &           zah(i,2)) ** 2

c     absorption on nuclei
            fnp = MAX (zah (i,1)-20.0_rk, 0.0_rk)
            fnp = MIN (fnp , 8.0_rk)
            fnh = MAX (40.0_rk - 
     &           (zah (i,2) - zah (i,1)) , 0.0_rk)
            fnh = MIN (fnh , 6.0_rk)
            scr3 = fnp * fnh

            cnu_m = ceq(i)+cpq(i)-cnq(i) - wc_mq  
            cel_m = ceq(i)
            
            eargn = tem1(i)*(emid(i)      -cnu_m)
            earge = tem1(i)*(emid(i)+qprim(i)-cel_m)
            
            IF (eargn .GT. 35.0_rk .AND. earge .GT. 35.0_rk) THEN
               scr2 = EXP (MIN (eargn-earge,bigeexp))
            ELSE
               scr2 = (1.0_rk + EXP (MIN (eargn,bigeexp))) / 
     &              (1.0_rk + EXP (MIN (earge,bigeexp)))
            END IF
            
            IF (emid(i)+qprim(i) .LE. wc_me) THEN
               scr2 = 0.0_rk
            ELSE
               scr2 = scr2 * wc_s0 * wc_ga2 / 14.0_rk * 
     &              (emid (i) + qprim (i)) / wc_me ** 2 *
     &              SQRT((emid (i) + qprim (i)) ** 2 - wc_me ** 2)
            END IF

            kapa (i,1) = kapa(i,1) + scr2 * scr3 * dnuc (i)
            
         END DO
#endif

c     Bremsstrahlung (secundum Burrows, Reddy & Thompson)
         DO i = 1, m
            kapa (i,3) = cbrems (i,ie) 
         END DO
c     energy exchange in scattering reactions
         DO i = 1, m
            kapa (i,3) = kapa (i,3) +
#undef MORE_HEATING
#ifdef MORE_HEATING
     &           3.0_rk*
#endif
     &           (ABS(emid (i) - 3.0_rk * tmev(i))/ wc_mnc) * kaps (i,3)
         END DO

         DO i = 1, m
            IF (tem(i).lt.1e8_rk) THEN
               kapa (i,:) = 0.0_rk
               kaps (i,:) = 0.0_rk
            END IF
         END DO
         
c         IF (ANY(kaps.NE.kaps .OR. kapa.ne.kapa)) THEN
c     PRINT *,'Problem with mu_e (cpo)',j,k
c            DO i = 1, m
c               PRINT*,'fmt',i,j,k,den(i),tem(i),xiq(i,config%qn),
c     &              cpot(i,j,k,2),kaps(i,:),kapa(i,:)
c            END DO
c         END IF

         tau (m,:) = 0.0_rk
         DO i = m-1, 1, -1
            tau (i,:) = tau (i+1,:) +
     &           (kapa (i,:) + kaps (i,:)) * 
c     &           SQRT ((kapa (i,:) + kaps (i,:)) * kapa (i,:)) * 
     &           delta_r_if (i) * phi_r (i,j,k) ** 2 
         END DO

c     -----------------------------------------------------------------
c     Formal solution of the transport equation for two rays
c     (outgoing/ingoing) to obtain an estimate for the flux factor
c     -----------------------------------------------------------------

c     ingoing ray
         f_nu (m,1,:) = 0.0_rk
         f_nu (0,0,:) = 0.0_rk
         f_eq (0,:) = f_eq (1,:)
         f_eq (m+1,:) = f_eq (m,:)
c     For computing the flux factor, it could make sense to
c     use the total opacity (absorption + scattering)
         kapeff (1:m,:) = kaps (1:m,:) + kapa (1:m,:)
         DO is = 1, 3
            DO i = 1, m 
c     "anti-scattering"
            kapeff (i,is) = MAX (kapeff (i,is) -
     &              0.4_rk / (r (i) * phi (i,j,k) ** 2), 0.0_rk)
c            kapeff (i,is) = MAX (kapeff (i,is) -
c     &              1.0_rk / (r (i) * phi (i,j,k) ** 2), 0.0_rk)
            END DO
         END DO
         kapeff (0,:) = kapeff (1,:)
         kapeff (m+1,:) = kapeff (m,:)


cC$$$     formal solution along a single radial ray
C$$$c     ingoing ray
C$$$         DO i = m-1, 0, -1
C$$$            f_nu (i,1,:) = 
C$$$     &           (f_nu (i+1,1,:) +
C$$$     &           delta_r (i+1) * phi (i+1,j,k) ** 2 *
C$$$     &           0.5_rk * (1*kapeff (i,:) + 1*kapeff (i+1,:)) * 
C$$$     &           0.5_rk * (1*f_eq (i,:) + 1*f_eq (i+1,:))) /
C$$$     &           (1.0_rk + delta_r (i+1) * phi (i+1,j,k) ** 2 *
C$$$     &           0.5_rk * (1*kapeff (i,:) + 1*kapeff (i+1,:)))
C$$$         END DO

C$$$c     outgoing ray
C$$$c         f_nu (0,0,:) = f_eq (1,:)
C$$$         f_nu (0,0,:) = f_nu (0,1,:)
C$$$         DO i = 1, m
C$$$            f_nu (i,0,:) = 
C$$$     &           (f_nu (i-1,0,:) +
C$$$     &           delta_r (i) * phi (i,j,k) ** 2 *
C$$$     &           0.5_rk * (1*kapeff (i,:) + 1*kapeff (i+1,:)) * 
C$$$     &           0.5_rk * (1*f_eq (i,:) + 1*f_eq (i+1,:))) /
C$$$     &           (1.0_rk + delta_r (i) * phi (i,j,k) ** 2 *
C$$$     &           0.5_rk * (1*kapeff (i,:) + 1*kapeff (i+1,:)))
C$$$         END DO


c     Solution including scattering using Riccati transformation
         DO is = 1, 3
            chi (m) = 1.0_rk
            psi (m) = 0.0_rk
            DO i = m - 1, 0, -1
               scre = 0.5_rk * (kapa (i,is) + kapa (i+1,is)) *
     &              delta_r_if (i) * phi_r (i,j,k) ** 2
               scr2 = chi (i+1) + scre
               scr3 = (kapeff (i,is) + kapeff (i+1,is)) *
     &              delta_r_if (i) * phi_r (i,j,k) ** 2
               chi (i) = 2.0_rk * scr2 /
     &              (1.0_rk + SQRT (1.0_rk + 2.0_rk * scr3 * scr2))
               psi (i) = (psi (i+1) - scre * 
     &              0.5_rk * (f_eq (i,is) + f_eq (i+1,is))) /
     &              (1.0_rk + 0.5_rk * scr3 * chi (i))
            END DO
            f_nu (0,0,is) = - psi (0) / chi (0)
            DO i = 1, m
               scr3 = 0.5_rk  * (kapeff (i,is) + kapeff (i+1,is)) *
     &              delta_r_if (i-1) * phi_r (i-1,j,k) ** 2
               f_nu (i,0,is) = (f_nu (i-1,0,is) - scr3 * psi (i)) /
     &              (1.0_rk + scr3 * chi (i))
               f_nu (i,1,is) = psi (i) + chi (i) * f_nu (i,0,is)
            END DO
            DO i = 1, m
               scr2 = f_nu (i,0,is) + f_nu (i,1,is)
               scr3 = f_nu (i,0,is) - f_nu (i,1,is)
               f_nu (i,0,is) = scr2
               f_nu (i,1,is) = scr3
            END DO
         END DO

c     flux factor
         DO is = 1, 3 
            f_h (0,is) = 0.0_rk 
c            rns (is) = 0.0_rk
c            tns (is) = tem (i)
c            nden (is) = 1.0_rk

            DO i = 1, m
               scre = f_nu (i,1,is) / (f_nu (i,0,is) + 1e-80_rk)
               scre = MAX (scre, 1e-80_rk)
               IF (ABS (scre - 1.0_rk) .GT. 0.005_rk) THEN
                  f_h (i,is) = 1.0_rk + 2.0_rk * scre / 
     &                 (1.0 - scre) + 2.0_rk / LOG (scre)
               ELSE
                  scre = scre - 1.0_rk
                  f_h (i,is) = 
     &                 - (((((863.0_rk * scre - 1134.0_rk)) * 
     &                 scre + 1596.0_rk) * scre - 2520.0_rk) *
     &                 scre + 5040.0_rk) * scre / 30240.0_rk
               END IF
            END DO
            END DO

c     -----------------------------------------------------------------
c     Solve the stationary 0th moment equation now that the flux
c     factor has been determined:
c
c     The ODE r^-2 d(H r^2)/dr = kappa_star * (J_eq - H/f_h)
c     has critical points where f_h =0. We therefore solve it
c     on different patches by right/left integration using 
c     H(r_crit)=0 as boundary condition.
c     -----------------------------------------------------------------

c     new: approximately account for Doppler shift by means of
c     an effective lab-frame opacity for electron neutrinos and
c     antineutrinos
         DO i = 1, m
            scr2 = MAX(MIN (v_1 (i,j,k) / phi (i,j,k) ** 2, 0.5_rk),
     &           -0.5_rk)
            kapa (i,1) = kapa (i,1) * (w (i,j,k) * (1.0_rk - 
     &           0.5_rk * (f_h (i-1,1)  + f_h (i,1)) * scr2)) ** 3
            kapa (i,2) = kapa (i,2) * (w (i,j,k) * (1.0_rk - 
     &           0.5_rk * (f_h (i-1,2)  + f_h (i,2)) * scr2)) ** 3
         END DO

         f_nu (0,0,:) = f_eq (1,:)
         f_eq (0,:) = f_eq (1,:)
         f_eq (m+1,:) = f_eq (m,:)
         kapa (0,:) = kapa (1,:)
         kapa (m+1,:) = kapa (m,:)
         f_nu (0,0,:) = f_nu (0,1,:)
         
         DO is = 1, 3

            f_nu (1:m,:,is) = 0.0_rk

c     positive flux factor
            DO i = 1, m
               IF (f_h (i-1,is) .GT. 0.51_rk) THEN
c     two-moment closure
                  scr2 = (1.0 - 2.0_rk * f_h (i-1,is) +
     &                 4.0_rk * f_h (i-1,is) ** 2) / 3.0_rk
                  scr3 =  (1.0_rk - 4.0_rk * f_h (i-1,is) ** 2) /
     &                 3.0_rk
                  f_h (i,is) = MIN(f_h (i,is),
     &                 f_h (i-1,is) + delta_r (i) *
     &                 (f_h (i-1,is) * (scr2 - 1.0_rk) * r_inv (i) +
     &                 kapa (i,is) * (f_h (i-1,is) ** 2 -
     &                 scr2 + f_eq (i,is) / f_nu (i-1,0,is) *
     &                 f_h (i-1,is) * scr2) +
     &                 kaps (i,is) * f_h (i-1,is) ** 2) / scr3)
                  f_h (i,is) = MIN (f_h (i,is), 1.0_rk)
                  f_h (i,is) = MAX (f_h (i,is), 0.51_rk)
               END IF
               IF (f_h (i,is) .GT. 0.0_rk) THEN
                  f_nu (i,0,is) = 
     &                 (f_nu (i-1,0,is) * dare (i-1) + dvol (i) *
     &                 0.5_rk * (kapa (i,is) + kapa (i+1,is)) * 
     &                 0.5_rk * (f_eq (i,is) + f_eq (i+1,is))) /
     &                 (dare (i) + dvol (i) * 0.5_rk *
     &                 (kapa (i,is) + kapa (i+1,is)) / f_h (i,is)) 
               END IF
            END DO

c     negative flux factor
            DO i = m-1, 1, -1
               IF (f_h (i,is) .LT. 0.0_rk) THEN
                  f_nu (i,1,is) = 
     &                 (f_nu (i+1,1,is) * dare (i+1) - dvol (i+1) * 
     &                 0.5_rk * (kapa (i,is) + kapa (i+1,is)) * 
     &                 0.5_rk * (f_eq (i,is) + f_eq (i+1,is))) /
     &                 (dare (i) - dvol (i+1) * 0.5_rk *
     &                 (kapa (i,is) + kapa (i+1,is)) / f_h (i,is)) 
               END IF
               f_nu (i,0,is) = f_nu (i,0,is) + f_nu (i,1,is)
            END DO
         END DO

c     calculate contribution to source term
         DO i = 0, m 
            lum (i,:) = f_nu (i,0,:) * dare (i)
c            lum (i,3) = 0.0_rk
         END DO
         DO i = 1, m
            scre = (epl (i) ** 3 - emn (i) ** 3) / 3.0_rk
            scr1 = scre * f_nu (i,0,:) 
            fnutot (i,j,k,:) = fnutot (i,j,k,:) + scr1 * emid0 (ie) /
     &           alpha_r (i,j,k)
            scr1 = ABS (scr1) / MAX (ABS (f_h (i,:)), 1e-30_rk)
            enutot (i,j,k,:) = enutot (i,j,k,:) + scr1 * emid0 (ie) /
     &           alpha_r (i,j,k)
            dnutot (i,j,k,:) = dnutot (i,j,k,:) + scr1 

            scre = SUM (lum (i,:) * g_nu (:) -
     &           lum (i-1,:) * g_nu (:)) /
     &           (dvol (i) * alpha (i,j,k) ** 5) * 
c     &          (dvol (i) * alpha (i,j,k) ** 3) * 
     &           (eran (ie) ** 4 - eran (ie-1) ** 4) / 4.0_rk
            qentot (i,j,k)   = qentot (i,j,k) - scre * pc_meverg 
            scre = SUM (lum (i,:) * siglep (:) - 
     &           lum (i-1,:) * siglep (:)) /
     &           (dvol (i) * alpha (i,j,k) ** 4) * 
c     &           dvol (i) * alpha (i,j,k) ** 2 * 
     &           (eran (ie) ** 3 - eran (ie-1) ** 3) / 3.0_rk
            qyetot (i,j,k,1) = qyetot (i,j,k,1) - scre * pc_mb 
         END DO

      END DO

      DO is = 1, 3
         DO i = 1, m
            fnutot (i,j,k,is) = fnutot (i,j,k,is) * pc_meverg *
     &           sqrt_gamma (i,j,k)
            enutot (i,j,k,is) = enutot (i,j,k,is) * pc_meverg / pc_cl *
     &           sqrt_gamma (i,j,k)
            dnutot (i,j,k,is) = dnutot (i,j,k,is) / pc_cl *
     &           sqrt_gamma (i,j,k)
         END DO
      END DO


      IF (bounce) THEN
         DO i = 1, m
            qentot (i,j,k  ) = qentot (i,j,k  ) * w (i,j,k) ** 2
            qyetot (i,j,k,1) = qyetot (i,j,k,1) * w (i,j,k) 
            IF (xnnu(i,j,k,config%qn) .GT. 0.7_rk .OR.
     &           entropy(i,j,k) .LT. 0.2_rk) THEN
               qyetot(i,j,k,1)=0.0_rk
               qentot(i,j,k) = 0.0_rk
            END IF
            IF (xnnu(i,j,k,config%qn) .GT. 0.57_rk) THEN
               qyetot(i,j,k,1)=min(qyetot(i,j,k,1),0.0_rk)
            END IF
         END DO
      ELSE 
         DO i = 1, m
            qentot (i,j,k  ) = 0.0_rk
            qyetot (i,j,k,1) = 0.0_rk
         END DO
      END IF

c      print *,g_nu
c      DO i = 30, 70
c         WRITE (32,*),r(i)/1e5_rk,lum(i,1),
c         PRINT '(2(x,i3),7(x,e16.7))',i,j,r(i)/1e5_rk,
c     &        4.0_rk*pc_pi*fnutot(i,j,k,1)*r_if(i)**2,
c     &        4.0_rk*pc_pi*fnutot(i,j,k,2)*r_if(i)**2,
c     &        4.0_rk*pc_pi*fnutot(i,j,k,3)*r_if(i)**2,
c     &        (qentot(i,j,k)*pc_ergmev)/(rho(i,j,k)*pc_geog/pc_mb),
c     &        qyetot(i,j,k,1)/(rho(i,j,k)*pc_geog),
c     &        sum (kapeff (i:m,1)*(delta_r (i:m)))
c         PRINT '(2(x,i3),7(x,e16.7))',i,j,r(i)/1e5_rk,fnutot(i,j,k,:),
c     &        enutot(i,j,k,:)
c      END DO


#ifndef CFC_TRANSPORT2
c     Newtonian hydro: reset alpha to unity.
      IF (config%i_grtr .NE. 0) THEN
         alpha (:, j, k) = 1.0_rk
         alpha_r (:, j, k) = 1.0_rk
      ENDIF
#endif /* CFC_TRANSPORT2 */         
      
      END DO !j,k-loop
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END DO
#endif


c      STOP
      RETURN

      END SUBROUTINE fast_neutrino

c     ==================================================================




      END MODULE coconut

c     ==================================================================
