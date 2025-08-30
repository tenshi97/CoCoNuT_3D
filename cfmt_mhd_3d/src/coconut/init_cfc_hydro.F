c     ==================================================================

      MODULE gr_hyd_init

c     ==================================================================

      IMPLICIT NONE
      
      CONTAINS



c     ==================================================================

      SUBROUTINE init_cfc_hydro(restart)
      
c     ==================================================================

      USE precision
      
!      USE gfloat_hy, ONLY: cfl
!      USE intgrs_hy, ONLY: isym
      USE totare_hy
      USE phycon
      
      USE size_cfc
      USE parameters_cfc
      USE conserved_cfc
      USE hydro_primitives_cfc
      USE boundary_cfc
      USE grid_cfc
      USE metric_cfc
      USE perm_aux_cfc

      USE poisson
#ifdef CFC_MHD
      USE induction
#endif
      
      USE eos3d_routine

      USE gr_parameters
      USE gr_initial_setup
      
#ifdef EXTRACT_GW
      USE gw3d, ONLY: qijq1, qijq2
#endif
      
      use configure

      IMPLICIT NONE
      
      logical, intent(in) :: restart
      real(kind=rk)       :: eos3ds_self(2), eos3ds_children(2)
      real(kind=rk) :: alpha_grid,dr_tmp,r_tmp
      real (kind=rk) :: xnu_sum (1 : m)
      
      integer(kind=ik) :: i, j, k, jk, l
      
      logical ler
      
      print*
      print*, 'Initialization of module CoCoNuT.'

      call read_parameters
      metric_step=1
      
c     same value for physical constants in VERTEX and CoCoNuT
      c_light=pc_cl
      g_grav=pc_gc
      
c     calculate some constants
      pi = asin (1.0_rk) * 2.0_rk
      p_geom_factor = pc_gc / pc_cl ** 4
      p_geom_factor_inv = 1.0_rk / p_geom_factor
      rho_geom_factor = pc_gc / pc_cl ** 2
      rho_geom_factor_inv = 1.0_rk / rho_geom_factor

      if (config%nsdim .ge. 2) then
         equatorial_symmetry = (config%isym .eq. 1)
         if (equatorial_symmetry) then
            delta_theta = pi / real (2*n,kind=rk)
            equatorial_symmetry_factor = 2.0_rk
            symm_v_theta_lower (:) = symm_v_theta_lower_eq (:)
         else
            delta_theta = pi / real (n,kind=rk)
            equatorial_symmetry_factor = 1.0_rk
            symm_v_theta_lower (:) = symm_v_theta_lower_noneq (:)
         endif
         delta_theta_if = delta_theta

c     compute theta grid vectors
         do j = - 3, n + 4
            theta (j) = (0.5_rk + real(j - 1, kind = rk)) /
     &           real (n,kind=rk) * pi
         enddo
         do j = - 1, n + 1
            theta_if (j) = real (j, kind = rk) / 
     &           real (n, kind = rk) * pi
         enddo
         if (equatorial_symmetry) then
            theta (:) = 0.5_rk * theta (:)
            theta_if (:) = 0.5_rk * theta_if (:)
         end if
      else
         equatorial_symmetry = .false.
         delta_theta = pi
         delta_theta_if = delta_theta
         equatorial_symmetry_factor = 1.0_rk

c     compute theta grid vectors
         do j = - 3, n + 4
            theta (j) = 0.5_rk * pi
         enddo
         do j = - 1, n + 1
            theta_if (j) = 0.5_rk * pi
         enddo
      endif

      delta_varphi = 2.0_rk * pi / real(o, kind=rk)
      delta_varphi_if = delta_varphi
      
      if (config%nsdim .eq. 3) then
         do k = - 3, o + 4
            varphi (k) = (0.5_rk + real(k - 1, kind = rk)) /
     &        real (o, kind=rk) * 2.0_rk * pi
         end do
         do k = - 1, o + 1
            varphi_if (k) = real (k, kind = rk) / 
     &           real (o, kind = rk) * 2.0_rk * pi
         enddo
      else
         varphi (:) = 0.0_rk
         varphi_if (:) = 0.0_rk
      end if

c     copy radial grid used in PROMETHEUS/VERTEX
      r   (1:m)   =xzntot(1:m)
      r_if(1:m)   =xzrtot(1:m)

c     extrapolate outwards:
      r_tmp=r_if(m)
      dr_tmp=xzrtot(m)-xzltot(m)
      alpha_grid=dr_tmp/(xzrtot(m-1)-xzltot(m-1))
      do i=m+1,m+4
         dr_tmp=alpha_grid*dr_tmp
         r_tmp=r_tmp+dr_tmp
         if (i.le.m+3) r_if(i)=r_tmp
         r(i)=0.5_rk*(r_if(i-1)+r_tmp)
      end do

      r (0) = - r (1)
      r (- 1) = - r (2)
      r (- 2) = - r (3)
      r (- 3) = - r (4)

      do i = 0, m + 3
         delta_r_if (i) = r (i + 1) - r (i)
      enddo
      delta_r_if (- 1) = delta_r_if (1)
      delta_r_if (- 2) = delta_r_if (2)
      delta_r_if (- 3) = delta_r_if (3)

      r_if (0) = 0.0_rk
      
      r_if (- 1) = - r_if (1)
      r_if (- 2) = - r_if (2)
      r_if (- 3) = - r_if (3)
      
      do i = 1, m + 3
         delta_r (i) = r_if (i) - r_if (i - 1)
      enddo
      
      delta_r (0) = delta_r (1)
      delta_r (- 1) = delta_r (2)
      delta_r (- 2) = delta_r (3)

      do i = 1, m + 3
         amp (i) = delta_r_if (i) / delta_r_if (i - 1)
      enddo

      call calculate_perm_aux_values

c     copy data from arrays used in PROMETHEUS/VERTEX:
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
!$OMP PARALLEL
!$OMP& PRIVATE(i,j,k,jk,l)

!$OMP DO
!$OMP& SCHEDULE(static)
#endif    
      do jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         do i = 1, m
            rho(i,j,k)=dentot(i,j,k)*rho_geom_factor
            v_1(i,j,k)=vextot(i,j,k)/pc_cl
            v_2(i,j,k)=veytot(i,j,k)/pc_cl*r(i) 
            v_3(i,j,k)=veztot(i,j,k)/pc_cl*r(i)*sinus_theta(j)
            
#ifdef CFC_TRANSPORT2
            w  (i,j,k)=1.0_rk/  
     &           sqrt(1.0_rk-(vextot(i,j,k)**2+veytot(i,j,k)**2+ 
     &           veztot(i,j,k)**2)/pc_cl**2)
#else
            w  (i,j,k)=1.0_rk
#endif /* CFC_TRANSPORT2 */
            t   (i,j,k)=temtot(i,j,k)
         end do

         do l = 1, config%qn
            do i = 1, m
               xnnu(i,j,k,l)=xnutot(i,j,k,l)
            end do
         end do
         
c     Make sure mass fractions sum up to unity
         xnu_sum (:) = 0.0_rk
           
         do l = 1, config%qn-1
            do i = 1, m
               xnu_sum(i) = xnu_sum (i) + xnnu (i,j,k,l)
            end do
         end do

         DO i = 1, m
            xnu_sum(i) = 1.0_rk / xnu_sum(i)
         END DO

         do l = 1, config%qn-1
            do i = 1, m
               xnnu (i, j, k, l) =
     &              xnnu (i, j, k, l) * xnu_sum (i)
            end do
         enddo
         
         
         call eos3ds_CFC(1,j,k,ler, eos3ds_self, eos3ds_children)

         do i = 1, m
#ifdef CFC_TRANSPORT2
            v_squared (i, j, k) = 1.0_rk - 1.0_rk / w (i, j, k) ** 2
#else
            v_squared (i, j, k) =
     &           (vextot(i,j,k)**2+veytot(i,j,k)**2+veztot(i,j,k)**2)/
     &           pc_cl**2
#endif /* CFC_TRANSPORT2 */            
         end do

      end do !j,k-loop

      if (.not. restart) then

         rho_nuc_thr_passed = .false.
         rho_nuc_thr_passed_old = .false.
         first_iteration = 1
         iteration = 0
         t_total = 0.0_rk
      
c     metric quantities:
#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
!$OMP PARALLEL DO
!$OMP& SCHEDULE(static)
!$OMP& PRIVATE(i,j,k,jk,l) 
#endif    
         do jk = 1, n_loc * o_loc

            k = int((jk + n_loc - 1) / n_loc )
            j = (n_s - 1) + (jk - (k - 1) * n_loc)
            k = k + o_s - 1
            
            phi       (:,j,k)=1.0_rk
            alpha     (:,j,k)=1.0_rk
            sqrt_gamma(:,j,k)=1.0_rk
            
         end do
      end if

#ifdef CFC_MHD
      if (.not. restart) call setup_b_field (1e10, 1e10, 2)
c      if (.not. restart) call setup_b_field (1e12, 1e12, 2)
c      if (.not. restart) call setup_b_field (1e-4, 1e-4, 3)
#endif /* CFC_MHD */         
      call initialize_conserved_quantities

#if  defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined (OPEN_MP_3D))
!$OMP END PARALLEL
#endif

      cfl_factor = config%cfl

      call comp_legendre
      call hilfsgroessen_r

#ifdef EXTRACT_GW
      qijq1(:,:) = 0.0_rk
      qijq2(:,:) = 0.0_rk
#endif /* EXTRACT_GW */

      return
      
      END SUBROUTINE init_cfc_hydro

c     ==================================================================



      END MODULE gr_hyd_init

c     ==================================================================
