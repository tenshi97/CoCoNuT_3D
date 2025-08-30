c     ==================================================================

      MODULE gr_input

c     ==================================================================

      IMPLICIT NONE

      CONTAINS



c     ==================================================================

      SUBROUTINE read_coconut_data

c     ==================================================================

      USE precision

      USE size_cfc
      USE grid_cfc
      USE parameters_cfc
      USE hydro_primitives_cfc
      USE conserved_cfc
      USE metric_cfc
      USE boundary_cfc
      USE perm_aux_cfc
      USE temp_aux_cfc
      USE nutra_cfc
      USE induction
      
#ifdef WRITE_BINARY_OUTPUT
      USE dataformat
#endif
      USE charac

#ifndef PROGRAM_remap
      USE metric
      USE cons_check
      USE recover_prim_vars
      USE renormalize
      USE gr_hyd_init
      USE gr_initial_setup
#endif
#ifdef COARSENING      
      USE mesh_coarsening
#endif
      
      use configure 
      USE gfloat_hy, ONLY: time

      IMPLICIT NONE

      type(datafile) :: fd

      integer (kind=ik) :: i, j, k, l
      logical :: perturb

#ifdef READ_OLD_RESTART_FILES 
      real (kind=ik) ::
     &     phi_last_metric (0:m+1,0:n+1,1:o,1:3),
     &     alpha_last_metric (0:m+1,0:n+1,1:o,1:3),
     &     beta_up_1_last_metric (0:m+1,0:n+1,1:o,1:3),
     &     beta_up_2_last_metric (0:m+1,0:n+1,1:o,1:3),
     &     beta_up_3_last_metric (0:m+1,0:n+1,1:o,1:3)
#endif /* READ_OLD_RESTART_FILES */

      logical :: ler

#if 1
      real (kind=rk), allocatable :: vpx (:,:,:), vpy (:,:,:)
#endif

#ifdef CFC_TRANSPORT
      call init_cfc_hydro(.true.)
#endif

      pi = asin (1.0_rk) * 2.0_rk

      print*
      print*, 'Reading in restart data.'

#ifndef READ_OLD_RESTART_FILES

      print *,'Opening file'
      call open_data_file(fd, trim(rstfil_cfc))

      print *,'d_cap_hat'      
      call readq(fd,"d_cap_hat",d_cap_hat (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'s_1_hat'
      call readq(fd,"s_1_hat",s_1_hat (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'s_2_hat'
      call readq(fd,"s_2_hat",s_2_hat (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'s_3_hat'
      call readq(fd,"s_3_hat",s_3_hat (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'tau_hat'
      call readq(fd,"tau_hat",tau_hat (1:m,n_s:n_e,o_s:o_e),    
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'xnnu'
      call readq(fd,"xnnu",xnnu (1:m,n_s:n_e,o_s:o_e,1:config%qn),
     &     (/1,m, n_s,n_e, o_s,o_e, 1,config%qn/), 
     &     (/1,m, 1,n, 1,o, 1,config%qn/))
      print *,'phi'
      call readq(fd,"phi",phi (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'alpha'
      call readq(fd,"alpha",alpha (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'beta_up_1'
      call readq(fd,"beta_up_1",beta_up_1 (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'beta_up_2'
      call readq(fd,"beta_up_2",beta_up_2 (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'beta_up_3'
      call readq(fd,"beta_up_3",beta_up_3 (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'enu'
      call readq(fd,"enu",enu(1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'fnu'
      call readq(fd,"fnu",fnu(1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'m_grav_1_ini'
      call readq(fd,"m_grav_1_ini",m_grav_1_ini)
      print *,'m_rest_ini'
      call readq(fd,"m_rest_ini",m_rest_ini)

#ifdef CFC_MHD
      print *,'b_cap_1_hat'
      call readq(fd,"b_cap_1_hat",b_cap_1_hat (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'b_cap_2_hat'
      call readq(fd,"b_cap_2_hat",b_cap_2_hat (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'b_cap_3_hat'
      call readq(fd,"b_cap_3_hat",b_cap_3_hat (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))

      print *,'b_1'
      call readq(fd,"b_1",b_cap_1_hat_old (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'b_2'
      call readq(fd,"b_2",b_cap_2_hat_old (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      print *,'b_3'
      call readq(fd,"b_3",b_cap_3_hat_old (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
c      print *,'psi'
      psi_old (1 : m, n_s : n_e, o_s : o_e) = 0.0_rk
      call readq(fd,"psi",psi_old (1:m,n_s:n_e,o_s:o_e),
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))

      b_1 (1 : m, n_s : n_e, o_s : o_e) =
     &     b_cap_1_hat_old (1 : m, n_s : n_e, o_s : o_e)
      b_2 (1 : m, n_s : n_e, o_s : o_e) =
     &     b_cap_2_hat_old (1 : m, n_s : n_e, o_s : o_e)
      b_3 (1 : m, n_s : n_e, o_s : o_e) =
     &  b_cap_3_hat_old (1 : m, n_s : n_e, o_s : o_e)

      psi (1 : m, n_s : n_e, o_s : o_e) =
     &  psi_old (1 : m, n_s : n_e, o_s : o_e)
      
#endif /* CFC_MHD */


      print *,'Closing file'
      call close_data_file(fd)

#else /* READ_OLD_RESTART_FILES */

      open (11, file = rstfil_cfc, form = 'unformatted')
 
      read (11) d_cap_hat (1:m,1:n,1:o),
     &     s_1_hat (1:m,1:n,1:o),
     &     s_2_hat (1:m,1:n,1:o),
     &     s_3_hat (1:m,1:n,1:o),
     &     tau_hat (1:m,1:n,1:o),
     &     xnnu (1:m,1:n,1:o,1:config%qn),
     &     phi (0:m+1,0:n+1,1:o),
     &     alpha (0:m+1,0:n+1,1:o),
     &     beta_up_1 (0:m+1,0:n+1,1:o),
     &     beta_up_2 (0:m+1,0:n+1,1:o),
     &     beta_up_3 (0:m+1,0:n+1,1:o),
     &     phi_last_metric (0:m+1,0:n+1,1:o,1:3),
     &     alpha_last_metric (0:m+1,0:n+1,1:o,1:3),
     &     beta_up_1_last_metric (0:m+1,0:n+1,1:o,1:3),
     &     beta_up_2_last_metric (0:m+1,0:n+1,1:o,1:3),
     &     beta_up_3_last_metric (0:m+1,0:n+1,1:o,1:3)

      read (11) enu(1:m,1:n,1:o),fnu(1:m,1:n,1:o)

      read (11) m_grav_1_ini,m_rest_ini

      close (11)

#endif /* READ_OLD_RESTART_FILES */

      do i=1,config%qn
         d_cap_xnu_hat (1:m,n_s:n_e,o_s:o_e,i)=
     &   d_cap_hat     (1:m,n_s:n_e,o_s:o_e)*
     &   xnnu          (1:m,n_s:n_e,o_s:o_e,i)
      end do
      sqrt_gamma (1:m,n_s:n_e,o_s:o_e) =
     &     phi (1:m,n_s:n_e,o_s:o_e) ** 6


#if 0
      do k = o_s, o_e
         do j = n_s, n_e
            tau_hat (1:m,j,k) = tau_hat (1:m,j,k) +
     &           0.5_rk * psi (1:m,j,k) ** 2
         end do
      end do
#endif
#if 0
c      if (time .LE. 1e-8_rk) then
         PRINT *,'Restarting from multi-D initial model'
         call init_cfc_hydro(.false.)
c      end if
#endif
      call renormalize_mass_fractions

#ifdef ORIGIN
      s_1_hat(1,:,:)=0.0_rk
#endif


#ifdef COARSENING         
      if (config%nsdim .eq. 3) then
         print *,'Setup: mesh coarsening...'
         call setup_coarse_grid
         print *,'Done!'
         call coarse_grid_average (0)
         call prolongation
c         call initialize_conserved_quantities
      endif
#endif

      
#ifndef PROGRAM_remap
      print *,'perm_aux'
      call calculate_perm_aux_values

#if defined(OPEN_MP_2D) || defined (OPEN_MP_3D)
!$OMP PARALLEL
#endif
      print *,'metric'
      call calculate_metric
      print *,'metric calculated'
#if defined(OPEN_MP_2D) || defined (OPEN_MP_3D)
!$OMP END PARALLEL
#endif

      call do_conservation_check

      print *,'recovery'
      call recover_primitives

      
#if 0
      call setup_b_field (1e6, 1e6, 1)
      call initialize_conserved_quantities
      call recover_primitives
#endif

      
#if 0
c     perturbations

      INQUIRE(file='perturb',exist=perturb)
      IF (perturb) THEN
         print *,'perturbations'
         allocate (vpx(m,n_s:n_e,o_s:o_e))
c     allocate (vpy(m,n_s:n_e,o_s:o_e))
c      open (20,file='perturbation.i3e',form='unformatted',
c     &     action='read')
c      read (20) vpx, vpy
c      close(20)
         do k = 1, o_s
         CALL RANDOM_SEED
         end do
         call RANDOM_NUMBER(vpx)
c     call RANDOM_NUMBER(vpy)
         
         do k = o_s, o_e
            do j = n_s, n_e
               s_1_hat (1:m,j,k) = s_1_hat (1:m,j,k) *
     &              (1.0_rk + vpx(1:m,j,k) * 5e-3_rk)
            end do
         end do
c     print *,'vpx',d_cap_hat(80,n_s,:)
c      print *,'vpx2',vpx(80,n_s,:)
      
         deallocate(vpx)

         call recover_primitives

c      call initialize_conserved_quantities 
         END IF
         
#endif
#endif /* PROGRAM_remap */
      
      print *,'Done!'

      return

      END SUBROUTINE read_coconut_data

c     ==================================================================



      END MODULE gr_input

c     ==================================================================
