#undef MC_LIMITER
#define VAN_LEER_LIMITER
#undef SUPERBEE_LIMITER

#define THIRD_ORDER_PROLONGATION
c#undef THIRD_ORDER_PROLONGATION
      
c     ==================================================================

      MODULE mesh_coarsening

c     ==================================================================
c     Adaptive coarsening of the VARPHI grid near the polar axes,
c     designed to avoid overly severe restrictions on the CFL time-step.
c     Based on the idea of P. Cerda-Duran, parallel implementation and
c     treatment of mass fractions by B. Mueller.
c     ------------------------------------------------------------------
      
      CONTAINS

c     ==================================================================

      SUBROUTINE setup_coarse_grid

c     ==================================================================
c     - Determines coarsening levels for different latitudes such that
c       the aspect ratio of the super-cells in the theta-varphi-plane
c       is smaller than 1:2
c     - Sets up communicators for super-cells exceeding the local
c       domain in the MPI case.
c     ------------------------------------------------------------------

      USE precision
      USE conserved_cfc
      USE grid_cfc
      USE mesh_cfc
      USE perm_aux_cfc

      USE mo_mpi

      IMPLICIT NONE

      INTEGER (KIND=ik) :: i, j, k, level, region, istat, jj, kk, ksc
      INTEGER (KIND=ik) :: nps, npe, ops, ope

c#ifdef MPI_HYDRO
      INTEGER (KIND=ik) :: ranks (0 : o / 2 - 1) 
      INTEGER (KIND=ik) :: coords (2)
      INTEGER (KIND=ik) :: mycoords (2)
      INTEGER (KIND=ik) :: dims (2)
      LOGICAL :: periods (2)
      INTEGER (KIND=ik) :: comm_size
      INTEGER (KIND=ik) :: dummy_comm, dummy_group
      INTEGER (KIND=ik) :: ierr
c#endif      

c      ALLOCATE (clevel (n_s : n_e), stat = istat)
      ALLOCATE (clevel (1 : n), stat = istat)
      
      DO j = 1, n
         clevel (j) = INT (- LOG (sinus_theta (j) * delta_varphi /
     &        delta_theta) / LOG (2.0_rk))
         clevel (j) = MAX (clevel (j), 0_ik)
      END DO
      PRINT *,n_s,n_e,'clevel',clevel(:)

c     Count regions with different coarsening levels, distinguish
c     patches with theta coordinate of the relevant MPI process
      region = 1
      region_s = 1
      region_e = 1
      level = clevel (1)
      DO j = 2, n
         IF (clevel (j) .NE. level .OR. MOD (j, n_loc) .EQ. 1) THEN
            level = clevel (j)
            region = region + 1
         END IF
         IF (j .EQ. n_s) region_s = region
         IF (j .EQ. n_e) region_e = region
      END DO

      nregions = region
      PRINT *,'nregions',nregions,region_s,region_e
      
      ALLOCATE (j_i (1 : nregions), stat = istat)
      ALLOCATE (j_f (1 : nregions), stat = istat)
      ALLOCATE (clevel_region (1 : nregions), stat = istat)

c#ifdef MPI_HYDRO
      ALLOCATE (left_root  (region_s : region_e), stat = istat)
      ALLOCATE (right_root (region_s : region_e), stat = istat)
      ALLOCATE (need_intracell  (region_s : region_e), stat = istat)
      ALLOCATE (comm_super_cell (region_s : region_e), stat = istat)
      ALLOCATE (group_super_cell (region_s : region_e), stat = istat)
c#endif


c     Determine boundaries of coarsening regions in THETA direction
      region = 1
      level = clevel (1)
      clevel_region (region) = level
      j_i (region) = 1
      j_f (region) = 1

      DO j = 2, n
         IF (clevel (j) .NE. level .OR. MOD (j, n_loc) .EQ. 1) THEN
            region = region + 1
            level = clevel (j)
            clevel_region (region) = level
            j_i (region) = j
            j_f (region) = j
         ELSE
            j_f (region) = j
         END IF
         PRINT *,'j',j,level,region
      END DO

      PRINT *,'clevel_region',clevel_region(:)
      PRINT *,'clevel_region1',j_i(:)
      PRINT *,'clevel_region2',j_f(:)

      CALL MPI_Cart_get (cart_comm, 2_ik, dims, periods,
     &     mycoords, ierr)

c     loop over processes
      DO region = 1, nregions

         jj = (j_i (region) - 1) / n_loc !1st process coordinate
         
         PRINT *,'region1',region,2**clevel_region(region),
     &        o_loc,nregions,jj
         
         IF (2 ** clevel_region (region) .GT. o_loc) THEN
            
            IF (region .GE. region_s .AND. region .LE. region_e)
     &           need_intracell (region) = .TRUE.

c     loop over super-cells
            DO ksc = 1, o / (2 ** clevel_region (region) / o_loc)

               PRINT *,'b',ksc

               CALL MPI_Cart_coords (cart_comm, myproc, 2_ik, coords,
     &           ierr)
               PRINT *,'c'

               k = (o_s - 1) / 2 ** clevel_region (region)
               comm_size = 
     &              MAX (2 ** clevel_region (region) / o_loc, 1_ik)
               
               DO i = 0, comm_size - 1
                  PRINT *,'subcell',i,o/2-1
                  coords (1) = jj
                  coords (2) = (ksc-1) * comm_size + i
                  CALL MPI_Cart_rank (cart_comm, coords, ranks (i), 
     &                 ierr)
                  PRINT *,'subcell',i,coords(:),k,ranks(i)
               END DO
c               CALL  MPI_Comm_group (cart_comm, cart_group, ierr)

               PRINT *,'level',region,clevel_region(region)
               PRINT *,'Group_incl',comm_size,ranks(0:comm_size-1)
               PRINT *,'cart_group',cart_group, cart_comm

               IF (mycoords (1) .EQ. coords (1) .AND.
     &              mycoords (2) / comm_size .EQ. ksc-1) THEN

c     this process is included in the super-cell
                  CALL MPI_Group_incl (cart_group, comm_size,
     &                 ranks (0:comm_size-1), group_super_cell
     &                 (region), ierr)
                  PRINT *,'comm_create',group_super_cell(region)
                  CALL MPI_Comm_create (cart_comm,
     &                 group_super_cell (region),
     &                 comm_super_cell (region), ierr)

                  left_root  (region) = 0_ik
                  right_root (region) = comm_size - 1_ik
                  PRINT *,'roots',left_root(region),right_root(region)

               ELSE

c     only particapte in the creating the communicator, which does
c     not have to be stored

                  CALL MPI_Group_incl (cart_group, comm_size,
     &                 ranks (0:comm_size-1), dummy_group, ierr)
                  PRINT *,'comm_create',dummy_group
                  CALL MPI_Comm_create (cart_comm, dummy_group,
     &                 dummy_comm, ierr)

               END IF

            END DO !loop over super-cells

         ELSE
          
            IF (region .GE. region_s .AND. region .LE. region_e)
     &           need_intracell (region) = .FALSE.
            
         END IF
                        
c         PRINT *,'need_intracell',need_intracell(region)
      END DO                    !loop over regions
c#endif

      PRINT *,'Mesh coarsening: setup complete!'

      DEALLOCATE (clevel, stat = istat)

c      ALLOCATE (fft_comm (1:dims(1)))


      RETURN

      END SUBROUTINE setup_coarse_grid

c     ==================================================================




c     ==================================================================

      SUBROUTINE coarse_grid_average (mode)

c     ==================================================================
c     - Computes super-cell averages of conserved quantities and mass
c       fractions.
c     - Reconstructs gradients of conserved quantities and mass
c       fractions to recover higher-order accuracy.
c     ------------------------------------------------------------------

      USE precision


      USE conserved_cfc
#ifdef CFC_MHD
      USE hydro_primitives_cfc, ONLY: xnnu, psi
#else
      USE hydro_primitives_cfc, ONLY: xnnu
#endif      
      USE mesh_cfc

      USE mo_mpi

      USE configure

      USE, INTRINSIC :: iso_c_binding

C      INCLUDE 'fftw3-mpi.f'

      IMPLICIT NONE

      INTEGER (KIND=ik) :: mode
      
      REAL (KIND=rk), ALLOCATABLE :: q (:,:,:), qxnu (:,:,:)
      REAL (KIND=rk), ALLOCATABLE :: qr (:,:,:), qrxnu (:,:,:)
      REAL (KIND=rk), ALLOCATABLE :: ql (:,:,:), qlxnu (:,:,:)
      REAL (KIND=rk) :: dxdphi_plus (1:m), dxdphi_minus (1:m),
     &     dxdphi (1:m), cflat (1:m), corr_fac (1:m)
      REAL (KIND=rk) :: dl, dr, slope, scal, coeff

#ifdef CFC_MHD
      INTEGER (KIND=ik), PARAMETER :: nvars = 10
#else
      INTEGER (KIND=ik), PARAMETER :: nvars = 5
#endif      
      
      INTEGER (KIND=ik) :: i, j, k, kk, i_nuc, ji1, jf1, region, istat,
     &     level, ki, kf, km1, kp1, l, nn, nd, width, width0

c#ifdef MPI_HYDRO
      INTEGER (KIND=ik) :: src, dest, ierr, mpistat (MPI_STATUS_SIZE)
      LOGICAL :: at_boundary
c#endif

c      fftw_mpi_local_size_many_1d (o,m*qy_loc, cart_comm, 
c     &     FFTW_FORWARD, FFTW_SCRAMBLED_OUT, ni, i_start, no, o_start)

c      PRINT *,'b_cap_hat (vorher)',b_cap_1_hat(549,n_s,o_s),
c     &     b_cap_2_hat(549,n_s,o_s), b_cap_3_hat(549,n_s,o_s),
c     &     psi(549,n_s,o_s)
c      PRINT *,'tau_hat (vorher)',tau_hat(12,n_s,o_s)

      DO region = region_s, region_e
c         print *,'region2',region,2**clevel_region(region),
c     &        need_intracell(region),o_loc
         IF (clevel_region (region) .NE. 0) THEN
            
            ji1 = j_i (region)
            jf1 = j_f (region)
            !print *,'ji1,jf1',j_i(region),j_f(region)
            
            ALLOCATE (q (1:m, ji1:jf1, 1:nvars), stat = istat)
            ALLOCATE (qxnu (1:m, ji1:jf1, 1:config%qn), stat = istat)
            ALLOCATE (ql (1:m, ji1:jf1, 1:nvars), stat = istat)
            ALLOCATE (qlxnu (1:m, ji1:jf1, 1:config%qn), stat = istat)
            ALLOCATE (qr (1:m, ji1:jf1, 1:nvars), stat = istat)
            ALLOCATE (qrxnu (1:m, ji1:jf1, 1:config%qn), stat = istat)
            
c     Determine number of super-cells on local domain
            width0 = 2 ** clevel_region (region)
            width  = MIN (width0, qz_proc)
            scal = 1.0_rk / REAL (width0, kind = rk)
            nn = o_loc / width 

c     Number of domains sharing the same super-cell
            nd = MAX (1_ik, width0 / o_loc) 

c     Loop over super-cells (computation of cell averages)
            DO kk = 0, nn - 1
               !print *,'supercell kk',kk
               q    (:,:,:) = 0.0_rk
               qxnu (:,:,:) = 0.0_rk

c     Sum conserved variables over interior cells and average
c     Note: we assume that delta_varphi is constant!
               DO k = o_s + width * kk, o_s + width * (kk + 1) - 1
                  !print *,'k',k,ji1,jf1
                  DO j = ji1, jf1
                     q (:,j,1) = q (:,j,1) + d_cap_hat (:,j,k)
                     q (:,j,2) = q (:,j,2) + s_1_hat   (:,j,k)
                     q (:,j,3) = q (:,j,3) + s_2_hat   (:,j,k)
                     q (:,j,4) = q (:,j,4) + s_3_hat   (:,j,k)
                     q (:,j,5) = q (:,j,5) + tau_hat   (:,j,k)
#ifdef CFC_MHD
                     q (:,j,6) = q (:,j,6) + b_cap_1_hat   (:,j,k)
                     q (:,j,7) = q (:,j,7) + b_cap_2_hat   (:,j,k)
                     q (:,j,8) = q (:,j,8) + b_cap_3_hat   (:,j,k)
                     q (:,j,9) = q (:,j,9) + psi         (1:m,j,k)
                     q (:,j,10) = q (:,j,10) + eps_aux     (:,j,k)
#endif /* CFC_MHD */                     
                     DO i_nuc = 1, config%qn
                        qxnu (:,j,i_nuc) = qxnu (:,j,i_nuc) +
     &                       d_cap_xnu_hat (:,j,k,i_nuc)
                     END DO
                  END DO
               END DO           !summation loop

c#ifdef MPI_HYDRO
               IF (need_intracell (region)) THEN
                  !print *,'MPI_Allreduce'
                  !print *,'q(80,ji1,3)',ji1,q(80,ji1,3)
                  CALL MPI_AllReduce (q, ql,
     &            m * (jf1 - ji1 + 1) * nvars, MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, comm_super_cell (region), ierr)
                  q = ql
                  CALL MPI_AllReduce (qxnu, qlxnu,
     &            m * (jf1 - ji1 + 1) * config%qn, MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, comm_super_cell (region), ierr)
                  qxnu = qlxnu
               END IF
c#endif
               q    (:,:,:) = q    (:,:,:) / 
     &              REAL (width0, KIND=rk)
               qxnu (:,:,:) = qxnu (:,:,:) / 
     &              REAL (width0, KIND=rk)
               !print *,q(80,ji1,3),qxnu(80,ji1,config%qn),'reduce'

               DO k = o_s + width * kk, o_s + width * (kk + 1) - 1

                  d_cap_hat (:,ji1:jf1,k) = q (:,:,1)
                  s_1_hat   (:,ji1:jf1,k) = q (:,:,2)
                  s_2_hat   (:,ji1:jf1,k) = q (:,:,3)
                  s_3_hat   (:,ji1:jf1,k) = q (:,:,4)
                  tau_hat   (:,ji1:jf1,k) = q (:,:,5)
                  d_cap_xnu_hat (:,ji1:jf1,k,1:config%qn) =
     &                 qxnu (:,ji1:jf1,1:config%qn)
#ifdef CFC_MHD
                  b_cap_1_hat   (:,ji1:jf1,k) = q (:,:,6)
                  b_cap_2_hat   (:,ji1:jf1,k) = q (:,:,7)
                  b_cap_3_hat   (:,ji1:jf1,k) = q (:,:,8)
                  psi         (1:m,ji1:jf1,k) = q (:,:,9)
                  eps_aux       (:,ji1:jf1,k) = q (:,:,10)
#endif /* CFC_MHD */                  
                  DO i_nuc = 1, config%qn
                     xnnu (1:m,ji1:jf1,k,i_nuc) =
     &                    d_cap_xnu_hat (1:m,ji1:jf1,k,i_nuc) /
     &                    d_cap_hat (1:m,ji1:jf1,k)
                  END DO

               END DO
                  
            END DO              !loop over super-cells

            IF (mode .EQ. 0) GOTO 2000
            
c     2nd loop over super-cells (pre-reconstruction)
c            GO TO 2000

            DO kk = 0, nn - 1
               !print *,'2nd loop, kk',kk

c#ifdef MPI_HYDRO
               at_boundary = 
     &              ((o_s - 1) / 2 ** clevel_region (region) .NE.
     &              (o_s - 2) / 2 ** clevel_region (region)) .OR.
     &              o_s .EQ. 1
               !print *,at_boundary

               IF (kk. EQ. 0_ik .AND. at_boundary) THEN !first local super-cell
c     send to right neighbour, receive from left neighbour
c     (note of caution: sending data from super-cell nn, i.e. not the one
c     we are currently working on!)
                  q (:,:,1) = d_cap_hat (:,ji1:jf1,o_e)
                  q (:,:,2) = s_1_hat   (:,ji1:jf1,o_e) / q (:,:,1)
                  q (:,:,3) = s_2_hat   (:,ji1:jf1,o_e) / q (:,:,1)
                  q (:,:,4) = s_3_hat   (:,ji1:jf1,o_e) / q (:,:,1)
                  q (:,:,5) = tau_hat   (:,ji1:jf1,o_e) / q (:,:,1)
#ifdef CFC_MHD
                  q (:,:,6) = b_cap_1_hat   (:,ji1:jf1,o_e)
                  q (:,:,7) = b_cap_2_hat   (:,ji1:jf1,o_e)
                  q (:,:,8) = b_cap_3_hat   (:,ji1:jf1,o_e)
                  q (:,:,9) = psi         (1:m,ji1:jf1,o_e)
                  q (:,:,10) = eps_aux      (:,ji1:jf1,o_e) / q (:,:,1)
#endif /* CFC_MHD */
                  qxnu (:,:,1:config%qn) = 
c     &                 d_cap_hat_xnu (1:m,ji1:jf1,o_e,1:config%qn)
     &                 xnnu (1:m,ji1:jf1,o_e,1:config%qn)
                  CALL MPI_Cart_shift (cart_comm, 1_ik, nd,
     &                 src, dest, ierr)
                  !print *,'MPI_Sendrecv',src,dest,left_root(region)
                  IF (MOD (src / nd, 2_ik) .EQ. 1_ik) THEN
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5000_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5001_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Recv (ql, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5000_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qlxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5001_ik,
     &                    cart_comm, mpistat, ierr)
                  ELSE
                     CALL MPI_Recv (ql, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5000_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qlxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5001_ik,
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5000_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5001_ik,
     &                    cart_comm, ierr)
                  END IF

                  IF (need_intracell (region)) THEN
c                     !print *,'MPI_Bcast'
                     CALL MPI_Bcast (ql,    m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region),
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qlxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                  END IF

               ELSE IF (kk .EQ. 0_ik .AND. .NOT. at_boundary) THEN

                  IF (need_intracell (region)) THEN
c                     !print *,'MPI_Bcast'
                     CALL MPI_Bcast (ql,    m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region),
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qlxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region),
     &                    comm_super_cell (region), ierr)
                  END IF

               ELSE             !data is process-local
c#endif /* MPI_HYDRO */

                  k = o_s + width * kk
                  km1 = k - 1
                  IF (km1 .lt. o_s) km1 = o_e
                  
                  ql (:,:,1) = d_cap_hat (:,ji1:jf1,km1)
                  ql (:,:,2) = s_1_hat   (:,ji1:jf1,km1) / ql (:,:,1)
                  ql (:,:,3) = s_2_hat   (:,ji1:jf1,km1) / ql (:,:,1)
                  ql (:,:,4) = s_3_hat   (:,ji1:jf1,km1) / ql (:,:,1)
                  ql (:,:,5) = tau_hat   (:,ji1:jf1,km1) / ql (:,:,1)
#ifdef CFC_MHD
                  ql (:,:,6) = b_cap_1_hat   (:,ji1:jf1,km1)
                  ql (:,:,7) = b_cap_2_hat   (:,ji1:jf1,km1)
                  ql (:,:,8) = b_cap_3_hat   (:,ji1:jf1,km1)
                  ql (:,:,9) = psi        (1:m,ji1:jf1,km1)
                  ql (:,:,10) = eps_aux      (:,ji1:jf1,km1) / q (:,:,1)
#endif /* CFC_MHD */
                  qlxnu (:,:,1:config%qn) = 
c     &                 d_cap_hat_xnu (1:m,ji1:jf1,km1,1:config%qn)
     &                 xnnu (1:m,ji1:jf1,km1,1:config%qn)

c#ifdef MPI_HYDRO 
               END IF

               IF (kk. EQ. nn - 1 .AND.  at_boundary) THEN !last local super-cell

c     send to left neighbour, receive from right neighbour
c     (note of caution: sending data from super-cell 1, i.e. no the one
c     we are currently working on!)
                  q (:,:,1) = d_cap_hat (:,ji1:jf1,o_s)
                  q (:,:,2) = s_1_hat   (:,ji1:jf1,o_s) / q (:,:,1)
                  q (:,:,3) = s_2_hat   (:,ji1:jf1,o_s) / q (:,:,1)
                  q (:,:,4) = s_3_hat   (:,ji1:jf1,o_s) / q (:,:,1)
                  q (:,:,5) = tau_hat   (:,ji1:jf1,o_s) / q (:,:,1)
#ifdef CFC_MHD
                  q (:,:,6) = b_cap_1_hat   (:,ji1:jf1,o_s)
                  q (:,:,7) = b_cap_2_hat   (:,ji1:jf1,o_s)
                  q (:,:,8) = b_cap_3_hat   (:,ji1:jf1,o_s)
                  q (:,:,9) = psi         (1:m,ji1:jf1,o_s)
                  q (:,:,10) = eps_aux      (:,ji1:jf1,o_s) / q (:,:,1)
#endif /* CFC_MHD */
                  qxnu (:,:,1:config%qn) = 
c     &                 d_cap_hat_xnu (1:m,ji1:jf1,o_s,1:config%qn)
     &                 xnnu (1:m,ji1:jf1,o_s,1:config%qn)
                  CALL MPI_Cart_shift (cart_comm, 1_ik, -nd,
     &                 src, dest, ierr)
c                  print *,'MPI_sendrecv',src,dest,right_root(region),
c     &                 MOD (src / nd, 2_ik)
                  IF (MOD (src / nd, 2_ik) .EQ. 1_ik) THEN
c                     !print *,'send'
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5002_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5003_ik,
     &                    cart_comm, ierr)
c                     !print *,'recv'
                     CALL MPI_Recv (qr, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5002_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qrxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5003_ik,
     &                    cart_comm, mpistat, ierr)
                  ELSE
c                     !print *,'recv'
                     CALL MPI_Recv (qr, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5002_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qrxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5003_ik,
     &                    cart_comm, mpistat, ierr)
c                     !print *,'send'
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5002_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5003_ik,
     &                    cart_comm, ierr)
                  END IF

                  IF (need_intracell (region)) THEN
c                     !print *,'MPI_bcast'
                     CALL MPI_Bcast (qr, m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qrxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                  END IF

c                  !print *,'communication finished'

               ELSE IF (kk .EQ. nn - 1 .AND. .NOT. at_boundary) THEN

                  IF (need_intracell (region)) THEN
c                     !print *,'MPI_bcast'
                     CALL MPI_Bcast (qr,    m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qrxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                  END IF
c                  !print *,'communication finished'
               ELSE             !data is process-local
c#endif /* MPI_HYDRO */

                  k = o_s + width * (kk + 1) - 1
                  kp1 = k + 1
                  IF (kp1 .gt. o_e) kp1 = o_s
                  
                  qr    (:,:,1) = d_cap_hat (:,ji1:jf1,kp1) 
                  qr    (:,:,2) = s_1_hat   (:,ji1:jf1,kp1) / qr (:,:,1)
                  qr    (:,:,3) = s_2_hat   (:,ji1:jf1,kp1) / qr (:,:,1)
                  qr    (:,:,4) = s_3_hat   (:,ji1:jf1,kp1) / qr (:,:,1)
                  qr    (:,:,5) = tau_hat   (:,ji1:jf1,kp1) / qr (:,:,1)
#ifdef CFC_MHD
                  qr    (:,:,6) = b_cap_1_hat(:,ji1:jf1,kp1)
                  qr    (:,:,7) = b_cap_2_hat(:,ji1:jf1,kp1)
                  qr    (:,:,8) = b_cap_3_hat(:,ji1:jf1,kp1)
                  qr    (:,:,9) = psi      (1:m,ji1:jf1,kp1)
                  qr    (:,:,10) = eps_aux   (:,ji1:jf1,kp1) / q (:,:,1)
#endif /* CFC_MHD */                  
                  qrxnu (:,:,1:config%qn) = 
c     &                 d_cap_hat_xnu    (1:m,ji1:jf1,kp1,1:config%qn)
     &                 xnnu    (1:m,ji1:jf1,kp1,1:config%qn)

c#ifdef MPI_HYDRO
               END IF
c#endif /* MPI_HYDRO */

c     Copy conserved variables and mass fractions into temporary arrays
c     (simplifies reconstruction loop)
               q (:,:,1) = d_cap_hat (:,ji1:jf1,o_s+width*kk)
               q (:,:,2) = s_1_hat   (:,ji1:jf1,o_s+width*kk)/q (:,:,1)
               q (:,:,3) = s_2_hat   (:,ji1:jf1,o_s+width*kk)/q (:,:,1)
               q (:,:,4) = s_3_hat   (:,ji1:jf1,o_s+width*kk)/q (:,:,1)
               q (:,:,5) = tau_hat   (:,ji1:jf1,o_s+width*kk)/q (:,:,1)
#ifdef CFC_MHD               
               q (:,:,6) = b_cap_1_hat(:,ji1:jf1,o_s+width*kk)
               q (:,:,7) = b_cap_2_hat(:,ji1:jf1,o_s+width*kk)
               q (:,:,8) = b_cap_3_hat(:,ji1:jf1,o_s+width*kk)
               q (:,:,9) = psi      (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,10) = eps_aux   (:,ji1:jf1,o_s+width*kk)/q (:,:,1)
#endif /* CFC_MHD */               
               qxnu (:,:,1:config%qn) =
c     &              d_cap_hat_xnu (1:m,ji1:jf1,o_s+width*kk,1:config%qn)
     &              xnnu (1:m,ji1:jf1,o_s+width*kk,1:config%qn)

c     ------------------------------------------------------------------
c     Determine slope in super-cell and reconstruct sub-cells (q and
c     qxnu are overwritten in the process and contain the slopes in the end)

               DO j = ji1, jf1
                  !print *,'slopes',j,ji1,jf1
                  DO l = 1, nvars   !conserved variables
                     DO i = 1, m
                        dl = q  (i, j, l) - ql (i, j, l)
                        dr = qr (i, j, l) - q  (i, j, l)
                        ql (i, j, l) = dl
                        qr (i, j, l) = dr
#ifdef MC_LIMITER
c     Monotonized central limiter
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = MIN (ABS (0.5_rk * (dl + dr)),
     &                          2.0_rk * MIN (ABS (dl), ABS (dr)))
                           slope  = SIGN (slope, dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef VAN_LEER_LIMITER
c     van Leer's smooth slope limiter 
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = 2.0 * dl * dr / (dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef SUPERBEE_LIMITER
                        slope =
     &                       (sign (1.0_rk, dl) + sign (1.0_rk, dr)) *
     &                       min(abs (dl), abs (dr),
     &                       0.5_rk * max (abs (dl), abs (dr)))
#endif

                        q (i, j, l) = slope !* scal
                     END DO
                  END DO !l-loop

                  corr_fac (1:m) = q (1 :m, j, 1) / !needed later (-> mass fractions)
     &                 (12.0_rk * d_cap_hat (1:m,j,o_s+width*kk))

c     Additional slope limiting
                  DO l = 2, 5
                     DO i = 1, m
                        q (i, j, l) = q (i, j, l) /
     &                       (1.0_rk + 2.0_rk * ABS (corr_fac (i)))
                     END DO
                  END DO

                  DO i_nuc = 1, config%qn !mass fractions and Y_e
                     DO i = 1, m

                        dl = qxnu  (i, j, i_nuc) - qlxnu (i, j, i_nuc)
                        dr = qrxnu (i, j, i_nuc) - qxnu  (i, j, i_nuc)
#ifdef MC_LIMITER
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = MIN (ABS (0.5_rk * (dl + dr)),
     &                          2.0_rk * MIN (ABS (dl), ABS (dr)))
                           slope  = SIGN (slope, dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef VAN_LEER_LIMITER
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = 2.0 * dl * dr / (dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef SUPERBEE_LIMITER
                        slope =
     &                       (sign (1.0_rk, dl) + sign (1.0_rk, dr)) *
     &                       min(abs (dl), abs (dr),
     &                       0.5_rk * max (abs (dl), abs (dr)))
#endif

c                        qxnu (i, j, i_nuc) = slope !* scal
                        
c     Additional slope limiting to ensure positivity of X_i
c                        qxnu (i, j, i_nuc) = SIGN (MIN (ABS (slope),
c     &                       qxnu (i, j, i_nuc) /
c     &                       (0.5_rk + ABS (corr_fac (i)))),
c     &                       slope)
                        qxnu (i, j, i_nuc) = slope /
     &                       (1.0_rk + 2.0_rk * ABS (corr_fac (i))) 
c                        qxnu (i, j, i_nuc) = slope /
c     &                       MAX (1.0_rk,
c     &                       ABS (slope) / qxnu (i, j, i_nuc) *
c     &                       (0.5_rk + ABS (corr_fac (i))))
c                        qxnu (i, j, i_nuc) = 0.0_rk !debugging
c                        qxnu (i, j, i_nuc) = 0.0_rk

                     END DO
                  END DO !species-loop

c     ------------------------------------------------------------------
c     Linear reconstruction of conserved variables

                  DO k = o_s + width * kk, o_s + width * (kk + 1) - 1_ik

                     l = MOD (k - 1, width0)

c     The integral of the pre-reconstructed variable U over sub-cells
c     can be expressed analytically as <U>_supercell + coeff * dU/dvarphi:
                     coeff = (1_ik + 2_ik * l - width0) /
     &                    REAL (2_ik * width0, KIND=rk)

                     DO i = 1, m

                        s_1_hat   (i, j, k) = s_1_hat   (i, j, k) /
     &                       d_cap_hat (i, j, k)
                        s_2_hat   (i, j, k) = s_2_hat   (i, j, k) /
     &                       d_cap_hat (i, j, k)
                        s_3_hat   (i, j, k) = s_3_hat   (i, j, k) /
     &                       d_cap_hat (i, j, k)
                        tau_hat   (i, j, k) = tau_hat   (i, j, k) /
     &                       d_cap_hat (i, j, k)

                        d_cap_hat (i, j, k) = d_cap_hat (i, j, k) +
     &                       coeff * q (i, j, 1)

                        s_1_hat   (i, j, k) = (s_1_hat   (i, j, k) +
     &                       (coeff - corr_fac (i)) * q (i, j, 2)) *
     &                       d_cap_hat (i, j, k)
                        s_2_hat   (i, j, k) = (s_2_hat   (i, j, k) +
     &                       (coeff - corr_fac (i)) * q (i, j, 3)) *
     &                       d_cap_hat (i, j, k)
                        s_3_hat   (i, j, k) = (s_3_hat   (i, j, k) +
     &                       (coeff - corr_fac (i)) * q (i, j, 4)) *
     &                       d_cap_hat (i, j, k)
                        tau_hat   (i, j, k) = (tau_hat   (i, j, k) +
     &                       (coeff - corr_fac (i)) * q (i, j, 5)) *
     &                       d_cap_hat (i, j, k)

#ifdef CFC_MHD
                        b_cap_1_hat (i, j, k) = b_cap_1_hat (i, j, k) +
     &                       coeff * q (i, j, 6)
                        b_cap_2_hat (i, j, k) = b_cap_2_hat (i, j, k) +
     &                       coeff * q (i, j, 7)
                        b_cap_3_hat (i, j, k) = b_cap_3_hat (i, j, k) +
     &                       coeff * q (i, j, 8)
                        psi          (i, j, k) = psi         (i, j, k) +
     &                       coeff * q (i, j, 9)
                        eps_aux (i, j, k) = (eps_aux   (i, j, k) /
     &                       d_cap_hat (i, j, k) +
     &                       (coeff - corr_fac (i)) * q (i, j, 10)) *
     &                       d_cap_hat (i, j, k)
#endif /* CFC_MHD */                        
                     END DO

                  END DO        !loop over sub-cells

c     ------------------------------------------------------------------
c     Mass fractions:
c     Step 1: Ensure that  \sum_i dX_i/dvarphi = 1 by rescaling
                  dxdphi_plus  (1:m) = 0.0_rk
                  dxdphi_minus (1:m) = 0.0_rk
                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        IF (qxnu (i, j, i_nuc) .GT. 0.0_rk) THEN
                           dxdphi_plus (i) = dxdphi_plus (i) +
     &                          qxnu (i, j, i_nuc)
                        ELSE
                           dxdphi_minus (i) = dxdphi_minus (i) +
     &                          qxnu (i, j, i_nuc)
                        END IF
                     END DO
                  END DO
                  
                  dxdphi (:) = dxdphi_plus (:) + dxdphi_minus (:)
                  WHERE (dxdphi .GT. 0.0_rk)
                     cflat = - dxdphi_minus / dxdphi_plus
                  ELSE WHERE (dxdphi .LT. 0.0_rk)
                     cflat = - dxdphi_plus  / dxdphi_minus
                  ELSE WHERE
                     cflat = 1.0_rk
                  END WHERE

                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        IF (dxdphi (i) .GT. 0.0_rk) THEN
                           qxnu (i, j, i_nuc) = MAX (cflat (i) *
     &                          qxnu (i, j, i_nuc), 0.0_rk) +
     &                          MIN (qxnu (i, j, i_nuc), 0.0_rk) 
                        ELSE
                           qxnu (i, j, i_nuc) = MIN (cflat (i) *
     &                          qxnu (i, j, i_nuc), 0.0_rk) +
     &                          MAX (qxnu (i, j, i_nuc), 0.0_rk) 
                        END IF
                     END DO
                  END DO

c     Step 2: Linear reconstruction of mass fractions
                  DO k = o_s + width * kk, o_s + width * (kk + 1) - 1_ik

                     l = MOD (k - 1, width0)

                     coeff = (1_ik + 2_ik * l - width0) /
     &                    REAL (2_ik * width0, KIND=rk)

                     DO i_nuc = 1, config%qn
                        DO i = 1, m
                           xnnu (i, j, k, i_nuc) =
     &                          xnnu  (i, j, k, i_nuc) +
     &                          (coeff - corr_fac (i)) *
     &                          qxnu (i, j, i_nuc)
                           d_cap_xnu_hat (i, j, k, i_nuc) =
     &                          xnnu  (i, j, k, i_nuc) *
     &                          d_cap_hat (i, j, k) 
                        END DO
                     END DO

                  END DO        !loop over sub-cells (varphi direction)
               END DO           !loop over theta-index


            END DO              !loop over super-cells
            
c            CALL MPI_Barrier (ierr)

 2000       CONTINUE

            DEALLOCATE (q, qxnu, ql, qr, qlxnu, qrxnu)
c            !print *,'deallocation'

         END IF

      END DO !loop over regions

c      PRINT *,'b_cap_hat (nachher)',b_cap_1_hat(549,n_s,o_s),
c     &     b_cap_2_hat(549,n_s,o_s), b_cap_3_hat(549,n_s,o_s),
c     &     psi(549,n_s,o_s)
c      PRINT *,'tau_hat (nachher)',tau_hat(12,n_s,o_s)

      RETURN

      END SUBROUTINE coarse_grid_average

c     ==================================================================


      
c     ==================================================================

      SUBROUTINE prolongation

c     ==================================================================
c     - Reconstructs gradients of primitive variables
c       to recover higher-order accuracy.
c     ------------------------------------------------------------------

      USE precision
      USE hydro_primitives_cfc
      USE mesh_cfc
      USE mo_mpi

      USE configure

      USE, INTRINSIC :: iso_c_binding

C      INCLUDE 'fftw3-mpi.f'

      IMPLICIT NONE

#ifdef CFC_MHD
      INTEGER (KIND=ik), PARAMETER :: nvars = 11
#else
      INTEGER (KIND=ik), PARAMETER :: nvars = 7
#endif      
      REAL (KIND=rk), ALLOCATABLE :: q (:,:,:), qxnu (:,:,:)
      REAL (KIND=rk), ALLOCATABLE :: qr (:,:,:), qrxnu (:,:,:)
      REAL (KIND=rk), ALLOCATABLE :: ql (:,:,:), qlxnu (:,:,:)
      REAL (KIND=rk) :: splus (1:m), sminus (1:m),
     &     ssum (1:m), cflat (1:m)
      REAL (KIND=rk) :: a0 (1:nvars), a1 (1:nvars), a2 (1:nvars)
      REAL (KIND=rk) :: dl, dr, slope, scal, coeff, coeff2
      
      INTEGER (KIND=ik) :: i, j, k, jk, kk, i_nuc, ji1, jf1, region,
     &     istat, level, ki, kf, km1, kp1, l, nn, nd, width, width0

c#ifdef MPI_HYDRO
      INTEGER (KIND=ik) :: src, dest, ierr, mpistat (MPI_STATUS_SIZE)
      LOGICAL :: at_boundary
c#endif

      DO region = region_s, region_e
         IF (clevel_region (region) .NE. 0) THEN
            
            ji1 = j_i (region)
            jf1 = j_f (region)
            
            ALLOCATE (q (1:m, ji1:jf1, 1:nvars), stat = istat)
            ALLOCATE (qxnu (1:m, ji1:jf1, 1:config%qn), stat = istat)
            ALLOCATE (ql (1:m, ji1:jf1, 1:nvars), stat = istat)
            ALLOCATE (qlxnu (1:m, ji1:jf1, 1:config%qn), stat = istat)
            ALLOCATE (qr (1:m, ji1:jf1, 1:nvars), stat = istat)
            ALLOCATE (qrxnu (1:m, ji1:jf1, 1:config%qn), stat = istat)
            
c     Determine number of super-cells on local domain
            width0 = 2 ** clevel_region (region)
            width  = MIN (width0, qz_proc)
            scal = 1.0_rk / REAL (width0, kind = rk)
            nn = o_loc / width 

c     Number of domains sharing the same super-cell
            nd = MAX (1_ik, width0 / o_loc) 

c     loop over super-cells (pre-reconstruction)

            DO kk = 0, nn - 1

c#ifdef MPI_HYDRO
               at_boundary = 
     &              ((o_s - 1) / 2 ** clevel_region (region) .NE.
     &              (o_s - 2) / 2 ** clevel_region (region)) .OR.
     &              o_s .EQ. 1

               IF (kk. EQ. 0_ik .AND. at_boundary) THEN !first local super-cell
c     send to right neighbour, receive from left neighbour
c     (note of caution: sending data from super-cell nn, i.e. not the one
c     we are currently working on!)                  
                  q (:,:,1) = rho   (1:m,ji1:jf1,o_e)
                  q (:,:,2) = v_1   (1:m,ji1:jf1,o_e)
                  q (:,:,3) = v_2   (1:m,ji1:jf1,o_e)
                  q (:,:,4) = v_3   (1:m,ji1:jf1,o_e)
                  q (:,:,5) = eps   (1:m,ji1:jf1,o_e)
                  q (:,:,6) = p     (1:m,ji1:jf1,o_e)
                  q (:,:,7) = gamm  (1:m,ji1:jf1,o_e)
#ifdef CFC_MHD
                  q (:,:,8)  = b_1 (1:m,ji1:jf1,o_e)
                  q (:,:,9)  = b_2 (1:m,ji1:jf1,o_e)
                  q (:,:,10) = b_3 (1:m,ji1:jf1,o_e)
                  q (:,:,11) = psi (1:m,ji1:jf1,o_e)
#endif /* CFC_MHD */
                  qxnu (:,:,1:config%qn) = 
     &                 xnnu (1:m,ji1:jf1,o_e,1:config%qn)
                  CALL MPI_Cart_shift (cart_comm, 1_ik, nd,
     &                 src, dest, ierr)
                  IF (MOD (src / nd, 2_ik) .EQ. 1_ik) THEN
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5000_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5001_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Recv (ql, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5000_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qlxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5001_ik,
     &                    cart_comm, mpistat, ierr)
                  ELSE
                     CALL MPI_Recv (ql, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5000_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qlxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5001_ik,
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5000_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5001_ik,
     &                    cart_comm, ierr)
                  END IF

                  IF (need_intracell (region)) THEN
                     CALL MPI_Bcast (ql,    m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region),
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qlxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                  END IF

               ELSE IF (kk .EQ. 0_ik .AND. .NOT. at_boundary) THEN

                  IF (need_intracell (region)) THEN
                     CALL MPI_Bcast (ql,    m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region),
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qlxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region),
     &                    comm_super_cell (region), ierr)
                  END IF

               ELSE             !data is process-local
c#endif /* MPI_HYDRO */

                  k = o_s + width * kk
                  km1 = k - 1
                  IF (km1 .lt. o_s) km1 = o_e
                  
                  ql (:,:,1) = rho   (1:m,ji1:jf1,km1)
                  ql (:,:,2) = v_1   (1:m,ji1:jf1,km1)
                  ql (:,:,3) = v_2   (1:m,ji1:jf1,km1)
                  ql (:,:,4) = v_3   (1:m,ji1:jf1,km1)
                  ql (:,:,5) = eps   (1:m,ji1:jf1,km1)
                  ql (:,:,6) = p     (1:m,ji1:jf1,km1)
                  ql (:,:,7) = gamm  (1:m,ji1:jf1,km1)
#ifdef CFC_MHD
                  ql (:,:,8)  = b_1 (1:m,ji1:jf1,km1)
                  ql (:,:,9)  = b_2 (1:m,ji1:jf1,km1)
                  ql (:,:,10) = b_3 (1:m,ji1:jf1,km1)
                  ql (:,:,11) = psi (1:m,ji1:jf1,km1)
#endif /* CFC_MHD */
                  qlxnu (:,:,1:config%qn) = 
     &                 xnnu (1:m,ji1:jf1,km1,1:config%qn)

c#ifdef MPI_HYDRO 
               END IF

               IF (kk. EQ. nn - 1 .AND.  at_boundary) THEN !last local super-cell

c     send to left neighbour, receive from right neighbour
c     (note of caution: sending data from super-cell 1, i.e. no the one
c     we are currently working on!)
                  q (:,:,1) = rho (1:m,ji1:jf1,o_s)
                  q (:,:,2) = v_1 (1:m,ji1:jf1,o_s)
                  q (:,:,3) = v_2 (1:m,ji1:jf1,o_s)
                  q (:,:,4) = v_3 (1:m,ji1:jf1,o_s)
                  q (:,:,5) = eps (1:m,ji1:jf1,o_s)
                  q (:,:,6) = p   (1:m,ji1:jf1,o_s)
                  q (:,:,7) = gamm(1:m,ji1:jf1,o_s)
#ifdef CFC_MHD
                  q (:,:,8)  = b_1 (1:m,ji1:jf1,o_s)
                  q (:,:,9)  = b_2 (1:m,ji1:jf1,o_s)
                  q (:,:,10) = b_3 (1:m,ji1:jf1,o_s)
                  q (:,:,11) = psi (1:m,ji1:jf1,o_s)
#endif /* CFC_MHD */
                  qxnu (:,:,1:config%qn) = 
     &                 xnnu (1:m,ji1:jf1,o_s,1:config%qn)
                  CALL MPI_Cart_shift (cart_comm, 1_ik, -nd,
     &                 src, dest, ierr)
                  IF (MOD (src / nd, 2_ik) .EQ. 1_ik) THEN
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5002_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5003_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Recv (qr, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5002_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qrxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5003_ik,
     &                    cart_comm, mpistat, ierr)
                  ELSE
                     CALL MPI_Recv (qr, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5002_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qrxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5003_ik,
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5002_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5003_ik,
     &                    cart_comm, ierr)
                  END IF

                  IF (need_intracell (region)) THEN
                     CALL MPI_Bcast (qr, m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qrxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                  END IF


               ELSE IF (kk .EQ. nn - 1 .AND. .NOT. at_boundary) THEN

                  IF (need_intracell (region)) THEN
                     CALL MPI_Bcast (qr,    m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qrxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                  END IF

               ELSE             !data is process-local
c#endif /* MPI_HYDRO */

                  k = o_s + width * (kk + 1) - 1
                  kp1 = k + 1
                  IF (kp1 .gt. o_e) kp1 = o_s
                  
                  qr    (:,:,1) = rho (1:m,ji1:jf1,kp1) 
                  qr    (:,:,2) = v_1 (1:m,ji1:jf1,kp1)
                  qr    (:,:,3) = v_2 (1:m,ji1:jf1,kp1)
                  qr    (:,:,4) = v_3 (1:m,ji1:jf1,kp1)
                  qr    (:,:,5) = eps (1:m,ji1:jf1,kp1)
                  qr    (:,:,6) = p   (1:m,ji1:jf1,kp1)
                  qr    (:,:,7) = gamm(1:m,ji1:jf1,kp1)
#ifdef CFC_MHD
                  qr    (:,:,8)  = b_1 (1:m,ji1:jf1,kp1)
                  qr    (:,:,9)  = b_2 (1:m,ji1:jf1,kp1)
                  qr    (:,:,10) = b_3 (1:m,ji1:jf1,kp1)
                  qr    (:,:,11) = psi (1:m,ji1:jf1,kp1)
#endif /* CFC_MHD */                  
                  qrxnu (:,:,1:config%qn) = 
     &                 xnnu    (1:m,ji1:jf1,kp1,1:config%qn)

c#ifdef MPI_HYDRO
               END IF
c#endif /* MPI_HYDRO */

c     Copy conserved variables and mass fractions into temporary arrays
c     (simplifies reconstruction loop)
               q (:,:,1) = rho (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,2) = v_1 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,3) = v_2 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,4) = v_3 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,5) = eps (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,6) = p   (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,7) = gamm(1:m,ji1:jf1,o_s+width*kk)
#ifdef CFC_MHD               
               q (:,:,8)  = b_1 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,9)  = b_2 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,10) = b_3 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,11) = psi (1:m,ji1:jf1,o_s+width*kk)
#endif /* CFC_MHD */               
               qxnu (:,:,1:config%qn) =
     &              xnnu (1:m,ji1:jf1,o_s+width*kk,1:config%qn)


#ifdef THIRD_ORDER_PROLONGATION
c Following Engwirda & Kelley (2016)
               DO j = ji1, jf1

                  DO l = 1, nvars   
                     DO i = 1, m
                        dl =
     &                      (5.0_rk * q  (i, j, l) +
     &                       2.0_rk * ql (i, j, l) - 
     &                                qr (i, j, l)) / 6.0_rk
                        dr =
     &                      (5.0_rk * q  (i, j, l) -
     &                                ql (i, j, l) + 
     &                       2.0_rk * qr (i, j, l)) / 6.0_rk
                        slope = dr - dl
                        
                        IF ((q  (i, j, l) - ql (i, j, l)) *
     &                      (qr (i, j, l) - q  (i, j, l))
     &                       .LE. 0.0_rk) THEN
c     flattening at local extremum
                           dl = q (i, j, l)
                           dr = q (i, j, l)
                           slope = 0.0_rk
                        ELSE
                           slope = SIGN (
     &                          MIN (ABS (slope), 2.0_rk * MIN (
     &                          ABS (q  (i, j, l) - ql (i, j, l)),
     &                          ABS (qr (i, j, l) - q  (i, j, l)))),
     &                          slope)
c     ensure interface values are bounded
                           IF ((q  (i, j, l) - dl) * (dl - ql (i, j, l))
     &                          .LT. 0.0_rk)
     &                          dl = q (i, j, l) - 0.5_rk * slope
                           IF ((qr (i, j, l) - dr) * (dr - q  (i, j, l))
     &                          .LT. 0.0_rk)
     &                          dr = q (i, j, l) + 0.5_rk * slope
                           IF ((dr - dl) * (dl - 3.0_rk * q (i, j, l)
     &                          + 2.0_rk * dr) .LT. 0.0_rk) THEN
                              dl = 3.0_rk * q (i, j, l) - 2.0_rk * dr
                           END IF
                           IF ((dr - dl) * (- dr + 3.0_rk * q (i, j, l)
     &                          - 2.0_rk * dl) .LT. 0.0_rk) then
                              dr = 3.0_rk * q (i, j, l) - 2.0_rk * dl
                           END IF
                        END IF
                        ql (i, j, l) = dl
                        qr (i, j, l) = dr
c     qr and ql now contain interpolated values at super-cell interfaces                        
                     END DO
                  END DO

                  DO l = 1, config%qn !mass fractions and Y_e
                     DO i = 1, m
                        dl =
     &                      (5.0_rk * qxnu  (i, j, l) +
     &                       2.0_rk * qlxnu (i, j, l) - 
     &                                qrxnu (i, j, l)) / 6.0_rk
                        dr =
     &                      (5.0_rk * qxnu  (i, j, l) -
     &                                qlxnu (i, j, l) + 
     &                       2.0_rk * qrxnu (i, j, l)) / 6.0_rk
                        slope = dr - dl
                        
                        IF ((qxnu  (i, j, l) - qlxnu (i, j, l)) *
     &                      (qrxnu (i, j, l) - qxnu  (i, j, l))
     &                       .LE. 0.0_rk) THEN
c     flattening at local extremum
                           dl = qxnu (i, j, l)
                           dr = qxnu (i, j, l)
                           slope = 0.0_rk
                        ELSE
                           slope = SIGN (
     &                          MIN (ABS (slope), 2.0_rk * MIN (
     &                          ABS (qxnu (i,j,l) - qlxnu(i,j,l)),
     &                          ABS (qrxnu(i,j,l) - qxnu (i,j,l)))),
     &                          slope)
c     ensure interface values are bounded
                           IF ((qxnu (i, j, l) - dl) *
     &                          (dl - qlxnu (i, j, l)) .LT. 0.0_rk)
     &                          dl = qxnu (i, j, l) - 0.5_rk * slope
                           IF ((qrxnu (i, j, l) - dr) *
     &                          (dr - qxnu (i, j, l)) .LT. 0.0_rk)
     &                          dr = qxnu (i, j, l) + 0.5_rk * slope
                           IF ((dr - dl) *
     &                          (dl - 3.0_rk * qxnu (i, j, l) +
     &                          2.0_rk * dr) .LT. 0.0_rk) THEN
                              dl = 3.0_rk * qxnu (i, j, l) - 2.0_rk * dr
                           END IF
                           IF ((dr - dl) *
     &                          (- dr + 3.0_rk * qxnu (i, j, l) -
     &                          2.0_rk * dl) .LT. 0.0_rk) then
                              dr = 3.0_rk * qxnu (i, j, l) - 2.0_rk * dl
                           END IF
                        END IF
                        qlxnu (i, j, l) = dl
                        qrxnu (i, j, l) = dr
c     qrxnu and qlxnu now contain interpolated values at super-cell interfaces                        
                     END DO
                  END DO                  


c     ------------------------------------------------------------------
c     Mass fractions: Ensure that edge values sum to unity                 
c     left interface
                  splus  (1:m) = 0.0_rk
                  sminus (1:m) = 0.0_rk
                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        IF (qlxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
     &                       .GT. 0.0_rk) THEN
                           splus (i) = splus (i) +
     &                          qlxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        ELSE
                           sminus (i) = sminus (i) +
     &                          qlxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        END IF
                     END DO
                  END DO
                  
                  ssum (:) = splus (:) + sminus (:)
                  WHERE (ssum .GT. 1e-15_rk)
                     cflat = - sminus / splus
                  ELSE WHERE (ssum .LT. -1e-15_rk)
                     cflat = - splus  / sminus
                  ELSE WHERE
                     cflat = 1.0_rk
                  END WHERE

                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        slope = qlxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        IF (ssum (i) .GT. 0.0_rk) THEN
                           qlxnu (i,j,i_nuc) = qxnu (i,j,i_nuc) +
     &                          MAX (cflat (i) * slope, 0.0_rk) +
     &                          MIN (slope, 0.0_rk)
                        ELSE
                           qlxnu (i,j,i_nuc) = qxnu (i,j,i_nuc) +
     &                          MIN (cflat (i) * slope, 0.0_rk) +
     &                          MAX (slope, 0.0_rk)
                        END IF
                     END DO
                  END DO

c     right interface
                  splus  (1:m) = 0.0_rk
                  sminus (1:m) = 0.0_rk
                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        IF (qrxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
     &                       .GT. 0.0_rk) THEN
                           splus (i) = splus (i) +
     &                          qrxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        ELSE
                           sminus (i) = sminus (i) +
     &                          qrxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        END IF
                     END DO
                  END DO
                  
                  ssum (:) = splus (:) + sminus (:)
                  WHERE (ssum .GT. 0.0_rk)
                     cflat = - sminus / splus
                  ELSE WHERE (ssum .LT. 0.0_rk)
                     cflat = - splus  / sminus
                  ELSE WHERE
                     cflat = 1.0_rk
                  END WHERE

                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        slope = qrxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        IF (ssum (i) .GT. 0.0_rk) THEN
                           qrxnu (i,j,i_nuc) = qxnu (i,j,i_nuc) +
     &                          MAX (cflat (i) * slope, 0.0_rk) +
     &                          MIN (slope, 0.0_rk)
                        ELSE
                           qrxnu (i,j,i_nuc) = qxnu (i,j,i_nuc) +
     &                          MIN (cflat (i) * slope, 0.0_rk) +
     &                          MAX (slope, 0.0_rk)
                        END IF
                     END DO
                  END DO
                  
c     ------------------------------------------------------------------
c     reconstruct values in sub-cells
                  DO k = o_s + width * kk, o_s + width * (kk + 1) - 1_ik
                     
                     l = MOD (k - 1, width0)

                     coeff = REAL (1_ik + 2_ik * l - width0, KIND=rk)
                     coeff2 =
     &                    REAL (1_ik + 3_ik * coeff ** 2, KIND=rk) /
     &                    REAL (12_ik * width0 ** 2)
                     coeff = coeff /
     &                    REAL (2_ik * width0, KIND=rk)

                     DO i = 1, m
                        a0(:) = 0.25_rk *
     &                       (6.0_rk*q(i,j,:)-ql(i,j,:)-qr(i,j,:))
                        a1(:) = qr(i,j,:) - ql(i,j,:)
                        a2(:) = 3.0_rk *
     &                       (ql(i,j,:)+qr(i,j,:)-2.0_rk*q(i,j,:))

                        rho (i, j, k) = a0 (1) +
     &                       coeff * a1 (1) + coeff2 * a2 (1)
                        v_1 (i, j, k) = a0 (2) +
     &                       coeff * a1 (2) + coeff2 * a2 (2)
                        v_2 (i, j, k) = a0 (3) +
     &                       coeff * a1 (3) + coeff2 * a2 (3)
                        v_3 (i, j, k) = a0 (4) +
     &                       coeff * a1 (4) + coeff2 * a2 (4)
                        eps (i, j, k) = a0 (5) +
     &                       coeff * a1 (5) + coeff2 * a2 (5)
                        p   (i, j, k) = a0 (6) +
     &                       coeff * a1 (6) + coeff2 * a2 (6)
                        gamm(i, j, k) = a0 (7) +
     &                       coeff * a1 (7) + coeff2 * a2 (7)
#ifdef CFC_MHD
                        b_1 (i, j, k) = a0 (8) +
     &                       coeff * a1 (8) + coeff2 * a2 (8)
                        b_2 (i, j, k) = a0 (9) +
     &                       coeff * a1 (9) + coeff2 * a2 (9)
                        b_3 (i, j, k) = a0 (10) +
     &                       coeff * a1 (10) + coeff2 * a2 (10)
                        psi (i, j, k) = a0 (11) +
     &                       coeff * a1 (11) + coeff2 * a2 (11)
#endif /* CFC_MHD */                        
                     END DO

                     DO i_nuc = 1, config%qn
                        DO i = 1, n
                           a0 (1) = 0.25_rk *
     &                          (6.0_rk * qxnu(i, j, i_nuc)
     &                          - qlxnu (i, j, i_nuc)
     &                          - qrxnu (i, j, i_nuc))
                           a1 (1) = qrxnu (i, j, i_nuc)
     &                          - qlxnu (i, j, i_nuc)
                           a2 (1) = 3.0_rk *
     &                          (qlxnu (i, j, i_nuc)
     &                          + qrxnu (i, j, i_nuc)-
     &                          2.0_rk * qxnu(i, j, i_nuc))
                           xnnu (i, j, k, i_nuc) = a0 (1) +
     &                          coeff * a1 (1) + coeff2 * a2 (1)
                        END DO
                     END DO
                     
                  END DO
               END DO
                           
#else /* THIRD_ORDER_PROLONGATION */
c     piecewise linear reconstruction
               
c     ------------------------------------------------------------------
c     Determine slope in super-cell and reconstruct sub-cells (q and
c     qxnu are overwritten in the process and contain the slopes in the end)

               DO j = ji1, jf1
                  !print *,'slopes',j,ji1,jf1
                  DO l = 1, nvars   !conserved variables
                     DO i = 1, m
                        dl = q  (i, j, l) - ql (i, j, l)
                        dr = qr (i, j, l) - q  (i, j, l)
                        ql (i, j, l) = dl
                        qr (i, j, l) = dr
#ifdef MC_LIMITER
c     Monotonized central limiter
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = MIN (ABS (0.5_rk * (dl + dr)),
     &                          2.0_rk * MIN (ABS (dl), ABS (dr)))
                           slope  = SIGN (slope, dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef VAN_LEER_LIMITER
c     van Leer's smooth slope limiter 
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = 2.0 * dl * dr / (dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef SUPERBEE_LIMITER
                        slope =
     &                       (sign (1.0_rk, dl) + sign (1.0_rk, dr)) *
     &                       min(abs (dl), abs (dr),
     &                       0.5_rk * max (abs (dl), abs (dr)))
#endif
                        q (i, j, l) = slope !* scal

                     END DO
                  END DO !l-loop

                  DO i_nuc = 1, config%qn !mass fractions and Y_e
                     DO i = 1, m

                        dl = qxnu  (i, j, i_nuc) - qlxnu (i, j, i_nuc)
                        dr = qrxnu (i, j, i_nuc) - qxnu  (i, j, i_nuc)
#ifdef MC_LIMITER
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = MIN (ABS (0.5_rk * (dl + dr)),
     &                          2.0_rk * MIN (ABS (dl), ABS (dr)))
                           slope  = SIGN (slope, dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef VAN_LEER_LIMITER
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = 2.0 * dl * dr / (dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef SUPERBEE_LIMITER
                        slope =
     &                       (sign (1.0_rk, dl) + sign (1.0_rk, dr)) *
     &                       min(abs (dl), abs (dr),
     &                       0.5_rk * max (abs (dl), abs (dr)))
#endif
                        qxnu (i, j, i_nuc) = slope 
                        
                     END DO
                  END DO !species-loop

c     ------------------------------------------------------------------
c     Linear reconstruction of conserved variables

                  DO k = o_s + width * kk, o_s + width * (kk + 1) - 1_ik

                     l = MOD (k - 1, width0)

c     The integral of the pre-reconstructed variable U over sub-cells
c     can be expressed analytically as <U>_supercell + coeff * dU/dvarphi:
                     coeff = (1_ik + 2_ik * l - width0) /
     &                    REAL (2_ik * width0, KIND=rk)

                     DO i = 1, m

                        rho (i, j, k) = rho (i, j, k) +
     &                       coeff * q (i, j, 1)
                        v_1 (i, j, k) = v_1   (i, j, k) +
     &                       coeff * q (i, j, 2)
                        v_2 (i, j, k) = v_2   (i, j, k) +
     &                       coeff * q (i, j, 3)
                        v_3 (i, j, k) = v_3   (i, j, k) +
     &                       coeff * q (i, j, 4)
                        eps (i, j, k) = eps   (i, j, k) +
     &                       coeff * q (i, j, 5)
                        p   (i, j, k) = p     (i, j, k) +
     &                       coeff * q (i, j, 6)
                        gamm(i, j, k) = gamm  (i, j, k) +
     &                       coeff * q (i, j, 7)
#ifdef CFC_MHD
                        b_1 (i, j, k) = b_1 (i, j, k) +
     &                       coeff * q (i, j, 8)
                        b_2 (i, j, k) = b_2 (i, j, k) +
     &                       coeff * q (i, j, 9)
                        b_3 (i, j, k) = b_3 (i, j, k) +
     &                       coeff * q (i, j, 10)
                        psi (i, j, k) = psi (i, j, k) +
     &                       coeff * q (i, j, 11)
#endif /* CFC_MHD */                        
                     END DO

                  END DO        !loop over sub-cells

c     ------------------------------------------------------------------
c     Mass fractions:
c     Step 1: Ensure that  \sum_i dX_i/dvarphi = 1 by rescaling
                  splus  (1:m) = 0.0_rk
                  sminus (1:m) = 0.0_rk
                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        IF (qxnu (i, j, i_nuc) .GT. 0.0_rk) THEN
                           splus (i) = splus (i) +
     &                          qxnu (i, j, i_nuc)
                        ELSE
                           sminus (i) = sminus (i) +
     &                          qxnu (i, j, i_nuc)
                        END IF
                     END DO
                  END DO
                  
                  ssum (:) = splus (:) + sminus (:)
                  WHERE (ssum .GT. 0.0_rk)
                     cflat = - sminus / splus
                  ELSE WHERE (ssum .LT. 0.0_rk)
                     cflat = - splus  / sminus
                  ELSE WHERE
                     cflat = 1.0_rk
                  END WHERE

                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        IF (ssum (i) .GT. 0.0_rk) THEN
                           qxnu (i, j, i_nuc) = MAX (cflat (i) *
     &                          qxnu (i, j, i_nuc), 0.0_rk) +
     &                          MIN (qxnu (i, j, i_nuc), 0.0_rk) 
                        ELSE
                           qxnu (i, j, i_nuc) = MIN (cflat (i) *
     &                          qxnu (i, j, i_nuc), 0.0_rk) +
     &                          MAX (qxnu (i, j, i_nuc), 0.0_rk) 
                        END IF
                     END DO
                  END DO

c     Step 2: Linear reconstruction of mass fractions
                  DO k = o_s + width * kk, o_s + width * (kk + 1) - 1_ik

                     l = MOD (k - 1, width0)

                     coeff = (1_ik + 2_ik * l - width0) /
     &                    REAL (2_ik * width0, KIND=rk)

                     DO i_nuc = 1, config%qn
                        DO i = 1, m
                           xnnu (i, j, k, i_nuc) =
     &                          xnnu  (i, j, k, i_nuc) +
     &                          coeff * qxnu (i, j, i_nuc)
                        END DO
                     END DO

                  END DO        !loop over sub-cells (varphi direction)
               END DO           !loop over theta-index
#endif /* THIRD_ORDER_PROLONGATION */


            END DO              !loop over super-cells

            
c            CALL MPI_Barrier (ierr)


            DEALLOCATE (q, qxnu, ql, qr, qlxnu, qrxnu)

         END IF

      END DO !loop over regions

      
      RETURN

      END SUBROUTINE prolongation

c     ==================================================================
      


c     ==================================================================

      SUBROUTINE prolongation_output

c     ==================================================================
c     - Reconstructs gradients of primitive variables
c       to recover higher-order accuracy.
c     ------------------------------------------------------------------

      USE precision
      USE conserved_cfc, ONLY: d_cap_hat, tau_hat
      USE hydro_primitives_cfc
      USE mesh_cfc
      USE mo_mpi

      USE configure

      USE, INTRINSIC :: iso_c_binding

C      INCLUDE 'fftw3-mpi.f'

      IMPLICIT NONE

#ifdef CFC_MHD
      INTEGER (KIND=ik), PARAMETER :: nvars = 11 + 9
#else
      INTEGER (KIND=ik), PARAMETER :: nvars = 7 + 9
#endif      
      REAL (KIND=rk), ALLOCATABLE :: q (:,:,:), qxnu (:,:,:)
      REAL (KIND=rk), ALLOCATABLE :: qr (:,:,:), qrxnu (:,:,:)
      REAL (KIND=rk), ALLOCATABLE :: ql (:,:,:), qlxnu (:,:,:)
      REAL (KIND=rk) :: splus (1:m), sminus (1:m),
     &     ssum (1:m), cflat (1:m)
      REAL (KIND=rk) :: 
     &     a0 (-8:nvars-9), a1 (-8:nvars-9), a2 (-8:nvars-9)
      REAL (KIND=rk) :: dl, dr, slope, scal, coeff, coeff2
      
      INTEGER (KIND=ik) :: i, j, k, jk, kk, i_nuc, ji1, jf1, region,
     &     istat, level, ki, kf, km1, kp1, l, nn, nd, width, width0

c#ifdef MPI_HYDRO
      INTEGER (KIND=ik) :: src, dest, ierr, mpistat (MPI_STATUS_SIZE)
      LOGICAL :: at_boundary
c#endif

      DO region = region_s, region_e
         IF (clevel_region (region) .NE. 0) THEN
            
            ji1 = j_i (region)
            jf1 = j_f (region)
            
            ALLOCATE (q (1:m, ji1:jf1, -8:nvars-9), stat = istat)
            ALLOCATE (qxnu (1:m, ji1:jf1, 1:config%qn), stat = istat)
            ALLOCATE (ql (1:m, ji1:jf1, -8:nvars-9), stat = istat)
            ALLOCATE (qlxnu (1:m, ji1:jf1, 1:config%qn), stat = istat)
            ALLOCATE (qr (1:m, ji1:jf1, -8:nvars-9), stat = istat)
            ALLOCATE (qrxnu (1:m, ji1:jf1, 1:config%qn), stat = istat)
            
c     Determine number of super-cells on local domain
            width0 = 2 ** clevel_region (region)
            width  = MIN (width0, qz_proc)
            scal = 1.0_rk / REAL (width0, kind = rk)
            nn = o_loc / width 

c     Number of domains sharing the same super-cell
            nd = MAX (1_ik, width0 / o_loc) 

c     loop over super-cells (pre-reconstruction)

            DO kk = 0, nn - 1

c#ifdef MPI_HYDRO
               at_boundary = 
     &              ((o_s - 1) / 2 ** clevel_region (region) .NE.
     &              (o_s - 2) / 2 ** clevel_region (region)) .OR.
     &              o_s .EQ. 1

               IF (kk. EQ. 0_ik .AND. at_boundary) THEN !first local super-cell
c     send to right neighbour, receive from left neighbour
c     (note of caution: sending data from super-cell nn, i.e. not the one
c     we are currently working on!)                  
                  q (:,:, 0) = entropy(1:m,ji1:jf1,o_e)
                  q (:,:,-1) = cpot   (1:m,ji1:jf1,o_e,1)
                  q (:,:,-2) = cpot   (1:m,ji1:jf1,o_e,2)
                  q (:,:,-3) = cpot   (1:m,ji1:jf1,o_e,3)
                  q (:,:,-4) = cpot   (1:m,ji1:jf1,o_e,4)
                  q (:,:,-5) = t      (1:m,ji1:jf1,o_e)
                  q (:,:,-6) = w      (1:m,ji1:jf1,o_e)
                  q (:,:,-7) = d_cap_hat (1:m,ji1:jf1,o_e)
                  q (:,:,-8) = tau_hat   (1:m,ji1:jf1,o_e)

                  q (:,:,1) = rho   (1:m,ji1:jf1,o_e)
                  q (:,:,2) = v_1   (1:m,ji1:jf1,o_e)
                  q (:,:,3) = v_2   (1:m,ji1:jf1,o_e)
                  q (:,:,4) = v_3   (1:m,ji1:jf1,o_e)
                  q (:,:,5) = eps   (1:m,ji1:jf1,o_e)
                  q (:,:,6) = p     (1:m,ji1:jf1,o_e)
                  q (:,:,7) = gamm  (1:m,ji1:jf1,o_e)
#ifdef CFC_MHD
                  q (:,:,8)  = b_1 (1:m,ji1:jf1,o_e)
                  q (:,:,9)  = b_2 (1:m,ji1:jf1,o_e)
                  q (:,:,10) = b_3 (1:m,ji1:jf1,o_e)
                  q (:,:,11) = psi (1:m,ji1:jf1,o_e)
#endif /* CFC_MHD */
                  qxnu (:,:,1:config%qn) = 
     &                 xnnu (1:m,ji1:jf1,o_e,1:config%qn)
                  CALL MPI_Cart_shift (cart_comm, 1_ik, nd,
     &                 src, dest, ierr)
                  IF (MOD (src / nd, 2_ik) .EQ. 1_ik) THEN
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5000_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5001_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Recv (ql, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5000_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qlxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5001_ik,
     &                    cart_comm, mpistat, ierr)
                  ELSE
                     CALL MPI_Recv (ql, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5000_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qlxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5001_ik,
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5000_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5001_ik,
     &                    cart_comm, ierr)
                  END IF

                  IF (need_intracell (region)) THEN
                     CALL MPI_Bcast (ql,    m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region),
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qlxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                  END IF

               ELSE IF (kk .EQ. 0_ik .AND. .NOT. at_boundary) THEN

                  IF (need_intracell (region)) THEN
                     CALL MPI_Bcast (ql,    m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region),
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qlxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region),
     &                    comm_super_cell (region), ierr)
                  END IF

               ELSE             !data is process-local
c#endif /* MPI_HYDRO */

                  k = o_s + width * kk
                  km1 = k - 1
                  IF (km1 .lt. o_s) km1 = o_e
                  
                  ql (:,:, 0) = entropy(1:m,ji1:jf1,km1)
                  ql (:,:,-1) = cpot   (1:m,ji1:jf1,km1,1)
                  ql (:,:,-2) = cpot   (1:m,ji1:jf1,km1,2)
                  ql (:,:,-3) = cpot   (1:m,ji1:jf1,km1,3)
                  ql (:,:,-4) = cpot   (1:m,ji1:jf1,km1,4)
                  ql (:,:,-5) = t      (1:m,ji1:jf1,km1)
                  ql (:,:,-6) = w      (1:m,ji1:jf1,km1)
                  ql (:,:,-7) = d_cap_hat (1:m,ji1:jf1,km1)
                  ql (:,:,-8) = tau_hat   (1:m,ji1:jf1,km1)

                  ql (:,:,1) = rho   (1:m,ji1:jf1,km1)
                  ql (:,:,2) = v_1   (1:m,ji1:jf1,km1)
                  ql (:,:,3) = v_2   (1:m,ji1:jf1,km1)
                  ql (:,:,4) = v_3   (1:m,ji1:jf1,km1)
                  ql (:,:,5) = eps   (1:m,ji1:jf1,km1)
                  ql (:,:,6) = p     (1:m,ji1:jf1,km1)
                  ql (:,:,7) = gamm  (1:m,ji1:jf1,km1)
#ifdef CFC_MHD
                  ql (:,:,8)  = b_1 (1:m,ji1:jf1,km1)
                  ql (:,:,9)  = b_2 (1:m,ji1:jf1,km1)
                  ql (:,:,10) = b_3 (1:m,ji1:jf1,km1)
                  ql (:,:,11) = psi (1:m,ji1:jf1,km1)
#endif /* CFC_MHD */
                  qlxnu (:,:,1:config%qn) = 
     &                 xnnu (1:m,ji1:jf1,km1,1:config%qn)

c#ifdef MPI_HYDRO 
               END IF

               IF (kk. EQ. nn - 1 .AND.  at_boundary) THEN !last local super-cell

c     send to left neighbour, receive from right neighbour
c     (note of caution: sending data from super-cell 1, i.e. no the one
c     we are currently working on!)
                  q (:,:, 0) = entropy(1:m,ji1:jf1,o_s)
                  q (:,:,-1) = cpot   (1:m,ji1:jf1,o_s,1)
                  q (:,:,-2) = cpot   (1:m,ji1:jf1,o_s,2)
                  q (:,:,-3) = cpot   (1:m,ji1:jf1,o_s,3)
                  q (:,:,-4) = cpot   (1:m,ji1:jf1,o_s,4)
                  q (:,:,-5) = t      (1:m,ji1:jf1,o_s)
                  q (:,:,-6) = w      (1:m,ji1:jf1,o_s)
                  q (:,:,-7) = d_cap_hat (1:m,ji1:jf1,o_s)
                  q (:,:,-8) = tau_hat   (1:m,ji1:jf1,o_s)

                  q (:,:,1) = rho (1:m,ji1:jf1,o_s)
                  q (:,:,2) = v_1 (1:m,ji1:jf1,o_s)
                  q (:,:,3) = v_2 (1:m,ji1:jf1,o_s)
                  q (:,:,4) = v_3 (1:m,ji1:jf1,o_s)
                  q (:,:,5) = eps (1:m,ji1:jf1,o_s)
                  q (:,:,6) = p   (1:m,ji1:jf1,o_s)
                  q (:,:,7) = gamm(1:m,ji1:jf1,o_s)
#ifdef CFC_MHD
                  q (:,:,8)  = b_1 (1:m,ji1:jf1,o_s)
                  q (:,:,9)  = b_2 (1:m,ji1:jf1,o_s)
                  q (:,:,10) = b_3 (1:m,ji1:jf1,o_s)
                  q (:,:,11) = psi (1:m,ji1:jf1,o_s)
#endif /* CFC_MHD */
                  qxnu (:,:,1:config%qn) = 
     &                 xnnu (1:m,ji1:jf1,o_s,1:config%qn)
                  CALL MPI_Cart_shift (cart_comm, 1_ik, -nd,
     &                 src, dest, ierr)
                  IF (MOD (src / nd, 2_ik) .EQ. 1_ik) THEN
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5002_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5003_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Recv (qr, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5002_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qrxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5003_ik,
     &                    cart_comm, mpistat, ierr)
                  ELSE
                     CALL MPI_Recv (qr, m*(jf1-ji1+1)*nvars, 
     &                    MPI_DOUBLE_PRECISION, src, 5002_ik, 
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Recv (qrxnu, m*(jf1-ji1+1_ik)*config%qn, 
     &                    MPI_DOUBLE_PRECISION, src, 5003_ik,
     &                    cart_comm, mpistat, ierr)
                     CALL MPI_Send (q, m*(jf1-ji1+1)*nvars,
     &                    MPI_DOUBLE_PRECISION, dest, 5002_ik,
     &                    cart_comm, ierr)
                     CALL MPI_Send (qxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, dest, 5003_ik,
     &                    cart_comm, ierr)
                  END IF

                  IF (need_intracell (region)) THEN
                     CALL MPI_Bcast (qr, m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qrxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                  END IF


               ELSE IF (kk .EQ. nn - 1 .AND. .NOT. at_boundary) THEN

                  IF (need_intracell (region)) THEN
                     CALL MPI_Bcast (qr,    m*(jf1-ji1+1_ik)*nvars,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                     CALL MPI_Bcast (qrxnu, m*(jf1-ji1+1_ik)*config%qn,
     &                    MPI_DOUBLE_PRECISION, left_root (region), 
     &                    comm_super_cell (region), ierr)
                  END IF

               ELSE             !data is process-local
c#endif /* MPI_HYDRO */

                  k = o_s + width * (kk + 1) - 1
                  kp1 = k + 1
                  IF (kp1 .gt. o_e) kp1 = o_s
                  
                  qr    (:,:, 0) = entropy(1:m,ji1:jf1,kp1) 
                  qr    (:,:,-1) = cpot   (1:m,ji1:jf1,kp1,1) 
                  qr    (:,:,-2) = cpot   (1:m,ji1:jf1,kp1,2) 
                  qr    (:,:,-3) = cpot   (1:m,ji1:jf1,kp1,3) 
                  qr    (:,:,-4) = cpot   (1:m,ji1:jf1,kp1,4) 
                  qr    (:,:,-5) = t      (1:m,ji1:jf1,kp1) 
                  qr    (:,:,-6) = w      (1:m,ji1:jf1,kp1) 
                  qr    (:,:,-7) = d_cap_hat (1:m,ji1:jf1,kp1) 
                  qr    (:,:,-8) = tau_hat   (1:m,ji1:jf1,kp1) 

                  qr    (:,:,1) = rho (1:m,ji1:jf1,kp1) 
                  qr    (:,:,2) = v_1 (1:m,ji1:jf1,kp1)
                  qr    (:,:,3) = v_2 (1:m,ji1:jf1,kp1)
                  qr    (:,:,4) = v_3 (1:m,ji1:jf1,kp1)
                  qr    (:,:,5) = eps (1:m,ji1:jf1,kp1)
                  qr    (:,:,6) = p   (1:m,ji1:jf1,kp1)
                  qr    (:,:,7) = gamm(1:m,ji1:jf1,kp1)
#ifdef CFC_MHD
                  qr    (:,:,8)  = b_1 (1:m,ji1:jf1,kp1)
                  qr    (:,:,9)  = b_2 (1:m,ji1:jf1,kp1)
                  qr    (:,:,10) = b_3 (1:m,ji1:jf1,kp1)
                  qr    (:,:,11) = psi (1:m,ji1:jf1,kp1)
#endif /* CFC_MHD */                  
                  qrxnu (:,:,1:config%qn) = 
     &                 xnnu    (1:m,ji1:jf1,kp1,1:config%qn)

c#ifdef MPI_HYDRO
               END IF
c#endif /* MPI_HYDRO */

c     Copy conserved variables and mass fractions into temporary arrays
c     (simplifies reconstruction loop)
               q (:,:, 0) = entropy(1:m,ji1:jf1,o_s+width*kk)
               q (:,:,-1) = cpot   (1:m,ji1:jf1,o_s+width*kk,1)
               q (:,:,-2) = cpot   (1:m,ji1:jf1,o_s+width*kk,2)
               q (:,:,-3) = cpot   (1:m,ji1:jf1,o_s+width*kk,3)
               q (:,:,-4) = cpot   (1:m,ji1:jf1,o_s+width*kk,4)
               q (:,:,-5) = t      (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,-6) = w      (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,-7) = d_cap_hat (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,-8) = tau_hat   (1:m,ji1:jf1,o_s+width*kk)

               q (:,:,1) = rho (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,2) = v_1 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,3) = v_2 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,4) = v_3 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,5) = eps (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,6) = p   (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,7) = gamm(1:m,ji1:jf1,o_s+width*kk)
#ifdef CFC_MHD               
               q (:,:,8)  = b_1 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,9)  = b_2 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,10) = b_3 (1:m,ji1:jf1,o_s+width*kk)
               q (:,:,11) = psi (1:m,ji1:jf1,o_s+width*kk)
#endif /* CFC_MHD */               
               qxnu (:,:,1:config%qn) =
     &              xnnu (1:m,ji1:jf1,o_s+width*kk,1:config%qn)


#ifdef THIRD_ORDER_PROLONGATION
c Following Engwirda & Kelley (2016)
               DO j = ji1, jf1

                  DO l = - 8, nvars - 9
                     DO i = 1, m
                        dl =
     &                      (5.0_rk * q  (i, j, l) +
     &                       2.0_rk * ql (i, j, l) - 
     &                                qr (i, j, l)) / 6.0_rk
                        dr =
     &                      (5.0_rk * q  (i, j, l) -
     &                                ql (i, j, l) + 
     &                       2.0_rk * qr (i, j, l)) / 6.0_rk
                        slope = dr - dl
                        
                        IF ((q  (i, j, l) - ql (i, j, l)) *
     &                      (qr (i, j, l) - q  (i, j, l))
     &                       .LE. 0.0_rk) THEN
c     flattening at local extremum
                           dl = q (i, j, l)
                           dr = q (i, j, l)
                           slope = 0.0_rk
                        ELSE
                           slope = SIGN (
     &                          MIN (ABS (slope), 2.0_rk * MIN (
     &                          ABS (q  (i, j, l) - ql (i, j, l)),
     &                          ABS (qr (i, j, l) - q  (i, j, l)))),
     &                          slope)
c     ensure interface values are bounded
                           IF ((q  (i, j, l) - dl) * (dl - ql (i, j, l))
     &                          .LT. 0.0_rk)
     &                          dl = q (i, j, l) - 0.5_rk * slope
                           IF ((qr (i, j, l) - dr) * (dr - q  (i, j, l))
     &                          .LT. 0.0_rk)
     &                          dr = q (i, j, l) + 0.5_rk * slope
                           IF ((dr - dl) * (dl - 3.0_rk * q (i, j, l)
     &                          + 2.0_rk * dr) .LT. 0.0_rk) THEN
                              dl = 3.0_rk * q (i, j, l) - 2.0_rk * dr
                           END IF
                           IF ((dr - dl) * (- dr + 3.0_rk * q (i, j, l)
     &                          - 2.0_rk * dl) .LT. 0.0_rk) then
                              dr = 3.0_rk * q (i, j, l) - 2.0_rk * dl
                           END IF
                        END IF
                        ql (i, j, l) = dl
                        qr (i, j, l) = dr
c     qr and ql now contain interpolated values at super-cell interfaces                        
                     END DO
                  END DO

                  DO l = 1, config%qn !mass fractions and Y_e
                     DO i = 1, m
                        dl =
     &                      (5.0_rk * qxnu  (i, j, l) +
     &                       2.0_rk * qlxnu (i, j, l) - 
     &                                qrxnu (i, j, l)) / 6.0_rk
                        dr =
     &                      (5.0_rk * qxnu  (i, j, l) -
     &                                qlxnu (i, j, l) + 
     &                       2.0_rk * qrxnu (i, j, l)) / 6.0_rk
                        slope = dr - dl
                        
                        IF ((qxnu  (i, j, l) - qlxnu (i, j, l)) *
     &                      (qrxnu (i, j, l) - qxnu  (i, j, l))
     &                       .LE. 0.0_rk) THEN
c     flattening at local extremum
                           dl = qxnu (i, j, l)
                           dr = qxnu (i, j, l)
                           slope = 0.0_rk
                        ELSE
                           slope = SIGN (
     &                          MIN (ABS (slope), 2.0_rk * MIN (
     &                          ABS (qxnu (i,j,l) - qlxnu(i,j,l)),
     &                          ABS (qrxnu(i,j,l) - qxnu (i,j,l)))),
     &                          slope)
c     ensure interface values are bounded
                           IF ((qxnu (i, j, l) - dl) *
     &                          (dl - qlxnu (i, j, l)) .LT. 0.0_rk)
     &                          dl = qxnu (i, j, l) - 0.5_rk * slope
                           IF ((qrxnu (i, j, l) - dr) *
     &                          (dr - qxnu (i, j, l)) .LT. 0.0_rk)
     &                          dr = qxnu (i, j, l) + 0.5_rk * slope
                           IF ((dr - dl) *
     &                          (dl - 3.0_rk * qxnu (i, j, l) +
     &                          2.0_rk * dr) .LT. 0.0_rk) THEN
                              dl = 3.0_rk * qxnu (i, j, l) - 2.0_rk * dr
                           END IF
                           IF ((dr - dl) *
     &                          (- dr + 3.0_rk * qxnu (i, j, l) -
     &                          2.0_rk * dl) .LT. 0.0_rk) then
                              dr = 3.0_rk * qxnu (i, j, l) - 2.0_rk * dl
                           END IF
                        END IF
                        qlxnu (i, j, l) = dl
                        qrxnu (i, j, l) = dr
c     qrxnu and qlxnu now contain interpolated values at super-cell interfaces                        
                     END DO
                  END DO                  


c     ------------------------------------------------------------------
c     Mass fractions: Ensure that edge values sum to unity                 
c     left interface
                  splus  (1:m) = 0.0_rk
                  sminus (1:m) = 0.0_rk
                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        IF (qlxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
     &                       .GT. 0.0_rk) THEN
                           splus (i) = splus (i) +
     &                          qlxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        ELSE
                           sminus (i) = sminus (i) +
     &                          qlxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        END IF
                     END DO
                  END DO
                  
                  ssum (:) = splus (:) + sminus (:)
                  WHERE (ssum .GT. 1e-15_rk)
                     cflat = - sminus / splus
                  ELSE WHERE (ssum .LT. -1e-15_rk)
                     cflat = - splus  / sminus
                  ELSE WHERE
                     cflat = 1.0_rk
                  END WHERE

                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        slope = qlxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        IF (ssum (i) .GT. 0.0_rk) THEN
                           qlxnu (i,j,i_nuc) = qxnu (i,j,i_nuc) +
     &                          MAX (cflat (i) * slope, 0.0_rk) +
     &                          MIN (slope, 0.0_rk)
                        ELSE
                           qlxnu (i,j,i_nuc) = qxnu (i,j,i_nuc) +
     &                          MIN (cflat (i) * slope, 0.0_rk) +
     &                          MAX (slope, 0.0_rk)
                        END IF
                     END DO
                  END DO

c     right interface
                  splus  (1:m) = 0.0_rk
                  sminus (1:m) = 0.0_rk
                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        IF (qrxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
     &                       .GT. 0.0_rk) THEN
                           splus (i) = splus (i) +
     &                          qrxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        ELSE
                           sminus (i) = sminus (i) +
     &                          qrxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        END IF
                     END DO
                  END DO
                  
                  ssum (:) = splus (:) + sminus (:)
                  WHERE (ssum .GT. 0.0_rk)
                     cflat = - sminus / splus
                  ELSE WHERE (ssum .LT. 0.0_rk)
                     cflat = - splus  / sminus
                  ELSE WHERE
                     cflat = 1.0_rk
                  END WHERE

                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        slope = qrxnu (i, j, i_nuc) - qxnu (i, j, i_nuc)
                        IF (ssum (i) .GT. 0.0_rk) THEN
                           qrxnu (i,j,i_nuc) = qxnu (i,j,i_nuc) +
     &                          MAX (cflat (i) * slope, 0.0_rk) +
     &                          MIN (slope, 0.0_rk)
                        ELSE
                           qrxnu (i,j,i_nuc) = qxnu (i,j,i_nuc) +
     &                          MIN (cflat (i) * slope, 0.0_rk) +
     &                          MAX (slope, 0.0_rk)
                        END IF
                     END DO
                  END DO
                  
c     ------------------------------------------------------------------
c     reconstruct values in sub-cells
                  DO k = o_s + width * kk, o_s + width * (kk + 1) - 1_ik
                     
                     l = MOD (k - 1, width0)

                     coeff = REAL (1_ik + 2_ik * l - width0, KIND=rk) /
     &                    REAL (2_ik * width0, KIND=rk)
                     coeff2 =
     &                    REAL (1_ik+3_ik*(2_ik*l-width0)**2, KIND=rk) /
     &                    REAL (12_ik * width0 ** 2)

                     DO i = 1, m
                        a0(:) = 0.25_rk *
     &                       (6.0_rk*q(i,j,:)-ql(i,j,:)-qr(i,j,:))
                        a1(:) = qr(i,j,:) - ql(i,j,:)
                        a2(:) = 3.0_rk *
     &                       (ql(i,j,:)+qr(i,j,:)-2.0_rk*q(i,j,:))


                        entropy (i, j, k) = a0 (0) +
     &                       coeff * a1 (0) + coeff2 * a2 (0)
                        cpot    (i, j, k, 1) = a0 (-1) +
     &                       coeff * a1 (-1) + coeff2 * a2 (-1)
                        cpot    (i, j, k, 2) = a0 (-2) +
     &                       coeff * a1 (-2) + coeff2 * a2 (-2)
                        cpot    (i, j, k, 3) = a0 (-3) +
     &                       coeff * a1 (-3) + coeff2 * a2 (-3)
                        cpot    (i, j, k, 4) = a0 (-4) +
     &                       coeff * a1 (-4) + coeff2 * a2 (-4)
                        t       (i, j, k)    = a0 (-5) +
     &                       coeff * a1 (-5) + coeff2 * a2 (-5)
                        w       (i, j, k)    = a0 (-6) +
     &                       coeff * a1 (-6) + coeff2 * a2 (-6)
                        d_cap_hat (i, j, k)  = a0 (-7) +
     &                       coeff * a1 (-7) + coeff2 * a2 (-7)
                        tau_hat   (i, j, k)  = a0 (-8) +
     &                       coeff * a1 (-8) + coeff2 * a2 (-8)


                        rho (i, j, k) = a0 (1) +
     &                       coeff * a1 (1) + coeff2 * a2 (1)
                        v_1 (i, j, k) = a0 (2) +
     &                       coeff * a1 (2) + coeff2 * a2 (2)
                        v_2 (i, j, k) = a0 (3) +
     &                       coeff * a1 (3) + coeff2 * a2 (3)
                        v_3 (i, j, k) = a0 (4) +
     &                       coeff * a1 (4) + coeff2 * a2 (4)
                        eps (i, j, k) = a0 (5) +
     &                       coeff * a1 (5) + coeff2 * a2 (5)
                        p   (i, j, k) = a0 (6) +
     &                       coeff * a1 (6) + coeff2 * a2 (6)
                        gamm(i, j, k) = a0 (7) +
     &                       coeff * a1 (7) + coeff2 * a2 (7)
#ifdef CFC_MHD
                        b_1 (i, j, k) = a0 (8) +
     &                       coeff * a1 (8) + coeff2 * a2 (8)
                        b_2 (i, j, k) = a0 (9) +
     &                       coeff * a1 (9) + coeff2 * a2 (9)
                        b_3 (i, j, k) = a0 (10) +
     &                       coeff * a1 (10) + coeff2 * a2 (10)
                        psi (i, j, k) = a0 (11) +
     &                       coeff * a1 (11) + coeff2 * a2 (11)
#endif /* CFC_MHD */                        
                     END DO

                     DO i_nuc = 1, config%qn
                        DO i = 1, n
                           a0 (1) = 0.25_rk *
     &                          (6.0_rk * qxnu(i, j, i_nuc)
     &                          - qlxnu (i, j, i_nuc)
     &                          - qrxnu (i, j, i_nuc))
                           a1 (1) = qrxnu (i, j, i_nuc)
     &                          - qlxnu (i, j, i_nuc)
                           a2 (1) = 3.0_rk *
     &                          (qlxnu (i, j, i_nuc)
     &                          + qrxnu (i, j, i_nuc)-
     &                          2.0_rk * qxnu(i, j, i_nuc))
                           xnnu (i, j, k, i_nuc) = a0 (1) +
     &                          coeff * a1 (1) + coeff2 * a2 (1)
                        END DO
                     END DO
                     
                  END DO
               END DO
                           
#else /* THIRD_ORDER_PROLONGATION */
c     piecewise linear reconstruction
               
c     ------------------------------------------------------------------
c     Determine slope in super-cell and reconstruct sub-cells (q and
c     qxnu are overwritten in the process and contain the slopes in the end)

               DO j = ji1, jf1
                  !print *,'slopes',j,ji1,jf1
                  DO l = - 8, nvars - 9  !conserved variables
                     DO i = 1, m
                        dl = q  (i, j, l) - ql (i, j, l)
                        dr = qr (i, j, l) - q  (i, j, l)
                        ql (i, j, l) = dl
                        qr (i, j, l) = dr
#ifdef MC_LIMITER
c     Monotonized central limiter
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = MIN (ABS (0.5_rk * (dl + dr)),
     &                          2.0_rk * MIN (ABS (dl), ABS (dr)))
                           slope  = SIGN (slope, dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef VAN_LEER_LIMITER
c     van Leer's smooth slope limiter 
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = 2.0 * dl * dr / (dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef SUPERBEE_LIMITER
                        slope =
     &                       (sign (1.0_rk, dl) + sign (1.0_rk, dr)) *
     &                       min(abs (dl), abs (dr),
     &                       0.5_rk * max (abs (dl), abs (dr)))
#endif
                        q (i, j, l) = slope !* scal

                     END DO
                  END DO !l-loop

                  DO i_nuc = 1, config%qn !mass fractions and Y_e
                     DO i = 1, m

                        dl = qxnu  (i, j, i_nuc) - qlxnu (i, j, i_nuc)
                        dr = qrxnu (i, j, i_nuc) - qxnu  (i, j, i_nuc)
#ifdef MC_LIMITER
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = MIN (ABS (0.5_rk * (dl + dr)),
     &                          2.0_rk * MIN (ABS (dl), ABS (dr)))
                           slope  = SIGN (slope, dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef VAN_LEER_LIMITER
                        IF (dl * dr .GT. 0.0_rk) THEN
                           slope = 2.0 * dl * dr / (dl + dr)
                        ELSE
                           slope = 0.0_rk
                        END IF
#endif
#ifdef SUPERBEE_LIMITER
                        slope =
     &                       (sign (1.0_rk, dl) + sign (1.0_rk, dr)) *
     &                       min(abs (dl), abs (dr),
     &                       0.5_rk * max (abs (dl), abs (dr)))
#endif
                        qxnu (i, j, i_nuc) = slope 
                        
                     END DO
                  END DO !species-loop

c     ------------------------------------------------------------------
c     Linear reconstruction of conserved variables

                  DO k = o_s + width * kk, o_s + width * (kk + 1) - 1_ik

                     l = MOD (k - 1, width0)

c     The integral of the pre-reconstructed variable U over sub-cells
c     can be expressed analytically as <U>_supercell + coeff * dU/dvarphi:
                     coeff = (1_ik + 2_ik * l - width0) /
     &                    REAL (2_ik * width0, KIND=rk)

                     DO i = 1, m

                        entropy(i, j, k)    = entropy(i, j, k) +
     &                       coeff * q (i, j,  0)
                        cpot   (i, j, k, 1) = cpot   (i, j, k, 1) +
     &                       coeff * q (i, j, -1)
                        cpot   (i, j, k, 2) = cpot   (i, j, k, 2) +
     &                       coeff * q (i, j, -2)
                        cpot   (i, j, k, 3) = cpot   (i, j, k, 3) +
     &                       coeff * q (i, j, -3)
                        cpot   (i, j, k, 4) = cpot   (i, j, k, 4) +
     &                       coeff * q (i, j, -4)
                        t      (i, j, k)    = t      (i, j, k) +
     &                       coeff * q (i, j, -5)
                        w      (i, j, k)    = w      (i, j, k) +
     &                       coeff * q (i, j, -6)
                        d_cap_hat (i, j, k) = d_cap_hat (i, j, k) +
     &                       coeff * q (i, j, -7)
                        tau_hat   (i, j, k) = tau_hat   (i, j, k) +
     &                       coeff * q (i, j, -8)

                        rho (i, j, k) = rho (i, j, k) +
     &                       coeff * q (i, j, 1)
                        v_1 (i, j, k) = v_1   (i, j, k) +
     &                       coeff * q (i, j, 2)
                        v_2 (i, j, k) = v_2   (i, j, k) +
     &                       coeff * q (i, j, 3)
                        v_3 (i, j, k) = v_3   (i, j, k) +
     &                       coeff * q (i, j, 4)
                        eps (i, j, k) = eps   (i, j, k) +
     &                       coeff * q (i, j, 5)
                        p   (i, j, k) = p     (i, j, k) +
     &                       coeff * q (i, j, 6)
                        gamm(i, j, k) = gamm  (i, j, k) +
     &                       coeff * q (i, j, 7)
#ifdef CFC_MHD
                        b_1 (i, j, k) = b_1 (i, j, k) +
     &                       coeff * q (i, j, 8)
                        b_2 (i, j, k) = b_2 (i, j, k) +
     &                       coeff * q (i, j, 9)
                        b_3 (i, j, k) = b_3 (i, j, k) +
     &                       coeff * q (i, j, 10)
                        psi (i, j, k) = psi (i, j, k) +
     &                       coeff * q (i, j, 11)
#endif /* CFC_MHD */                        
                     END DO

                  END DO        !loop over sub-cells

c     ------------------------------------------------------------------
c     Mass fractions:
c     Step 1: Ensure that  \sum_i dX_i/dvarphi = 1 by rescaling
                  splus  (1:m) = 0.0_rk
                  sminus (1:m) = 0.0_rk
                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        IF (qxnu (i, j, i_nuc) .GT. 0.0_rk) THEN
                           splus (i) = splus (i) +
     &                          qxnu (i, j, i_nuc)
                        ELSE
                           sminus (i) = sminus (i) +
     &                          qxnu (i, j, i_nuc)
                        END IF
                     END DO
                  END DO
                  
                  ssum (:) = splus (:) + sminus (:)
                  WHERE (ssum .GT. 0.0_rk)
                     cflat = - sminus / splus
                  ELSE WHERE (ssum .LT. 0.0_rk)
                     cflat = - splus  / sminus
                  ELSE WHERE
                     cflat = 1.0_rk
                  END WHERE

                  DO i_nuc = 1, config%qn - 1
                     DO i = 1, m
                        IF (ssum (i) .GT. 0.0_rk) THEN
                           qxnu (i, j, i_nuc) = MAX (cflat (i) *
     &                          qxnu (i, j, i_nuc), 0.0_rk) +
     &                          MIN (qxnu (i, j, i_nuc), 0.0_rk) 
                        ELSE
                           qxnu (i, j, i_nuc) = MIN (cflat (i) *
     &                          qxnu (i, j, i_nuc), 0.0_rk) +
     &                          MAX (qxnu (i, j, i_nuc), 0.0_rk) 
                        END IF
                     END DO
                  END DO

c     Step 2: Linear reconstruction of mass fractions
                  DO k = o_s + width * kk, o_s + width * (kk + 1) - 1_ik

                     l = MOD (k - 1, width0)

                     coeff = (1_ik + 2_ik * l - width0) /
     &                    REAL (2_ik * width0, KIND=rk)

                     DO i_nuc = 1, config%qn
                        DO i = 1, m
                           xnnu (i, j, k, i_nuc) =
     &                          xnnu  (i, j, k, i_nuc) +
     &                          coeff * qxnu (i, j, i_nuc)
                        END DO
                     END DO

                  END DO        !loop over sub-cells (varphi direction)
               END DO           !loop over theta-index
#endif /* THIRD_ORDER_PROLONGATION */


            END DO              !loop over super-cells

            
c            CALL MPI_Barrier (ierr)


            DEALLOCATE (q, qxnu, ql, qr, qlxnu, qrxnu)

         END IF

      END DO !loop over regions

      
      RETURN

      END SUBROUTINE prolongation_output

c     ==================================================================
      


      END MODULE mesh_coarsening

c     ==================================================================



