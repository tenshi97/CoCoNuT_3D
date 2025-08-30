c     =================================================================

      MODULE fft_mpi

c     =================================================================

      IMPLICIT NONE

      CONTAINS
     

c     =================================================================

      SUBROUTINE fft (x_in, y_out, len, mode, index)

c     =================================================================

      USE precision
      USE configure
      USE mo_mpi

      IMPLICIT NONE

      INTEGER (KIND=ik), INTENT(IN) :: len
      REAL (KIND=rk), INTENT (IN) ::
     &     x_in  (1:len,qy_s:qy_e,qz_s:qz_e)
      REAL (KIND=rk), INTENT (OUT) ::
     &     y_out (1:len,qy_s:qy_e,qz_s:qz_e)
      INTEGER (KIND=ik), INTENT (OUT) ::  index (0:config%qz-1)

      INTEGER (KIND=ik), INTENT (in) :: mode
      
      INTEGER (KIND=ik) :: logn

      DOUBLE PRECISION :: pi
      DOUBLE PRECISION :: arg

      DOUBLE COMPLEX :: y  (1:len/2,qy_s:qy_e,qz_s-1:qz_e-1)
      DOUBLE COMPLEX :: y0 (1:len/2,qy_s:qy_e,qz_s-1:qz_e-1)

      DOUBLE COMPLEX :: om, omi
c     we never have more than 2**21 cells in varphi
      INTEGER (KIND=ik) :: binary (0:20)

      INTEGER (KIND=ik) :: i, j, jj, ji, jf, k, kk, jk, inx, nk
      INTEGER (KIND=ik) :: di, di0, di1
      INTEGER (KIND=ik) :: n_loc, o_loc

      INTEGER (KIND=ik) :: ndat, ierr, src, dest
      INTEGER (KIND=ik) :: mpistat (MPI_STATUS_SIZE)

      LOGICAL :: send_right


      IF (MOD (len, 2_ik) .NE. 0) THEN
         STOP 'Error: FFT needs even number of radial grids points'
      END IF

      pi = ACOS (-1.0d0)
     
      n_loc = qy_e - qy_s + 1
      o_loc = qz_e - qz_s + 1
      ndat = (len/2) * n_loc * o_loc

c     Determine log_2 (qz), integer arithmetic
      j = config%qz
      logn = 0
      DO i = 1, config%qz
         j = j / 2
         IF (j .EQ. 0) EXIT
         logn = logn + 1
         !PRINT *,'log_2 (qz) =',logn
      END DO

      DO i = 0, config%qz-1
         j = i
         DO k = logn-1, 0, -1
            binary (k) = j / 2 ** k
            j = j - binary (k) * 2 ** k
         END DO

         j=0
         DO k = 0, logn-1
            j = j + binary (logn-1-k) * 2 ** k
         END DO
         index (i) = j
      END DO


c     Pack the array that we transform into a complex array
      DO jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (qy_s - 1) + (jk - (k - 1) * n_loc)
         k = k + qz_s - 1

         DO i = 1, len / 2
            y0 (i,j,k-1) = DCMPLX (x_in (2*i-1,j,k), x_in (2*i,j,k))
c            y0 (i,j,k-1) = DCMPLX (x_in (i,j,k), 0.0_rk)
         END DO

      END DO

      IF (mode .EQ. 1) THEN !forward transform
      
         DO k = 1, logn
         
            di = 2 ** (logn - k)
            nk = 2 ** (k-1)
            arg = - DBLE (2.0d0*pi/DBLE (2*di))
            om = DCMPLX (COS (arg), SIN (arg))
            !PRINT *,'lev',k,di

            IF (2*di .LE. o_loc) THEN !no MPI communication needed

               ji = (qz_s - 1) / (2*di)
               jf = (qz_e - 1) / (2*di)

               DO j = ji, jf
               !PRINT *,'j=',j
                  di0 = 2*di*j
                  di1 = 2*di*j+di
                  DO i = 0, di-1
                     omi = DCMPLX (COS (i * arg), SIN (i * arg))
                     y (:,:,di0+i) =
     &                    y0 (:,:,di0+i) + y0 (:,:,di1+i)
                     y (:,:,di1+i) =
     &                    (y0 (:,:,di0+i) - y0 (:,:,di1+i)) * omi
                  END DO
               END DO

            ELSE !MPI communication needed

#ifdef MPI_HYDRO
               IF (MOD ((qz_s-1)/di, 2) .EQ. 0) THEN

                  send_right = .TRUE.
                  CALL MPI_Cart_shift (cart_comm, 1_ik,  di / o_loc,
     &              src, dest, ierr)
                  !PRINT *,'send',src, dest,di,o_loc
                  CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &                 dest, 6000_ik, cart_comm, ierr)
                  !PRINT *,'sent'
                  CALL MPI_Recv (y,  ndat, MPI_DOUBLE_COMPLEX,
     &                 dest, 6001_ik, cart_comm, mpistat, ierr) 
                  !PRINT *,'received'

                  DO i = qz_s-1, qz_e-1
                     inx = MOD (i, di)
                     omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                     y (:,:,i) =
     &                    y0 (:,:,i) + y (:,:,i)
                  END DO

               ELSE

                  send_right = .FALSE.
                  CALL MPI_Cart_shift (cart_comm, 1_ik, di / o_loc,
     &              src, dest, ierr)
!                  !PRINT *,'send (b)',src,dest
                  CALL MPI_Recv (y, ndat, MPI_DOUBLE_COMPLEX,
     &                 src, 6000_ik, cart_comm, mpistat, ierr)
c                  !PRINT *,'received (b)'
                  CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &                 src, 6001_ik, cart_comm, ierr) 
c                  !PRINT *,'sent (b)'
                  DO i = qz_s-1, qz_e-1
                     inx = MOD (i, di)
                     omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                     y (:,:,i) =
     &                    (y (:,:,i) - y0 (:,:,i)) * omi !om ** i
                  END DO
               END IF
#endif /* MPI_HYDRO */

            END IF

#if 0
            DO jk = 1, n_loc * o_loc
               kk = int((jk + n_loc - 1) / n_loc )
               jj = (qy_s - 1) + (jk - (kk - 1) * n_loc)
               kk = kk + qz_s - 2
               y0 (:,jj,kk) = y (:,jj,kk)
            END DO
#endif
            y0 = y
            
         END DO


c     -----------------------------------------------------------------
      ELSE IF (mode .EQ. -1) THEN ! backward transform

         DO jk = 1, n_loc * o_loc
            k = int((jk + n_loc - 1) / n_loc)
            j = (qy_s - 1) + (jk - (k - 1) * n_loc)
            k = k + qz_s - 2
            y0 (:,j,k) = y0 (:,j,k) / DBLE (config%qz)
         END DO

         DO k = logn, 1, -1

            di = 2 ** (logn - k)
            nk = 2 ** (k-1)
            arg = DBLE (2.0d0*pi/DBLE (2*di))
            om = DCMPLX (COS (arg), SIN (arg))
            !PRINT *,'lev',k,di

            IF (2*di .LE. o_loc) THEN !no MPI communication needed

#ifdef MPI_HYDRO
               ji = (qz_s - 1) / (2*di)
               jf = (qz_e - 1) / (2*di)

               DO j = ji, jf
               !PRINT *,'j=',j
                  
                  di0 = 2*di*j
                  di1 = 2*di*j+di
                  DO i = 0, di-1
                     omi = DCMPLX (COS (i * arg), SIN (i * arg))
                     y (:,:,di0+i) = !0.5_rk * 
     &                    (y0 (:,:,di0+i) + y0 (:,:,di1+i) * omi)
                     y (:,:,di1+i) = !0.5_rk *
     &                    (y0 (:,:,di0+i) - y0 (:,:,di1+i) * omi)
                  END DO
               END DO

            ELSE !MPI communication needed

               IF (MOD ((qz_s-1)/di, 2) .EQ. 0) THEN

                  send_right = .TRUE.
                  CALL MPI_Cart_shift (cart_comm, 1_ik,  di / o_loc,
     &              src, dest, ierr)
                  CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &                 dest, 6000_ik, cart_comm, ierr)
                  CALL MPI_Recv (y, ndat, MPI_DOUBLE_COMPLEX,
     &                 dest, 6001_ik, cart_comm, mpistat, ierr)

                  DO i = qz_s-1, qz_e-1
                     inx = MOD (i, di)
                     omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                     y (:,:,i) = !0.5_rk *
     &                    (y0 (:,:,i) + y (:,:,i) * omi)
                  END DO

               ELSE

                  send_right = .FALSE.
                  CALL MPI_Cart_shift (cart_comm, 1_ik, di / o_loc,
     &              src, dest, ierr)
                  CALL MPI_Recv (y, ndat, MPI_DOUBLE_COMPLEX,
     &                 src, 6000_ik, cart_comm, mpistat, ierr)
                  CALL MPI_Send (y0, ndat, MPI_DOUBLE_COMPLEX,
     &                 src, 6001_ik, cart_comm, ierr) 
                  DO i = qz_s-1, qz_e-1
                     inx = MOD (i, di)
                     omi = DCMPLX (COS (inx * arg), SIN (inx * arg))
                     y (:,:,i) = !0.5_rk *
     &                    (y (:,:,i) - y0 (:,:,i) * omi)
                  END DO
               END IF

#endif /* MPI_HYDRO */

            END IF

#if 0
            DO jk = 1, n_loc * o_loc
               kk = int((jk + n_loc - 1) / n_loc )
               jj = (qy_s - 1) + (jk - (kk - 1) * n_loc)
               kk = kk + qz_s - 2
               y0 (:,jj,kk) = y (:,jj,kk)
            END DO
#endif
            y0 = y
            
         END DO

      ELSE

         STOP 'FFT: mode not implemented!'

      END IF


c     Unpack array:
      DO jk = 1, n_loc * o_loc
         k = int((jk + n_loc - 1) / n_loc )
         j = (qy_s - 1) + (jk - (k - 1) * n_loc)
         k = k + qz_s - 1

         DO i = 1, len / 2
            y_out (2*i-1,j,k) = DBLE  (y (i,j,k-1))
            y_out (2*i  ,j,k) = AIMAG (y (i,j,k-1))
c            y_out (i  ,j,k) = DBLE (y (i,j,k-1))
         END DO
         
      END DO


      RETURN

      END SUBROUTINE fft
      
c     =================================================================
      


      END MODULE fft_mpi
      
c     =================================================================
