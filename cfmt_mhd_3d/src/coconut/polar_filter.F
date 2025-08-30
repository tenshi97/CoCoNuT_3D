#if 1

#undef SMOOTH_ENERGY
#undef SMOOTH_MASS_FRACTION
      
c     ==================================================================

      MODULE filtering

c     ==================================================================
      
      CONTAINS


c     ==================================================================

      SUBROUTINE polar_filter (mode)

c     ==================================================================
c     - Filter fluxes and pressure on varphi interfaces.
c     ------------------------------------------------------------------

      USE precision
      
      USE conserved_cfc
      USE fluxes_cfc
      USE grid_cfc
      USE hydro_areas_mod
#ifdef CFC_MHD
      USE hydro_primitives_cfc, ONLY: p, psi
#else
      USE hydro_primitives_cfc, ONLY: p
#endif /* CFC_MHD */
c      USE conserved_cfc, ONLY: b_cap_1_hat, b_cap_2_hat, b_cap_3_hat
      USE perm_aux_cfc, ONLY: delta_theta_inv, sinus_theta,
     &     sinus_theta_if
      USE phycon, ONLY: pc_pi
      USE size_cfc

      USE fft_mpi
      USE mo_mpi

      USE configure

      USE, INTRINSIC :: iso_c_binding

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: mode

#ifdef CFC_MHD
      REAL (KIND=rk) :: xin (1:m*(10+config%qn), n_s:n_e, o_s:o_e)
      REAL (KIND=rk) :: xout (1:m*(10+config%qn), n_s:n_e, o_s:o_e)
      REAL (KIND=rk) :: sbufz (1:m, n_s:n_e, (10+config%qn))
      REAL (KIND=rk) :: rbufz (1:m, n_s:n_e, (10+config%qn))
#else /* CFC_MHD */
      REAL (KIND=rk) :: xin (1:m*(6+config%qn), n_s:n_e, o_s:o_e)
      REAL (KIND=rk) :: xout (1:m*(6+config%qn), n_s:n_e, o_s:o_e)
      REAL (KIND=rk) :: sbufz (1:m, n_s:n_e, (6+config%qn))
      REAL (KIND=rk) :: rbufz (1:m, n_s:n_e, (6+config%qn))
#endif /* CFC_MHD */

      INTEGER (KIND=ik) :: src, dest, ierr, mpistat (MPI_STATUS_SIZE)
      INTEGER (KIND=ik), parameter :: tag_ysndl=1500_ik
      INTEGER (KIND=ik), parameter :: tag_ysndr=1501_ik
      INTEGER (KIND=ik), parameter :: tag_zsndl=1502_ik
      INTEGER (KIND=ik), parameter :: tag_zsndr=1503_ik

      REAL (KIND=rk), PARAMETER :: lg2 = LOG (2.0_rk)
      REAL (KIND=rk) :: wn, wn0, rspher, scr, sth, sfac
      REAL (KIND=rk) :: dsmooth (1:m)
      INTEGER (KIND=ik) :: index (1:o)

      INTEGER (KIND=ik) :: i, ii, j, k, jk, l, ixf, ioy, nvars
      INTEGER (KIND=ik) :: msmooth (n_s:n_e)

#ifdef CFC_MHD      
      nvars = 10+config%qn
#else
      nvars = 6+config%qn
#endif /* CFC_MHD */
      
      ixf    = areas%ix_are(1, 2)
      ioy    = areas%ix_are(1, 6)
      if (ioy.eq.n) then
         rspher = r (ixf + 1)
      else
         rspher = 1e6
      end if

      DO jk = 1, n_loc * o_loc
         k = INT((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         IF (mode .EQ. 0) THEN
            DO l = 1, 5
               xin (1+(l-1)*m:l*m,j,k)=f_num_if(1:m,j,k,l)
            END DO
            DO l = 1, config%qn
               xin (1+(l+4)*m:(l+5)*m,j,k)=f_num_xnu_if(1:m,j,k,l)
            END DO
            xin (1+(config%qn+5)*m:(config%qn+6)*m,j,k)=
     &           pav_varphi(1:m,j,k)
#ifdef CFC_MHD
            do l = 1, 4
               xin (1+(config%qn+5+l)*m:(config%qn+6+l)*m,j,k)=
     &              f_num_if(1:m,j,k,5+l)
            end do
#endif            
         ELSE
            xin (1+0*m:1*m,j,k)=d_cap_hat(1:m,j,k)
            xin (1+1*m:2*m,j,k)=s_1_hat  (1:m,j,k)
            xin (1+2*m:3*m,j,k)=s_2_hat  (1:m,j,k)
            xin (1+3*m:4*m,j,k)=s_3_hat  (1:m,j,k)
            xin (1+4*m:5*m,j,k)=tau_hat  (1:m,j,k)
#ifdef CFC_MHD
     &           - 0.5_rk *
     &          (b_cap_1_hat (1:m, j, k) ** 2 +
     &           b_cap_2_hat (1:m, j, k) ** 2 + 
     &           b_cap_3_hat (1:m, j, k) ** 2)            
#endif /* CFC_MHD */            
c            print*,'tau_hat(vorher)',j,k,tau_hat(177,j,k)
            DO l = 1, config%qn
               xin (1+(l+4)*m:(l+5)*m,j,k)=d_cap_xnu_hat(1:m,j,k,l)
#ifdef SMOOTH_MASS_FRACTIONS
     &              / d_cap_hat(1:m,j,k)
#endif               
            END DO
#ifdef CFC_MHD
            xin (1+(config%qn+5)*m:(config%qn+6)*m,j,k)=
     &           eps_aux(1:m,j,k)
            xin (1+(config%qn+6)*m:(config%qn+7)*m,j,k)=
     &           b_cap_1_hat  (1:m,j,k)
            xin (1+(config%qn+7)*m:(config%qn+8)*m,j,k)=
     &           b_cap_2_hat  (1:m,j,k)
            xin (1+(config%qn+8)*m:(config%qn+9)*m,j,k)=
     &           b_cap_3_hat  (1:m,j,k)
            xin (1+(config%qn+9)*m:(config%qn+10)*m,j,k)=
     &           psi  (1:m,j,k)
#endif            
         END IF

      END DO

      CALL FFT (xin (1:m*nvars, n_s:n_e, o_s:o_e),
     &     xout (1:m*nvars, n_s:n_e, o_s:o_e), m*nvars, 
     &     1, index)
         
      DO jk = 1, n_loc * o_loc
         k = INT((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         wn = REAL (MIN (index (k), o - index (k)), KIND=rk)
         wn = MAX (wn, 1e-35_rk)
         sth = MIN (sinus_theta_if (j-1), sinus_theta_if (j)) 
c         sth = MIN (sinus_theta (j-1), sinus_theta (j+1)) 
         sth = MAX (sth, sinus_theta (1))

         IF (mode .LE. 1) THEN
            msmooth (j) = m
            DO i = 1, m
               IF (sth * delta_varphi * 
     &              delta_theta_inv *
     &              r (i) / rspher * 2.0_rk
     &              * MAX(0.7_rk / config%cfl, 1.0_rk)
     &              .LT. 1.0_rk)
     &              msmooth (j) = i
            END DO

            DO ii = 0, (nvars-1)*m, m
               DO i = 1, m
                  xin (i+ii, j, k) =
     &                 xout (i+ii, j, k) *
     &                 MIN (sth * delta_varphi * 
     &                 delta_theta_inv *
     &                 r (i) / rspher /
     &                 SIN (0.5_rk * wn * delta_varphi) * 2.0_rk
c     &                 SIN (0.5_rk * wn * delta_varphi) * 0.5_rk
c     &                 * MAX(0.7_rk / config%cfl, 1.0_rk)
     &                 , 1.0_rk) ** 2
c                  IF (i.eq.130) print*,'fft',j,wn,ii/m,xin(i+ii,j,k),
c     &                 xout(i+ii,j,k)
               END DO
               IF ((j.EQ.1 .OR. j.EQ.n) .AND. wn.GE.0.9_rk)
     &              xin (ii+1:ii+msmooth(j), j, k) = 0.0_rk
c     &              xin (ii+1:ii+m, j, k) = 0.0_rk
            END DO
 
         ELSE

c     Include hyperviscosity for stable modes to damp the
c     Gibbs phenomenon

            msmooth (j) = m
            DO i = 1, m
               IF (sth * delta_varphi * 
     &              delta_theta_inv *
     &              r (i) / rspher !* 2.0_rk
c     &              * MAX(0.7_rk / config%cfl, 1.0_rk)
     &              .LT. 1.0_rk)
     &              msmooth (j) = i
            END DO
c            PRINT *,'msmooth',j,k,msmooth(j)

            DO ii = 0, (nvars-1)*m, m
               DO i = 1, msmooth(j)
                  sfac =  MIN ((sth * delta_varphi * 
     &                 delta_theta_inv *
     &                 r (i) / rspher /
c     &                 (wn * delta_varphi) * 4.0_rk) ** 2 * 0.5_rk,
     &                 (wn * delta_varphi)) ** 2 * 0.5_rk,
     &                 1.0_rk)
                  wn0 = sth * delta_theta_inv *
     &                 r (i) / rspher !* 0.25_rk
                  IF (wn .LE. wn0) sfac = MIN (sfac,
     &                 EXP (-(wn/wn0) ** 6 * lg2))
                  xin (i+ii, j, k) =
     &                 xout (i+ii, j, k) * sfac
               END DO
               DO i = msmooth(j) + 1, m
                  xin (i+ii, j, k) =
     &                 xout (i+ii, j, k)
               END DO
               IF ((j.EQ.1 .OR. j.EQ.n) .AND. wn.GE.0.9_rk)
     &              xin (ii+1:ii+msmooth(j), j, k) = 0.0_rk
c     &              xin (ii+1:ii+m, j, k) = 0.0_rk
            END DO
            
         END IF

      END DO

      CALL FFT (xin (1:m*nvars, n_s:n_e, o_s:o_e),
     &     xout (1:m*nvars, n_s:n_e, o_s:o_e), m*nvars,
     &     -1, index)
      
      DO jk = 1, n_loc * o_loc
         k = INT((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         IF (mode .EQ. 0) THEN
            DO l = 1, 5
               f_num_if(1:m,j,k,l) = xout (1+(l-1)*m:l*m,j,k)
            END DO
            DO l = 1, config%qn
               f_num_xnu_if(1:m,j,k,l) = xout (1+(l+4)*m:(l+5)*m,j,k)
            END DO
            pav_varphi(1:m,j,k) = 
     &           xout (1+(config%qn+5)*m:(config%qn+6)*m,j,k)
#ifdef CFC_MHD            
            DO l = 1, 4
               f_num_if(1:m,j,k,5+l) =
     &              xout (1+(config%qn+5+l)*m:(config%qn+6+l)*m,j,k)
            END DO
#endif /* CFC_MHD */            
         ELSE
            d_cap_hat(1:m,j,k)=xout (1+0*m:1*m,j,k)
            s_1_hat  (1:m,j,k)=xout (1+1*m:2*m,j,k)
            s_2_hat  (1:m,j,k)=xout (1+2*m:3*m,j,k)
            s_3_hat  (1:m,j,k)=xout (1+3*m:4*m,j,k)
            tau_hat  (1:m,j,k)=xout (1+4*m:5*m,j,k)

            DO l = 1, config%qn
               d_cap_xnu_hat(1:m,j,k,l)=xout (1+(l+4)*m:(l+5)*m,j,k)
#ifdef SMOOTH_MASS_FRACTIONS
     &              * d_cap_hat(1:m,j,k)
#endif               
            END DO

#ifdef CFC_MHD
            eps_aux(1:m,j,k) =
     &           xout (1+(config%qn+5)*m:(config%qn+6)*m,j,k)
            b_cap_1_hat  (1:m,j,k)=
     &           xout (1+(config%qn+6)*m:(config%qn+7)*m,j,k)
            b_cap_2_hat  (1:m,j,k)=
     &           xout (1+(config%qn+7)*m:(config%qn+8)*m,j,k)
            b_cap_3_hat  (1:m,j,k)=
     &           xout (1+(config%qn+8)*m:(config%qn+9)*m,j,k)
            tau_hat      (1:m,j,k) = tau_hat  (1:m,j,k) + 0.5_rk *
     &          (b_cap_1_hat (1:m, j, k) ** 2 +
     &           b_cap_2_hat (1:m, j, k) ** 2 + 
     &           b_cap_3_hat (1:m, j, k) ** 2)            
            psi  (1:m,j,k)=
     &           xout (1+(config%qn+9)*m:(config%qn+10)*m,j,k)
c            print*,'tau_hat(nachher)',j,k,tau_hat(177,j,k)
#endif /* CFC_MHD */
         END IF
      END DO
         


      IF (mode .EQ. 0) THEN
c     Communicate required fluxes to right neighbour
         DO l = 1, 5
            sbufz(1:m,n_s:n_e,l)=
     &        f_num_if (1:m,n_s:n_e,o_e,l)
         END DO
         DO l = 1, config%qn
            sbufz(1:m,n_s:n_e,l+5)=
     &           f_num_if (1:m,n_s:n_e,o_e,l)
         END DO
         sbufz(1:m,n_s:n_e,config%qn+6) = pav_varphi (1:m,n_s:n_e,o_s-1)
#ifdef CFC_MHD
         DO l = 1, 4
            sbufz(1:m,n_s:n_e,config%qn+6+l)=
     &        f_num_if (1:m,n_s:n_e,o_e,5+l)
         END DO
#endif /* CFC_MHD */         

c     determine source and destination ranks and send data
         CALL MPI_Cart_shift(cart_comm,1,1,src,dest,ierr)
            
         IF (MOD (src, 2_ik) .EQ. 1_ik) THEN
            CALL MPI_Send (sbufz,m*n_loc*nvars,
     &           MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &           cart_comm,ierr)
            CALL MPI_Recv (rbufz,m*n_loc*nvars,
     &           MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &           mpistat,ierr)
         ELSE
            CALL MPI_Recv (rbufz,m*n_loc*nvars,
     &           MPI_DOUBLE_PRECISION,src,tag_zsndr,cart_comm,
     &           mpistat,ierr)
            CALL MPI_Send (sbufz,m*n_loc*nvars,
     &           MPI_DOUBLE_PRECISION,dest,tag_zsndr,
     &           cart_comm,ierr)
         END IF
         
         DO l = 1, 5
            f_num_if (1:m,n_s:n_e,o_s-1,l) =
     &           rbufz (1:m,n_s:n_e,l)
         END DO
         DO l = 1, config%qn
            f_num_if (1:m,n_s:n_e,o_s-1,l) =
     &           rbufz (1:m,n_s:n_e,l+5)
         END DO
         pav_varphi (1:m,n_s:n_e,o_s-1) = rbufz(1:m,n_s:n_e,config%qn+6)
#ifdef CFC_MHD
         DO l = 1, 4
            f_num_if (1:m,n_s:n_e,o_e,5+l)=
     &           rbufz(1:m,n_s:n_e,config%qn+6+l)
         END DO
#endif /* CFC_MHD */         

      END IF
      
C      STOP
      RETURN

      END SUBROUTINE polar_filter

c     ==================================================================



      END MODULE filtering

c     ==================================================================
#endif


