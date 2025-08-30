c     ==================================================================

      MODULE gw3d

c     ==================================================================

      USE precision

      IMPLICIT NONE

      REAL (KIND=rk) :: qijq1(1:3,1:3)
      REAL (KIND=rk) :: qijq2(1:3,1:3)

      SAVE


      CONTAINS
      

c     ==================================================================

      SUBROUTINE extract_gw (dt)
      
c     ==================================================================

      USE precision

      USE conserved_cfc, ONLY: s_1_hat, s_2_hat, s_3_hat,
     &     s_1_hat_old, s_2_hat_old, s_3_hat_old
      USE grid_cfc, ONLY: r
      USE hydro_primitives_cfc, ONLY: v_1, v_2, v_3
      USE metric_cfc, ONLY: alpha, beta_up_1, beta_up_2, beta_up_3,
     &     g_up_11, g_up_22, g_up_33
      USE parameters_cfc
      USE perm_aux_cfc, ONLY: r_inv, sinus_theta, sinus_theta_inv,
     &     cosin_theta, sinus_varphi, cosin_varphi
      USE size_cfc
      USE sources_cfc, ONLY: s_total

      USE totare_hy, ONLY: dvxtot,dvytot,dvztot
      USE gfloat_hy, ONLY: vlfrac

#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      IMPLICIT NONE

      REAL (KIND=rk), INTENT (IN) :: dt

      INTEGER (KIND=ik) :: i, j, k, jk
      REAL (KIND=rk) :: scr1, scr2, scr3
      REAL (KIND=rk) :: sx, sy, sz 
      REAL (KIND=rk) :: sox, soy, soz 
      REAL (KIND=rk) :: vx, vy, vz
      REAL (KIND=rk) :: fx, fy, fz
      REAL (KIND=rk) :: xx, yy, zz
      REAL (KIND=rk) :: cth, sth, cph, sph, dv, dvom

#ifdef MPI_HYDRO
      REAL (KIND=rk) :: qijtmp (1:3,1:3,2)
      REAL (KIND=rk) :: qijtmp2(1:3,1:3,2)
      INTEGER (KIND=ik) :: ierr
#endif

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
      qijq1 (:,:) = 0.0_rk
      qijq2 (:,:) = 0.0_rk
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE nowait
#endif


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO REDUCTION (+: qijq1,qijq2)
#endif
      DO jk = 1, n_loc * o_loc

         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1

         sth = sinus_theta (j)
         cth = cosin_theta (j)
         cph = cosin_varphi (k)
         sph = sinus_varphi (k)

         dvom = equatorial_symmetry_factor * vlfrac *
     &        dvytot (j) * dvztot (k)

         DO i = 1, m

c     Transform \hat{S}_i and \hat{v}^i to Cartesian basis 
            scr1 = s_1_hat (i, j, k)
            scr2 = s_2_hat (i, j, k) *
     &           r_inv (i)
            scr3 = s_3_hat (i, j, k) *
     &           r_inv (i) * sinus_theta_inv (j)

            sx =   scr1  * sth * cph
     &           + scr2  * cth * cph
     &           - scr3        * sph
            sy =   scr1  * sth * sph
     &           + scr2  * cth * sph
     &           + scr3        * cph
            sz =   scr1  * cth
     &           - scr2  * sth

c     Same for  source terms
            scr1 = s_total (i, j, k, 2)
            scr2 = s_total (i, j, k, 3) *
     &           r_inv (i)
            scr3 = s_total (i, j, k, 4) *
     &           r_inv (i) * sinus_theta_inv (j)

            fx = alpha (i, j, k) * 
     &           ( scr1  * sth * cph
     &           + scr2  * cth * cph
     &           - scr3        * sph)
            fy = alpha (i, j, k) *
     &           ( scr1  * sth * sph
     &           + scr2  * cth * sph
     &           + scr3        * cph)
            fz = alpha (i, j, k) * 
     &           ( scr1  * cth
     &           - scr2  * sth)
            
c     This assumes that we're calling this routine after
c     rescale_velocties (2)
            scr1 = g_up_11 (i, j, k) * v_1 (i, j, k) -
     &           beta_up_1 (i, j, k) / alpha (i, j, k)
            scr2 = g_up_22 (i, j, k) * v_2 (i, j, k) -
     &           beta_up_2 (i, j, k) / alpha (i, j, k)
            scr3 = g_up_33 (i, j, k) * v_3 (i, j, k) -
     &           beta_up_3 (i, j, k) / alpha (i, j, k)

            scr2 = scr2 * r (i)
            scr3 = scr3 * r (i) * sth

            vx =   scr1  * sth * cph
     &           + scr2  * cth * cph
     &           - scr3        * sph
            vy =   scr1  * sth * sph
     &           + scr2  * cth * sph
     &           + scr3        * cph
            vz =   scr1  * cth
     &           - scr2  * sth

            xx = r (i) * sth * cph
            yy = r (i) * sth * sph
            zz = r (i) * cth

c     S_hat at last time step -- for time-integrated quadrupole formula
            scr1 = s_1_hat_old (i, j, k)
            scr2 = s_2_hat_old (i, j, k) *
     &           r_inv (i)
            scr3 = s_3_hat_old (i, j, k) *
     &           r_inv (i) * sinus_theta_inv (j)

            sox =  scr1  * sth * cph
     &           + scr2  * cth * cph
     &           - scr3        * sph
            soy =  scr1  * sth * sph
     &           + scr2  * cth * sph
     &           + scr3        * cph
            soz =  scr1  * cth
     &           - scr2  * sth


            dv = dvxtot (i) * dvom
            
c     Contribution to integral in stress formula
c     Let us not symmetrize for the time being
            qijq2 (1, 1) = qijq2 (1, 1) +
     &           (sx * vx + xx * fx) * dv
            qijq2 (1, 2) = qijq2 (1, 2) +
     &           (sx * vy + xx * fy) * dv
            qijq2 (1, 3) = qijq2 (1, 3) +
     &           (sx * vz + xx * fz) * dv
            qijq2 (2, 1) = qijq2 (2, 1) +
     &           (sy * vx + yy * fx) * dv
            qijq2 (2, 2) = qijq2 (2, 2) +
     &           (sy * vy + yy * fy) * dv
            qijq2 (2, 3) = qijq2 (2, 3) +
     &           (sy * vz + yy * fz) * dv
            qijq2 (3, 1) = qijq2 (3, 1) +
     &           (sz * vx + zz * fx) * dv
            qijq2 (3, 2) = qijq2 (3, 2) +
     &           (sz * vy + zz * fy) * dv
            qijq2 (3, 3) = qijq2 (3, 3) +
     &           (sz * vz + zz * fz) * dv

            qijq1 (1, 1) = qijq1 (1, 1) +
     &           xx * (sx - sox) * dv
            qijq1 (1, 2) = qijq1 (1, 2) +
     &           xx * (sy - soy) * dv
            qijq1 (1, 3) = qijq1 (1, 3) +
     &           xx * (sz - soz) * dv
            qijq1 (2, 1) = qijq1 (2, 1) +
     &           yy * (sx - sox) * dv
            qijq1 (2, 2) = qijq1 (2, 2) +
     &           yy * (sy - soy) * dv
            qijq1 (2, 3) = qijq1 (2, 3) +
     &           yy * (sz - soz) * dv
            qijq1 (3, 1) = qijq1 (3, 1) +
     &           zz * (sx - sox) * dv
            qijq1 (3, 2) = qijq1 (3, 2) +
     &           zz * (sy - soy) * dv
            qijq1 (3, 3) = qijq1 (3, 3) +
     &           zz * (sz - soz) * dv

         END DO

      END DO


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif

      qijq1 = qijq1/dt

#ifdef MPI_HYDRO
      qijtmp (1:3,1:3,1) = qijq1 (1:3,1:3)
      qijtmp (1:3,1:3,2) = qijq2 (1:3,1:3)

c      CALL MPI_Allreduce (qijtmp, qijtmp2, 18,
c     &     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_Reduce (qijtmp, qijtmp2, 18,
     &     MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      IF (myproc .EQ. 0) qijtmp = qijtmp2
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif
      qijq1 (1:3,1:3) = qijtmp (1:3,1:3,1)
      qijq2 (1:3,1:3) = qijtmp (1:3,1:3,2)
#endif

c     Since all our variables are in cm (geometric units),
c     no conversion to cgs is necessary
      qijq1 (1:3,1:3) = 4.0_rk * qijq1 (1:3,1:3)
      qijq2 (1:3,1:3) = 4.0_rk * qijq2 (1:3,1:3)
c     Check that the pre-factor is indeed 4 (not 2?)

      RETURN

      END SUBROUTINE extract_gw
      
c     ==================================================================


      END MODULE gw3d
      
c     ==================================================================
