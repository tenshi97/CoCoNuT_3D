c     ==================================================================

      MODULE cons_check

c     ==================================================================

      IMPLICIT NONE

      CONTAINS

      

c     ==================================================================

      SUBROUTINE do_conservation_check

c     ==================================================================

      USE precision

      USE size_cfc
      USE parameters_cfc
      USE metric_cfc
      USE conserved_cfc
      USE hydro_primitives_cfc
      USE sources_cfc
      USE nutra_cfc

      USE totare_hy, ONLY: dvxtot,dvytot,dvztot
      USE gfloat_hy, ONLY: vlfrac

#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, jk

#ifdef MPI_HYDRO
      real (kind=rk) :: masses (2), masses2 (2)
      integer (kind=ik) :: ierr
#endif

#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
      m_rest = 0.0_rk
      m_grav_1 = 0.0_rk
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE nowait
#endif


#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP DO REDUCTION (+: m_rest,m_grav_1)
#endif
      do jk = 1, n_loc * o_loc
            
         k = int((jk + n_loc - 1) / n_loc )
         j = (n_s - 1) + (jk - (k - 1) * n_loc)
         k = k + o_s - 1
            
         do i = 1, m

            m_rest = m_rest +
     &           equatorial_symmetry_factor *
     &           d_cap_hat (i, j, k) *
     &           dvxtot (i) * dvytot(j) * dvztot (k) *
     &           vlfrac
            m_grav_1 = m_grav_1 +
     &           equatorial_symmetry_factor *
     &           (2.0_rk / phi (i, j, k) *
     &           (tau_hat (i, j, k) + d_cap_hat (i, j, k) +
     &           enu (i, j, k)) +
     &           0.125_rk / pi / phi (i, j, k) ** 7 *
     &           k_cap_hat_squared (i, j, k)) *
     &           dvxtot (i) * dvytot (j) * dvztot (k) *
     &           vlfrac

         enddo
      enddo

#ifdef MPI_HYDRO
      masses(1)=m_grav_1
      masses(2)=m_rest
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP SINGLE
#endif
      call MPI_Allreduce (masses, masses2, 2,
     &     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      masses=masses2
#if defined(OPENMP_CFC) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$OMP END SINGLE
#endif
      m_grav_1=masses(1)
      m_rest  =masses(2)
#endif

      print *,'rest mass',m_rest
      
ccc   store initial values
      if (iteration .eq. 1) then
         m_grav_1_ini = m_grav_1
         m_rest_ini = m_rest
      endif
         
      return

      END SUBROUTINE do_conservation_check

c     ==================================================================



      END MODULE cons_check

c     ==================================================================
