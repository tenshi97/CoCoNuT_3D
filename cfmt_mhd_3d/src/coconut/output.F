c     ==================================================================

      MODULE gr_output

c     ==================================================================

      IMPLICIT NONE

      CONTAINS



c     ==================================================================

      SUBROUTINE write_coconut_data

c     ==================================================================

      USE precision

      USE size_cfc
      USE parameters_cfc
      USE hydro_primitives_cfc
      USE conserved_cfc
      USE metric_cfc
      USE nutra_cfc
      
      USE nucparam, only : name_xnuc

      USE charac
#ifdef WRITE_BINARY_OUTPUT
      USE dataformat_vertex
#endif

#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      use configure
      IMPLICIT NONE

      integer (kind=ik) :: i, j, k, l
#ifdef MPI_HYDRO
      integer (kind=ik) :: fh, ierr
#endif
      type(datafile) :: fd

      call create_vertex_file(fd, rstfil_cfc, "CoCoNuT restart file")
      
      call writeq(fd,"d_cap_hat",d_cap_hat (1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian baryon mass density', 'geometric',
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"s_1_hat",s_1_hat (1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian r-momentum density', 'geometric', 
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"s_2_hat",s_2_hat (1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian theta-momentum density', 'geometric', 
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"s_3_hat",s_3_hat (1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian varphi-momentum density', 'geometric',
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"tau_hat",tau_hat (1:m,n_s:n_e,o_s:o_e),    
     &     'Eulerian energy density', 'geometric', 'radius:xzn',
     &     'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"xnnu",xnnu (1:m,n_s:n_e,o_s:o_e,1:config%qn),
     &     'mass fractions', 'geometric', 'radius:xzn',
     &     "species:" // stringlist(name_xnuc),
     &     'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e, 1,config%qn/), 
     &     (/1,m, 1,n, 1,o, 1,config%qn/))
      call writeq(fd,"phi",phi (1:m,n_s:n_e,o_s:o_e),
     &     'conformal factor', 'geometric', 'radius:xzn', 
     &     'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"alpha",alpha (1:m,n_s:n_e,o_s:o_e),
     &     'lapse function', 'geometric', 'radius:xzn', 
     &     'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"beta_up_1",beta_up_1 (1:m,n_s:n_e,o_s:o_e),
     &     'radial shift', 'geometric', 'radius:xzn', 
     &     'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"beta_up_2",beta_up_2 (1:m,n_s:n_e,o_s:o_e),
     &     'meridional shift', 'geometric',
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"beta_up_3",beta_up_3 (1:m,n_s:n_e,o_s:o_e),
     &     'zonal shift', 'geometric', 
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"enu",enu(1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian neutrino energy density', 'geometric',
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"fnu",fnu(1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian neutrino flux density', 'geometric', 
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"m_grav_1_ini",m_grav_1_ini,
     &     'gravitational mass', 'geometric')
      call writeq(fd,"m_rest_ini",m_rest_ini,
     &     'baryonic mass', 'geometric')

#ifdef CFC_MHD
      call writeq(fd,"b_cap_1_hat",b_cap_1_hat (1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian r-momentum density', 'geometric', 
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"b_cap_2_hat",b_cap_2_hat (1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian r-momentum density', 'geometric', 
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"b_cap_3_hat",b_cap_3_hat (1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian r-momentum density', 'geometric', 
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))

      b_cap_1_hat_old (1 : m, n_s : n_e, o_s : o_e) =
     &     b_1 (1 : m, n_s : n_e, o_s : o_e)
      b_cap_2_hat_old (1 : m, n_s : n_e, o_s : o_e) =
     &     b_2 (1 : m, n_s : n_e, o_s : o_e)
      b_cap_3_hat_old (1 : m, n_s : n_e, o_s : o_e) =
     &     b_3 (1 : m, n_s : n_e, o_s : o_e)

      call writeq(fd,"b_1",b_cap_1_hat_old (1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian r-momentum density', 'geometric', 
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"b_2",b_cap_2_hat_old (1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian r-momentum density', 'geometric', 
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"b_3",b_cap_3_hat_old (1:m,n_s:n_e,o_s:o_e),
     &     'Eulerian r-momentum density', 'geometric', 
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
      call writeq(fd,"psi",psi (1:m,n_s:n_e,o_s:o_e),
     &     'divergence cleaning auxilliary field', 'geometric', 
     &     'radius:xzn', 'theta:yzn', 'phi:zzn',
     &     (/1,m, n_s,n_e, o_s,o_e/), 
     &     (/1,m, 1,n, 1,o /))
#endif /*CFC_MHD */
    
      call close_data_file(fd)

      return

      END SUBROUTINE write_coconut_data

c     ==================================================================


      END MODULE gr_output

c     ==================================================================
