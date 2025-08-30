c     ==================================================================

      MODULE gr_alloc

c     ==================================================================

      IMPLICIT NONE

      CONTAINS



c     ==================================================================

      SUBROUTINE allocate_coconut_arrays

c     ==================================================================

      USE precision
      USE size_cfc, ONLY: m, n, o, n_s, n_e, o_s, o_e, n_loc, o_loc

      USE conserved_cfc
      USE enecons_cfc
      USE fluxes_cfc
      USE hydro_primitives_cfc
      USE shock_cfc
      USE interface_hydro_cfc
      USE metric_cfc
      USE nutra_cfc
      USE ppm_cfc
      USE sources_cfc
      USE temp_aux_cfc

#ifdef CFC_TRANSPORT2
      USE metric_aux_cfc
      USE mod_poisson_aux
#endif

#ifndef CFC_TRANSPORT2
      USE modpot_cfc
#endif

#ifdef MPI_HYDRO
      USE mo_mpi
#endif

      USE abort

      USE coconut, ONLY : qtot, xtot, sumdvy, sumdvy2
      USE gr_tstep, ONLY: allocate_gr_tstep
      USE legendre, ONLY: allocate_legendre
      USE mod_hilfsgroessen_r, ONLY: allocate_hilfsgroessen_r
      USE grid_cfc, ONLY: allocate_grid_cfc
      USE perm_aux_cfc, ONLY: allocate_perm_aux_cfc
      USE ppm_cfc, ONLY: allocate_ppm_cfc

      USE print_stdout_mod, only : print_memory_alloc

      USE configure
      USE meminfo

      IMPLICIT NONE

      real (kind=rk) :: mem_local, mem_global

      integer (kind=ik) :: istat

#ifdef CHECK_MEMORY
      call meminfo_start("GR allocate")
#endif

      mem_global = 0.0_rk

      ! allocate arrays which were originally static in subroutine
      ! avare_hy

#ifdef CFC_MHD
      allocate(qtot(1:m,11),xtot(1:m,config%qn), 
     &         sumdvy(1:m), sumdvy2(1:m), stat=istat)
#else /* CFC_MHD */      
      allocate(qtot(1:m,5),xtot(1:m,config%qn), 
     &         sumdvy(1:m), sumdvy2(1:m), stat=istat)
#endif /* CFC_MHD */
      
      if (istat .ne. 0) then
         raise_abort("allocation of module cococnut failed")
      end if

      mem_local = 8*(m*5+m*config%qn+m+m)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 1")

#ifdef MPI_HYDRO
      n_s = qy_s
      n_e = qy_e
      n_loc = n_e - n_s + 1
c      o_s = qz_s
c      o_e = qy_e
      o_s = qz_s
      o_e = qz_e
      o_loc = o_e - o_s + 1
#else
      n_s = 1
      n_e = config%qy
      n_loc = config%qy
      o_s = 1
      o_e = config%qz
      o_loc = config%qz
#endif

c     ------------------------------------------------------------------
c     conserved variables (module conserved_cfc)

      allocate (
     &     d_cap_hat (1 : m, n_s : n_e, o_s : o_e),
     &     s_1_hat   (1 : m, n_s : n_e, o_s : o_e),
     &     s_2_hat   (1 : m, n_s : n_e, o_s : o_e),
     &     s_3_hat   (1 : m, n_s : n_e, o_s : o_e),
     &     tau_hat   (1 : m, n_s : n_e, o_s : o_e),
     &     d_cap_xnu_hat (1 : m, n_s : n_e, o_s : o_e, 1 :  config%qn),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module conserved_cfc (1) failed")
      end if

      mem_local = 8*(m*(n_e-n_s+1)*(o_e-o_s+1)*(5+config%qn))
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 2")
      
      allocate (
     &     d_cap_hat_old (1 : m, n_s : n_e, o_s : o_e),
     &     s_1_hat_old   (1 : m, n_s : n_e, o_s : o_e),
     &     s_2_hat_old   (1 : m, n_s : n_e, o_s : o_e),
     &     s_3_hat_old   (1 : m, n_s : n_e, o_s : o_e),
     &     tau_hat_old   (1 : m, n_s : n_e, o_s : o_e),
     &     d_cap_xnu_hat_old (1 : m, n_s : n_e, o_s : o_e, 
     &                                              1 : config%qn),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module conserved_cfc (2) failed")
      end if

      mem_local = 8*(m*(n_e-n_s+1)*(o_e-o_s+1)*(5+config%qn))
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 3")


c     ------------------------------------------------------------------
c     old gravitational potential/lapse function (module enecons_cfc)

      allocate (
     &     phi_potential_alt (1 : m, n_s : n_e, o_s : o_e),
     &     alpha_alt (1 : m , n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module enecons_cfc failed")
      end if

      mem_local = 8*(m*(n_e-n_s+1)*(o_e-o_s+1)*2)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 4")

c     ------------------------------------------------------------------
c     interface fluxes & pressure (module fluxes_cfc)

      allocate (
#ifndef CFC_MHD      
     &     f_num_if (0 : m, n_s - 1 : n_e, o_s-1 : o_e, 1 : 5),
#else      
     &     f_num_if (0 : m, n_s - 1 : n_e, o_s-1 : o_e, 1 : 11),
#endif
     &     f_num_xnu_if (0 : m, n_s - 1 : n_e, o_s - 1 : o_e,
     &     1 : config%qn),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module fluxes_cfc failed")
      end if

      mem_local = 8*((m+1)*(n_e-n_s+2)*(o_e-o_s+2)*(5+config%qn))
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 5")

      allocate (
#ifndef CFC_MHD      
     &     f_diff (1 : m, n_s : n_e, o_s : o_e, 1 : 5),
#else
     &     f_diff (1 : m, n_s : n_e, o_s : o_e, 1 : 11),
#endif
     &     f_diff_xnu (1 : m, n_s : n_e, o_s : o_e, 1 : config%qn),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module fluxes_cfc failed")
      end if

      mem_local = 8*(m*(n_e-n_s+1)*(o_e-o_s+1)*(5+config%qn))
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 6")


      allocate (
     &     pav_r (0 : m, n_s : n_e, o_s : o_e),
     &     pav_theta (1 : m, n_s - 1 : n_e, o_s : o_e),
     &     pav_varphi (1 : m, n_s : n_e, o_s - 1 : o_e),
     &     stat = istat)

      mem_local = 8*((m+1)*(n_e-n_s+1)*(o_e-o_s+1)+
     &     m*(n_e-n_s+2)*(o_e-o_s+1)+
     &     m*(n_e-n_s+1)*(o_e-o_s+2))
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 7")

      if (istat .ne. 0) then
         raise_abort("allocation of module fluxes_cfc failed")
      end if


c     ------------------------------------------------------------------
c     hydro variables and shock switch (modules hydro_primitives_cfc and shock_cfc)

      allocate (
     &     rho (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     v_1 (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     v_2 (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     v_3 (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     eps (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
#ifdef CFC_MHD
     &     eps_aux (1 : m, n_s : n_e, o_s : o_e),
     &     eps_old (1 : m, n_s : n_e, o_s : o_e),
#endif      
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module hydro_primitives_cfc failed")
      end if


      mem_local = 8*((m+8)*(n_e-n_s+9)*(o_e-o_s+9)*5)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 8")


      allocate (
     &     p     (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     w     (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     h     (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     t     (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     gamm  (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     cpot    (1 : m, n_s : n_e, o_s : o_e, 4),
     &     entropy (1 : m, n_s : n_e, o_s : o_e),
     &     v_squared       (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     c_sound_squared (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     xnnu  (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4, 
     &                                               1 : config%qn),
     &     stat = istat)

      mem_local = 8*((m+8)*(n_e-n_s+9)*(o_e-o_s+9)*(7+config%qn)+
     &     m*(n_e-n_s+1)*(o_e-o_s+1)*5)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNuT 9")

      if (istat .ne. 0) then
         raise_abort("allocation of module hydro_primitives_cfc failed")
      end if

      allocate (lshock (0 : m+1, n_s-1 : n_e+1, o_s - 1 : o_e +1),
     &     stat = istat)
      
      if (istat .ne. 0) then
         raise_abort("allocation of module shock_cfc failed")
      end if

      mem_local = 8*((m+1)*(n_e-n_s+3)*(o_e-o_s+3))
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 10")


c     ------------------------------------------------------------------
c     hydro variables at cell interfaces (module interface_hydro_cfc)

      allocate (
     &    rho_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &    v_1_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &    v_2_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &    v_3_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &    eps_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &    xnu_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1,
     &                                                   1 : config%qn),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module interface_hydro_cfc failed")
      end if

      mem_local = 8*((m+4)*(n_e-n_s+5)*(o_e-o_s+5)*(5+config%qn))
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 11")


      allocate (
     &     p_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &     h_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &     gamm_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &     v_squared_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &     c_sound_squared_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2,
     &     0 : 1),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module interface_hydro_cfc failed")
      end if

      mem_local = 8*((m+4)*(n_e-n_s+5)*(o_e-o_s+5)*5)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 11")



#ifdef CFC_MHD
      allocate (
     &     b_1 (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     b_2 (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     b_3 (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     e_1 (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     e_2 (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     e_3 (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     psi (-3 : m+4, n_s-4 : n_e+4, o_s-4 : o_e+4),
     &     psi_old (1:m,  n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of primmitives (MHD) failed")
      end if

      mem_local = 8*((m+8)*(n_e-n_s+9)*(o_e-o_s+9)*(7+config%qn)+
     &     m*(n_e-n_s+1)*(o_e-o_s+1)*3)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNuT-MHD 1")


      allocate (
     &     b_cap_1_hat   (1 : m, n_s : n_e, o_s : o_e),
     &     b_cap_2_hat   (1 : m, n_s : n_e, o_s : o_e),
     &     b_cap_3_hat   (1 : m, n_s : n_e, o_s : o_e),
     &     b_cap_1_hat_old   (1 : m, n_s : n_e, o_s : o_e),
     &     b_cap_2_hat_old   (1 : m, n_s : n_e, o_s : o_e),
     &     b_cap_3_hat_old   (1 : m, n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of conserved MHD variables failed")
      end if

      mem_local = 8*(m*(n_e-n_s+1)*(o_e-o_s+1)*6)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNuT-MHD 2")

      allocate (
     &    b_1_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &    b_2_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &    b_3_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &    e_1_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &    e_2_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &    e_3_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &    psi_if (-1 : m+2, n_s-2 : n_e+2, o_s-2 : o_e+2, 0 : 1),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of MHD variables (interaces) failed")
      end if

      mem_local = 8*((m+4)*(n_e-n_s+5)*(o_e-o_s+5)*3)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNuT-MHD 3")

      allocate (
     &     f_diff_b (1 : m, n_s : n_e, o_s : o_e, 1 : 3),
     &     emf_face (0 : m+1, n_s-1 : n_e+1, o_s-1 : o_e+1, 2, 1 : 3),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of MHD flux variables failed")
      end if

      mem_local = 8*(m*(n_e-n_s+1)*(o_e-o_s+1)*3+
     &     (m+1)*(n_e-n_s+2)*(o_e-o_s+2)*3)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNuT-MHD 4")

#endif /* CFC_MHD */


c     ------------------------------------------------------------------
c     metric functions, derivatives, potentials etc. (modules
c     metric_cfc and modpot_cfc)

      allocate (
     &     phi       (0 : m+1, n_s-1 : n_e+1, o_s-1 : o_e+1),
     &     alpha     (0 : m+1, n_s-1 : n_e+1, o_s-1 : o_e+1),
     &     beta_up_1 (0 : m+1, n_s-1 : n_e+1, o_s-1 : o_e+1),
     &     beta_up_2 (0 : m+1, n_s-1 : n_e+1, o_s-1 : o_e+1),
     &     beta_up_3 (0 : m+1, n_s-1 : n_e+1, o_s-1 : o_e+1),
     &     g_11 (1 : m, n_s : n_e, o_s : o_e),
     &     g_22 (1 : m, n_s : n_e, o_s : o_e),
     &     g_33 (1 : m, n_s : n_e, o_s : o_e),
     &     g_up_11 (1 : m, n_s : n_e, o_s : o_e),
     &     g_up_22 (1 : m, n_s : n_e, o_s : o_e),
     &     g_up_33 (1 : m, n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module metric_cfc failed")
      end if

      mem_local = 8*((m+2)*(n_e-n_s+3)*(o_e-o_s+3)*5+
     &     m*(n_e-n_s+1)*(o_e-o_s+1)*6)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 12")


      allocate (
     &     phi_r       (0 : m, n_s : n_e, o_s : o_e),
     &     alpha_r     (0 : m, n_s : n_e, o_s : o_e),
     &     beta_up_1_r (0 : m, n_s : n_e, o_s : o_e),
     &     beta_up_2_r (0 : m, n_s : n_e, o_s : o_e),
     &     beta_up_3_r (0 : m, n_s : n_e, o_s : o_e),
     &     g_11_r (0 : m, n_s : n_e, o_s : o_e),
     &     g_22_r (0 : m, n_s : n_e, o_s : o_e),
     &     g_33_r (0 : m, n_s : n_e, o_s : o_e),
     &     g_up_11_r (0 : m, n_s : n_e, o_s : o_e),
     &     g_up_22_r (0 : m, n_s : n_e, o_s : o_e),
     &     g_up_33_r (0 : m, n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module metric_cfc failed")
      end if

      mem_local = 8*((m+1)*(n_e-n_s+1)*(o_e-o_s+1)*11)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 13")


      allocate (
     &     phi_theta       (1 : m, n_s-1 : n_e, o_s : o_e),
     &     alpha_theta     (1 : m, n_s-1 : n_e, o_s : o_e),
     &     beta_up_1_theta (1 : m, n_s-1 : n_e, o_s : o_e),
     &     beta_up_2_theta (1 : m, n_s-1 : n_e, o_s : o_e),
     &     beta_up_3_theta (1 : m, n_s-1 : n_e, o_s : o_e),
     &     g_11_theta (1 : m, n_s-1 : n_e, o_s : o_e),
     &     g_22_theta (1 : m, n_s-1 : n_e, o_s : o_e),
     &     g_33_theta (1 : m, n_s-1 : n_e, o_s : o_e),
     &     g_up_11_theta (1 : m, n_s-1 : n_e, o_s : o_e),
     &     g_up_22_theta (1 : m, n_s-1 : n_e, o_s : o_e),
     &     g_up_33_theta (1 : m, n_s-1 : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module metric_cfc failed")
      end if

      mem_local = 8*(m*(n_e-n_s+2)*(o_e-o_s+1)*11)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 14")


      allocate (
     &     phi_varphi       (1 : m, n_s : n_e, o_s-1 : o_e),
     &     alpha_varphi     (1 : m, n_s : n_e, o_s-1 : o_e),
     &     beta_up_1_varphi (1 : m, n_s : n_e, o_s-1 : o_e),
     &     beta_up_2_varphi (1 : m, n_s : n_e, o_s-1 : o_e),
     &     beta_up_3_varphi (1 : m, n_s : n_e, o_s-1 : o_e),
     &     g_11_varphi (1 : m, n_s : n_e, o_s-1 : o_e),
     &     g_22_varphi (1 : m, n_s : n_e, o_s-1 : o_e),
     &     g_33_varphi (1 : m, n_s : n_e, o_s-1 : o_e),
     &     g_up_11_varphi (1 : m, n_s : n_e, o_s-1 : o_e),
     &     g_up_22_varphi (1 : m, n_s : n_e, o_s-1 : o_e),
     &     g_up_33_varphi (1 : m, n_s : n_e, o_s-1 : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module metric_cfc failed")
      end if

      mem_local = 8*(m*(n_e-n_s+1)*(o_e-o_s+2)*11)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 15")


      allocate (
     &     k_cap_hat (1 : m, n_s : n_e, o_s : o_e, 1 : 3, 1 : 3),
     &     k_cap_hat_squared (1 : m, n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module metric_cfc failed")
      end if

      mem_local = 8*(m*(n_e-n_s+1)*(o_e-o_s+1)*10)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 16")


      allocate (
     &     sqrt_gamma (0 : m+1, n_s-1 : n_e+1, o_s-1 : o_e+1),
     &     sqrt_gamma_inv (0 : m+1, n_s-1 : n_e+1, o_s-1 : o_e+1),
     &     sqrt_gamma_r (0 : m, n_s : n_e, o_s : o_e),
     &     sqrt_gamma_theta (1 : m, n_s-1 : n_e, o_s : o_e),
     &     sqrt_gamma_varphi (1 : m, n_s : n_e, o_s-1 : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module metric_cfc failed")
      end if

      mem_local = 8*((m+2)*(n_e-n_s+3)*(o_e-o_s+3)*2+
     &     (m+1)*(n_e-n_s+1)*(o_e-o_s+1)+
     &     m*(n_e-n_s+2)*(o_e-o_s+1)+
     &     m*(n_e-n_s+1)*(o_e-o_s+1))
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 17")


      allocate (
     &     beta (0 : m + 1, n_s - 1 : n_e + 1,o_s - 1 : o + 1, 3),
     &     w_tmp (0 : m + 1, n_s-1 : n_e + 1, o_s-1 : o_e+1, 3),
     &     dw_tmp_dr (1 : m, n_s : n_e, o_s : o_e, 3),
     &     dw_tmp_dtheta (1 : m, n_s : n_e, o_s : o_e, 3),
     &     stat = istat)
      
      if (istat .ne. 0) then
         raise_abort("allocation of module metric_cfc failed")
      end if

      mem_local = 8*((m+2)*(n_e-n_s+3)*(o_e-o_s+3)*6+
     &     m*(n_e-n_s+1)*(o_e-o_s+1)*6)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 18")


      allocate (
     &     scalar_tmp      (1 : m, n_s : n_e, o_s : o_e),
     &     s_cap_star      (1 : m, n_s : n_e, o_s : o_e, 3),
     &     s_cap           (1 : m, n_s : n_e, o_s : o_e, 3),
     &     rho_hat_phi     (1 : m, n_s : n_e, o_s : o_e),
     &     rho_hat_phi_1   (1 : m, n_s : n_e, o_s : o_e),
     &     rho_hat_phi_2   (1 : m, n_s : n_e, o_s : o_e),
     &     rho_hat_alpha   (1 : m, n_s : n_e, o_s : o_e),
     &     rho_hat_alpha_1 (1 : m, n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module metric_cfc failed")
      end if

      mem_local = 8*(m*(n_e-n_s+3)*(o_e-o_s+3)*12)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 19")

#ifndef CFC_TRANSPORT2
      allocate (
     &     phi_potential (0 : m+1, n_s-1 : n_e+1, o_s-1 : o_e+1),
     &     dphi_potential_dr (1 : m, n_s : n_e, o_s : o_e),
     &     dphi_potential_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     dphi_potential_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     dphi_potential_alt_dr (1 : m, n_s : n_e, o_s : o_e),
     &     dphi_potential_alt_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     dphi_potential_alt_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module metric_cfc failed")
      end if

      allocate (rho_grav (1 : m, n_s : n_e, o_s : o_e),
     &     rho_tangential_1_spherical (1 : m + 1),
     &     rho_tangential_2_spherical (1 : m + 1),
     &     rho_mod_spherical (1 : m + 1),
     &     eps_spherical (1 : m + 1),
     &     p_spherical (1 : m + 1),
     &     h_spherical (1 : m + 1),
     &     v_1_spherical (1 : m + 1),
     &     pnu_spherical (1 : m + 1),
     &     m_spherical_newton (1 : m),
     &     m_spherical_newton_if (0 : m),
     &     m_spherical_tov (1 : m),
     &     m_spherical_tov_if (0 : m),
     &     potential_spherical_newton (0 : m + 1),
     &     potential_spherical_tov (0 : m + 1),
     &     delta_potential_spherical (0 : m + 1),
     &     stat = istat)


      if (istat .ne. 0) then
         raise_abort("allocation of module modpot_cfc failed")
      end if
#endif


c     ------------------------------------------------------------------
c     Components of neutrino stress-energy tensor on hydro grid
c     (module nutra_cfc)

      allocate (
     &     enu (1 : m, n_s : n_e, o_s : o_e),
     &     fnu (1 : m, n_s : n_e, o_s : o_e),
     &     pnu (1 : m, n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module nutra_cfc failed")
      end if

      mem_local = 8*(m*(n_e-n_s+3)*(o_e-o_s+3)*3)
      mem_local = mem_local/1024._rk/1024._rk
      mem_global = mem_global + mem_local
      call print_memory_alloc(mem_local, mem_global, "CoCoNut 22")

#ifdef CFC_TRANSPORT2
c     ------------------------------------------------------------------
c     auxiliary variables for Poisson solver (module mod_poisson_aux)

      allocate (
     &     beta_x (0 : m + 1, n_s - 2 : n_e + 2, o_s - 1 : o_e + 1),
     &     beta_y (0 : m + 1, n_s - 2 : n_e + 2, o_s - 1 : o_e + 1),
     &     beta_z (0 : m + 1, n_s - 2 : n_e + 2, o_s - 1 : o_e + 1),
     &     chi    (0 : m + 1, n_s - 2 : n_e + 2, o_s - 1 : o_e + 1),
     &     s_x   (1 : m, n_s : n_e, o_s : o_e),
     &     s_y   (1 : m, n_s : n_e, o_s : o_e),
     &     s_z   (1 : m, n_s : n_e, o_s : o_e),
     &     s_tmp (1 : m, n_s : n_e, o_s : o_e),
     &     stat=istat)

c      mem_local = 8*(m*(n_e-n_s+3)*(o_e-o_s+3)*3+
c      mem_local = mem_local/1024._rk/1024._rk
c      mem_global = mem_global + mem_local
c      call print_memory_alloc(mem_local, mem_global, "CoCoNut 23")
#endif /* CFC_TRANSPORT2 */

c     ------------------------------------------------------------------
c     Source terms and some more metric derivatives (module sources_cfc)


#ifdef CFC_MHD
      allocate (s_total (1 : m, n_s : n_e, o_s : o_e, 1 : 9),
     &     stat = istat)
#else
      allocate (s_total (1 : m, n_s : n_e, o_s : o_e, 1 : 5),
     &     stat = istat)
#endif /* CFC_MHD */

      if (istat .ne. 0) then
         raise_abort("allocation of module sources_cfc failed")
      end if

      allocate (
     &     dln_sqrt_gamma_dr (1 : m, n_s : n_e, o_s : o_e),
     &     dln_sqrt_gamma_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     dln_sqrt_gamma_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     beta_up_k_k (1 : m, n_s : n_e, o_s : o_e),
     &     dphi_dt_o_phi (1 : m, n_s : n_e, o_s : o_e),
     &     alpha_dalpha_dr (1 : m, n_s : n_e, o_s : o_e),
     &     alpha_dalpha_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     alpha_dalpha_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     dalpha_dr (1 : m, n_s : n_e, o_s : o_e),
     &     dalpha_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     dalpha_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module sources_cfc failed")
      end if

      allocate (
     &     dg_00_dr (1 : m, n_s : n_e, o_s : o_e),
     &     dg_11_dr (1 : m, n_s : n_e, o_s : o_e),
     &     dg_22_dr (1 : m, n_s : n_e, o_s : o_e),
     &     dg_33_dr (1 : m, n_s : n_e, o_s : o_e),
     &     dg_00_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     dg_11_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     dg_22_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     dg_33_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     dg_00_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     dg_11_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     dg_22_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     dg_33_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     dg_01_dr (1 : m, n_s : n_e, o_s : o_e),
     &     dg_02_dr (1 : m, n_s : n_e, o_s : o_e),
     &     dg_03_dr (1 : m, n_s : n_e, o_s : o_e),
     &     dg_01_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     dg_02_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     dg_03_dtheta (1 : m, n_s : n_e, o_s : o_e),
     &     dg_01_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     dg_02_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     dg_03_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     stat = istat)

      allocate (
     &     excrv (1 : m, n_s : n_e, o_s : o_e, 1 : 3, 1 : 3),
     &     excrv2 (1 : m, n_s : n_e, o_s : o_e),
     &     excrv_beta (1 : m, n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module sources_cfc failed")
      end if


c     ------------------------------------------------------------------
c     some auxiliary arrays (module temp_aux_cfc)

      allocate (
     &     rho_gw (1 : m, n_s : n_e, o_s : o_e),
     &     v_gw_r (1 : m, n_s : n_e, o_s : o_e),
     &     v_gw_theta (1 : m, n_s : n_e, o_s : o_e),
     &     v_gw_varphi (1 : m, n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module temp_aux_cfc failed")
      end if


#ifdef CFC_TRANSPORT2
c     ------------------------------------------------------------------
c     auxiliary variables for metric solver (module metric_aux_cfc)

      allocate (
     &     alpha_phi6 (0 : m + 1, n_s - 1 : n_e + 1, o_s - 1 : o_e + 1),
     &     da_phi_dr      (1 : m, n_s : n_e, o_s : o_e),
     &     da_phi_dtheta  (1 : m, n_s : n_e, o_s : o_e),
     &     da_phi_dvarphi (1 : m, n_s : n_e, o_s : o_e),
     &     drel_tmp (n_s : n_e, o_s : o_e),
     &     stat = istat)

      if (istat .ne. 0) then
         raise_abort("allocation of module metric_aux_cfc failed")
      end if
#endif /* CFC_TRANSPORT2 */

      call allocate_gr_tstep
      call allocate_legendre
      call allocate_hilfsgroessen_r
      call allocate_grid_cfc
      call allocate_perm_aux_cfc
      call allocate_ppm_cfc

#ifdef CHECK_MEMORY
      call meminfo_stop("GR allocate")
#endif

      print *, 'Allocation of arrays for module CoCoNuT successful.'

      return

      END SUBROUTINE allocate_coconut_arrays

c     ==================================================================



      END MODULE gr_alloc

c     ==================================================================
