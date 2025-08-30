module hydro_memory

  use precision

  implicit none

  private

  public allocate_hydro_memory, deallocate_hydro_memory, alloc_hyd_mem_glob

  real(kind=rk) :: alloc_hyd_mem_glob, alloc_hyd_mem_loc

contains

!>
!> \par This subroutine allocates one after the other all the arrays which
!>      are defined in modules and used for the hydro part of the code
!>
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information
!>   $Revision:$
!>   $Date:$
!>
!> \endverbatim
!>
subroutine allocate_hydro_memory

  use precision
  use abort

  use configure
  use vnew_hy
  use vold_hy
  use mesh_hy
  use grd_hy
  use hydro_hy
  use physcs_hy
  use reman_hy
  use intrp_hy
  use spez_hy

  use totare_hy
  use nutrio_hy
  use totgrq_hy
!  use arecon_hy
#ifdef CFC_TRANSPORT
  use rshlat_hy
#endif

  use  gfloat_hy

  use marker_hy

  use massio_hy

  use revsho_hy
  use bndinf_hy
  use ancient_hy
#ifdef CFC_TRANSPORT
  use ancient_cfc
#endif
#ifdef HTCL
  use htcl_hy
#endif
  use vnuw_hy

  use mapare_proc

#ifdef CHECK_MEMORY
  use meminfo
#endif

#ifdef CHECK_THREAD_AFFINITY
  use thread_affinity, only : init_thread_affinity
#endif

  use nusource_data
  use print_stdout_mod
  use hydro_areas_mod
  use netw_thiel

  implicit none

  integer(kind=ik) :: istat

  alloc_hyd_mem_glob = 0._rk
  alloc_hyd_mem_loc  = 0._rk

#ifdef CHECK_MEMORY
  if (meminfo_flag) call meminfo_start("Hydro allocate")
#endif

  call allocate_gfloat_hy(alloc_hyd_mem_glob)

! from mod_ppm.F90 module vnew_old
  call allocate_vnew_hy(alloc_hyd_mem_glob)

! from mod_ppm.F90 module vold_old
  call allocate_vold_hy(alloc_hyd_mem_glob)

! from mod_ppm.F module mesh_hy
  call allocate_mesh_hy(alloc_hyd_mem_glob)

! still missing gfloat_hy ! not bit-identical ????

! from mod_ppm.F module grd_hy
  call allocate_grd_hy(alloc_hyd_mem_glob)

! from mod_ppm.F module hydro_hy
  call allocate_hydro_hy(alloc_hyd_mem_glob)

! from mod_ppm.F module physcs_hy
  call allocate_physcs_hy(alloc_hyd_mem_glob)

! from mod_ppm.F module reman_hy
  call allocate_reman_hy(alloc_hyd_mem_glob)

! from mod_ppm.F module spez_hy
  call allocate_spez_hy(alloc_hyd_mem_glob)

! from mod_ppm.F module intrp_hy
  call allocate_intrp_hy(alloc_hyd_mem_glob)


!from mod_hydro.F90 module totare_hy
  call allocate_totare_hy(alloc_hyd_mem_glob)

! from mod_hydro.F90 module rshlat_hy
#ifdef CFC_TRANSPORT
  call allocate_rshlat_hy
#endif /* CFC_TRANSPORT */

! from mod_hydro.F90 module nutrio_hy
  call allocate_nutrio_hy(alloc_hyd_mem_glob)

! from mod_hydro.F90 totgrq_hy
  call allocate_totgrq_hy(alloc_hyd_mem_glob)

 ! from mod_hydro.F90 module arecon_hy
  call allocate_hydro_areas(alloc_hyd_mem_glob)

! from mod_hydro.F90 module marker_hy
  call allocate_marker_hy(alloc_hyd_mem_glob)


! from mod_hydro.F90 module massio_hy
  call allocate_massio_hy(alloc_hyd_mem_glob)

! from mod_hydro.F90 module revsho_hy
  call allocate_revsho_hy(alloc_hyd_mem_glob)

! from mod_hydro.F90 module bndinf_hy
  call allocate_bndinf_hy(alloc_hyd_mem_glob)

! still missing bndinf_hy


! from mod_hydro.F90 module ancient_hy
  call allocate_ancient_hy(alloc_hyd_mem_glob)

#ifdef CFC_TRANSPORT
! from mod_ancient_cfc.F90 module ancient_cfc
  call allocate_ancient_cfc !No global memory accounting yet!
#endif

#ifdef HTCL
! from mod_htcl.F90 module htcl_hy
  call allocate_htcl_hy(alloc_hyd_mem_glob)
#endif

! from mod_vnuw.F90 module vnuw_hy
  call allocate_vnuw_hy(alloc_hyd_mem_glob)

! from mapare.F90 mpi buffers
  call allocate_mapare(alloc_hyd_mem_glob)

  ! from mod_nusource.F90
  call allocate_nusource_data(alloc_hyd_mem_glob)


#if 0
  ! from mod_network_thiel.F90
  if (config%use_network) then
     call allocate_network(alloc_hyd_mem_glob)
  endif
#endif

#ifdef CHECK_MEMORY
  if (meminfo_flag) call meminfo_stop("Hydro allocate")
#endif

#ifdef CHECK_THREAD_AFFINITY
  call init_thread_affinity
#endif


  call printit_taskX(0,"----------------------------------------------------- ")
  call printit_taskX(0," Totally  allocated (HYDRO) [Mb]",alloc_hyd_mem_glob)
  call printit_taskX(0,"----------------------------------------------------- ")

end subroutine allocate_hydro_memory


!>
!> \par This subroutine deallocates one after the other all the arrays which
!>      are defined in modules and used for the hydro part of the code
!>
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 2122 $
!>   $Date: 2011-12-09 21:58:41 +0100 (Fri, 09 Dec 2011) $
!>
!> \endverbatim
!>
subroutine deallocate_hydro_memory

  use precision
  use abort

  use vnew_hy
  use vold_hy
  use mesh_hy
  use grd_hy
  use hydro_hy
  use physcs_hy
  use reman_hy
  use intrp_hy
  use spez_hy

  use totare_hy
  use nutrio_hy
  use totgrq_hy
!  use arecon_hy
#ifdef CFC_TRANSPORT
  use rshlat_hy
#endif

  use gfloat_hy

  use marker_hy

  use massio_hy

  use revsho_hy
  use bndinf_hy
  use ancient_hy
#ifdef CFC_TRANSPORT
  use ancient_cfc
#endif
#ifdef HTCL
  use htcl_hy
#endif
  use vnuw_hy

  use mapare_proc

#ifdef CHECK_MEMORY
  use meminfo
#endif

#ifdef CHECK_THREAD_AFFINITY
  use thread_affinity, only : cleanup_thread_affinity
#endif

  use hydro_areas_mod
  implicit none

  integer(kind=ik) :: istat


  call deallocate_gfloat_hy
! from mod_ppm.F90 module vnew_old
  call deallocate_vnew_hy

! from mod_ppm.F90 module vold_old
  call deallocate_vold_hy

! from mod_ppm.F module mesh_hy
  call deallocate_mesh_hy

! still missing gfloat_hy ! not bit-identical ????

! from mod_ppm.F module grd_hy
  call deallocate_grd_hy

! from mod_ppm.F module hydro_hy
  call deallocate_hydro_hy

! from mod_ppm.F module physcs_hy
  call deallocate_physcs_hy

! from mod_ppm.F module reman_hy
  call deallocate_reman_hy

! from mod_ppm.F module spez_hy
  call deallocate_spez_hy

! from mod_ppm.F module intrp_hy
  call deallocate_intrp_hy


!from mod_hydro.F90 module totare_hy
  call deallocate_totare_hy

! from mod_hydro.F90 module rshlat_hy
#ifdef CFC_TRANSPORT
  call deallocate_rshlat_hy
#endif /* CFC_TRANSPORT */

! from mod_hydro.F90 module nutrio_hy
  call deallocate_nutrio_hy

! from mod_hydro.F90 totgrq_hy
  call deallocate_totgrq_hy

! from mod_hydro.F90 module marker_hy
  call deallocate_marker_hy

! from mod_hydro.F90 module massio_hy
  call deallocate_massio_hy

! from mod_hydro.F90 module revsho_hy
  call deallocate_revsho_hy

! from mod_hydro.F90 module bndinf_hy
  call deallocate_bndinf_hy

! still missing bndinf_hy


! from mod_hydro.F90 module ancient_hy
  call deallocate_ancient_hy

#ifdef CFC_TRANSPORT
! from mod_ancient_cfc.F90 module ancient_cfc
  call deallocate_ancient_cfc
#endif

#ifdef HTCL
! from mod_htcl.F90 module htcl_hy
  call deallocate_htcl_hy
#endif

! from mod_vnuw.F90 module vnuw_hy
  call deallocate_vnuw_hy

! from mapare.F90 mpi buffers
  call deallocate_mapare

#ifdef CHECK_MEMORY
  call cleanup_meminfo
#endif


#ifdef CHECK_THREAD_AFFINITY
  call cleanup_thread_affinity
#endif


end subroutine deallocate_hydro_memory

end module hydro_memory
