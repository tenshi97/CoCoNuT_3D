module timings

  private
  public :: timings_out

  contains

!> \par prints the timing information
!>
!> \author M. Rampp
!>
!> \detail 
!> print the timing information of the first six cycles
!> for the transport and hydro part
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!> \endverbatim
subroutine timings_out
  use precision
  use cputim
  use intgrs_hy
  !use arecon_hy
  use param_rt

  use mo_mpi

  use abort

  use hydro_areas_mod
  use configure
  implicit none
  integer(kind=ik) :: i
  real(kind=rk), dimension(2) :: parcpu, parcpu1


#ifndef DEBUG_TIMINGS

  if (myproc .eq. 0) then
     print *," Summary of timings: "
     print '(" CPU/WC time for one cycle:                 ",2f8.2," s")',timer%cycle_tot
     if (config%p_ntr .ne. 0) then
        print '(" Ratio transport/total:                     " ,2f8.2)',timer%transp_tot/timer%cycle_tot
        print '(" time for total transport cycle:            " ,2f8.2," [s]")',timer%transp_tot
     endif
     print '(" time for total hydro cycles:               " ,2f8.2,    " [s]")',timer%hydro_total(:)
     print '(" + # of H. steps:                             " , i4         )',areas%ix_are(1,11)
     print '("   hydro communication:                     " ,2f8.2,    " [%]")',timer%hydro_comm/timer%cycle_tot*100._rk
     print '("   transport Communication:                 " ,2f8.2," [%]")',timer%transp_comm/timer%cycle_tot*100._rk
     print '("   OpenMP parallelisation:                  " ,2f8.2," [%]")',timer%omp_par/timer%cycle_tot*100._rk
     print *," "
     print *,"detailed report: "
     print *," "
     print *,"hydro part (in order of called routines):    "
     print '(" hydro_step / children:                      " ,2f8.2," / " ,2f8.2," [s]")', &
        timer%rady_hydrostep(:), timer%rady_hydrostep_children(:)
     parcpu(:) = timer%hydro_setinf(:,1) + timer%hydro_setinf(:,2)
     parcpu1(:) = timer%hydro_setinf_children(:,1) + timer%hydro_setinf_children(:,2) 
#ifndef CFC_TRANSPORT
     print '(" |-> setinf / children:                      " ,2f8.2," / " ,2f8.2," [s]")',parcpu,parcpu1
     print '(" |-> detect_sh_are / children:               " ,2f8.2," / " ,2f8.2," [s]")', &
        timer%hydro_detectshare, timer%hydro_detectshare_children
     parcpu(:) = timer%hydro_fluxcor(:,1) + timer%hydro_fluxcor(:,2)
     parcpu1(:) = timer%hydro_fluxcor_children(:,1) + timer%hydro_fluxcor_children(:,2)
     print '(" |-> fluxcor / children:                     " ,2f8.2," / " ,2f8.2," [s]")',parcpu,parcpu1
     parcpu(:) = timer%hydro_filare(:,1) + timer%hydro_filare(:,2)
     parcpu1(:) = timer%hydro_filare_children(:,1) + timer%hydro_filare_children(:,2)
     print '(" |-> filare / children:                      " ,2f8.2," / " ,2f8.2," [s]")',parcpu,parcpu1
     parcpu(:) = timer%hydro_tstep(:,1) + timer%hydro_tstep(:,2)
     parcpu1(:) = timer%hydro_tstep_children(:,1) + timer%hydro_tstep_children(:,2)
     print '(" |->Hydro_are 1 / children:                  " ,2f8.2," / "  ,2f8.2," [s]")', &
        timer%hydro_hydroare(:,1), timer%hydro_hydroare(:,1)
     print '("    |-> Sweeps 1 / children:                 " ,2f8.2," / "   ,2f8.2," [s]")', &
        timer%hydroare_sweep(2,:,1),timer%hydroare_sweep_children(2,:,1)
     print '("        |-> hydrow loop / children:          " ,2f8.2," / "   ,2f8.2," [s]")', &
        timer%sweep_hydrow(2,:,1,1),timer%sweep_hydrow_children(2,:,1,1)
     print '("        |-> accelx loop / children:          " ,2f8.2," / "   ,2f8.2," [s]")', &
        timer%sweep_accel(2,:,1,1),timer%sweep_accel_children(2,:,1,1)
     print '("        |-> burning / children:              " ,2f8.2," / "   ,2f8.2," [s]")', &
        timer%sweep_burn(2,:,1,1),timer%sweep_burn_children(2,:,1,1)
     print '("        |-> eos3d / children:                " ,2f8.2," / "   ,2f8.2," [s]")', &
        timer%sweep_eos3d(2,:,1,1),timer%sweep_eos3d_children(2,:,1,1)

     if (config%nsdim .gt. 1) then
        print '("     |-> Sweeps 2  / children:               " ,2f8.2," / "    ,2f8.2," [s]")', &
                timer%hydroare_sweep(2,:,2), timer%hydroare_sweep_children(2,:,2)
        print '("         |-> hydrow loop / children:         ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_hydrow(2,:,2,1),timer%sweep_hydrow_children(2,:,2,1)
        print '("         |-> accelx loop / children:         ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_accel(2,:,2,1),timer%sweep_accel_children(2,:,2,1)
        print '("         |-> burning / children:             ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_burn(2,:,2,1),timer%sweep_burn_children(2,:,2,1)
        print '("         |-> eos3d / children:               ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_eos3d(2,:,2,1),timer%sweep_eos3d_children(2,:,2,1)
     endif

     if (config%nsdim .gt. 2) then
        print '("     |-> Sweeps 3  / children:               " ,2f8.2," / "    ,2f8.2," [s]")', &
                timer%hydroare_sweep(2,:,3), timer%hydroare_sweep_children(2,:,3)
        print '("         |-> hydrow loop / children:         ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_hydrow(2,:,2,1),timer%sweep_hydrow_children(2,:,3,1)
        print '("         |-> accelx loop / children:         ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_accel(2,:,3,1),timer%sweep_accel_children(2,:,3,1)
        print '("         |-> burning / children:             ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_burn(2,:,3,1),timer%sweep_burn_children(2,:,3,1)
        print '("         |-> eos3d / children:               ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_eos3d(2,:,3,1),timer%sweep_eos3d_children(2,:,3,1)
     endif

     print '(" |-> tstep: / children:                      " ,2f8.2," / "  ,2f8.2," [s]")',parcpu,parcpu1

     print '(" |-> poisson / children:                     " ,2f8.2," / "  ,2f8.2," [s]")', & 
        timer%hydro_poisson, timer%hydro_poisson_children

     print '(" |-> accelare / children:                    " ,2f8.2," / "   ,2f8.2," [s]")', &
        timer%hydro_accelare, timer%hydro_accelare_children
     print '("     |-> Sweeps 1  / children:               " ,2f8.2," / "    ,2f8.2," [s]")', &
        timer%accelare_sweep(1,:,1), timer%accelare_sweep_children(1,:,1)
     print '("         |-> hydrow loop / children:         ",2f8.2," / "     ,2f8.2," [s]")', &
        timer%sweep_hydrow(1,:,1,2),timer%sweep_hydrow_children(1,:,1,2)
     print '("         |-> accelx loop / children:         ",2f8.2," / "     ,2f8.2," [s]")', &
        timer%sweep_accel(1,:,1,2),timer%sweep_accel_children(1,:,1,2)
     print '("         |-> burning / children:             ",2f8.2," / "     ,2f8.2," [s]")', &
        timer%sweep_burn(1,:,1,2),timer%sweep_burn_children(1,:,1,2)
     print '("         |-> eos3d / children:               ",2f8.2," / "     ,2f8.2," [s]")', &
        timer%sweep_eos3d(1,:,1,2),timer%sweep_eos3d_children(1,:,1,2)

     if (config%nsdim .gt. 1) then
        print '("     |-> Sweeps 2  / children:               " ,2f8.2," / "    ,2f8.2," [s]")', &
                timer%accelare_sweep(1,:,2), timer%accelare_sweep_children(1,:,2)
        print '("         |-> hydrow loop / children:         ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_hydrow(1,:,2,2),timer%sweep_hydrow_children(1,:,2,2)
        print '("         |-> accelx loop / children:         ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_accel(1,:,2,2),timer%sweep_accel_children(1,:,2,2)
        print '("         |-> burning / children:             ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_burn(1,:,2,2),timer%sweep_burn_children(1,:,2,2)
        print '("         |-> eos3d / children:               ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_eos3d(1,:,2,2),timer%sweep_eos3d_children(1,:,2,2)
     endif

     if (config%nsdim .gt. 2) then
        print '("     |-> Sweeps 3  / children:               " ,2f8.2," / "    ,2f8.2," [s]")', &
                timer%accelare_sweep(1,:,3), timer%accelare_sweep_children(1,:,3)
        print '("         |-> hydrow loop / children:         ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_hydrow(1,:,3,2),timer%sweep_hydrow_children(1,:,3,2)
        print '("         |-> accelx loop / children:         ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_accel(1,:,3,2),timer%sweep_accel_children(1,:,3,2)
        print '("         |-> burning / children:             ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_burn(1,:,3,2),timer%sweep_burn_children(1,:,3,2)
        print '("         |-> eos3d / children:               ",2f8.2," / "     ,2f8.2," [s]")', &
                timer%sweep_eos3d(1,:,3,2),timer%sweep_eos3d_children(1,:,3,2)
     endif
     print '(" |-> nusource / children:                    " ,2f8.2," / "   ,2f8.2," [s]")', &
                timer%hydro_nusource, timer%hydro_nusource_children
#endif /* CFC_TRANSPORT */

     if (config%p_ntr .gt. 0) then
        print *," "
        print *,"detailed report: "
        print *," "
        print *,"transport part (in order of called routines): "
        print '(" neutrino / children:                        " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%rady_neutrino(:), timer%rady_neutrino_children(:)
        print '(" |-> nusource / children:                    " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%neutrino_nusource(:), timer%neutrino_nusource_children(:)
        print '(" |-> savare / children:                      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%neutrino_savare(:), timer%neutrino_savare_children(:)
        print '(" |-> map_ra2hyd / children:                  " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%neutrino_mapra2hyd(:), timer%neutrino_mapra2hyd_children(:)
        print '(" |-> qterms / children:                      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%neutrino_qterms(:), timer%neutrino_qterms_children(:)
        print '("     |-> avecoeff_mass / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_avecoeffmass,timer%qterms_avecoeffmass_children
        print '("     |-> rst_sects:                          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_rstsects,timer%qterms_rstsects_children

     if (config%use_spherical_eddington_factor) then
        print '("     |-> Aver. Edd-Fact:                     " ,2f8.2," [s]")',timer%average_eddfac
        print '("         |-> get_sect / children:            " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_getsect(:,1),timer%qterms_getsect_children(:,1)
        print '("         |-> map_hyd2ra / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_maphyd2ra(:,1), timer%qterms_maphyd2ra_children(:,1)
        print '("         |-> transp / children:              " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_transp(:,1),timer%qterms_transp_children(:,1)
        print '("           |-> himodyn / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_himodyn(:,1),timer%transp_himodyn_children(:,1) 
        print '("           |-> flagspec / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_flagspec(:,1),timer%transp_flagspec_children(:,1) 
        print '("           |-> matkoeff / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,1,1),timer%transp_matkoeff_children(:,1,1) 
        print '("               |-> eos / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,3,1),timer%transp_matkoeff_children(:,3,1) 
        print '("               |-> emislte / children:       " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,4,1),timer%transp_matkoeff_children(:,4,1) 
        print '("               |-> absonucl / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,5,1),timer%transp_matkoeff_children(:,5,1) 
        print '("               |-> scatnucl / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,6,1),timer%transp_matkoeff_children(:,6,1) 
        print '("               |-> scations / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,7,1),timer%transp_matkoeff_children(:,7,1) 
        print '("               |-> absonuci / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,8,1),timer%transp_matkoeff_children(:,8,1)
        print '("               |-> brems / children:         " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,9,1),timer%transp_matkoeff_children(:,9,1)
        print '("               |-> nunu / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,10,1),timer%transp_matkoeff_children(:,10,1)
        print '("               |-> scatele / children:       " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,11,1),timer%transp_matkoeff_children(:,11,1)
        print '("               |-> abscnucl / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,12,1),timer%transp_matkoeff_children(:,12,1)
        print '("               |-> isn / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,13,1),timer%transp_matkoeff_children(:,13,1)
        print '("           |-> arvisco / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_arvisco(:,1),timer%transp_arvisco_children(:,1) 
        print '("           |-> mebte loop:                   " ,2f8.2," [s]")',timer%transp_mebteloop(:,1)
        if (use_mpi) then
           print '("                # iterations:                "       , i4   )',imebte(nymoms,nzmoms)-1
        else
           print '("                # iterations:                "       , i4   )',imebte(config%nystrt,1)-1
        endif
        print '("              |-> eddfak / children:         " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_eddfak(:,1),timer%transp_eddfak_children(:,1) 
        print '("              |-> momeq / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_momeq_mebte(:,1),timer%transp_momeq_mebte_children(:,1) 
        print '("                |-> momeqsub / children:     " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%momeq_momeqsub(:,1,1),timer%momeq_momeqsub_children(:,1,1) 
        print '("           |-> full problem:                 " ,2f8.2," [s]")',timer%transp_fullproblem(:,1)
        print '("              |-> momeq / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_momeq(:,1),timer%transp_momeq_children(:,1)
        print '("                |-> momeqsub / children:     " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%momeq_momeqsub(:,1,2),timer%momeq_momeqsub_children(:,1,2) 
     endif ! config%use_spherical_eddington_factor

     print '("         |-> put_sect / children:            " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_putsect(:,1),timer%qterms_putsect_children(:,1)
     print '("         |-> sav_eddi / children:            " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_saveddi(:),timer%qterms_saveddi_children(:)

     print '("     |-> Neutrino advection / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_advec(:),timer%qterms_advec_children(:) 
     if (.not.config%use_spherical_eddington_factor) then
        print '("     |-> without 1D-Edd fact       ")'
        print '("         |-> get_sect / children:            " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_getsect(:,2),timer%qterms_getsect_children(:,2)
        print '("         |-> map_hyd2ra / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_maphyd2ra(:,2),timer%qterms_maphyd2ra_children(:,2)
        print '("         |-> map_eddi / children:            " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_mapeddi(:,2),timer%qterms_mapeddi_children(:,2)
        print '("         |-> transp / children:              " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_transp(:,2),timer%qterms_transp_children(:,2)
        print '("           |-> himodyn / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_himodyn(:,2),timer%transp_himodyn_children(:,2) 
        print '("           |-> flagspec / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_flagspec(:,2),timer%transp_flagspec_children(:,2) 
        print '("           |-> matkoeff / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,1,2),timer%transp_matkoeff_children(:,1,2) 
        print '("               |-> eos / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,3,2),timer%transp_matkoeff_children(:,3,2) 
        print '("               |-> emislte / children:       " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,4,2),timer%transp_matkoeff_children(:,4,2) 
        print '("               |-> absonucl / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,5,2),timer%transp_matkoeff_children(:,5,2) 
        print '("               |-> scatnucl / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,6,2),timer%transp_matkoeff_children(:,6,2) 
        print '("               |-> scations / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,7,2),timer%transp_matkoeff_children(:,7,2) 
        print '("               |-> absonuci / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,8,2),timer%transp_matkoeff_children(:,8,2)
        print '("               |-> brems / children:         " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,9,2),timer%transp_matkoeff_children(:,9,2)
        print '("               |-> nunu / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,10,2),timer%transp_matkoeff_children(:,10,2)
        print '("               |-> scatele / children:       " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,11,2),timer%transp_matkoeff_children(:,11,2)
        print '("               |-> abscnucl / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,12,2),timer%transp_matkoeff_children(:,12,2)
        print '("               |-> isn / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,13,2),timer%transp_matkoeff_children(:,13,2)
        print '("           |-> arvisco / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_arvisco(:,2),timer%transp_arvisco_children(:,2) 
        print '("           |-> mebte loop:                   " ,2f8.2," [s]")',timer%transp_mebteloop(:,2)
        if (use_mpi) then
           print '("                # iterations:                "       , i4   )',imebte(nymoms,nzmoms)-1
        else
           print '("                # iterations:                "       , i4   )',imebte(config%nystrt,1)-1
        endif
        print '("              |-> eddfak / children:         " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_eddfak(:,2),timer%transp_eddfak_children(:,2) 
        print '("              |-> momeq / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_momeq_mebte(:,2),timer%transp_momeq_mebte_children(:,2)
        print '("                |-> momeqsub / children:     " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%momeq_momeqsub(:,2,1),timer%momeq_momeqsub_children(:,2,1) 
        print '("           |-> full problem:                 " ,2f8.2," [s]")',timer%transp_fullproblem(:,2)
        print '("              |-> momeq / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_momeq(:,2),timer%transp_momeq_children(:,2)
        print '("                |-> momeqsub / children:     " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%momeq_momeqsub(:,2,2),timer%momeq_momeqsub_children(:,2,2) 
        print '("         |-> put_sect / children:            " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_putsect(:,2),timer%qterms_putsect_children(:,2)
        print '("         |-> sourcet / children:             " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_sourcet(:,2),timer%qterms_sourcet_children(:,2)
        print '("         |-> enemom / children:              " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_enemom(:,2),timer%qterms_enemom_children(:,2)
        if (config%nymom .eq. 1) then
           print '("         |-> conserve_1d / children   :      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_conserve1d(:,2),timer%qterms_conserve1d_children(:,2)
        endif
     endif ! .not. sperhical_eddington_factor
     if (config%use_spherical_eddington_factor) then
        print '("     |-> with 1D-Edd-fact                    ")'
        print '("         |-> get_sect / children:            " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_getsect(:,3),timer%qterms_getsect_children(:,3)
        print '("         |-> map_hyd2ra / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_maphyd2ra(:,3),timer%qterms_maphyd2ra_children(:,3)
        print '("         |-> map_eddi / children:            " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_mapeddi(:,3),timer%qterms_mapeddi_children(:,3)
        print '("         |-> transp / children:              " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_transp(:,3), timer%qterms_transp_children(:,3)
        print '("           |-> himodyn / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_himodyn(:,3),timer%transp_himodyn_children(:,3) 
        print '("           |-> flagspec / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_flagspec(:,3),timer%transp_flagspec_children(:,3) 
        print '("           |-> matkoeff / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,1,3),timer%transp_matkoeff_children(:,1,3) 
        print '("               |-> eos / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,3,3),timer%transp_matkoeff_children(:,3,3) 
        print '("               |-> emislte / children:       " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,4,3),timer%transp_matkoeff_children(:,4,3) 
        print '("               |-> absonucl / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,5,3),timer%transp_matkoeff_children(:,5,3) 
        print '("               |-> scatnucl / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,6,3),timer%transp_matkoeff_children(:,6,3) 
        print '("               |-> scations / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,7,3),timer%transp_matkoeff_children(:,7,3) 
        print '("               |-> absonuci / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,8,3),timer%transp_matkoeff_children(:,8,3)
        print '("               |-> brems / children:         " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,9,3),timer%transp_matkoeff_children(:,9,3)
        print '("               |-> nunu / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,10,3),timer%transp_matkoeff_children(:,10,3)
        print '("               |-> scatele / children:       " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,11,3),timer%transp_matkoeff_children(:,11,3)
        print '("               |-> abscnucl / children:      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,12,3),timer%transp_matkoeff_children(:,12,3)
        print '("               |-> isn / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_matkoeff(:,13,3),timer%transp_matkoeff_children(:,13,3)
        print '("           |-> arvisco / children:           " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_arvisco(:,3),timer%transp_arvisco_children(:,3) 
        print '("           |-> mebte loop:                   " ,2f8.2," [s]")',timer%transp_mebteloop(:,3)
        if (use_mpi) then
           print '("                # iterations:                "       , i4   )',imebte(nymoms,nzmoms)-1
        else
           print '("                # iterations:                "       , i4   )',imebte(config%nystrt,1)-1
        endif
        print '("              |-> eddfak / children:         " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_eddfak(:,3),timer%transp_eddfak_children(:,3) 
        print '("              |-> momeq / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_momeq_mebte(:,3),timer%transp_momeq_mebte_children(:,3)
        print '("                |-> momeqsub / children:     " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%momeq_momeqsub(:,3,1),timer%momeq_momeqsub_children(:,3,1) 
        print '("           |-> full problem:                 " ,2f8.2," [s]")',timer%transp_fullproblem(:,3)
        print '("              |-> momeq / children:          " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%transp_momeq(:,3),timer%transp_momeq_children(:,3)
        print '("                |-> momeqsub / children:     " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%momeq_momeqsub(:,3,2),timer%momeq_momeqsub_children(:,3,2) 
        print '("         |-> put_sect / children:            " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_putsect(:,3),timer%qterms_putsect_children(:,3)
        print '("         |-> sourcet / children:             " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_sourcet(:,3),timer%qterms_sourcet_children(:,3)
        print '("         |-> enemom / children:              " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_enemom(:,3),timer%qterms_enemom_children(:,3)
        if (config%nymom .eq. 1) then
           print '("         |-> conserve_1d / children   :      " ,2f8.2," / " ,2f8.2," [s]")', &
                timer%qterms_conserve1d(:,3),timer%qterms_conserve1d_children(:,3)
        endif
     endif ! use_spherical


  endif ! p_ntr

!        print '("               |-> ME-BTE it:           " ,2f8.2," [s]")',timer%mebte_iter
!
!        print '("                   |-> BTE:              " ,2f8.2," [s]")',timer%bte
!        print '("                   |-> ME:               " ,2f8.2," [s]")',timer%me
!        print '("                       |-> Setup , large:  " ,2f8.2," [s]")',timer%mebte_setup(:,1,1)
!        print '("                           |-> matkoeff:     " ,2f8.2," [s]")',timer%mebte_setup(:,1,2)
!        print '("                              |-> eos:        " ,2f8.2," [s]")',timer%mebte_setup(:,1,3)
!        print '("                              |-> emislte:    " ,2f8.2," [s]")',timer%mebte_setup(:,1,4)
!        print '("                              |-> absonucl:   " ,2f8.2," [s]")',timer%mebte_setup(:,1,5)
!        print '("                              |-> scatnucl:   " ,2f8.2," [s]")',timer%mebte_setup(:,1,6)
!        print '("                              |-> scations:   " ,2f8.2," [s]")',timer%mebte_setup(:,1,7)
!        print '("                              |-> absonuci:   " ,2f8.2," [s]")',timer%mebte_setup(:,1,8)
!        print '("                              |-> brems:      " ,2f8.2," [s]")',timer%mebte_setup(:,1,9)
!        print '("                              |-> nunu:       " ,2f8.2," [s]")',timer%mebte_setup(:,1,10)
!        print '("                              |-> scatele:    " ,2f8.2," [s]")',timer%mebte_setup(:,1,11)
!        print '("                              |-> abscnucl:   " ,2f8.2," [s]")',timer%mebte_setup(:,1,12)
!        print '("                              |-> isn:        " ,2f8.2," [s]")',timer%mebte_setup(:,1,13)
!
!        print '("                   |-> Matrix, large:  " ,2f8.2," [s]")',timer%mebte_matrix(:,1,1)
!        if (use_mpi) then
!        print '("                     # it.:"     , f5.2       )', &
!                real(timer%num_mebte_matrix(1,nymoms,nzmoms),kind=rk) !/real(imebte(config%nystrt,1),kind=rk)! we want absolute 
                                                            ! numbers here; the average can be computed easily
!        else
!           print '("                 # it.:"     , f5.2       )', &
!                real(timer%num_mebte_matrix(1,config%nystrt,1),kind=rk) !/real(imebte(config%nystrt,1),kind=rk) ! we want absolute 
                                                            ! numbers here; the average can be computed easily
!        endif
!        print '("                      |-> libs:        " ,2f8.2," [s]")',timer%mebte_matrix(:,1,2)
! 
!        print '("                       |-> Setup,  small:  " ,2f8.2," [s]")',timer%mebte_setup(:,2,1)
!        print '("                           |-> matkoeff:     " ,2f8.2," [s]")',timer%mebte_setup(:,2,2)
!        print '("                              |-> eos:        " ,2f8.2," [s]")',timer%mebte_setup(:,2,3)
!        print '("                              |-> emislte:    " ,2f8.2," [s]")',timer%mebte_setup(:,2,4)
!        print '("                              |-> absonucl:   " ,2f8.2," [s]")',timer%mebte_setup(:,2,5)
!        print '("                              |-> scatnucl:   " ,2f8.2," [s]")',timer%mebte_setup(:,2,6)
!        print '("                              |-> scations:   " ,2f8.2," [s]")',timer%mebte_setup(:,2,7)
!        print '("                              |-> absonuci:   " ,2f8.2," [s]")',timer%mebte_setup(:,2,8)
!        print '("                              |-> brems:      " ,2f8.2," [s]")',timer%mebte_setup(:,2,9)
!        print '("                              |-> nunu:       " ,2f8.2," [s]")',timer%mebte_setup(:,2,10)
!        print '("                              |-> scatele:    " ,2f8.2," [s]")',timer%mebte_setup(:,2,11)
!        print '("                              |-> abscnucl:   " ,2f8.2," [s]")',timer%mebte_setup(:,2,12)
!        print '("                              |-> isn:        " ,2f8.2," [s]")',timer%mebte_setup(:,2,13)
!        
!        print '("                   |-> Matrix, small:  " ,2f8.2," [s]")',timer%mebte_matrix(:,2,1)
!        if (use_mpi) then
!        print '("                     # it.:"     , f5.2       )', &
!                   real(timer%num_mebte_matrix(2,nymoms,nzmoms),kind=rk) !/real(imebte(config%nystrt,1),kind=rk) ! we want absolute 
                                                            ! numbers here; the average can be computed easily
!        else
!        print '("                     # it.:"     , f5.2       )', &
!                real(timer%num_mebte_matrix(2,config%nystrt,1),kind=rk) !/real(imebte(config%nystrt,1),kind=rk) ! we want absolute 
                                                            ! numbers here; the average can be computed easily
!        endif
!        print '("                      |-> libs:        " ,2f8.2," [s]")',timer%mebte_matrix(:,2,2)
!        print '("               |-> momeq (excl.MEBTE): " ,2f8.2," [s]")',timer%full_system - timer%mebte_iter     
!        print '("               |-> ME, full :           " ,2f8.2," [s]")',timer%me_full
!        print '("                 |->  Setup,  large:  " ,2f8.2," [s]")',timer%full_system_setup(:,1,1)
!        print '("                      |-> matkoeff   :  " ,2f8.2," [s]")',timer%full_system_setup(:,1,2)
!        print '("                          |-> eos      :  " ,2f8.2," [s]")',timer%full_system_setup(:,1,3)
!        print '("                          |-> emislte:    " ,2f8.2," [s]")',timer%full_system_setup(:,1,4)
!        print '("                          |-> absonucl:   " ,2f8.2," [s]")',timer%full_system_setup(:,1,5)
!        print '("                          |-> scatnucl:   " ,2f8.2," [s]")',timer%full_system_setup(:,1,6)
!        print '("                          |-> scations:   " ,2f8.2," [s]")',timer%full_system_setup(:,1,7)
!        print '("                          |-> absonuci:   " ,2f8.2," [s]")',timer%full_system_setup(:,1,8)
!        print '("                          |-> brems:      " ,2f8.2," [s]")',timer%full_system_setup(:,1,9)
!        print '("                          |-> nunu:       " ,2f8.2," [s]")',timer%full_system_setup(:,1,10)
!        print '("                          |-> scatele:    " ,2f8.2," [s]")',timer%full_system_setup(:,1,11)
!        print '("                          |-> abscnucl:   " ,2f8.2," [s]")',timer%full_system_setup(:,1,12)
!        print '("                          |-> isn:        " ,2f8.2," [s]")',timer%full_system_setup(:,1,13)
!
!        print '("                |-> Matrix, large:  " ,2f8.2," [s]")',timer%full_system_matrix(:,1,1)
!        if (use_mpi) then
!        print '("                   # it.:"     , f5.2       )', &
!             real(SUM(timer%num_full_system_matrix(1,nymoms:nymome,nzmoms:nzmome)),kind=rk) / &
!               (real(nymom_proc,kind=rk)*real(nzmom_proc,kind=rk)) ! we want absolute values
!        else
!        print '("                   # it.:"     , f5.2       )', &
!             real(SUM(timer%num_full_system_matrix(1,1:config%nymom,1:config%nztra)),kind=rk) / &
!               (real(config%nymom,kind=rk)*real(config%nztra,kind=rk)) ! we want absolute values
!        endif
!
!        print '("                  |-> libs:        " ,2f8.2," [s]")',timer%full_system_matrix(:,1,2)
!        print '("      ++ Setup,  small:  " ,2f8.2," [s]")',timer%full_system_setup(:,2,1)
!        print '("         + matkoeff   :  " ,2f8.2," [s]")',timer%full_system_setup(:,2,2)
!        print '("           + eos      :  " ,2f8.2," [s]")',timer%full_system_setup(:,2,3)
!        print '("           + emislte:    " ,2f8.2," [s]")',timer%full_system_setup(:,2,4)
!        print '("           + absonucl:   " ,2f8.2," [s]")',timer%full_system_setup(:,2,5)
!        print '("           + scatnucl:   " ,2f8.2," [s]")',timer%full_system_setup(:,2,6)
!        print '("           + scations:   " ,2f8.2," [s]")',timer%full_system_setup(:,2,7)
!        print '("           + absonuci:   " ,2f8.2," [s]")',timer%full_system_setup(:,2,8)
!        print '("           + brems:      " ,2f8.2," [s]")',timer%full_system_setup(:,2,9)
!        print '("           + nunu:       " ,2f8.2," [s]")',timer%full_system_setup(:,2,10)
!        print '("           + scatele:    " ,2f8.2," [s]")',timer%full_system_setup(:,2,11)
!        print '("           + abscnucl:   " ,2f8.2," [s]")',timer%full_system_setup(:,2,12)
!        print '("           + isn:        " ,2f8.2," [s]")',timer%full_system_setup(:,2,13)
!
!
!        print '("      ++ Matrix, small:  " ,2f8.2," [s]")',timer%full_system_matrix(:,2,1)
!        if (use_mpi) then
!        print '("             # it.:"     , f5.2       )', &
!                   real(SUM(timer%num_full_system_matrix(2,nymoms:nymome,nzmoms:nzmome)),kind=rk) / &
!                     (real(nymom_proc,kind=rk)*real(nzmom_proc,kind=rk)) ! we want absolute values
!        else
!        print '("             # it.:"     , f5.2       )', &    
!                   real(SUM(timer%num_full_system_matrix(2,1:config%nymom,1:config%nztra)),kind=rk) / &
!                     (real(config%nymom,kind=rk)*real(config%nztra,kind=rk)) ! we want absolute values
!        endif
!
!        print '("         ++ libs:        " ,2f8.2," [s]")',timer%full_system_matrix(:,2,2)
!
!
!
!
!        print '("           |-> put_sect:         " ,2f8.2," [s]")',timer%transp_putsect(:,1)
!        print '("           |-> sav_eddi:         " ,2f8.2," [s]")',timer%transp_saveddi(:)
!        print '("           |-> put_sect:             " ,2f8.2," [s]")',timer%transp_putsect(:,3)
!        print '("           |-> sourcet:              " ,2f8.2," [s]")',timer%transp_sourcet(:,3)
!        print '("           |-> enemom:               " ,2f8.2," [s]")',timer%transp_enemom(:,3)
!        if (config%nymom .eq. 1) then
!        print '("           |-> conserve_1d:          " ,2f8.2," [s]")',timer%transp_conserve1d(:,3)
!        endif
!        endif
!
!        print '("   |-> nusource:               " ,2f8.2," [s]")',timer%transp_nusource
!        print '("   |-> map_ra2hyd:             " ,2f8.2," [s]")',timer%transp_mapra2hyd
!        print '("   |-> savare:                 " ,2f8.2," [s]")',timer%transp_savare
!
!
!
!
!
!
!
!        
!     endif ! p_ntr


endif ! myproc


#endif /* DEBUG_TIMINGS */

end subroutine timings_out

end module timings
