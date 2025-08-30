!-------------------------------------------------------------------
!>
!> \par This module provides some variables needed for timing the run-time conditions of the VERTEX code
!>
!> \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module cputim
  use precision
  use mo_mpi

  implicit none
! LOCAL varibales that are nor in modules

  SAVE

  private 
  public timer, timer_val, allocate_cputim, deallocate_cputim, second_v, imebte, cputime_flag, initialize_subtiming_counters


  logical :: cputime_flag=.true.

  type timer_val
     ! timers are named in the form: calling subroutine -> subroutine called
     ! e.g. in transp() call himodyn => transp_himodyn

     real(kind=rk), allocatable, dimension(:)       :: hydro_total
     real(kind=rk), allocatable, dimension(:)       :: rady_hydrostep
     real(kind=rk), allocatable, dimension(:)       :: rady_hydrostep_children
     real(kind=rk), allocatable, dimension(:,:)     :: hydro_setinf
     real(kind=rk), allocatable, dimension(:,:)     :: hydro_setinf_children
     real(kind=rk), allocatable, dimension(:)       :: hydro_detectshare
     real(kind=rk), allocatable, dimension(:)       :: hydro_detectshare_children
     real(kind=rk), allocatable, dimension(:,:)     :: hydro_hydroare
     real(kind=rk), allocatable, dimension(:,:)     :: hydro_hydroare_children
     real(kind=rk), allocatable, dimension(:,:)     :: hydro_fluxcor
     real(kind=rk), allocatable, dimension(:,:)     :: hydro_fluxcor_children
     real(kind=rk), allocatable, dimension(:,:)     :: hydro_filare
     real(kind=rk), allocatable, dimension(:,:)     :: hydro_filare_children
     real(kind=rk), allocatable, dimension(:,:)     :: hydro_tstep
     real(kind=rk), allocatable, dimension(:,:)     :: hydro_tstep_children
     real(kind=rk), allocatable, dimension(:)       :: hydro_poisson
     real(kind=rk), allocatable, dimension(:)       :: hydro_poisson_children
     real(kind=rk), allocatable, dimension(:)       :: hydro_accelare
     real(kind=rk), allocatable, dimension(:)       :: hydro_accelare_children
     real(kind=rk), allocatable, dimension(:)       :: hydro_nusource
     real(kind=rk), allocatable, dimension(:)       :: hydro_nusource_children
     real(kind=rk), allocatable, dimension(:,:,:)   :: accelare_sweep
     real(kind=rk), allocatable, dimension(:,:,:)   :: accelare_sweep_children
     real(kind=rk), allocatable, dimension(:,:,:,:) :: sweep_hydrow
     real(kind=rk), allocatable, dimension(:,:,:,:) :: sweep_hydrow_children
     real(kind=rk), allocatable, dimension(:,:,:)   :: hydroare_sweep
     real(kind=rk), allocatable, dimension(:,:,:)   :: hydroare_sweep_children
     real(kind=rk), allocatable, dimension(:,:,:,:) :: sweep_accel
     real(kind=rk), allocatable, dimension(:,:,:,:) :: sweep_accel_children
     real(kind=rk), allocatable, dimension(:,:,:,:) :: sweep_eos3d
     real(kind=rk), allocatable, dimension(:,:,:,:) :: sweep_eos3d_children
     real(kind=rk), allocatable, dimension(:,:,:,:) :: sweep_burn
     real(kind=rk), allocatable, dimension(:,:,:,:) :: sweep_burn_children

     ! transport part
     real(kind=rk), allocatable, dimension(:)       :: rady_neutrino
     real(kind=rk), allocatable, dimension(:)       :: rady_neutrino_children
     real(kind=rk), allocatable, dimension(:)       :: neutrino_nusource
     real(kind=rk), allocatable, dimension(:)       :: neutrino_nusource_children
     real(kind=rk), allocatable, dimension(:)       :: neutrino_qterms
     real(kind=rk), allocatable, dimension(:)       :: neutrino_qterms_children
     real(kind=rk), allocatable, dimension(:)       :: neutrino_mapra2hyd
     real(kind=rk), allocatable, dimension(:)       :: neutrino_mapra2hyd_children

     real(kind=rk), allocatable, dimension(:)       :: neutrino_savare
     real(kind=rk), allocatable, dimension(:)       :: neutrino_savare_children
     real(kind=rk), allocatable, dimension(:)       :: qterms_avecoeffmass
     real(kind=rk), allocatable, dimension(:)       :: qterms_avecoeffmass_children
     real(kind=rk), allocatable, dimension(:)       :: qterms_rstsects
     real(kind=rk), allocatable, dimension(:)       :: qterms_rstsects_children
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_getsect
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_getsect_children
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_maphyd2ra
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_maphyd2ra_children
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_putsect
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_putsect_children
     real(kind=rk), allocatable, dimension(:)       :: qterms_saveddi
     real(kind=rk), allocatable, dimension(:)       :: qterms_saveddi_children
     real(kind=rk), allocatable, dimension(:)       :: qterms_advec
     real(kind=rk), allocatable, dimension(:)       :: qterms_advec_children
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_mapeddi
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_mapeddi_children
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_transp
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_transp_children
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_sourcet
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_sourcet_children
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_enemom
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_enemom_children
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_conserve1d
     real(kind=rk), allocatable, dimension(:,:)     :: qterms_conserve1d_children

     real(kind=rk), allocatable, dimension(:,:)     :: transp_himodyn
     real(kind=rk), allocatable, dimension(:,:)     :: transp_himodyn_children
     real(kind=rk), allocatable, dimension(:,:)     :: transp_flagspec
     real(kind=rk), allocatable, dimension(:,:)     :: transp_flagspec_children
     real(kind=rk), allocatable, dimension(:,:,:)   :: transp_matkoeff
     real(kind=rk), allocatable, dimension(:,:,:)   :: transp_matkoeff_children
     real(kind=rk), allocatable, dimension(:,:)     :: transp_arvisco
     real(kind=rk), allocatable, dimension(:,:)     :: transp_arvisco_children
     real(kind=rk), allocatable, dimension(:,:)     :: transp_eddfak
     real(kind=rk), allocatable, dimension(:,:)     :: transp_eddfak_children
     real(kind=rk), allocatable, dimension(:,:)     :: transp_momeq
     real(kind=rk), allocatable, dimension(:,:)     :: transp_momeq_children
     real(kind=rk), allocatable, dimension(:,:)     :: transp_momeq_mebte
     real(kind=rk), allocatable, dimension(:,:)     :: transp_momeq_mebte_children
     real(kind=rk), allocatable, dimension(:,:)     :: transp_mebteloop
     real(kind=rk), allocatable, dimension(:,:)     :: transp_fullproblem
     real(kind=rk), allocatable, dimension(:,:,:)   :: momeq_momeqsub
     real(kind=rk), allocatable, dimension(:,:,:)   :: momeq_momeqsub_children


     real(kind=rk), allocatable, dimension(:)       :: matkoeff_eos
     real(kind=rk), allocatable, dimension(:)       :: matkoeff_eos_children
     real(kind=rk), allocatable, dimension(:)       :: emislte
     real(kind=rk), allocatable, dimension(:)       :: emislte_children
     real(kind=rk), allocatable, dimension(:)       :: abscnucl
     real(kind=rk), allocatable, dimension(:)       :: abscnucl_children
     real(kind=rk), allocatable, dimension(:)       :: scatnucl
     real(kind=rk), allocatable, dimension(:)       :: scatnucl_children
     real(kind=rk), allocatable, dimension(:)       :: scations
     real(kind=rk), allocatable, dimension(:)       :: scations_children
     real(kind=rk), allocatable, dimension(:)       :: absonuci
     real(kind=rk), allocatable, dimension(:)       :: absonuci_children
     real(kind=rk), allocatable, dimension(:)       :: pairs
     real(kind=rk), allocatable, dimension(:)       :: pairs_children
     real(kind=rk), allocatable, dimension(:)       :: brems
     real(kind=rk), allocatable, dimension(:)       :: brems_children
     real(kind=rk), allocatable, dimension(:)       :: nunu
     real(kind=rk), allocatable, dimension(:)       :: nunu_children
     real(kind=rk), allocatable, dimension(:)       :: absonucl
     real(kind=rk), allocatable, dimension(:)       :: absonucl_children
     real(kind=rk), allocatable, dimension(:)       :: isn
     real(kind=rk), allocatable, dimension(:)       :: isn_children
     real(kind=rk), allocatable, dimension(:)       :: scatele
     real(kind=rk), allocatable, dimension(:)       :: scatele_children

     real(kind=rk), allocatable, dimension(:)       :: cycle_tot
     real(kind=rk), allocatable, dimension(:)       :: transp_tot

     real(kind=rk), allocatable, dimension(:)       :: transp_comm
     real(kind=rk), allocatable, dimension(:)       :: omp_par
     real(kind=rk), allocatable, dimension(:)       :: average_eddfac
     real(kind=rk), allocatable, dimension(:)       :: hydro_comm
     integer(kind=ik)                               :: hydro_sweep_mode
     integer(kind=ik)                               :: transp_call_mode
     integer(kind=ik)                               :: momeq_call_mode

  end type timer_val

  type(timer_val) :: timer


  integer(kind=ik), allocatable :: imebte(:,:) ! numer of mebte iterations

contains

!>
!> \par This subroutine allocates the arrays from module cputim
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine allocate_cputim
    use precision
    use abort

    use configure
    implicit none
    integer(kind=ik) :: istat


    allocate(timer%rady_hydrostep(2), timer%rady_hydrostep_children(2),    &   
             timer%transp_himodyn(2,3), timer%hydro_total(2),              &
              timer%transp_himodyn_children(2,3),                          &
             timer%hydro_setinf(2,2), timer%hydro_setinf_children(2,2),    &
             timer%hydro_detectshare(2),                                   &
             timer%hydro_detectshare_children(2),                          &
             timer%hydro_hydroare(2,3), timer%hydro_hydroare_children(2,3),&
             timer%hydro_fluxcor(2,3), timer%hydro_fluxcor_children(2,3),  &
             timer%hydro_filare(2,3),  timer%hydro_filare_children(2,3),   &
             timer%hydro_tstep(2,3),  timer%hydro_tstep_children(2,3),     &
             timer%hydro_poisson(2), timer%hydro_poisson_children(2),      &
             timer%hydro_accelare(2), timer%hydro_accelare_children(2),    &
             timer%hydro_nusource(2),  timer%hydro_nusource_children(2),   &
             timer%accelare_sweep(2,2,3),                                  &
             timer%accelare_sweep_children(2,2,3),                         &
             timer%sweep_hydrow(2,2,3,2),                                  &
             timer%sweep_hydrow_children(2,2,3,2),                         &
             timer%hydroare_sweep(2,2,3),                                  &
             timer%hydroare_sweep_children(2,2,3),                         &
             timer%sweep_accel(2,2,3,2),                                   &
             timer%sweep_accel_children(2,2,3,2),                          &
             timer%sweep_eos3d(2,2,3,2),                                   &
             timer%sweep_eos3d_children(2,2,3,2),                          &
             timer%sweep_burn(2,2,3,2),                                    &
             timer%sweep_burn_children(2,2,3,2),                           &
             timer%neutrino_nusource(2),                                   &
             timer%neutrino_nusource_children(2) ,                         &
             timer%neutrino_qterms(2), timer%neutrino_qterms_children(2),  &
             timer%neutrino_mapra2hyd(2),                                  &
             timer%neutrino_mapra2hyd_children(2),                         &
             timer%rady_neutrino(2),                                       &
             timer%rady_neutrino_children(2),                              &
             timer%neutrino_savare(2), timer%neutrino_savare_children(2),  &
             timer%qterms_avecoeffmass(2),                                 &
             timer%qterms_avecoeffmass_children(2),                        &
             timer%qterms_rstsects(2), timer%qterms_rstsects_children(2),  &  
             timer%qterms_getsect(2,3), timer%qterms_getsect_children(2,3),  &
             timer%qterms_maphyd2ra(2,3),                                    &
             timer%qterms_maphyd2ra_children(2,3),                           &
             timer%qterms_putsect(2,3), timer%qterms_putsect_children(2,3),  &  
             timer%qterms_saveddi(2),  timer%qterms_saveddi_children(2),     &
             timer%qterms_advec(2), timer%qterms_advec_children(2),          &
             timer%qterms_mapeddi(2,4), timer%qterms_mapeddi_children(2,4),  &
             timer%qterms_transp(2,3), timer%qterms_transp_children(2,3),    &
             timer%qterms_sourcet(2,3), timer%qterms_sourcet_children(2,3),  &
             timer%qterms_enemom(2,3), timer%qterms_enemom_children(2,3),    &
             timer%qterms_conserve1d(2,3),                                   &
             timer%qterms_conserve1d_children(2,3),                          &
             timer%transp_flagspec(2,3), timer%transp_flagspec_children(2,3),&
             timer%transp_matkoeff(2,13,3),                                  &
             timer%transp_matkoeff_children(2,13,3),                         &
             timer%transp_arvisco(2,3), timer%transp_arvisco_children(2,3),  &
             timer%transp_eddfak(2,3), timer%transp_eddfak_children(2,3),    &
             timer%transp_momeq(2,3), timer%transp_momeq_children(2,3),      &
             timer%momeq_momeqsub(2,3,2),                                    &
             timer%momeq_momeqsub_children(2,3,2),                           &
             timer%transp_momeq_mebte(2,3),                                  &
             timer%transp_momeq_mebte_children(2,3),                         &
             timer%transp_mebteloop(2,3),                                    &
             timer%transp_fullproblem(2,3),                                  &
             timer%matkoeff_eos(2), timer%matkoeff_eos_children(2),          &
             timer%emislte(2), timer%emislte_children(2),                    & 
             timer%abscnucl(2), timer%abscnucl_children(2),                  &
             timer%scatnucl(2), timer%scatnucl_children(2),                  & 
             timer%scations(2), timer%scations_children(2),                  &
             timer%absonuci(2), timer%absonuci_children(2),                  &
             timer%pairs(2), timer%pairs_children(2),                        &
             timer%brems(2), timer%brems_children(2),                        &
             timer%nunu(2), timer%nunu_children(2),                          &
             timer%absonucl(2), timer%absonucl_children(2),                  &
             timer%isn(2), timer%isn_children(2),                            &
             timer%scatele(2), timer%scatele_children(2),                    &
             timer%omp_par(2), timer%average_eddfac(2),                      &

stat = istat)


    if (use_mpi) then
       allocate(imebte(nymoms:nymome,nzmoms:nzmome), stat=istat)
    else
       allocate(imebte(config%nystrt:config%nymom,config%nztra), stat=istat)

    endif

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of modul cputim 1 failed")
    end if

    imebte(:,:)                           = 0
 


    allocate(timer%cycle_tot(2), timer%transp_tot(2), timer%transp_comm(2), &
             timer%hydro_comm(2) )
 


    timer%hydro_total                   = 0._rk
    timer%rady_hydrostep                = 0._rk
    timer%rady_hydrostep_children       = 0._rk
    timer%hydro_setinf(:,:)             = 0._rk
    timer%hydro_setinf_children(:,:)    = 0._rk
    timer%hydro_detectshare(:)          = 0._rk
    timer%hydro_detectshare_children(:) = 0._rk
    timer%hydro_hydroare(:,:)           = 0._rk
    timer%hydro_hydroare_children(:,:)  = 0._rk
    timer%hydro_fluxcor(:,:)            = 0._rk
    timer%hydro_fluxcor_children(:,:)   = 0._rk
    timer%hydro_filare(:,:)             = 0._rk
    timer%hydro_filare_children(:,:)    = 0._rk
    timer%hydro_tstep(:,:)              = 0._rk
    timer%hydro_tstep_children(:,:)     = 0._rk
    timer%hydro_poisson(:)              = 0._rk
    timer%hydro_poisson_children(:)     = 0._rk
    timer%hydro_accelare(:)             = 0._rk
    timer%hydro_accelare_children(:)    = 0._rk
    timer%hydro_nusource(:)              = 0._rk
    timer%hydro_nusource_children(:)     = 0._rk
    timer%accelare_sweep(:,:,:)          = 0._rk
    timer%accelare_sweep_children(:,:,:) = 0._rk
    timer%sweep_hydrow(:,:,:,:)          = 0._rk
    timer%sweep_hydrow_children(:,:,:,:) = 0._rk
    timer%hydroare_sweep(:,:,:)          = 0._rk
    timer%hydroare_sweep_children(:,:,:) = 0._rk
    timer%sweep_accel(:,:,:,:)           = 0._rk
    timer%sweep_accel_children(:,:,:,:)  = 0._rk
    timer%sweep_eos3d(:,:,:,:)           = 0._rk
    timer%sweep_eos3d_children(:,:,:,:)  = 0._rk
    timer%sweep_burn(:,:,:,:)            = 0._rk
    timer%sweep_burn_children(:,:,:,:)   = 0._rk

    ! transport part
    timer%rady_neutrino(:)               = 0._rk
    timer%rady_neutrino_children(:)      = 0._rk
    timer%neutrino_nusource(:)           = 0._rk
    timer%neutrino_nusource_children(:)  = 0._rk
    timer%neutrino_qterms(:)             = 0._rk
    timer%neutrino_qterms_children(:)    = 0._rk
    timer%neutrino_mapra2hyd(:)          = 0._rk
    timer%neutrino_mapra2hyd_children(:) = 0._rk
    timer%neutrino_savare(:)             = 0._rk
    timer%neutrino_savare_children(:)     = 0._rk            
    timer%qterms_avecoeffmass(:)          = 0._rk
    timer%qterms_avecoeffmass_children(:) = 0._rk
    timer%qterms_rstsects(:)              = 0._rk
    timer%qterms_rstsects_children(:)     = 0._rk
    timer%qterms_getsect(:,:)             = 0._rk
    timer%qterms_getsect_children(:,:)    = 0._rk
    timer%qterms_maphyd2ra(:,:)           = 0._rk
    timer%qterms_maphyd2ra_children(:,:)  = 0._rk
    timer%qterms_putsect(:,:)             = 0._rk
    timer%qterms_putsect_children(:,:)    = 0._rk
    timer%qterms_saveddi(:)               = 0._rk
    timer%qterms_saveddi_children(:)      = 0._rk
    timer%qterms_advec(:)                 = 0._rk
    timer%qterms_advec_children(:)        = 0._rk
    timer%qterms_mapeddi(:,:)             = 0._rk
    timer%qterms_mapeddi_children(:,:)    = 0._rk
    timer%qterms_transp(:,:)              = 0._rk
    timer%qterms_transp_children(:,:)     = 0._rk
    timer%qterms_sourcet(:,:)             = 0._rk
    timer%qterms_sourcet_children(:,:)    = 0._rk
    timer%qterms_enemom(:,:)              = 0._rk
    timer%qterms_enemom_children(:,:)     = 0._rk
    timer%qterms_conserve1d(:,:)          = 0._rk
    timer%qterms_conserve1d_children(:,:) = 0._rk
    timer%transp_himodyn(:,:)             = 0._rk
    timer%transp_himodyn_children(:,:)    = 0._rk
    timer%transp_flagspec(:,:)            = 0._rk
    timer%transp_flagspec_children(:,:)   = 0._rk
    timer%transp_matkoeff(:,:,:)          = 0._rk
    timer%transp_matkoeff_children(:,:,:) = 0._rk
    timer%transp_arvisco(:,:)             = 0._rk
    timer%transp_arvisco_children(:,:)    = 0._rk
    timer%transp_eddfak(:,:)              = 0._rk
    timer%transp_eddfak_children(:,:)     = 0._rk
    timer%transp_momeq(:,:)               = 0._rk
    timer%transp_momeq_children(:,:)      = 0._rk
    timer%momeq_momeqsub(:,:,:)           = 0._rk
    timer%momeq_momeqsub_children(:,:,:)  = 0._rk
    timer%transp_momeq_mebte(:,:)         = 0._rk
    timer%transp_momeq_mebte_children(:,:)= 0._rk
    timer%transp_mebteloop(:,:)           = 0._rk
    timer%transp_fullproblem(:,:)         = 0._rk

    timer%matkoeff_eos(:)                 = 0._rk
    timer%matkoeff_eos_children(:)        = 0._rk
    timer%emislte(:)                      = 0._rk
    timer%emislte_children(:)             = 0._rk
    timer%abscnucl(:)                     = 0._rk
    timer%abscnucl_children(:)            = 0._rk
    timer%scatnucl(:)                     = 0._rk
    timer%scatnucl_children(:)            = 0._rk
    timer%scations(:)                     = 0._rk
    timer%scations_children(:)            = 0._rk
    timer%absonuci(:)                     = 0._rk
    timer%absonuci_children(:)            = 0._rk
    timer%pairs(:)                        = 0._rk
    timer%pairs_children(:)               = 0._rk
    timer%brems(:)                        = 0._rk
    timer%brems_children(:)               = 0._rk
    timer%nunu(:)                         = 0._rk
    timer%nunu_children(:)                = 0._rk
    timer%absonucl(:)                     = 0._rk
    timer%absonucl_children(:)            = 0._rk
    timer%isn(:)                          = 0._rk
    timer%isn_children(:)                 = 0._rk
    timer%scatele(:)                      = 0._rk
    timer%scatele_children(:)             = 0._rk
 
    timer%omp_par(:)                      = 0._rk
    timer%cycle_tot(:)                    = 0._rk
    timer%transp_tot(:)                   = 0._rk
   
    timer%transp_comm(:)                  = 0._rk
    timer%hydro_comm(:)                   = 0._rk
    timer%average_eddfac(:)               = 0._rk
  end subroutine allocate_cputim

!>
!> \par This subroutine deallocates the arrays from module cputim
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine deallocate_cputim
    use precision
    use abort

    implicit none
    integer(kind=ik) :: istat


    deallocate(timer%rady_hydrostep, timer%rady_hydrostep_children, &
               timer%hydro_total, timer%hydro_setinf,               &
               timer%hydro_setinf_children, timer%hydro_detectshare,&
               timer%hydro_detectshare_children,                    &
               timer%hydro_hydroare, timer%hydro_hydroare_children, &
               timer%hydro_fluxcor_children, timer%hydro_filare,    &
               timer%hydro_filare_children, timer%hydro_tstep,      &
               timer%hydro_tstep_children, timer%hydro_poisson,     &
               timer%hydro_poisson_children, timer%hydro_accelare,  &
               timer%hydro_accelare_children, timer%hydro_nusource, &
               timer%hydro_nusource_children, timer%accelare_sweep, &
               timer%accelare_sweep_children, timer%sweep_hydrow,   &
               timer%sweep_hydrow_children, timer%hydroare_sweep,   &
               timer%hydroare_sweep_children, timer%sweep_accel,    &
               timer%sweep_accel_children, timer%sweep_eos3d,       &
               timer%sweep_eos3d_children, timer%sweep_burn,        &
               timer%sweep_burn_children,timer%neutrino_qterms,     &
               timer%neutrino_qterms_children,                      &
               timer%neutrino_mapra2hyd,                            &
               timer%neutrino_mapra2hyd_children,                   &
               timer%rady_neutrino, timer%rady_neutrino_children,   &
               timer%neutrino_savare, timer%neutrino_savare_children, &
               timer%qterms_avecoeffmass,                             &
               timer%qterms_avecoeffmass_children,                    &
               timer%qterms_rstsects, timer%qterms_rstsects_children, &
               timer%qterms_getsect, timer%qterms_getsect_children,   &
               timer%qterms_maphyd2ra, timer%qterms_maphyd2ra_children, &
               timer%qterms_putsect,timer%qterms_putsect_children,      &
               timer%qterms_saveddi, timer%qterms_saveddi_children,     &
               timer%qterms_advec,  timer%qterms_advec_children,        &
               timer%qterms_mapeddi, timer%qterms_mapeddi_children,     &
               timer%qterms_transp, timer%qterms_transp_children,       &
               timer%qterms_sourcet, timer%qterms_sourcet_children,     &
               timer%qterms_enemom, timer%qterms_enemom_children,       &
               timer%qterms_conserve1d,                                 &
               timer%qterms_conserve1d_children,                        &
               timer%transp_himodyn, timer%transp_himodyn_children,     &
               timer%transp_flagspec, timer%transp_flagspec_children,   &
               timer%transp_matkoeff, timer%transp_matkoeff_children,   &
               timer%transp_arvisco, timer%transp_arvisco_children,     &
               timer%transp_eddfak, timer%transp_eddfak_children,       &
               timer%transp_momeq, timer%transp_momeq_children,         &
               timer%momeq_momeqsub, timer%momeq_momeqsub_children,     &
               timer%transp_momeq_mebte,                                &
               timer%transp_momeq_mebte_children,                       &
               timer%transp_mebteloop, timer%transp_fullproblem,        &
               timer%matkoeff_eos, timer%matkoeff_eos_children,         &
               timer%abscnucl, timer%abscnucl_children,                 &
               timer%scatnucl, timer%scatnucl_children,                 &
               timer%scations, timer%scations_children,                 &
               timer%absonuci, timer%absonuci_children, timer%pairs,    &
               timer%pairs_children, timer%brems, timer%brems_children, &
               timer%nunu, timer%nunu_children, timer%isn,              &
               timer%isn_children, timer%scatele,                       &
               timer%scatele_children, timer%omp_par,                   &
               timer%average_eddfac)


    deallocate( imebte, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of modul cputim 1 failed")
    end if

    deallocate(timer%cycle_tot, timer%transp_tot,                      &
               timer%transp_comm, timer%hydro_comm,                    &
               timer%emislte, timer%absonucl, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of modul cputim 2 failed")
    end if

  end subroutine deallocate_cputim


!>
!> \par This subroutine calls machine-dependent the respective clock functions for meassuring the cpu and wall clock time
!>
!>  return the CPU (t(1)) and wallclock-time (t(2)) in [s]
!>
!>  \author M. Rampp
!>
!> \param t  time-array(2)
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine second_v(t)
    use precision
    use setup_c
    implicit none

    real(kind=rk), intent(out) :: t(2)    ! t(1) = cpu time
                                          ! t(2) = wall clock time
    real(kind=rk) :: cpu

    if (cputime_flag) then
       call cpu_time(cpu)
       t(1) = cpu
       t(2) = real(ms_since_epoch(), kind=rk) / 1e6_rk
    endif
    
  end subroutine second_v

  subroutine initialize_subtiming_counters
    use configure
    implicit none


    ! timers which are needed fot timing the hydro part
    timer%hydro_total(:)                 = 0._rk
    timer%rady_hydrostep(:)              = 0._rk
    timer%rady_hydrostep_children(:)     = 0._rk
    timer%transp_comm(:)                 = 0._rk
    timer%hydro_setinf(:,:)              = 0.0_rk
    timer%hydro_setinf_children(:,:)     = 0.0_rk  
    timer%hydro_detectshare(:)           = 0.0_rk
    timer%hydro_detectshare_children(:)  = 0.0_rk
    timer%hydro_hydroare(:,:)            = 0.0_rk  
    timer%hydro_hydroare_children(:,:)   = 0.0_rk 
    timer%hydro_fluxcor(:,:)             = 0.0_rk
    timer%hydro_fluxcor_children(:,:)    = 0.0_rk
    timer%hydro_filare(:,:)              = 0.0_rk
    timer%hydro_filare_children(:,:)     = 0.0_rk
    timer%hydro_tstep(:,:)               = 0.0_rk
    timer%hydro_tstep_children(:,:)      = 0.0_rk
    timer%hydro_poisson(:)               = 0.0_rk
    timer%hydro_poisson_children(:)      = 0.0_rk
    timer%hydro_accelare(:)              = 0.0_rk
    timer%hydro_accelare_children(:)     = 0.0_rk
    timer%hydro_nusource(:)              = 0.0_rk
    timer%hydro_nusource_children(:)     = 0.0_rk
    timer%accelare_sweep(:,:,:)          = 0._rk
    timer%accelare_sweep_children(:,:,:) = 0._rk
    timer%sweep_hydrow(:,:,:,:)          = 0._rk
    timer%sweep_hydrow_children(:,:,:,:) = 0._rk
    timer%hydroare_sweep(:,:,:)          = 0._rk
    timer%hydroare_sweep_children(:,:,:) = 0._rk
    timer%sweep_accel(:,:,:,:)           = 0._rk
    timer%sweep_accel_children(:,:,:,:)  = 0._rk
    timer%sweep_eos3d(:,:,:,:)           = 0._rk
    timer%sweep_eos3d_children(:,:,:,:)  = 0._rk
    timer%sweep_burn(:,:,:,:)            = 0._rk
    timer%sweep_burn_children(:,:,:,:)   = 0._rk

    ! transport part
    timer%rady_neutrino(:)               = 0.0_rk
    timer%rady_neutrino_children(:)      = 0.0_rk
    timer%neutrino_nusource(:)           = 0.0_rk
    timer%neutrino_nusource_children(:)  = 0.0_rk
    timer%neutrino_qterms(:)             = 0.0_rk
    timer%neutrino_qterms_children(:)    = 0.0_rk
    timer%neutrino_mapra2hyd(:)          = 0.0_rk
    timer%neutrino_mapra2hyd_children(:) = 0.0_rk
    timer%neutrino_savare(:)             = 0.0_rk
    timer%neutrino_savare_children(:)    = 0.0_rk
    timer%qterms_avecoeffmass(:)          = 0.0_rk
    timer%qterms_avecoeffmass_children(:) = 0.0_rk
    timer%qterms_rstsects(:)              = 0.0_rk
    timer%qterms_rstsects_children(:)     = 0.0_rk
    timer%qterms_getsect(:,:)             = 0.0_rk
    timer%qterms_getsect_children(:,:)    = 0.0_rk
    timer%qterms_maphyd2ra(:,:)           = 0.0_rk
    timer%qterms_maphyd2ra_children(:,:)  = 0.0_rk
    timer%qterms_putsect(:,:)             = 0.0_rk
    timer%qterms_putsect_children(:,:)    = 0.0_rk
    timer%qterms_saveddi(:)               = 0.0_rk
    timer%qterms_saveddi_children(:)      = 0.0_rk
    timer%qterms_advec(:)                 = 0.0_rk
    timer%qterms_advec_children(:)        = 0.0_rk
    timer%qterms_mapeddi(:,:)             = 0.0_rk
    timer%qterms_mapeddi_children(:,:)    = 0.0_rk
    timer%qterms_transp(:,:)              = 0.0_rk
    timer%qterms_transp_children(:,:)     = 0.0_rk
    timer%qterms_sourcet(:,:)             = 0.0_rk
    timer%qterms_sourcet_children(:,:)    = 0.0_rk
    timer%qterms_enemom(:,:)              = 0.0_rk
    timer%qterms_enemom_children(:,:)     = 0.0_rk
    timer%qterms_conserve1d(:,:)          = 0.0_rk
    timer%qterms_conserve1d_children(:,:) = 0.0_rk
    timer%transp_himodyn(:,:)             = 0._rk
    timer%transp_himodyn_children(:,:)    = 0._rk
    timer%transp_flagspec(:,:)            = 0._rk
    timer%transp_flagspec_children(:,:)   = 0._rk

    timer%transp_matkoeff(:,:,:)          = 0._rk
    timer%transp_matkoeff_children(:,:,:) = 0._rk
    timer%transp_arvisco(:,:)             = 0._rk
    timer%transp_arvisco_children(:,:)    = 0._rk
    timer%transp_eddfak(:,:)              = 0._rk
    timer%transp_eddfak_children(:,: )    = 0._rk
    timer%transp_momeq(:,:)               = 0._rk
    timer%transp_momeq_children(:,: )     = 0._rk
    timer%momeq_momeqsub(:,:,:)           = 0._rk
    timer%momeq_momeqsub_children(:,:,:)  = 0._rk
    timer%transp_momeq_mebte(:,:)         = 0._rk
    timer%transp_momeq_mebte_children(:,:)= 0._rk
    timer%transp_mebteloop(:,:)           = 0._rk
    timer%transp_fullproblem(:,:)         = 0._rk

    timer%matkoeff_eos(:)                 = 0._rk
    timer%matkoeff_eos_children(:)        = 0._rk
    timer%abscnucl(:)                     = 0._rk
    timer%abscnucl_children(:)            = 0._rk




    timer%transp_tot(:)         = 0._rk
    timer%transp_comm(:)        = 0._rk

    timer%hydro_comm(:)         = 0.0_rk   

    timer%omp_par(:)            = 0.0_rk

    ! transport timers 
    timer%transp_tot(:)          = 0.0_rk
    timer%transp_comm(:)         = 0.0_rk
    timer%average_eddfac(:)      = 0.0_rk
  end subroutine initialize_subtiming_counters

end module cputim
