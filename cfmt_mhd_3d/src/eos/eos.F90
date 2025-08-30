!>
!> \verbatim
!> this module provides the EoS interface for the VERTEX code
!>
!>  Author: M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module eos_sn2
  use precision
  use abort

  implicit none
! LOCAL variables that are not in modules
  private
  public init_eos,            &
         eos,                 &
         eos_network_wrapper, &
#ifndef ONEMG_EOS
         lsrolo
#else
         lsrolo, lsrolo2
#endif         

  real(kind=rk), save :: lsrolo=-1._rk
#ifdef ONEMG_EOS
  real(kind=rk), save :: lsrolo2=1._rk
#endif

contains
!> \verbatim
!>
!>  Initialize the EoS
!>
!> Author : M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine init_eos
    use precision
    use abort
    use nucparam
    use eos_data_type, only : nullify_pointer, zero_data_variables, &
                              high_den_eos, intermediate_eos,       &
                              highden_table,       &
                              intermediate_eos, lept_radi_table

    use data_highden_eos
    use data_wolffeos_low
    use lepton_radiation_low_den_eos, only : lepton_radiation_eos_init
    use data_low_density_nse

#ifdef ONEMG_EOS
      USE burn_lim, ONLY: thigh2
#endif

!#ifdef PROGRAM_remap
!    use restart, only : use_temp_restart
!#endif

    use mpi_vertex, only : myproc
    use print_stdout_mod


!    use hydro_memory, only : allocate_hydro_memory


    implicit none
    real(kind=rk)    :: global_memory
    integer(kind=ik) :: eos_number, i

    logical, save :: first_call = .true.
#ifdef PROGRAM_remap
    if (config%use_temp_restart) then
#endif

       call nullify_pointer
       call zero_data_variables

       if (high_den_eos .eq. 1) then
#ifndef USE_MEMORY_MAPPING
! L&S-EoS table is used for high density regime
          call load_eos_table    ('./tables/eos_ls_220.i3e',1, 50 ,1,50 ,30, 50)      ! filename of L&S-EoS-table
#else /* USE_MEMORY_MAPPING */

#ifdef LINUX
          call load_eos_table    ('./tables/eos_ls.i3e_little_endian_units',1,50,1,50,30,50)
#else
          raise_abort("not yet tested")
#endif /* LINUX */

#endif /* USE_MEMORY_MAPPING */

       endif

#ifdef SI_EOS
    ! SI-EoS table is used for high density regime
       call load_eos_table    ('./tables/eos_si.i3e',1,50,1,50,30,50)      ! filename of Si-EoS-table
#endif

#ifdef SKA_EOS
    ! Ska-EoS table is used for high density regime
       call load_eos_table    ('./tables/eos_ska.i3e',1,50,1,50,30,50)      ! filename of Ska-EoS-table
#endif

#ifdef SKM_EOS
    ! SkM*-EoS table is used for high density regime
       call  load_eos_table   ('./tables/eos_skm.i3e',1,50,1,50,30,50)      ! filename of SkM*-EoS-table
#endif

#ifdef SIII_EOS
    ! SIII-EoS table is used for high density regime
       call load_eos_table    ('./tables/eos_siii.i3e',1,50,1,50,30,50)     ! filename of SIII-EoS-table
#endif

#ifdef SFHO_EOS
    ! SFHO-EoS table is used for high density regime
       call load_eos_table    ('./tables/eos_sfho.i3e',1,50,1,50,30,50)      ! filename of SFHO-EoS-table
#endif

#ifdef SFHX_EOS
    ! SFHO-EoS table is used for high density regime
       call load_eos_table    ('./tables/eos_sfhx.i3e',1,50,1,50,30,50)      ! filename of SFHX-EoS-table
#endif
       
       if (high_den_eos .eq. 2) then
       ! WOLFF-EoS table is used for high density regime
          call load_eos_table ('./tables/eos_wolff.i3e',1,50,1,50,30,50)    ! filename of WOLFF-EoS-table
       endif

       if (high_den_eos .eq. 3) then
       ! SHEN-EoS table is used fo high density regime
          call load_eos_table ('./tables/eos_shen.i3e',1,50,1,50,30,50)     ! filename of SHEN-EoS-table
       endif

       if (intermediate_eos .eq. 4) then
          call LOADTB_WOLFFEOS_LOW('./tables/eos_wolfflow.i3e')
       endif

    ! check whether an HIGH denity EoS was specified
       if (high_den_eos .eq. 0) then
          raise_abort("eos.F90(): You did not specify an high density EoS")
       endif

#if !defined(LATTIMER_EOS) && !defined(SI_EOS) && !defined(SKA_EOS) && !defined(SKM_EOS) && !defined(SIII_EOS)  && !defined(WOLFF_EOS) && !defined(SHEN_EOS) && !defined(SFHO_EOS) && !defined(SFHX_EOS) && defined(PROGRAM_vertex)
#error ================================================
#error  YOU MUST SPECIFY AN EoS FOR HIGH DENSITY REGIME
#error ================================================
#endif /* no eos defined */

       call lepton_radiation_eos_init('./tables/eos_tj.i3e')      ! filename of ThJ-EoS-table

#ifdef PROGRAM_remap
    endif ! config%use_temp_restart
#endif
    global_memory = 0.0_rk 
    call get_nuc_indices(global_memory)

#ifdef PROGRAM_remap
    if (config%use_temp_restart) then
#endif
       if (config%low_den_nse_eos .eq. 4) then
          call init_low_density_nse('./tables/eos_kk.i3e')
       endif

       if (config%low_den_nse_eos .eq. 17) then
          call init_low_density_nse('./tables/eostab_kok_nse_17.i3e')
       endif

       if (config%low_den_nse_eos .eq. 22) then
          call init_low_density_nse('./tables/eostab_kok_nse_22.i3e')
       endif

       if (config%low_den_nse_eos .eq. 23) then
          call init_low_density_nse('./tables/eostab_kok_nse_23.i3e')
       endif

       ! initialize eos boundaries between
       ! highdensity - (wolff lowdensity) - lowdensity EoS

       ! determine eos_number stop if more than one EoS is loaded
       ! (Wofflow_eos does not count here)
       eos_number=0

       do i=1, size(highden_table(:))
          if ((i .ne. 4) .and. (i .eq. highden_table(i)%eos_number)) then
             if (eos_number .ne. 0) then
                raise_abort("eos.f90(): more than one eos-table loaded")
             endif
             eos_number=i
          endif
       enddo

       if (intermediate_eos .eq. 4 .and. config%tj_ls_rho .eq. -1._rk) then
          ! use wolff-low-den eos but with predefined transition densities
          ! corresponds to the Vertex-case: intermediate_eos = 4 and 
          ! config%tj_wolff_rho == -1

          config%tj_wolff_rho=-1._rk
       endif

       if (intermediate_eos .eq. 4 .and. config%tj_ls_rho .ne. -1._rk) then
          ! use wolff-low-den eos and with choosen transition densities
          ! corresponds to the Vertex-case: intermediate_eos = 4 and 
          ! tj_wolff_rho == value and wolff_ls_rho=value

          raise_abort("THIS EOS CASE DOES NOT WORK")

          config%tj_wolff_rho=config%tj_ls_rho
          config%wolff_ls_rho=config%tj_ls_rho
       endif

       if (intermediate_eos .ne. 4 .and. config%tj_ls_rho .eq. -1._rk) then
          ! do not use wolff-low-den eos and use predefined value
          ! for transition between high-density eos and low-density eos
          ! corresponds to the Vertex-case: intermediate_eos != 4 and 
          ! tj_ls_rho == -1. and wolff_ls_rho=value

          config%tj_ls_rho=-1._rk
       endif

       if (intermediate_eos .ne. 4 .and. config%tj_ls_rho .ne. -1._rk) then
          ! do not use wolff-low-den eos and use chosen value
          ! for transition between high-density eos and low-density eos
          ! corresponds to the Vertex-case: intermediate_eos != 4 and 
          ! tj_ls_rho == value and wolff_ls_rho=value

       endif


       call printit_taskX(0,"EoS Table partitioning: ")
       call printit_taskX(0,"This partitioning might be a change")
       call printit_taskX(0,"compared to a previous start")

       if (config%eos_sw .ne. 0 .and. .not.(first_call)) then
          ! if eos_sw == 1 then you want to SWITCH lsrolo during the run-time
          ! _AFTER_ the first timestep to another value in an energy conserving way
          ! i.e. in the first timestep lsrolo etc. should be tj_ls_rho in the SECOND
          ! timestep lsrolo = tj_ls_rho
! 
          if (intermediate_eos .ne. 4) then
             highden_table(eos_number)%romin = config%tj_ls_rho
             lept_radi_table(1)%romax = config%tj_ls_rho
             lsrolo = config%tj_ls_rho

          else
             raise_abort("Not yet implemented")

             highden_table(eos_number)%romin = config%tj_wolff_rho
             lept_radi_table(1)%romax = config%tj_wolff_rho
             lsrolo  = config%tj_wolff_rho
          endif
       endif

       if (first_call) then
          call printit_taskX(0,"EoS Table partitioning: ")

          if (intermediate_eos .eq. 4) then
             if (config%tj_wolff_rho == -1) then
              ! this is the standard setup for the High-Density EOS
                ! SWromin is taken from the EoS-table

                highden_table(eos_number)%romin = 3.e9_rk
                highden_table(4)%romax = 3.e9_rk 

                lsrolo            = 9.e7_rk
                lept_radi_table(1)%romax           = lsrolo
                highden_table(4)%romin = lsrolo
                config%tj_wolff_rho      = 9.e7_rk
                config%wolff_ls_rho      = 3.e9_rk

             else
                highden_table(eos_number)%romin = config%wolff_ls_rho
                highden_table(4)%romax = config%wolff_ls_rho

                highden_table(4)%romin = config%tj_wolff_rho
                lept_radi_table(1)%romax     = config%tj_wolff_rho
                lsrolo      = config%tj_wolff_rho
             endif ! tj_wolff_rh == -1

             call printit_taskX(0,"TJ:                   rho < ", config%tj_wolff_rho)
             call printit_taskX(0,highden_table(4)%name, config%tj_wolff_rho)
             call printit_taskX(0,"                   <= rho < ", config%wolff_ls_rho)
             call printit_taskX(0,highden_table(eos_number)%name, config%wolff_ls_rho)
             call printit_taskX(0,"                   <= rho < ", highden_table(eos_number)%romax)

          else ! intermediate_eos != 4

#ifndef ONEMG_EOS
             ! Standard EoS setup
             if (config%tj_ls_rho == -1) then
                ! this is the standard setup for the High-Density EOS
                ! lsromin is taken from the EoS-table
                lsrolo=1.07 * 10**highden_table(eos_number)%lro_hd(2)
             endif

             if (config%tj_ls_rho /= -1 .and. config%eos_sw .eq. 1) then

                ! Careful: if tj_ls_rho is set, and we want to switch
                ! the EOS boundary, i.e. eos_sw = 1, then you have
                ! to start in the FIRST timestep nevertheless with
                ! the original value !!

                lsrolo=1.07 * 10**highden_table(eos_number)%lro_hd(2)
             endif

             if (config%tj_ls_rho /= -1 .and. config%eos_sw .eq. 0) then

                ! Careful: if tj_ls_rho is set, and we want to switch
                ! the EOS boundary, i.e. eos_sw = 1, then you have
                ! to start in the FIRST timestep nevertheless with
                ! the original value !!

                lsrolo=config%tj_ls_rho
             endif

             highden_table(eos_number)%romin = lsrolo
#ifdef DYNAMIC_EOS
             highden_table(eos_number)%romin_lc = lsrolo
#endif
             lept_radi_table(1)%romax        = lsrolo


             if (myproc .eq. 0) then
                write(*,'(a,1pe12.4)')           "TJ:                 rho < ", lept_radi_table(1)%romax
                write(*,'(a,1pe12.4,a,1pe12.4)') "LS: ", lsrolo, &
                                " <= rho < ", highden_table(eos_number)%romax
             endif


#else /* ONEMG_EOS */
     ! Special case for ONeMg progenitors
     !
     !          T
     !          ^         |
     !          |         |  high. den. EoS
     !          |         |
     ! tem_cut >|         \-----\
     !          |               |
     !          | low. den. EoS |
     !          |               |
     !          0---------------|-----------> rho
     !                    ^     ^           
     !                    |     |_ romin_high_hd
     !                    |_ romin_low_hd

             lsrolo=1.07_rk * 10.0_rk**highden_table(eos_number)%lro_hd(2)
             lsrolo2=1.0e11_rk
             config%tj_ls_rho=lsrolo
             highden_table(eos_number)%romin_low_hd  = config%tj_ls_rho
             highden_table(eos_number)%romin_high_hd = lsrolo2
             highden_table(eos_number)%tem_cut       = thigh2

             lept_radi_table(1)%romax_low_ld  = highden_table(eos_number)%romin_low_hd
             lept_radi_table(1)%romax_high_ld = highden_table(eos_number)%romin_high_hd
             lept_radi_table(1)%tem_cut       = highden_table(eos_number)%tem_cut
             if (myproc .eq. 0) then
                write(*,'(a,f5.3,a)') "For T >= ", highden_table(eos_number)%tem_cut * pc_kmev, " MeV:"
                write(*,'(a,1pe12.4)')           "TJ:                 rho < ", lept_radi_table(1)%romax_low_ld
                write(*,'(a,1pe12.4,a,1pe12.4)') "LS: ", highden_table(eos_number)%romin_low_hd, &
                  " <= rho < ", highden_table(eos_number)%romax

                write(*,'(a,f5.3,a)') "For T < ", highden_table(eos_number)%tem_cut * pc_kmev, " MeV:"
                write(*,'(a,1pe12.4)')           "TJ:                 rho < ", lept_radi_table(1)%romax_high_ld
                write(*,'(a,1pe12.4,a,1pe12.4)') "LS: ", highden_table(eos_number)%romin_high_hd, &
                  " <= rho < ", highden_table(eos_number)%romax
             endif


#endif /* ONEMG_EOS */

          endif ! intermediate-eos == 4

          if (config%eos_sw .ne. 0) then

             call printit_taskX(0," ")
             if (intermediate_eos .eq. 4) then

                call printit_taskX(0,"LSROLO,TJROMAX,and WOLFF-EOS-LOW-ROMIN will be changed to ", config%tj_wolff_rho)
                call printit_taskX(0,"in the second timestep")
             else
                call printit_taskX(0,"LSROLO,TJROMAX,and High-density-ROMIN will be changed to ",config%tj_ls_rho)
                call printit_taskX(0,"in the second timestep")
             endif ! intermediate eos eq 4
             call printit_taskX(0," ")
          endif
          first_call = .false.
       endif

       call printit_taskX(0," ")
#ifdef PROGRAM_remap
    endif ! config%use_temp_restart
#endif /* defined(PROGRAM_remap) */

  end subroutine init_eos

    !> SUBROUTINE eos
    !>
    !> Author : L. Huedepohl, A. Marek, originally M. Rampp
    !>
    !>
    !> \verbatim
    !>   SVN - Information  
    !>   $Revision$
    !>   $Date$
    !>   
    !> \endverbatim
    !>
    !>     task         :  provides an interface between PPM and 
    !>                     two EoS: EOS_LS_TB and THJEOS
    !>
    !>     input/output :
    !>
    !> \param rho [in]  The "density" in [g/cc], meaning the baryon number density times atomic unit mass
    !>              
    !>  \param tem [in,out] The temperature in [K]
    !>
    !> \param ede [in,out] internal energy density [erg/cc] (see below for the exact definitions concerning rest-mass)
    !>
    !> \param  p [out] total pressure in [erg/cc]
    !>
    !>  \param gam [out] gam = d(ln p) / d(ln rho) at constant entropy and ye
    !>
    !>  \param s [out] total entropy in [1/by/k_b]
    !>
    !>  \param  xi(*,1:nuc-1) [in,out]
    !>    mass fractions:
    !>    xi(*,1)       = neutron fraction
    !>    xi(*,2)       = proton fraction
    !>    xi(*,3)       = alpha fraction
    !>    xi(*,4:nuc-1) = heavy nuclei fractions
    !>
    !> xi(*,nuc) [in] electron fraction
    !>  \param  xhrep(*) [out] mass fraction of the representative heavy nucleus
    !> za(*,1) [out] Z of representative nucleus
    !> za(*,2) [out] A of representative nucleus
    !>
    !>  \param ccu [out] Neutrino chemical potential without Nucleon rest mass in [MeV]
    !>  \param ccn [out] Neutron chemical potential without Nucleon rest mass in [MeV]
    !>  \param ccp [out] Proton chemical potential without Nucleon rest mass in [MeV]
    !>  \param cce [out] Electron chemical potential with rest mass in [MeV]
    !>
    !>  \param mode [in]
    !>    1  <==>  rho, tem, ye  ==>  ede,  xi, za, p, gam, s
    !>    2  <==>  rho, ede, ye  ==>  tem,          p, gam, s
    !>    3  <==>  rho, ede, ye  ==>  tem,  xi, za, p, gam, s
    !>
    !>  \param nsemode [in]
    !>    0  <==>  the NSE-composition is MAPPED to the 
    !>             ensemble of nuclei being advected
    !>    1  <==>  the NSE-composition is not MAPPED to the 
    !>             ensemble of nuclei being advected
    !>    nsemode only affects the high-density-regime (NSE)
    !>
    !>   Normalization of energy density:
    !>
    !>    *****************************************************************
    !>    *  e_rel = ede + n_b*(m_n*c^2) - n_b*(8.8 MeV) + n_e*(m_e*c^2)  *
    !>    *****************************************************************
    !>
    !>             where: e_rel is the relativistic energy density 
    !>                                            (including rest mass)
    !>                    n_b is the barion number-density
    !>                    m_n is the rest-mass of the neutron
    !>                    m_e is the rest-mass of the electron
    !>
    !>          the value 8.8 MeV (chosen as recommended by Swesty) 
    !>              and used in producing the L&S-tables is approx. equal
    !>              the binding energy (per barion) 
    !>              of Iron, which implies (assuming m_n=m_p) 
    !>              ede ~> 0. for T->0
    !>
    !>   the parameter lsrolo is used to switch (according to value 
    !>         of rho) between L&S-EoS and ThJ-Eos
    !>
    !>----------------------------------------
    subroutine eos(rho, tem, xi, xhrep, za,    ede,      p, gam,  s, ccu, cce, ccn, ccp, selftime, childrentime, &
    !              g/cc   K   1      1   1  erg/cc  erg/cc    1 k_b  MeV  MeV  MeV  MeV         s             s
                  mode, nsemode,  ler, use_ls_eos, renorm)
    !                1        1  bool        bool,    bool


  use precision
  use abort

  use eos_data_type, only : highden_table, lept_radi_table, high_den_eos, &
                            intermediate_eos


  use high_density_eos, only : in_highdensity_eos,                     &
#ifdef DYNAMIC_EOS
                               acquire_sync_locks, release_sync_locks, &
#endif /* DYNAMIC_EOS */
                               highdensity_eos

  use driver_low_density_eos, only : low_density_eos, &
                                      in_low_density_eos
  use eostools, only: my_pack, my_pack_xi
  use phycon
  use nucparam
  use param_rt

#ifdef PROGRAM_remap
  use eos_table_remaper
#endif

  use cputim
  use configure
  use print_stdout_mod
  implicit none

  ! arguments
  real(kind=rk), intent(out)      :: selftime(2), childrentime(2)
  real(kind=rk)                   :: selftime_start(2)
  integer(kind=ik), intent(in)    :: mode,nsemode
  real(kind=rk),    intent(in)    :: rho(:)
  real(kind=rk),    intent(out)   :: p(:), gam(:), s(:), za(:,:), &
                                     xhrep(:), ccu(:), cce(:),    &
                                     ccn(:), ccp(:)

  real(kind=rk),    intent(inout) :: tem(:), xi(:,:), ede(:)
  logical, intent(out)            :: ler

  logical, optional, intent(in)   :: use_ls_eos(:)
  logical, intent(in), optional   :: renorm

  ! variables
  logical :: renorm_flag
  logical, dimension(size(rho))     :: ls_zone,tj_zone,wolff_zone

  integer(kind=ik) :: nuc,i,k,nn
  integer(kind=ik) :: eos_number

  logical :: highdensity_error,tj_error,boundary_error,overlap_error

  logical, save :: first_call=.true.
  logical, save :: first_eos=.true.

  real(kind=rk) :: offs

#ifndef BLOCKWISE_EOS
  integer(kind=ik) :: ls_count, tj_count, wolff_count
  real(kind=rk), dimension(size(rho)) :: rh_wrk,tt_wrk,ed_wrk,xh_wrk,pr_wrk,  &
                                         cu_wrk,ce_wrk,cn_wrk,cp_wrk,gam_wrk,s_wrk
  real(kind=rk), dimension(size(za,dim=1),size(za,dim=2)) :: za_wrk
  real(kind=rk), dimension(size(xi,dim=1),size(xi,dim=2)) :: xi_wrk
#else
  integer(kind=ik) :: kstart,kend
#endif /* BLOCKWISE_EOS */

  selftime     = 0._rk
  childrentime = 0._rk
#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

! the first to do is to check whether we still have a correct
! sum of mass fractions. It turned out that it helps if this
! is checked at every timestep

  nuc = size(xi,dim=2)

  do i=1,size(rho)
     if ( abs(sum(xi(i,1:nuc - 1)) - 1.0_rk) .gt. 1e-6_rk) then
        call printit_taskX(0,"ERROR: Mass fractions do not sum to 1.0 in zone ", i)
        do nn = 1, nuc
           call printit_taskX(0," ",nn, name_xnuc(nn), xi(i, nn))
        enddo
        call printit_taskX(0,"SUM - 1. = ",(sum(xi(i,1:nuc-1)) - 1.0_rk))

        raise_abort("eos.f90(): sum error in massfractions")
     endif
  end do

! determine eos_number stop if more than one EoS is loaded
! (Wofflow_eos does not count here)
  eos_number=0

  do i=1, size(highden_table(:))
     if ((i .ne. 4) .and. (i .eq. highden_table(i)%eos_number)) then
        if (eos_number .ne. 0) then
           raise_abort("eos.f90(): more than one eos-table loaded")
        endif
        eos_number=i
     endif
  enddo

#ifdef CFC_TRANSPORT 
  ! Quick fix: the following procedure might have unexpected consequences
  !            if one forgets about it and it should be removed once 
  !            the root cause for unreasonable temperatures is 
  !            understood and removed       
  tem=max(tem,lept_radi_table(1)%ttmin*1.05_rk)
#endif

  highdensity_error=.false.
  tj_error=.false.
  boundary_error=.false.
  overlap_error=.false.

  ls_zone    = in_highdensity_eos(rho, tem, xi(:,nuc),eos_number)

  if (intermediate_eos .eq. 4) then
     wolff_zone = in_highdensity_eos(rho, tem, xi(:,nuc),4)
  endif

  tj_zone    = in_low_density_eos(rho, tem, xi(:,nuc))

#if defined(CFC_TRANSPORT)
  if (intermediate_eos .ne. 4) then
     if (present(use_ls_eos)) then
        ls_zone=use_ls_eos
        tj_zone=.not. use_ls_eos
     endif
  endif
#endif

  if (ANY(ls_zone .and. tj_zone)) overlap_error=.TRUE.

  if (intermediate_eos .eq. 4) then
     if (ANY(ls_zone .and. wolff_zone)) overlap_error=.TRUE.
     if (ANY(wolff_zone .and. tj_zone)) overlap_error=.TRUE.
  endif


  if (intermediate_eos .eq. 4) then
     if (.not. ALL(ls_zone .or. wolff_zone .or. tj_zone)) boundary_error = .true.
  else
     if (.not. ALL(ls_zone .or. tj_zone)) boundary_error = .true.
  endif

  if (overlap_error) then
    write(*,*) "tj_ls_rho = ", config%tj_ls_rho
    write(*,*) "lsrolo = ", lsrolo
    write(*,*) "highden_table(eos_number)%romin = ", highden_table(eos_number)%romin
    write(*,*) "Overlap error: rho, tem, ye, in LS, in TJ"
    do i = 1, size(ls_zone)
      write(*,*) rho(i), tem(i), xi(i,nuc), ls_zone(i), tj_zone(i)
    end do
    raise_abort('EoS overlap error')
  endif
  if (boundary_error) goto 1000


  if (.not.(present(renorm))) then
    renorm_flag = .true.
  else
    renorm_flag = renorm
  endif

  if(config%restmass_version .gt. 0) then

     if (renorm_flag .eqv. .true.) then

  ! Renormalize
        do k=1,size(ede)

           if (config%restmass_version .eq. 1) then
              offs = moffs + SUM(xi(k,1:n_he4)*mbar(1:n_he4)) + (1._rk-SUM(xi(k,1:n_he4)))*mbar(n_rep)
           endif

           if (config%restmass_version .eq. 2) then

              offs = moffs + SUM(xi(k,1:config%qn-1)*mbar(1:config%qn-1))
           endif


           if (config%restmass_version .eq. 3) then
              offs = moffs + SUM(xi(k,1:config%qn)*mbar(1:config%qn))
           endif

           ede(k) = ede(k) + offs*rho(k)
        enddo

     endif

  endif ! restmass_version > 0


#ifdef DYNAMIC_EOS
  call acquire_sync_locks
#endif


#ifdef BLOCKWISE_EOS
      ! ! ! ! ! ! ! ! ! ! 
      ! High density EoS!
      ! ! ! ! ! ! ! ! ! !

  kstart = 1
  kend = 1


  do while (kend <= size(rho))
     ! Find start of block if possible
     do while (kstart <= size(rho))
        if (.not. ls_zone(kstart)) then
           kstart = kstart + 1
        else
           exit
        endif
     enddo

     if(kstart <= size(rho)) then
        ! Find end of block, this is always possible
        kend = kstart + 1
        do while (kend <= size(rho))
           if (ls_zone(kend)) then
              kend = kend + 1
           else
              exit
           endif
        enddo

          ! Calculate this block.
        call highdensity_eos(eos_number,rho(kstart:kend - 1),  &
                     tem(kstart:kend - 1),  &
                     xi(kstart:kend - 1,:), &
                     xhrep(kstart:kend - 1),&
                     za(kstart:kend - 1,:), &
                     ede(kstart:kend - 1),  &
                     p(kstart:kend - 1),    &
                     gam(kstart:kend - 1),  &
                     s(kstart:kend - 1),    &
                     ccu(kstart:kend - 1),  &
                     cce(kstart:kend - 1),  &
                     ccn(kstart:kend - 1),  &
                     ccp(kstart:kend - 1),  &
                     mode,              &
                     nsemode,           &
                     highdensity_error)
        if (highdensity_error) goto 1000
        kstart = kend
     else
        kend = size(rho) + 1
     endif
  enddo

  if (intermediate_eos .eq. 4) then
      ! ! ! ! ! ! ! ! ! ! ! ! ! !
      ! Intermediate density EoS!
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

     eos_number = intermediate_eos

     kstart = 1
     kend = 1
     do while (kend <= size(rho))
        ! Find start of block if possible
        do while (kstart <= size(rho))
           if (.not. wolff_zone(kstart)) then
              kstart = kstart + 1
           else
              exit
           endif
        enddo

        if(kstart <= size(rho)) then
           ! Find end of block, this is always possible
           kend = kstart + 1
           do while (kend <= size(rho))
              if (wolff_zone(kend)) then
                 kend = kend + 1
              else
                 exit
              endif
           enddo

           ! Calculate this block.
           call highdensity_eos(eos_number,rho(kstart:kend - 1),       &
                                tem(kstart:kend - 1),       &
                                xi(kstart:kend - 1,:),      &
                                xhrep(kstart:kend - 1),     &
                                za(kstart:kend - 1,:),      &
                                ede(kstart:kend - 1),       &
                                p(kstart:kend - 1),         &
                                gam(kstart:kend - 1),       &
                                s(kstart:kend - 1),         &
                                ccu(kstart:kend - 1),       &
                                cce(kstart:kend - 1),       &
                                ccn(kstart:kend - 1),       &
                                ccp(kstart:kend - 1),       &
                                mode,                   &
                                nsemode,                &
                                highdensity_error)
           if (highdensity_error) goto 1000
           kstart = kend
        else
           kend = size(rho) + 1
        endif
     enddo
  endif ! intermediate_eos .eq. 4

      ! ! ! ! ! ! ! ! ! ! 
      ! Low density EoS !
      ! ! ! ! ! ! ! ! ! !
  kstart = 1
  kend = 1
  do while (kend <= size(rho))
     ! Find start of block if possible
     do while (kstart <= size(rho))
        if (.not. tj_zone(kstart)) then
           kstart = kstart + 1
        else
           exit
        endif
     enddo

     if(kstart <= size(rho)) then
        ! Find end of block, this is always possible
        kend = kstart + 1
        do while (kend <= size(rho))
           if (tj_zone(kend)) then
              kend = kend + 1
           else
              exit
           endif
        enddo
          ! Calculate this block.
        call low_density_eos(rho(kstart:kend - 1),       &
                     tem(kstart:kend - 1),       &
                     xi(kstart:kend - 1,:),      &
                     xhrep(kstart:kend - 1),     &
                     za(kstart:kend - 1,:),      &
                     ede(kstart:kend - 1),       &
                     p(kstart:kend - 1),         &
                     gam(kstart:kend - 1),       &
                     s(kstart:kend - 1),         &
                     ccu(kstart:kend - 1),       &
                     cce(kstart:kend - 1),       &
                     ccn(kstart:kend - 1),       &
                     ccp(kstart:kend - 1),       &
                     mode,                   &
                     nsemode,                &
                     tj_error)
        if (tj_error) goto 1000
        kstart = kend
     else
        kend = size(rho) + 1
     endif
  enddo
#else /* BLOCKWISE_EOS */
      ! ! ! ! ! ! ! ! ! ! 
      ! High density EoS!
      ! ! ! ! ! ! ! ! ! !




  ls_count = COUNT(ls_zone)

  if (ls_count .gt. 0) then
     ! gather
     rh_wrk = my_pack(rho, ls_zone)
     tt_wrk = my_pack(tem, ls_zone)
     xi_wrk(:,:) = my_pack_xi(xi(:,:), ls_zone)
     if (mode.eq.2 .or. mode.eq.3) ed_wrk = my_pack(ede, ls_zone)
     if (mode.eq.4) s_wrk = my_pack(s, ls_zone)
     call highdensity_eos(eos_number,rh_wrk(1:ls_count),  &
                   tt_wrk(1:ls_count),  &
                   xi_wrk(1:ls_count,:),&
                   xh_wrk(1:ls_count),  &
                   za_wrk(1:ls_count,:),&
                   ed_wrk(1:ls_count),  &
                   pr_wrk(1:ls_count),  &
                   gam_wrk(1:ls_count), &
                   s_wrk(1:ls_count),   &
                   cu_wrk(1:ls_count),  &
                   ce_wrk(1:ls_count),  &
                   cn_wrk(1:ls_count),  &
                   cp_wrk(1:ls_count),  &
                   mode,                &
                   nsemode,             &
                   highdensity_error )
     if (highdensity_error) goto 1000


     ! scatter
     if (mode.eq.1.or. mode.eq.4) ede = unpack(ed_wrk, ls_zone, ede)
     if (mode.eq.2 .or. mode.eq.3.or. mode.eq.4) tem = unpack(tt_wrk, ls_zone, tem)
     p   = unpack(pr_wrk,  ls_zone, p)
     gam = unpack(gam_wrk, ls_zone, gam)
     s   = unpack(s_wrk,   ls_zone, s)
     ccu = unpack(cu_wrk,  ls_zone, ccu)
     cce = unpack(ce_wrk,  ls_zone, cce)
     ccn = unpack(cn_wrk,  ls_zone, ccn)
     ccp = unpack(cp_wrk,  ls_zone, ccp)

     do nn = 1, nuc - 1 
        xi(:,nn) = unpack(xi_wrk(:,nn), ls_zone, xi(:,nn))
     enddo

     if (nsemode .eq. 1) then
        xhrep   = unpack(xh_wrk,      ls_zone, xhrep)
        za(:,1) = unpack(za_wrk(:,1), ls_zone, za(:,1))
        za(:,2) = unpack(za_wrk(:,2), ls_zone, za(:,2))
     endif
  endif



  if (intermediate_eos .eq. 4) then
     ! ! ! ! ! ! ! ! ! ! ! ! ! !
     ! Intermediate density EoS!
     ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

     eos_number = intermediate_eos

     wolff_count = COUNT(wolff_zone)


     if (wolff_count .gt. 0) then
        ! gather
        rh_wrk = my_pack(rho,  wolff_zone)
        tt_wrk = my_pack(tem, wolff_zone)
        xi_wrk(:,:) = my_pack_xi(xi(:,:),     wolff_zone)
        if (mode.eq.2 .or. mode.eq.3) ed_wrk = my_pack(ede, wolff_zone)
        if (mode.eq.4) s_wrk = my_pack(s, wolff_zone)
        call highdensity_eos(eos_number,rh_wrk(1:wolff_count), &
                             tt_wrk(1:wolff_count), &
                             xi_wrk(1:wolff_count,:),&
                             xh_wrk(1:wolff_count), &
                             za_wrk(1:wolff_count,:),&
                             ed_wrk(1:wolff_count), &
                             pr_wrk(1:wolff_count), &
                             gam_wrk(1:wolff_count),&
                             s_wrk(1:wolff_count),  &
                             cu_wrk(1:wolff_count), &
                             ce_wrk(1:wolff_count), &
                             cn_wrk(1:wolff_count), &
                             cp_wrk(1:wolff_count), & 
                             mode,                  &
                             nsemode,               &
                             highdensity_error)
        if (highdensity_error) goto 1000
        ! scatter
        if (mode .eq. 1 .or. mode .eq. 4) ede = unpack(ed_wrk, wolff_zone, ede)

        if (mode .eq. 2 .or. mode .eq. 3 .or. mode .eq. 4) &
             tem = unpack(tt_wrk, wolff_zone, tem)

        p   = unpack(pr_wrk,  wolff_zone, p)
        gam = unpack(gam_wrk, wolff_zone, gam)
        s   = unpack(s_wrk,   wolff_zone, s)
        ccu = unpack(cu_wrk,  wolff_zone, ccu)
        cce = unpack(ce_wrk,  wolff_zone, cce)
        ccn = unpack(cn_wrk,  wolff_zone, ccn)
        ccp = unpack(cp_wrk,  wolff_zone, ccp)

        do nn = 1, nuc - 1
           xi(:,nn) = unpack(xi_wrk(1:wolff_count,nn), wolff_zone, xi(:,nn))
        enddo

        if (nsemode .eq. 1) then
           xhrep   = unpack(xh_wrk,      wolff_zone, xhrep)
           za(:,1) = unpack(za_wrk(:,1), wolff_zone, za(:,1))
           za(:,2) = unpack(za_wrk(:,2), wolff_zone, za(:,2))
        endif
     endif

  endif ! intermediate_eos = 4



!----NON-NSE EoS

      ! ! ! ! ! ! ! ! ! ! 
      ! Low density EoS !
      ! ! ! ! ! ! ! ! ! !
  tj_count = COUNT(tj_zone)



  if (tj_count .gt. 0) then

     ! gather
     rh_wrk = my_pack(rho, tj_zone)
     tt_wrk = my_pack(tem, tj_zone)
     xi_wrk = my_pack_xi(xi, tj_zone)

     if (mode .eq. 2 .or. mode .eq. 3) ed_wrk = my_pack(ede, tj_zone)
     if        (mode .eq. 4)            s_wrk = my_pack(s,   tj_zone)

     call low_density_eos(rh_wrk(1:tj_count),  &
                   tt_wrk(1:tj_count),  &
                   xi_wrk(1:tj_count,1:nuc),&
                   xh_wrk(1:tj_count),  &
                   za_wrk(1:tj_count,1:2),&
                   ed_wrk(1:tj_count),  &
                   pr_wrk(1:tj_count),  &
                   gam_wrk(1:tj_count), &
                   s_wrk(1:tj_count),   &
                   cu_wrk(1:tj_count),  &
                   ce_wrk(1:tj_count),  &
                   cn_wrk(1:tj_count),  &
                   cp_wrk(1:tj_count),  & 
                   mode,                &
                   nsemode,             &
                   tj_error)

     if (tj_error) goto 1000

        ! scatter
     if (mode .eq. 1 .or. mode .eq. 4) ede = unpack(ed_wrk, tj_zone, ede)

     if (mode.eq.2 .or. mode.eq.3.or. mode.eq.4) &
          tem = unpack(tt_wrk, tj_zone, tem)

     p   = unpack(pr_wrk,  tj_zone, p)
     gam = unpack(gam_wrk, tj_zone, gam)
     s   = unpack(s_wrk,   tj_zone, s)
     ccu = unpack(cu_wrk,  tj_zone, ccu)
     cce = unpack(ce_wrk,  tj_zone, cce)
     ccn = unpack(cn_wrk,  tj_zone, ccn)
     ccp = unpack(cp_wrk,  tj_zone, ccp)

     do nn = 1, nuc - 1
        xi(:,nn) = unpack(xi_wrk(1:tj_count,nn), tj_zone, xi(:,nn))
     enddo
     if (nsemode .eq. 1) then 
        where(tj_zone)
           xhrep(:)=0._rk
           za(:,1)=0._rk
           za(:,2)=0._rk
        endwhere


!           xhrep   = unpack(xh_wrk,      tj_zone, xhrep)
!           za(:,1) = unpack(za_wrk(:,1), tj_zone, za(:,1))
!           za(:,2) = unpack(za_wrk(:,2), tj_zone, za(:,2))
     endif
  endif


#endif /* BLOCKWISE_EOS */

#ifdef DYNAMIC_EOS
  call release_sync_locks
#endif

! correct neutrino chemical potential  
  ccu (:) = cce (:) + ccp (:) -ccn (:) - wc_mq
  
  if (config%restmass_version .gt. 0) then
     if (renorm_flag .eqv. .true.) then

      ! Denormalize
        do k=1,size(ede)

           if (config%restmass_version .eq. 1) then
              offs = moffs + SUM(xi(k,1:n_he4)*mbar(1:n_he4)) + (1._rk-SUM(xi(k,1:n_he4)))*mbar(n_rep)
           endif

           if (config%restmass_version .eq. 2) then
              offs = moffs + SUM(xi(k,1:config%qn-1)*mbar(1:config%qn-1))
           endif

           if (config%restmass_version .eq. 3) then
              offs = moffs + SUM(xi(k,1:config%qn)*mbar(1:config%qn))
           endif

           ede(k) = ede(k) - offs*rho(k)
           if (ede(k) < 0.) then
              write (*,*) 'Negative energy density at index ', k
              raise_abort('Negative energy density')
           endif
        enddo
     endif
  endif ! restmass_version > 0


  ler=.false.
#ifndef DEBUG_TIMINGS
 call second_v(selftime)
 selftime = selftime - selftime_start
#endif
  return

!-- error handling

1000 continue

  if (boundary_error) then
     write(*,*)
     write(*,'(''eos> EoS boundaries exceeded: '', ''[k/rho/t/ede/ye]'')')
     write(*,*) 'mode=',mode
     write(*,*)
     do k=1,size(rho)
#ifdef WOLFFEOS_LOW
        if (.not.(ls_zone(k) .or. wolff_zone(k) .or. tj_zone(k))) then
#else
           if (.not.(ls_zone(k) .or. tj_zone(k))) then
#endif
              write(*,'(I4,4(1x,1pe12.5))')  k, rho(k), tem(k), ede(k), xi(k,nuc)
          endif
       enddo
    elseif (highdensity_error) then
       write(*,*) 'eos> error in ',highden_table(eos_number)%name,"  mode = ",mode
    elseif (tj_error) then
       write(*,*) 'eos> error in ThJ EoS; mode = ',mode
    else
       write(*,*) 'eos> PANIC; mode = ',mode
       raise_abort("eos> PANIC")
    endif
    ler=.true.
  end subroutine eos

!> this function is a wrapper that is called in the c-network code.
!> for a given density, temperature, and composition it computes
!> the energy density and dedtem
!>
!> But be carefully : The network code does NOT use the dummy nuclei (i.e.
!> only xnuc(1:config%qn_network) instead of xnuc(1:config%qn-1) AND this subset of nuclei
!> was renomalized to unity. This has to be redone here before calling the
!> eos and afterwards
!>
!> Author A. Marek, Aug. 2009
!>
!>  \param mode
!>  \param rho_in  density
!>  \param tem_im  temperature
!>  \param xi_in   composition
!>  \param ede_out energy density
!>  \param dedtem_out dEnergy_density/dtemperature
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
!> available in c as:
!> extern void eos_network_wrapper(int *mode, double *rho, double *tem,double *xi, double *e, double *dedt);
!>
  subroutine eos_network_wrapper(mode,rho_in,tem_in,xi_in,ede_out,dedtem_out) bind(C)
    use abort
    use precision
    use iso_c_binding

    use netw_thiel , only : ye_netw, xi_dummy_netw, l_renorm_xnuc
!    use nucparam, only : netw_species

    use configure
    implicit none

! Attention: The network-code uses confiq%qn_network species, the hydro / transport
!            code uses config%qn species! config%qn = config%qn_network + 2 "dummy nuclei" + Ye

    real(kind=rk)                           :: eos_self(2), eos_children(2)  
    real(kind=C_DOUBLE),  intent(in)        :: rho_in, tem_in, &
                                               xi_in(1:config%qn_network)
    integer(kind=C_INT),  intent(in)        :: mode
    real(kind=C_DOUBLE), intent(out)        :: ede_out,dedtem_out

    real(kind=rk), dimension(1)             :: rho,tem,ede,dummy,ede_up,tem_up
    real(kind=rk), dimension(1,1:config%qn) :: xi
    real(kind=rk), dimension(1,2)           :: za
    real(kind=rk)                           ::h, summe

    real(kind=rk), parameter :: eps = 1.e-7_rk

    integer(kind=ik), parameter :: nsemode=1
    integer(kind=ik) :: i

    logical :: error

    rho(1)=real(rho_in, kind=rk)
    tem(1)=real(tem_in, kind=rk)
! copy composition parsed from the network to the array needed in EoS
    xi(1,1:config%qn_network)=real(xi_in(1:config%qn_network), kind=rk)

! renormalize the mass fractions

    if (l_renorm_xnuc .eq. 1) then
        do i=1,config%qn_network
           xi(1,i)=xi(1,i)*(sum(xi_in(1:config%qn_network)) - sum(xi_dummy_netw(1:2)))
        enddo
    endif

! copy composition of the dummy nuclei
    do i=1,2
       xi(1,config%qn_network+i)= xi_dummy_netw(i)
    enddo

! copy Ye to EoS array
    xi(1,config%qn) = ye_netw


    dummy(1)=0._rk
    za(1,:) =0._rk

    if (abs(sum(xi(1,1:config%qn -1)) -1._rk) .gt. 1.e-11) then
       summe =sum(xi(1,1:config%qn -1))
       do i=1,config%qn-1
          xi(1,i)=xi(1,i)/summe
       enddo
    endif

! compute energy density
    call eos(rho,tem,xi,dummy,za,ede,dummy,dummy,dummy,dummy,dummy, &
                        dummy,dummy,eos_self, eos_children,int(mode, kind=ik),nsemode,error)

    if (error) raise_abort("eos_network_wrapper(): error 1")


! network code wants erg/g
    ede(1)= ede(1)/rho(1)

    ede_out = real(ede(1), kind=C_DOUBLE)


! compute ede + epsilon

    tem_up=tem    + eps*tem
    h     =tem_up(1) - tem(1)

! compute energy density
    call eos(rho,tem_up,xi,dummy,za,ede_up,dummy,dummy,dummy,dummy,dummy, &
                        dummy,dummy,eos_self, eos_children,int(mode, kind=ik),nsemode,error)

    if (error) raise_abort("eos_network_wrapper(): error 2")

! network code wants erg/g
    ede_up(1) = ede_up(1) / rho(1)
! compute dedtem

    dedtem_out = real((ede_up(1) - ede(1))/h, kind=C_DOUBLE)

  end subroutine eos_network_wrapper

!-----------------------------------------------------------------------
end module eos_sn2


