!>
!> \verbatim
!> this module provides all the high-density EoS tables
!>
!>  Author: A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module eos_data_type
  use precision
  use iso_c_binding

  implicit none
#if !(defined(PROGRAM_remap))
#else /* PROGRAM_remap */
  public
#endif

  save

#ifndef PROGRAM_remap
  integer(kind=ik), parameter :: high_den_eos     = HIGH_DEN_EOS_FLAG
  integer(kind=ik), parameter :: intermediate_eos = INTERMEDIATE_EOS_FLAG
#else /* PROGRAM_remap */
  integer(kind=ik)            :: high_den_eos     = 1  != eosfile
  integer(kind=ik)            :: intermediate_eos = 3  != lowden +3 
                                                       ! careful lwoden = 1 or 0 
                                                       ! but intermediate eos 
                                                       ! should be 4 or 0
#endif /* PROGRAM_remap */

#ifdef DYNAMIC_EOS
  integer(kind=ik) :: thread_counter = 0, reload_counter = 0
#endif


#if defined(USE_MEMORY_MAPPING) || defined(DYNAMIC_EOS)

   integer(kind=c_int), bind(c ,name='hd_eos_timestamp') :: timestamp_c_f
   integer(kind=c_int), bind(c ,name='hd_eos_nrho') :: nro_c_f
   integer(kind=c_int), bind(c, name='hd_eos_ntt')  :: ntt_c_f
   integer(kind=c_int), bind(c, name='hd_eos_nye')  :: nye_c_f

   ! the "global" values of the whole table
   real(kind=c_float), bind(c ,name='den_min') :: romin_c_f
   real(kind=c_float), bind(c ,name='den_max') :: romax_c_f
   real(kind=c_float), bind(c ,name='tem_min') :: ttmin_c_f
   real(kind=c_float), bind(c ,name='tem_max') :: ttmax_c_f
   real(kind=c_float), bind(c ,name='ye_min')  :: yemin_c_f
   real(kind=c_float), bind(c ,name='ye_max')  :: yemax_c_f

#ifdef DYNAMIC_EOS
   ! the "local" values of the actually loaded part of the table

   real(kind=c_float), bind(c ,name='den_min_lc') :: romin_lc_c_f
   real(kind=c_float), bind(c ,name='den_max_lc') :: romax_lc_c_f
   real(kind=c_float), bind(c ,name='tem_min_lc') :: ttmin_lc_c_f
   real(kind=c_float), bind(c ,name='tem_max_lc') :: ttmax_lc_c_f
   real(kind=c_float), bind(c ,name='ye_min_lc')  :: yemin_lc_c_f
   real(kind=c_float), bind(c ,name='ye_max_lc')  :: yemax_lc_c_f
#endif /* DYNAMIC_EOS */
#endif /*  USE_MEMORY_MAPPING */

  type eos_tables

!> this data type is used for the High Density EoS and
!> for the Low Density NSE -EoS. Variables not needed in either
!> of both will be nullified at the corresponding loading routine


! variables that are used in all tables
   character(len=30) :: name
   character(120)    :: tablefile

   integer(kind=ik)  :: nro,ntt,nye
   real(kind=rk)     :: romax, romin, ttmax, ttmin, yemax, yemin

#ifdef DYNAMIC_EOS
   ! local values
   integer(kind=ik)  :: nro_lc,ntt_lc,nye_lc
   real(kind=rk)     :: romax_lc, romin_lc, ttmax_lc, ttmin_lc, &
                        yemax_lc, yemin_lc
   integer(kind=ik)  :: nro_max_lc, nro_min_lc, ntt_max_lc, ntt_min_lc, &
                        nye_max_lc, nye_min_lc
#endif /* DYNAMIC_EOS */


   real(kind=rk)     :: royemax, royemin
#ifdef ONEMG_EOS
   real(kind=rk)     :: tem_cut
#endif

   logical           :: eos_energy_offs ! will be set if the
                                        ! engery offset in the EoS
                                        ! is NOT 8.8 (i.e. enullm)
   real(kind=rk)     :: enullm
  
! variables only used in the lepton-radiation table
! lro_ld is rho*ye here
   real(kind=rk), dimension(:,:), pointer     :: let, le, lp, ls, gam
   real(kind=rk)                              :: dlttab, dldtab
   real(kind=rk), dimension(:,:,:,:), pointer :: be, bp

! variables only used in Low Density NSE EoS
   integer(kind=ik) :: nspec
   real(kind=rk)    :: aa,zz
   integer(kind=ik) :: n_hv

   real(kind=rk), dimension(:), pointer :: lro_ld, ltt_ld, ye_ld
   real(kind=rk), dimension(:,:,:), pointer :: eltab,hetab,xhtab
   real(kind=rk), dimension(:,:,:,:), pointer :: xitab
   real(kind=rk) :: eltab_min

#ifdef ONEMG_EOS
   real(kind=rk) :: romax_low_ld, romax_high_ld
#endif

! variables only used in High Density EoS
   integer(kind=ik)  :: eos_number
! the following numbers will be used (this vallue is set by HIGH_DEN_EOS_FLAG
! and /or INTERMEDIATE_EOS_FLAG in make.inc.logic)
! 1 = LS-EoS       = HIGH_DEN_EOS_FLAG
! 2 = WOLFF-EoS    = HIGH_DEN_EOS_FLAG
! 3 = SHEN-EOS     = HIGH_DEN_EOS_FLAG
! 4 = WOLFFEOS-LOW = INTERMEDIATE_EOS_FLAG

   integer(kind=ik) :: timestamp

#ifdef ONEMG_EOS
   real(kind=rk) :: romin_low_hd, romin_high_hd
#endif

   real(kind=rk_eos), dimension(:),     pointer :: lro_hd, ltt_hd, ye_hd
#ifdef DYNAMIC_EOS
   real(kind=rk_eos), dimension(:),     pointer :: lro_hd_lc, ltt_hd_lc, ye_hd_lc
#endif /* DYNAMIC_EOS */
   real(kind=rk_eos), dimension(:,:,:), pointer :: lpr, led, cei, cni, cpi, &
                                                   cui, gai, sti, xxp, xxn, &
                                                   xxa, xxh, xhz, xha

#ifdef NSMIX
   real(kind=rk_eos), dimension(:,:,:), pointer :: stn
#endif /* NSMIX */

  end type eos_tables

!#ifdef DYNAMIC_EOS
!  type pointers
!   real(kind=rk_eos), dimension(:),     pointer :: lro_hd, ltt_hd, ye_hd
!   real(kind=rk_eos), dimension(:,:,:), pointer :: lpr, led, cei, cni, cpi, &
!                                                   cui, gai, sti, xxp, xxn, &
!                                                   xxa, xxh, xhz, xha
!  end type pointers
!
!  type(pointers) :: dynamic_eos(10)
!#endif /* DYNAMIC_EOS */

#if defined(USE_MEMORY_MAPPING) || defined(DYNAMIC_EOS)
  type(c_ptr), BIND(C,name='hd_eos_lro_ptr') :: lro_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_ltt_ptr') :: ltt_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_ye_ptr')  :: ye_hd_c_f
#ifdef DYNAMIC_EOS
  type(c_ptr), BIND(C,name='hd_eos_lro_lc_ptr') :: lro_hd_lc_c_f
  type(c_ptr), BIND(C,name='hd_eos_ltt_lc_ptr') :: ltt_hd_lc_c_f
  type(c_ptr), BIND(C,name='hd_eos_ye_lc_ptr')  :: ye_hd_lc_c_f
#endif
  type(c_ptr), BIND(C,name='hd_eos_lpr_ptr') :: lpr_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_led_ptr') :: led_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_cei_ptr') :: cei_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_cni_ptr') :: cni_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_cpi_ptr') :: cpi_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_cui_ptr') :: cui_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_gai_ptr') :: gai_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_sti_ptr') :: sti_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_xxn_ptr') :: xxn_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_xxp_ptr') :: xxp_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_xxa_ptr') :: xxa_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_xxh_ptr') :: xxh_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_xha_ptr') :: xha_hd_c_f
  type(c_ptr), BIND(C,name='hd_eos_xhz_ptr') :: xhz_hd_c_f
#endif /* USE_MEMORY_MAPPING */

  type(eos_tables) :: highden_table(10)
  type(eos_tables) :: lowden_table(3)
  type(eos_tables) :: lept_radi_table(1)

contains 
!>
!> \verbatim
!> This functions nullifies all pointers. It should be called BEFORE
!> the pointers are assigned to targets
!>
!>  Author: A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>   
  subroutine nullify_pointer

  use precision
  use abort

    implicit none
      
    integer(kind=ik) :: i
 
    do i=1,10

       nullify(highden_table(i)%let,highden_table(i)%le,      &
               highden_table(i)%lp,highden_table(i)%ls,       &
               highden_table(i)%gam,highden_table(i)%be,      &
               highden_table(i)%bp,highden_table(i)%lro_ld,   &
               highden_table(i)%ltt_ld,highden_table(i)%ye_ld,&
               highden_table(i)%eltab,highden_table(i)%hetab, &
               highden_table(i)%xhtab,highden_table(i)%xitab)

       nullify(highden_table(i)%lro_hd,highden_table(i)%ltt_hd, &
               highden_table(i)%ye_hd,highden_table(i)%lpr,     &
               highden_table(i)%led,highden_table(i)%cei,       &
               highden_table(i)%cni,highden_table(i)%cpi,       &
               highden_table(i)%cui,highden_table(i)%gai,       &
               highden_table(i)%sti,highden_table(i)%xxp,       &
               highden_table(i)%xxn,highden_table(i)%xxa,       &
               highden_table(i)%xxh,highden_table(i)%xhz,       &
               highden_table(i)%xha)

#ifdef DYNAMIC_EOS
       nullify(highden_table(i)%lro_hd_lc,highden_table(i)%ltt_hd_lc, &
               highden_table(i)%ye_hd_lc)

#endif /* DYNAMIC_EOS */

!#ifdef DYNAMIC_EOS
!       nullify(dynamic_eos(i)%lro_hd, dynamic_eos(i)%ltt_hd,   &
!               dynamic_eos(i)%ye_hd, dynamic_eos(i)%lpr,       &
!               dynamic_eos(i)%led, dynamic_eos(i)%cei,         &
!               dynamic_eos(i)%cni, dynamic_eos(i)%cpi,         &
!               dynamic_eos(i)%cui, dynamic_eos(i)%gai,         &
!               dynamic_eos(i)%sti, dynamic_eos(i)%xxp,         &
!               dynamic_eos(i)%xxn, dynamic_eos(i)%xxa,         &
!               dynamic_eos(i)%xxh, dynamic_eos(i)%xhz,         &
!               dynamic_eos(i)%xha)
!
!#endif /* DYNAMIC_EOS */

    enddo

    do i=1,3
     
       nullify(lowden_table(i)%let,lowden_table(i)%le,      &
               lowden_table(i)%lp,lowden_table(i)%ls,       &
               lowden_table(i)%gam,lowden_table(i)%be,      &
               lowden_table(i)%bp,lowden_table(i)%lro_ld,   &
               lowden_table(i)%ltt_ld,lowden_table(i)%ye_ld,&
               lowden_table(i)%eltab,lowden_table(i)%hetab, &
               lowden_table(i)%xhtab,lowden_table(i)%xitab)

       nullify(lowden_table(i)%lro_hd,lowden_table(i)%ltt_hd, &
               lowden_table(i)%ye_hd,lowden_table(i)%lpr,     &
               lowden_table(i)%led,lowden_table(i)%cei,       &
               lowden_table(i)%cni,lowden_table(i)%cpi,       &
               lowden_table(i)%cui,lowden_table(i)%gai,       &
               lowden_table(i)%sti,lowden_table(i)%xxp,       &
               lowden_table(i)%xxn,lowden_table(i)%xxa,       &
               lowden_table(i)%xxh,lowden_table(i)%xhz,       &
               lowden_table(i)%xha)

#ifdef DYNAMIC_EOS
       nullify(lowden_table(i)%lro_hd_lc,lowden_table(i)%ltt_hd_lc, &
               lowden_table(i)%ye_hd_lc)

#endif /* DYNAMIC_EOS */

    enddo


    do i=1,1
    
       nullify(lept_radi_table(i)%let,lept_radi_table(i)%le,      &
               lept_radi_table(i)%lp,lept_radi_table(i)%ls,       &
               lept_radi_table(i)%gam,lept_radi_table(i)%be,      &
               lept_radi_table(i)%bp,lept_radi_table(i)%lro_ld,   &
               lept_radi_table(i)%ltt_ld,lept_radi_table(i)%ye_ld,&
               lept_radi_table(i)%eltab,lept_radi_table(i)%hetab, &
               lept_radi_table(i)%xhtab,lept_radi_table(i)%xitab)

       nullify(lept_radi_table(i)%lro_hd,lept_radi_table(i)%ltt_hd,  &
               lept_radi_table(i)%ye_hd,lept_radi_table(i)%lpr,      &
               lept_radi_table(i)%led,lept_radi_table(i)%cei,        &
               lept_radi_table(i)%cni,lept_radi_table(i)%cpi,        &
               lept_radi_table(i)%cui,lept_radi_table(i)%gai,        &
               lept_radi_table(i)%sti,lept_radi_table(i)%xxp,        &
               lept_radi_table(i)%xxn,lept_radi_table(i)%xxa,        &
               lept_radi_table(i)%xxh,lept_radi_table(i)%xhz,        &
               lept_radi_table(i)%xha)
#ifdef DYNAMIC_EOS
       nullify(lept_radi_table(i)%lro_hd_lc,lept_radi_table(i)%ltt_hd_lc,  &
               lept_radi_table(i)%ye_hd_lc)
#endif /* DYNAMIC_EOS */
    enddo

  end subroutine  nullify_pointer

!>
!> \verbatim
!> This functions "zeros" all variables. 
!>
!>  Author: A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>  
  subroutine zero_data_variables

    use precision

    implicit none

! this subroutine sets all variables that are defined in the eos types to zero

    integer(kind=ik) :: i
 
    do i=1,10

       ! global values
       highden_table(i)%name        =" "
       highden_table(i)%nro         =0
       highden_table(i)%ntt         =0
       highden_table(i)%nye         =0
       highden_table(i)%romax       =0._rk
       highden_table(i)%romin       =0._rk
       highden_table(i)%ttmax       =0._rk
       highden_table(i)%ttmin       =0._rk
       highden_table(i)%yemax       =0._rk
       highden_table(i)%yemin       =0._rk

#ifdef DYNAMIC_EOS
       ! local values
       highden_table(i)%nro_lc      =0
       highden_table(i)%ntt_lc      =0
       highden_table(i)%nye_lc      =0
       highden_table(i)%romax_lc    =0._rk
       highden_table(i)%romin_lc    =0._rk
       highden_table(i)%ttmax_lc    =0._rk
       highden_table(i)%ttmin_lc    =0._rk
       highden_table(i)%yemax_lc    =0._rk
       highden_table(i)%yemin_lc    =0._rk
#endif
       highden_table(i)%royemax     =0._rk
       highden_table(i)%royemin     =0._rk
#ifdef ONEMG_EOS
       highden_table(i)%tem_cut     =0._rk
#endif
       highden_table(i)%nspec       =0
       highden_table(i)%aa          =0._rk
       highden_table(i)%zz          =0._rk
       highden_table(i)%n_hv        =0

#ifdef ONEMG_EOS
       highden_table(i)%romax_low_ld =0._rk
       highden_table(i)%romax_high_ld=0._rk
#endif
       highden_table(i)%eos_number   =0

    enddo
    
    do i=1,3

       lowden_table(i)%name        =" "
       ! global values
       lowden_table(i)%nro         =0
       lowden_table(i)%ntt         =0
       lowden_table(i)%nye         =0
       lowden_table(i)%romax       =0._rk
       lowden_table(i)%romin       =0._rk
       lowden_table(i)%ttmax       =0._rk
       lowden_table(i)%ttmin       =0._rk
       lowden_table(i)%yemax       =0._rk
       lowden_table(i)%yemin       =0._rk

#ifdef DYNAMIC_EOS
       ! local values
       lowden_table(i)%nro_lc      =0
       lowden_table(i)%ntt_lc      =0
       lowden_table(i)%nye_lc      =0
       lowden_table(i)%romax_lc    =0._rk
       lowden_table(i)%romin_lc    =0._rk
       lowden_table(i)%ttmax_lc    =0._rk
       lowden_table(i)%ttmin_lc    =0._rk
       lowden_table(i)%yemax_lc    =0._rk
       lowden_table(i)%yemin_lc    =0._rk
#endif
       lowden_table(i)%royemax     =0._rk
       lowden_table(i)%royemin     =0._rk
#ifdef ONEMG_EOS
       lowden_table(i)%tem_cut     =0._rk
#endif
       lowden_table(i)%nspec       =0
       lowden_table(i)%aa          =0._rk
       lowden_table(i)%zz          =0._rk
       lowden_table(i)%n_hv        =0

#ifdef ONEMG_EOS
       lowden_table(i)%romax_low_ld =0._rk
       lowden_table(i)%romax_high_ld=0._rk
#endif
       lowden_table(i)%eos_number   =0
     
    enddo
 
    do i=1,1

       lept_radi_table(i)%name       =" "
       ! global values
       lept_radi_table(i)%nro        =0
       lept_radi_table(i)%ntt        =0
       lept_radi_table(i)%nye        =0
       lept_radi_table(i)%romax      =0._rk
       lept_radi_table(i)%romin      =0._rk
       lept_radi_table(i)%ttmax      =0._rk
       lept_radi_table(i)%ttmin      =0._rk
       lept_radi_table(i)%yemax      =0._rk
       lept_radi_table(i)%yemin      =0._rk
       
#ifdef DYNAMIC_EOS
       ! local values
       lept_radi_table(i)%nro_lc=0
       lept_radi_table(i)%ntt_lc=0
       lept_radi_table(i)%nye_lc=0
       lept_radi_table(i)%romax_lc=0._rk
       lept_radi_table(i)%romin_lc=0._rk
       lept_radi_table(i)%ttmax_lc=0._rk
       lept_radi_table(i)%ttmin_lc=0._rk
       lept_radi_table(i)%yemax_lc=0._rk
       lept_radi_table(i)%yemin_lc=0._rk
#endif
       lept_radi_table(i)%royemax =0._rk
       lept_radi_table(i)%royemin =0._rk
#ifdef ONEMG_EOS
       lept_radi_table(i)%tem_cut =0._rk
#endif
       lept_radi_table(i)%nspec   =0
       lept_radi_table(i)%aa      =0._rk
       lept_radi_table(i)%zz      =0._rk
       lept_radi_table(i)%n_hv     =0

#ifdef ONEMG_EOS
       lept_radi_table(i)%romax_low_ld=0._rk
       lept_radi_table(i)%romax_high_ld=0._rk
#endif
       lept_radi_table(i)%eos_number=0
       
    enddo

  end subroutine zero_data_variables
#if defined(USE_MEMORY_MAPPING) || defined(DYNAMIC_EOS)
!>
!> \verbatim
!> this module provides a mapping between the c pointers
!> and fortran pointers if memorz mapping is used
!>
!>  Author: A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine associate_hd_eos_c_to_fortran(eos_number, nro_start, nro_end, &
                                           ntt_start, ntt_end, nye_start,  &
                                           nye_end)

    use precision
    use abort
    use phycon
    use iso_c_binding

    implicit none 

    integer(kind=ik), intent(in) :: eos_number, nro_start, nro_end, &
                                    ntt_start, ntt_end, nye_start,  &
                                    nye_end

    integer(kind=ik) :: nro, ntt, nye, i, istat
    integer(kind=ik) :: nro_shape, ntt_shape, nye_shape
    
    logical, save :: first_call = .true. 

    ! global values 
    highden_table(eos_number)%timestamp = timestamp_c_f                            
    highden_table(eos_number)%nro = nro_c_f
    highden_table(eos_number)%ntt = ntt_c_f
    highden_table(eos_number)%nye = nye_c_f

    nro = nro_c_f
    ntt = ntt_c_f
    nye = nye_c_f

    if (first_call) then
       highden_table(eos_number)%romin    = 10.0**real(romin_c_f,kind=rk)
#ifdef DYNAMIC_EOS
       highden_table(eos_number)%romin_lc = 10.0**real(romin_lc_c_f,kind=rk)
#endif
       first_call = .false. 
    endif

    highden_table(eos_number)%romax = 10.0**real(romax_c_f,kind=rk)
    highden_table(eos_number)%ttmin = 10.0**(real(ttmin_c_f,kind=rk))
    highden_table(eos_number)%ttmax = 10.0**(real(ttmax_c_f,kind=rk))
    highden_table(eos_number)%yemin = real(yemin_c_f,kind=rk)
    highden_table(eos_number)%yemax = real(yemax_c_f,kind=rk)

    nro_shape = nro
    ntt_shape = ntt
    nye_shape = nye

#ifdef DYNAMIC_EOS
    ! local values                                
    ! nro_lc, ntt_lc, and nye_lc contain the number of the
    ! loaded table points 
    highden_table(eos_number)%nro_lc     = nro_end - nro_start + 1
    highden_table(eos_number)%ntt_lc     = ntt_end - ntt_start + 1
    highden_table(eos_number)%nye_lc     = nye_end - nye_start + 1

    nro_shape = nro_end - nro_start + 1
    ntt_shape = ntt_end - ntt_start + 1
    nye_shape = nye_end - nye_start + 1

    highden_table(eos_number)%nro_max_lc = nro_end
    highden_table(eos_number)%nro_min_lc = nro_start
    highden_table(eos_number)%ntt_max_lc = ntt_end
    highden_table(eos_number)%ntt_min_lc = ntt_start
    highden_table(eos_number)%nye_max_lc = nye_end

    highden_table(eos_number)%nye_min_lc = nye_start

    highden_table(eos_number)%romax_lc = 10.0**real(romax_lc_c_f,kind=rk)
    highden_table(eos_number)%ttmin_lc = 10.0**(real(ttmin_lc_c_f,kind=rk))
    highden_table(eos_number)%ttmax_lc = 10.0**(real(ttmax_lc_c_f,kind=rk))
    highden_table(eos_number)%yemin_lc = real(yemin_lc_c_f,kind=rk)
    highden_table(eos_number)%yemax_lc = real(yemax_lc_c_f,kind=rk)
#endif /* DYNAMIC_EOS */

    call c_f_pointer(lro_hd_c_f, highden_table(eos_number)%lro_hd, (/highden_table(eos_number)%nro/))

    call c_f_pointer(ltt_hd_c_f, highden_table(eos_number)%ltt_hd, (/highden_table(eos_number)%ntt/))

    call c_f_pointer(ye_hd_c_f,  highden_table(eos_number)%ye_hd,  (/highden_table(eos_number)%nye/))

#ifdef DYNAMIC_EOS
    call c_f_pointer(lro_hd_lc_c_f, highden_table(eos_number)%lro_hd_lc, (/highden_table(eos_number)%nro_lc/))

    call c_f_pointer(ltt_hd_lc_c_f, highden_table(eos_number)%ltt_hd_lc, (/highden_table(eos_number)%ntt_lc/))

    call c_f_pointer(ye_hd_lc_c_f,  highden_table(eos_number)%ye_hd_lc,  (/highden_table(eos_number)%nye_lc/))
#endif /* DYNAMIC_EOS */


    call c_f_pointer(lpr_hd_c_f, highden_table(eos_number)%lpr, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/))

    call c_f_pointer(led_hd_c_f, highden_table(eos_number)%led, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/)) 

    call c_f_pointer(cei_hd_c_f, highden_table(eos_number)%cei, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/)) 

    call c_f_pointer(cpi_hd_c_f, highden_table(eos_number)%cpi, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/)) 

    call c_f_pointer(cni_hd_c_f, highden_table(eos_number)%cni, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/))

    call c_f_pointer(cui_hd_c_f, highden_table(eos_number)%cui, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/))

    call c_f_pointer(gai_hd_c_f, highden_table(eos_number)%gai, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/))

    call c_f_pointer(sti_hd_c_f, highden_table(eos_number)%sti, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/))

    call c_f_pointer(xxn_hd_c_f, highden_table(eos_number)%xxn, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/))

    call c_f_pointer(xxp_hd_c_f, highden_table(eos_number)%xxp, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/))


    call c_f_pointer(xxa_hd_c_f, highden_table(eos_number)%xxa, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/))

    call c_f_pointer(xxh_hd_c_f, highden_table(eos_number)%xxh, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/))
    call c_f_pointer(xha_hd_c_f, highden_table(eos_number)%xha, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/))
    call c_f_pointer(xhz_hd_c_f, highden_table(eos_number)%xhz, (/nro_shape, &
                                                                  ntt_shape, &
                                                                  nye_shape/))

  end subroutine associate_hd_eos_c_to_fortran
#endif /* USE_MEMORY_MAPPING */

end module eos_data_type
!>
!> \verbatim
!> this module provides the lepton-radiation EoS table
!>
!>  Author: A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module data_lep_rad_table
  use precision
  use mpi_vertex, only : myproc
!  use eos_data_type, only : lept_radi_table

  implicit none
  private

  public load_lepton_radiation_table, be, bp

#ifdef CFC_TRANSPORT
  public tjttmin
#endif

!- boundaries of the ThJ table
  real(kind=rk), save :: tjroyemin,tjroyemax,tjttmin,tjttmax

!- dimensions of the ThJ table
  integer(kind=ik), parameter ::  ntstep=20, ltmin=5,   ltmax=12
  integer(kind=ik), parameter ::  ntt = (ltmax-ltmin)*ntstep+1
  integer(kind=ik), parameter ::  nrstep=20, lrmin=-10, lrmax=12
  integer(kind=ik), parameter ::  nro = (lrmax-lrmin)*nrstep+1

!- data of the ThJ table
  real(kind=rk), target :: ltt_ld(1:ntt),lro_ld(1:nro),      &
                let(1:ntt,1:nro),le(1:ntt,1:nro),       &
                lp(1:ntt,1:nro),ls(1:ntt,1:nro),gam(1:ntt,1:nro)


!- helpers for bicubic interpolation
  real(kind=rk), target :: be(4,4,ntt,nro)
  real(kind=rk), target :: bp(4,4,ntt,nro)

 contains
!>
!> \verbatim
!> load the lepton-radiation EoS table
!>
!>  Author: A. Marek
!> \endverbatim
!>
!> \param tbfile
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine load_lepton_radiation_table(tbfile)
    use precision
    use abort
    use eos_data_type, only : lept_radi_table
    use print_stdout_mod
    use mpi_vertex, only : myproc
    implicit none

    character*(*)    ::  tbfile
    integer(kind=ik) :: n1, n2, table_size
    real(kind=rk)    :: mem_local, mem_global

!> \todo mem_global should be known here
    mem_global = 0._rk
    mem_local  = 0._rk

    open (14, file = tbfile,  &
         form = 'unformatted')
    read (14) n1,n2,ltt_ld,lro_ld,let,le,lp,ls,gam
    close (14)

    table_size = 8 + (ntt+nro)*8 + ntt*nro*5*8

    if   ((n1 .ne. ntt) .or. (n2 .ne. nro))   then
       write (*,*) "Task ",myproc
       write (*,'(2(a,i4))') ' n1 = ',n1,' , n2 = ',n2
       raise_abort("load_lepton_radiation_table(): wrong boundaries")
    end if

    tjroyemin=10.0_rk**lro_ld(2)
    tjroyemax=10.0_rk**lro_ld(nro-1)
    tjttmin=10.0_rk**ltt_ld(2)
    tjttmax=10.0_rk**ltt_ld(ntt-1)

! set pointer to eos data

    lept_radi_table(1)%name       = tbfile
    lept_radi_table(1)%tablefile  = tbfile
    lept_radi_table(1)%eos_number = 0
    lept_radi_table(1)%nro        = nro
    lept_radi_table(1)%ntt        = ntt
    lept_radi_table(1)%nye        = 0
    lept_radi_table(1)%nspec      = 0
    lept_radi_table(1)%zz         = 0._rk
    lept_radi_table(1)%aa         = 0._rk


    lept_radi_table(1)%royemax      = tjroyemax
    lept_radi_table(1)%royemin      = tjroyemin
    lept_radi_table(1)%ttmax      = tjttmax
    lept_radi_table(1)%ttmin      = tjttmin     
    lept_radi_table(1)%yemax      = 0._rk
    lept_radi_table(1)%yemin      = 0._rk

    lept_radi_table(1)%lro_ld    => lro_ld
    lept_radi_table(1)%ltt_ld    => ltt_ld

    lept_radi_table(1)%let       => let
    lept_radi_table(1)%le        => le
    lept_radi_table(1)%lp        => lp 
    lept_radi_table(1)%ls        => ls
    lept_radi_table(1)%gam       => gam

    lept_radi_table(1)%dlttab    =  real((ntt-1),kind=rk)/(ltt_ld(ntt)-ltt_ld(1))
    lept_radi_table(1)%dldtab    =  real((nro-1),kind=rk)/(lro_ld(nro)-lro_ld(1))


    nullify(lept_radi_table(1)%be,lept_radi_table(1)%bp)

    nullify(lept_radi_table(1)%lro_hd, &
            lept_radi_table(1)%ltt_hd, &
            lept_radi_table(1)%ye_hd,  & 
            lept_radi_table(1)%lpr,    &
            lept_radi_table(1)%led,    &
            lept_radi_table(1)%cei,    &
            lept_radi_table(1)%cni,    &
            lept_radi_table(1)%cpi,    &
            lept_radi_table(1)%cui,    &
            lept_radi_table(1)%gai,    &
            lept_radi_table(1)%sti,    &
            lept_radi_table(1)%xxp,    &
            lept_radi_table(1)%xxn,    &
            lept_radi_table(1)%xxa,    &
            lept_radi_table(1)%xxh,    &
            lept_radi_table(1)%xhz,    &
#ifdef NSMIX
            lept_radi_table(1)%stn,    &
#endif
            lept_radi_table(1)%xha)

    nullify(lept_radi_table(1)%hetab,  &
            lept_radi_table(1)%xhtab,  &
            lept_radi_table(1)%xitab)


    call printit_taskX(0," ")
    call printit_taskX(0,"================================================================================")
    call printit_taskX(0," ")
    call printit_taskX(0,"  load_lepton_radiation_table> eos table --> ",trim(tbfile)," <-- successfully installed!")
    call printit_taskX(0,"      nroye = ", lept_radi_table(1)%nro)
    call printit_taskX(0,"      ntt   = ", lept_radi_table(1)%ntt)
    call printit_taskX(0," ")
    call printit_taskX(0,"      royemin = ", lept_radi_table(1)%royemin)
    call printit_taskX(0,"      royemax = ", lept_radi_table(1)%royemax)
    call printit_taskX(0,"      ttmin   = ", lept_radi_table(1)%ttmin)
    call printit_taskX(0,"      ttmax   = ", lept_radi_table(1)%ttmax)
    call printit_taskX(0," ")
    call printit_taskX(0,"   Size of table: [MB] ", table_size/1024/1024)
    call printit_taskX(0,"================================================================================")


    mem_local = table_size/1024._rk/1024._rk

      call print_memory_alloc(mem_local, mem_global, "load_lepton_radiation_table")
       
  end subroutine load_lepton_radiation_table


end module data_lep_rad_table


!>
!> \verbatim
!> this module provides the low-density NSE EoS table
!>
!>  Author: A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module data_low_density_nse
  use precision

  implicit none

  
  private
  public init_low_density_nse

  integer(kind=ik), target :: nt, nr, ne, nq
  
  real(kind=rk), allocatable, target :: &
                       tltab(:),dltab(:),yetab(:), eltab(:,:,:), &
                       xitab(:,:,:,:), hetab(:,:,:), xhtab(:,:,:)
  real(kind=rk), target :: zz, aa
  integer(kind=ik) :: n_hv

  contains
!>
!> \verbatim
!> load the low-density NSE table
!>
!>  Author: A. Marek
!> \endverbatim
!>
!> \param tbfile
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
    subroutine init_low_density_nse(tbfile)
      use precision
      use abort
      use nucparam, only :  n_fe56
      use eos_data_type
      use print_stdout_mod
      use mpi_vertex, only : myproc
      use configure
      use phycon, only : pc_kmev
      implicit none

      character*(*),intent(in) :: tbfile
      character(10)           :: tab_name
      character*5             :: nuc_name
      integer(kind=ik)        :: i, j, k, istat
      real(kind=rk)           :: mem_local, mem_global


!> \todo global memory should be known here 
      mem_local   = 0._rk
      mem_global = 0._rk
      
      if (config%use_network) then
! When using the Thielemann network, we can only switch to
! the NSE table, once NSE is FULLY established
         config%tkok = 0.8e10_rk*pc_kmev
      else
! When "flashing" is used to describe nuclear transmutations
! the NSE treshold temperature must be lower (T \approx 0.5 MeV):
! While full NSE may not be reached at this point, QSE already
! applies and is well approximated by the NSE composition.
! Retaining everything as 56Ni up to 8*10^9 K and neglecting
! dissociation into alphas can have fatal consequences.
         config%tkok = 0.5_rk
      endif

      if (config%low_den_nse_eos .gt. 2) then

         open (15, file = tbfile,form = 'unformatted')

         call printit_taskX(0," ")
         call printit_taskX(0,"***********************************************")
         call printit_taskX(0," ")
         call printit_taskX(0," Reading Low-den-NSE-table ",tbfile)
         call printit_taskX(0," Species:")

         if (config%low_den_nse_eos .eq. 4) then
            read (15) nt, nr, ne
         else
            read (15) nt, nr, ne, nq
         endif
         !      read(15) tab_name

         if (config%low_den_nse_eos .eq. 4) then
            allocate(tltab(nt),dltab(nr),yetab(ne),eltab(nt,nr,ne), &
                 hetab(nt,nr,ne),xhtab(nt,nr,ne),stat=istat)
            if(istat .ne. 0) then
               raise_abort("init_low_density_nse(): error in allocating arrays")
            endif

            mem_local = nt*8._rk + nr*8._rk + ne*8._rk +            &
                        nt*nr*ne*3*8._rk

            mem_local = mem_local / 1024._rk/1024._rk

            call print_memory_alloc(mem_local, mem_global, "eos_data: low_den_nse_eos = 4")
            
         else
            allocate(tltab(nt),dltab(nr),yetab(ne),eltab(nt,nr,ne), &
                 xitab(nt,nr,ne,nq),stat=istat)
            if(istat .ne. 0) then
               raise_abort("init_low_density_nse(): error in allocating arrays")
            endif

         
            mem_local = nt*8._rk + nr*8._rk + ne*8._rk +            &
                        nt*nr*ne*8._rk + nt*nr*ne*nq*8._rk

            mem_local = mem_local / 1024._rk/1024._rk   

            call print_memory_alloc(mem_local, mem_global, "eos_data: low_den_nse_eos")

         endif ! low_den_nse_eos = 4


         ! tab_name should be read
         tab_name = "dummy"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! this if clause should be removed, if eostab_kok_nse_17.i3e has
!!!!!!!!! the same structure as the other nse-tables !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (config%low_den_nse_eos .ne. 17) then
            do i=1,nq
              read (15) nuc_name
              call printit_taskX(0," ",nuc_name)
            enddo
            call printit_taskX(0," ")
         endif

         if (config%low_den_nse_eos .eq. 4) then
            read (15) tltab(:),dltab(:),yetab(:)
            read (15) zz,aa
            read (15) eltab(:,:,:),hetab(:,:,:),xhtab(:,:,:)
         else

            read (15) tltab(:),dltab(:),yetab(:)
            read (15) eltab(:,:,:)
            read (15) xitab(:,:,:,:)

            aa=0._rk
            zz=0._rk
         endif ! low_den_nse_eos = 4
         
         close (15)

         if (config%low_den_nse_eos .eq. 4) then
! the composition in the table is given for mn_54, charge neutrality
! is obtained using mn_54, but the result for X_heavy is put into
! Fe_56!

            call printit_taskX(0," ")
            call printit_taskX(0,"init_low_density_nse Warning: using Fe^56 !")
            call printit_taskX(0," ")

            n_hv=n_fe56
         else
            n_hv=0
         endif ! low_den_nse_eos = 4


         ! check whether mass fractions sum up to unit
         if (config%low_den_nse_eos .ne. 4) then
            do i=1,size(tltab)
               do j=1,size(dltab)
                  do k=1,size(yetab)
                     
                     if ( abs( sum(xitab(i,j,k,:)) -1.) .gt. 1e-6) then
                        call printit_taskX(0,"ERROR: Mass fractions do not sum to 1.0 in init_low_density_nse")
                        call printit_taskX(0,"SUM - 1. = ",sum(xitab(i,j,k,:)) -1.)
                        raise_abort("eos_data.f90(): sum error in massfractions")
                        
                     endif

                  enddo
               enddo
            enddo
         else
            call printit_taskX(0," ")
            call printit_taskX(0,"Cannot check sum of mass fractions in low_den_nse_eos = 4")
            call printit_taskX(0," ")
         endif
! set pointer to eos data

         lowden_table(1)%name   = tab_name
         lowden_table(1)%tablefile = tbfile
         lowden_table(1)%nro    = nr
         lowden_table(1)%ntt    = nt
         lowden_table(1)%nye    = ne

         lowden_table(1)%nspec  = nq
         lowden_table(1)%aa     = aa
         lowden_table(1)%zz     = zz
         
         lowden_table(1)%romax  = (dltab(nr))
         lowden_table(1)%romin  = (dltab(1))
         lowden_table(1)%ttmax  = (tltab(nt))
         lowden_table(1)%ttmin  = (tltab(1))
         lowden_table(1)%yemax  = (yetab(ne))
         lowden_table(1)%yemin  = (yetab(1))
         
         lowden_table(1)%lro_ld => dltab
         lowden_table(1)%ltt_ld => tltab
         lowden_table(1)%ye_ld  => yetab
         
         if (config%low_den_nse_eos .eq. 4) then
            lowden_table(1)%hetab  => hetab
            lowden_table(1)%xhtab  => xhtab
         
            nullify(lowden_table(1)%xitab)
         else
            nullify(lowden_table(1)%hetab)
            nullify(lowden_table(1)%xhtab)
         
            lowden_table(1)%xitab   =>  xitab
         endif ! low_den_nse_eos = 4     
         
         lowden_table(1)%eltab   =>  eltab
         lowden_table(1)%eltab_min = minval(eltab)
         
         nullify(lowden_table(1)%lro_hd, &
                 lowden_table(1)%ltt_hd, &
                 lowden_table(1)%ye_hd,  & 
                 lowden_table(1)%lpr   , &
                 lowden_table(1)%led,    &
                 lowden_table(1)%cei,    &
                 lowden_table(1)%cni,    &
                 lowden_table(1)%cpi,    &
                 lowden_table(1)%cui,    &
                 lowden_table(1)%gai,    &
                 lowden_table(1)%sti,    &
                 lowden_table(1)%xxp,    &
                 lowden_table(1)%xxn,    &
                 lowden_table(1)%xxa,    &
                 lowden_table(1)%xxh,    &
                 lowden_table(1)%xhz,    &
#ifdef NSMIX
                 lowden_table(1)%stn,    &
#endif
                 lowden_table(1)%xha)

         call printit_taskX(0," ")
         call printit_taskX(0," ",10._rk**lowden_table(1)%romin, " < rho < ",10._rk**lowden_table(1)%romax)
         call printit_taskX(0," ",10._rk**lowden_table(1)%ttmin, " < tem < ",10._rk**lowden_table(1)%ttmax)
         call printit_taskX(0,"Used for T > ", config%tkok)
         call printit_taskX(0," ",lowden_table(1)%yemin, " < ye < " ,lowden_table(1)%yemax)
         call printit_taskX(0," ")
         call printit_taskX(0," init_low_density_nse> read ", lowden_table(1)%nspec,"  species nucl.comp for thj")
         call printit_taskX(0," ")
         call printit_taskX(0,"*********************************************** ")
         call printit_taskX(0," ")

      endif ! low_den_nse_eos > 2
               
  end subroutine init_low_density_nse
end module data_low_density_nse


!>
!> \verbatim
!> this module loads the L&S-EoS table
!>
!>  Author: A. Marek
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module data_highden_eos

  use precision
  

  implicit none

  private
  public load_eos_table


#if !defined(USE_MEMORY_MAPPING) && !defined(DYNAMIC_EOS)
  integer(kind=ik), target :: nro, ntt, nye, timestamp

  real(kind=rk)    ,target :: romax,romin,ttmax,ttmin,yemax,yemin


  real(kind=rk_eos), allocatable, target :: lpr(:,:,:),led(:,:,:), &
                                            cei(:,:,:),cni(:,:,:), &
                                            cpi(:,:,:),cui(:,:,:), &
                                            gai(:,:,:),sti(:,:,:) 

  real(kind=rk_eos), allocatable, target :: lro(:), ltt(:), yei(:)



#ifdef NSMIX
  real(kind=rk_eos), allocatable, target :: stn(:,:,:)
#endif /* NSMIX */

!  Chemical composition
!     xxp: Protons
!     xxn: Neutrons
!     xxa: Alphas
!     xxh: Heavy
!     xhz: Z of heavy
!     xha: A of heavy

  real(kind=rk_eos), allocatable, target :: xxp(:,:,:), &
                                            xxn(:,:,:), &
                                            xxa(:,:,:), &
                                            xxh(:,:,:), &
                                            xhz(:,:,:), &
                                            xha(:,:,:)
#endif /* !USE_MEMORY_MAPPING && !DYNAMIC_EOS*/

contains

!> load the EoS table
!>
!> \author Original-Version: Wolfgang Keil, Markus Rampp, this version: A. Marek, MPA
!>
!> \param tbfile  file-name of table         
!>
  subroutine load_eos_table(tbfile, nrho_min, nrho_max, ntem_min, ntem_max, nye_min, nye_max)

    use precision
    use abort
#if defined(USE_MEMORY_MAPPING) || defined(DYNAMIC_EOS)
    use eos_c
#endif

    use phycon
    use nucparam, only: enullm
    use eos_data_type
    use print_stdout_mod
    use mpi_vertex, only : myproc
    implicit none

    character*(*)     :: tbfile
    character(72)    :: comment
    real(kind=rk)    :: enullm_in_table
    integer(kind=ik) :: file_error, table_size, i, j, k
    real(kind=rk)    :: mem_local, mem_global

#ifndef PROGRAM_remap
    integer(kind=ik), parameter :: eos_number= HIGH_DEN_EOS_FLAG
#else
    integer(kind=ik), parameter :: eos_number= 1
#endif

#ifdef NSMIX
    real(kind=rk) :: den(nro),tem(ntt)
#endif

    integer(kind=ik), intent(in) :: nrho_min, nrho_max, ntem_min, ntem_max, &
                                    nye_min, nye_max

#if defined(USE_MEMORY_MAPPING) || defined(DYNAMIC_EOS)
    integer(kind=ik) :: err
#endif


!----------------------------------------------------------------------
!     Dimensions of the EOS-table:                                     
!     nro = number of RO values                                        
!     ntt = number of TT values                                        
!     nye = number of YE values                                        
!     (calculation of log. rest value:                                 
!         INT(log(xxmax) - log(lastdecade)) * GPPD)                    
!
!     EOS-table:  eos_ls.i3e (REAL*4 - Version !!!!)           
!     lro         = log_10 of density in [g/cm^3]                      
!     ltt         = log_10 of temperature in [MeV]                     
!     yei         = electron fraction 
!                                 
!     lpr             = log_10 of pressure in [MeV/cm^3] 
!     led             = log_10 of energy density in [MeV/cm^3] 
!     cei,cni,cpi,cui = chemical potential of electrons, neutrons,      
!                       protons and Neutrinos in [MeV]
!                        cni is without restmass
!                        cei includes restmass 
!                        cpi is the chemical potential of protons 
!                            without restmass but - Q (see below)
!                        cui is never used in the code
!
!     gai             = dln(p)/dln(\rho)               
!     sti             = entropy
!     xxn,xxp,xxa,xxh = chemical composition
!     xhz,xha         
!----------------------------------------------------------------------

  highden_table(eos_number)%eos_energy_offs =.false.


  mem_local  = 0._rk
  mem_global = 0._rk

#if !defined(USE_MEMORY_MAPPING) && ! defined(DYNAMIC_EOS)

    open(52,file=tbfile,form = 'unformatted')
    read(52) timestamp, nro, ntt, nye
    print *,'EoS:',timestamp, nro, ntt, nye
    if (timestamp .ne. 20100907) then
        raise_abort("You have to use the EoS table with converted units")
    endif

    print *,'Allocating...'
    allocate(lro(nro), ltt(ntt),  yei(nye), lpr(nro,ntt,nye),        &
             led(nro,ntt,nye), cei(nro,ntt,nye), cni(nro,ntt,nye),   &
             cpi(nro,ntt,nye), cui(nro,ntt,nye), gai(nro,ntt,nye),   &
             sti(nro,ntt,nye), xxp(nro,ntt,nye), xxn(nro,ntt,nye),   &
             xxa(nro,ntt,nye), xxh(nro,ntt,nye), xhz(nro,ntt,nye),   &
             xha(nro,ntt,nye))
    print *,'Success!'

#ifdef DOUBLE_PRECISION_EOS
    table_size =  nro*8 + ntt*8 + nye*8 + nro*ntt*nye*8 * 14 
#else
    table_size =  nro*4 + ntt*4 + nye*4 + nro*ntt*nye*4 * 14 
#endif /* DOUBLE_PRECISION */
    
#ifdef NSMIX
    allocate(stn(nro,ntt,nye))

#ifdef DOUBLE_PRECISION_EOS
    table_size =  table_size + nro*ntt+nye*8
#else
    table_size =  table_size + nro*ntt+nye*4  
#endif /* DOUBLE_PRECISION */
#endif /* NSMIX */

    print *,'table_size =',table_size
!    rewind(52)
    close(52)
    open(52,file=tbfile,form = 'unformatted')
    print *,'rewind'

    read(52) timestamp, nro,ntt,nye,lro,ltt,yei,lpr,led,cei,cni,cpi,cui,gai, &
             sti,xxn,xxp,xxa,xxh,xhz,xha
    print *,lro
    print *,ltt
    print *,yei
    print *,'EoS: table read'
    read(52, iostat=file_error) enullm_in_table, comment
    print *,'EoS: enullm',enullm_in_table
    close(52)
    
    highden_table(eos_number)%timestamp = timestamp
    
! \todo decide what to do with different enullms

!      enullm = 8.8

  
     if (file_error /= 0) then

        ! enullm_in_table not read we assume the energy normalization
        ! to be 8.8 MeV then
        comment = "WARNING! Old format EOS, using default enullm = 8.8"
     else
        
        if (enullm_in_table .ne. enullm) then
           ! we have to do a renormalization

           highden_table(eos_number)%eos_energy_offs =.true.
           highden_table(eos_number)%enullm = enullm_in_table

        endif

     endif


     ! check whether mass fractions sum up to unit

     do i=1,size(lro)
        do j=1,size(ltt)
           do k=1,size(yei)

              if ( abs( (xxn(i,j,k) + xxp(i,j,k) + xxa(i,j,k) + xxh(i,j,k)) -1.) .gt. 1e-6) then
                 call printit_taskX(0,"ERROR: Mass fractions do not sum to 1.0 in load_eos_table")
                 call printit_taskX(0,"SUM - 1. = ",real(((xxn(i,j,k) + xxp(i,j,k) + xxa(i,j,k) + xxh(i,j,k)) - 1.0), kind=rk))
                 raise_abort("eos_data.f90(): sum error in massfractions")
                 
              endif

           enddo
        enddo
     enddo


#endif /* !USE_MEMORY_MAPPING && !DYNAMIC_EOS */

! nullify pointers which are not used in the HIGH DENSITY EOS
    nullify(highden_table(eos_number)%lro_ld, &
            highden_table(eos_number)%ltt_ld, &
            highden_table(eos_number)%ye_ld,  &
            highden_table(eos_number)%eltab,  &
            highden_table(eos_number)%hetab,  &
            highden_table(eos_number)%xhtab,  &
            highden_table(eos_number)%xitab)

    ! set variables to zero which are not used in the HIGH DENSITY EOS
    highden_table(eos_number)%nspec=0
    highden_table(eos_number)%n_hv=0
    highden_table(eos_number)%aa=0._rk
    highden_table(eos_number)%zz=0._rk

#if defined(USE_MEMORY_MAPPING)

     if(test_byte_order() .eq. 1) then
        err = load_mmap_hd_eos(tbfile // '_little_endian',1)
     else
        err = load_mmap_hd_eos(tbfile // '_big_endian',1)
     endif

     if (err == 1) then
        raise_abort("error in c-function load_mmap_hd_eos()")
     endif
 
     
     call associate_hd_eos_c_to_fortran(eos_number, nrho_min, nrho_max, &
                                        ntem_min, ntem_max, nye_min, nye_max)


#ifdef DOUBLE_PRECISION_EOS
    table_size =  (nrho_max-nrho_min+1)*8 + (ntem_max-ntem_min+1)*8 + (nye_max-nye_min+1)*8 + &
                  (nrho_max-nrho_min+1)*(ntem_max-ntem_min+1)*(nye_max-nye_min+1)*8 * 14 
#else
       table_size =  (nrho_max-nrho_min+1)*4 + (ntem_max-ntem_min+1)*4 + (nye_max-nye_min+1)*4 + &
                  (nrho_max-nrho_min+1)*(ntem_max-ntem_min+1)*(nye_max-nye_min+1)*4 * 14  
#endif /* DOUBLE_PRECISION_EOS */
     
#endif /* USE_MEMORY_MAPPING */

#if defined(DYNAMIC_EOS)

!    if (.not.(present(nrho_min))) then
       ! at first call set arbitarily chosen values for table dimensions
!       nrho_min = 1
!       nrho_max = 50
!       ntem_min = 1
!       ntem_max = 50
!       nye_min = 30
!       nye_max = 50
!    endif

      err = load_dynamic_hd_eos(tbfile, nrho_min, nrho_max, ntem_min, &
                                       ntem_max, nye_min, nye_max)
     
     if (err == 1) then
        raise_abort("error in c-function load_dynamic_hd_eos()")
     endif

     call associate_hd_eos_c_to_fortran(eos_number, nrho_min, nrho_max, ntem_min, &
                                        ntem_max, nye_min, nye_max)
     ! enullm_in_table not read we assume the energy normalization
     ! to be 8.8 MeV then
     comment = "WARNING! Old format EOS, using default enullm = 8.8"

#ifdef DOUBLE_PRECISION_EOS
    table_size =  (nrho_max-nrho_min+1)*8 + (ntem_max-ntem_min+1)*8 + (nye_max-nye_min+1)*8 + &
                  (nrho_max-nrho_min+1)*(ntem_max-ntem_min+1)*(nye_max-nye_min+1)*8 * 14 
#else
       table_size =  (nrho_max-nrho_min+1)*4 + (ntem_max-ntem_min+1)*4 + (nye_max-nye_min+1)*4 + &
                  (nrho_max-nrho_min+1)*(ntem_max-ntem_min+1)*(nye_max-nye_min+1)*4 * 14  
#endif /* DOUBLE_PRECISION_EOS */

#endif /* DYNAMIC_EOS */

#if !defined(USE_MEMORY_MAPPING) && !defined(DYNAMIC_EOS)

! define boundaries of the eos-table
!     (maybe it is better not to choose the outermost boundaries,      
!      because there the interpolation works slightly different)

!  if memory mapping is used, then this assignment was already done
!  in subroutine associate_c_to_fortan

    romax = 10._rk**real(lro(nro-1),kind=rk)
    romin = 10._rk**real(lro(2),kind=rk  )
    ttmax = 10._rk**real(ltt(ntt-1),kind=rk)
    ttmin = 10._rk**real(ltt(2),kind=rk  )
    yemax = real(yei(nye-1),kind=rk)
    yemin = real(yei(1),kind=rk)

    !  if memory mapping or dynamic_eos is used, then this assignment was already done
    !  in subroutine associate_c_to_fortan

    ! global values
    highden_table(eos_number)%nro        = nro
    highden_table(eos_number)%ntt        = ntt
    highden_table(eos_number)%nye        = nye

    highden_table(eos_number)%romax      = romax
    highden_table(eos_number)%romin      = romin
    highden_table(eos_number)%ttmax      = ttmax
    highden_table(eos_number)%ttmin      = ttmin     
    highden_table(eos_number)%yemax      = yemax
    highden_table(eos_number)%yemin      = yemin

!  if memory mapping is used, then this assignment was already done
!  in subroutine associate_c_to_fortan

    highden_table(eos_number)%lro_hd     => lro
    highden_table(eos_number)%ltt_hd     => ltt
    highden_table(eos_number)%ye_hd      => yei

    highden_table(eos_number)%lpr       => lpr
    highden_table(eos_number)%led       => led

    highden_table(eos_number)%cei       => cei
    highden_table(eos_number)%cni       => cni
    highden_table(eos_number)%cpi       => cpi
    highden_table(eos_number)%cui       => cui

    highden_table(eos_number)%gai       => gai
    highden_table(eos_number)%sti       => sti

    highden_table(eos_number)%xxp       => xxp
    highden_table(eos_number)%xxn       => xxn
    highden_table(eos_number)%xxa       => xxa

    highden_table(eos_number)%xxh       => xxh
    highden_table(eos_number)%xhz       => xhz
    highden_table(eos_number)%xha       => xha

#endif /* !definded(USE_MEMORY_MAPPING) && !defined(DYNAMIC_EOS) */ 

    if (highden_table(eos_number)%timestamp .ne. 20100907) then
        raise_abort("You have to use the EoS table with converted units")
    endif


! set pointer to eos data

    if (high_den_eos .eq. 1) then
       highden_table(eos_number)%name       = "LS-EoS"
    else if (high_den_eos .eq. 2) then
       highden_table(eos_number)%name       = "Wolff-EoS"
    else if (high_den_eos .eq. 3) then
       highden_table(eos_number)%name       = "SHEN-EoS"
    endif
    highden_table(eos_number)%tablefile     = tbfile
    highden_table(eos_number)%eos_number = eos_number



    ! almost done, now check the energy monotonicity
#ifndef DYNAMIC_EOS
    do i=1,highden_table(eos_number)%nro
       do k=1,highden_table(eos_number)%nye
          do j=2,highden_table(eos_number)%ntt

             if(     highden_table(eos_number)%led(i,j,k) &
                .le. highden_table(eos_number)%led(i,j-1,k)) then
!                print *,i,j,k,10**highden_table(eos_number)%lro_hd(i), &
!                              10**highden_table(eos_number)%ltt_hd(j), &
!                                  highden_table(eos_number)%ye_hd(k),  &
!                                  highden_table(eos_number)%led(i,j,k), &
!                                  highden_table(eos_number)%led(i,j-1,k)
!
!                raise_abort("loadtb_highden(): energy density not monoton")
             endif
          enddo
       enddo
    enddo
#else
    do i=1,highden_table(eos_number)%nro_lc
       do k=1,highden_table(eos_number)%nye_lc
          do j=2,highden_table(eos_number)%ntt_lc

         !    if(     highden_table(eos_number)%led(i,j,k) &
         !       .le. highden_table(eos_number)%led(i,j-1,k)) then
         !       print *,i,j,k,10**highden_table(eos_number)%lro_hd_lc(i), &
         !                     10**highden_table(eos_number)%ltt_hd_lc(j), &
         !                         highden_table(eos_number)%ye_hd_lc(k),  &
         !                         highden_table(eos_number)%led(i,j,k), &
         !                         highden_table(eos_number)%led(i,j-1,k)
         !       raise_abort("loadtb_highden(): energy density not monoton")
         !    endif
          enddo
       enddo
    enddo
#endif

    call printit_taskX(0," ")
    call printit_taskX(0,"================================================================================")
    call printit_taskX(0," ")
    call printit_taskX(0,"   LOADTB> eos table --> ",trim(tbfile),"  <-- successfully installed!")
    call printit_taskX(0,"   ",comment)
    call printit_taskX(0,"           enullm = [MeV]",highden_table(eos_number)%enullm)
    call printit_taskX(0,"           Global dimensions of the (whole) table:")
    call printit_taskX(0,"           nro = ",highden_table(eos_number)%nro)
    call printit_taskX(0,"           ntt = ",highden_table(eos_number)%ntt)
    call printit_taskX(0,"           nye = ",highden_table(eos_number)%nye)
    call printit_taskX(0," ")
    call printit_taskX(0,"           romin = [g/cc] ",  highden_table(eos_number)%romin)
    call printit_taskX(0,"           romax = [g/cc] ",  highden_table(eos_number)%romax)
    call printit_taskX(0,"           ttmin = [K]    ",  highden_table(eos_number)%ttmin)
    call printit_taskX(0,"           ttmax = [K]    ",  highden_table(eos_number)%ttmax)
    call printit_taskX(0,"           yemin = [1/by] ",  highden_table(eos_number)%yemin)
    call printit_taskX(0,"           yemax = [1/by] ",  highden_table(eos_number)%yemax)
    call printit_taskX(0," ")
#ifdef DYNAMIC_EOS
    call printit_taskX(0,"           Local (momentarily) dimensions of the loaded table:")
    call printit_taskX(0,"           nro = ",highden_table(eos_number)%nro_lc)
    call printit_taskX(0,"           from  ",highden_table(eos_number)%nro_min_lc)
    call printit_taskX(0,"           to    ",highden_table(eos_number)%nro_max_lc 
    call printit_taskX(0,"           ntt = ",highden_table(eos_number)%ntt_lc)
    call printit_taskX(0,"           from  ",highden_table(eos_number)%ntt_min_lc)
    call printit_taskX(0,"           to    ",highden_table(eos_number)%ntt_max_lc 
    call printit_taskX(0,"           ntt = ",highden_table(eos_number)%nye_lc)
    call printit_taskX(0,"           from  ",highden_table(eos_number)%nye_min_lc)
    call printit_taskX(0,"           to    ",highden_table(eos_number)%nye_max_lc 
    call printit_taskX(0," ")
    call printit_taskX(0,"           romin = [g/cc] ",  highden_table(eos_number)%romin_lc)
    call printit_taskX(0,"           romax = [g/cc] ",  highden_table(eos_number)%romax_lc)
    call printit_taskX(0,"           ttmin = [K]    ",  highden_table(eos_number)%ttmin_lc)
    call printit_taskX(0,"           ttmax = [K]    ",  highden_table(eos_number)%ttmax_lc)
    call printit_taskX(0,"           yemin = [1/by] ",  highden_table(eos_number)%yemin_lc)
    call printit_taskX(0,"           yemax = [1/by] ",  highden_table(eos_number)%yemax_lc)
    call printit_taskX(0," ")

#endif    
    call printit_taskX(0,"   Size of table: [MB] ",table_size/1024/1024)
    call printit_taskX(0,"================================================================================")


    mem_local = table_size/1024._rk/1024._rk

      call print_memory_alloc(mem_local, mem_global, "load_eos_table")

#ifdef NSMIX
!NUSTO QUICKFIX

#error "The units have changed, is this still correct?"

    tem(:) = 10._rk**real(ltt(:),kind=rk)
    den(:) = 10._rk**real(lro(:),kind=rk)
    do i=1,nye
       do j=1,ntt
          tema = tem(j)
          stn(:,j,i) = sti(:,j,i) + 1._rk/6._rk * ( tema * 1e-5_rk / wc_hc)**3 *      &
                       (3._rk * 7._rk*pc_pi**2/15._rk + (cui(:,j,i)/tema)**2) / den(:) * (pc_mb * 1e15_rk)
       enddo
    enddo

! pointer
    highden_table(eos_number)%stn => stn

#endif /* NSMIX */

  end subroutine load_eos_table

end module data_highden_eos

module data_wolffeos_low

  use precision
  

  implicit none

  private
  public loadtb_wolffeos_low


  integer(kind=ik), target :: nro, ntt, nye
  real(kind=rk),    target :: romax,romin,ttmax,ttmin,yemax,yemin

  real(kind=rk_eos), allocatable, target :: lpr(:,:,:),led(:,:,:), &
                                            cei(:,:,:),cni(:,:,:), &
                                            cpi(:,:,:),cui(:,:,:), &
                                            gai(:,:,:),sti(:,:,:) 

  real(kind=rk_eos), allocatable, target :: lro(:), ltt(:), yei(:)


#ifdef NSMIX
  real(kind=rk_eos), allocatable, target :: stn(:,:,:)
#endif /* NSMIX */

!  Chemical composition
!     xxp: Protons
!     xxn: Neutrons
!     xxa: Alphas
!     xxh: Heavy
!     xhz: Z of heavy
!     xha: A of heavy

  real(kind=rk_eos), allocatable, target :: xxp(:,:,:), &
                                            xxn(:,:,:), &
                                            xxa(:,:,:), &
                                            xxh(:,:,:), &
                                            xhz(:,:,:), &
                                            xha(:,:,:)

!----------------------------------------------------------------------
contains

  subroutine loadtb_wolffeos_low(tbfile)
                                                                      
!     loads table of Lattimer&Swesty EOS 
!                                                                      
!     Author       :  Original-Version :Wolfgang Keil, Markus Rampp
!                     this version : A. Marek, MPA
!
!     input:                                                           
!             tbfile = file-name of table         
!
!     important: it is assumed that the table entries are single 
!                precision and the numeric storage unit is 4 BYTES

    use precision
    use abort
    
    use phycon
    use nucparam, only: enullm
    use eos_data_type

    implicit none
    character(*)     :: tbfile
    character(72)    :: comment
    real(kind=rk)    :: enullm_in_table
    integer(kind=ik) :: file_error

    integer(kind=ik), parameter :: nnro=180, nntt=120, nnye=50
#ifndef PROGRAM_remap
    integer(kind=ik), parameter :: eos_number=INTERMEDIATE_EOS_FLAG
#else
    !> \todo What a mess with this numbers!
    integer(kind=ik) :: eos_number = 3
#endif

#ifdef NSMIX
    real(kind=rk)  :: den(nro), tem(ntt)
#endif

!----------------------------------------------------------------------
!     Dimensions of the EOS-table:                                     
!     nro = number of RO values                                        
!     ntt = number of TT values                                        
!     nye = number of YE values                                        
!     (calculation of log. rest value:                                 
!         INT(log(xxmax) - log(lastdecade)) * GPPD)                    
!
!     EOS-table:  eos_ls.i3e (REAL*4 - Version !!!!)           
!     lro         = log_10 of density in [g/cm^3]                      
!     ltt         = log_10 of temperature in [MeV]                     
!     yei         = electron fraction 
!                                 
!     lpr             = log_10 of pressure in [MeV/cm^3] 
!     led             = log_10 of energy density in [MeV/cm^3] 
!     cei,cni,cpi,cui = chemical potential of electrons, neutrons,      
!                       protons and Neutrinos in [MeV]
!                        cni is without restmass
!                        cei includes restmass 
!                        cpi is the chemical potential of protons 
!                            without restmass but - Q (see below)
!                        cui is never used in the code
!
!     gai             = dln(p)/dln(\rho)               
!     sti             = entropy
!     xxn,xxp,xxa,xxh = chemical composition
!     xhz,xha         
!----------------------------------------------------------------------
    if (intermediate_eos .eq. 4) then
       open(52,file=tbfile,form = 'unformatted')
       read(52) nro, ntt, nye

       ! check with grid points of table

       if (nro .ne. nnro .or. ntt .ne. nntt .or. nye .ne. nnye) then
          raise_abort("loadtb_wolffeos_low(): dimensions incorrect")
       endif

       allocate(lro(nro), ltt(ntt), yei(nye), lpr(nro,ntt,nye),        &
                led(nro,ntt,nye), cei(nro,ntt,nye), cni(nro,ntt,nye),  &
                cpi(nro,ntt,nye), cui(nro,ntt,nye), gai(nro,ntt,nye),  &
                sti(nro,ntt,nye), xxp(nro,ntt,nye), xxn(nro,ntt,nye),  &
                xxa(nro,ntt,nye), xxh(nro,ntt,nye),  xhz(nro,ntt,nye), &
                xha(nro,ntt,nye))

#ifdef NSMIX
       allocate(stn(nro,ntt,nye))
#endif

       rewind(52)
       print *,'rewind'
       read(52) nro,ntt,nye, lro,ltt,yei,lpr,led,cei,cni,cpi,cui,gai, &
                sti,xxn,xxp,xxa,xxh,xhz,xha
       read(52, iostat=file_error) enullm_in_table, comment
       close(52)

!
! FIXME
! decide what to do with different enullms
!

    !      enullm = 8.8
       comment = "WARNING! Using enullm = 8.8 regardles of EoS setting"
!     if (file_error /= 0) then
!       enullm = 8.8
!       comment = "WARNING! Old format EOS, using default enullm = 8.8"
!     else
!       enullm = enullm_in_table
!     endif



! check accuracy of mass-fractions ("warn of 'systematic' errors")
!      eps=2.0*EPSILON(xxn(1,1,1))
!      if (ANY(abs(1.0-xxn-xxp-xxa-xxh).gt.eps)) 
!     &     stop 'LOADTB> mass-fractions inconsistent'



      !---
      ! Convert into calculcation units
      !

      ! Atomic charge
      ! It actually makes quite a difference to multiply first and
      ! then interpolate or the other way round, e.g. "mapnuc"ed 
      ! mass ratios change by more than 10%
      ! xhz = xhz * xha

! define boundaries of the eos-table
!     (maybe it is better not to choose the outermost boundaries,      
!      because there the interpolation works slightly different)

       romax = 10._rk**real(lro(nro-1),kind=rk)
       romin = 10._rk**real(lro(2),kind=rk  )
       ttmax = 10._rk**real(ltt(ntt-1),kind=rk)
       ttmin = 10._rk**real(ltt(2),kind=rk  )
       yemax = real(yei(nye-1),kind=rk)
       yemin = real(yei(1),kind=rk)


! set pointer to eos data

       highden_table(eos_number)%name       = "WOLFFLOW-EoS"
       highden_table(eos_number)%tablefile  = tbfile
       highden_table(eos_number)%eos_number = eos_number
       highden_table(eos_number)%nro        = nro
       highden_table(eos_number)%ntt        = ntt
       highden_table(eos_number)%nye        = nye

       highden_table(eos_number)%romax      = romax
       highden_table(eos_number)%romin      = romin
       highden_table(eos_number)%ttmax      = ttmax
       highden_table(eos_number)%ttmin      = ttmin     
       highden_table(eos_number)%yemax      = yemax
       highden_table(eos_number)%yemin      = yemin


       nullify(highden_table(eos_number)%lro_ld, &
               highden_table(eos_number)%ltt_ld, &
               highden_table(eos_number)%ye_ld,  &
               highden_table(eos_number)%eltab,  &
               highden_table(eos_number)%hetab,  &
               highden_table(eos_number)%xhtab,  &
               highden_table(eos_number)%xitab)
       
       highden_table(eos_number)%nspec=0
       highden_table(eos_number)%n_hv=0
       highden_table(eos_number)%aa=0._rk
       highden_table(eos_number)%zz=0._rk
       
       highden_table(eos_number)%lro_hd       => lro
       highden_table(eos_number)%ltt_hd       => ltt
       highden_table(eos_number)%ye_hd        => yei
    
       highden_table(eos_number)%lpr       => lpr
       highden_table(eos_number)%led       => led
       
       highden_table(eos_number)%cei       => cei
       highden_table(eos_number)%cni       => cni
       highden_table(eos_number)%cpi       => cpi
       highden_table(eos_number)%cui       => cui
       
       highden_table(eos_number)%gai       => gai
       highden_table(eos_number)%sti       => sti
       
       highden_table(eos_number)%xxp       => xxp
       highden_table(eos_number)%xxn       => xxn
       highden_table(eos_number)%xxa       => xxa
    
       highden_table(eos_number)%xxh       => xxh
       highden_table(eos_number)%xhz       => xhz
       highden_table(eos_number)%xha       => xha
       

       write(*,'(80("="))')
       write(*,*)
       write(*,'("   LOADTB> eos table --> ",a," <-- successfully installed!")') trim(tbfile)
       write(*,'("           ",a,"")') comment
       write(*,'("           enullm = ",f5.2," MeV")') enullm
       write(*,'("           nro = ",i4)') nro
       write(*,'("           ntt = ",i4)') ntt
       write(*,'("           nye = ",i4)') nye
       write(*,*)
       write(*,'("           romin = ",1pe12.5," romax = ",1pe12.5," [g/cc]")') romin, romax
       write(*,'("           ttmin = ",1pe12.5," ttmax = ",1pe12.5," [K]")')    ttmin, ttmax
       write(*,'("           yemin = ",1pe12.5," yemax = ",1pe12.5," [1/by]")') yemin, yemax
       write(*,*)
       write(*,'(80("="))')
       
#ifdef NSMIX
!NUSTO QUICKFIX
       
#error "The units have changed, is this still correct?"

       tem(:) = 10._rk**real(ltt(:),kind=rk)
       den(:) = 10._rk**real(lro(:),kind=rk)
       do i=1,nye
          do j=1,ntt
             tema = tem(j)
             stn(:,j,i) = sti(:,j,i) + 1._rk/6._rk * ( tema * 1e-5_rk / wc_hc)**3 *      &
               (3._rk * 7._rk*pc_pi**2/15._rk + (cui(:,j,i)/tema)**2) / den(:) * (pc_mb * 1e15_rk)
          enddo
       enddo

       ! pointer
       highden_table(eos_number)%stn => stn

#endif

    endif ! intermediate_eos == 4

  end subroutine LOADTB_WOLFFEOS_LOW


end module data_wolffeos_low

