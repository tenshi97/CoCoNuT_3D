#ifdef DOUBLE_PRECISION_EOS
#define eos_cast(x) (x)
#else
#define eos_cast(x) real(x,kind=rk)
#endif

#undef DEBUG

!> \verbatim
!> this module provides the EoS interface for the high-density EoS in the
!> VERTEX code
!>
!>  Author: L. Huedepohl, A. Marek, originally M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
module high_density_eos

use precision
use abort
use error_screen

use eos_data_type, only : highden_table,                          &
#ifdef DYNAMIC_EOS
                          associate_hd_eos_c_to_fortran,          &
#endif /* DYNAMIC_EOS */
                          high_den_eos
implicit none

#ifdef DYNAMIC_EOS
real(kind=rk)    :: entry_lock_counter, exit_lock_counter
real(kind=rk)    :: reload_code_block_lock
real(kind=rk)    :: reload_count_entry_lock, reload_count_exit_lock
real(kind=rk)    :: check_reload_lock
#endif /* DYNAMIC_EOS */


private
public in_highdensity_eos,  &
#ifdef DYNAMIC_EOS
       acquire_sync_locks,  &
       release_sync_locks,  &
#endif /* DYNAMIC_EOS */
       highdensity_eos

contains
!> \verbatim
!>
!>  Determine whether grid point is in HIGH-density EoS    
!>
!> Author : L. Huedepohl
!>=======================================================================
!> \endverbatim
!>
!> \param rho density
!> \param tem temperature
!> \param Yee Y_e
!> \param eos_mode determine EoS table
!> \param is_in logical
!> 
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  function in_highdensity_eos(rho, tem, yee, eos_mode) result (intable)
    use precision
    implicit none
    
    integer(kind=ik), intent(in) :: eos_mode


    real(kind=rk),  intent(in) :: rho(:), tem(:), yee(:)
    logical           :: intable(size(rho))
#ifndef ONEMG_EOS
    intable = (rho < highden_table(eos_mode)%romax .and. rho >= highden_table(eos_mode)%romin) .and. &
              (tem < highden_table(eos_mode)%ttmax .and. tem >= highden_table(eos_mode)%ttmin) .and. &
              (yee < highden_table(eos_mode)%yemax .and. yee >= highden_table(eos_mode)%yemin)
#else /*  ONEMG_EOS */
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
     
    intable = (rho < highden_table(eos_mode)%romax                                                               .and. &
                ((tem >= highden_table(eos_mode)%tem_cut .and. rho >= highden_table(eos_mode)%romin_low_hd ) .or.      &
                 (tem  < highden_table(eos_mode)%tem_cut .and. rho >= highden_table(eos_mode)%romin_high_hd))  ) .and. &
              (tem < highden_table(eos_mode)%ttmax .and. tem >= highden_table(eos_mode)%ttmin)                   .and. &
              (yee < highden_table(eos_mode)%yemax .and. yee >= highden_table(eos_mode)%yemin)
#endif /* ONEMG_EOS */
  end function in_highdensity_eos

#ifdef DYNAMIC_EOS
!> \verbatim
!>
!>  Determine whether grid point is in the local cube of the highdensity EoS    
!>
!> Author : A. Marek
!>=======================================================================
!> \endverbatim
!>
!> \param rho density
!> \param tem temperature
!> \param Yee Y_e
!> \param eos_mode determine EoS table
!> \param is_in logical
!> 
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  function outof_local_highdensity_eos(rho, tem, yee, eos_mode) result (outoftable)
    use precision
    implicit none
    
    integer(kind=ik), intent(in) :: eos_mode

    real(kind=rk),  intent(in) :: rho(:), tem(:), yee(:)
    logical           :: outoftable(6)  ! 1 density lower bound
                                        ! 2 density upper bound
                                        ! 3 temperature lower bound
                                        ! 4 temperature upper bound
                                        ! 5 ye lower bound
                                        ! 6 ye upper bound

    outoftable(1) = .not.(minval(rho) > highden_table(eos_mode)%romin_lc)
    outoftable(2) = .not.(maxval(rho) < highden_table(eos_mode)%romax_lc)

    outoftable(3) = .not.(minval(tem) > highden_table(eos_mode)%ttmin_lc)
    outoftable(4) = .not.(maxval(tem) < highden_table(eos_mode)%ttmax_lc)

    outoftable(5) = .not.(minval(yee) > highden_table(eos_mode)%yemin_lc)
    outoftable(6) = .not.(maxval(yee) < highden_table(eos_mode)%yemax_lc)

  end function outof_local_highdensity_eos

#endif /* DYNAMIC_EOS */

! TODO: force/suggest inlining in a portable way
  
!> \verbatim
!>
!>  Extact from EoS-table for each quantity the eight points on a cube
!>  which are needed for interpolation
!>  THIS SUBROUTINE SHOULD BE INLINED
!>
!> Author : A. Marek
!>=======================================================================
!> \endverbatim
!>
!> \param nzones      number of zones on which to perform extaction
!> \param eos_mode    determine EoS table
!> \param xtrct_data  array of extracted points
!> \param i_rho0      index array for lower density index
!> \param i_rho1      index array for upper density index
!> \param i_tem0      index array for lower temperature index
!> \param i_tem1      index array for lower temperature index
!> \param i_ye0       index array for lower Y_e index
!> \param i_ye1       index array for upper Y_e index 
!> 
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine xtrct_intrpl_data_eos(nzones,eos_mode,xtrct_data,i_rho0,i_rho1,i_tem0,i_tem1,i_ye0,i_ye1 )
    use precision
    use phycon

    use eos_data_type, only : highden_table
    implicit none

    integer(kind=ik), intent(in) :: nzones,i_rho0(:),i_rho1(:),i_tem0(:), &
                                           i_tem1(:),i_ye0(:),i_ye1(:),eos_mode
    integer(kind=ik) :: i
    real(kind=rk_eos), intent(inout) :: xtrct_data(:,:,:)


! define some macros
#define DRU   1
#define ENER  2
#define MUE   3
#define MUN   4
#define MUP   5
#define MUNE  6
#define ADIA  7
#define ENTR  8
#define MASSN 9
#define MASSP 10
#define MASSA 11
#define MASSH 12
#define CHRG  13
#define ATMAS 14

    do i = 1, nzones

! the pressure values
       xtrct_data(i,DRU,1)=highden_table(eos_mode)%LPR(i_rho0(i),i_tem0(i),i_ye0(i)) 
       xtrct_data(i,DRU,2)=highden_table(eos_mode)%LPR(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,DRU,3)=highden_table(eos_mode)%LPR(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,DRU,4)=highden_table(eos_mode)%LPR(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,DRU,5)=highden_table(eos_mode)%LPR(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,DRU,6)=highden_table(eos_mode)%LPR(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,DRU,7)=highden_table(eos_mode)%LPR(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,DRU,8)=highden_table(eos_mode)%LPR(i_rho1(i),i_tem1(i),i_ye1(i))

! the energy values
       xtrct_data(i,ENER,1)=highden_table(eos_mode)%LED(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,ENER,2)=highden_table(eos_mode)%LED(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,ENER,3)=highden_table(eos_mode)%LED(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,ENER,4)=highden_table(eos_mode)%LED(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,ENER,5)=highden_table(eos_mode)%LED(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,ENER,6)=highden_table(eos_mode)%LED(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,ENER,7)=highden_table(eos_mode)%LED(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,ENER,8)=highden_table(eos_mode)%LED(i_rho1(i),i_tem1(i),i_ye1(i))


! the electron chemical potential 
       xtrct_data(i,MUE,1)=highden_table(eos_mode)%CEI(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MUE,2)=highden_table(eos_mode)%CEI(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MUE,3)=highden_table(eos_mode)%CEI(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MUE,4)=highden_table(eos_mode)%CEI(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,MUE,5)=highden_table(eos_mode)%CEI(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MUE,6)=highden_table(eos_mode)%CEI(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MUE,7)=highden_table(eos_mode)%CEI(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MUE,8)=highden_table(eos_mode)%CEI(i_rho1(i),i_tem1(i),i_ye1(i))

! the neutron  chemical potential 
       xtrct_data(i,MUN,1)=highden_table(eos_mode)%CNI(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MUN,2)=highden_table(eos_mode)%CNI(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MUN,3)=highden_table(eos_mode)%CNI(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MUN,4)=highden_table(eos_mode)%CNI(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,MUN,5)=highden_table(eos_mode)%CNI(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MUN,6)=highden_table(eos_mode)%CNI(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MUN,7)=highden_table(eos_mode)%CNI(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MUN,8)=highden_table(eos_mode)%CNI(i_rho1(i),i_tem1(i),i_ye1(i))


! the proton  chemical potential 
       xtrct_data(i,MUP,1)=highden_table(eos_mode)%CPI(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MUP,2)=highden_table(eos_mode)%CPI(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MUP,3)=highden_table(eos_mode)%CPI(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MUP,4)=highden_table(eos_mode)%CPI(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,MUP,5)=highden_table(eos_mode)%CPI(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MUP,6)=highden_table(eos_mode)%CPI(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MUP,7)=highden_table(eos_mode)%CPI(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MUP,8)=highden_table(eos_mode)%CPI(i_rho1(i),i_tem1(i),i_ye1(i))

! the neutrino chemical potential 
       xtrct_data(i,MUNE,1)=highden_table(eos_mode)%CUI(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MUNE,2)=highden_table(eos_mode)%CUI(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MUNE,3)=highden_table(eos_mode)%CUI(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MUNE,4)=highden_table(eos_mode)%CUI(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,MUNE,5)=highden_table(eos_mode)%CUI(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MUNE,6)=highden_table(eos_mode)%CUI(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MUNE,7)=highden_table(eos_mode)%CUI(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MUNE,8)=highden_table(eos_mode)%CUI(i_rho1(i),i_tem1(i),i_ye1(i))

       
! the adiabatic index
       xtrct_data(i,ADIA,1)=highden_table(eos_mode)%GAI(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,ADIA,2)=highden_table(eos_mode)%GAI(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,ADIA,3)=highden_table(eos_mode)%GAI(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,ADIA,4)=highden_table(eos_mode)%GAI(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,ADIA,5)=highden_table(eos_mode)%GAI(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,ADIA,6)=highden_table(eos_mode)%GAI(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,ADIA,7)=highden_table(eos_mode)%GAI(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,ADIA,8)=highden_table(eos_mode)%GAI(i_rho1(i),i_tem1(i),i_ye1(i))

! the entropy
       xtrct_data(i,ENTR,1)=highden_table(eos_mode)%STI(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,ENTR,2)=highden_table(eos_mode)%STI(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,ENTR,3)=highden_table(eos_mode)%STI(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,ENTR,4)=highden_table(eos_mode)%STI(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,ENTR,5)=highden_table(eos_mode)%STI(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,ENTR,6)=highden_table(eos_mode)%STI(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,ENTR,7)=highden_table(eos_mode)%STI(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,ENTR,8)=highden_table(eos_mode)%STI(i_rho1(i),i_tem1(i),i_ye1(i))

! the neutron mass fraction
       xtrct_data(i,MASSN,1)=highden_table(eos_mode)%XXN(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MASSN,2)=highden_table(eos_mode)%XXN(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MASSN,3)=highden_table(eos_mode)%XXN(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MASSN,4)=highden_table(eos_mode)%XXN(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,MASSN,5)=highden_table(eos_mode)%XXN(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MASSN,6)=highden_table(eos_mode)%XXN(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MASSN,7)=highden_table(eos_mode)%XXN(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MASSN,8)=highden_table(eos_mode)%XXN(i_rho1(i),i_tem1(i),i_ye1(i))
       

! the proton mass fraction
       xtrct_data(i,MASSP,1)=highden_table(eos_mode)%XXP(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MASSP,2)=highden_table(eos_mode)%XXP(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MASSP,3)=highden_table(eos_mode)%XXP(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MASSP,4)=highden_table(eos_mode)%XXP(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,MASSP,5)=highden_table(eos_mode)%XXP(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MASSP,6)=highden_table(eos_mode)%XXP(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MASSP,7)=highden_table(eos_mode)%XXP(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MASSP,8)=highden_table(eos_mode)%XXP(i_rho1(i),i_tem1(i),i_ye1(i))

! the alpha mass fraction
       xtrct_data(i,MASSA,1)=highden_table(eos_mode)%XXA(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MASSA,2)=highden_table(eos_mode)%XXA(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MASSA,3)=highden_table(eos_mode)%XXA(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MASSA,4)=highden_table(eos_mode)%XXA(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,MASSA,5)=highden_table(eos_mode)%XXA(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MASSA,6)=highden_table(eos_mode)%XXA(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MASSA,7)=highden_table(eos_mode)%XXA(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MASSA,8)=highden_table(eos_mode)%XXA(i_rho1(i),i_tem1(i),i_ye1(i))
       
! the heavy mass fraction
       xtrct_data(i,MASSH,1)=highden_table(eos_mode)%XXH(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MASSH,2)=highden_table(eos_mode)%XXH(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MASSH,3)=highden_table(eos_mode)%XXH(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MASSH,4)=highden_table(eos_mode)%XXH(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,MASSH,5)=highden_table(eos_mode)%XXH(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,MASSH,6)=highden_table(eos_mode)%XXH(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,MASSH,7)=highden_table(eos_mode)%XXH(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,MASSH,8)=highden_table(eos_mode)%XXH(i_rho1(i),i_tem1(i),i_ye1(i))

! the charge number
       xtrct_data(i,CHRG,1)=highden_table(eos_mode)%XHZ(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,CHRG,2)=highden_table(eos_mode)%XHZ(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,CHRG,3)=highden_table(eos_mode)%XHZ(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,CHRG,4)=highden_table(eos_mode)%XHZ(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,CHRG,5)=highden_table(eos_mode)%XHZ(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,CHRG,6)=highden_table(eos_mode)%XHZ(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,CHRG,7)=highden_table(eos_mode)%XHZ(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,CHRG,8)=highden_table(eos_mode)%XHZ(i_rho1(i),i_tem1(i),i_ye1(i))
       
! the atomic number
       xtrct_data(i,ATMAS,1)=highden_table(eos_mode)%XHA(i_rho0(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,ATMAS,2)=highden_table(eos_mode)%XHA(i_rho0(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,ATMAS,3)=highden_table(eos_mode)%XHA(i_rho0(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,ATMAS,4)=highden_table(eos_mode)%XHA(i_rho0(i),i_tem1(i),i_ye1(i))
       xtrct_data(i,ATMAS,5)=highden_table(eos_mode)%XHA(i_rho1(i),i_tem0(i),i_ye0(i))
       xtrct_data(i,ATMAS,6)=highden_table(eos_mode)%XHA(i_rho1(i),i_tem0(i),i_ye1(i))
       xtrct_data(i,ATMAS,7)=highden_table(eos_mode)%XHA(i_rho1(i),i_tem1(i),i_ye0(i))
       xtrct_data(i,ATMAS,8)=highden_table(eos_mode)%XHA(i_rho1(i),i_tem1(i),i_ye1(i))

    enddo

         
  end subroutine xtrct_intrpl_data_eos

#ifdef DYNAMIC_EOS

 subroutine determine_subcube_dimenions(border_flag, eos_mode, den_min, den_max, tem_min,    &
                                        tem_max, ye_min, ye_max, ede_min, ede_max, nrho_min, &
                                        nrho_max, ntem_min, ntem_max, nye_min, nye_max,mode, ler)

    use precision
    use print_stdout_mod
    implicit none

    
    logical, intent(in) :: ler
    logical, intent(in) :: border_flag(6) ! if true then this border is not in table
                                          ! 1 den_min, 2 den_max, 3 tem_min, 4 tem_max
                                          ! 5 ye_min, 6 ye_max

    integer(kind=ik), intent(in)  :: eos_mode, mode
    real(kind=rk), intent(in)     :: den_min, den_max, ye_min, ye_max     
                                 
    real(kind=rk), intent(in)     :: ede_min, ede_max
    real(kind=rk), intent(in)     :: tem_min, tem_max 
    real(kind=rk)                 :: temmin, temmax, ede_min_lg, ede_max_lg,denmin, denmax
    integer(kind=ik), intent(out) :: nrho_min, nrho_max, nye_min, nye_max

    integer(kind=ik)              :: tem_indices(5)
    integer(kind=ik), intent(out) :: ntem_min, ntem_max
    integer(kind=ik)              :: nro, ntt, nye,iteration_index
    real(kind=rk)                 :: frac_value


    denmin = log10(den_min)
    denmax = log10(den_max)

    ! assume that temperature changes might be 30 %

    temmax = log10(tem_max * 1.3_rk)
    temmin = log10(tem_min * 0.7_rk)

    ! CAREFUL THIS IT THE GLOBAL I.E. "MAXIMUM" GRID
    nro=highden_table(eos_mode)%nro
    ntt=highden_table(eos_mode)%ntt
    nye=highden_table(eos_mode)%nye

    if (border_flag(1)) then
       iteration_index = 1
       do while(denmin .ge. highden_table(eos_mode)%lro_hd(iteration_index) .and. iteration_index .le. nro)
          iteration_index = iteration_index + 1 
       end do
       
       nrho_min = iteration_index - 1
         
       if (nrho_min .lt. 1) then
          raise_abort("determine_subcub_dimensions(): nrho_min < 1 !")  
       endif
       if (nrho_min .gt. nro) then
          raise_abort("determine_subcub_dimensions(): nrho_min > nro !")  
       endif

       call printit_taskX(0," ")
       call printit_taskX(0," Reloading EoS-table...")
       call printit_taskX(0, Lower density border recomputed!")
       call printit_taskX(0," ")
 
    else
       nrho_min = highden_table(eos_mode)%nro_min_lc
    endif

    if (border_flag(2)) then
       iteration_index = highden_table(eos_mode)%nro
       
       do while(denmax .le. highden_table(eos_mode)%lro_hd(iteration_index) .and. iteration_index .ge. 1)
          iteration_index = iteration_index - 1 
       end do
       nrho_max = iteration_index + 1
       
       if (nrho_max .lt. 1 ) then
          raise_abort("determine_subcub_dimensions(): nrho_max < 1 !")  
       endif
       if (nrho_max .gt. highden_table(eos_mode)%nro ) then
          raise_abort("determine_subcub_dimensions(): nrho_max > nro !")  
       endif

       call printit_taskX(0," ")
       call printit_taskX(0," Reloading EoS-table...")
       call printit_taskX(0," Upper density border recomputed!")
       call printit_taskX(0," ")

    else
       nrho_max = highden_table(eos_mode)%nro_max_lc
    endif

    if (border_flag(5)) then
       iteration_index = 1
       do while(ye_min .ge. highden_table(eos_mode)%ye_hd(iteration_index) .and. iteration_index .le. nye)
          iteration_index = iteration_index + 1 
       end do
       nye_min = iteration_index - 1

       if (nye_min .lt. 1) then
          raise_abort("determine_subcub_dimensions(): nye_min < 1 !")  
       endif
       if (nye_min .gt. nye) then
          raise_abort("determine_subcub_dimensions(): nye_min > nye !")  
       endif

       call printit_taskX(0," ")
       call printit_taskX(0," Reloading EoS-table...")
       call printit_taskX(0," Lower Ye border recomputed!")
       call printit_taskX(0," ")

    else
       nye_min = highden_table(eos_mode)%nye_min_lc
    endif

    if (border_flag(6)) then
       iteration_index = highden_table(eos_mode)%nye
       do while(ye_max .le. highden_table(eos_mode)%ye_hd(iteration_index) .and. iteration_index .ge. 1)
          iteration_index = iteration_index - 1 
       end do
       nye_max = iteration_index + 1
       if (nye_max .lt. 1) then
          raise_abort("determine_subcub_dimensions(): nye_max < 1")  
       endif
       if (nye_max .gt. highden_table(eos_mode)%nye) then
          raise_abort("determine_subcub_dimensions(): nye_max > nye")  
       endif

       call printit_taskX(0," ")
       call printit_taskX(0," Reloading EoS-table...")
       call printit_taskX(0," Upper Ye border recomputed!")
       call printit_taskX(0," ")


    else
       nye_max = highden_table(eos_mode)%nye_max_lc
    endif

    ! now check boundaries for the temperature grid this is done such:
    ! determine index according to the temperature (this is correct in mode 1,
    !                                              but the temperature from the
    !                                              previous timestep for mode > 1)
    ! determine index according to energy          (this is correct if mode > 1)
    ! determine the index by taking the min/max of temperature and energy index

    ! the temperature grid

    if (border_flag(3)) then
       
       iteration_index = 1
       do while(temmax .ge. highden_table(eos_mode)%ltt_hd(iteration_index) .and. iteration_index .lt. ntt)
          iteration_index = iteration_index + 1 
       end do
       tem_indices(1) = iteration_index

       call printit_taskX(0," ")
       call printit_taskX(0," Reloading EoS-table...")
       call printit_taskX(0," Lower temperature border recomputed!")
       call printit_taskX(0," ")


    else
       tem_indices(1) = highden_table(eos_mode)%ntt_min_lc
    endif

       if (tem_indices(1) .lt. 1) then
          raise_abort("determine_subcub_dimensions(): tem_indices(1) < 1")  
       endif
       if (tem_indices(1) .gt. highden_table(eos_mode)%ntt) then
          raise_abort("determine_subcub_dimensions(): tem_indices(1) > ntt")  
       endif

    if (mode .ge. 2) then ! only check energy if energy inversion has to be done
       ! check lower energy bound

       ! start value
       iteration_index = tem_indices(1)
       ede_min_lg = log10(ede_min)
       !
       do while(ede_min_lg .le. highden_table(eos_mode)%led(nrho_min,iteration_index,nye_min) .and. iteration_index .gt. 1)
          iteration_index = iteration_index - 1 
       end do
       tem_indices(2) = iteration_index

       if (tem_indices(2) .lt. 1) then
          raise_abort("determine_subcub_dimensions(): tem_indices(2) < 1")  
       endif
       if (tem_indices(2) .gt. highden_table(eos_mode)%ntt) then
          raise_abort("determine_subcub_dimensions(): tem_indices(2) > ntt")  
       endif

       iteration_index = tem_indices(2)      
       !
       do while(ede_min_lg .le. highden_table(eos_mode)%led(nrho_min,iteration_index,nye_max) .and. iteration_index .gt. 1)
          iteration_index = iteration_index - 1 
       end do
       tem_indices(3) = iteration_index

       if (tem_indices(3) .lt. 1) then
          raise_abort("determine_subcub_dimensions(): tem_indices(3) < 1")  
       endif
       if (tem_indices(3) .gt. highden_table(eos_mode)%ntt) then
          raise_abort("determine_subcub_dimensions(): tem_indices(3) > ntt")  
       endif

       iteration_index = tem_indices(3)      

       do while(ede_min_lg .le. highden_table(eos_mode)%led(nrho_max,iteration_index,nye_min) .and. iteration_index .gt. 1)
          iteration_index = iteration_index - 1 
       end do
       tem_indices(4) = iteration_index

       if (tem_indices(4) .lt. 1) then
          raise_abort("determine_subcub_dimensions(): tem_indices(4) < 1")  
       endif
       if (tem_indices(4) .gt. highden_table(eos_mode)%ntt) then
          raise_abort("determine_subcub_dimensions(): tem_indices(4) > ntt")  
       endif

       iteration_index = tem_indices(4)      
       !
       do while(ede_min_lg .le. highden_table(eos_mode)%led(nrho_max,iteration_index,nye_max) .and. iteration_index .gt. 1)
          iteration_index = iteration_index - 1 
       end do

       tem_indices(5) = iteration_index

       if (tem_indices(5) .lt. 1) then
          raise_abort("determine_subcub_dimensions(): tem_indices(5) < 1")  
       endif
       if (tem_indices(5) .gt. highden_table(eos_mode)%ntt) then
          raise_abort("determine_subcub_dimensions(): tem_indices(5) > ntt")  
       endif

       ntem_min = max(1,minval(tem_indices(:)))
    else
       ntem_min = tem_indices(1)
    endif

    if (border_flag(4)) then
       ! ntem_max is not given externaly, thus compute
       iteration_index = highden_table(eos_mode)%ntt
       do while(temmax .le. highden_table(eos_mode)%ltt_hd(iteration_index) .and. iteration_index .gt. 1)
          iteration_index = iteration_index - 1 
       end do
       tem_indices(1) = iteration_index + 1

       call printit_taskX(0," ")
       call printit_taskX(0," Reloading EoS-table...")
       call printit_taskX(0," Upper temperature border recomputed!")
       call printit_taskX(0," ")

    else
       tem_indices(1) = highden_table(eos_mode)%ntt_max_lc
    endif
    
       if (tem_indices(1) .lt. 1) then
          raise_abort("determine_subcub_dimensions(): tem_indices(1) < 1")  
       endif
       if (tem_indices(1) .gt. highden_table(eos_mode)%ntt) then
          raise_abort("determine_subcub_dimensions(): tem_indices(1) > ntt")  
       endif


    if (mode .ge. 2) then
       ! check upper energy bound

       ! start value
       iteration_index = tem_indices(1)
       ede_max_lg = log10(ede_max)
       !
       do while(ede_max_lg .ge. highden_table(eos_mode)%led(nrho_min,iteration_index,nye_min) .and. iteration_index .lt. ntt)
          iteration_index = iteration_index + 1 
       end do

       tem_indices(2) = iteration_index

       if (tem_indices(2) .lt. 1) then
          raise_abort("determine_subcub_dimensions(): tem_indices(2) < 1")  
       endif
       if (tem_indices(2) .gt. highden_table(eos_mode)%ntt) then
          raise_abort("determine_subcub_dimensions(): tem_indices(2) > ntt")  
       endif

       iteration_index = tem_indices(2)      
       !
       do while(ede_max_lg .ge. highden_table(eos_mode)%led(nrho_min,iteration_index,nye_max) .and. iteration_index .lt. ntt)
          iteration_index = iteration_index + 1 
       end do

       tem_indices(3) = iteration_index

       if (tem_indices(3) .lt. 1) then
          raise_abort("determine_subcub_dimensions(): tem_indices(3) < 1")  
       endif
       if (tem_indices(3) .gt. highden_table(eos_mode)%ntt) then
          raise_abort("determine_subcub_dimensions(): tem_indices(3) > ntt")  
       endif

       iteration_index = tem_indices(3)      
       !
       do while(ede_max_lg .ge. highden_table(eos_mode)%led(nrho_max,iteration_index,nye_min) .and. iteration_index .lt. ntt)
          iteration_index = iteration_index + 1 
       end do

       tem_indices(4) = iteration_index

       if (tem_indices(4) .lt. 1) then
          raise_abort("determine_subcub_dimensions(): tem_indices(4) < 1")  
       endif
       if (tem_indices(4) .gt. highden_table(eos_mode)%ntt) then
          raise_abort("determine_subcub_dimensions(): tem_indices(4) > ntt")  
       endif

       iteration_index = tem_indices(4)      
       !
       do while(ede_max_lg .ge. highden_table(eos_mode)%led(nrho_max,iteration_index,nye_max) .and. iteration_index .lt. ntt)
          iteration_index = iteration_index + 1 
       end do

       tem_indices(5) = iteration_index

       if (tem_indices(5) .lt. 1) then
          raise_abort("determine_subcub_dimensions(): tem_indices(5) < 1")  
       endif
       if (tem_indices(5) .gt. highden_table(eos_mode)%ntt) then
          raise_abort("determine_subcub_dimensions(): tem_indices(5) > ntt")  
       endif

       ntem_max = min(highden_table(eos_mode)%ntt_lc,maxval(tem_indices(:)))
    else
       ntem_max = tem_indices(1)
    endif

    if (ler) then
       ! energy inversion did not work since eos table bounds were too small

       if (border_flag(3)) then
          ! lower bound was too high

          if (ntem_min .gt. 1) then
             ntem_min = ntem_min - 1
          endif

          if (ntem_min .eq. 1) then
             
             if (nrho_min .gt. 1) then
                nrho_min = nrho_min - 1
             endif

             if (nye_min .gt. 1) then
                nye_min = nye_min - 1
             endif
          endif

          call printit_taskX(0," ")
          call printit_taskX(0," Reloading EoS-table since energy inversion did not work..." )
          call printit_taskX(0," Lower energy border was too high!")
          call printit_taskX(0," ")


       endif
  
       if (border_flag(4)) then
          ! upper bound was too low

          ! first increase temperature grid
          if (ntem_max .lt. highden_table(eos_mode)%ntt) then
             ntem_max = ntem_max + (highden_table(eos_mode)%ntt - ntem_max)/2
          endif

          if (ntem_max .eq. highden_table(eos_mode)%ntt-1) then
             ntem_max = highden_table(eos_mode)%ntt
          endif

          if (ntem_max .eq. highden_table(eos_mode)%ntt) then
             if (nrho_max .lt. highden_table(eos_mode)%nro) then
                nrho_max = nrho_max + (highden_table(eos_mode)%nro - nrho_max)/2
             endif

             if (nrho_max .eq. highden_table(eos_mode)%nro-1) then
                nrho_max = highden_table(eos_mode)%nro
             endif

             if (nye_max .lt. highden_table(eos_mode)%nye) then
                nye_max = nye_max + (highden_table(eos_mode)%nye - nye_max)/2
             endif

             if (nye_max .eq. highden_table(eos_mode)%nye-1) then
                nye_max = highden_table(eos_mode)%nye
             endif
          endif

          call printit_taskX(0," ")
          call printit_taskX(0," Reloading EoS-table since energy inversion did not work..." )
          call printit_taskX(0," Upper energy border was too low!")
          call printit_taskX(0," ")

       endif

    endif

  end subroutine determine_subcube_dimenions

#endif /* DYNAMIC_EOS */


  subroutine compute_table_indices(eos_mode, den, tem, ye, irho0, irho1, &
                                   item0, item1, iye0, iye1, wr0, wr1,   &
                                   wt0, wt1, wy0, wy1)

    use precision
    use eos_data_type
    use specfun, only : fastlog10
    implicit none

    real(kind=rk), intent(in) :: den(:), tem(:), ye(:)

    integer(kind=ik), intent(in) :: eos_mode
    integer(kind=ik), dimension(size(den)), intent(out) :: irho0, irho1,    &
                                                           iye0,  iye1,     &
                                                           item0, item1
    real(kind=rk), dimension(size(den)), intent(out) :: wr1, wr0, wt0, wt1, &
                                                        wy0, wy1

    real(kind=rk), dimension(size(den)) :: rho0, tem0, ye0, logden, logtem, yetab
    real(kind=rk) :: rho_max, rho_min, tem_max, tem_min, ye_max, ye_min
    integer(kind=ik) :: i, nos, nro, ntt, nye

    integer(kind=ik) :: omp_get_thread_num,omp_get_num_threads,j

    nos=size(den)


#ifndef DYNAMIC_EOS
    ! we work on the "global" grid,
    nro=highden_table(eos_mode)%nro
    ntt=highden_table(eos_mode)%ntt
    nye=highden_table(eos_mode)%nye

    rho_max = highden_table(eos_mode)%lro_hd(nro)
    rho_min = highden_table(eos_mode)%lro_hd(1)

    tem_max = highden_table(eos_mode)%ltt_hd(ntt)
    tem_min = highden_table(eos_mode)%ltt_hd(1)

    ye_max  = highden_table(eos_mode)%ye_hd(nye)
    ye_min  = highden_table(eos_mode)%ye_hd(1)
#else /*DYNAMIC_EOS*/
    ! only part of the table is loaded, thus we work on the "local" grid

    nro=highden_table(eos_mode)%nro_lc
    ntt=highden_table(eos_mode)%ntt_lc
    nye=highden_table(eos_mode)%nye_lc

    rho_max = highden_table(eos_mode)%lro_hd_lc(nro)
    rho_min = highden_table(eos_mode)%lro_hd_lc(1)

    tem_max = highden_table(eos_mode)%ltt_hd_lc(ntt)
    tem_min = highden_table(eos_mode)%ltt_hd_lc(1)

    ye_max  = highden_table(eos_mode)%ye_hd_lc(nye)
    ye_min  = highden_table(eos_mode)%ye_hd_lc(1)
#endif /*DYNAMIC_EOS*/

    logden(:) = den(:)
    logtem(:) = tem(:)

    call fastlog10(logden,nos)
    call fastlog10(logtem,nos)

    do i = 1, nos
        yetab(i) = ye(i)

       rho0(i) = max( 1.0_rk, min(real((NRO-1),kind=rk),       &
            1.0_rk+real((NRO-1),kind=rk)*(logden(i)-rho_min)   &
            / (rho_max-rho_min)))
       
       tem0(i) = max( 1.0_rk, min( real((NTT-1),kind=rk),      &
            1.0_rk+real((NTT-1),kind=rk)*(logtem(i)-tem_min)   &
            / (tem_max-tem_min)))
       

       ye0 (i) = max( 1.0_rk, min( real((NYE-1),kind=rk),        &
            1.0_rk+real((NYE-1),kind=rk)*(    yetab(i) - ye_min) &
            / (ye_max-ye_min))) 
       
#ifdef IBM_OPT
       ! split some loops
    enddo

    do i = 1, nos
#endif /* IBM_OPT */
       
       irho0(i) = int( rho0(i) )        ! calculate index of rho grid
       item0(i) = int( tem0(i) )        ! calculate index of tem grid
       iye0(i) = int(  ye0(i) )         ! calculate index of ye  grid
       
       irho1(i) = irho0(i) + 1
       item1(i) = item0(i) + 1
       iye1(i) =  iye0(i)  + 1

#ifdef DYNAMIC_EOS       
       if (highden_table(eos_mode)%lro_hd_lc(irho0(i)) .gt. logden(i) ) then
          if (irho0(i) .gt. 1) then
             irho0(i) = irho0(i) - 1
          else
             raise_abort("compute_table_indices(): lower bi-section density > density!")  
          endif
       endif

       if (highden_table(eos_mode)%lro_hd_lc(irho1(i)) .lt. logden(i) ) then
          if (irho1(i) .lt. nro) then
             irho1(i) = irho1(i) + 1
          else
             raise_abort("compute_table_indices(): upper bi-section density < density!")
          endif
       endif

       if (highden_table(eos_mode)%ltt_hd_lc(item0(i)) .gt. logtem(i) ) then
          if (item0(i) .gt. 1) then
             item0(i) = item0(i) - 1
          else
             raise_abort("compute_table_indices(): lower bi-section temperature > temperature!")
          endif
       endif

       if (highden_table(eos_mode)%ltt_hd_lc(item1(i)) .lt. logtem(i) ) then
          if (item1(i) .lt. ntt) then
             item1(i) = item1(i) + 1
          else
             raise_abort("compute_table_indices(): upper bi-section temperature < temperature!")
          endif
       endif

       if (highden_table(eos_mode)%ye_hd_lc(iye0(i)) .gt. (yetab(i))) then
          if (iye0(i) .gt. 1) then
             iye0(i) = iye0(i) - 1
          else
             raise_abort("compute_table_indices(): lower bi-section Ye > Ye!")
          endif
       endif

       if (highden_table(eos_mode)%ye_hd_lc(iye1(i)) .lt. (yetab(i))) then
          if (iye1(i) .lt. nye) then
             iye1(i) = iye1(i) + 1
          else
             raise_abort("compute_table_indices(): upper bi-section Ye < Ye!")
          endif
       endif
#endif /* DYNAMIC_EOS */

#ifdef IBM_OPT
       ! split some loops
    enddo
    
    do i = 1, nos
#endif /* IBM_OPT */

       ! compute the interpolation wheight factors
#ifndef DYNAMIC_EOS
       wr1(i) = (logden(i)-highden_table(eos_mode)%lro_hd(irho0(i))) / &
                (highden_table(eos_mode)%lro_hd(irho1(i)) - highden_table(eos_mode)%lro_hd(irho0(i)))
       wr0(i) = 1.0_rk - wr1(i)

       wt1(i) = (logtem(i)-highden_table(eos_mode)%ltt_hd(item0(i))) / &
                (highden_table(eos_mode)%ltt_hd(item1(i)) - highden_table(eos_mode)%ltt_hd(item0(i)))
       wt0(i) = 1.0_rk - wt1(i)

       wy1(i) = (yetab(i)-highden_table(eos_mode)%ye_hd(iye0(i))) / &
                (highden_table(eos_mode)%ye_hd(iye1(i)) - highden_table(eos_mode)%ye_hd(iye0(i)))
       wy0(i) = 1.0_rk - wy1(i)
#else
       wr1(i) = (logden(i)-highden_table(eos_mode)%lro_hd_lc(irho0(i))) / &
                (highden_table(eos_mode)%lro_hd_lc(irho1(i)) - highden_table(eos_mode)%lro_hd_lc(irho0(i)))
       wr0(i) = 1.0_rk - wr1(i)

       wt1(i) = (logtem(i)-highden_table(eos_mode)%ltt_hd_lc(item0(i))) / &
                (highden_table(eos_mode)%ltt_hd_lc(item1(i)) - highden_table(eos_mode)%ltt_hd_lc(item0(i)))
       wt0(i) = 1.0_rk - wt1(i)

       wy1(i) = (yetab(i)-highden_table(eos_mode)%ye_hd_lc(iye0(i))) / &
                (highden_table(eos_mode)%ye_hd_lc(iye1(i)) - highden_table(eos_mode)%ye_hd_lc(iye0(i)))
       wy0(i) = 1.0_rk - wy1(i)
#endif /* DYNAMIC_EOS */

    enddo

#ifdef DYNAMIC_EOS
    do i = 1, nos

#if DEBUG
       if (highden_table(eos_mode)%lro_hd_lc(irho0(i)) .gt. logden(i) ) then
          raise_abort("compute_table_indices(): lower bi-section density > density!")
       endif

       if (highden_table(eos_mode)%lro_hd_lc(irho1(i)) .lt. logden(i) ) then
          raise_abort("compute_table_indices(): upper bi-section density < density!")
       endif

       if (highden_table(eos_mode)%ltt_hd_lc(item0(i)) .gt. logtem(i) ) then
          raise_abort("compute_table_indices(): lower bi-section temperature > temperature!")
       endif

       if (highden_table(eos_mode)%ltt_hd_lc(item1(i)) .lt. logtem(i) ) then
          raise_abort("compute_table_indices(): upper bi-section temperature < temperature!")
       endif

       if (highden_table(eos_mode)%ye_hd_lc(iye0(i)) .gt. (yetab(i))) then
          raise_abort("compute_table_indices(): lower bi-section Ye > Ye!")
       endif

       if (highden_table(eos_mode)%ye_hd_lc(iye1(i)) .lt. (yetab(i))) then
          raise_abort("compute_table_indices(): upper bi-section Ye < Ye!")
       endif
#endif /* DEBUG */
    enddo

#endif /*DYNAMIC_EOS*/

  end subroutine compute_table_indices

#ifdef DYNAMIC_EOS
  subroutine acquire_sync_locks

    implicit none

    call omp_init_nest_lock(reload_count_entry_lock)
    call omp_init_nest_lock (check_reload_lock)
    call omp_init_nest_lock(reload_count_exit_lock)
    call omp_init_nest_lock(reload_code_block_lock)
    call omp_init_nest_lock(entry_lock_counter)
    call omp_init_nest_lock(exit_lock_counter)

  end subroutine acquire_sync_locks

  subroutine release_sync_locks

    implicit none
     
    call omp_destroy_nest_lock(reload_count_entry_lock)
    call omp_destroy_nest_lock(check_reload_lock)   
    call omp_destroy_nest_lock(reload_code_block_lock)
    call omp_destroy_nest_lock(reload_count_exit_lock)
    call omp_destroy_nest_lock(entry_lock_counter)    
    call omp_destroy_nest_lock(exit_lock_counter)

  end subroutine release_sync_locks
#endif /* DYNAMIC_EOS */


!> \verbatim
!>-----------------------------------------
!>
!> Author            : Wolfgang Keil, Markus Rampp (MPA)
!>                     Vectorized bisection due to Maximilian Ruffert 
!>                     Optimized version for IBM Power6 by A. Marek
!> 
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
!> Purpose           :
!>                     
!>                                                                     
!>               E q u a t i o n   o f   s t a t e                     
!>                                                                     
!>     Nuclear eos from tables, Lattimer & Swesty                      
!>                                                                     
!>     The state variables are calculated at nos grid points         
!>                                                                   
!>                                                                   
!>     input :                                                       
!>                                                                   
!>     rho   =    mass density   (g / cm^3)                        
!>     tem   =    temperature      (K)                             
!>     yee   =    electron fraction(1/by)                            
!>     input =   1  <==>  input  :  RHO , T                          
!>               2  <==>  input  :  RHO , E                          
!>                                        T (as guess value)         
!>     output:                                                       
!>                                                                   
!>     pr    =    pressure         (erg/cm^3)                           
!>     ce    =    chemical potential of electrons (MeV, including 
!>                                                      restmass)
!>     cn,cp =    chemical potential of neutrons,protons (MeV, without 
!>                                                             restmass)
!>     cu    =    chemical potential of neutrinos 
!>                                                                      
!>     input =   1  <==>  output  :  ED (erg/cm^3) 
!>               2  <==>  output  :  T  (K)
!>                                   ED (unchanged)
!>
!> IMPORTANT: the bisection works only for E(T) being MONOTONICALLY
!>           INCREASING, this property has to be guaranteed prior to call
!>            there is no check in these routines 
!>
!>-----------------------------------------   
!> \endverbatim
  subroutine highdensity_eos(eos_mode, rho, tem, xi, xhrep, za, ed, pr, ga, &
                                       st, cu, ce, cn, cp, mode, nsemode, ler)
    use eostools
    use precision
    use abort

    use nucparam, only : n_n,n_p,n_d,n_he3,n_he4, enullm

    use phycon, only : pc_mb, pc_meverg, pc_mevk, pc_kmev
#ifdef DYNAMIC_EOS
    use data_highden_eos, only : load_eos_table
    use  eos_data_type, only : thread_counter, reload_counter
#endif
    use specfun, only : fastlog10

    use mpi_vertex, only : myproc
    use configure
    use print_stdout_mod
    implicit none

      integer(kind=ik), intent(in)    :: eos_mode                    
      integer(kind=ik), intent(in)    :: mode, nsemode
      real(kind=rk),    intent(in)    :: rho(:)
      real(kind=rk),    intent(out)   :: pr(:), ga(:), za(:,:),  &
                                         xhrep(:), cu(:), ce(:), cn(:), &
                                         cp(:)

      real(kind=rk),    intent(inout) :: tem(:), xi(:,:), ed(:), st(:)
      logical, intent(out)   :: ler

      !> \var wr_0  interpolation weight (low) in density
      !> \var wr_1  interpolation weight (high) in density
      !> \var wt_0  interpolation weight (low) in temperature
      !> \var wt_1  interpolation weight (high) in temperature
      !> \var wy_0  interpolation weight (low) in Y_e
      !> \var wy_1  interpolation weight (high) in Y_e
      real(kind=rk), dimension(size(rho)) :: rho0, tem0, ye0,           &
                                    wr1, wr0, wy1, wy0, wt1, wt0,       &
                                    scl, sch, fun
      real(kind=rk), dimension(size(rho)) :: xn, xp, xa, xh
      real(kind=rk), dimension(size(rho)) :: a, z_a
      integer(kind=ik), dimension(size(rho)) :: irho0, irho1,    &
                                       iye0,  iye1,     &
                                       item0, item1,    &
                                       ittm

      real(kind=rk), dimension(size(rho)) :: isect ! isect is now a float
      real(kind=rk) :: sumisec
      real(kind=rk) scratch(4), scm, summ, z
      integer(kind=ik) nos,nuc,i,kk, err, omp_get_thread_num, &
                       omp_get_num_threads
#ifdef DYNAMIC_EOS
      logical :: out_subcube(6)
#endif

#ifdef DYNAMIC_EOS
      integer(kind=ik) :: nrho_min, nrho_max, ntem_min, ntem_max, nye_min, &
                          nye_max
      integer(kind=ik) :: retry_reload_counter 
#endif /* DYNAMIC_EOS */


!> new data array that will contain for every zone, for the 14 interpolation
!> quantities of the eos (ugly that this is hard coded - improve! ) the
!> eight interpolation points for tri-linear interpolation with

!> intrpl_data(:,:,1) = rho_low ,tem_low, ye_low
!> intrpl_data(:,:,2) = rho_low ,tem_low, ye_high 
!> intrpl_data(:,:,3) = rho_low ,tem_high, ye_low
!> intrpl_data(:,:,4) = rho_low ,tem_high, ye_high
!> intrpl_data(:,:,5) = rho_thigh ,tem_low, ye_low
!> intrpl_data(:,:,6) = rho_thigh ,tem_low, ye_high 
!> intrpl_data(:,:,7) = rho_thigh ,tem_high, ye_low
!> intrpl_data(:,:,8) = rho_thigh ,tem_high, ye_high

! a hash would be nice here
       
      real(kind=rk_eos), dimension(size(rho),14,8) :: intrpl_data

      integer(kind=ik) :: nro, ntt, nye

#ifdef DYNAMIC_EOS
    interface
       integer(kind=c_int) FUNCTION flush_table_subcube() bind(C)
         use, intrinsic :: iso_c_binding
         implicit none
 
       end FUNCTION flush_table_subcube
    end interface
    
    retry_reload_counter = 0
#endif /* DYNAMIC_EOS */

#ifndef DYNAMIC_EOS
    ! use the global, i.e. maximum possible grid dimensions
      nro=highden_table(eos_mode)%nro
      ntt=highden_table(eos_mode)%ntt
      nye=highden_table(eos_mode)%nye
#else /* DYNAMIC_EOS */
    ! use the local, i.e. minimal possible grid dimensions
      nro=highden_table(eos_mode)%nro_lc
      ntt=highden_table(eos_mode)%ntt_lc
      nye=highden_table(eos_mode)%nye_lc
      out_subcube(:)   =.false.
#endif /* DYNAMIC_EOS */

      nos=size(rho)
      nuc=size(xi, dim=2)

      ler=.false.

      if (highden_table(eos_mode)%eos_energy_offs) then

         ! High density EoS is tabulated with a different
         ! energy offset than the "normal" 8.8 MeV
         ! we have to correct for this here

         ed(:) = ed(:) +                                   &
                 (highden_table(eos_mode)%enullm - enullm) &
                 * rho(:) * pc_meverg / pc_mb

      endif

#ifdef DYNAMIC_EOS

      ! check whether the given points of rho, tem, and, ye are
      ! still located in the currently loaded sub-cube of the EoS

      out_subcube(:) = outof_local_highdensity_eos(rho(:), tem(:), xi(:,nuc), eos_mode)

      if (any(out_subcube)) then
         ! determine how many thread would like to reload the EoS table
         call omp_set_nest_lock  (reload_count_entry_lock)
         reload_counter = reload_counter + 1
         call omp_unset_nest_lock(reload_count_entry_lock)
         

         do while (reload_counter .gt. 1)
            ! maybe one thread increased table enough, thus check again to
            ! be sure
            call omp_set_nest_lock(check_reload_lock)
            if (.not.(any(outof_local_highdensity_eos(rho(:), tem(:), &
                          xi(:,nuc), eos_mode) )) ) then
               reload_counter = reload_counter - 1

            endif
            call omp_unset_nest_lock(check_reload_lock)        

         enddo
      endif

      out_subcube(:) = outof_local_highdensity_eos(rho(:), tem(:), xi(:,nuc), eos_mode)

1111  continue

      if (any(out_subcube)) then

         call omp_set_nest_lock(reload_code_block_lock)
         do while (thread_counter .gt. 0) 

         enddo

         call determine_subcube_dimenions(out_subcube, eos_mode, minval(rho(:)),maxval(rho(:)),      &
                                          minval(tem(:)), maxval(tem(:)), minval(xi(:,nuc)),         &
                                          maxval(xi(:,nuc)), minval(ed(:)), maxval(ed(:)), nrho_min, &
                                          nrho_max, ntem_min, ntem_max, nye_min, nye_max, mode, ler)

#ifdef DEBUG
         write (*,*) "NEW EOS indices: ",nrho_min, nrho_max, ntem_min, ntem_max,  &
                                     nye_min, nye_max
#endif /* DYNAMIC_EOS */
  
         ! deallocate previos subcube
         err = flush_table_subcube()

         if (err .ne. 0) then
            raise_abort("highdensity():error when freeing sub-table memory")
         endif

         ! reload the table
#ifdef LATTIMER_EOS
         call load_eos_table("./tables/eos_ls.i3e",nrho_min, nrho_max,    &
                             ntem_min, ntem_max, nye_min, nye_max )
#endif

#ifdef WOLFF_EOS
         call load_eos_table("./tables/eos_wolff.i3e",nrho_min, nrho_max,    &
                             ntem_min, ntem_max, nye_min, nye_max )
#endif

#ifdef SHEN_EOS
         call load_eos_table("./tables/eos_shen.i3e",nrho_min, nrho_max,    &
                             ntem_min, ntem_max, nye_min, nye_max )
#endif
         call associate_hd_eos_c_to_fortran(eos_mode, nrho_min, nrho_max, &
                                            ntem_min, ntem_max, nye_min,  &
                                            nye_max)
         call omp_unset_nest_lock(reload_code_block_lock)

      endif
   
      ! a "resize" .ie. the loaded table is bigger than the needed one has
      ! still to be implemented

      ! reloading finished of actual thread => reduce counter 
      call omp_set_nest_lock  (reload_count_exit_lock)
      if (reload_counter .gt. 0) then
         reload_counter = reload_counter - 1
      endif
      call omp_unset_nest_lock(reload_count_exit_lock)
#endif /* DYNAMIC_EOS */
 
#ifdef DYNAMIC_EOS
      ! counter-variable: update number of threads actually in the EoS
      ! evaluation
      call omp_set_nest_lock(entry_lock_counter)
      thread_counter = thread_counter + 1
      call omp_unset_nest_lock(entry_lock_counter)
#endif /* DYNAMIC_EOS */  
  

      ! this has to be called for every mode anyway
      ! the only difference is that mode = 1 uses irho0, irho1, item0, item1,
      ! _AND_ iye0, iye1, and the other modes just irho0, irho1, item0, item1
      call compute_table_indices(eos_mode, rho(:), tem(:), xi(:,nuc),      &
                                 irho0(:), irho1(:), item0(:), item1(:),   &
                                 iye0(:), iye1(:), wr0(:), wr1(:), wt0(:), &
                                 wt1(:), wy0(:), wy1(:) )
#if defined(DYNAMIC_EOS) && defined(DEBUG)
      ! check energy bounds

      if (mode .ge. 2) then

      fun(:) = ed(:)

      call fastlog10(fun, size(fun) )

      do i=1,nos
         
         if(eos_cast(highden_table(eos_mode)%LED(irho0(i),1,iye0(i))) .gt. &
             fun(i) ) then
            
            write (*,*) i,eos_cast(highden_table(eos_mode)%LED(irho0(i),1,iye0(i))), &
                    fun(i)
            raise_abort("A")
         endif

         if(eos_cast(highden_table(eos_mode)%LED(irho0(i),1,iye1(i))) .gt. &
             fun(i) ) then
            
            write (*,*) i,eos_cast(highden_table(eos_mode)%LED(irho0(i),1,iye1(i))), &
                    fun(i)
            raise_abort("B")
         endif

         if(eos_cast(highden_table(eos_mode)%LED(irho1(i),1,iye0(i))) .gt. &
             fun(i) ) then
            
            write (*,*) i,eos_cast(highden_table(eos_mode)%LED(irho1(i),1,iye0(i))), &
                    fun(i)
            raise_abort("C")
         endif

         if(eos_cast(highden_table(eos_mode)%LED(irho1(i),1,iye1(i))) .gt. &
             fun(i) ) then
            
            write (*,*) i,eos_cast(highden_table(eos_mode)%LED(irho1(i),1,iye1(i))), &
                    fun(i)
            raise_abort("D")
         endif

         ntt = highden_table(eos_mode)%ntt_lc  ! this can either be to global index or local one. CAREFUL

         
         if(eos_cast(highden_table(eos_mode)%LED(irho0(i),ntt,iye0(i))) .gt. &
             fun(i) ) then
            
            write (*,*) i,eos_cast(highden_table(eos_mode)%LED(irho0(i),ntt,iye0(i))), &
                    fun(i)
            raise_abort("A1")
         endif

         if(eos_cast(highden_table(eos_mode)%LED(irho0(i),ntt,iye1(i))) .gt. &
             fun(i) ) then
            
            write (*,*) i,eos_cast(highden_table(eos_mode)%LED(irho0(i),ntt,iye1(i))), &
                    fun(i)
            raise_abort("B1")
         endif

         if(eos_cast(highden_table(eos_mode)%LED(irho1(i),ntt,iye0(i))) .gt. &
             fun(i) ) then
            
            write (*,*) i,eos_cast(highden_table(eos_mode)%LED(irho1(i),ntt,iye0(i))), &
                    fun(i)
            raise_abort("C1")
         endif

         if(eos_cast(highden_table(eos_mode)%LED(irho1(i),ntt,iye1(i))) .gt. &
             fun(i) ) then
            
            write (*,*) i,eos_cast(highden_table(eos_mode)%LED(irho1(i),ntt,iye1(i))), &
                    fun(i)
            raise_abort("D1")
         endif
      enddo
   endif
#endif /* DYNAMIC_EOS && DEBUG */

      select case (mode)
        case(1)   ! Density, ye and  t e m p e r a t u r e  given

          call xtrct_intrpl_data_eos(nos,eos_mode,intrpl_data,irho0,irho1, &
                                     item0,item1,iye0,iye1)

#ifdef IBM_OPT
! on IBM Power6 it is better to split the loops, however on
! Vector-Machines you do not want to do this

!IBM* PREFETCH_FOR_LOAD (highden_table(eos_mode)%LPR  ,wr0, wr1, wt0, wt1, wy0, wy1, irho0, irho1, item0, item1, iye0, iye1, highden_table(eos_mode)%LED, highden_table(eos_mode)%CEI, highden_table(eos_mode)%CNI,highden_table(eos_mode)%CPI,highden_table(eos_mode)%CUI,highden_table(eos_mode)%GAI,highden_table(eos_mode)%STI)

#endif /* IBM_OPT */
          do i = 1, nos
#define      index DRU
#define      result pr
#include     "high_density_eos.trilinear.X90"
             !debug(i)
             !debug(wr0(i))
             !debug(wr1(i))
             !debug(wt0(i))
             !debug(wt1(i))
             !debug(wy0(i))
             !debug(wy1(i))
             pr(i) = 10.0_rk**pr(i)

#define      index ENER
#define      result ed
#include     "high_density_eos.trilinear.X90"
             ed(i) = 10.0_rk**ed(i)
             
#define      index MUE
#define      result ce
#include     "high_density_eos.trilinear.X90"

#define      index MUN
#define      result cn
#include     "high_density_eos.trilinear.X90"

#define      index MUP
#define      result cp
#include     "high_density_eos.trilinear.X90"

#define      index MUNE
#define      result cu
#include     "high_density_eos.trilinear.X90"

#define      index ADIA
#define      result ga
#include     "high_density_eos.trilinear.X90"

    if(high_den_eos .eq. 2) then
  ! test if we're out of WOLFF boundaries


             if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
                if ((eos_cast(highden_table(eos_mode)%GAI(irho0(i),item0(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho0(i),item0(i),iye1(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho0(i),item1(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho0(i),item1(i),iye1(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item0(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item0(i),iye1(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item1(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item1(i),iye1(i))) .lt. 0.)) then

                   raise_abort("out of wolff boundaries2")
                end if
             endif
    endif ! high_den_eos == 2

#define      index ENTR
#define      result st
#include     "high_density_eos.trilinear.X90"

#ifdef IBM_OPT
         enddo

!IBM* PREFETCH_FOR_LOAD (highden_table(eos_mode)%XXN  ,wr0, wr1, wt0, wt1, wy0, wy1, irho0, irho1, item0, item1, iye0, iye1, highden_table(eos_mode)%XXP, highden_table(eos_mode)%XXA, highden_table(eos_mode)%CNI,highden_table(eos_mode)%XXH,highden_table(eos_mode)%XHZ,highden_table(eos_mode)%XHA)
          do i = 1, nos
#endif /* IBM_OPT */

#define      index MASSN
#define      result xn
#include     "high_density_eos.trilinear.X90"

    if(high_den_eos .eq. 2) then
             if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
                if ((eos_cast(highden_table(eos_mode)%XXN(irho0(i),item0(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%XXN(irho0(i),item0(i),iye1(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%XXN(irho0(i),item1(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%XXN(irho0(i),item1(i),iye1(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item0(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item0(i),iye1(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item1(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
                   raise_abort("eos_ls(): xn negative; pos 1")
                end if
             endif
    endif ! high_den_eos == 2

#define      index MASSP
#define      result xp
#include     "high_density_eos.trilinear.X90"

    if(high_den_eos .eq. 2) then
           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if ((eos_cast(highden_table(eos_mode)%XXP(irho0(i),item0(i),iye0(i))) .lt. 0.) .or.        &
                   (eos_cast(highden_table(eos_mode)%XXP(irho0(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho0(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho0(i),item1(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
!                 raise_abort("eos_ls(): xp negative; pos 1")
         
              end if
           endif
    endif ! high_den_eos == 2

#define      index MASSA
#define      result xa
#include     "high_density_eos.trilinear.X90"
 
    if(high_den_eos .eq. 2) then

           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if ((eos_cast(highden_table(eos_mode)%XXA(irho0(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho0(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho0(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho0(i),item1(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
                 raise_abort("eos_ls(): xa negative; pos 1")
              end if
           endif
    endif ! high_den_eos == 2

#define      index MASSH
#define      result xh
#include     "high_density_eos.trilinear.X90"

    if(high_den_eos .eq. 2) then


           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if ((eos_cast(highden_table(eos_mode)%XXH(irho0(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho0(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho0(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho0(i),item1(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
                 raise_abort("eos_ls(): xh negative; pos 1")
              endif

           end if
    endif ! high_den-eos == 2

#define      index CHRG
#define      result z_a
#include     "high_density_eos.trilinear.X90"

#define      index ATMAS
#define      result a
#include     "high_density_eos.trilinear.X90"

          enddo

       case(2:3)  ! search for temperature
          fun(:) = ed(:)

          call fastlog10(fun, size(fun))

         do i = 1, nos

           item0(i) = 1
#ifndef DYNAMIC_EOS
           item1(i) = highden_table(eos_mode)%ntt  ! this can either be to global index or local one. CAREFUL

#else
           item1(i) = highden_table(eos_mode)%ntt_lc  ! this can either be to global index or local one. CAREFUL
#endif
   
           scl(i) =     &
               wr0(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%LED(irho0(i),item0(i),iye0(i)))  &
                         + wy1(i)*eos_cast(highden_table(eos_mode)%LED(irho0(i),item0(i),iye1(i))) )        &
             + wr1(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%LED(irho1(i),item0(i),iye0(i)))  &
                         + wy1(i)*eos_cast(highden_table(eos_mode)%LED(irho1(i),item0(i),iye1(i))) )

           sch(i) =     &
               wr0(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%LED(irho0(i),item1(i),iye0(i)))  &
                         + wy1(i)*eos_cast(highden_table(eos_mode)%LED(irho0(i),item1(i),iye1(i))) )        &
             + wr1(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%LED(irho1(i),item1(i),iye0(i)))  &
                         + wy1(i)*eos_cast(highden_table(eos_mode)%LED(irho1(i),item1(i),iye1(i))) )

         enddo
   
         if (ANY(fun(:).lt.scl(:)) ) then
            ler=.TRUE. 
#ifdef DYNAMIC_EOS
            retry_reload_counter = retry_reload_counter + 1
            out_subcube(:) = .false.
            out_subcube(3) = .true.
            if (retry_reload_counter .gt. 50) then
               raise_abort("too many reloads of EOS")  
            endif
            goto 1111
#endif /* DYNAMIC_EOS */
         endif


         if (ANY(fun(:).gt.sch(:)) ) then
            ler=.TRUE.
#ifdef DYNAMIC_EOS
            retry_reload_counter = retry_reload_counter + 1
            out_subcube(:) = .false.
            out_subcube(4) = .true.
            if (retry_reload_counter .gt. 50) then
               raise_abort("too many reloads of EOS")  
            endif
            goto 1111
#endif /* DYNAMIC_EOS */
         endif



         if (ANY(fun(:).lt.scl(:).or.fun(:).gt.sch(:))) then 
            ler=.TRUE.
         else
            ler=.false.
         endif


    !-- report an error if boundaries exceeded
        if(ler) then

           do i = 1, nos
              if (fun(i).lt.scl(i).or.fun(i).gt.sch(i)) then
                 call show_error_screen("eos_ls","energy inversion", &
                      i,scl(i),fun(i),sch(i))
              endif   
           enddo

         call show_error_screen("eos_ls","EOS boundaries exceeded:")

         call show_error_screen("eos_ls","\[k,rho,ed,ye\]")
         call show_error_screen(nos,rho,ed,xi(:,nuc))
         return
      endif

  ! --- initalis. for bisection
        do i = 1, nos
           if ( item1(i) .eq. item0(i)+1 ) then
              isect(i) = 0._rk
           else
              isect(i) = 1._rk
           endif
        enddo

  ! --- bisection iteration
    100 continue

         do i = 1, nos
           if ( isect(i) .eq. 1._rk ) then

              ittm(i) = (item1(i) + item0(i)) / 2

              scm =     &
                  wr0(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%LED(irho0(i),ittm(i),iye0(i)))        &
                            + wy1(i)*eos_cast(highden_table(eos_mode)%LED(irho0(i),ittm(i),iye1(i))) )      &
                + wr1(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%LED(irho1(i),ittm(i),iye0(i)))        &
                            + wy1(i)*eos_cast(highden_table(eos_mode)%LED(irho1(i),ittm(i),iye1(i))) )

              if ( fun(i) .gt. scm ) then
                 scl(i)  = scm
                 item0(i) = ittm(i)
              else
                 sch(i)  = scm
                 item1(i) = ittm(i)
              endif
           endif
         enddo

        do i = 1, nos
           if ( item1(i) .gt. item0(i)+1 ) then
              isect(i) = 1._rk
           else
              isect(i) = 0._rk
           endif
        enddo

        sumisec = 0._rk
        do i = 1, nos
           sumisec = sumisec + isect(i)
        enddo

        if ( sumisec .ne. 0._rk ) goto 100
  ! -- end bisection


        do i = 1, nos

           wt1(i) = ( fun(i) - scl(i) ) / ( sch(i) - scl(i) )
           wt0(i) = 1.0_rk - wt1(i)

           tem(i) = 10.0_rk**( wt0(i)*eos_cast(highden_table(eos_mode)%LTT_HD(item0(i))) +        &
                wt1(i)*eos_cast(highden_table(eos_mode)%LTT_HD(item1(i))) )


        enddo
        
          call xtrct_intrpl_data_eos(nos,eos_mode,intrpl_data,irho0,irho1, &
                                     item0,item1,iye0,iye1)
      
          do i = 1, nos

  ! ... pressure
#define    index DRU
#define    result pr
#include   "high_density_eos.trilinear.X90"
           pr(i) = 10.0_rk**pr(i)

  ! ... internal energy density
#define    index ENER
#define    result ed
#include   "high_density_eos.trilinear.X90"
           ed(i) = 10.0_rk**ed(i)

  ! chemical potentials
#define    index MUE
#define    result ce
#include   "high_density_eos.trilinear.X90"

#define    index MUN
#define    result cn
#include   "high_density_eos.trilinear.X90"

#define    index MUP
#define    result cp
#include   "high_density_eos.trilinear.X90"

#define    index MUNE
#define    result cu
#include   "high_density_eos.trilinear.X90"

#ifdef IBM_OPT
        enddo

        do i=1, nos
#endif /* IBM_OPT */

#define    index ADIA
#define    result ga
#include   "high_density_eos.trilinear.X90"

    if(high_den_eos .eq. 2) then
  ! test if we're out of WOLFF boundaries

             if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
                if ((eos_cast(highden_table(eos_mode)%GAI(irho0(i),item0(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho0(i),item0(i),iye1(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho0(i),item1(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho0(i),item1(i),iye1(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item0(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item0(i),iye1(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item1(i),iye0(i))) .lt. 0.) .or.       &
                     (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
                   
                   raise_abort("eos_ls(): Out of Wolff boundaries 2")
                end if
             endif
    endif ! high_den_eos == 2

  ! ... entropy
#define    index ENTR
#define    result st
#include   "high_density_eos.trilinear.X90"

  ! .. composition
#define    index MASSN
#define    result xn
#include   "high_density_eos.trilinear.X90"

    if(high_den_eos .eq. 2) then
           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if((eos_cast(highden_table(eos_mode)%XXN(irho0(i),item0(i),iye0(i))) .lt. 0.) .or.        &
                   (eos_cast(highden_table(eos_mode)%XXN(irho0(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho0(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho0(i),item1(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
                 raise_abort("eos_ls(): xn negative; pos 1")
              end if
           endif
    endif ! high_den_eos

#define    index MASSP
#define    result xp
#include   "high_density_eos.trilinear.X90"

    if(high_den_eos .eq. 2) then
           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if ((eos_cast(highden_table(eos_mode)%XXP(irho0(i),item0(i),iye0(i))) .lt. 0.) .or.        &
                   (eos_cast(highden_table(eos_mode)%XXP(irho0(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho0(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho0(i),item1(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
                 raise_abort("eos_ls(): xp negative; pos 1")
              end if
           endif
    endif ! high_den_eos == 2
#define    index MASSA
#define    result xa
#include   "high_density_eos.trilinear.X90"

    if(high_den_eos .eq. 2) then

           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if ((eos_cast(highden_table(eos_mode)%XXA(irho0(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho0(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho0(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho0(i),item1(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
                 raise_abort("eos_ls(): xa negative; pos 1")
              end if
           endif
    endif ! high_den_eos == 2

#define    index MASSH
#define    result xh
#include   "high_density_eos.trilinear.X90"


    if(high_den_eos .eq. 2) then


           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if ((eos_cast(highden_table(eos_mode)%XXH(irho0(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho0(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho0(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho0(i),item1(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
                 raise_abort("eos_ls(): xh negative; pos 1")
              endif
           end if
    endif ! high_den-eos == 2

#define    index CHRG
#define    result z_a
#include   "high_density_eos.trilinear.X90"

#define    index ATMAS
#define    result a
#include   "high_density_eos.trilinear.X90"

        enddo

        case(4)
           fun(:) = st(:)


        do i = 1, nos
           item0(i) = 1
           item1(i) = NTT ! CAREFUL this can either be the global or local index
   
           scl(i) =     &
               wr0(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%STI(irho0(i),item0(i),iye0(i)))  &
                         + wy1(i)*eos_cast(highden_table(eos_mode)%STI(irho0(i),item0(i),iye1(i))) )&
             + wr1(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%STI(irho1(i),item0(i),iye0(i)))  &
                         + wy1(i)*eos_cast(highden_table(eos_mode)%STI(irho1(i),item0(i),iye1(i))) )


           sch(i) =     &
               wr0(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%STI(irho0(i),item1(i),iye0(i)))  &
                         + wy1(i)*eos_cast(highden_table(eos_mode)%STI(irho0(i),item1(i),iye1(i))) )&
             + wr1(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%STI(irho1(i),item1(i),iye0(i)))  &
                         + wy1(i)*eos_cast(highden_table(eos_mode)%STI(irho1(i),item1(i),iye1(i))) )


        enddo
        if (ANY(fun(:).lt.scl(:).or.fun(:).gt.sch(:))) ler=.TRUE.
        

  !-- report an error if boundaries exceeded
        if(ler) then
           do i = 1, nos
              if (fun(i).lt.scl(i).or.fun(i).gt.sch(i)) then
                call show_error_screen("eos_ls","entropy inversion", &
                     i,scl(i),fun(i),sch(i))
              endif   
           enddo

         call show_error_screen("eos_ls","EOS boundaries exceeded:")

          call show_error_screen("eos_ls","\[k,rho,sto,ye\]")         
          call show_error_screen(nos,rho,st,xi(:,nuc))

           return
        endif 

  ! --- initalis. for bisection
        do i = 1, nos
           if ( item1(i) .eq. item0(i)+1 ) then
              isect(i) = 0._rk
           else
              isect(i) = 1._rk
           endif
        enddo

  ! --- bisection iteration
   101  continue

        do i = 1, nos
           if ( isect(i) .eq. 1._rk ) then

              ittm(i) = (item1(i) + item0(i)) / 2

              scm =     &
                  wr0(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%STI(irho0(i),ittm(i),iye0(i)))        &
                            + wy1(i)*eos_cast(highden_table(eos_mode)%STI(irho0(i),ittm(i),iye1(i))) )      &
                + wr1(i) * (  wy0(i)*eos_cast(highden_table(eos_mode)%STI(irho1(i),ittm(i),iye0(i)))        &
                            + wy1(i)*eos_cast(highden_table(eos_mode)%STI(irho1(i),ittm(i),iye1(i))) )

              if ( fun(i) .gt. scm ) then
                 scl(i)  = scm
                 item0(i) = ittm(i)
              else
                 sch(i)  = scm
                 item1(i) = ittm(i)
              endif
           endif
        enddo

        do i = 1, nos
           if ( item1(i) .gt. item0(i)+1 ) then
              isect(i) = 1._rk
           else
              isect(i) = 0._rk
           endif
        enddo

        sumisec = 0._rk
        do i = 1, nos
           sumisec = sumisec + isect(i)
        enddo


        if ( sumisec .ne. 0._rk ) goto 101
  ! -- end bisection

        do i = 1, nos

  !         write(44,*) item0(i),10**ltt_hd(item0(i)),fun(i),
  !     &        (10**sch(i)-10**scl(i))/10**fun(i)

           wt1(i) = ( fun(i) - scl(i) ) /       &
                ( sch(i) - scl(i) )
           wt0(i) = 1.0_rk- wt1(i)

           tem(i) = 10.0_rk**( wt0(i)*eos_cast(highden_table(eos_mode)%ltt_hd(item0(i))) +        &
                wt1(i)*eos_cast(highden_table(eos_mode)%ltt_hd(item1(i))) )


        enddo

          call xtrct_intrpl_data_eos(nos,eos_mode,intrpl_data,irho0,irho1, &
                                     item0,item1,iye0,iye1)


        do i = 1, nos

  ! ... pressure
#define    index DRU
#define    result pr
#include   "high_density_eos.trilinear.X90"
           pr(i) = 10.0_rk**pr(i)
           
  ! ... internal energy density
#define    index ENER
#define    result ed
#include   "high_density_eos.trilinear.X90"
           ed(i) = 10.0_rk**ed(i)
 
  ! chemical potentials
#define    index MUE
#define    result ce
#include   "high_density_eos.trilinear.X90"

#define    index MUN
#define    result cn
#include   "high_density_eos.trilinear.X90"

#define    index MUP
#define    result cp
#include   "high_density_eos.trilinear.X90"

#define    index MUNE
#define    result cu
#include   "high_density_eos.trilinear.X90"

  ! .. gamma
#define    index ADIA
#define    result ga
#include   "high_density_eos.trilinear.X90"

    if (high_den_eos .eq. 2) then
  ! test if we're out of WOLFF boundaries

           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if ((eos_cast(highden_table(eos_mode)%GAI(irho0(i),item0(i),iye0(i))) .lt. 0.) .or.       &
                   (eos_cast(highden_table(eos_mode)%GAI(irho0(i),item0(i),iye1(i))) .lt. 0.) .or.       &
                   (eos_cast(highden_table(eos_mode)%GAI(irho0(i),item1(i),iye0(i))) .lt. 0.) .or.       &
                   (eos_cast(highden_table(eos_mode)%GAI(irho0(i),item1(i),iye1(i))) .lt. 0.) .or.       &
                   (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item0(i),iye0(i))) .lt. 0.) .or.       &
                   (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item0(i),iye1(i))) .lt. 0.) .or.       &
                   (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item1(i),iye0(i))) .lt. 0.) .or.       &
                   (eos_cast(highden_table(eos_mode)%GAI(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
                 
                 raise_abort("eos_ls(): Out of Wolff boundaries 2")
              end if
           endif
    endif ! high_den_eos == 2

  ! ... entropy
#define    index ENTR
#define    result st
#include   "high_density_eos.trilinear.X90"

  ! .. composition
#define    index MASSN
#define    result xn
#include   "high_density_eos.trilinear.X90"

      if (high_den_eos .eq. 2) then
           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if ((eos_cast(highden_table(eos_mode)%XXN(irho0(i),item0(i),iye0(i))) .lt. 0.) .or.        &
                   (eos_cast(highden_table(eos_mode)%XXN(irho0(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho0(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho0(i),item1(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXN(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
                 
                 raise_abort("eos_ls(): xn negative; pos 1")
              end if
           endif
      endif ! high_den_eos == 2
#define    index MASSP
#define    result xp
#include   "high_density_eos.trilinear.X90"

      if (high_den_eos .eq. 2) then
           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if ((eos_cast(highden_table(eos_mode)%XXP(irho0(i),item0(i),iye0(i))) .lt. 0.) .or.        &
                   (eos_cast(highden_table(eos_mode)%XXP(irho0(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho0(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho0(i),item1(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXP(irho1(i),item1(i),iye1(i))) .lt. 0.)) then

                 raise_abort("eos_ls(): xp negative; pos 1")
              end if
           endif

      endif ! high_den_eos == 2
#define    index MASSA
#define    result xa
#include   "high_density_eos.trilinear.X90"

      if (high_den_eos .eq. 2) then
           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if ((eos_cast(highden_table(eos_mode)%XXA(irho0(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho0(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho0(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho0(i),item1(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXA(irho1(i),item1(i),iye1(i))) .lt. 0.)) then

                 raise_abort("eos_ls(): xa negative; pos 1")
          end if
       endif

      endif ! high_den_eos == 2
#define    index MASSH
#define    result xh
#include   "high_density_eos.trilinear.X90"

      if (high_den_eos .eq. 2) then
           if (eos_mode .eq. 4 .or. highden_table(eos_mode)%name .eq. "wolffeos_low") then 
              if ((eos_cast(highden_table(eos_mode)%XXH(irho0(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho0(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho0(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho0(i),item1(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item0(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item0(i),iye1(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item1(i),iye0(i))) .lt. 0.) .or. &
                   (eos_cast(highden_table(eos_mode)%XXH(irho1(i),item1(i),iye1(i))) .lt. 0.)) then
       
                 
                 raise_abort("eos_ls(): xh negative; pos 1")
              end if
           endif
      endif ! high_den_eos == 2
#define    index CHRG
#define    result z_a
#include   "high_density_eos.trilinear.X90"

#define    index ATMAS
#define    result a
#include   "high_density_eos.trilinear.X90"

        enddo

        case default
           raise_abort('something wrong in eos_ls')
      end select

      ! counter-variable
#ifdef DYNAMIC_EOS
      call omp_set_nest_lock(exit_lock_counter)
      thread_counter = thread_counter - 1
      call omp_unset_nest_lock(exit_lock_counter)

#endif /* DYNAMIC_EOS */

! this enforces X_n+X_p+X_alfa+X_heavy=1.0 to machine accuracy
!  which might be greater than the accuracy of the table-entries
!   by correcting the largest of the four summands
!   
      do i = 1, nos
         scratch(1)=max(0._rk,xn(i))
         scratch(2)=max(0._rk,xp(i))
         scratch(3)=max(0._rk,xa(i))
         scratch(4)=max(0._rk,xh(i))

         kk=SUM(MAXLOC(scratch))
         summ=SUM(scratch)-scratch(kk)
         scratch(kk)=1.0_rk-summ

         xn(i)=scratch(1)
         xp(i)=scratch(2)
         xa(i)=scratch(3)
         xh(i)=scratch(4)
      enddo

      select case(nsemode)
        case(1) ! use the NSE-composition calculated by LS-EoS
          xi(:,n_n) = xn
          xi(:,n_p) = xp
          xi(:,n_he4) = xa
          xi(:,n_he4+1:nuc-1) = 0.0
          xhrep = xh
          if (config%qn .gt. 21) then
             xi(:,n_d:n_he3)=0.0_rk
          endif
          za(:,1) = z_a * a
          za(:,2) = a
        case(0) ! use the NSE-composition MAPPED to the ensemble of nuclei being advected
!CDIR IEXPAND (mapnuc)
!DIR$ PREFERVECTOR
          do i = 1, nos
            z = z_a(i) * a(i)
            call mapnuc(xi(i,1:nuc),rho(i),xn(i),xp(i),xa(i),xh(i),z,a(i))

!           write(*,'(I4,6E13.5,1e21.13)') k,xi(k,1),xi(k,2), &
!           xi(k,3),xi(k,nuc-3),xi(k,nuc-2),xi(k,nuc-1),      &
!           SUM(xi(k,1:nuc-1))
          enddo
        case default
          raise_abort("lseos(): nsemode not implemented")
      end select


      if (highden_table(eos_mode)%eos_energy_offs) then

         ! High density EoS is tabulated with a different
         ! energy offset than the "normal" 8.8 MeV
         ! we have to correct for this here

         ed(:) = ed(:) -                                   &
                 (highden_table(eos_mode)%enullm - enullm) &
                 * rho(:) * pc_meverg / pc_mb

      endif

      
      end subroutine highdensity_eos

end module high_density_eos

