! unused defines:
#undef BURN_NETW_NOA


!> \verbatim
!>
!> This module provides the nuclear parameters (what nuclei are used ?
!> their masses, their index in the xnuc-array...)
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
module nucparam

  use precision

  use phycon

  use configure

  implicit none
! LOCAL variables that are not in modules

  public
  private nc

  SAVE
      
! When using the Thielemann network, we can only switch to
! the NSE table, once NSE is FULLY established
!
! When "flashing" is used to describe nuclear transmutations
! the NSE treshold temperature must be lower (T \approx 0.5 MeV):
! While full NSE may not be reached at this point, QSE already
! applies and is well approximated by the NSE composition.
! Retaining everything as 56Ni up to 8*10^9 K and neglecting
! dissociation into alphas can have fatal consequences.

real(kind=rk) :: moffs
  real(kind=rk), allocatable :: pc_nuc(:,:), mbar(:)
!
  integer(kind=ik) :: nucset

  real(kind=rk), allocatable :: ebind(:) 


#ifdef PROGRAM_remap
  integer(kind=ik) :: restmass_version = 0
#else
!  integer(kind=ik), parameter :: restmass_version = RESTMASS_VERS
#endif

  integer(kind=ik)       ::  n_rep

  character(len=5), allocatable :: name_xnuc(:)

! add other symbols here      
  integer(kind=ik)        :: n_n, n_p, n_d, n_t, n_he3, n_he4, n_c12,    &
                             n_n14, n_o16, n_ne20, n_mg24, n_si28,       &
                             n_s32, n_ar36, n_ca40, n_ti44, n_cr48,      &
                             n_ni56, n_mn54, n_fe52, n_fe56, n_fe60,     &
                             n_ni70, n_ni120, n_zr200 


  integer(kind=ik)         :: nc
  real(kind=rk), parameter :: enullm=8.8_rk, massn=939.5731_rk

!.. pc_nuc(:,1):    Z
!.. pc_nuc(:,2):    A
!.. pc_nuc(:,3):    A mass excess [in MeV]
!.. pc_nuc(:,4:7):  fit parameters for partition function (stat. weight)
! the statistical fit function is: g=a*(1+exp(b/t9+c+d*t9)



contains
!> \verbatim
!>
!> This subroutine identifies the used nuclear species and assigns
!> an unique index to each nucleus
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
  subroutine get_nuc_indices(mem_global)

    use precision
    use abort

! purpose: finds the index values n_n, n_p ,... for a selection
!          of nuclear species (defined by (A,Z) at runtime
!          this is to avoid that changes in the above table for 
!          pc_nuc affect the identification of a particular species
!          somewhere elese in the code (e.g. subroutine burn)

    use print_stdout_mod
    use mpi_vertex

    implicit none
! LOCAL variables that are not in modules
    real(kind=rk), intent(inout)  :: mem_global
    character(len=60)             :: data_set
    integer(kind=ik)              :: ndat
    integer(kind=ik)              :: ndat_rd, i
    real(kind=rk)                 :: rd1, rd2, rd3, rd4, rd5, rd6, rd7, &
                                     mem_local
    integer(kind=ik)              :: istat

    mem_local = 0._rk

    allocate(mbar(config%qn), pc_nuc(config%qn-1,7), name_xnuc(config%qn), &
             stat=istat)

    if (istat .ne. 0) then
        raise_abort("get_nuc_indicies(): allocation of failed")
    end if

    mem_local = config%qn*2*8._rk + (config%qn-1)*7*8._rk

    
    allocate(ebind(config%qn), stat=istat)

    if (istat .ne. 0) then
        raise_abort("get_nuc_indicies(): allocation of failed 2")
    end if

    mem_local = mem_local + config%qn*8._rk


    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "totare_hy")

    select case (config%qn)
      case (21)
       if (nucset .eq. 0) then
          data_set="./tables/20_species_std.data"
       else
          data_set="./tables/20_species_thiel.data"
       endif

      case (25)
       data_set="./tables/24_species_thiel.data"

      case (26)
       data_set="./tables/25_species_thiel.data"

      case default
        raise_abort("Invalid value for parameter 'qn', valid are 21, 25, and 26")

    end select

    ndat=config%qn - 1

    ! Load nuclear data sets

    call printit_taskX(0,"Using nuclear data set from ",trim(data_set))

    open(67, file=data_set, form='unformatted', status='old')
    read(67) ndat_rd
    if (ndat_rd .ne. ndat) then
       raise_abort("mod_nucparam.F90(): wrong dimension of data_set "//data_set)
    endif

    do i=1,ndat
       read(67) rd1, rd2, rd3, rd4, rd5, rd6, rd7
       pc_nuc(i,1)=rd1
       pc_nuc(i,2)=rd2
       pc_nuc(i,3)=rd3
       pc_nuc(i,4)=rd4
       pc_nuc(i,5)=rd5
       pc_nuc(i,6)=rd6
       pc_nuc(i,7)=rd7
    enddo

    close(67)


    n_n =-99999
    n_p =-99999
    n_d =-99999
    n_t =-99999
    n_he3=-99999
    n_he4=-99999
    n_c12 =-99999
    n_n14=-99999
    n_o16 =-99999
    n_ne20=-99999
    n_mg24=-99999
    n_si28=-99999
    n_s32 =-99999
    n_ar36=-99999
    n_ca40=-99999
    n_ti44=-99999
    n_cr48=-99999
    n_ni56=-99999
    n_mn54=-99999
    n_fe52=-99999
    n_fe56=-99999
    n_fe60=-99999
    n_ni70=-99999
    n_ni120=-99999
    n_zr200=-99999

    call printit_taskX(0," ")
    call printit_taskX(0,"get_nuc_indices:")

    name_xnuc(config%qn) = "   Ye"


    do nc=1, config%qn-1

       if (pc_nuc(nc,1) .eq. 0._rk .and. pc_nuc(nc,2) .eq. 1._rk) then
          n_n=nc

          call printit_taskX(0,"  n: ",n_n, pc_nuc(nc,3)*pc_meverg/pc_mb)
          name_xnuc(n_n) = "    n"
       endif
       if (pc_nuc(nc,1) .eq. 1._rk .and. pc_nuc(nc,2) .eq. 1._rk) then
          n_p=nc

          call printit_taskX(0,"  p: ",n_p, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_p) = "    p"
       endif
       if (pc_nuc(nc,1) .eq. 1._rk .and. pc_nuc(nc,2) .eq. 2._rk) then
          n_d=nc
          call printit_taskX(0," 2H: ",n_d, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_d) = "   H2"
       endif
       if (pc_nuc(nc,1) .eq. 1._rk .and. pc_nuc(nc,2) .eq. 3._rk) then
          n_t=nc
          call printit_taskX(0,"  3H: ",n_t, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_t) = "   H3"
       endif
       if (pc_nuc(nc,1) .eq. 2._rk .and. pc_nuc(nc,2) .eq. 3._rk) then
          n_he3=nc
          call printit_taskX(0," 3He: ",n_he3, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_he3) = "  He3"
       endif
       if (pc_nuc(nc,1) .eq. 2._rk .and. pc_nuc(nc,2) .eq. 4._rk) then
          n_he4=nc
          call printit_taskX(0," 4He: ",n_he4, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_he4) = "  He4"
       endif
       if (pc_nuc(nc,1) .eq. 6._rk .and. pc_nuc(nc,2) .eq.12._rk) then
          n_c12=nc
          call printit_taskX(0,"12C : ",n_c12, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_c12) = "  C12"
       endif
       if (pc_nuc(nc,1) .eq. 7._rk .and. pc_nuc(nc,2) .eq. 14._rk) then
          n_n14=nc
          call printit_taskX(0,"14N : ",n_n14, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_n14) = "  N14"
       endif
       if (pc_nuc(nc,1) .eq. 8._rk .and. pc_nuc(nc,2) .eq. 16._rk) then
          n_o16=nc
          call printit_taskX(0,"16O : ",n_o16, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_o16) = "  O16"
       endif
       if (pc_nuc(nc,1) .eq. 10._rk .and. pc_nuc(nc,2) .eq. 20._rk) then
          n_ne20=nc
         call printit_taskX(0,"20Ne: ",n_ne20, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_ne20) = " Ne20"
       endif
       if (pc_nuc(nc,1) .eq. 12._rk .and. pc_nuc(nc,2) .eq. 24._rk) then
          n_mg24=nc
          call printit_taskX(0,"24Mg: ",n_mg24, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_mg24) = " Mg24"
       endif
       if (pc_nuc(nc,1) .eq. 14._rk .and. pc_nuc(nc,2) .eq. 28._rk) then
          n_si28=nc
          call printit_taskX(0,"28Si: ",n_si28, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_si28) = " Si28"
       endif
       if (pc_nuc(nc,1) .eq. 16._rk .and. pc_nuc(nc,2) .eq. 32._rk) then
          n_s32=nc
          call printit_taskX(0,"32S : ",n_s32, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_s32) = "  S32"
       endif
       if (pc_nuc(nc,1) .eq. 18._rk .and. pc_nuc(nc,2) .eq. 36._rk) then
          n_ar36=nc
          call printit_taskX(0,"36Ar: ",n_ar36, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_ar36) = " Ar36"
       endif
       if (pc_nuc(nc,1) .eq. 20._rk .and. pc_nuc(nc,2) .eq. 40._rk) then
          n_ca40=nc
          call printit_taskX(0,"40Ca: ",n_ca40, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_ca40) = " Ca40"
       endif
       if (pc_nuc(nc,1) .eq. 22._rk .and. pc_nuc(nc,2) .eq. 44._rk) then
          n_ti44=nc
          call printit_taskX(0,"44Ti: ",n_ti44, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_ti44) = " Ti44"
       endif
       if (pc_nuc(nc,1) .eq.  24._rk .and. pc_nuc(nc,2) .eq.  48._rk) then
          n_cr48=nc
          call printit_taskX(0,"48Cr: ",n_cr48, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_cr48) = " Cr48"
       endif
       if (pc_nuc(nc,1) .eq.  26._rk .and. pc_nuc(nc,2) .eq. 52._rk) then
          n_fe52=nc
          call printit_taskX(0,"52Fe: ",n_fe52, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_fe52) = " Fe52"
       endif
       if (pc_nuc(nc,1) .eq. 26._rk .and. pc_nuc(nc,2) .eq. 56._rk) then
          n_fe56=nc
          call printit_taskX(0,"56Fe: ",n_fe56, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_fe56) = " Fe56"
       endif

       if (pc_nuc(nc,1) .eq. 28._rk .and. pc_nuc(nc,2) .eq. 56._rk) then
          n_ni56=nc
          call printit_taskX(0,"56Ni: ",n_ni56, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_ni56) = " Ni56"
       endif
       if (pc_nuc(nc,1) .eq. 25._rk .and. pc_nuc(nc,2) .eq. 54._rk) then
          n_mn54=nc
          call printit_taskX(0,"54Mn: ",n_mn54, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_mn54) = " Mn54"
       endif
       if (pc_nuc(nc,1) .eq. 26._rk .and. pc_nuc(nc,2) .eq. 60._rk) then
          n_fe60=nc
          call printit_taskX(0,"60Fe: ",n_fe60, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_fe60) = " Fe60"
       endif
       if (pc_nuc(nc,1) .eq. 28._rk .and. pc_nuc(nc,2) .eq. 70._rk) then
          n_ni70=nc
          call printit_taskX(0,"70Ni: ",n_ni70, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_ni70) = " Ni70"
       endif
       if (pc_nuc(nc,1) .eq. 28._rk .and. pc_nuc(nc,2) .eq. 120._rk) then
          n_ni120=nc
          call printit_taskX(0,"120Ni: ",n_ni120, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_ni120) = "Ni120"
       endif
       if (pc_nuc(nc,1) .eq. 40._rk .and. pc_nuc(nc,2) .eq. 200._rk) then
          n_zr200=nc
          call printit_taskX(0,"200Zr: ",n_zr200, pc_nuc(nc,3)*pc_meverg/pc_mb)

          name_xnuc(n_zr200) = "Zr200"
       endif
    enddo
    call printit_taskX(0,"=====================================================================")

#if !(defined(PROGRAM_remap))

    if (config%use_flash_c) then
       if (n_si28.le.0 .or. n_mg24.le.0 .or. n_o16.le.0 &
                       .or. n_ne20.le.0) then
          raise_abort("get_nuc_indices: Carbon burning not possible")
       else
          call printit_taskX(0,"get_nuc_indices: Carbon burning switched on")
       endif
    endif

    if (config%use_flash_o) then
       if (n_si28.le.0 .or. n_mg24.le.0 .or. n_o16.le.0 .or. n_ne20.le.0) then
          raise_abort("get_nuc_indices: Oxygen burning not possible")
       else
          call printit_taskX(0,"get_nuc_indices: Oxygen burning switched on")
       endif
    endif

    if (config%use_flash_si) then
       if (n_si28.le.0 .or. n_ni56.le.0) then
          raise_abort("get_nuc_indices: Silicon burning not possible")
       else
          call printit_taskX(0,"get_nuc_indices: Silicon burning switched on")
       endif
    endif

#ifdef BURN_SS
    if (n_fe56.le.0 .or. n_si28.le.0 .or. n_ni56.le.0 .or. n_c12.le.0 &
         .or. n_mg24.le.0 .or. n_o16.le.0 .or. n_ne20.le.0) then
       raise_abort("get_nuc_indices(): Supersonic burning not possible")
    else
       call printit_taskX(0,"get_nuc_indices: Supersonic burning switched on")
    endif
#endif

#ifdef BURN_NETW_NOA
    if (n_he4.le.0 .or. n_c12.le.0 .or.n_ne20.le.0 .or.  n_o16.le.0 .or. n_mg24.le.0) then
       raise_abort("get_nuc_indices(): Burning Network No Alphas not possible")
    else
       call printit_taskX(0,"get_nuc_indices: Burning Network No Alphas switched on")
    endif
#endif

#ifdef TAK_RATES
    call printit_taskX(0,"Takahara electron capture rates switched on")
#endif

#ifdef BURN_NETWORK
    if (n_he4.le.0 .or. n_c12.le.0 .or.n_ne20.le.0 .or.  n_o16.le.0 .or. n_mg24.le.0) then
       raise_abort("get_nuc_indices(): Burn network not possible")
    else
       call printit_taskX(0,"get_nuc_indices: Burn network switched on")
    endif
#endif

#endif /* not remaper */

    ebind(:)=0.0_rk
    do nc=1,config%qn-1
       ebind(nc)=pc_nuc(nc,3)*pc_meverg/pc_mb
    enddo



    if (config%restmass_version .gt. 0) then
#ifndef PROGRAM_remap
       call printit_taskX(0,"Running with No-Rest-Mass version ",config%restmass_version)
#endif
       if (config%restmass_version .eq. 1) n_rep=n_fe56

       do nc=1,config%qn-1
          mbar(nc) = (pc_nuc(nc,3)+pc_nuc(nc,2)*wc_mb) &
            /(pc_nuc(nc,2)*pc_mb)*pc_meverg
       enddo


       if (config%restmass_version .eq. 3) then
          !c electron rest mass is subtracted
          mbar(config%qn) = wc_me/pc_mb*pc_meverg
          raise_abort("restmass_version 3 is not yet energy conservative!")
       else
          mbar(config%qn) = 0.0_rk
       endif

 endif ! restmass_version > 0

    moffs = (-massn + enullm)/pc_mb*pc_meverg
  end subroutine get_nuc_indices

end module nucparam
