module multigrid_rt

  use precision


  implicit none
! LOCAL variables that are not in modules

  SAVE

  real(kind=rk), allocatable :: ralag_1d(:), ralagold_1d(:)
  real(kind=rk), allocatable :: ralag_3d(:,:,:), ralagold_3d(:,:,:)

#ifdef CFC_TRANSPORT2
  real(kind=rk), allocatable ::  thphmean(:,:,:,:),thphmass(:,:,:,:)
#else
  real(kind=rk), allocatable :: thphmean(:,:,:),thphmass(:,:,:,:)
#endif /*CFC_TRANSPORT2*/


#ifdef CFC_TRANSPORT2
  real(kind=rk), allocatable :: sthe(:)
#endif

  integer(kind=ik), allocatable :: jmin(:),jmax(:)

  integer(kind=ik), allocatable :: kmin(:),kmax(:)
contains

!>
!> \par This subroutine allocates the arrays from module multigrid_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine allocate_multigrid_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use mo_mpi

    use configure

    implicit none
    integer(kind=ik)       :: istat
    real(kind=rk)          :: mem_global, mem_local

    allocate(jmin(config%nystrt:config%nymom),jmax(config%nystrt:config%nymom), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module multigrid_rt 3 failed")
    end if

    mem_local = (config%nymom-config%nystrt+1)*4._rk * 2._rk

    allocate(kmin(config%nztra),kmax(config%nztra), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module multigrid_rt 4 failed")
    end if

    mem_local = mem_local + config%nztra*4._rk *2._rk

#ifdef CFC_TRANSPORT2
    allocate(sthe(config%qy), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module multigrid_rt 4a failed")
    end if

    mem_local = mem_local + config%qy*8._rk
#endif

    if (config%use_multid_collapse) then
       allocate(ralag_3d (-1:config%imaxp +1,config%nystrt:config%nymom,config%nztra), &
                ralagold_3d(-1:config%imaxp +1,config%nystrt:config%nymom,config%nztra), stat=istat)
    else

       allocate(ralag_1d   (-1:config%imaxp +1), &
                ralagold_1d(-1:config%imaxp +1), stat=istat)

    endif ! multid_collapse

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module multigrid_rt 1 failed")
    end if

    if (config%use_multid_collapse) then
       mem_local = mem_local + (config%imaxp +1+2)*(config%nymom-config%nystrt)*config%nztra*8._rk*2._rk
    else
       mem_local = mem_local + (config%imaxp +1+2)*8._rk*2._rk
    endif

#ifdef CFC_TRANSPORT2
    allocate(thphmean(config%qx,qy_s:qy_e,qz_s:qz_e,0:1),    &
             thphmass(config%qx,qy_s:qy_e,qz_s:qz_e,0:1),    &
             stat=istat)

    mem_local = mem_local + config%qx*(qy_e-qy_s+1)*        &
                           (qy_e-qz_s+1)*2*8._rk * 2._rk
#else /* CFC_TRANSPORT2 */
    allocate(thphmean(qy_s:qy_e,qz_s:qz_e,0:1),              &
             thphmass(config%qx,qy_s:qy_e,qz_s:qz_e,0:1),    &
             stat=istat)
    mem_local = mem_local + config%qx*(qy_e-qy_s+1)*(qy_e-qz_s+1)*2*8._rk  + &
                               (qy_e-qy_s+1)*(qy_e-qz_s+1)*2*8._rk
#endif /* CFC_TRANSPORT2 */
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module multigrid_rt 2 failed")
    end if

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

     call print_memory_alloc(mem_local, mem_global, "multigrid_rt")

  end subroutine allocate_multigrid_rt
!>
!> \par This subroutine deallocates the arrays from module multigrid_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine deallocate_multigrid_rt
    use precision
    use abort

    use mo_mpi
    use configure

    implicit none
    integer(kind=ik) :: istat

    deallocate(jmin, jmax, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module multigrid_rt 3 failed")
    end if

    deallocate(kmin, kmax, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module multigrid_rt 4 failed")
    end if

#ifdef CFC_TRANSPORT2
    deallocate(sthe, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module multigrid_rt 4a failed")
    end if
#endif

    if (config%use_multid_collapse) then
       deallocate(ralag_3d, ralagold_3d, stat=istat)
    else
       deallocate(ralag_1d, ralagold_1d, stat=istat)
    endif
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module multigrid_rt 1 failed")
    end if

#ifdef CFC_TRANSPORT2
    deallocate(thphmean, thphmass, stat=istat)
#else /* CFC_TRANSPORT2 */
    deallocate(thphmean, thphmass, stat=istat)
#endif /* CFC_TRANSPORT2 */
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module multigrid_rt 2 failed")
    end if


  end subroutine deallocate_multigrid_rt
end module multigrid_rt

! -----------------------------------------------------------------
module averagrid_rt
  use precision

  implicit none
! LOCAL variables that are not in modules

  SAVE

  real(kind=rk), allocatable ::  rqlag(:), gelag(:,:,:), xhlag(:,:,:), &
                                 felag(:,:,:), xjlag(:,:,:), fql_rd(:,:)

contains 
!>
!> \par This subroutine allocates the arrays from module averagrid_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine allocate_averagrid_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(rqlag(0:config%imaxp +1), gelag(-1:config%imaxp +1  ,config%iemax,config%isma), &
             xhlag(-1:config%imaxp +1  ,config%iemax,config%isma),                  &
             felag( 0:config%imaxp +1  ,config%iemax,config%isma),                  &
             xjlag( 0:config%imaxp +1+1,config%iemax,config%isma),                  &
             fql_rd(config%iemax,config%isma), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module averagrid_rt 1 failed")
    end if

    mem_local = (config%imaxp +1+1) *8._rk     + &
                (config%imaxp +1+2)*config%iemax*config%isma*8._rk*3._rk + &
                (config%imaxp +1+1)*config%iemax*config%isma*8._rk       + &
                config%iemax*config%isma*8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "averagrid_rt")

  end subroutine allocate_averagrid_rt
!>
!> \par This subroutine deallocates the arrays from module averagrid_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine deallocate_averagrid_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

    deallocate(rqlag, gelag, xhlag, felag, xjlag, fql_rd, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module averagrid_rt 1 failed")
    end if
  end subroutine deallocate_averagrid_rt

end module averagrid_rt

! -----------------------------------------------------------------
module sourceterms_rt
  use precision

  implicit none
! LOCAL variables that are not in modules

  SAVE
  real(kind=rk), allocatable ::  selag(:,:,:),   &
                                 smlag(:,:,:),   &
                                 sylag(:,:,:,:), &
                                 sylag_rd(:,:,:)


#ifdef FCNC_CALC
  real(kind=rk), allocatable :: selag_n(:,:,:), &
                                smlag_n(:,:,:), &
                                sylag_n(:,:,:)
#endif /* FCNC_CALC */

contains
!>
!> \par This subroutine allocates the arrays from module sourceterms_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine allocate_sourceterms_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use mo_mpi

    use configure
    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(selag(-4:config%imaxp +1+4,nymoms:nymome,nzmoms:nzmome),   &
             smlag(-4:config%imaxp +1+4,nymoms:nymome,nzmoms:nzmome),   &
             sylag(-4:config%imaxp +1+4,nymoms:nymome,nzmoms:nzmome,5), &
             sylag_rd(-4:config%imaxp +1+4,nymoms:nymome,nzmoms:nzmome), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module sourceterms_rt 1 failed")
    end if

    mem_local = (config%imaxp +1+8+1)*(nymome-nymoms)*(nzmome-nzmoms)*8._rk*8._rk


#ifdef FCNC_CALC

    allocate(selag_n(-4:config%imaxp +1+4,nymoms:nymome,nzmoms:nzmome), &
             smlag_n(-4:config%imaxp +1+4,nymoms:nymome,nzmoms:nzmome), &
             sylag_n(-4:config%imaxp +1+4,nymoms:nymome,nzmoms:nzmome), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module sourceterms_rt 2 failed")
    end if

    mem_local = mem_local + (config%imaxp +1+8+1)*(nymome-nymoms)*(nzmome-nzmoms)*8._rk*3._rk

#endif /* FCNC_CALC */

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "sourceterms_rt")

  end subroutine allocate_sourceterms_rt
!>
!> \par This subroutine deallocates the arrays from module sourceterms_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine deallocate_sourceterms_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

    deallocate(selag, smlag, sylag, sylag_rd, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module sourceterms_rt 1 failed")
    end if


#ifdef FCNC_CALC
    deallocate(selag_n, smlag_n,  sylag_n, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module sourceterms_rt 2 failed")
    end if
#endif /* FCNC_CALC */

  end subroutine deallocate_sourceterms_rt

end module sourceterms_rt
! -----------------------------------------------------------------
module stressten_rt
  use precision

  implicit none
! LOCAL variables that are not in modules

  SAVE
  real(kind=rk), allocatable :: dnulag(:,:,:,:), &
                                gnulag(:,:,:,:), &
                                enulag(:,:,:,:), &
                                fnulag(:,:,:,:), &
                                pnulag(:,:,:,:)

contains
!>
!> \par This subroutine allocates the arrays from module stressten_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine allocate_stressten_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use mo_mpi

    use configure
    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(dnulag(-4:config%imaxp +1+4,config%isma,nymoms:nymome,nzmoms:nzmome), &
             gnulag(-4:config%imaxp +1+4,config%isma,nymoms:nymome,nzmoms:nzmome), &
             enulag(-4:config%imaxp +1+4,config%isma,nymoms:nymome,nzmoms:nzmome), &
             fnulag(-4:config%imaxp +1+4,config%isma,nymoms:nymome,nzmoms:nzmome), &
!             fnulag(-4:config%imaxp +1+4,config%isma,nymom,nzmoms:nzmome), &
             pnulag(-4:config%imaxp +1+4,config%isma,nymoms-1:nymome+1,nzmoms:nzmome), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module stressten_rt failed")
    end if

    mem_local = (config%imaxp +1+8+1)*config%isma*(nymome-nymoms+1)*(nzmome-nzmoms*1)*8._rk*4._rk +&
                (config%imaxp +1+8+1)*config%isma*(nymome-nymoms+3)*(nzmome-nzmoms*1)*8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "stressten_rt")

  end subroutine allocate_stressten_rt
!>
!> \par This subroutine deallocates the arrays from module stressten_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine deallocate_stressten_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

    deallocate(dnulag, gnulag, enulag, fnulag, pnulag, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module stressten_rt failed")
    end if

  end subroutine deallocate_stressten_rt

end module stressten_rt
! -----------------------------------------------------------------
module lagradq_rt
  use precision

  implicit none
! LOCAL variables that are not in modules

  SAVE

  real(kind=rk), allocatable, dimension(:,:,:,:,:) :: xlag, wlag, elag
#ifdef BLUEGENE
  real(kind=rk), allocatable, target, dimension(:,:,:,:,:,:) :: riplag, rimlag
#else
  real(kind=rk), allocatable,         dimension(:,:,:,:,:,:) :: riplag, rimlag
#endif
  real(kind=rk), allocatable :: beythq_2d(:,:,:), beyqth_2d(:,:,:) 

contains
!>
!> \par This subroutine allocates the arrays from module lagradq_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine allocate_lagradq_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use mo_mpi

    use configure
    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(beythq_2d(-1:config%imaxp +1,config%nystrt:config%nymom,nzmoms:nzmome), &
             beyqth_2d( 0:config%imaxp +1,config%nystrt:config%nymom,nzmoms:nzmome), stat=istat)
  if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module lagradq_rt 1 failed")
    end if

    mem_local = (config%imaxp +1+2)*(config%nymom-config%nystrt+1)*(nzmome-nzmoms+1)*8._rk + &
                (config%imaxp +1+1)*(config%nymom-config%nystrt+1)*(nzmome-nzmoms+1)*8._rk


    if ( (.not.(use_mpi))) then
       allocate(xlag(-1:config%imaxp*2+3,config%iemax,config%isma,config%nystrt:config%nymom,config%nztra), &
                wlag(-1:config%imaxp*2+3,config%iemax,config%isma,config%nystrt:config%nymom,config%nztra), stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_transport_arrays(): allocation of module lagradq_rt failed")
       end if

       mem_local = mem_local + (config%imaxp*2+3+2)*config%iemax*config%isma*(config%nymom-config%nystrt+1)* &
                   config%nztra*8._rk * 2._rk

       allocate(elag(-1:config%imaxp*2+3,config%iemax,2*config%isma,config%nystrt:config%nytra,config%nztra), &
                riplag(-1:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma, &
                config%nystrt:config%nytra,config%nztra), stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_transport_arrays(): allocation of module lagradq_rt 2 failed")
       end if

       mem_local = mem_local + (config%imaxp*2+3+2)*config%iemax*2*config%isma*(config%nytra-config%nystrt+1)* &
                   config%nztra*8._rk + &
                   (config%imaxp +1+2)*(config%imaxp +1-config%cmin+1)*config%iemax*config%isma*(config%nytra-config%nystrt+1)* &
                   config%nztra*8._rk


       allocate(rimlag(-1:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma,config%nystrt:config%nytra, &
                config%nztra), stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_transport_arrays(): allocation of module lagradq_rt 3 failed")
       end if
           mem_local = mem_local + &
                   (config%imaxp +1+2)*(config%imaxp +1-config%cmin+1)*config%iemax*config%isma*(config%nytra-config%nystrt+1)* &
                   config%nztra*8._rk

    else ! USE_MPI

    !> \todo Attention! As soon as adveclat is called in the third dimension, one has to change
    !>       the dimension of xlag, wlag, xlags and wlags  to nzmoms-1:nzmome+1
      allocate(xlag(-1:config%imaxp*2+3,config%iemax,config%isma,nymoms-1:nymome+1,nzmoms:nzmome),    &
               wlag(-1:config%imaxp*2+3,config%iemax,config%isma,nymoms-1:nymome+1,nzmoms:nzmome),         &
               elag(-1:config%imaxp*2+3,config%iemax,2*config%isma,nymoms:nymome,nzmoms:nzmome),  &
!> \todo riplag and rimlag should be running from a new variable (mod_mpi!!!) nytras:nytrae 
!>       in order to be consistent with NON-MPI case 
               riplag(-1:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma,nymoms:nymome,nzmoms:nzmome),&
               rimlag(-1:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma,nymoms:nymome,nzmoms:nzmome),&
               stat=istat)
      if (istat .ne. 0) then
         raise_abort("allocate_transport_arrays(): allocation of module lagradq_rt failed")
      end if

           mem_local = mem_local + &
                   (config%imaxp*2+3*2)*config%iemax*config%isma*(nymome-nymoms+3)*(nzmome-nzmoms+1)*8._rk*2._rk + &  
                (config%imaxp*2+3+2)*config%iemax*2*config%isma*(nymome-nymoms+1)*(nymome-nymoms+1)*8._rk +        &
                (config%imaxp +1+2)*(config%imaxp +1-config%cmin+1)*config%iemax*config%isma*(nymome-nymoms+1) *                   &
                (nzmome-nzmoms+1)*8._rk*2._rk


   endif ! use_mpi


    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "lagradq_rt")

  end subroutine allocate_lagradq_rt
!>
!> \par This subroutine deallocates the arrays from module lagradq_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine deallocate_lagradq_rt
    use precision
    use abort

    use mo_mpi

    implicit none
    integer(kind=ik) :: istat


    deallocate(beythq_2d,  beyqth_2d, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module lagradq_rt 1 failed")
    end if

    deallocate(xlag, wlag, elag, riplag, rimlag, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module lagradq_rt failed")
    end if

  end subroutine deallocate_lagradq_rt

end module lagradq_rt
! -----------------------------------------------------------------      

module lagback_rt
!
!    arrays containing quantities on all transport-grids
!c
!c   wwlag: material-coefficients
!c   cplag: LTE-Composition 
!c   hylag: helper for output, only for debugging
  use precision

  implicit none
! LOCAL variables that are not in modules

  SAVE

  real(kind=rk), allocatable, dimension(:,:,:,:,:,:) :: wwlag 



! dimensions: (see also interactions.F90 and module matkoeff_rt)
! wwlag(:,:,:,:,:,0) : neutrino energy equilibrium density
! wwlag(:,:,:,:,:,1) : absorption (summ of all rates in the code)
! wwlag(:,:,:,:,:,2) : scattering (summ of all rates in the code)
! wwlag(:,:,:,:,:,3) : inelastic scattering (summ of all rates in the code)
! wwlag(:,:,:,:,:,4) : first moment inelastic scattering 
! wwlag(:,:,:,:,:,5) : second moment inelastic scattering 
! wwlag(:,:,:,:,:,6) : absorption on nuclei (only LMS rates)
! wwlag(:,:,:,:,:,7) : absorption on nuclei (only BRUENN rates)
! wwlag(:,:,:,:,:,8) : absorption on nucleons (only FFN rates)
! wwlag(:,:,:,:,:,9) : absorption on nucleons (only KJT rates)

!     &     hylag(-1:config%imaxp +1+1,config%nystrt:config%nymom,config%nztra,0:4)
!     &     cplag(0:config%imaxp +1,config%nystrt:config%nymom,config%nztra,qn),
!     &     rhlag(0:config%imaxp +1,config%nystrt:config%nymom,config%nztra),
!     &     ttlag(0:config%imaxp +1,config%nystrt:config%nymom,config%nztra),
!     &     xhlag(0:config%imaxp +1,config%nystrt:config%nymom,config%nztra),
!     &     zalag(0:config%imaxp +1,config%nystrt:config%nymom,config%nztra,2)


contains 
!>
!> \par This subroutine allocates the arrays from module lagback_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine allocate_lagback_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use mo_mpi
    use configure
    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    if ( (.not.(use_mpi))) then
       allocate(wwlag(0:config%imaxp +1,config%iemax,config%isma,config%nystrt:config%nymom,config%nztra,0:9), stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_transport_arrays(): allocation of module lgback_rt failed")
       end if
       mem_local = (config%imaxp +1+1)*config%iemax*config%isma*(config%nymom-config%nystrt+1)*config%nztra*10*8._rk

    else
       allocate(wwlag(0:config%imaxp +1,config%iemax,config%isma,nymoms:nymome,nzmoms:nzmome,0:9), stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_transport_arrays(): allocation of module lgback_rt failed")
       end if
       mem_local = (config%imaxp +1+1)*config%iemax*config%isma*(nymome-nymoms+1)*(nzmome-nzmoms+1)*10*8._rk

    endif ! use_mpi

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "lagback_rt")

  end subroutine allocate_lagback_rt
!>
!> \par This subroutine deallocates the arrays from module lagback_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine deallocate_lagback_rt
    use precision
    use abort

    use mo_mpi

    implicit none
    integer(kind=ik) :: istat

    deallocate(wwlag, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module lgback_rt failed")
    end if

  end subroutine deallocate_lagback_rt
end module lagback_rt
! -----------------------------------------------------------------      
module lagrold_rt

!
! arrays
!   contain the radiation field on all transport-grids at the 
!     beginning of the transport-timestep.
!   this has to be restored on all grids if the transport-timestep is
!     not successful on any of the grids !!
!
! used in restsec, copysec

  use precision

  implicit none
! LOCAL variables that are not in modules

  SAVE

  real(kind=rk),allocatable :: xlags(:,:,:,:,:),wlags(:,:,:,:,:), &
                               riplags(:,:,:,:,:,:),rimlags(:,:,:,:,:,:)

contains
!>
!> \par This subroutine allocates the arrays from module lagrold_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine allocate_lagrold_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use mo_mpi

    use configure

    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    if ( (.not.(use_mpi))) then


       allocate(xlags(-1:config%imaxp*2+3,config%iemax,config%isma,config%nystrt:config%nymom,config%nztra), &
                wlags(-1:config%imaxp*2+3,config%iemax,config%isma,config%nystrt:config%nymom,config%nztra), &
                riplags(-1:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma, &
                config%nystrt:config%nytra,config%nztra), stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_transport_arrays(): allocation of module lagrold_rt failed")
       end if

       mem_local = (config%imaxp*2+3+2)*config%iemax*config%isma*(config%nymom-config%nystrt+1)*config%nztra*8._rk*2._rk + &
                    (config%imaxp +1+2)*(config%imaxp +1-config%cmin+1)*config%iemax*config%isma* &
                    (config%nytra-config%nystrt+1)*config%nztra*8._rk

       allocate(rimlags(-1:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma, &
                        config%nystrt:config%nytra,config%nztra),stat=istat)
       if (istat .ne. 0) then
          raise_abort("allocate_transport_arrays(): allocation of module lagrold_rt 2 failed")
       end if

       mem_local = mem_local + (config%imaxp +1+2)*(config%imaxp +1-config%cmin+1)*config%iemax*config%isma* &
                   (config%nytra-config%nystrt+1)*config%nztra*8._rk


    else ! USE_MPI

       allocate( xlags(-1:config%imaxp*2+3,config%iemax,config%isma,nymoms-1:nymome+1,nzmoms:nzmome),    &
                 wlags(-1:config%imaxp*2+3,config%iemax,config%isma,nymoms-1:nymome+1,nzmoms:nzmome),         &
                 riplags(-1:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma,nymoms:nymome,nzmoms:nzmome),&
                 rimlags(-1:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma,nymoms:nymome,nzmoms:nzmome),&
                 stat=istat)

       if (istat .ne. 0) then
          raise_abort("allocate_transport_arrays(): allocation of module lagrold_rt failed")
       end if

       mem_local = (config%imaxp*2+3+2)*config%iemax*config%isma*(nymome-nymoms+3)*(nzmome-nzmoms+1)*8._rk*2._rk + &
                   (config%imaxp +1+2)*(config%imaxp +1-config%cmin+1)*config%iemax*config%isma*(nymome-nymoms+1) * &
                   (nzmome-nzmoms+1)*8._rk*2._rk

    endif ! use_mpi

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local
    
    call print_memory_alloc(mem_local, mem_global, "lagrold_rt")

  end subroutine allocate_lagrold_rt
!>
!> \par This subroutine allocates the arrays from module lagrold_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
  subroutine deallocate_lagrold_rt
    use precision
    use abort

    implicit none
    integer(kind=ik) :: istat

    deallocate(xlags, wlags, riplags, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module lagrold_rt failed")
    end if

    deallocate(rimlags,stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module lagrold_rt 2 failed")
    end if


  end subroutine deallocate_lagrold_rt

end module lagrold_rt
