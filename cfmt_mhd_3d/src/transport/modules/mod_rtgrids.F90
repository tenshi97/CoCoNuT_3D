!>
!> \par This module provides the arrays that are needed for the tangent ray grid
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>
module tanray_grids_rt
  use precision

  implicit none
! LOCAL variables that are not in modules
  SAVE
  real(kind=rk), allocatable :: amue(:,:), amueq(:,:)
  real(kind=rk), allocatable :: amueold(:,:), amueqold(:,:)

  real(kind=rk), allocatable :: drmue(:,:), drmueq(:,:)

  real(kind=rk), allocatable :: wika_rd(:),wikc_rd(:), &
                                wikb(:,:), &
                                wikd(:,:), &
                                wikaq(:,:), &
                                wikcq(:,:)

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && ((defined(EDD_FACT_2D) || defined(EDD_FACT_3D) ) && !(defined(OPEN_MP_1D)))
!$omp threadprivate ( amue,amueold,drmue,wika_rd,wikc_rd,wikb, &
!$omp                 wikd,wikaq,wikcq)
#endif /* EDD_FACT_3D ||  EDD_FACT_3D */
#endif

contains
!>
!> \par This subroutine allocates the arrays from module tanray_grids_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>

  subroutine allocate_tanray_grids_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

    allocate(amueq(0:config%imaxp +1,config%cmin:config%imaxp +1), amueqold( 0:config%imaxp +1,config%cmin:config%imaxp +1), &
             drmueq(-1:config%imaxp +1,config%cmin:config%imaxp), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module tanray_grids_r 0 failed")
    end if 

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && ((defined(EDD_FACT_2D) || defined(EDD_FACT_3D) ) && !(defined(OPEN_MP_1D)))
!$omp parallel private(istat)
#endif
#endif
    allocate(amue(-1:config%imaxp +1,config%cmin:config%imaxp +1),  &
             amueold(-1:config%imaxp +1,config%cmin:config%imaxp +1),  &
             wika_rd(config%cmin:config%imaxp +1),wikc_rd(config%cmin:config%imaxp +1), &
             wikb(-1:config%imaxp +1,config%cmin:config%imaxp +1), wikd(-1:config%imaxp +1,config%cmin:config%imaxp +1), &
             wikaq(0:config%imaxp +1,config%cmin:config%imaxp +1), wikcq(0:config%imaxp +1,config%cmin:config%imaxp +1), &
             drmue(0:config%imaxp +1,config%cmin:config%imaxp +1),  &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module tanray_grids_r 1 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && ((defined(EDD_FACT_2D) || defined(EDD_FACT_3D) ) && !(defined(OPEN_MP_1D)))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local =(config%imaxp +1+2)*(config%imaxp +1-config%cmin+1)*8._rk*4._rk + &
               (config%imaxp +1+1)*(config%imaxp +1-config%cmin+1)*8._rk*3._rk + &
               (config%imaxp +1-config%cmin+1)*8._rk*2._rk
    mem_local = mem_local * thread_num

    mem_local = mem_local + (config%imaxp +1+1)*(config%imaxp +1-config%cmin+1)*8._rk*2._rk + &
                (config%imaxp +1+2)*(config%imaxp +1-config%cmin+1)*8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "tanray_grids_rt")
    
  end subroutine allocate_tanray_grids_rt
!>
!> \par This subroutine deallocates the arrays from module tanray_grids_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>

  subroutine deallocate_tanray_grids_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

    deallocate(amueq, amueqold,  drmueq, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module tanray_grids_r 0 failed")
    end if
#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && ((defined(EDD_FACT_2D) || defined(EDD_FACT_3D) ) && !(defined(OPEN_MP_1D)))
!$omp parallel private(istat)
#endif
#endif
    deallocate(amue, amueold, wika_rd,wikc_rd, &
               wikb, wikd,  wikaq, wikcq,  drmue, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module tanray_grids_r 1 failed")
    end if
#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && ((defined(EDD_FACT_2D) || defined(EDD_FACT_3D) ) && !(defined(OPEN_MP_1D)))
!$omp end parallel
#endif
#endif
  end subroutine deallocate_tanray_grids_rt

end module tanray_grids_rt

!-----------------------------------------------------------------------

!>
!> \par This module provides the arrays that are needed for the radiation field
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>
module radfield_rt

  use precision

  implicit none
! LOCAL variables that are not in modules

  save

  real(kind=rk), allocatable :: x(:,:,:), w(:,:,:),       &
                                xalt(:,:,:), walt(:,:,:), &
                                ripoldi(:,:,:,:),         &
                                rimoldi(:,:,:,:)
#ifdef BLUEGENE
  real(kind=rk), pointer     :: rip(:,:,:,:), rim(:,:,:,:)
#else
  real(kind=rk), allocatable :: rip(:,:,:,:), rim(:,:,:,:)
#endif
  real(kind=rk), allocatable :: fe(:,:,:), fq(:,:,:), ge(:,:,:), &
                                xjit(:,:,:), xhit(:,:,:), &
                                feold(:,:,:)

  real(kind=rk), allocatable :: xhinnen(:,:), xm1rd(:,:)

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp threadprivate(x,w,xalt,walt,fe,feold,fq,ge,xjit,xhit, &
!$omp               xhinnen,xm1rd)

#if defined(EDD_FACT_2D) || defined(EDD_FACT_3D)
!$omp threadprivate(rip,rim,ripoldi,rimoldi)
#endif /* EDD_FACT_2D  || EDD_FACT_3D */

#endif /* OPEN_MP_2D || OPEN_MP_3D_TRA */
#endif

contains
#ifdef BLUEGENE
  subroutine associate_with_lbound(ptr,array,u1,u2,u3,u4)
    implicit none
    integer(kind=ik), intent(in) :: u1,u2,u3,u4
    real(kind=rk), dimension(u1:,u2:,u3:,u4:), intent(in), target :: array
    real(kind=rk), dimension(:,:,:,:), pointer, intent(inout) :: ptr
    ptr => array
  end subroutine associate_with_lbound
#endif
!>
!> \par This subroutine allocates the arrays from module radfield_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  

  subroutine allocate_radfield_rt(mem_global)

    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik, thread_num_2=1_ik
    real(kind=rk)    :: mem_global, mem_local, mem_local_2

    mem_local = 0._rk
    mem_local_2= 0._rk

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    allocate(x(-1:config%imaxp*2+3,config%iemax,config%isma), w(-1:config%imaxp*2+3,config%iemax,config%isma),       &
             xalt(-1:config%imaxp*2+3,config%iemax,config%isma), walt(-1:config%imaxp*2+3,config%iemax,config%isma), &
             fe(0:config%imaxp +1,config%iemax,config%isma), fq(-1:config%imaxp +1,config%iemax,config%isma),        &
             ge(-1:config%imaxp +1,config%iemax,config%isma), xjit(0:config%imaxp +1+1,config%iemax,config%isma),    &
             xhit(-1:config%imaxp +1,config%iemax,config%isma), feold(0:config%imaxp +1,config%iemax,config%isma),   &
             xhinnen(config%iemax,config%isma), xm1rd(config%iemax,config%isma), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module radfield_rt 1 failed")
    end if 

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif
    mem_local = (config%imaxp*2+3+2)*config%iemax*config%isma*8._rk * 4._rk   + &
                (config%imaxp +1+1)*config%iemax*config%isma*8._rk  * 2._rk   + &
                (config%imaxp +1+2)*config%iemax*config%isma*8._rk  * 4._rk   + &
                config%iemax*config%isma*8._rk * 2._rk

    mem_local = mem_local * thread_num

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && ((defined(OPEN_MP_2D) && defined(EDD_FACT_2D)) || (defined(OPEN_MP_3D_TRA) && defined(EDD_FACT_3D)))
!$omp parallel private(istat)
#endif
#endif

    allocate(ripoldi( 0:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma), &
             rimoldi(-1:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma), &
#ifndef BLUEGENE
             rip( 0:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma), &
             rim(-1:config%imaxp +1,config%cmin:config%imaxp +1,config%iemax,config%isma), & 
#endif
             stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module radfield_rt 2 failed")
    end if

#ifdef BLUEGENE
    nullify(rip)
    nullify(rim)
#endif

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && ((defined(OPEN_MP_2D) && defined(EDD_FACT_2D)) || (defined(OPEN_MP_3D_TRA) && defined(EDD_FACT_3D)))
!$omp single
    thread_num_2 = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

#ifdef BLUEGENE
    mem_local_2 = (config%imaxp +1+1)*(config%imaxp +1-config%cmin+1)*config%iemax*config%isma*8._rk   + &
                  (config%imaxp +1+2)*(config%imaxp +1-config%cmin+1)*config%iemax*config%isma*8._rk
#else
    mem_local_2 = (config%imaxp +1+1)*(config%imaxp +1-config%cmin+1)*config%iemax*config%isma*8._rk * 2._rk  + &
                  (config%imaxp +1+2)*(config%imaxp +1-config%cmin+1)*config%iemax*config%isma*8._rk * 2._rk
#endif

    mem_local = mem_local + mem_local_2

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "radfield_rt")
    
  end subroutine allocate_radfield_rt
!>
!> \par This subroutine allocates the arrays from module radfield_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  

  subroutine deallocate_radfield_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat
#if !(defined(PROGRAM_remap))
#if defined(OPENMP_TRANSPORT) && ( defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    deallocate(x, w, xalt, walt, fe, fq, ge, xjit, xhit, feold, &
               xhinnen, xm1rd, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module radfield_rt 1 failed")
    end if 
#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp end parallel
#endif
#endif

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && ((defined(OPEN_MP_2D) && defined(EDD_FACT_2D)) || (defined(OPEN_MP_3D_TRA) && defined(EDD_FACT_3D)))
!$omp parallel private(istat)
#endif
#endif
    deallocate(ripoldi, rimoldi, rip, rim, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module radfield_rt 2 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && ((defined(OPEN_MP_2D) && defined(EDD_FACT_2D)) || (defined(OPEN_MP_3D_TRA) && defined(EDD_FACT_3D)))
!$omp end parallel
#endif
#endif
  end subroutine deallocate_radfield_rt
end module radfield_rt
!-----------------------------------------------------------------------
!>
!> \par This module provides the quantities which are needed for the cutting of the neutrino spectra
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>
module speccut_rt
  use precision

  use configure
  implicit none
! LOCAL variables that are not in modules

  save

#ifdef SPECCUT_SEVEN
  real(kind=rk), parameter :: specran=1.e-7_rk
#else
  real(kind=rk), parameter :: specran=1.e-12_rk
#endif
!  logical :: lowebin(0:config%imaxp +1,iemax+1,isma,2)
  logical, allocatable :: lowebin(:,:,:,:)

!  integer :: ierandb(0:config%imaxp +1,isma,2)
  integer(kind=ik), allocatable :: ierandb(:,:,:)

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp threadprivate(lowebin,ierandb)
#endif /* OPEN_MP_2D || OPEN_MP_3D_TRA */
#endif

contains

!>
!> \par This subroutine allocates the arrays from module speccut_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  

  subroutine allocate_speccut_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif

    allocate(lowebin(0:config%imaxp +1,config%iemax+1,config%isma,2), ierandb(0:config%imaxp +1,config%isma,2), &
             stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module  speccut_rt1 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local = (config%imaxp +1+1)*(config%iemax+1)*config%isma*2*2._rk + &
                (config%imaxp +1+1)*config%isma*2*4._rk
    mem_local = mem_local * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "speccut_rt")
    
  end subroutine allocate_speccut_rt

!>
!> \par This subroutine deallocates the arrays from module speccut_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  

  subroutine deallocate_speccut_rt
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat


#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
#endif
    deallocate(lowebin, ierandb, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module  speccut_rt1 failed")
    end if

#if !(defined(PROGRAM_remap))
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp end parallel
#endif
#endif
  end subroutine deallocate_speccut_rt

end module speccut_rt
