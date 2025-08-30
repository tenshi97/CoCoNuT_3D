module ephhlp_hy
  use precision

  implicit none
!c LOCAL variables that are not in modules
  save

!  real(kind=rk) :: ephi(0:q)
  real(kind=rk), allocatable :: ephi(:)

#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp threadprivate (ephi)
#endif

contains

!>
!> \par This subroutine allocates the arrays from module ephhlp_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>
!> \endverbatim
!>
  subroutine allocate_ephhlp_hy(mem_global)
    use precision
    use print_stdout_mod, only : print_memory_alloc
    use abort

    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
    allocate(ephi(0:config%q), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module ephhlp_hy 1 failed")
    end if

#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
    mem_local = (config%q+1)*8._rk * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "ephhlp_hy")
    
  end subroutine allocate_ephhlp_hy

!>
!> \par This subroutine deallocates the arrays from module ephhlp_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  
  subroutine deallocate_ephhlp_hy
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat

#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp parallel private(istat)
#endif
    deallocate(ephi, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of module ephhlp_hy 1 failed")
    end if 
#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D_TRA))
!$omp end parallel
#endif
  end subroutine deallocate_ephhlp_hy

end module ephhlp_hy
