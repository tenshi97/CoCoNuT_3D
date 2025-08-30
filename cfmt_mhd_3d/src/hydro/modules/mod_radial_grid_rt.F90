!>
!> \par This module provides some variables needed for the radial transport grid
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module radial_grid_rt
  use precision

  implicit none
! LOCAL variables that are not in modules

  SAVE

  real(kind=rk), allocatable, dimension(:) :: polrq,polrq1,  &
                                              polrqq

  real(kind=rk), allocatable, dimension(:) :: r, rq, dv,dvq,  &
                                              ss, sq,dvalt,   &
                                              dvqalt, dvw,    &
                                              dvwq, dv_raw, dvq_raw
  real(kind=rk), allocatable :: rr(:)

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp threadprivate(polrq,polrq1,polrqq,r,rq,dv,dvq,ss, &
!$omp               sq,dvalt,dvqalt,dvw,dvwq,dv_raw,dvq_raw,rr)
#endif
#endif

contains 
!>
!> \par This subroutine allocates the arrays from module radial_grid_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine allocate_radial_grid_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif

    allocate(polrq(-1:config%imaxp +1),polrq1(0:config%imaxp +1),  &
             polrqq(0:config%imaxp +1), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module radial_grid_rt failed")
    end if

    allocate(r(config%cmin:config%imaxp +1), rq(0:config%imaxp +1), dv(-1:config%imaxp +1), dvq(0:config%imaxp +1), &
             ss(-1:config%imaxp +1), sq(0:config%imaxp +1), dvalt(-1:config%imaxp +1),              &
             dvqalt(0:config%imaxp +1), dvw(-1:config%imaxp +1), dvwq(0:config%imaxp +1),           &
             dv_raw(-1:config%imaxp +1), dvq_raw(0:config%imaxp +1), stat=istat)



    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module radial_grid_rt 2 failed")
    end if

    allocate(rr(-2:config%imaxp +1+1), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module radial_grid_rt 3 failed")
    end if

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif
#endif

    mem_local = (config%imaxp +1+2)*8._rk*6._rk + (config%imaxp +1+1)*8._rk*8._rk + &
                (config%imaxp +1-config%cmin +1)*8._rk + (config%imaxp +1+4)*8._rk * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "radial_grid_rt")

  end subroutine allocate_radial_grid_rt
!>
!> \par This subroutine deallocates the arrays from module radial_grid_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine deallocate_radial_grid_rt
    use precision
    use abort
    use configure
    implicit none

    integer(kind=ik) :: istat

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
#endif
    deallocate(polrq,polrq1,polrqq, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module radial_grid_rt failed")
    end if

    deallocate(r, rq, dv ,dvq, ss, sq, dvalt, dvqalt, dvw, dvwq, dv_raw, &
               dvq_raw, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module radial_grid_rt 2 failed")
    end if

    deallocate(rr, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module radial_grid_rt 3 failed")
    end if

#if !(defined(PROGRAM_remap))
#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp end parallel
#endif
#endif

  end subroutine deallocate_radial_grid_rt

end module radial_grid_rt
