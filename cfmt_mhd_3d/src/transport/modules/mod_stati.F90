!>
!> \par This module provides some variables needed for the transport
!> \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
module stati_rt

  use precision

  implicit none
!  LOCAL variables thata are not in modules

  save

! for timestep-control 
  real(kind=rk) :: sigsy,sigse
  real(kind=rk), allocatable, dimension(:,:) :: sigmaj, sigmah, &
                                                sigman, sigmaf

  integer(kind=ik) :: irkritj,irkrith,irkritn,irkritf, &
                      iekritj,iskritj,iekritn,iskritn, &
                      iekrith,iskrith,iekritf,iskritf, &
                      ikritsy,ikritse

! for ME-BTE Iteration

  real(kind=rk), allocatable :: sigit(:,:)

  integer(kind=ik) :: ieit,isit
  integer(kind=ik), allocatable :: ikritit(:,:)

! for Newton-Iteration

  real(kind=rk), allocatable :: relxj(:,:,:), relxh(:,:,:), relso(:)

  integer(kind=ik), allocatable :: iirelj(:,:,:),iirelh(:,:,:), irelso(:), &
                                   ierelj(:), isrelj(:),ierelh(:),  &
                                   isrelh(:)

#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp threadprivate(                                         &
!$omp      irkritj,irkrith,irkritn,irkritf,                  &
!$omp      iekritj,iskritj,iekritn,iskritn,                  &
!$omp      iekrith,iskrith,iekritf,iskritf,                  &
!$omp      ikritsy,ikritse,                                  &
!$omp      ikritit,ieit,isit,                                &
!$omp      iirelj,iirelh,irelso,ierelj,isrelj,ierelh,isrelh) 
!$omp threadprivate(                                         &
!$omp      sigmaj,sigmah,sigman,sigmaf,sigsy,sigse,sigit,    &
!$omp      relxj,relxh,relso)
#endif

contains
!>
!> \par This subroutine allocates the arrays from module stati_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
  subroutine allocate_stati_rt(mem_global)

    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none
    integer(kind=ik) :: istat, omp_get_num_threads, thread_num=1_ik
    real(kind=rk)    :: mem_global, mem_local

#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif

    allocate(sigmaj(config%iemax,config%isma), sigmah(config%iemax,config%isma), sigman(config%iemax,config%isma), &
             sigmaf(config%iemax,config%isma), sigit(config%iemax,config%isma), ikritit(config%iemax,config%isma), &
             stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of modul stati_rt 1 failed")
    end if


    allocate(relxj(config%iemax,config%isma,2),relxh(config%iemax,config%isma,2), relso(2), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of modul stati_rt 2 failed")
    end if

    allocate(iirelj(config%iemax,config%isma,2), iirelh(config%iemax,config%isma,2), irelso(2), &
             ierelj(2),isrelj(2),ierelh(2), isrelh(2), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of modul stati_rt 3 failed")
    end if

#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp single
    thread_num = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif

    mem_local = config%iemax*config%isma*8._rk *5._rk + config%iemax*config%isma*4._rk + &
                config%iemax*config%isma*2*8._rk*2._rk +2*8._rk          + &
                config%iemax*config%isma*2*4._rk*2._rk +2*4._rk*5._rk 
    mem_local = mem_local * thread_num

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "stati_hy")

  end subroutine allocate_stati_rt
!>
!> \par This subroutine deallocates the arrays from module stati_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
  subroutine deallocate_stati_rt
    use precision
    use abort

    implicit none
    integer(kind=ik) :: istat

#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp parallel private(istat)
#endif
    deallocate(sigmaj, sigmah, sigman, sigmaf, sigit, ikritit, stat=istat)

    if (istat .ne. 0) then
       raise_abort("deallocate_transport_arrays(): deallocation of modul stati_rt 1 failed")
    end if

    deallocate(relxj, relxh, relso, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of modul stati_rt 2 failed")
    end if

    deallocate(iirelj, iirelh, irelso, ierelj,isrelj,ierelh, isrelh, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): deallocation of modul stati_rt 3 failed")
    end if

#if  defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp end parallel
#endif

  end subroutine deallocate_stati_rt
    


end module stati_rt
