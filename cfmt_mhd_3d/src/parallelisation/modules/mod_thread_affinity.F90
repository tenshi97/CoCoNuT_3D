module thread_affinity

  use precision
  use setup_c

  implicit none
  public :: check_thread_affinity, print_thread_affinity, &
            init_thread_affinity, cleanup_thread_affinity

  private
! integer(kind=ik) :: thread_num
  integer(kind=ik) :: thread_max
  integer(kind=ik) :: process_cpu_id

  integer(kind=ik), allocatable :: cpu_ids(:)

contains

  subroutine init_thread_affinity
    
    use precision
    use abort

    use print_stdout_mod
    implicit none


    integer(kind=ik), external :: OMP_GET_MAX_THREADS
    integer(kind=ik) :: istat

#if defined(OPENMP_HYD)
    thread_max = OMP_GET_MAX_THREADS()
#else
    thread_max = 1
#endif    
    if(.not.(allocated(cpu_ids))) then
       allocate(cpu_ids(thread_max), stat=istat)

       if (istat .ne. 0) then
          raise_abort("Error when allocating init_thread_affinity")
       endif
    endif

  end subroutine init_thread_affinity

  subroutine cleanup_thread_affinity
    
    use precision
    use abort
    implicit none

    integer(kind=ik) :: istat

    if((allocated(cpu_ids))) then
       deallocate(cpu_ids, stat=istat)

       if (istat .ne. 0) then
          raise_abort("Error when deallocating init_thread_affinity")
       endif
    endif

  end subroutine cleanup_thread_affinity


  subroutine check_thread_affinity
    
    use precision
    implicit none
    
    integer(kind=ik) :: thread_cpu_id
    
    integer(kind=ik), external :: OMP_GET_THREAD_NUM
    
    integer(kind=ik) :: actuall_num, i
       

    call get_process_affinity(process_cpu_id)

#ifdef CHECK_THREAD_AFFINITY
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(NONE) &
!$OMP   PRIVATE(i,thread_cpu_id,actuall_num) &
!$OMP   SHARED(thread_max,cpu_ids) &
!$OMP   SCHEDULE(STATIC)
    do i=1,thread_max
       call get_thread_affinity(thread_cpu_id)
       actuall_num=omp_get_thread_num()
       cpu_ids(actuall_num+1)=thread_cpu_id
    enddo
#endif

  end subroutine check_thread_affinity

  subroutine print_thread_affinity
  
    use precision
    use mpi_vertex, only : myproc
    use print_stdout_mod
    implicit none

    integer(kind=ik) :: i
    call printit_taskX(0," ")
    call printit_taskX(0," =================================================")
    call printit_taskX(0,"              Thread affinity ")
    call printit_taskX(0," ")

    if (myproc .eq. 0) then
       write(*,'("In total ",i4," threads on this node or run")') thread_max
       do i=1,thread_max

          write(*,'("Thread ",i4," is running on logical CPU-ID ",i4)') i,cpu_ids(i)
       enddo
    endif
    call printit_taskX(0," ")


    write(97,*) " "
    write(97,*) " ======================================================="
    write(97,*) "           Thread affinity            "
    write(97,*) " "
    write(97,'("In total ",i4," threads on this node or run")') thread_max
    do i=1,thread_max

       write(97,'("Thread ",i4," is running on logical CPU-ID ",i4)') i,cpu_ids(i)
    enddo
    write(97,*) " "

  end subroutine print_thread_affinity


end module thread_affinity
