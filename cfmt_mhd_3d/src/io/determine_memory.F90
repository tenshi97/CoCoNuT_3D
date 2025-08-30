module meminfo

  use precision
  use memory_c

  implicit none

  private

  public :: init_memory_control, meminfo_start, meminfo_stop, &
            cleanup_meminfo, meminfo_flag

#ifdef SGE_QUEUE
  public :: query_sge_memory  
#endif  

  integer(kind=ik) :: memtot, conversion_factor, pagesize

  logical :: meminfo_flag
  
  integer(kind=ik), save :: instances_start=0
  integer(kind=ik), save :: instances_stop=0

  logical, save :: init_flag = .true.

  type meminfo_vals
     
  character(len=255)  :: descriptor
  integer(kind=ik)    :: memtotal, memfree
  integer(kind=ik)    :: max_mem, heap_mem, heap_avail, heap_max
  integer(kind=ik)    :: stack_mem, stack_avail, mempercore
  integer(kind=ik)    :: unit_conv, memused

  end type meminfo_vals

  type(meminfo_vals), allocatable :: val_of_instance_start(:)
  type(meminfo_vals), allocatable :: val_of_instance_stop(:)
  type(meminfo_vals), allocatable :: tmp(:)


contains

  subroutine cleanup_meminfo
    use precision
    use abort
    use print_stdout_mod

    implicit none

    integer(kind=ik) :: istat

    deallocate(val_of_instance_start, val_of_instance_stop, stat=istat)

    if (istat .ne. 0) then
       raise_abort("error when deallocating  meminfo")
    endif

    if (allocated(tmp)) then
       deallocate(tmp, stat=istat)
       if (istat .ne. 0) then
          raise_abort("error when deallocating  meminfo")
       endif
    endif
  end subroutine cleanup_meminfo


#ifdef SGE_QUEUE
  subroutine query_sge_memory    
    use mo_mpi
    use print_stdout_mod
    implicit none

    call printit_taskX(0,"HOSTNAME                ARCH         NCPU  LOAD  MEMTOT  MEMUSE  SWAPTO  SWAPUS")
    call system("qhost | grep "// hostnames(myproc)%procnameg(1:procname_lengthg(myproc)) )

  end subroutine query_sge_memory
#endif

  subroutine init_memory_control
    use print_stdout_mod   
    use precision
    
    implicit none
    
    integer(kind=ik) :: memused

    integer(kind=ik) :: rusage_max_mem, rusage_heap_mem, rusage_heap_avail, &
                        rusage_heap_max
    integer(kind=ik) :: rusage_stack_mem, rusage_stack_avail, rusage_mempercore

    rusage_max_mem     = 0_ik
    rusage_heap_mem    = 0_ik
    rusage_heap_avail  = 0_ik
    rusage_heap_max    = 0_ik
    rusage_stack_mem   = 0_ik
    rusage_stack_avail = 0_ik
    rusage_mempercore  = 0_ik

    call get_pagesize(pagesize)
    
    memtot=memavail()
    
#ifdef LINUX
    ! to get used memory one can read /proc/self/stat
    call read_proc_self_stat(memused)
#endif


    call printit_taskX(0," ")
    call printit_taskX(0,"Memory information")
    call printit_taskX(0," ")
    call printit_taskX(0,"Pagesize of the system: [kB]",pagesize)
    call printit_taskX(0,"Total memory available: [MB]",memtot)

    call printit_taskX(0," ")
#ifndef IBM
    call printit_taskX(0,"Information from malloc_stat:")
    call print_malloc_stat("init_memory_control")
#endif

#ifdef IBM

    call get_rusage_info(rusage_max_mem, rusage_heap_mem, rusage_heap_avail,   &
                        rusage_heap_max, rusage_stack_mem, rusage_stack_avail, &
                        rusage_mempercore)

    
    call printit_taskX(0," ")
    call printit_taskX(0,"Information from rusage:")
    call printit_taskX(0," MAX_MEM      [MB]: ",rusage_max_mem)
    call printit_taskX(0," HEAP_MEM     [MB]: ",rusage_heap_mem)
#ifdef BLUEGENE
    call printit_taskX(0," HEAP_AVAIL   [MB]: ",rusage_heap_avail)
#endif
    call printit_taskX(0," STACK_MEM    [MB]: ",rusage_stack_mem)
#ifdef BLUEGENE
    call printit_taskX(0," STACK AVAIL  [MB]: ",rusage_stack_avail)
    call printit_taskX(0," MEM PER CORE [MB]: ",rusage_mempercore)
#endif
    call printit_taskX(0," ")
#endif /* IBM */

  end subroutine init_memory_control


  function memavail( ) result(memtot)

    use precision
    implicit none
    
    character(len=50) :: input_string, input_string2
    integer(kind=ik)  :: input
    integer(kind=ik)  :: memtot
    
#if defined(LINUX) | defined(BLUEGENE)
!    integer(kind=ik)  :: i,label_length
    integer(kind=ik)  :: conversion_factor
#endif /* LINUX */


#if defined(LINUX) || defined(BLUEGENE)

  ! on an LINUX system one can obtain the total memory from
  ! the system /proc/meminfo file
    open(11,file='/proc/meminfo',status='old',form='formatted')
!    do i=1,1
       read(11,*,end=120) input_string, input, input_string2
       
       if (input_string2 .eq. 'kB') then
          conversion_factor = 1024
       else if (input_string2 .eq. 'MB') then
          conversion_factor = 1
       else 
          conversion_factor = 1 
       end if

!    enddo
    
    memtot = input / conversion_factor
120 close(11)
  
  
#endif /* LINUX */


#if defined(IBM) && !(defined(BLUEGENE))
    call ibm_mem_free(memtot)
#endif /* IBM && !BLUEGENE */
     
#ifdef NEC
    memtot = 0
#endif


end function memavail

#if (defined(LINUX))
  subroutine read_proc_self_stat(memused)

    use print_stdout_mod

    use precision

    implicit none


    integer(kind=ik) :: ind, j, nmem
    character(len=300) :: dir, dir2, file


    integer(kind=ik), intent(out) :: memused


    file='/proc/self/stat'

    open(unit=11,file=file,form='formatted')
    read(11,'(A300)') dir
    close(11)

    ind=300
    j=0
    do while(j<23)
       ind=index(dir,' ')
       dir2=dir(ind+1:300)
       j=j+1
       dir=dir2
    enddo
    
    ind=index(dir,' ')
    dir2=dir(1:ind)
    read(dir2,'(I12)') nmem   

    memused = (nmem*pagesize)/1024/1024 ! MB

  end subroutine read_proc_self_stat

#endif /* LINUX */



  subroutine meminfo_start(label)
    use print_stdout_mod
    use precision
    implicit none
    
    character*(*),intent(in) :: label
!    character(len=50) :: input_string, input_string2
    
!    integer(kind=ik) :: input
    integer(kind=ik) :: label_length
 
    

    instances_start=instances_start+1
    
    ! ugly implementation of a push operation
    if(instances_start .eq. 1) then
       allocate(val_of_instance_start(instances_start))
    else
       allocate(tmp(instances_start-1))
       tmp(:)=val_of_instance_start(1:instances_start-1)
       deallocate(val_of_instance_start)
       allocate(val_of_instance_start(instances_start))
       val_of_instance_start(1:instances_start-1)=tmp(:)
       deallocate(tmp)
    endif

    label_length=len_trim(label)      

    val_of_instance_start(instances_start)%descriptor(1:label_length)=label(1:label_length)
    val_of_instance_start(instances_start)%descriptor(label_length+1:)=" "
    
    if (init_flag) then
       val_of_instance_start(1)%memtotal = memtot ! from module
       init_flag = .false.
    else
       val_of_instance_start(instances_start)%memtotal = val_of_instance_start(1)%memtotal 
    endif

  
    val_of_instance_start(instances_start)%memused     = 0_ik
    val_of_instance_start(instances_start)%max_mem     = 0_ik
    val_of_instance_start(instances_start)%heap_mem    = 0_ik
    val_of_instance_start(instances_start)%heap_avail  = 0_ik
    val_of_instance_start(instances_start)%heap_max    = 0_ik
    val_of_instance_start(instances_start)%stack_mem   = 0_ik
    val_of_instance_start(instances_start)%stack_avail = 0_ik
    val_of_instance_start(instances_start)%mempercore  = 0_ik
    val_of_instance_start(instances_start)%memfree     = 0_ik



#if ((defined(LINUX)) && !defined(BLUEGENE))
    call read_proc_self_stat(val_of_instance_start(instances_start)%max_mem)
#endif


#if (defined(IBM)) || (defined(BLUEGENE))
    call get_rusage_info(val_of_instance_start(instances_start)%max_mem,    &
                         val_of_instance_start(instances_start)%heap_mem,   &
                         val_of_instance_start(instances_start)%heap_avail, &
                         val_of_instance_start(instances_start)%heap_max,   &
                         val_of_instance_start(instances_start)%stack_mem,  &
                         val_of_instance_start(instances_start)%stack_avail,&
                         val_of_instance_start(instances_start)%mempercore)
#endif



  end subroutine meminfo_start

  subroutine meminfo_stop(label)
    use print_stdout_mod
    use precision
    implicit none

    character*(*),intent(in) :: label
    character(len=50) :: input_string, input_string2
    integer(kind=ik) :: input
    integer(kind=ik) :: label_length
    integer(kind=ik) :: val,i, nref
    



    instances_stop=instances_stop+1


    ! ugly implementation of a push operation
    if(instances_stop .eq. 1) then
       allocate(val_of_instance_stop(instances_stop))
    else
       allocate(tmp(instances_stop-1))
       tmp(:)=val_of_instance_stop(1:instances_stop-1)
       
       deallocate(val_of_instance_stop)
       allocate(val_of_instance_stop(instances_stop))
       val_of_instance_stop(1:instances_stop-1)=tmp(:)
       deallocate(tmp)
    endif
    
    label_length=len_trim(label) 


    val_of_instance_stop(instances_stop)%descriptor(1:label_length)=label(1:label_length)
    val_of_instance_stop(instances_stop)%descriptor(label_length+1:)=" "
    

    val_of_instance_stop(instances_stop)%memtotal = val_of_instance_start(1)%memtotal 


    val_of_instance_stop(instances_stop)%memused     = 0_ik
    val_of_instance_stop(instances_stop)%max_mem     = 0_ik
    val_of_instance_stop(instances_stop)%heap_mem    = 0_ik
    val_of_instance_stop(instances_stop)%heap_avail  = 0_ik 
    val_of_instance_stop(instances_stop)%heap_max    = 0_ik
    val_of_instance_stop(instances_stop)%stack_mem   = 0_ik
    val_of_instance_stop(instances_stop)%stack_avail = 0_ik
    val_of_instance_stop(instances_stop)%mempercore  = 0_ik
    val_of_instance_stop(instances_stop)%memfree     = 0_ik

#if ((defined(LINUX)) && !defined(BLUEGENE))
    call read_proc_self_stat(val_of_instance_stop(instances_stop)%max_mem)
#endif
     



#if (defined(IBM)) || (defined(BLUEGENE)) 
    call get_rusage_info(val_of_instance_stop(instances_stop)%max_mem,    &
                         val_of_instance_stop(instances_stop)%heap_mem,   &
                         val_of_instance_stop(instances_stop)%heap_avail, &
                         val_of_instance_stop(instances_stop)%heap_max,   &
                         val_of_instance_stop(instances_stop)%stack_mem,  &
                         val_of_instance_stop(instances_stop)%stack_avail,&
                         val_of_instance_stop(instances_stop)%mempercore)

#endif



    do i=1,size(val_of_instance_start(:))
       if (val_of_instance_stop(instances_stop)%descriptor .eq. val_of_instance_start(i)%descriptor) then
          nref=i
       endif
    enddo

#ifdef LINUX
    call print_malloc_stat(val_of_instance_stop(instances_stop)%descriptor(1:label_length))
#endif


    call printit_taskX(0,"Measured block                  ",    &
                         val_of_instance_stop(instances_stop)%descriptor(1:label_length))

    val=val_of_instance_stop(instances_stop)%max_mem
    call printit_taskX(0,"Maximal memory used :              [MB]", val)
#ifdef BLUEGENE
    val=val_of_instance_stop(instances_stop)%heap_max
    call printit_taskX(0,"Maximum heap :                     [MB]", val)
    val=val_of_instance_stop(instances_stop)%heap_avail
    call printit_taskX(0,"Heap available :                   [MB]", val)
#endif
    val=val_of_instance_stop(instances_stop)%heap_mem
    call printit_taskX(0,"Heap memory used :                 [MB]", val)
    val=val_of_instance_stop(instances_stop)%stack_mem
    call printit_taskX(0,"Stack memory used :                [MB]", val)
#ifdef BLUEGENE
    val=val_of_instance_stop(instances_stop)%stack_avail
    call printit_taskX(0,"Stack memory available :           [MB]", val)
    val=val_of_instance_stop(instances_stop)%mempercore
    call printit_taskX(0,"Memory per core :                  [MB]", val)
#endif

    call printit_taskX(0," ")

#ifdef SGE_QUEUE
!    call query_sge_memory
#endif

  end subroutine meminfo_stop



end module meminfo
