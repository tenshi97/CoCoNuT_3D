module memory_c
  use precision

  private
#ifdef IBM
  public get_rusage_info
#endif
#ifndef IBM
  public print_malloc_stat
#endif
  public get_pagesize


#ifndef IBM
  interface
    subroutine print_malloc_stat_c(name) bind(C, name="print_malloc_stat")
      use iso_c_binding
      character(kind=C_CHAR), dimension(*), intent(in) :: name
    end subroutine
  end interface
#endif
  
  interface
    subroutine get_pagesize_c(pagesize) bind(C, name="get_pagesize")
      use iso_c_binding
      implicit none
      integer(kind=C_INT), intent(out) :: pagesize
    end subroutine
  end interface

#ifdef IBM
  interface
    subroutine get_rusage_info_c(max_mem, heap_mem, heap_avail, heap_max , &
                                 stack_mem, stack_avail, mempercore)       &
                                 bind(C, name="get_rusage_info")
      use iso_c_binding
      implicit none
      integer(kind=C_INT) :: max_mem
      integer(kind=C_INT) :: heap_mem
      integer(kind=C_INT) :: heap_avail
      integer(kind=C_INT) :: heap_max
      integer(kind=C_INT) :: stack_mem
      integer(kind=C_INT) :: stack_avail
      integer(kind=C_INT) :: mempercore
    end subroutine get_rusage_info_c
  end interface
#endif

#if defined(IBM) || defined(BLUEGENE)
  interface
    subroutine ibm_mem_free_c(mem) bind(C, name="ibm_mem_free")
      use iso_c_binding
      implicit none
      integer(kind=C_INT) :: mem
    end subroutine
  end interface
#endif

  contains

#ifndef IBM
    subroutine print_malloc_stat(name)
      use iso_c_binding
      implicit none
      character(*), intent(in) :: name
      call print_malloc_stat_c(trim(name) // C_NULL_CHAR)
    end subroutine
#endif

    subroutine get_pagesize(pagesize)
      use iso_c_binding
      implicit none
      integer(kind=ik), intent(out) :: pagesize
      integer(kind=C_INT) :: pagesize_c
      call get_pagesize_c(pagesize_c)
      pagesize = int(pagesize_c, kind=ik)
    end subroutine
    
#ifdef IBM
    subroutine get_rusage_info(max_mem, heap_mem, heap_avail,   & 
                              heap_max, stack_mem, stack_avail, &
                              mempercore)
      use iso_c_binding
      implicit none
      integer(kind=ik), intent(out) :: max_mem, heap_mem, heap_avail,    &
                                       heap_max, stack_mem, stack_avail, &
                                       mempercore

      integer(kind=C_INT) :: max_mem_c, heap_mem_c, heap_avail_c,    &
                             heap_max_c, stack_mem_c, stack_avail_c, &
                             mempercore_c

      call get_rusage_info_c(max_mem_c, heap_mem_c, heap_avail_c, heap_max_c, &
                           stack_mem_c, stack_avail_c, mempercore_c)
      max_mem       = int(max_mem_c,   kind=ik)
      heap_mem      = int(heap_mem_c,  kind=ik)
      heap_avail    = int(heap_avail_c, kind=ik)
      heap_max      = int(heap_max_c, kind=ik)
      stack_mem     = int(stack_mem_c, kind=ik)
      stack_avail   = int(stack_avail_c, kind=ik)
      mempercore    = int(mempercore_c, kind=ik)
    end subroutine get_rusage_info
#endif

end module memory_c
