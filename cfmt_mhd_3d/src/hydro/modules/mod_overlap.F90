!-----------------------------------------------------------------------
!>
!> \par This module provides some variables needed for the lateral ransport grid
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
module overlap
  use precision

  implicit none
! LOCAL varibales that are not in modules
  
  SAVE

  integer(kind=ik), allocatable ::  nxlb(:),nxub(:),novl(:)

contains

!>
!> \par This subroutine allocates the arrays from module overlap
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine allocate_overlap(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc
    use configure
    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local
    
    allocate(nxlb(config%nystrt:config%nymom), nxub(config%nystrt:config%nymom), novl(config%nystrt:config%nymom), &
             stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module overlap failed")
    end if

    mem_local = (config%nymom-config%nystrt+1)*4._rk * 3._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "overlap")

  end subroutine allocate_overlap
!>
!> \par This subroutine deallocates the arrays from module overlap
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine deallocate_overlap
    use precision
    use abort
    implicit none
    integer(kind=ik) :: istat
      
    deallocate(nxlb,nxub,novl, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module overlap failed")
    end if

  end subroutine deallocate_overlap

end module overlap
