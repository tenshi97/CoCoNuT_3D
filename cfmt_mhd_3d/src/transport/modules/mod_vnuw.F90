!>
!> \par This module provides the neutrino sourceterms
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
module vnuw_hy
  use precision

  implicit none
! LOCAL variables that are not in modules
  save

  real(kind=rk), allocatable, dimension(:,:,:) :: qen, qmo, qmy
  real(kind=rk), allocatable :: qye(:,:,:,:)

contains

!>
!> \par This subroutine allocates the arrays from module vnuw_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 726 $
!>   $Date: 2010-02-02 20:36:46 +0100 (Tue, 02 Feb 2010) $
!>   
!> \endverbatim
!>   
  subroutine allocate_vnuw_hy(mem_global)

    use precision
    use abort

    use mo_mpi
    use print_stdout_mod, only : print_memory_alloc
 
    use configure
   implicit none

    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(qye(config%qx,qy_s:qy_e,qz_s:qz_e,5), qen(config%qx,qy_s:qy_e,qz_s:qz_e), & 
             qmo(config%qx,qy_s:qy_e,qz_s:qz_e), qmy(config%qx,qy_s:qy_e,qz_s:qz_e), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module vnuw_hy 1 failed")
    end if

    mem_local = config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk * 8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "vnuw_hy")
    
  end subroutine allocate_vnuw_hy

!>
!> \par This subroutine allocates the arrays from module vnuw_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 726 $
!>   $Date: 2010-02-02 20:36:46 +0100 (Tue, 02 Feb 2010) $
!>   
!> \endverbatim
!>   
  subroutine deallocate_vnuw_hy

    use precision
    use abort

    implicit none
    integer(kind=ik) :: istat

    deallocate(qye, qen, qmo, qmy, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module vnuw_hy 1 failed")
    end if

  end subroutine deallocate_vnuw_hy

end module vnuw_hy
