module nusource_data
  use precision
  use mo_mpi

  implicit none
  
  save

  integer(kind=ik):: qqq

  real(kind=rk), allocatable, dimension(:,:,:)  :: syeold,senold,smoold, &
                                                   smyold
  real(kind=rk), allocatable, dimension(:,:,:)  :: syeold1,syeold2, &
                                                   syeold3,syeold4

!  data syeold /qqq*0._rk/, senold /qqq*0._rk/, smoold /qqq*0._rk/, &
!                           smyold /qqq*0._rk/, syeold1 /qqq*0._rk/
!

contains

!>
!> \par This subroutine allocates the arrays from module nusource_data
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>  

  subroutine allocate_nusource_data(mem_global)
    use precision
    use abort
    use mo_mpi
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none
    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local


    qqq=config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)

    allocate(syeold(config%qx,qy_s:qy_e,qz_s:qz_e), senold(config%qx,qy_s:qy_e,qz_s:qz_e), &
             smoold(config%qx,qy_s:qy_e,qz_s:qz_e), smyold(config%qx,qy_s:qy_e,qz_s:qz_e), &
             syeold1(config%qx,qy_s:qy_e,qz_s:qz_e),syeold2(config%qx,qy_s:qy_e,qz_s:qz_e), &
             syeold3(config%qx,qy_s:qy_e,qz_s:qz_e),syeold4(config%qx,qy_s:qy_e,qz_s:qz_e), &
             stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_transport_arrays(): allocation of module arvisco_rt failed")
    end if


    mem_local = config%qx*(qy_e-qy_s+1)*(qz_e-qz_s+1)*8._rk*8._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local
    
    call print_memory_alloc(mem_local, mem_global, "nusource_rt")
    
  end subroutine allocate_nusource_data

end module nusource_data
