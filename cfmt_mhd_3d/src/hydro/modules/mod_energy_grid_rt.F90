!-----------------------------------------------------------------------
!>
!> \par This module provides some variables needed for the energy transport grid
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
module energy_grid_rt
  use precision

  implicit none
! LOCAL varibales that are not in modules

  SAVE

  real(kind=rk), allocatable, dimension(:) :: ener, emid, eh3, eh4, eh2de, &
                                              eh3de, pole, pole1, expol,   &
                                              expol1, e3lft, e3rgt, e4lft, &
                                              e4rgt, exe3lft, exe3rgt,     &
                                              exe4lft, exe4rgt

contains

!>
!> \par This subroutine allocates the arrays from module energy_grid_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine allocate_energy_grid_rt(mem_global)
    use precision
    use abort
    use print_stdout_mod, only : print_memory_alloc

    use configure
    implicit none

    integer(kind=ik) :: istat
    real(kind=rk)    :: mem_global, mem_local

    allocate(ener(config%iemax+2), emid(config%iemax+1),       &
              eh3(config%iemax+1),  eh4(config%iemax+1),       &
            eh2de(0:config%iemax+1), eh3de(0:config%iemax+1),  &
             pole(config%iemax), pole1(config%iemax),          &
            expol(config%iemax), expol1(config%iemax), stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of module energy_grid_rt failed")
    end if
    mem_local = (config%iemax+2)*8._rk *3._rk + (config%iemax+1)*8._rk * 3._rk + &
                (config%iemax) * 8._rk * 4._rk

    mem_local = mem_local/1024._rk/1024._rk

    mem_global = mem_global + mem_local

    call print_memory_alloc(mem_local, mem_global, "energy_grid_rt")

  end subroutine allocate_energy_grid_rt
!>
!> \par This subroutine deallocates the arrays from module energy_grid_rt
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine deallocate_energy_grid_rt
    use precision
    use abort
    implicit none

    integer(kind=ik) :: istat

    deallocate(ener,emid,eh3,eh4, eh2de,eh3de, &
               pole,pole1,expol,expol1, stat=istat)
    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of module energy_grid_rt failed")
    end if
  end subroutine deallocate_energy_grid_rt


end module energy_grid_rt
