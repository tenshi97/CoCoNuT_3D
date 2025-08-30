MODULE rshlat_hy

!=======================================================================

  USE precision

  USE configure

  IMPLICIT NONE
  
  SAVE
  
  real(kind=rk), allocatable :: rshlat(:,:,:)
  
contains
!>
!> \par This subroutine allocates the arrays from module rshlat_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine allocate_rshlat_hy
    use precision
    use abort
    use configure
    implicit none

    integer(kind=ik) :: istat
    print *, " Module rshlat_hy_hy allocating arrays..."

    allocate(rshlat(1:config%qx,1:config%qy,1:config%qz), stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): allocation of modul rshlat_hy failed")
    end if
  end subroutine allocate_rshlat_hy
!>
!> \par This subroutine deallocates the arrays from module rshlat_hy
!>  \author A. Marek
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!> 
  subroutine deallocate_rshlat_hy
    use precision
    use abort
    implicit none

    integer(kind=ik) :: istat

    deallocate(rshlat, stat=istat)

    if (istat .ne. 0) then
       raise_abort("allocate_hydro_arrays(): deallocation of modul rshlat_hy failed")
    end if
  end subroutine deallocate_rshlat_hy

END MODULE rshlat_hy
