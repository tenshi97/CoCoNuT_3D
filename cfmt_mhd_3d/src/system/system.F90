!> \verbatim
!> This module provides machine dependent values, such as the largest
!> possible exponent
!>
!> \endverbatim
!>
!> \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date$
!>   
!> \endverbatim
!>
module machcons
  use precision

  implicit none
! LOCAL variables that are not in modules

  SAVE

  real(kind=rk) :: bigeexp,big10exp,emach,big


contains
!> \verbatim
!> This subroutine  sets the machine dependent values, such as the largest
!> possible exponent
!>
!> \endverbatim
!>
!> \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 527 $
!>   $Date$
!>   
!> \endverbatim
!>
  subroutine set_machcons
    use precision
    use print_stdout_mod

    implicit none
    emach = max(1.e-14_rk,5.0_rk*EPSILON(1.0_rk))
    big = HUGE(1.0_rk)
    bigeexp=log(big)-1.0_rk    ! safe side
    big10exp=log10(big)-1.0_rk    ! safe side

! this write-statement is switched off, since set_machcons is called before init_mpi
!    call printit_taskX(0,"Machine constants: ",emach,big,bigeexp,big10exp)

  end subroutine set_machcons
  
end module machcons
