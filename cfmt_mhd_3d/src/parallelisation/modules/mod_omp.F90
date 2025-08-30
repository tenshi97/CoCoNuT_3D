!>
!> \verbatim
!> This module provides the OpenMP functionality in the VERTEX
!> code
!>
!>  Author: M. Rampp
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module mod_2d

  use precision
  implicit none

  SAVE


  !$    external OMP_GET_THREAD_NUM
  !$    integer OMP_GET_THREAD_NUM
  integer(kind=ik) :: j_2d,k_2d,jknum,ithrd_omp
#if !(defined(PROGRAM_remap))
#if defined(OPENMP_TRANSPORT) && (defined(OPEN_MP_2D) || defined(OPEN_MP_3D))
!$omp threadprivate (j_2d,k_2d,jknum)
#endif
#endif

end module mod_2d



