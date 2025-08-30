!>
!> \par This module provides the real-kind, integer-kind definitions that are used in VERTEX
!>
!> \param rk    double precision floating point
!> \param rk4   single precision floating point
!> \param ik    the system integer
!>
!> \author A. Marek, MPA, Jan. 2009
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module precision
  use iso_c_binding, only : C_FLOAT, C_DOUBLE, C_INT32_T

  implicit none
      
  integer, parameter :: rk  = C_DOUBLE
  integer, parameter :: rk4 = C_FLOAT
      
  integer, parameter :: ik = C_INT32_T

#ifdef DOUBLE_PRECISION_EOS
  integer, parameter :: rk_eos = C_DOUBLE
#else
  integer, parameter :: rk_eos = C_FLOAT
#endif

end module precision
