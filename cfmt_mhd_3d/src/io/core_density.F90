module coredensity

private
public :: core_density

contains



!> \par Simple funciton to get the "core" density
!> 
!> \author Lorenz
!> 
!> \detail
!>   In 1D this is just dentot(5,1,1), in multi-D
!>   it is maxval(dentot(5,:,:))
!> 
!> \verbatim
!>   SVN - Information  
!>   $Revision: 850 $
!>   $Date: 2010-03-05 17:21:18 +0100 (Fri, 05 Mar 2010) $
!> \endverbatim



function core_density() result(rho_c)
  use precision
  use totare_hy, only : dentot

  use mo_mpi

  implicit none
  real(kind=rk) :: rho_c ,rho_rcv
  integer :: ierr

  
  if (use_mpi) then

     rho_c = maxval(dentot(5,:,:))
     rho_rcv = 0._rk
     call MPI_Allreduce(rho_c, rho_rcv, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
     rho_c = rho_rcv

  else

     rho_c = maxval(dentot(5,:,:))
  endif
end function core_density

end module coredensity
