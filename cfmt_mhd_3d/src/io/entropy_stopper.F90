module stopentropy

private

public :: stop_entropy

contains

!> \par simple function to stop the program at bounce
!> 
!> \author Lorenz
!> 
!> 
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!> \endverbatim
function stop_entropy() result(stopit)
  use precision
  use abort
  use totare_hy, only : stotot, xzntot

  use configure 
  implicit none
  logical :: stopit

#ifdef MPI_HYDRO
  raise_abort("Entropy stopper not yet implemented in MPI Case!")
#endif

  stopit =.false.

  if (maxval(stotot(1:config%qx,1,1), xzntot(1:config%qx) .lt. 1e7) .gt. 3.2) then
     stopit=.true.
  endif

end function stop_entropy


end module stopentropy
