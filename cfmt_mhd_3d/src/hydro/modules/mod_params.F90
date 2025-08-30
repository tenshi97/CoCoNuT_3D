!> \verbatim
!>
!>     model parameters:
!>
!>
!>     parameters of the neutrino transport:
!>
!>     p_ntr   = key to switch transport on ( = 1) and off ( = 0 )
!>     p_nbk   = key to switch backreaction on ( = 1) and off ( = 0 )
!>     i_grtr  = key to switch ART-Transport on ( = 1) and off ( = 0 )
!>     ieul    = key to switch Lagrangian grid ( = 0) or
!>                             Eulerian grid ( = 1)
!>               ieul is set by Makefile-option !!!
!>     ibc = inner boundary condition
!>
!> \endverbatim
!> \author M. Rampp
!> 
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module param_rt
  use precision


  implicit none
! LOCAL variables that are not in modules

  SAVE

  integer(kind=ik) ,parameter :: ibc=-22

  integer(kind=ik)            ::  p_tau , p_raa , &
                                 p_rbb , p_sot

  integer(kind=ik)            :: isw3, &
                                     iswit,irat
!  integer(kind=ik)            :: geoen
  

  logical                     :: lisw1,lfixibc

  real(kind=rk)               :: zsf1,zsf2,zsf3, &
                                     dximax

  integer(kind=ik) :: labels, irvers_r
     

end module param_rt
