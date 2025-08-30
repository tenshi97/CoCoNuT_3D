module fcnc_calc

  private
  public :: initialize_fcnc_calc

  contains


!> \verbatim
!> This subroutine initialises the fcnc neutrino reactions
!>
!>  Author: A. Marek, MPA, March 2009
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
subroutine initialize_fcnc_calc
  use precision

  use neutrinotypes

  implicit none
! LOCAL variables that are not in modules

  real(kind=rk) :: GTMR_delta

  open(15,file='GTMR_delta')
  read(15,*) GTMR_delta
  read(15,*) fcnc_fac
  close(15)

  fcnc_delta(1)=GTMR_delta
  fcnc_delta(2)=GTMR_delta
  fcnc_delta(3)=0._rk

  if (fcnc_fac.eq.2.) then
     fcnc_delta(1:2) = 0.0_rk
     fcnc_delta(3)   = GTMR_delta
     fcnc_fac = 0.0_rk
  endif

  write (*,*) 'GTMR, mode ',fcnc_fac,', delta = ',fcnc_delta(1:3)
  
  if (isma.eq.1) fcnc_fac = 0._rk

end subroutine initialize_fcnc_calc

end module fcnc_calc
