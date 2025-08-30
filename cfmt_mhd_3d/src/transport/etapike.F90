!===============================================================

MODULE mod_degeneracy_parameter

!===============================================================

  IMPLICIT NONE
  
  PUBLIC etapike
  
CONTAINS

! --------------------------------------------------------------
  
  FUNCTION etapike(dn0)

! --------------------------------------------------------------
! Autor             : Markus Rampp
! modul             : $RCSfile: interactions.F,v $
! version           : $Revision: 1.12 $
! date of creation  : $Date: 2003/06/11 08:31:18 $
!
! Purpose: Approximately calculate eta from given (normalized)
!          Number density dn0
!          dn0 = fac * dn
!            where dn is the number density and
!                  fac=3*pi**2*(hquer*c)**3/(2*m*c**2*kT)**1.5
!
! Reference: Th. Hecht, Diplomarbeit, TUM 1998 Eq.(2.22)
!            ===> BEWARE OF MISPRINT in Hechts Eq.(2.22) <===
    USE precision
    
    USE phycon, ONLY : pc_pi
    
    IMPLICIT NONE

    ! LOCAL variables that are not in modules
    REAL(KIND=rk) etapike,dn83,dn23,dn43,dn0
    
    REAL(KIND=rk) ,parameter :: pi=pc_pi
    REAL(KIND=rk) ,parameter :: r=0.9235_rk,s=0.7665_rk,u=1.2970_rk
    
    dn23=dn0**(2._rk/3._rk)
    dn43=dn23*dn23         ! = dn0**(4./3.)
    dn83=dn43*dn43         ! = dn0**(8./3.)
    
    etapike=dn23*(1._rk-pi**2/(12._rk*dn43))*(dn83/(u+dn83))+ &
         log(4._rk*dn0/(3._rk*sqrt(pi)))/(1._rk+r*dn0**2+s*dn0**4)
    
    
  end function etapike
  
  

END MODULE mod_degeneracy_parameter

!===============================================================
