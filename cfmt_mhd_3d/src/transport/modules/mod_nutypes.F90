!>
!> \par This module provides the definition of the neutrino flavours used in Vertex
!>
!>  \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision: 731 $
!>   $Date: 2010-02-03 18:27:12 +0100 (Wed, 03 Feb 2010) $
!>   
!> \endverbatim
!>
module neutrinotypes      

  use precision

  use abort
  implicit none
! LOCAL variables that are not in modules

  SAVE

  character (13) ,parameter, dimension(3) :: &
          neustr=(/'electron     ','anti-electron', &
                   'mu/tau       '/)

  integer(kind=ik), parameter, dimension(3) :: neu=(/1,2,3/)
  integer(kind=ik), parameter, dimension(3) :: neu_a=(/2,1,3/)

  real(kind=rk), parameter, dimension(3) :: neu_qw=(/1.0_rk,1.0_rk,4.0_rk /)

  real(kind=rk), dimension(3) :: siglep
#ifdef FCNC_CALC
  real(kind=rk) :: fcnc_delta(3),fcnc_fac
#endif
!   use the array "neu" to attribute a specific neutrino-type to
!    its logical number in the code:
!      e.g. if only electron-antineutrinos are to be transported,
!    (isma=1 and) neu(1)=2 

! Table:
! 
!  1 = electron
!  2 = anti-electron
!  3 = muon
!  4 = anti-muon
!  5 = tauon
!  6 = anti-tauon

contains 


  subroutine nutypes(inu)
    use precision
    use print_stdout_mod
    use mpi_vertex

    implicit none
! LOCAL variables that are not in modules

    integer(kind=ik) ,intent(in) :: inu
    
    real(kind=rk), parameter, dimension(3) :: sigl=(/+1.0_rk,-1.0_rk,0.0_rk/)

    integer(kind=ik) :: is,nu
      
    if (inu .gt. 3) then
       raise_abort("nutypes(): do not use more than three neutrino types")
    endif


    call printit_taskX(0,"nutypes>")
    call printit_taskX(0," ")
    call printit_taskX(0,"--> transport of ....")

    do is=1,inu

       nu=neu(is)
       siglep(is)=sigl(nu)
       if (myproc .eq. 0) then
          write(*,'(I1,A2,A13,1x,A20,1f4.1)') is,': ', &
              neustr(neu(is)),'neutrinos; siglep = ',siglep(is)
       endif
    enddo
    call printit_taskX(0," ")
    call printit_taskX(0,"end nutypes>")

  end subroutine nutypes


end module neutrinotypes
