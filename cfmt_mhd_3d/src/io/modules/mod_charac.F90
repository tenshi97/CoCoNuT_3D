!>
!> \verbatim
!> This module provides the variables which contain information
!> about the modell name , name of the restart files and so on... 
!> with the MPI_HYDRO or the purley OPENMP version
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
module charac
  use precision

  implicit none
! LOCAL variables that are not in modules

  SAVE

  !      character(8)  label
  character(200) rstfil, outfil, rstfil_ra, outfil_ra
 ! character(80) rst_mode

  character(200) rstfil_cfc,outfil_cfc
 ! character(200) basenm
 ! character(200) calculation
  character(200)  filpos
#if defined(MPI_HYDRO) && !(defined(PROGRAM_remap))
 ! character(7)  suffix
#else
 ! character(2)  suffix
#endif

#if defined(NEW_OUTPUT_LABELS) || defined(PROGRAM_remap)
  integer(kind=ik) :: file_number
 ! character (len=8):: suffix_nr
#endif

  character (len=8):: myproc_character

! modinp,
#if defined(PROGRAM_remap)
  character(200)    :: gridfile,param_fil,gridfile_hydro
  character(200)    ::  basenm_i,basenm_o
#endif

end module charac
