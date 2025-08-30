!-----------------------------------------
!
! Subroutines for MPI-Parallelization
!
!-----------------------------------------
!>
!> \verbatim
!> This module provides some variables wich are needed for the 
!> MPI-version of the VERTEX-code
!>
!>  Author: K. Benkert and B. Mueller
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module mpi_vertex

  use precision

#ifndef MPI_HYDRO
  use mpi_stubs
#else
#ifdef USE_MPIF_H
  implicit none
  include "mpif.h"
#else
  use mpi
  implicit none
#endif
#endif
  
! LOCAL variables that are not in modules
  save


  public  & ! MPI-setup
          use_mpi,  myproc, hostnames, procname_lengthg, &
          pidg, ppidg, cart_comm, printit, string_array, &
          use_1neighbour_comm, use_2neighbour_comm,      &
          use_4neighbour_comm, cart_group


  public & ! "MPI-DIMENSIONS"
          nprocs, nymomsg,nymomeg, qy_sg, qy_eg, nzmomsg, nzmomeg,     &
          qz_sg, qz_eg, nymoms, nymome, nzmoms, nzmome, qy_s, qy_e,    &
          qz_s, qz_e, qy_proc, qz_proc, nymom_proc, nzmom_proc
  

  public &  ! message_tags 
          tag_sweep1, tag_sweep2, tag_sweep3,                          &
          tag_sweep31, tag_sweep32, tag_sweep4, tag_sweep41,           &
          tag_sweep42, tag_restrt1, tag_restrt2, tag_detectsh1,        &
          tag_detectsh2, tag_lag1, tag_lag2, tag_lag3,                 &
          tag_lag4, tag_lag5, tag_lag6

  public &  ! remap dimensions
          nymom_isg,nymom_ieg,nzmom_isg,nzmom_ieg, &
          nymom_osg,nymom_oeg,nzmom_osg,nzmom_oeg, &
          qy_isg,qy_ieg,qz_isg,qz_ieg, &
          qy_osg,qy_oeg,qz_osg,qz_oeg, &
          nymom_is,nymom_ie,nzmom_is,nzmom_ie, &
          qy_is,qy_ie,qz_is,qz_ie, &
          qy_os,qy_oe,qz_os,qz_oe, &
          nymom_os,nymom_oe,nzmom_os,nzmom_oe

  private



  ! parameters for domain decomposition
  !   qy_s: start
  !   qy_e: end
  !   qy_proc: local number

  INTEGER(kind=ik) :: myproc, nprocs, nymoms, nymome, nymom_proc, qy_s, &
                      qy_e, qy_proc, nzmoms, nzmome, nzmom_proc, qz_s,  &
                      qz_e, qz_proc
  INTEGER(kind=ik), DIMENSION(:), allocatable:: nymomsg, nymomeg, qy_sg, &
                                                qy_eg, nzmomsg, nzmomeg, &
                                                qz_sg, qz_eg, pidg, ppidg

  INTEGER(kind=ik) :: cart_comm, cart_group
  
  integer(kind=ik), allocatable, dimension(:) :: procname_lengthg
 
  type string_array
     character(LEN=MPI_MAX_PROCESSOR_NAME) ::  procnameg
  end type string_array

  type(string_array), allocatable :: hostnames(:)

! for remaper
  integer(kind=ik), dimension(:), allocatable:: nymom_isg,nymom_ieg,nzmom_isg,nzmom_ieg, &
                                                nymom_osg,nymom_oeg,nzmom_osg,nzmom_oeg, &
                                                qy_isg,qy_ieg,qz_isg,qz_ieg, &
                                                qy_osg,qy_oeg,qz_osg,qz_oeg

  integer(kind=ik) :: nymom_is,nymom_ie,nzmom_is,nzmom_ie, &
                      qy_is,qy_ie,qz_is,qz_ie, &
                      qy_os,qy_oe,qz_os,qz_oe, &
                      nymom_os,nymom_oe,nzmom_os,nzmom_oe

  LOGICAL ::  printit

  ! tags for mpi communication, numbers are arbitrary
  integer(kind=ik), parameter :: tag_sweep1 = 13
  integer(kind=ik), parameter :: tag_sweep2 = 14
  integer(kind=ik), parameter :: tag_sweep3 = 15
  integer(kind=ik), parameter :: tag_sweep31 = 151
  integer(kind=ik), parameter :: tag_sweep32 = 152
  integer(kind=ik), parameter :: tag_sweep4 = 16
  integer(kind=ik), parameter :: tag_sweep41 = 161
  integer(kind=ik), parameter :: tag_sweep42 = 162
      
  integer(kind=ik), parameter :: tag_restrt1 = 17
  integer(kind=ik), parameter :: tag_restrt2 = 18

  integer(kind=ik), parameter :: tag_detectsh1 = 24
  integer(kind=ik), parameter :: tag_detectsh2 = 25

  integer(kind=ik), parameter :: tag_lag1 = 63
  integer(kind=ik), parameter :: tag_lag2 = 64
  integer(kind=ik), parameter :: tag_lag3 = 65
  integer(kind=ik), parameter :: tag_lag4 = 66
  integer(kind=ik), parameter :: tag_lag5 = 67
  integer(kind=ik), parameter :: tag_lag6 = 68

  logical :: use_mpi, use_1neighbour_comm, use_2neighbour_comm, &
             use_4neighbour_comm

end module mpi_vertex

module mo_mpi

  use mpi_vertex

#ifndef MPI_HYDRO
  use mpi_stubs
#else

#ifdef USE_MPIF_H
  implicit none
  include "mpif.h"
#else
  use mpi
  implicit none
#endif

#endif

end module mo_mpi
