!> \verbatim
!> This module provides the subroutine qave_are which is overloaded 
!> with the purley OpenMP or the hybdird MPI/OpenMP version
!>
!> Here is a list of subroutines which are overloaded depending on which
!> parallelisation model is used
!>
!>  interface-name    real-name          compiled when?
!>  qave_are          qave_are_MPI_HYDRO   MPI_HYDRO
!>  qave_are          qave_are_OPENMO    if not MPI_HYDRO
!>
!>  The purpose of qave_are is to compute an average of the different
!>  areas
!>
!>  Author: A. Marek, MPA, Nov. 2009
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module qave_overload

  use precision

  public qave_are


  interface qave_are

     module procedure qave_are_MPI_HYDRO

  end interface

contains



!=======================================================================
!> \verbatim
!>=======================================================================
!>      MPI version
!>=======================================================================
!>
!>
!>     task:  perform an average
!>
!> Author : B. Mueller
!>=======================================================================
!> \endverbatim
!>
!> \param  isw  iarea  identify the different areas by their number 
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
SUBROUTINE qave_are_MPI_HYDRO(iarea)

  use precision
  use abort


  use vnuw_hy ! all 
  use totare_hy, only: dvytot, dvztot
!  use arecon_hy, only: ix_are 
  use nutrio_hy, only: qentot, qyetot, qmotot, qmytot 

  use mo_mpi
#ifndef DEBUG_TIMINGS
  use cputim
#endif

  use configure
  use hydro_areas_mod
  IMPLICIT NONE
  
  real(kind=rk) :: qtmp (3*config%qx), qtmp_buf(3*config%qx)
  
  integer(kind=ik), intent(in) :: iarea
  integer(kind=ik) :: i,j,k,m,l
  integer(kind=ik) :: ixi,ixf,iox,iyi,iyf,ioy,izi,izf,ioz,nx,ny,nz
  integer(kind=ik) :: ierr
      
  real(kind=rk) :: domega_total
  real(kind=rk) :: tim1(2), tim2(2)

!-----------------------------------------------------------------------
!     get control indices:
!-----------------------------------------------------------------------

  ixi    = areas%ix_are(iarea, 1)
  ixf    = areas%ix_are(iarea, 2)
  iox    = areas%ix_are(iarea, 3)
  iyi    = areas%ix_are(iarea, 4)
  iyf    = areas%ix_are(iarea, 5)
  ioy    = areas%ix_are(iarea, 6)
  izi    = areas%ix_are(iarea, 7)
  izf    = areas%ix_are(iarea, 8)
  ioz    = areas%ix_are(iarea, 9)
  
  nz = izf - izi
  ny = iyf - iyi
  nx = ixf - ixi
  
  nz = nz/ioz + 1
  ny = ny/ioy + 1
  nx = nx/iox + 1
  
  if (ioy .eq. config%qy .and. ioz .eq. config%qz) then
! i.e. the inner spherical symmetric area is treated here !!!

!-- due to nonlocal transport effects qyetot,qentot,qmotot  vary
!    within the unresolved hydro areas; therefore averages are 
!    appropriate
     domega_total=1.0_rk/(SUM(dvytot(1:config%qy))*SUM(dvztot(1:config%qz)))
     k=qz_s
     j=qy_s
     qen(ixi:ixf,j,k)=0.0_rk
     qye(ixi:ixf,j,k,:)=0.0_rk
     qmo(ixi:ixf,j,k)=0.0_rk
     qmy(ixi:ixf,j,k)=0.0_rk
     do m=qz_s,qz_e
        do l=qy_s,qy_e
           do i = 1,nx
              qen(i,j,k)  =qen(i,j,k)+ qentot(i,l,m)*dvytot(l)*dvztot(m)
              qye(i,j,k,1)=qye(i,j,k,1)+ qyetot(i,l,m,1)*dvytot(l)*dvztot(m)
!             qye(i,j,k,2)=qye(i,j,k,2)+! qyetot(i,l,k,2)*dvytot(l)*dvztot(m)
!             qye(i,j,k,3)=qye(i,j,k,3)+! qyetot(i,l,k,3)*dvytot(l)*dvztot(m)
!             qye(i,j,k,4)=qye(i,j,k,4)+! qyetot(i,l,k,4)*dvytot(l)*dvztot(m)
!             qye(i,j,k,5)=qye(i,j,k,5)+! qyetot(i,l,k,5)*dvytot(l)*dvztot(m)

              qmo(i,j,k)  =qmo(i,j,k)+qmotot(i,l,m)*dvytot(l)*dvztot(m)
              qmy(i,j,k)  =qmy(i,j,k)+ qmytot(i,l,m)*dvytot(l)*dvztot(m)
           end do
        end do
     end do

     qen(ixi:ixf,j,k)   = qen(ixi:ixf,j,k)  *domega_total
     qye(ixi:ixf,j,k,1) = qye(ixi:ixf,j,k,1)*domega_total
!            qye(ixi:ixf,j,k,2) = qye(ixi:ixf,j,k,2)*domega_total
!            qye(ixi:ixf,j,k,3) = qye(ixi:ixf,j,k,3)*domega_total
!            qye(ixi:ixf,j,k,4) = qye(ixi:ixf,j,k,4)*domega_total
!            qye(ixi:ixf,j,k,5) = qye(ixi:ixf,j,k,5)*domega_total


     qmo(ixi:ixf,j,k)   = qmo(ixi:ixf,j,k)  *domega_total
     qmy(ixi:ixf,j,k)   = 0.0_rk



     if (use_mpi) then
        !     MPI-Allreduce (Summe) fuer qen(ixi:ixf,j,k),qye(ixi:ixf,j,k,1),qmo(ixi:ixf,j,k)
        !     qye(:,:,:,2) is only a diagnostic quantity -> no averaging necessary


#ifndef DEBUG_TIMINGS
        call second_v(tim1)
#endif
        
        qtmp(     1:  nx)=qen(1:nx,qy_s,qz_s)
        qtmp(  nx+1:2*nx)=qye(1:nx,qy_s,qz_s,1)
        qtmp(2*nx+1:3*nx)=qmo(1:nx,qy_s,qz_s)


        call MPI_allreduce(qtmp(1:3*nx), qtmp_buf(1:3*nx), 3*nx,MPI_DOUBLE_PRECISION, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)

        qtmp(:)=qtmp_buf(:)
        
        qen(1:nx,qy_s,qz_s)  =qtmp_buf(     1:  nx)
        qye(1:nx,qy_s,qz_s,1)=qtmp_buf(  nx+1:2*nx)

        qmo(1:nx,qy_s,qz_s)  =qtmp_buf(2*nx+1:3*nx)
#ifndef DEBUG_TIMINGS
        call second_v(tim2)
        timer%transp_comm = timer%transp_comm + (tim2-tim1)
#endif
     endif ! use_mpi


     do k = qz_s,qz_e
        do j = qy_s,qy_e
            qentot(ixi:ixf,j,k)   = qen(1:nx,qy_s,qz_s)
            qyetot(ixi:ixf,j,k,1) = qye(1:nx,qy_s,qz_s,1)
!            qyetot(ixi:ixf,j,1,2) = qye(1:nx,qy_s,1,2)
!            qyetot(ixi:ixf,j,1,3) = qye(1:nx,qy_s,1,3) 
!            qyetot(ixi:ixf,j,1,4) = qye(1:nx,qy_s,1,4) 
!            qyetot(ixi:ixf,j,1,5) = qye(1:nx,qy_s,1,5)       
            qmotot(ixi:ixf,j,k)   = qmo(1:nx,qy_s,qz_s)
            qmytot(ixi:ixf,j,k)   = 0.0_rk
        end do
     end do

  else if (ioy.eq.1 .and. ioz.eq.1) then

     do k = qz_s,qz_e
        do j = qy_s,qy_e
          qen(1:nx,j,k)   = qentot(ixi:ixf,j,k)
          qye(1:nx,j,k,1) = qyetot(ixi:ixf,j,k,1)
          qye(1:nx,j,k,2) = qyetot(ixi:ixf,j,k,2)
          qye(1:nx,j,k,3) = qyetot(ixi:ixf,j,k,3)
          qye(1:nx,j,k,4) = qyetot(ixi:ixf,j,k,4)
          qye(1:nx,j,k,5) = qyetot(ixi:ixf,j,k,5)
          qmo(1:nx,j,k)   = qmotot(ixi:ixf,j,k)
          qmy(1:nx,j,k)   = qmytot(ixi:ixf,j,k)
        end do
     end do

  else if (ioy.ne.1) then
     raise_abort("qave_are(): error: ioy .ne. 1")
!         call stopit_mpi("qave_are: error: ioy .ne. 1")
  end if


END SUBROUTINE qave_are_MPI_HYDRO

!=======================================================================



end module qave_overload
