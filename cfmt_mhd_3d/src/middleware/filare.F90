!> \par provides the overloaded subroutine fileare
!>
!> \author A. Marek, MPA, Nov. 2009
!>
!> \detail This module provides the subroutine filare which are overloaded 
!> with the MPI_HYDRO or the purley OPENMP version
!>
!> Here is a list of subroutines which are overloaded depending on whether
!> purley OPENMP or hybrid MPI/OPENMP is used
!>
!>  interface-name    real-name        compiled when?
!>  filare            filare_OPENMP    if not MPI_HYDRO
!>  filare            filare_MPI_HYDRO        MPI_HYDRO
!>
!>  The purpose of filare is to fill area in low resolution regimes
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!> \endverbatim
!>
module filare_overload

  use precision

  public :: filare


  interface filare
     module procedure filare_MPI_HYDRO
  end interface


  contains


!> \par fills area in low resolution regimes - New version (MPI-parallelisation)
!>
!> \author B. Mueller and F. Hanke
!>
!> \param  iarea description of the aera by its number
!>
!> \detail
!> the purpose of filare is to fill area in low resolution regimes. This
!> routine is needed to handle different calculation areas. For a more
!> detailed description of these areas have look into the documentation
!> of setinf.f90. \n
!> If 2 areas are used (more are not supported), the first area is
!> sphericall symmetric. Therefore the calculation is done only in 1D
!> in the first angular zones and all other unresolved angular "tot"-
!> zones have to be filled with the calculated 1D grid point values. For this
!> task this values are just copied into all angular zones values of the 
!> "tot"-arrays.
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!> \endverbatim
!>
subroutine filare_MPI_HYDRO(iarea, selftime, childrentime)

  use precision
  use abort
  
  use totare_hy 
!  use arecon_hy, only: ix_are 
  use massio_hy, only: dflxtot, eflxtot, vxflxtot, vyflxtot, vzflxtot, xnflxtot 

  use mo_mpi

#ifdef HTCL
  use nutrio_hy, only: qentot,qyetot,cpotot
#else
  use nutrio_hy, only: cpotot
#endif

  use configure
  use hydro_areas_mod

#ifndef DEBUG_TIMINGS
  use cputim
#endif
  IMPLICIT NONE
!>-----------------------------------------------------------------------
!> \verbatim
!>     control indices:
!>              i*i = initial index of area in the total area
!>              i*f = final index of area in the total area
!>              io* = index offset -> full resolution = 1, etc.
!>              qn  = number of chemical species
!> \endverbatim
!>-----------------------------------------------------------------------
  real(kind=rk), intent(out)   :: selftime(2), childrentime(2)
  real(kind=rk)                :: selftime_start(2)
  integer(kind=ik), intent(in) :: iarea
  
  integer(kind=ik)             :: ixi, ixf, iox, iyi, iyf, ioy, izi, izf, ioz
  integer(kind=ik)             :: i, j, k, n, ny, nz


  selftime     = 0._rk
  childrentime = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

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

!-----------------------------------------------------------------------
!     calculate dimensions of the PPM-grid:
!-----------------------------------------------------------------------


!      write(*,'(''filare-test> ixi = '',i3,'' ixf = '',i3,
!     &          '' iox = '',i3,'' qn = '',i2)')
!     &         ixi, ixf, iox, qn
!      write(*,'(''filare-test> iyi = '',i3,'' iyf = '',i3,
!     &          '' ioy = '',i3)')
!     &         iyi, iyf, ioy
!      write(*,'(''filare-test> izi = '',i3,'' izf = '',i3,
!     &          '' ioz = '',i3)')
!     &         izi, izf, ioz

!-----------------------------------------------------------------------
!     copy a part of the total calculation area in the PPM-arrays:
!-----------------------------------------------------------------------

  ny = iyf - iyi
  ny = ny/ioy + 1

  nz = izf - izi
  nz = nz/ioz +1

  if (ioy .eq. config%qy .and. ioz .eq. config%qz) then
     
    do k=qz_s,qz_e
     do j=qy_s,qy_e
        do i=ixi,ixf 
           vextot(i,j,k)=vextot(i,qy_s,qz_s)
           veytot(i,j,k)=veytot(i,qy_s,qz_s)
           veztot(i,j,k)=veztot(i,qy_s,qz_s)
           dentot(i,j,k)=dentot(i,qy_s,qz_s)
           enetot(i,j,k)=enetot(i,qy_s,qz_s)
           pretot(i,j,k)=pretot(i,qy_s,qz_s)
           gaetot(i,j,k)=gaetot(i,qy_s,qz_s)
           gactot(i,j,k)=gactot(i,qy_s,qz_s)
           stotot(i,j,k)=stotot(i,qy_s,qz_s)
           temtot(i,j,k)=temtot(i,qy_s,qz_s)
           qgrtot(i,j,k)=qgrtot(i,qy_s,qz_s)
#ifdef HTCL
           qentot(i,j,k)=qentot(i,qy_s,qz_s)
           qyetot(i,j,k,1)=qyetot(i,qy_s,qz_s,1)
#endif /* HTCL */
        end do

        dflxtot (1,j,k,iarea)=dflxtot (1,qy_s,qz_s,iarea)
        dflxtot (2,j,k,iarea)=dflxtot (2,qy_s,qz_s,iarea)
        eflxtot (1,j,k,iarea)=eflxtot (1,qy_s,qz_s,iarea)
        eflxtot (2,j,k,iarea)=eflxtot (2,qy_s,qz_s,iarea)
        vxflxtot(1,j,k,iarea)=vxflxtot(1,qy_s,qz_s,iarea)
        vxflxtot(2,j,k,iarea)=vxflxtot(2,qy_s,qz_s,iarea)
        vyflxtot(1,j,k,iarea)=vyflxtot(1,qy_s,qz_s,iarea)
        vyflxtot(2,j,k,iarea)=vyflxtot(2,qy_s,qz_s,iarea)
        vzflxtot(1,j,k,iarea)=vzflxtot(1,qy_s,qz_s,iarea)
        vzflxtot(2,j,k,iarea)=vzflxtot(2,qy_s,qz_s,iarea)
        
        do n=1,config%qn
           xnflxtot(1,j,k,n,iarea)=xnflxtot(1,qy_s,qz_s,n,iarea)
           xnflxtot(2,j,k,n,iarea)=xnflxtot(2,qy_s,qz_s,n,iarea)
        end do
     end do

     do n=1,config%qn
        do j=qy_s,qy_e
           do i=ixi,ixf
              xnutot(i,j,k,n) = xnutot(i,qy_s,qz_s,n)
           enddo
        end do
     enddo

     do n=1,4
        do j=qy_s,qy_e
            do i=ixi,ixf
                cpotot(i,j,k,n) = cpotot(i,qy_s,qz_s,n)
            enddo
        enddo
     enddo   

    enddo
         
  else if (ioy .ne. 1 .and. ioz .ne. 1) then 
     raise_abort("filare(): ioy .ne. 1 .and. ioz .ne. 1")
!         call stopit_mpi("filare: ioy .ne. 1 .and. ioz .ne. 1")
  end if
#ifndef DEBUG_TIMINGS
  call second_v(selftime)
  selftime = selftime - selftime_start
#endif
  return
end subroutine filare_MPI_HYDRO

!=======================================================================


end module filare_overload
