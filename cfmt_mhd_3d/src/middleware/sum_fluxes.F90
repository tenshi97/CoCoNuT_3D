!>
!> \verbatim
!> This module provides the subroutine sumcq and sumcqflux which are needed 
!> to compute all fluxes leaving the grid 
!>
!> Both functions are declared public to be useable in other procedures
!>
!>
!>  Author: A. Marek, MPA, Nov. 2009
!>
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!> \endverbatim
!>
module sum_fluxes

public sumcq, sumcqflux

contains
!>
!> This function computes the sum of the neutrino sourceterms
!> \param dencq sourceterm (e.g. qentot, qyetot)
!> \param dti   timestep
!> 
!>
!>  Author: M. Rampp
!>
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!>
 function sumcq(dencq,dti)

   use precision
   use abort

   use gfloat_hy
   use totare_hy

   use mo_mpi

   implicit none
! LOCAL variables that are not in modules
   
   real(kind=rk) ,intent(in) :: dencq(:,:,:),dti
   real(kind=rk) :: sumcq,dsurf
   real(kind=rk) :: dvr(SIZE(dencq,DIM=1))
   integer(kind=ik) :: nnx,nny,nnz,j,k

   nnx=SIZE(dencq,DIM=1)
   nny=SIZE(dencq,DIM=2)
   nnz=SIZE(dencq,DIM=3)
   
   if (use_mpi) then
      ! MPI Kontrollabfrage if nny .ne. qy_proc then error
      if (nny .ne. qy_proc .or. nnz .ne. qz_proc) then
         raise_abort("sumcq(): nny .ne. qy_proc")
         !         call stopit_mpi ("sumcq: nny .ne. qy_proc")
      endif
   endif

   dvr = (xzrtot(1:nnx)**3 - xzltot(1:nnx)**3)/ 3._rk
   sumcq = 0._rk

   do k  = 1, qz_proc
      do j = 1, qy_proc
         dsurf = vlfrac * dvztot(k+qz_s-1) * dvytot(j+qy_s-1)
         sumcq = sumcq + dsurf  * SUM( dvr * dencq(1:nnx,j,k) ) * dti
      enddo
   enddo
   
 end function sumcq
    

!>
!> This function computes the sum of the neutrino fluxes over
!> the grid boundaries
!>
!> \param flucq neutrino flux
!> \param dti   timestep
!> 
!>
!>  Author: M. Rampp
!>
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>   
!>
 function sumcqflux(flucq,dti)

   use precision
   use abort
   
   use gfloat_hy
   use totare_hy

   use mo_mpi

   implicit none
! LOCAL variables that are not in modules

   real(kind=rk) ,intent(in) :: flucq(:,:),dti
   real(kind=rk) :: sumcqflux,dsurf
   integer(kind=ik) :: nny,nnz,j,k


   nny=SIZE(flucq,DIM=1)
   nnz=SIZE(flucq,DIM=2)


   if(use_mpi) then
      ! MPI Kontrollabfrage if nny .ne. qy_proc then error
      if (nny .ne. qy_proc .or. nnz .ne. qz_proc) then
         raise_abort("sumcqflux(): nny .ne. qy_proc")
         !         call stopit_mpi ("sumcqflux: nny .ne. qy_proc")
      endif
   endif
      
   sumcqflux = 0._rk

   do k  = 1, qz_proc
      do j = 1, qy_proc
         dsurf = vlfrac * dvztot(k+qz_s-1) * dvytot(j+qy_s-1)
         sumcqflux = sumcqflux + dsurf * flucq(j,k) * dti
      enddo
   enddo


 end function sumcqflux

end module sum_fluxes


