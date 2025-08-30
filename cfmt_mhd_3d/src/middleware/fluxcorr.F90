!> \par provides the overloaded subroutine fluxcorr
!>
!> \author A. Marek, MPA, Nov. 2009
!>
!> \detail
!> This module provides the subroutine fluxcorr and fluxcorr_are which are 
!> overloaded with the MPI_HYDRO or the purley OPENMP version
!>
!> Here is a list of subroutines which are overloaded depending on whether
!> purley OPENMP or hybrid MPI/OPENMP is used
!>
!>  interface-name    real-name               compiled when?
!>  fluxcorr          fluxcorr_PROM           if not CFC_TRANSPORT
!>  fluxcorr_are      fluxcorr_are_MPI_HYDRO         MPI_HYDRO, no CFC
!>  fluxcorr_are      fluxcorr_are_OPENMP     if not MPI_HYDRO, no CFC
!>
!>  The purpose of fluxcorr / fluxcorr_are is to fill area in low resolution 
!>  regimes with the fluxes
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!> \endverbatim
!>
module fluxcorr_overload

use precision

public fluxcorr, fluxcorr_are

interface fluxcorr
#ifndef CFC_TRANSPORT
   module procedure fluxcorr_PROM
#endif
end interface


interface fluxcorr_are
   module procedure fluxcorr_are_MPI_HYDRO
end interface


contains

!> \par Driver for fluxcorr_are
!>
!> \author M. Rampp
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!>
!> \endverbatim
!>
subroutine fluxcorr_PROM(selftime, childrentime)
  use precision
  
!  use arecon_hy
  use totare_hy
 
  use hydro_areas_mod
  use configure
#ifndef DEBUG_TIMINGS
  use cputim
#endif
  implicit none
! LOCAL variables that are not in modules
  
  integer(kind=ik)           :: ia
  real(kind=rk), intent(out) :: selftime(2), childrentime(2)
  real(kind=rk)              :: selftime_start(2)


  selftime     = 0._rk
  childrentime = 0._rk

#ifndef DEBUG_TIMINGS
  call second_v(selftime_start)
#endif

! calculate flux correction for all areas


  if (areas%are_nu .gt. 1) then 
     do ia = 1, areas%are_nu         
        call fluxcorr_are(ia)
     enddo
  endif


#ifndef DEBUG_TIMINGS
 call second_v(selftime)
 selftime = selftime - selftime_start
#endif
 return
end subroutine fluxcorr_PROM


!> \par corrects the fluxes at the interface of two areas - Originall version (NEW version - MPI parallelisation)
!>
!> \author B. Mueller and F. Hanke
!>
!> \param iarea description of the aera by its number
!>
!> \detail
!> The task of fluxcorr is to correct the fluxes at the interface of different
!> calculation areas. This routine is needed to handle different calculation
!> areas. For a more detailed description of this areas have a look into the
!> documentation of setinf.f90. \n
!> Only two areas are handled, a spherical symmetric area in the innermost
!> region and a multi-D outer region. Therefore the flux correction has just
!> to be corrected at the right boundary of the first area. This
!> flux-correction is necessary, since a multi-D area covers a certain (maybe
!> the whole) angular domain and the 1D inner area not.
!> At first the corrected flux at the right boundary is calculated, for
!> example by:
!> \f[ \mbox{dflxr}(1,1) = \frac{\sum_{j,k}
!> \mbox{dvytot}(j) \mbox{dvztot}(k) \left\{ \mbox{dflxtot}(1,j,k,2) -
!> \mbox{dflxtot}(2,j,k,1) \right\}}{\mbox{dvy}(1)\mbox{dvz}(1)\mbox{dvx}(nx)}
!> \f]
!> Thereby is \f$ \mbox{dvy}(1) = \cos(\mbox{yzltot}(1)) - \cos(\mbox{yzrtot}
!> (qy)) \f$ and \f$ \mbox{dvz}(1) = \mbox{zzrtot}(qz) - \mbox{zzltot}(1) \f$
!> the total resolved angular domain.
!> With this angular corrected fluxes the hydrodynamic states are adjusted.
!> Therefore the flux times the timestep is added to the hydro quantity,e.g. :
!> \f[ \mbox{energy}(nx,1,1) = \frac{ \mbox{densty}(1,1) \mbox{energy}(nx,
!> 1,1) - \mbox{dt\_are}(iarea)\mbox{eflxr}(1,1))}{ \mbox{densty}(nx,1,1)}
!> \f]
!> Thereby dt_are is the timestep of this inner area.
!>
!> \verbatim
!>   SVN - Information  
!>   $Revision$
!>   $Date$
!> \endverbatim
!>
!> \verbatim
!>   SVN - Information
!>   $Revision$
!>   $Date$
!> \endverbatim
!>
SUBROUTINE fluxcorr_are_MPI_HYDRO(iarea)

!=======================================================================

  use precision
  use abort
  
  use vnew_hy 
  use totare_hy 

  use massio_hy 
  use mesh_hy

  use mo_mpi
#ifndef DEBUG_TIMINGS
  use cputim
#endif  
 
  use hydro_areas_mod

  use configure
 IMPLICIT NONE
      
  integer(kind=ik), intent(in) :: iarea
  integer(kind=ik) :: n,i,j,k,ii,iii,iif,jji,jjf,kki,kkf
  integer(kind=ik) :: ixi,ixf,iox,iyi,iyf,ioy,izi,izf,ioz, nx, ny,nz

  integer(kind=ik) :: ial,iar,nyl,nyr,ifinl,ifinr, ierr 
!  real(kind=rk), dimension(qy,qz) :: scrch_1,scrch_2

  real(kind=rk) :: rbuf(config%qn+5), rbuf_buf(config%qn+5)
 
  real(kind=rk) :: tim1(2), tim2(2)

!-- get control indices:
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


!-- control indices of adjacent areas
  ial=max(iarea-1,1)
  iar=min(iarea+1,areas%are_nu)

  nyl = (areas%ix_are(ial,5)-areas%ix_are(ial,4))/areas%ix_are(ial, 6)+1
  nyr = (areas%ix_are(iar,5)-areas%ix_are(iar,4))/areas%ix_are(iar, 6)+1

!-- check if lower adjacent area has finer y-spacing
  if (nyl .gt. ny) then
     ifinl=1
  else
     ifinl=0
  endif

! check if upper adjacent area has finer y-spacing
  if (nyr.gt.ny) then
     ifinr=1
  else
     ifinr=0
  endif


  dvx(1:nx)  = dvxtot(ixi:ixf)

  if(ioy .eq. 1) then
     yznl(1:ny)   = yzltot(iyi:iyf:ioy)
     yzn(1:ny)    = yzntot(iyi:iyf:ioy)
     yznr(1:ny)   = yzrtot(iyi:iyf:ioy)
     dvy(1:ny)    = dvytot(iyi:iyf:ioy)
  else if (ioy .eq. config%qy) then
     dvy(qy_s) = cos(yzltot(1)) - cos(yzrtot(config%qy)) !igeomy = 4 is assumed!
  else
     raise_abort("fluxcorr_are(): error: ioy .ne. 1")
!          call stopit_mpi("fluxcorr_are: error: ioy .ne. 1")
  endif

  if(ioz .eq. 1) then
     zznl(1:nz)   = zzltot(izi:izf:ioz)
     zzn(1:nz)    = zzntot(izi:izf:ioz)
     zznr(1:nz)   = zzrtot(izi:izf:ioz)
     dvz(1:nz)    = dvztot(izi:izf:ioz)
  else if (ioz .eq. config%qz) then
     dvz(qz_s) = zzrtot(config%qz) - zzltot(1)
  else
     raise_abort("fluxcorr_are(): error: ioz .ne. 1")
  endif


  if (ioy .eq. config%qy .and. ioz .eq. config%qz .and. iarea .eq. 1) then
!     spherical symmetry only in innermost region -> only right boundary matters
!     (perhaps generalize treatment at some stage, but try to avoid unnecessary
!     reduction operations)

     densty(1:nx,qy_s,qz_s) = dentot(ixi:ixf,qy_s,qz_s)
     energy(1:nx,qy_s,qz_s) = enetot(ixi:ixf,qy_s,qz_s)
     xnuc(1:nx,qy_s,qz_s,1:config%qn) = xnutot(ixi:ixf,qy_s,qz_s,1:config%qn)
     velx  (1:nx,qy_s,qz_s) = vextot(ixi:ixf,qy_s,qz_s)
     vely  (1:nx,qy_s,qz_s) = veytot(ixi:ixf,qy_s,qz_s)
     velz  (1:nx,qy_s,qz_s) = veztot(ixi:ixf,qy_s,qz_s)


! sum up flux-corrections
      dflxr(qy_s,qz_s) = 0._rk
      eflxr(qy_s,qz_s) = 0._rk
      vxflxr(qy_s,qz_s) = 0._rk
      vyflxr(qy_s,qz_s) = 0._rk
      vzflxr(qy_s,qz_s) = 0._rk
      xnflxr(qy_s,qz_s,:) = 0._rk

      do k = qz_s,qz_e
        do j = qy_s,qy_e
!            dflxl(j,k)=SUM( spread(dvytot(jji:jjf),dim=2,ncopies=kkf-kki+1)* &
!                            spread(dvztot(kki:kkf),dim=1,ncopies=jjf-jji+1)*(&
!                            -dflxtot(1+ifinl,jji:jjf,kki:kkf,iarea-ifinl) &
!                            +dflxtot(1,      jji:jjf,kki:kkf,iarea)) ) / &
!                           (dvz(k)*dvy(j)*dvx(1))

            dflxr(qy_s,qz_s)= dflxr(qy_s,qz_s) + dvytot(j)*dvztot(k)* &
              (dflxtot(2-ifinr,j,k,iarea+ifinr)-dflxtot(2,j,k,iarea))

!            eflxl(j,k)=SUM( spread(dvytot(jji:jjf),dim=2,ncopies=kkf-kki+1)* &
!                            spread(dvztot(kki:kkf),dim=1,ncopies=jjf-jji+1)*( &
!                           -eflxtot(1+ifinl,jji:jjf,kki:kkf,iarea-ifinl) &
!                            +eflxtot(1,      jji:jjf,kki:kkf,iarea)) ) / &
!                            (dvz(k)*dvy(j)*dvx(1))

            eflxr(qy_s,qz_s)= eflxr(qy_s,qz_s) + dvytot(j)*dvztot(k)* &
              (eflxtot(2-ifinr,j,k,iarea+ifinr)-eflxtot(2,j,k,iarea))

!            vxflxl(j,k)=SUM( spread(dvytot(jji:jjf),dim=2,ncopies=kkf-kki+1)* &
!                             spread(dvztot(kki:kkf),dim=1,ncopies=jjf-jji+1)*( &
!                             -vxflxtot(1+ifinl,jji:jjf,kki:kkf,iarea-ifinl) &
!                             +vxflxtot(1,      jji:jjf,kki:kkf,iarea)) ) / &
!                             (dvz(k)*dvy(j)*dvx(1))

           vxflxr(qy_s,qz_s)= vxflxr(qy_s,qz_s) + dvytot(j)*dvztot(k)* &
              (vxflxtot(2-ifinr,j,k,iarea+ifinr)-vxflxtot(2,j,k,iarea))

!            vyflxl(j,k)=SUM( spread(dvytot(jji:jjf),dim=2,ncopies=kkf-kki+1)* &
!                        spread(dvztot(kki:kkf),dim=1,ncopies=jjf-jji+1)*( &
!                        -vyflxtot(1+ifinl,jji:jjf,kki:kkf,iarea-ifinl) &
!                        +vyflxtot(1,      jji:jjf,kki:kkf,iarea)) ) / &
!                        (dvz(k)*dvy(j)*dvx(1))

           vyflxr(qy_s,qz_s)= vyflxr(qy_s,qz_s) + dvytot(j)*dvztot(k)* &
              (vyflxtot(2-ifinr,j,k,iarea+ifinr)-vyflxtot(2,j,k,iarea))

!            vzflxl(j,k)=SUM(spread(dvytot(jji:jjf),dim=2,ncopies=kkf-kki+1)* &
!                            spread(dvztot(kki:kkf),dim=1,ncopies=jjf-jji+1)*(&
!                            -vzflxtot(1+ifinl,jji:jjf,kki:kkf,iarea-ifinl) &
!                            +vzflxtot(1,      jji:jjf,kki:kkf,iarea)) ) / &
!                            (dvz(k)*dvy(j)*dvx(1))

           vzflxr(qy_s,qz_s)= vzflxr(qy_s,qz_s) + dvytot(j)*dvztot(k)* &
              (vzflxtot(2-ifinr,j,k,iarea+ifinr)-vzflxtot(2,j,k,iarea))

            do n=1,config%qn
!               xnflxl(j,k,n)=SUM(spread(dvytot(jji:jjf),dim=2,ncopies=kkf-kki+1)* &
!                                 spread(dvztot(kki:kkf),dim=1,ncopies=jjf-jji+1)*( &
!                                 -xnflxtot(1+ifinl,jji:jjf,kki:kkf,n,iarea-ifinl) &
!                                 +xnflxtot(1,      jji:jjf,kki:kkf,n,iarea)) ) / &
!                                (dvz(k)*dvy(j)*dvx(1))

              xnflxr(qy_s,qz_s,n)= xnflxr(qy_s,qz_s,n) + dvytot(j)*dvztot(k)* &
                (xnflxtot(2-ifinr,j,k,n,iarea+ifinr)-xnflxtot(2,j,k,n,iarea))

            enddo

        enddo
      enddo

!renormalize
      dflxr(qy_s,qz_s) = dflxr(qy_s,qz_s) / (dvz(qz_s)*dvy(qy_s)*dvx(nx))
      eflxr(qy_s,qz_s) = eflxr(qy_s,qz_s) / (dvz(qz_s)*dvy(qy_s)*dvx(nx))
      vxflxr(qy_s,qz_s) = vxflxr(qy_s,qz_s) / (dvz(qz_s)*dvy(qy_s)*dvx(nx))
      vyflxr(qy_s,qz_s) = vyflxr(qy_s,qz_s) / (dvz(qz_s)*dvy(qy_s)*dvx(nx))
      vzflxr(qy_s,qz_s) = vzflxr(qy_s,qz_s) / (dvz(qz_s)*dvy(qy_s)*dvx(nx))

      do n=1,config%qn
        xnflxr(qy_s,qz_s,n) = xnflxr(qy_s,qz_s,n) / (dvz(qz_s)*dvy(qy_s)*dvx(nx))
      enddo



      if (use_mpi) then
#ifndef DEBUG_TIMINGS
         call second_v(tim1)
#endif
         !     MPI-Allreduce (Summe) fuer dflxl(j,k),dflxr(j,k),...xnflxl(j,k,1:config%qn),xnflxr(1:config%qn)
         rbuf(1) = dflxr(qy_s,qz_s)
         rbuf(2) = eflxr(qy_s,qz_s)
         rbuf(3) = vxflxr(qy_s,qz_s)
         rbuf(4) = vyflxr(qy_s,qz_s)
         rbuf(5) = vzflxr(qy_s,qz_s)
         rbuf(6:5+config%qn) = xnflxr(qy_s,qz_s,1:config%qn)
 
     

         call MPI_allreduce(rbuf, rbuf_buf, config%qn+5,MPI_DOUBLE_PRECISION, &
                            MPI_SUM, MPI_COMM_WORLD, ierr)
 
         dflxr(qy_s,qz_s)  = rbuf_buf(1) 
         eflxr(qy_s,qz_s)  = rbuf_buf(2)
         vxflxr(qy_s,qz_s) = rbuf_buf(3)
         vyflxr(qy_s,qz_s) = rbuf_buf(4)
         vzflxr(qy_s,qz_s) = rbuf_buf(5)
         xnflxr(qy_s,qz_s,1:config%qn) = rbuf_buf(6:5+config%qn) 
#ifndef DEBUG_TIMINGS
         call second_v(tim2)

         timer%transp_comm = timer%transp_comm + (tim2-tim1)
#endif
      endif ! use_mpi

! correct hydrodynamic state


!      scrch_1(qy_s,1)=densty(1 ,qy_s,1)
!     scrch_2(qy_s,qz_s)=densty(nx,qy_s,qz_s)

!      densty(1, qy_s,1)=densty(1 ,qy_s,1)
!     &     - areas%dt_are(iarea)*dflxl(qy_s,1)
     densty(nx,qy_s,qz_s)=densty(nx,qy_s,qz_s) - areas%dt_are(iarea)*dflxr(qy_s,qz_s)

!      energy(1, qy_s,1)=(scrch_1(qy_s,1)*energy(1 ,qy_s,1)
!     &     - areas%dt_are(iarea)*eflxl(qy_s,1))/densty(1, qy_s,1)
     energy(nx,qy_s,qz_s)=(densty(nx,qy_s,qz_s)*energy(nx,qy_s,qz_s) &
                        - areas%dt_are(iarea)*eflxr(qy_s,qz_s))/densty(nx,qy_s,qz_s)

!      velx(1, qy_s,1)=(scrch_1(qy_s,1)*velx(1 ,qy_s,1)
!     &     - areas%dt_are(iarea)*vxflxl(qy_s,1))/densty(1, qy_s,1)
     velx(nx,qy_s,qz_s)=(densty(nx,qy_s,qz_s)*velx(nx,qy_s,qz_s) &
                       - areas%dt_are(iarea)*vxflxr(qy_s,qz_s))/densty(nx,qy_s,qz_s)

!      vely(1, qy_s,1)=(scrch_1(qy_s,1)*vely(1 ,qy_s,1)
!     &     - areas%dt_are(iarea)*vyflxl(qy_s,1))/densty(1, qy_s,1)
     vely(nx,qy_s,qz_s)=(densty(nx,qy_s,qz_s)*vely(nx,qy_s,qz_s) &
                       - areas%dt_are(iarea)*vyflxr(qy_s,qz_s))/densty(nx,qy_s,qz_s)

!      velz(1, qy_s,1)=(scrch_1(qy_s,1)*velz(1 ,qy_s,1)
!     &     - areas%dt_are(iarea)*vzflxl(qy_s,1))/densty(1, qy_s,1)
     velz(nx,qy_s,qz_s)=(densty(nx,qy_s,qz_s)*velz(nx,qy_s,qz_s) &
                       - areas%dt_are(iarea)*vzflxr(qy_s,qz_s))/densty(nx,qy_s,qz_s)

     do n=1,config%qn
!         xnuc(1, qy_s,1,n)=(scrch_1(qy_s,1)*xnuc(1 ,qy_s,1,n)
!     &        - areas%dt_are(iarea)*xnflxl(qy_s,1,n))/densty(1, qy_s,1)
          xnuc(nx,qy_s,qz_s,n)=(densty(nx,qy_s,qz_s)*xnuc(nx,qy_s,qz_s,n) &
                       - areas%dt_are(iarea)*xnflxr(qy_s,qz_s,n))/densty(nx,qy_s,qz_s)
     enddo


! write back arrays
     do k=qz_s,qz_e
        do j=qy_s,qy_e
          dentot(ixi:ixf,j,k) = densty(1:nx,qy_s,qz_s)
          enetot(ixi:ixf,j,k) = energy(1:nx,qy_s,qz_s)
          xnutot(ixi:ixf,j,k,1:config%qn) =  xnuc(1:nx,qy_s,qz_s,1:config%qn)
          vextot(ixi:ixf,j,k) = velx  (1:nx,qy_s,qz_s)
          veytot(ixi:ixf,j,k) = vely  (1:nx,qy_s,qz_s)
          veztot(ixi:ixf,j,k) = velz  (1:nx,qy_s,qz_s)
        end do
     end do

  else if (ioy.ne.1) then
!          call stopit_mpi("fluxcorr_are: error: ioy .ne. 1")

     raise_abort("fluxcorr_are(): error: ioy .ne. 1")
  end if

  return

END SUBROUTINE fluxcorr_are_MPI_HYDRO

!=======================================================================



end module fluxcorr_overload
