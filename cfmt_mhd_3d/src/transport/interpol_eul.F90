module interpol_eul

implicit none

contains

#ifndef NOTRA
#if defined(EULGRID) && !defined(PROGRAM_remap)

subroutine intpint

!***********************************************************************
!
! Author     : Markus Rampp (MPA) 
! Modul      : $RCSfile: interpol_eul.F,v $
! Version    : $Revision$
! Date       : $Date$
!
! PURPOSE: eddfak_vec updates the Intensities by a 
!          lagrangian derivative of I
!
!           for an EULERIAN GRID
!              this is simulated by radial interpolation 
!              of the Intensities given at the old timelevel onto
!              the "backwardly moved" radial grid
!
!           in order to do the interpolation at fixed angle-cosines mue
!              the interpolation has to be a bilinear one 
!
! NOTE: everything can be done outside the ME-BTE iteration-loop 
!       rip,rim must contain Intensities given at the OLD TIMELEVEL
!
!***********************************************************************

  use precision
  
  use radial_grid_rt
  use tanray_grids_rt
  use backquants_rt
  use radfield_rt !, ripold => rip, rimold => rim
  use time_rt
  use abort

  use configure

! LOCAL variables that are not in modules

  implicit none
  integer(kind=ik) :: i,is,ie,j,k,krq,klq,kr,kl
  integer(kind=ik) :: iip1q,iiq,iip1,ii,niqp1,irq
  integer(kind=ik) :: niq,irqp1,nip1,ni,irp1,np,ir
  real(kind=rk)    ::  ralt(-1:config%imaxp+1),rqalt(0:config%imaxp+1)
  real(kind=rk)    :: bq1,bq2,b3,b4,b1,b2,bq3,bq4,aq1
  real(kind=rk)    :: aq2,a2,a1
  integer(kind=ik) :: klpos(-1:config%imaxp+1,config%cmin:config%imaxp+1),  &
                      krpos(-1:config%imaxp+1,config%cmin:config%imaxp+1),  &
                      klposq(0:config%imaxp+1,config%cmin:config%imaxp+1),  &
                      krposq(0:config%imaxp+1,config%cmin:config%imaxp+1),  &
                      ipos(-1:config%imaxp+1),iposq(0:config%imaxp+1)
!-->      k[l,r]pos,amue    must have the same dimensions (for sort_vec)
!-->      k[l,r]pos,amueq   must have the same dimensions (for sort_vec)

  real(kind=rk)   ::  pa(-1:config%imaxp+1),paq(0:config%imaxp+1),   &
                      pb1(-1:config%imaxp+1,config%cmin:config%imaxp+1),    &
                      pb1q(0:config%imaxp+1,config%cmin:config%imaxp+1),    &
                      pb2(-1:config%imaxp+1,config%cmin:config%imaxp+1),    &
                      pb2q(0:config%imaxp+1,config%cmin:config%imaxp+1)

  integer(kind=ik) :: imp2

  logical          :: lstop

  imp2 = config%imaxp+1+2
! ---------

!CDIR NEXPAND(sort_vec)

  lstop=.FALSE.

#ifdef CFC_TRANSPORT2
  do i=-1,config%imaxp+1
     ralt(i) =r(i)-be(i)*a_phi2(i)/dxi1ex(i)
  enddo
  do i= 0,config%imaxp+1
     rqalt(i)=rq(i)-beq(i)*a_phi2q(i)/dxi1exq(i)
  enddo
#else
  do i=-1,config%imaxp+1
     ralt(i) =r(i)-be(i)*ex(i)/dxi1ex(i)
  enddo
  do i= 0,config%imaxp+1
     rqalt(i)=rq(i)-beq(i)*exq(i)/dxi1exq(i)
  enddo
#endif

!      write(*,*) 'WARNING: TEST !!!'
!      do i= 0,config%imaxp+1
!         ralt(i)=r(i)+1.5*(r(i)-r(i-1))
!      enddo
!      do i= 0,imaxp
!         rqalt(i)=rq(i)+1.7*(rq(i+1)-r(i))
!      enddo
!      ralt(config%imaxp+1)=2.0*r(config%imaxp+1)
!      ralt(0)=0.2*r(0)
!      rqalt(config%imaxp+1)=1.5*rq(config%imaxp+1)


! -- find positions
  call sort_vec(r(-1:),config%imaxp+1+2,1,ralt(-1:),ipos(-1:),config%imaxp+1+2,1)
  call sort_vec(rq(0:),config%imaxp+1+1,1,rqalt(0:),iposq(0:),config%imaxp+1+1,1)



!      i=-1
!         np     = i-config%cmin+1
!         ir=max(ipos(i)-2,-1)
!         irp1=min(ir+1,config%imaxp+1)
!         ni=ir-config%cmin+1
!         nip1=irp1-config%cmin+1
!         call sort_vec(amue(ir,  config%cmin),ni,  imp2,
!     &                 amue(i,config%cmin),klpos(i,config%cmin),np ,imp2)
!         call sort_vec(amue(irp1,config%cmin),nip1,imp2,
!     &        amue(i,config%cmin),krpos(i,config%cmin),np ,imp2)

  do i=0,config%imaxp+1
     np     = i-config%cmin+1
! extrapolate
     if (ralt(i).lt.r(-1)) then
        lstop=.TRUE.
     elseif (ralt(i).ge.r(config%imaxp+1)) then
        ipos(i)  =config%imaxp+1-1
     else
! interpolate
        ipos(i)   =ipos(i)-2
     endif

! extrapolate
     if (rqalt(i).lt.rq(0)) then
        iposq(i) =0
     elseif (rqalt(i).ge.rq(config%imaxp+1)) then
        iposq(i) =config%imaxp+1-1
! interpolate
     else
        iposq(i) =iposq(i)-1
     endif

     ir=ipos(i)
     irp1=ir+1
     ni=ir-config%cmin+1
     nip1=irp1-config%cmin+1
     call sort_vec(amue(ir,  config%cmin:),ni,  1, amue(i,config%cmin:),klpos(i,config%cmin:),np,1)
     call sort_vec(amue(irp1,config%cmin:),nip1,1, amue(i,config%cmin:),krpos(i,config%cmin:),np,1)

     irq=iposq(i)
     irqp1=irq+1
     niq=irq-config%cmin
     niqp1=irqp1-config%cmin+1
     call sort_vec(amueq(irq,  config%cmin:),niq,1, amueq(i,config%cmin:),klposq(i,config%cmin:),np-1,1)
     call sort_vec(amueq(irqp1,config%cmin:),niq,1, amueq(i,config%cmin:),krposq(i,config%cmin:),np-1,1)

  enddo

  if (lstop)  raise_abort("intpint(): this should not happen for ibc=-22")

! -- compute interpolation-weights

!      i=-1
!         ii  =max(ipos(i)-2,-1)
!         iip1=ii+1
!         if (ralt(i).ge.r(config%imaxp+1)) then   !extrapolate
!            pa(i) = zero
!         elseif (ralt(i).lt.r(-1)) then   !extrapolate
!            pa(i) = zero
!         else                             !interpolate
!            pa(i) =(ralt(i)-r(ii))/(r(iip1)-r(ii))
!         endif
!         do k=config%cmin,i-1
!            kl=klpos(i,k)+config%cmin-1
!            kr=krpos(i,k)+config%cmin-1
!            pb1(i,k) =(amue(i,k)-amue(ii  ,kl))/
!     &           (amue(ii,kl+1)-amue(ii,kl))
!            pb2(i,k) =(amue(i,k)-amue(iip1,kr))/
!     &           (amue(iip1,kr+1)-amue(iip1,kr))
!         enddo
  do i=0,config%imaxp+1
     ii=ipos(i)
     iip1=ii+1
     
     iiq=iposq(i)
     iip1q=iiq+1
     
     pa(i) =(ralt(i)-r(ii))/(r(iip1)-r(ii))
     paq(i)=(rqalt(i)-rq(iiq))/(rq(iip1q)-rq(iiq))
        


     do k=config%cmin,i-1
        kl=klpos(i,k)+config%cmin-1
        kr=krpos(i,k)+config%cmin-1
        klq=klposq(i,k)+config%cmin-1
        krq=krposq(i,k)+config%cmin-1
        
        if (ii .lt. 0) then
           pb1(i,k) = 0.0_rk
        else
           pb1(i,k) =(amue(i,k)-amue(ii  ,kl))/ (amue(ii,kl+1)-amue(ii,kl))
        endif
        pb2(i,k) =(amue(i,k)-amue(iip1,kr))/ (amue(iip1,kr+1)-amue(iip1,kr))

        pb1q(i,k)=(amueq(i,k)-amueq(iiq  ,klq))/ (amueq(iiq,klq+1)-amueq(iiq,klq))
        pb2q(i,k)=(amueq(i,k)-amueq(iip1q,krq))/ (amueq(iip1q,krq+1)-amueq(iip1q,krq))
     enddo
  enddo

!      do i=-1,config%imaxp+1
!         ir=ipos(i)-2
!         write(*,'(2I4,3e13.5)') i,ir,r(ir),ralt(i),r(i+1)
!      enddo
!c      write(*,*) '-----------'
!c      do i=0,201,202
!c         write(*,'(1I4,1e18.10)') i,be(i)
!c         ir=ipos(i)-2
!c         irp1=min(ir+1,config%imaxp+1)
!c         do k=config%cmin,i
!c            kl=klpos(i,k)+config%cmin-1
!c            kr=krpos(i,k)+config%cmin-1
!c            write(*,'(2I4,3e18.10)') k,kl,amue(ir,kl),amue(i,k),
!c     &           amue(ir,kl+1)
!c            write(*,'(2I4,3e18.10)') k,kr,amue(irp1,kr),amue(i,k),
!c     &           amue(irp1,kr+1)
!c            write(*,*) ' '
!c         enddo
!c      enddo
!c      stop
!
!      open(63,form='unformatted')
!      write(63) pa,paq
!      write(63) pb1,pb2,pb1q,pb2q
!      close(63)


! -- do the linear interpolation
  do is=1,config%isma
     do ie=1,config%iemax
!         i=-1
!            ii  =max(ipos(i)-2,-1)
!            iip1=min(ii+1,config%imaxp+1)
!
!            a2 =pa(i)
!            a1 =one-a2
!            do k=config%cmin,i-1
!               kl=klpos(i,k)+config%cmin-1
!               kr=krpos(i,k)+config%cmin-1
!               
!               b2=pb1(i,k)
!               b1=one-b2
!               b4=pb2(i,k)
!               b3=one-b4
!              
!               rimoldi(i,k,ie,is)=a1*(b1*rimold(ii,  kl,ie,is)+
!     &                                b2*rimold(ii,  kl+1,ie,is))+
!     &                            a2*(b3*rimold(iip1,kr,ie,is)+
!     &                                b4*rimold(iip1,kr+1,ie,is))
!            enddo

        do i=0,config%imaxp+1
           ii=ipos(i)
           iip1=ii+1
           iiq  =iposq(i)
           iip1q =iiq+1

           a2 =pa(i)
           a1 =1.0_rk-a2
           aq2=paq(i)
           aq1=1.0_rk-aq2
           do k=config%cmin,i-1
              kl=klpos(i,k)+config%cmin-1
              kr=krpos(i,k)+config%cmin-1
              klq=klposq(i,k)+config%cmin-1
              krq=krposq(i,k)+config%cmin-1

              b2=pb1(i,k)
              b1=1.0_rk-b2
              b4=pb2(i,k)
              b3=1.0_rk-b4

              bq2=pb1q(i,k)
              bq1=1.0_rk-bq2
              bq4=pb2q(i,k)
              bq3=1.0_rk-bq4
              
              rimoldi(i,k,ie,is)=a1*(b1*rim(ii,  kl,ie,is)+         &
                                     b2*rim(ii,  kl+1,ie,is))+      &
                                 a2*(b3*rim(iip1,kr,ie,is)+         &
                                     b4*rim(iip1,kr+1,ie,is))
     
              ripoldi(i,k,ie,is)=aq1*(bq1*rip(iiq,  klq,ie,is)   +  &
                                      bq2*rip(iiq,  klq+1,ie,is))+  &
                                 aq2*(bq3*rip(iip1q,krq  ,ie,is) +  &
                                      bq4*rip(iip1q,krq+1,ie,is))

! (only) extrapolation could have produced negative values
              ripoldi(i,k,ie,is)=max(ripoldi(i,k,ie,is) , 0.0_rk)
           enddo
        enddo
        j=-1
        rimoldi(j,:,ie,is) = 1.e+100_rk
        do j=0,config%imaxp+1
           rimoldi(j,j,ie,is) = 1.e+100_rk
           ripoldi(j,j,ie,is) = 1.e+100_rk
        enddo
     enddo
  enddo

!      open(64,form='unformatted')
!      write(64) r(-1:config%imaxp+1),ralt(-1:config%imaxp+1),amue,
!     &     rimoldi(:,:,2,1)*1e-30,
!     &     rimold(:,:,2,1)*1e-30
!      write(64) rq(0:config%imaxp+1),rqalt(0:config%imaxp+1),amueq,
!     &     ripoldi(:,:,2,1)*1e-30,
!     &     ripold(:,:,2,1)*1e-30
!      close(64)
!      stop

end subroutine intpint
#endif /* EULGRID && !PROGRAM_remap */


!----------------------------------------------------------------------
#if defined(LAGGRID) && !defined(PROGRAM_remap)
subroutine intpint

!***********************************************************************
!
! Autor             : Markus Rampp (MPA) 
! Modul             : $RCSfile: interpol_eul.F,v $
! Version           : $Revision$
! Date of creation  : $Date$
! Time of creation  : 
! Date of extraction: %D%
! Time of extraction: %T%
!
! PURPOSE: eddfak_vec updates the Intensities by a 
!          lagrangian derivative of I
!
!           for a LAGRANGIAN GRID,
!             the angle-cosines mue of the tangent-ray grid
!             change during the timestep
!
!           in order to keep them fixed during the timestep
!              the Intensities given at the old timelevel (and therefore
!              defined onto the OLD mue-grid are interpolated in angle 
!              onto the NEW mue-grid 
!
! NOTE: everything can be done outside the ME-BTE iteration-loop 
!       rip,rim must contain Intensities given at the OLD TIMELEVEL
!***********************************************************************
  use precision
  
  use tanray_grids_rt
  use backquants_rt
  use radfield_rt, ripold => rip, rimold => rim
! LOCAL variables that are not in modules

  implicit none

  integer(kind=ik) ::    ktar(-1:config%imaxp+1,config%cmin:config%imaxp+1), &
                         ktarq(0:config%imaxp+1,config%cmin:config%imaxp+1)
!-->      ktar,amue    must have the same dimensions (for sort_vec)
!-->      ktarq,amueq  must have the same dimensions (for sort_vec)

  integer, parameter  :: imp2=config%imaxp+1+2
  logical             :: omp_in_parallel
!CDIR NEXPAND (linintp,monintp,sort_vec)
!---     do the sorting once for all energies      ---
!               timedependent !
!
!---                      
  do jk=0,config%imaxp+1
     np     = jk-config%cmin+1
!DIR$ INLINE
     call sort_vec(amueqold(jk,config%cmin:),np-1,1,amueq(jk,config%cmin:), ktarq(jk,config%cmin:),np-1,1)
  enddo

  do jk=-1,config%imaxp+1
     np     = jk-config%cmin+1
!DIR$ INLINE
     call sort_vec(amueold(jk,config%cmin:),np,1,amue(jk,config%cmin:), ktar(jk,config%cmin:),np,1)
  enddo



!---           evaluate interpolating polynomial
!                     to be done whenever grid has moved
! 
!---
#if defined(OPENMP_TRANSPORT) && defined(OPEN_MP_1D)
!OMP$ PARALLEL DO  IF (.not.OMP_IN_PARALLEL())
!OMP$&     SHARED(amueq,amueqold,amue,amueold,ripold,ripoldi,
!OMP$&            rimold,rimoldi,ktar,ktarq,imp2,
!OMP$&            config)
!OMP$&     PRIVATE(ie,is,jk,np,kk,kkk,aderq,ader)
#endif
  do ie=1,config%iemax
     do is=1,config%isma

        do jk=0,config%imaxp+1
           np=jk-config%cmin+1

!            do kk=config%cmin,jk-1
!               kkk=ktarq(jk,kk)+config%cmin-1
!               aderq = ( ripold(jk,kkk+1,ie,is)-
!     &              ripold(jk,kkk,ie,is) ) /
!     &              ( amueqold(jk,kkk+1)-amueqold(jk,kkk) )
!
!               ripoldi(jk,kk,ie,is)=ripold(jk,kkk,ie,is)+
!     &              aderq * (amueq(jk,kk)-amueqold(jk,kkk))   
!            enddo

!CDIR$ INLINE
!c -- monintp requires 3times more CPU than linintp
!            call monintp(amueqold(jk,config%cmin),ripold(jk,config%cmin,ie,is),
!     &               np,imp2-1, ktarq(jk,config%cmin), amueq(jk,config%cmin),
!     &           ripoldi(jk,config%cmin,ie,is),np-1,imp2-1,2,0)
!DIR$ INLINE
           call linintp(amueqold(jk,config%cmin),ripold(jk,config%cmin,ie,is),    &
                        np,imp2-1, ktarq(jk,config%cmin), amueq(jk,config%cmin),  &
                        ripoldi(jk,config%cmin,ie,is),np-1,imp2-1)

    
        enddo

        do jk=-1,config%imaxp+1
           np=jk-config%cmin+1

!            do kk=config%cmin,jk-1
!               kkk=ktar(jk,kk)+config%cmin-1
!               ader = ( rimold(jk,kkk+1,ie,is)-
!     &           rimold(jk,kkk,ie,is) ) /
!     &           ( amueold (jk,kkk+1)-amueold (jk,kkk) )
!               rimoldi(jk,kk,ie,is)=rimold(jk,kkk,ie,is)+
!     &              ader *(amue(jk,kk)- amueold(jk,kkk)) 
!            enddo
  
!CDIR$ INLINE
!c -- monintp requires 3times more CPU than linintp
!            call monintp(amueold(jk,config%cmin),rimold(jk,config%cmin,ie,is),
!     &           np,imp2, ktar(jk,config%cmin),amue(jk,config%cmin),
!     &           rimoldi(jk,config%cmin,ie,is),np, imp2,2,1)   
!DIR$ INLINE
           call linintp(amueold(jk,config%cmin),rimold(jk,config%cmin,ie,is),     &
                        np,imp2, ktar(jk,config%cmin),amue(jk,config%cmin),       &
                        rimoldi(jk,config%cmin,ie,is),np, imp2) 

        enddo

        do jk=-1,config%imaxp+1
           rimoldi(jk,jk,ie,is) = rimold(jk,jk,ie,is)
           if (jk .ge. 0) ripoldi(jk,jk,ie,is) = ripold(jk,jk,ie,is)
        enddo

     enddo
  enddo

  return
end subroutine intpint
#endif /* defined(LAGGRID) && !defined(PROGRAM_remap) */

#endif /* NOTRA */
! -------------------------------------------------------------------
SUBROUTINE sort_vec(xx,nx,incx,rr,ixx0,nr,incr)

!-----------------------------------------
! Author(s)         : Max Ruffert, Markus Rampp
! Modul             : $RCSfile: interpol_eul.F,v $
! Version           : $Revision$
! Date of creation  : $Date$
! Time of creation  : 
! Date of extraction: %D%
! Time of extraction: %T%
!
! Purpose :    search for the locations ix0 of a Vector r in a Vector x
!
!                   x has to be in  monotonically (strict) 
!                     as- or descending order    
!
! Input: x:      "basic grid"; section of xx       
!        nx:     Dimension of x
!        incx:   Increment between Elements of xx to be searched
! 
!        r:      Positions to be searched for; section of rr 
!        nr:     Dimension of r (and ix0)
!        incr:   Increment between Elements of rr (and ixx0)
!
! (Output): ix0  |  x(ix0(k)) .le. r(k) .lt. x(ix0(k)+1) if x ascending
!           ix0  |  x(ix0(k)) .ge. r(k) .gt. x(ix0(k)+1) if x descending
!
!         ix0(k) = 0    if r(k) .lt.(.gt.) x(1)  .AND. x as(des)cending
!         ix0(k) = nx   if r(k) .gt.(.lt.) x(nx) .AND. x as(des)cending
!
!  Output: the original Section of ixx0
!         
!-----------------------------------------
  use precision

! LOCAL variables that are not in modules
!  IMPLICIT REAL(KIND=RK) (A-H,O-Z),INTEGER(KIND=IK) (I-K,M-N), &
!           LOGICAL (L)

  implicit none


  integer(kind=ik), intent(in)    :: nx,nr,incx,incr

  real(kind=rk)                   :: x(nx),r(nr)
  real(kind=rk), intent(in)       :: xx(:),rr(:)
  integer(kind=ik), intent(out)   :: ixx0(:)
  Integer(kind=ik), DIMENSION(nr) :: ix1,ittm, ix0

  real(kind=rk)                   :: isect(nr),isecsum
  integer(kind=ik)                :: i,ii
  real(kind=rk)                   :: xm
  logical                         :: ldescnd

! -- set up the 1D-sections
  do i=1,nx
     ii=(i-1)*incx+1
     x(i)=xx(ii)
  enddo
  ldescnd = x(nx).lt.x(1)

  do i=1,nr
     ii=(i-1)*incr+1
     r(i)=rr(ii)
  enddo

! -- prepare for iteration
  do i=1,nr
     ix0(i) = 1
     ix1(i) = nx
  enddo

  do i = 1,nr

#ifdef IBM_OPT
     
     isect(i)=FSEL((-1._rk*(ix1(i)-(ix0(i)+1))),0._rk,1._rk)
     

#else /* IBM_OPT */
     if ( ix1(i) .eq. ix0(i)+1 ) then
        isect(i) = 0._rk
     else
        isect(i) = 1._rk
     endif

#endif /* IBM_OPT */
     if ((r(i) .gt.  x(1)).eqv.ldescnd) then
        isect(i) = 0._rk
        ix0(i)=0
     endif
     if ((r(i) .lt. x(nx)).eqv.ldescnd) then
        isect(i) = 0._rk
        ix0(i)=nx
     endif
  enddo
      
! --- bisection iteration

100 continue

  do i = 1,nr
     if ( isect(i) .eq. 1._rk ) then

        ittm(i) = (ix1(i) + ix0(i)) / 2

        xm = x(ittm(i))
        if (ldescnd) then


!#ifdef IBM_OPT
!           ix0(i)=int(FSEL((-1._rk*(r(i)-xm)),real(ittm(i),kind=rk),real(ix0(i),kind=rk) ))
!           ix1(i)=int(FSEL((xm-r(i)),real(ix1(i),kind=rk),real(ittm(i),kind=rk)))
!#else /* IBM_OPT */
           
           if ( r(i) .le. xm ) then 
              ix0(i) = ittm(i)
           else
              ix1(i) = ittm(i)
           endif
!#endif /* IBM_OPT */
        else

!#ifdef IBM_OPT
!           ix0(i)=int(FSEL((r(i)-xm),real(ittm(i),kind=rk),real(ix0(i),kind=rk)))

!           ix1(i)=int(FSEL((r(i)-xm),real(ix1(i),kind=rk),real(ittm(i),kind=rk)))
!#else /* IBM_OPT */
           if ( r(i) .ge. xm ) then 
              ix0(i) = ittm(i)
           else
              ix1(i) = ittm(i)
           endif
!#endif /* IBM_OPT */
        endif
     endif
  enddo
  


  do i = 1, nr
!#ifdef IBM_OPT
!     
!     ! ( ix1(i) .gt. ix0(i)+1 ) translates to ix1(i) .ge. ix0(i)+2
!     isect(i)=FSEL(real((ix1(i)-(ix0(i)+2)),kind=rk),1._rk,0._rk)
!
!#else /* IBM_OPT */

     if ( ix1(i) .gt. ix0(i)+1 ) then
        isect(i) = 1._rk
     else
        isect(i) = 0._rk
     endif

!#endif /* IBM_OPT */
  enddo

  isecsum = 0._rk

  do i = 1, nr
     isecsum = isecsum + isect(i)
  enddo

  if ( isecsum .ne. 0._rk ) goto 100

!-- copy back
  do i=1,nr
     ii=(i-1)*incr+1
     ixx0(ii)=ix0(i)
  enddo

  return
end SUBROUTINE sort_vec

! -------------------------------------------------------------------
SUBROUTINE linintp(xol,yol,no,inco,ipoxx,xne,yne,nn,incn)
!
!     Linear interpolation of an array-section
!     
!     Extrapolation of values at outer(inner)most gridpoints if
!       necessary
!
!
! LOCAL varibales that are not in modules

  use precision

!  IMPLICIT REAL(KIND=RK) (A-H,O-Z),INTEGER(KIND=IK) (I-K,M-N), &
!           LOGICAL (L)
  implicit none

  integer(kind=ik) :: i,ii
  integer(kind=ik), intent(in) :: no,nn,incn,inco
  real(kind=rk), intent(in) ::  xol(:),yol(:),xne(:)
  real(kind=rk), intent(inout) :: yne(:)
  integer(kind=ik),intent(in) :: ipoxx(:)
  

  real(kind=rk) :: xo(no),yo(no),xn(nn),yn(nn)
  integer(kind=ik) :: ipox(nn)
  real(kind=rk) ::  s(no)

       
! -- set up the 1D-sections
  do i=1,no
     ii=(i-1)*inco+1
     xo(i)=xol(ii)
     yo(i)=yol(ii)
  enddo
  do i=1,nn
     ii=(i-1)*incn+1
     xn(i)=xne(ii)
     ipox(i)=ipoxx(ii)
  enddo
!--determine slopes
  do i=1, no-1
     s(i) = ( yo(i+1) - yo(i) ) / (xo(i+1) -  xo(i))
  enddo
  s(no) = 0._rk

!-- calculate y(x):
  do i=1, nn
     ii=ipox(i)
     if (ii.gt.no) then
        yn(i) = yo(no)
     elseif (ii.le.0) then
        yn(i) = yo(1)
     else
        yn(i) = yo(ii) + s(ii)*(xn(i) - xo(ii))
     endif
  enddo

!-- copy back
  do i=1,nn
     ii=(i-1)*incn+1
     yne(ii)=yn(i)
  enddo
  
  RETURN
END SUBROUTINE linintp

!c -------------------------------------------------------------------
SUBROUTINE monintp(xol,yol,no,inco,ipoxx,xne,yne,nn,incn, iibc,iobc)
!
!     Second order polynomial interpolation due to 
!     M. Steffen, A&A, 239, p.443 (1990)
!     adapted from taz-Version
!
!
!     Input: iibc,iobc  inner and outer boundary condition for slopes
!                      =0: slope=0 
!                      =1: slope of next zone
!                      =2: assumes a functional dependence of
!                           f=A+B*x^2; (A;B) given by the points n-1,n
!                DEFAULT: as in Steffen
!

!
! LOCAL varibales that are not in modules
  use precision


!  IMPLICIT REAL(KIND=RK) (A-H,O-Z),INTEGER(KIND=IK) (I-K,M-N), &
!           LOGICAL (L)

  implicit none
  integer(kind=ik), intent(in) :: no,nn,ipoxx(:),inco,incn,iibc,iobc
  real(kind=rk), intent(in) :: xol(:),xne(:),yol(:)
  real(kind=rk), intent(out)::  yne(:)

  real(kind=rk) ::  xo(no),yo(no),xn(nn),yn(nn)
   integer(kind=ik) :: ipox(nn)
  real(kind=rk) ::  s(no),h(no),yl(no)
  integer(kind=ik) :: i,ii
  real(kind=rk) :: dx,b,a,p,div


       
! -- set up the 1D-sections
  do i=1,no
     ii=(i-1)*inco+1
     xo(i)=xol(ii)
     yo(i)=yol(ii)
  enddo
  do i=1,nn
     ii=(i-1)*incn+1
     xn(i)=xne(ii)
     ipox(i)=ipoxx(ii)
  enddo
  !
  DO i=1, no-1
     h(i) = xo(i+1) -  xo(i)
     s(i) = ( yo(i+1) - yo(i) ) / h(i)
  ENDDO
!-- boundary values  
  SELECT CASE(iibc)
  CASE(0)
     yl(1)=0.0_rk
  CASE(1)
     yl(1)=s(1)
  CASE(2)
     yl(1)=2.0_rk*xo(1)*(yo(2)-yo(1))/(xo(2)**2-xo(1)**2)
  CASE DEFAULT
     div    = h(1)/(h(1) + h(2))
     p      = s(1)* (1._rk+div) - s(2)*div
     yl(1) = ( sign(1._rk,s(1)) + sign(1._rk,p) ) *( min( abs(s(1)), 0.5_rk*abs(p) ) )
  END SELECT
  SELECT CASE(iobc)
  CASE(0)
     yl(no)=0._rk
  CASE(1)
     yl(no)=s(no-1)
  CASE(2)
     yl(no)=2._rk*xo(no)*(yo(no)-yo(no-1))/(xo(no)**2- xo(no-1)**2)
  CASE DEFAULT
     div    = h(no-1)/(h(no-1) + h(no-2))
     p      = s(no-1)* (1._rk+div) - s(no-2)*div
     yl(no) = ( sign(1._rk,s(no-1)) + sign(1._rk,p) ) * ( min( abs(s(no-1)), 0.5_rk*abs(p) ) )
  END SELECT
!-- determine slopes
  DO i=2, no-1
     p       = ( s(i-1)*h(i) + s(i)*h(i-1) ) / ( h(i-1) + h(i) )
     yl(i)   = ( sign(1._rk,s(i-1)) + sign(1._rk,s(i)) ) * ( min( abs(s(i-1)), abs(s(i)), 0.5_rk*abs(p) ) )
  ENDDO

!-- calculate y(x):
  DO i=1, nn
     ii=ipox(i)
     if (ii.ge.no) then
        yn(i) = yo(no)
     elseif(ii.le.0) then
        yn(i) = yo(1)
     else
        a       = ( yl(ii) + yl(ii+1) - 2._rk*s(ii) )/ (h(ii)*h(ii))
        b       = ( 3._rk*s(ii) - 2._rk*yl(ii) - yl(ii+1) ) / h(ii)
        dx      = xn(i) - xo(ii)
        yn(i) = ((( a*dx ) + b )*dx + yl(ii) )*dx + yo(ii)
     endif
  ENDDO

!-- copy back
  do i=1,nn
     ii=(i-1)*incn+1
     yne(ii)=yn(i)
  enddo

  RETURN
END SUBROUTINE monintp

end module interpol_eul
