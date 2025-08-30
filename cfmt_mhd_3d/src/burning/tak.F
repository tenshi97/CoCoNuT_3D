! TAK_RATES 

c=======================================================================

      MODULE mod_tak

c=======================================================================

      IMPLICIT NONE

      CONTAINS


      subroutine tak(rhoi,yei,temi,ecne,ecmg,enne,enmg)
c=======================================================================
c Autor           : Francisco Shu Kitaura
c changed by A.Marek
c Modul           :
c Version         :
c Date            : $Date: 2004/14/05 18:37:27 $
c
c purpose: electron capture rates
c=======================================================================
      use precision

      use rateselec
      use cap

      implicit none
c     LOCAL variables that are not in modules

      real(kind=rk),    intent(in) :: rhoi, yei,temi

      real(kind=rk),    intent(out):: ecne,ecmg,enne,enmg !,gane,gamg

      real(kind=rk) :: rhoye0, tmp0,
     &     wr1 , wr0 , wy1 , wy0 , wt1, wt0
      real(kind=rk) :: takroymax,takroymin,taktemax,taktemin
      integer(kind=ik) :: irhoye0, irhoye1,
     &     itmp0, itmp1
      integer(kind=ik) i,j,nos

      call get_cap

c     rhoye0 = max( 1._rk, min( real(NROYE-1),
c     &         1._rk+real(NROYE-1)*(log10(rhoi*yei)-LROYE(1)) /
c     &         (LROYE(NROYE)-LROYE(1) )))
c     tmp0 = max( 1._rk, min( real(NTEM-1),
c     &         1._rk+real(NTEM-1)*(log10(temi)-LTEM(1)) /
c     &         (LTEM(NTEM)-LTEM(1) )))

c     irhoye0 = int( rhoye0 )    ! calculate index of rho*ye grid
c     itmp0 = int( tmp0 )        ! calculate index of tem grid

      irhoye0=1
      do while (LROYE(irhoye0) .le. log10(rhoi*yei))
         irhoye0 = irhoye0 + 1
      end do
      irhoye0 = irhoye0 - 1

      itmp0=1
      do while (LTEM(itmp0) .le. log10(temi))
         itmp0 = itmp0 + 1
      end do
      itmp0 = itmp0 - 1

      irhoye1 = irhoye0 + 1
      itmp1 = itmp0 + 1

      wr1 = (log10(rhoi*yei)-LROYE(irhoye0)) /
     &     (LROYE(irhoye1)-LROYE(irhoye0) )
      wr0 = 1.e0_rk - wr1

      wt1 = (log10(temi)-LTEM(itmp0)) /
     &     (LTEM(itmp1)-LTEM(itmp0) )
      wt0 = 1.e0_rk - wt1

      ecne = 10.e0_rk**( wr0*( wt0*( real(recne(irhoye0,itmp0)) )
     &     + wt1*( real(recne(irhoye0,itmp1)) ) )
     &     + wr1*( wt0*( real(recne(irhoye1,itmp0)) )
     &     + wt1*( real(recne(irhoye1,itmp1)) ) ) )

      ecmg = 10.e0_rk**( wr0*( wt0*( real(recmg(irhoye0,itmp0)) )
     &     + wt1*( real(recmg(irhoye0,itmp1)) ) )
     &     + wr1*( wt0*( real(recmg(irhoye1,itmp0)) )
     &     + wt1*( real(recmg(irhoye1,itmp1)) ) ) )

      enne = 10.e0_rk**( wr0*( wt0*( real(enune(irhoye0,itmp0)) )
     &     + wt1*( real(enune(irhoye0,itmp1)) ) )
     &     + wr1*( wt0*( real(enune(irhoye1,itmp0)) )
     &     + wt1*( real(enune(irhoye1,itmp1)) ) ) )

      enmg = 10.e0_rk**( wr0*( wt0*( real(enumg(irhoye0,itmp0)) )
     &     + wt1*( real(enumg(irhoye0,itmp1)) ) )
     &     + wr1*( wt0*( real(enumg(irhoye1,itmp0)) )
     &     + wt1*( real(enumg(irhoye1,itmp1)) ) ) )

c     gane = 10**( wr0*( wt0*( real(gamne(irhoye0,itmp0)) )
c     &                + wt1*( real(gamne(irhoye0,itmp1)) ) )
c     &          + wr1*( wt0*( real(gamne(irhoye1,itmp0)) )
c     &                + wt1*( real(gamne(irhoye1,itmp1)) ) ) )

c     gamg = 10**( wr0*( wt0*( real(gammg(irhoye0,itmp0)) )
c     &                + wt1*( real(gammg(irhoye0,itmp1)) ) )
c     &          + wr1*( wt0*( real(gammg(irhoye1,itmp0)) )
c     &                + wt1*( real(gammg(irhoye1,itmp1)) ) ) )

c     write(*,*) irhoye0,itmp0,recmg(irhoye0,itmp0 ),
c     &    recne(irhoye0,itmp0),enune(irhoye0,itmp0),
c     &    enumg(irhoye0,itmp0)
c     write(*,*) irhoye1,itmp1,recmg(irhoye1,itmp1),
c     &    recne(irhoye1,itmp1),enune(irhoye1,itmp1),
c     &    enumg(irhoye1,itmp1)

      end subroutine tak



      subroutine tak_vektor(rhoi,yei,temi,ecne,ecmg,enne,enmg)
c=======================================================================
c Autor           : Francisco Shu Kitaura
c Modul           :
c Version         :
c Date            : $Date: 2004/14/05 18:37:27 $
c
c purpose: electron capture rates
c=======================================================================

      use precision

      use rateselec
      use cap

      implicit none

      real(kind=rk),    intent(in) :: rhoi(:), yei(:),temi(:)

      real(kind=rk),    intent(out):: ecne(:),ecmg(:),enne(:),enmg(:) !,gane,gamg

      real(kind=rk), dimension(size(rhoi)) :: rhoye0, tmp0,
     &                  wr1 , wr0 , wy1 , wy0 , wt1, wt0
      real(kind=rk) :: takroymax,takroymin,taktemax,taktemin
      integer(kind=ik), dimension(size(rhoi)) :: irhoye0, irhoye1,
     &                  itmp0, itmp1

      integer(kind=ik) :: n_tak
      integer(kind=ik) :: i,j,nos

      n_tak=size(rhoi)

#if 0
      irhoye0=1
      do while (LROYE(irhoye0) .le. log10(rhoi*yei))
         irhoye0 = irhoye0 + 1
      end do
      irhoye0 = irhoye0 - 1

      itmp0=1
      do while (LTEM(itmp0) .le. log10(temi))
         itmp0 = itmp0 + 1
      end do
      itmp0 = itmp0 - 1
#else
      irhoye0=1
      do j=1,nroye
         do i=1,n_tak
            if (LROYE(irhoye0(i)) .le. log10(rhoi(i)*yei(i)))
     +           irhoye0(i)= irhoye0(i)+1
         end do
      end do
      irhoye0=irhoye0-1

      itmp0=1
      do j=1,nroye
         do i=1,n_tak
            if  (LTEM(itmp0(i)) .le. log10(temi(i)))
     +           itmp0(i)=itmp0(i)+1
         end do
      end do
      itmp0=itmp0-1

#endif

      irhoye1 = min(irhoye0 + 1,nroye)
      itmp1 = min(itmp0+1,ntem)
      irhoye0 = max(irhoye0,1)
      itmp0 = max(itmp0,1)

      do i=1,n_tak

         if (irhoye1(i).eq.irhoye0(i)) then
            wr1(i)=1.0_rk
            wr0(i)=0.0_rk
         else
            wr1(i) = (log10(rhoi(i)*yei(i))-LROYE(irhoye0(i))) /
     &           (LROYE(irhoye1(i))-LROYE(irhoye0(i)) )
            wr0(i) = 1._rk - wr1(i)
         end if

         if (itmp1(i).eq.itmp0(i)) then
            wt1(i)=1.0_rk
            wt0(i)=0.0_rk
         else
            wt1(i) = (log10(temi(i))-LTEM(itmp0(i))) /
     &           (LTEM(itmp1(i))-LTEM(itmp0(i)) )
            wt0(i) = 1._rk - wt1(i)
         end if

         ecne(i)=10.0_rk**
     &        (wr0(i)*(wt0(i)*(dble(recne(irhoye0(i),itmp0(i))))
     &        +wt1(i)*(dble(recne(irhoye0(i),itmp1(i)))))
     &        +wr1(i)*(wt0(i)*(dble(recne(irhoye1(i),itmp0(i))))
     &        +wt1(i)*(dble(recne(irhoye1(i),itmp1(i))))))

         ecmg(i)=10.0_rk**
     &        (wr0(i)*(wt0(i)*(dble(recmg(irhoye0(i),itmp0(i))))
     &        +wt1(i)*(dble(recmg(irhoye0(i),itmp1(i)))))
     &        +wr1(i)*(wt0(i)*(dble(recmg(irhoye1(i),itmp0(i))))
     &        +wt1(i)*(dble(recmg(irhoye1(i),itmp1(i))))))

         enne(i)=10.0_rk**
     &        (wr0(i)*(wt0(i)*(dble(enune(irhoye0(i),itmp0(i))))
     &        +wt1(i)*(dble(enune(irhoye0(i),itmp1(i)))))
     &        +wr1(i)*(wt0(i)*(dble(enune(irhoye1(i),itmp0(i))))
     &        +wt1(i)*(dble(enune(irhoye1(i),itmp1(i))))))

         enmg(i)=10.0_rk**
     &        (wr0(i)*(wt0(i)*(dble(enumg(irhoye0(i),itmp0(i))))
     &        +wt1(i)*(dble(enumg(irhoye0(i),itmp1(i)))))
     &        +wr1(i)*(wt0(i)*(dble(enumg(irhoye1(i),itmp0(i))))
     &        +wt1(i)*(dble(enumg(irhoye1(i),itmp1(i))))))

c     gane(i)=10.0_rk**(wr0(i)*(wt0(i)*(dble(gamne(irhoye0(i),itmp0(i))))
c     &                +wt1(i)*(dble(gamne(irhoye0(i),itmp1(i)))))
c     &          +wr1(i)*(wt0(i)*(dble(gamne(irhoye1(i),itmp0(i))))
c     &                +wt1(i)*(dble(gamne(irhoye1(i),itmp1(i))))))

c     gamg(i)=10.0_rk**(wr0(i)*(wt0(i)*(dble(gammg(irhoye0(i),itmp0(i))))
c     &                +wt1(i)*(dble(gammg(irhoye0(i),itmp1(i)))))
c     &          +wr1(i)*(wt0(i)*(dble(gammg(irhoye1(i),itmp0(i))))
c     &                +wt1(i)*(dble(gammg(irhoye1(i),itmp1(i))))))

c     write(*,*) irhoye0(i),itmp0(i),recmg(irhoye0(i),itmp0(i)),
c     &    recne(irhoye0(i),itmp0(i)),enune(irhoye0(i),itmp0(i)),
c     &    enumg(irhoye0(i),itmp0(i))
c     write(*,*) irhoye1(i),itmp1(i),recmg(irhoye1(i),itmp1(i)),
c     &    recne(irhoye1(i),itmp1(i)),enune(irhoye1(i),itmp1(i)),
c     &    enumg(irhoye1(i),itmp1(i))

      enddo

      end subroutine tak_vektor



      END MODULE mod_tak
