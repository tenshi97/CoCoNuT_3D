      MODULE mod_get_rate

c=======================================================================
c
c     This modules contains scalar (get_rate) and vector (get_rate1,
c     get_rate2, get_rate3, get_rate4, get_rate5) routine to
c     compute the rates for the Kitaura burning network. Note that
c     calling get_rate from inner loops may be computationally
c     inefficient because of bad inlining properties.
c
c ----------------------------------------------------------------------

      IMPLICIT NONE

      CONTAINS


c=======================================================================

      subroutine get_rate (tmpr, dens, sigb, reacr, reacn, mode)

c=======================================================================
c Autor           : Shu Kitaura
c                   changed by A.Marek
c Modul           :
c Version         :
c Date            : $Date: 2004/14/05 18:37:27 $
c
c purpose: calculates the rate
c
c t913 = t9^(1./3.)
c t953 = t9^(5./3.)
c
c rate = exp(rate(i,1) + rate(i,2)/t9 + rate(i,3)/t913 + rate(i,4)*t913 +
c      rate(i,5)*t9 + rate(i,6)*t953 + rate(i,7)*ln(t9))
c
c     (for the reaction i)
c
c corresponds to ln(2)/t12 for decays, NA<ab> for two-body reactions and
c NA**2<abc> for three-body reactions.
c=======================================================================
      use precision
      use abort

      use burnrate

      implicit none
c LOCAL variables that are not in modules

      integer(kind=ik), intent(in) :: mode, reacr,reacn
      real(kind=rk),    intent(in) :: tmpr, dens
      real(kind=rk),    intent(out):: sigb

      real(kind=rk) ::tem9,temp9,tem9l,tem9lthd,tem9lfthd,t09(1:7),dummy

      integer(kind=ik) :: j

      real(kind=rk), parameter ::  thrd = 1.0_rk/3.0_rk

      call rates

         write(0,*) "tmpr = ", tmpr

         tem9      = tmpr*1.e-9_rk
         tem9l     = log(tem9)
         tem9lthd  = thrd*tem9l
         write(0,*) "tem9lthd = ", tem9lthd
         tem9lfthd = 5.0_rk*tem9lthd
         t09(1)  = 1.0_rk
         t09(2)  = 1.0_rk/tem9
         t09(3)  = exp( -tem9lthd )
         write(0,*) "t09(3) = ", t09(3)
         t09(4)  = 1.0_rk/t09(3)
         t09(5)  = tem9
         t09(6)  = exp( tem9lfthd )
         t09(7)  = tem9l
         temp9  = tem9



            sigb = 0.0_rk


      select case(mode)

c     rates due to reactions involving a single nuclei
      case(1)

         do j = 1,7
            sigb = sigb+t09(j)*rate(reacr,j)
         end do

         sigb =  exp ( sigb )

       case(2)

         do j = 1,7
            sigb = sigb+t09(j)*rate(reacn,j)
         end do

        sigb =   exp ( sigb )

        dummy = sigb

          do j = 1,7
            sigb = sigb+t09(j)*rate(reacr,j)
         end do

         sigb = sigb + dummy

c     rates due to reactions involving two nuclei
      case(3)

         do j = 1,7
            sigb = sigb+t09(j)*rate(reacr,j)
         end do

         sigb =  dens * exp ( sigb )


      case(4)

         do j = 1,7
            sigb = sigb+t09(j)*rate(reacn,j)
         end do

        sigb =  dens *  exp ( sigb )

        dummy = sigb

          do j = 1,7
            sigb = sigb+t09(j)*rate(reacr,j)
         end do

         sigb = dens * sigb + dummy

      case default
         raise_abort("get_rate(): wrong mode")
      end select

      end subroutine get_rate

c=======================================================================






c=======================================================================

      subroutine get_rate1 (tmpr, dens, sigb, reacr, reacn)

c=======================================================================
c Autor           : Shu Kitaura
c                   changed by A.Marek, vectorized by B.Mueller
c Modul           :
c Version         :
c Date            : $Date: 2004/14/05 18:37:27 $
c
c purpose: calculates the rate
c
c t913 = t9^(1./3.)
c t953 = t9^(5./3.)
c
c rate = exp(rate(i,1) + rate(i,2)/t9 + rate(i,3)/t913 + rate(i,4)*t913 +
c      rate(i,5)*t9 + rate(i,6)*t953 + rate(i,7)*ln(t9))
c
c     (for the reaction i)
c
c corresponds to ln(2)/t12 for decays, NA<ab> for two-body reactions and
c NA**2<abc> for three-body reactions.
c=======================================================================

      USE precision
      USE burnrate

      IMPLICIT NONE

      integer(kind=ik), intent(in) :: reacr,reacn
      real(kind=rk),    intent(in) :: tmpr(:), dens(:)
      real(kind=rk),    intent(out):: sigb(:)

      real(kind=rk), dimension(size(dens)) ::
     &     tem9,temp9,tem9l,tem9lthd,tem9lfthd,dummy
      real(kind=rk), dimension(size(dens),7) :: t09

      integer(kind=ik) :: j

      real(kind=rk), parameter ::  thrd = 1.0_rk/3.0_rk

      tem9      = tmpr*1.d-9
      tem9l     = log(tem9)
      tem9lthd  = thrd*tem9l
      tem9lfthd = 5.d0*tem9lthd
      t09(:,1)  = 1.d0
      t09(:,2)  = 1.d0/tem9
      t09(:,3)  = exp( -tem9lthd )
      t09(:,4)  = 1.d0/t09(:,3)
      t09(:,5)  = tem9
      t09(:,6)  = exp( tem9lfthd )
      t09(:,7)  = tem9l
      temp9  = tem9


      sigb = 0.0d0


c     rates due to reactions involving a single nuclei

      do j = 1,7
         sigb = sigb+t09(:,j)*rate(reacr,j)
      end do

      sigb =  exp ( sigb )

      return
      end subroutine get_rate1

c=======================================================================



c=======================================================================

      subroutine get_rate2 (tmpr, dens, sigb, reacr, reacn)

c=======================================================================
c Autor           : Shu Kitaura
c                   changed by A.Marek
c Modul           :
c Version         :
c Date            : $Date: 2004/14/05 18:37:27 $
c
c purpose: calculates the rate
c
c t913 = t9^(1./3.)
c t953 = t9^(5./3.)
c
c rate = exp(rate(i,1) + rate(i,2)/t9 + rate(i,3)/t913 + rate(i,4)*t913 +
c      rate(i,5)*t9 + rate(i,6)*t953 + rate(i,7)*ln(t9))
c
c     (for the reaction i)
c
c corresponds to ln(2)/t12 for decays, NA<ab> for two-body reactions and
c NA**2<abc> for three-body reactions.
c=======================================================================

      USE precision
      USE burnrate

      IMPLICIT NONE

      integer(kind=ik), intent(in) :: reacr,reacn
      real(kind=rk),    intent(in) :: tmpr(:), dens(:)
      real(kind=rk),    intent(out):: sigb(:)

      real(kind=rk), dimension(size(dens)) ::
     &     tem9,temp9,tem9l,tem9lthd,tem9lfthd,dummy
      real(kind=rk), dimension(size(dens),7) :: t09

      integer(kind=ik) :: j

      real(kind=rk), parameter ::  thrd = 1.0_rk/3.0_rk

      tem9      = tmpr*1.d-9
      tem9l     = log(tem9)
      tem9lthd  = thrd*tem9l
      tem9lfthd = 5.d0*tem9lthd
      t09(:,1)  = 1.d0
      t09(:,2)  = 1.d0/tem9
      t09(:,3)  = exp( -tem9lthd )
      t09(:,4)  = 1.d0/t09(:,3)
      t09(:,5)  = tem9
      t09(:,6)  = exp( tem9lfthd )
      t09(:,7)  = tem9l
      temp9  = tem9


      sigb = 0.0d0


      do j = 1,7
         sigb=sigb+t09(:,j)*rate(reacn,j)
      end do

      sigb=exp(sigb)

      dummy=sigb
c     neu
      sigb=0.0d0

      do j = 1,7
         sigb=sigb+t09(:,j)*rate(reacr,j)
      end do

      sigb=exp(sigb)
      sigb=sigb+dummy

      return
      end subroutine get_rate2

c=======================================================================



c=======================================================================

      subroutine get_rate3 (tmpr, dens, sigb, reacr, reacn)

c=======================================================================
c Autor           : Shu Kitaura
c                   changed by A.Marek
c Modul           :
c Version         :
c Date            : $Date: 2004/14/05 18:37:27 $
c
c purpose: calculates the rate
c
c t913 = t9^(1./3.)
c t953 = t9^(5./3.)
c
c rate = exp(rate(i,1) + rate(i,2)/t9 + rate(i,3)/t913 + rate(i,4)*t913 +
c      rate(i,5)*t9 + rate(i,6)*t953 + rate(i,7)*ln(t9))
c
c     (for the reaction i)
c
c corresponds to ln(2)/t12 for decays, NA<ab> for two-body reactions and
c NA**2<abc> for three-body reactions.
c=======================================================================

      USE precision
      USE burnrate

      IMPLICIT NONE

      integer(kind=ik), intent(in) :: reacr,reacn
      real(kind=rk),    intent(in) :: tmpr(:), dens(:)
      real(kind=rk),    intent(out):: sigb(:)

      real(kind=rk), dimension(size(dens)) ::
     &     tem9,temp9,tem9l,tem9lthd,tem9lfthd,dummy
      real(kind=rk), dimension(size(dens),7) :: t09

      integer(kind=ik) :: j

      real(kind=rk), parameter ::  thrd = 1.0_rk/3.0_rk

      tem9      = tmpr*1.d-9
      tem9l     = log(tem9)
      tem9lthd  = thrd*tem9l
      tem9lfthd = 5.d0*tem9lthd
      t09(:,1)  = 1.d0
      t09(:,2)  = 1.d0/tem9
      t09(:,3)  = exp( -tem9lthd )
      t09(:,4)  = 1.d0/t09(:,3)
      t09(:,5)  = tem9
      t09(:,6)  = exp( tem9lfthd )
      t09(:,7)  = tem9l
      temp9  = tem9


      sigb = 0.0d0

      do j = 1,7
         sigb = sigb+t09(:,j)*rate(reacr,j)
      end do

      sigb =  dens * exp ( sigb )

      return
      end subroutine get_rate3

c=======================================================================



c=======================================================================

      subroutine get_rate4 (tmpr, dens, sigb, reacr, reacn)

c=======================================================================
c Autor           : Shu Kitaura
c                   changed by A.Marek
c Modul           :
c Version         :
c Date            : $Date: 2004/14/05 18:37:27 $
c
c purpose: calculates the rate
c
c t913 = t9^(1./3.)
c t953 = t9^(5./3.)
c
c rate = exp(rate(i,1) + rate(i,2)/t9 + rate(i,3)/t913 + rate(i,4)*t913 +
c      rate(i,5)*t9 + rate(i,6)*t953 + rate(i,7)*ln(t9))
c
c     (for the reaction i)
c
c corresponds to ln(2)/t12 for decays, NA<ab> for two-body reactions and
c NA**2<abc> for three-body reactions.
c=======================================================================

      USE precision
      USE burnrate

      IMPLICIT NONE

      integer(kind=ik), intent(in) :: reacr,reacn
      real(kind=rk),    intent(in) :: tmpr(:), dens(:)
      real(kind=rk),    intent(out):: sigb(:)

      real(kind=rk), dimension(size(dens)) ::
     &     tem9,temp9,tem9l,tem9lthd,tem9lfthd,dummy
      real(kind=rk), dimension(size(dens),7) :: t09

      integer(kind=ik) :: j

      real(kind=rk), parameter ::  thrd = 1.0_rk/3.0_rk

      tem9      = tmpr*1.d-9
      tem9l     = log(tem9)
      tem9lthd  = thrd*tem9l
      tem9lfthd = 5.d0*tem9lthd
      t09(:,1)  = 1.d0
      t09(:,2)  = 1.d0/tem9
      t09(:,3)  = exp( -tem9lthd )
      t09(:,4)  = 1.d0/t09(:,3)
      t09(:,5)  = tem9
      t09(:,6)  = exp( tem9lfthd )
      t09(:,7)  = tem9l
      temp9  = tem9


      sigb = 0.0d0

      do j = 1,7
         sigb = sigb+t09(:,j)*rate(reacn,j)
      end do

      sigb =  dens *  exp ( sigb )
      dummy = sigb
c     neu
      sigb=0.0d0

      do j = 1,7
         sigb = sigb+t09(:,j)*rate(reacr,j)
      end do

c     neu
      sigb=exp(sigb)

      sigb = dens * sigb + dummy

      return
      end subroutine get_rate4

c=======================================================================



c=======================================================================

      subroutine get_rate5 (tmpr, dens, sigb, reacr, reacn)

c=======================================================================
c Autor           : Shu Kitaura
c                   changed by A.Marek
c Modul           :
c Version         :
c Date            : $Date: 2004/14/05 18:37:27 $
c
c purpose: calculates the rate
c
c t913 = t9^(1./3.)
c t953 = t9^(5./3.)
c
c rate = exp(rate(i,1) + rate(i,2)/t9 + rate(i,3)/t913 + rate(i,4)*t913 +
c      rate(i,5)*t9 + rate(i,6)*t953 + rate(i,7)*ln(t9))
c
c     (for the reaction i)
c
c corresponds to ln(2)/t12 for decays, NA<ab> for two-body reactions and
c NA**2<abc> for three-body reactions.
c=======================================================================

      USE precision
      USE burnrate

      IMPLICIT NONE

      integer(kind=ik), intent(in) :: reacr,reacn
      real(kind=rk),    intent(in) :: tmpr(:), dens(:)
      real(kind=rk),    intent(out):: sigb(:)

      real(kind=rk), dimension(size(dens)) ::
     &     tem9,temp9,tem9l,tem9lthd,tem9lfthd,dummy
      real(kind=rk), dimension(size(dens),7) :: t09

      integer(kind=ik) :: j

      real(kind=rk), parameter ::  thrd = 1.0_rk/3.0_rk

      tem9      = tmpr*1.d-9
      tem9l     = log(tem9)
      tem9lthd  = thrd*tem9l
      tem9lfthd = 5.d0*tem9lthd
      t09(:,1)  = 1.d0
      t09(:,2)  = 1.d0/tem9
      t09(:,3)  = exp( -tem9lthd )
      t09(:,4)  = 1.d0/t09(:,3)
      t09(:,5)  = tem9
      t09(:,6)  = exp( tem9lfthd )
      t09(:,7)  = tem9l
      temp9  = tem9


      sigb = 0.0d0

      do j = 1,7
         sigb = sigb+t09(:,j)*rate(reacn,j)
      end do

      sigb =  dens*dens *  exp ( sigb )

      dummy = sigb
      sigb = 0.e0

      do j = 1,7
         sigb = sigb+t09(:,j)*rate(reacr,j)
      end do

      sigb =   exp ( sigb )

      sigb = dens*dens * sigb + dummy

      return
      end subroutine get_rate5

c=======================================================================



      END MODULE mod_get_rate
