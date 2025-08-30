!> \verbatim
!> This module provides some machine dependent procedures
!> \endverbatim
!>
!> \author A. Marek
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
module specfun
  use precision

  public isrmax, isrmin_v, isamin_v,isrmax_v,isamax_v, whenfgt_v, &
         fastexp, fastlog, fastlog10, fastsqrt, fastcos, fastsin

  private

contains

!> \verbatim
!> Calls vector exponential in case fast vector library
!> is available (VMASS on IBM AIX, or MKL in case Intel's
!> Fortran Compiler is used).
!>
!> Calls regular exponential function otherwise
!> \endverbatim
!>
!> \author M. Rampp
!> \param nmax
!> \param xvec
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
  subroutine fastexp ( xvec,nmax )

    use precision
    implicit none

    real(kind=rk), intent(inout) :: xvec(:)
    integer(kind=ik), intent(in) :: nmax

#if WITH_VENDOR_VECLIBS

#if defined(IBM_COMPILER)
    call vexp( xvec,xvec,nmax ) !on IBM machines the fast vector routine
                                !should always be called, KJT rates would
                                !be slowed down considerably otherwise
#elif defined(LINUX)
    call vdexp( nmax,xvec,xvec )
#else
#error No implementation for fastexp found!
#endif

#else /* no WITH_VENDOR_VECLIBS */
    xvec(1:nmax) = exp(xvec(1:nmax))
#endif

    return
  end subroutine fastexp

!>
!> \verbatim
!> Calls vector natural logarithm in case fast vector library
!> is available (VMASS on IBM AIX, or MKL in case Intel's
!> Fortran Compiler is used).
!>
!> Calls regular log function otherwise
!> \endverbatim
!>
!> \author M. Rampp
!> \param nmax
!> \param xvec
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
  subroutine fastlog ( xvec,nmax )

    use precision
    implicit none

    real(kind=rk), intent(inout) :: xvec(:)
    integer(kind=ik), intent(in) :: nmax

#if WITH_VENDOR_VECLIBS

#if defined(IBM_COMPILER)
    call vlog( xvec,xvec,nmax )
#elif defined(LINUX)
    call vdln( nmax,xvec,xvec )
#else
#error No implementation for fastlog found!
#endif

#else /* no WITH_VENDOR_VECLIBS */
    xvec(1:nmax) = log(xvec(1:nmax))
#endif

    return
  end subroutine fastlog

!>
!> \verbatim
!> Calls vector logarithm (log10) in case fast vector library
!> is available (VMASS on IBM AIX, or MKL in case Intel's
!> Fortran Compiler is used).
!>
!> Calls regular log function otherwise
!> \endverbatim
!>
!> \author A. Marek
!> \param nmax
!> \param xvec
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
  subroutine fastlog10 ( xvec,nmax )

    use precision
    implicit none

    real(kind=rk), intent(inout) :: xvec(:)
    integer(kind=ik), intent(in) :: nmax

#if WITH_VENDOR_VECLIBS

#if defined(IBM_COMPILER)
    call vlog10( xvec,xvec,nmax )
#elif defined(LINUX)
    call vdlog10( nmax,xvec,xvec )
#else
#error No implementation for fastlog10 found!
#endif

#else /* no WITH_VENDOR_VECLIBS */
    xvec(1:nmax) = log10(xvec(1:nmax))
#endif

    return
  end subroutine fastlog10

!>
!> \verbatim
!> Calls vector sqrt in case fast vector library
!> is available (VMASS on IBM AIX, or MKL in case Intel's
!> Fortran Compiler is used).
!>
!> Calls regular log function otherwise
!> \endverbatim
!>
!> \author A. Marek
!> \param nmax
!> \param xvec
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
  subroutine fastsqrt ( xvec,nmax )
    use precision

    implicit none

    real(kind=rk), intent(inout) :: xvec(:)
    integer(kind=ik), intent(in) :: nmax

#if WITH_VENDOR_VECLIBS

#if defined(IBM_COMPILER)
    call vsqrt( xvec,xvec,nmax )
#elif defined(LINUX)
    call vdsqrt( nmax,xvec,xvec )
#else
#error No implementation for fastsqrt found!
#endif

#else /* no WITH_VENDOR_VECLIBS */
    xvec(1:nmax) = sqrt(xvec(1:nmax))
#endif

    return
  end subroutine fastsqrt

!>
!> \verbatim
!> Calls vector cos in case fast vector library
!> is available (VMASS on IBM AIX, or MKL in case Intel's
!> Fortran Compiler is used).
!>
!> Calls regular log function otherwise
!> \endverbatim
!>
!> \author A. Marek
!> \param nmax
!> \param xvec
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
  subroutine fastcos ( xvec,nmax )
    use precision

    implicit none

    real(kind=rk) :: xvec(:)
    integer(kind=ik) :: m,nmax
#if WITH_VENDOR_VECLIBS

#if defined(IBM_COMPILER)
    call vcos( xvec,xvec,nmax )
#elif defined(LINUX) && defined(WITH_VENDOR_VECLIBS)
    call vdcos( nmax,xvec,xvec )
#else
#error No implementation for fastcos found!
#endif

#else /* no WITH_VENDOR_VECLIBS */
    xvec(1:nmax) = cos(xvec(1:nmax))
#endif

    return
  end subroutine fastcos

!>
!> \verbatim
!> Calls vector sin in case fast vector library
!> is available (VMASS on IBM AIX, or MKL in case Intel's
!> Fortran Compiler is used).
!>
!> Calls regular log function otherwise
!> \endverbatim
!>
!> \author A. Marek
!> \param nmax
!> \param xvec
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
  subroutine fastsin ( xvec,nmax )
    use precision

    implicit none

    real(kind=rk), intent(inout) :: xvec(:)
    integer(kind=ik), intent(in) :: nmax

#if WITH_VENDOR_VECLIBS

#if defined(IBM_COMPILER)
    call vsin( xvec,xvec,nmax )
#elif defined(LINUX)
    call vdsin( nmax,xvec,xvec )
#else
#error No implementation for fastsin found!
#endif

#else /* no WITH_VENDOR_VECLIBS */
    xvec(1:nmax) = sin(xvec(1:nmax))
#endif

    return
  end subroutine fastsin


!>
!> \verbatim
!> searches a matrix for its maximum value
!>
!> \endverbatim
!>
!> \author M. Rampp
!> \param xx matrix(m,n)
!> \param m
!> \param n
!> \param jm position of maximum in column
!> \param jn position of maximum in row
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
  subroutine isrmax(xx,m,n,jm,jn)

    use precision

!DIR$ INLINEALWAYS isrmax
! searches a MATRIX for its maximum value
    implicit none
! LOCAL variables that are not in modules

!    integer(kind=ik), intrinsic :: maxloc ! forcheck
     integer(kind=ik) ,intent(in) :: m,n
    real(kind=rk) ,intent(in) :: xx(m,n)
    integer(kind=ik) ,intent(out) :: jm,jn

    integer(kind=ik) j(2)

#ifdef CRAY
    ii=ISMAX(m*n,xx,1)-1
    jm=mod(ii,m)+1
    jn=ii/m+1
#else
    j=MAXLOC(xx)
    jm=j(1)
    jn=j(2)
#endif

  end subroutine isrmax

! -------------------------------------------------------------------
!>
!> \verbatim
!> searches a vector for its first location of a minimum
!>
!> \endverbatim
!>
!> \author M. Rampp
!> \param n
!> \param xx
!> \param incx
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
  function isrmin_v(n,xx,incx) RESULT(min_return)
!DIR$ INLINEALWAYS isrmin_v


    use precision

    implicit none
!      integer(kind=ik), intrinsic :: minloc ! forcheck
      integer(kind=ik) :: n,il(1),incx,i,ii,min_return
      real(kind=rk) :: xx(n*incx),x(n)

#ifdef CRAY
      ISRMIN_V = ISMIN(n,xx,incx)
#else
      do i=1,n
         ii=(i-1)*incx+1
         x(i)=xx(ii)
      enddo

      il=MINLOC(x)
      min_return=il(1)
#endif
      return
    end function isrmin_v
! -------------------------------------------------------------------
!>
!> \verbatim
!> searches a vector for its first location of its absolute minimum value
!>
!> \endverbatim
!>
!> \author M. Rampp
!> \param n
!> \param xx
!> \param incx
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
    function isamin_v(n,xx,incx) result(min_return)
!DIR$ INLINEALWAYS isamin_v

      use precision

      implicit none
!      real(kind=rk), intrinsic :: abs ! forechk

      integer(kind=ik) :: n,incx,il(1), min_return,i,ii

      real(kind=rk) :: xx(n*incx),x(n)

#ifdef CRAY
      min_return = ISAMIN(n,xx,incx)
#else
      do i=1,n
         ii=(i-1)*incx+1
         x(i)=abs(xx(ii))
      enddo

      il=MINLOC(x)
      min_return=il(1)
#endif
      return
    end function isamin_v
!>
!> \verbatim
!> searches a vector for its first location of its maximum value
!>
!> \endverbatim
!>
!> \author M. Rampp
!> \param n
!> \param xx
!> \param incx
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
    function isrmax_v(n,xx,incx) result(max_return)
!DIR$ INLINEALWAYS isrmax_v

      use precision

      integer(kind=ik) :: n,incx,il(1), max_return,i,ii
      real(kind=rk) :: xx(n*incx),x(n)


#ifdef CRAY
      max_return = ISMAX(n,xx,incx)
#else
      do i=1,n
         ii=(i-1)*incx+1
         x(i)=xx(ii)
      enddo

      il=MAXLOC(x)
      max_return=il(1)
#endif
      return
    end function isrmax_v

! -------------------------------------------------------------------
!>
!> \verbatim
!> searches a vector for its first location of its maximum absolute value
!>
!> \endverbatim
!>
!> \author M. Rampp
!> \param n
!> \param xx
!> \param incx
!>
!> \verbatim
!>   SVN - Information
!>   $Revision: 527 $
!>   $Date: 2009-11-27 14:19:54 +0100 (Fri, 27 Nov 2009) $
!>
!> \endverbatim
!>
    function isamax_v(n,xx,incx) result(max_return)
!DIR$ INLINEALWAYS isamax_v

      use precision

      implicit none

      integer(kind=ik) :: il(1),n,incx,max_return,i,ii
      real(kind=rk) :: xx(n*incx),x(n)


#ifdef CRAY
      max_return = ISAMAX(n,xx,incx)
#else
      do i=1,n
         ii=(i-1)*incx+1
         x(i)=abs(xx(ii))
      enddo
      il=MAXLOC(x)
      max_return=il(1)
#endif
      return
    end function isamax_v

! -------------------------------------------------------------------
    subroutine whenfgt_v(n,x,incx,target,index,nn)
      use precision

      implicit none
! LOCAL variables that are not in modules


      integer(kind=ik) ,intent(in) :: n,incx
      real(kind=rk) ,intent(in) :: x(n*incx),target
      integer(kind=ik) ,intent(out) :: index(n),nn
      integer(kind=ik) :: i,ii
#ifdef CRAY
!DIR$ INLINE
!DIR$ INLINE
      call whenfgt(n,x,incx,target,index,nn)
#else
      nn=0
      do i=1,n
         ii=(i-1)*incx+1
         if (x(ii).gt.target) then
            nn=nn+1
            index(nn)=i
         endif
      enddo
#endif
    end subroutine whenfgt_v

end module specfun
