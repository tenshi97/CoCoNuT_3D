module lapack

 public dgemm, dgetrf, dgetrs, dgemv

 private


 contains

FUNCTION IEEECK( ISPEC, ZERO, ONE ) result(ieeeck_out)
  
  use precision
  implicit none

!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1998
!
!     .. Scalar Arguments ..
  INTEGER(kind=ik), intent(in) :: ISPEC
  INTEGER(kind=ik)             :: ieeeck_out
  REAL(kind=rk), intent(in)    :: ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
  REAL(kind=rk) ::  NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
                    NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
  IEEECK_OUT = 1
!
  POSINF = ONE / ZERO
  IF( POSINF.LE.ONE ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
!
  NEGINF = -ONE / ZERO
  IF( NEGINF.GE.ZERO ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
!
  NEGZRO = ONE / ( NEGINF+ONE )
  IF( NEGZRO.NE.ZERO ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
!
  NEGINF = ONE / NEGZRO
  IF( NEGINF.GE.ZERO ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
  !
  NEWZRO = NEGZRO + ZERO
  IF( NEWZRO.NE.ZERO ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
  !
  POSINF = ONE / NEWZRO
  IF( POSINF.LE.ONE ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
  !
  NEGINF = NEGINF*POSINF
  IF( NEGINF.GE.ZERO ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
  !
  POSINF = POSINF*POSINF
  IF( POSINF.LE.ONE ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
  IF( ISPEC.EQ.0 ) &
        RETURN
!
  NAN1 = POSINF + NEGINF
!
  NAN2 = POSINF / NEGINF
!
  NAN3 = POSINF / POSINF
!
  NAN4 = POSINF*ZERO
!
  NAN5 = NEGINF*NEGZRO
!
  NAN6 = NAN5*0.0
!
  IF( NAN1.EQ.NAN1 ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
!
  IF( NAN2.EQ.NAN2 ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
!
  IF( NAN3.EQ.NAN3 ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
!
  IF( NAN4.EQ.NAN4 ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
!
  IF( NAN5.EQ.NAN5 ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
!
  IF( NAN6.EQ.NAN6 ) THEN
     IEEECK_OUT = 0
     RETURN
  END IF
!yyy
  RETURN
END FUNCTION IEEECK

FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 ) result(ilaenv_out)
  
  use precision
  implicit none

  !
  !  -- LAPACK auxiliary routine (version 3.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     June 30, 1999
  !
  !     .. Scalar Arguments ..
  CHARACTER*( * ), intent(in)  :: NAME, OPTS
  INTEGER(kind=ik), intent(in) :: ISPEC, N1, N2, N3, N4
  integer(kind=ik)             :: ilaenv_out
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  ILAENV is called from the LAPACK routines to choose problem-dependent
  !  parameters for the local environment.  See ISPEC for a description of
  !  the parameters.
  !
  !  This version provides a set of parameters which should give good,
  !  but not optimal, performance on many of the currently available
  !  computers.  Users are encouraged to modify this subroutine to set
  !  the tuning parameters for their particular machine using the option
  !  and problem size information in the arguments.
  !
  !  This routine will not function correctly if it is converted to all
  !  lower case.  Converting it to all upper case is allowed.
  !
  !  Arguments
  !  =========
  !
  !  ISPEC   (input) INTEGER
  !          Specifies the parameter to be returned as the value of
  !          ILAENV.
  !          = 1: the optimal blocksize; if this value is 1, an unblocked
  !               algorithm will give the best performance.
  !          = 2: the minimum block size for which the block routine
  !               should be used; if the usable block size is less than
  !               this value, an unblocked routine should be used.
  !          = 3: the crossover point (in a block routine, for N less
  !               than this value, an unblocked routine should be used)
  !          = 4: the number of shifts, used in the nonsymmetric
  !               eigenvalue routines
  !          = 5: the minimum column dimension for blocking to be used;
  !               rectangular blocks must have dimension at least k by m,
  !               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
  !          = 6: the crossover point for the SVD (when reducing an m by n
  !               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
  !               this value, a QR factorization is used first to reduce
  !               the matrix to a triangular form.)
  !          = 7: the number of processors
  !          = 8: the crossover point for the multishift QR and QZ methods
  !               for nonsymmetric eigenvalue problems.
  !          = 9: maximum size of the subproblems at the bottom of the
  !               computation tree in the divide-and-conquer algorithm
  !               (used by xGELSD and xGESDD)
  !          =10: ieee NaN arithmetic can be trusted not to trap
  !          =11: infinity arithmetic can be trusted not to trap
  !
  !  NAME    (input) CHARACTER*(*)
  !          The name of the calling subroutine, in either upper case or
  !          lower case.
  !
  !  OPTS    (input) CHARACTER*(*)
  !          The character options to the subroutine NAME, concatenated
  !          into a single character string.  For example, UPLO = 'U',
  !          TRANS = 'T', and DIAG = 'N' for a triangular routine would
  !          be specified as OPTS = 'UTN'.
  !
  !  N1      (input) INTEGER
  !  N2      (input) INTEGER
  !  N3      (input) INTEGER
  !  N4      (input) INTEGER
  !          Problem dimensions for the subroutine NAME; these may not all
  !          be required.
  !
  ! (ILAENV) (output) INTEGER
  !          >= 0: the value of the parameter specified by ISPEC
  !          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
  !
  !  Further Details
  !  ===============
  !
  !  The following conventions have been used when calling ILAENV from the
  !  LAPACK routines:
  !  1)  OPTS is a concatenation of all of the character options to
  !      subroutine NAME, in the same order that they appear in the
  !      argument list for NAME, even if they are not used in determining
  !      the value of the parameter specified by ISPEC.
  !  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
  !      that they appear in the argument list for NAME.  N1 is used
  !      first, N2 second, and so on, and unused problem dimensions are
  !      passed a value of -1.
  !  3)  The parameter value returned by ILAENV is checked for validity in
  !      the calling subroutine.  For example, ILAENV is used to retrieve
  !      the optimal blocksize for STRTRI as follows:
  !
  !      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
  !      IF( NB.LE.1 ) NB = MAX( 1, N )
  !
  !  =====================================================================
  !
  !     .. Local Scalars ..
  LOGICAL          :: CNAME, SNAME
  CHARACTER*1      :: C1
  CHARACTER*2      :: C2, C4
  CHARACTER*3      :: C3
  CHARACTER*6      :: SUBNAM
  INTEGER(kind=ik) :: I, IC, IZ, NB, NBMIN, NX
  !     ..
  !     .. Intrinsic Functions ..
  !      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
  !     ..
  !     .. External Functions ..
  !      INTEGER            IEEECK
  !      EXTERNAL           IEEECK
  !     ..
  !     .. Executable Statements ..
  !
  GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000, &
          1100 ) ISPEC
  !
  !     Invalid value for ISPEC
  !
  ILAENV_OUT = -1
  RETURN
  !
100 CONTINUE
  !
  !     Convert NAME to upper case if the first character is lower case.
  !
  ILAENV_OUT = 1
  SUBNAM = NAME
  IC = ICHAR( SUBNAM( 1:1 ) )
  IZ = ICHAR( 'Z' )
  IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
     !
     !        ASCII character set
     !
     IF( IC.GE.97 .AND. IC.LE.122 ) THEN
        SUBNAM( 1:1 ) = CHAR( IC-32 )
        DO I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           IF( IC.GE.97 .AND. IC.LE.122 ) &
                SUBNAM( I:I ) = CHAR( IC-32 )
        enddo
     END IF
     !
  ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
     !
     !        EBCDIC character set
     !
     IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.  &
          ( IC.GE.145 .AND. IC.LE.153 ) .OR.  &
          ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
        SUBNAM( 1:1 ) = CHAR( IC+64 )
        DO I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                ( IC.GE.162 .AND. IC.LE.169 ) )    &
                SUBNAM( I:I ) = CHAR( IC+64 )
        enddo
     END IF
     !
  ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
     !
     !        Prime machines:  ASCII+128
     !
     IF( IC.GE.225 .AND. IC.LE.250 ) THEN
        SUBNAM( 1:1 ) = CHAR( IC-32 )
        DO I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           IF( IC.GE.225 .AND. IC.LE.250 ) &
                SUBNAM( I:I ) = CHAR( IC-32 )
        enddo
     END IF
  END IF
  !
  C1 = SUBNAM( 1:1 )
  SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
  CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
  IF( .NOT.( CNAME .OR. SNAME ) ) &
       RETURN
  C2 = SUBNAM( 2:3 )
  C3 = SUBNAM( 4:6 )
  C4 = C3( 2:3 )
  !
  GO TO ( 110, 200, 300 ) ISPEC
  !
110 CONTINUE
  !
  !     ISPEC = 1:  block size
  !
  !     In these examples, separate code is provided for setting NB for
  !     real and complex.  We assume that NB will take the same value in
  !     single or double precision.
  !
  NB = 1
  !
  IF( C2.EQ.'GE' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
          C3.EQ.'QLF' ) THEN
        IF( SNAME ) THEN
           NB = 32
        ELSE
           NB = 32
        END IF
     ELSE IF( C3.EQ.'HRD' ) THEN
        IF( SNAME ) THEN
           NB = 32
        ELSE
           NB = 32
        END IF
     ELSE IF( C3.EQ.'BRD' ) THEN
        IF( SNAME ) THEN
           NB = 32
        ELSE
           NB = 32
        END IF
     ELSE IF( C3.EQ.'TRI' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( C2.EQ.'PO' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( C2.EQ.'SY' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
        NB = 32
     ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
        NB = 64
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        NB = 64
     ELSE IF( C3.EQ.'TRD' ) THEN
        NB = 32
     ELSE IF( C3.EQ.'GST' ) THEN
        NB = 64
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     END IF
  ELSE IF( C2.EQ.'GB' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           IF( N4.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        ELSE
           IF( N4.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        END IF
     END IF
  ELSE IF( C2.EQ.'PB' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           IF( N2.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        ELSE
           IF( N2.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        END IF
     END IF
  ELSE IF( C2.EQ.'TR' ) THEN
     IF( C3.EQ.'TRI' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( C2.EQ.'LA' ) THEN
     IF( C3.EQ.'UUM' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
     IF( C3.EQ.'EBZ' ) THEN
        NB = 1
     END IF
  END IF
  ILAENV_OUT = NB
  RETURN
  !
200 CONTINUE
  !
  !     ISPEC = 2:  minimum block size
  !
  NBMIN = 2
  IF( C2.EQ.'GE' ) THEN
     IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
          C3.EQ.'QLF' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     ELSE IF( C3.EQ.'HRD' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     ELSE IF( C3.EQ.'BRD' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     ELSE IF( C3.EQ.'TRI' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     END IF
  ELSE IF( C2.EQ.'SY' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NBMIN = 8
        ELSE
           NBMIN = 8
        END IF
     ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
        NBMIN = 2
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
     IF( C3.EQ.'TRD' ) THEN
        NBMIN = 2
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     END IF
  END IF
  ILAENV_OUT = NBMIN
  RETURN
  !
300 CONTINUE
  !
  !     ISPEC = 3:  crossover point
  !
  NX = 0
  IF( C2.EQ.'GE' ) THEN
     IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
          C3.EQ.'QLF' ) THEN
        IF( SNAME ) THEN
           NX = 128
        ELSE
           NX = 128
        END IF
     ELSE IF( C3.EQ.'HRD' ) THEN
        IF( SNAME ) THEN
           NX = 128
        ELSE
           NX = 128
        END IF
     ELSE IF( C3.EQ.'BRD' ) THEN
        IF( SNAME ) THEN
           NX = 128
        ELSE
           NX = 128
        END IF
     END IF
  ELSE IF( C2.EQ.'SY' ) THEN
     IF( SNAME .AND. C3.EQ.'TRD' ) THEN
        NX = 32
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
     IF( C3.EQ.'TRD' ) THEN
        NX = 32
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             C4.EQ.'BR' ) THEN
           NX = 128
        END IF
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             C4.EQ.'BR' ) THEN
           NX = 128
        END IF
     END IF
  END IF
  ILAENV_OUT = NX
  RETURN
  !
400 CONTINUE
  !
  !     ISPEC = 4:  number of shifts (used by xHSEQR)
  !
  ILAENV_OUT = 6
  RETURN
  !
500 CONTINUE
  !
  !     ISPEC = 5:  minimum column dimension (not used)
  !
  ILAENV_OUT = 2
  RETURN
  !
600 CONTINUE 
  !
  !     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
  !
  ILAENV_OUT = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
  RETURN
  !
700 CONTINUE
  !
  !     ISPEC = 7:  number of processors (not used)
  !
  ILAENV_OUT = 1
  RETURN
  !
800 CONTINUE
  !
  !     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
  !
  ILAENV_OUT = 50
  RETURN
  !
900 CONTINUE
  !
  !     ISPEC = 9:  maximum size of the subproblems at the bottom of the
  !                 computation tree in the divide-and-conquer algorithm
  !                 (used by xGELSD and xGESDD)
  !
  ILAENV_OUT = 25
  RETURN
  !
1000 CONTINUE
  !
  !     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
  !
  !     ILAENV_OUT = 0
  ILAENV_OUT = 1
  IF( ILAENV_OUT.EQ.1 ) THEN
     ILAENV_OUT = IEEECK( 0, 0.0_rk, 1.0_rk) 
  END IF
  RETURN
  !
1100 CONTINUE
  !
  !     ISPEC = 11: infinity arithmetic can be trusted not to trap
  !
  !     ILAENV_OUT = 0
  ILAENV_OUT = 1
  IF( ILAENV_OUT.EQ.1 ) THEN
     ILAENV_OUT = IEEECK( 1, 0.0_rk, 1.0_rk) 
  END IF
  RETURN
  !
  !     End of ILAENV
  !
END FUNCTION ILAENV

 
FUNCTION LSAME(CA,CB) result(lsame_out)
  use precision
  
  implicit none

!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
  CHARACTER, intent(in) :: CA,CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
!  INTRINSIC ICHAR
!     ..
!     .. Local Scalars ..
  INTEGER(kind=ik) :: INTA,INTB,ZCODE
  logical          :: lsame_out
!     ..
!
!     Test if the characters are equal
!
  LSAME_OUT = CA .EQ. CB
  IF (LSAME_OUT) RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
  ZCODE = ICHAR('Z')
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
  INTA = ICHAR(CA)
  INTB = ICHAR(CB)
!
  IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
     IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
     IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
!
  ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
     IF (INTA.GE.129 .AND. INTA.LE.137 .OR. &
         INTA.GE.145 .AND. INTA.LE.153 .OR. &
         INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
     IF (INTB.GE.129 .AND. INTB.LE.137 .OR. &
         INTB.GE.145 .AND. INTB.LE.153 .OR. &
         INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
!
  ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
     IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
     IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
  END IF

  LSAME_OUT = INTA .EQ. INTB
!
!     RETURN
!
!     End of LSAME

END FUNCTION LSAME


SUBROUTINE XERBLA( SRNAME, INFO )
  use precision
  use abort

  implicit none


!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
  CHARACTER*6, intent(in)       :: SRNAME
  INTEGER(kind=ik), intent(in)  :: INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========


!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Executable Statements ..
!
  WRITE( *, FMT = 9999 )SRNAME, INFO
!
  raise_abort("Lapack error handler called")
!
9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
            'an illegal value' )
!
!     End of XERBLA
!
END SUBROUTINE XERBLA




SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

  use precision
  
  implicit none

!     .. Scalar Arguments ..
  real(kind=rk), intent(in)    :: ALPHA,BETA
  integer(kind=ik), intent(in) :: LDA, LDB, LDC

  integer(kind=ik)             :: m, n, k
  CHARACTER, intent(in)        ::  TRANSA, TRANSB
!     ..
!     .. Array Arguments ..
  real(kind=rk), intent(in)    :: A(LDA,*), B(LDB,*)
  real(kind=rk), intent(inout) :: C(LDC,*)
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!  =====================================================================
!
!     .. External Functions ..
!  LOGICAL, external :; LSAME
!     ..
!     .. External Subroutines ..
!  EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!  INTRINSIC MAX
!     ..
!     .. Local Scalars ..
  real(kind=rk)    :: TEMP
  integer(kind=ik) :: I,INFO,J,L,NCOLA,NROWA,NROWB
  LOGICAL          :: NOTA,NOTB
!     ..
!     .. Parameters ..
  real(kind=rk), parameter :: ONE=1.0_rk,ZERO=0.0_rk

!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
  NOTA = LSAME(TRANSA,'N')
  NOTB = LSAME(TRANSB,'N')
  IF (NOTA) THEN
     NROWA = M
     NCOLA = K
  ELSE
     NROWA = K
     NCOLA = M
  END IF
  IF (NOTB) THEN
     NROWB = K
  ELSE
     NROWB = N
  END IF
!
!     Test the input parameters.
!
  INFO = 0
  IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. &
      (.NOT.LSAME(TRANSA,'T'))) THEN
     INFO = 1
  ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. &
           (.NOT.LSAME(TRANSB,'T'))) THEN
     INFO = 2
  ELSE IF (M.LT.0) THEN
     INFO = 3
  ELSE IF (N.LT.0) THEN
     INFO = 4
  ELSE IF (K.LT.0) THEN
     INFO = 5
  ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
     INFO = 8
  ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
     INFO = 10
  ELSE IF (LDC.LT.MAX(1,M)) THEN
     INFO = 13
  END IF
  IF (INFO.NE.0) THEN
     CALL XERBLA('DGEMM ',INFO)
     RETURN
  END IF
  !
!     Quick return if possible.
!
  IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
      (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
!
!     And if  alpha.eq.zero.
!
  IF (ALPHA.EQ.ZERO) THEN
     IF (BETA.EQ.ZERO) THEN
        DO J = 1,N
           DO I = 1,M
              C(I,J) = ZERO
           enddo
        enddo
     ELSE
        DO J = 1,N
           DO I = 1,M
              C(I,J) = BETA*C(I,J)
           enddo
        enddo
     END IF
     RETURN
  END IF
!
!     Start the operations.
!
  IF (NOTB) THEN
     IF (NOTA) THEN
!
!           Form  C := alpha*A*B + beta*C.
!
        DO J = 1,N
           IF (BETA.EQ.ZERO) THEN
              DO I = 1,M
                 C(I,J) = ZERO
              enddo
           ELSE IF (BETA.NE.ONE) THEN
              DO I = 1,M
                 C(I,J) = BETA*C(I,J)
              enddo
           END IF
           DO L = 1,K
              IF (B(L,J).NE.ZERO) THEN
                 TEMP = ALPHA*B(L,J)
                 DO I = 1,M
                    C(I,J) = C(I,J) + TEMP*A(I,L)
                 enddo
              END IF
           enddo
        enddo
     ELSE
        !
        !           Form  C := alpha*A'*B + beta*C
        !
        DO J = 1,N
           DO I = 1,M
              TEMP = ZERO
              DO L = 1,K
                 TEMP = TEMP + A(L,I)*B(L,J)
              enddo
              IF (BETA.EQ.ZERO) THEN
                 C(I,J) = ALPHA*TEMP
              ELSE
                 C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
           enddo
        enddo
     END IF
  ELSE
     IF (NOTA) THEN
        !
        !           Form  C := alpha*A*B' + beta*C
        !
        DO J = 1,N
           IF (BETA.EQ.ZERO) THEN
              DO I = 1,M
                 C(I,J) = ZERO
              enddo
           ELSE IF (BETA.NE.ONE) THEN
              DO I = 1,M
                 C(I,J) = BETA*C(I,J)
              enddo
           END IF
           DO L = 1,K
              IF (B(J,L).NE.ZERO) THEN
                 TEMP = ALPHA*B(J,L)
                 DO I = 1,M
                    C(I,J) = C(I,J) + TEMP*A(I,L)
                 enddo
              END IF
           enddo
        enddo
     ELSE
        !
        !           Form  C := alpha*A'*B' + beta*C
        !
        DO J = 1,N
           DO I = 1,M
              TEMP = ZERO
              DO L = 1,K
                 TEMP = TEMP + A(L,I)*B(J,L)
              enddo
              IF (BETA.EQ.ZERO) THEN
                 C(I,J) = ALPHA*TEMP
              ELSE
                 C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
           enddo
        enddo
     END IF
  END IF
  !
  RETURN
  !
  !     End of DGEMM .
                                                  !
END SUBROUTINE DGEMM


SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)

  use precision

  implicit none
!     .. Scalar Arguments ..
  REAL(kind=rk), intent(in)    :: ALPHA
  INTEGER(kind=ik), intent(in) :: INCX, INCY, LDA, M, N
!     ..
!     .. Array Arguments ..
  REAL(kind=rk), intent(inout) :: A(LDA,*)
  REAL(kind=rk), intent(inout) :: X(*), Y(*)
!     ..
!
!  Purpose
!  =======
!
!  SGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Arguments
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL             array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - REAL             array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!  =====================================================================
!
!     .. Parameters ..
  REAL(kind=rk), parameter :: ZERO=0.0_rk
!     ..
!     .. Local Scalars ..
  REAL(kind=rk)            :: TEMP
  INTEGER(kind=ik)         :: I, INFO, IX, J, JY, KX
!     ..
!     .. External Subroutines ..
!  EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!  INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
  INFO = 0
  IF (M.LT.0) THEN
     INFO = 1
  ELSE IF (N.LT.0) THEN
     INFO = 2
  ELSE IF (INCX.EQ.0) THEN
     INFO = 5
  ELSE IF (INCY.EQ.0) THEN
     INFO = 7
  ELSE IF (LDA.LT.MAX(1,M)) THEN
     INFO = 9
  END IF
  IF (INFO.NE.0) THEN
     CALL XERBLA('SGER  ',INFO)
     RETURN
  END IF
!
!     Quick return if possible.
!
  IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
  IF (INCY.GT.0) THEN
     JY = 1
  ELSE
     JY = 1 - (N-1)*INCY
  END IF
  IF (INCX.EQ.1) THEN
     DO J = 1,N
        IF (Y(JY).NE.ZERO) THEN
           TEMP = ALPHA*Y(JY)
           DO I = 1,M
              A(I,J) = A(I,J) + X(I)*TEMP
           enddo
        END IF
        JY = JY + INCY
     enddo
  ELSE
     IF (INCX.GT.0) THEN
        KX = 1
     ELSE
        KX = 1 - (M-1)*INCX
     END IF
     DO J = 1,N
        IF (Y(JY).NE.ZERO) THEN
           TEMP = ALPHA*Y(JY)
           IX = KX
           DO I = 1,M
              A(I,J) = A(I,J) + X(IX)*TEMP
              IX = IX + INCX
           enddo
        END IF
        JY = JY + INCY
     enddo
  END IF
!
  RETURN
!
!     End of SGER  .
!
END SUBROUTINE SGER


SUBROUTINE SSCAL(N,SA,SX,INCX)

  use precision

  implicit none
!     .. Scalar Arguments ..
  REAL(kind=rk), intent(in)    :: SA
  INTEGER(kind=ik), intent(in) :: INCX, N
  !     ..
  !     .. Array Arguments ..
  REAL(kind=rk), intent(inout) :: SX(*)
!     ..
!
!  Purpose
!  =======
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to 1.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(:)
!
!  =====================================================================
!
!     .. Local Scalars ..
  INTEGER(kind=ik) :: I, M, MP1, NINCX
!     ..
!     .. Intrinsic Functions ..
!  INTRINSIC MOD
!     ..
  IF (N.LE.0 .OR. INCX.LE.0) RETURN
  IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
  NINCX = N*INCX
  DO I = 1,NINCX,INCX
     SX(I) = SA*SX(I)
  enddo
  RETURN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
20 M = MOD(N,5)
  IF (M.EQ.0) GO TO 40
  DO I = 1,M
     SX(I) = SA*SX(I)
  end DO
  IF (N.LT.5) RETURN
40 MP1 = M + 1
  DO I = MP1,N,5
     SX(I) = SA*SX(I)
     SX(I+1) = SA*SX(I+1)
     SX(I+2) = SA*SX(I+2)
     SX(I+3) = SA*SX(I+3)
     SX(I+4) = SA*SX(I+4)
  end DO
  RETURN
END SUBROUTINE SSCAL


SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)

  use precision

  implicit none

!     .. Scalar Arguments ..
  INTEGER(kind=ik), intent(in) :: INCX, INCY, N
!     ..
!     .. Array Arguments ..
  REAL(kind=rk), intent(inout) :: SX(*), SY(*)
!     ..
!
!  Purpose
!  =======
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal to 1.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(:)
!
!  =====================================================================
!
!     .. Local Scalars ..
  REAL(kind=rk)    :: STEMP
  INTEGER(kind=ik) :: I, IX, IY, M, MP1
!     ..
!     .. Intrinsic Functions ..
!  INTRINSIC MOD
!     ..
  IF (N.LE.0) RETURN
  IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
  IX = 1
  IY = 1
  IF (INCX.LT.0) IX = (-N+1)*INCX + 1
  IF (INCY.LT.0) IY = (-N+1)*INCY + 1
  DO I = 1,N
     STEMP = SX(IX)
     SX(IX) = SY(IY)
     SY(IY) = STEMP
     IX = IX + INCX
     IY = IY + INCY
  enddo
  RETURN
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
20 M = MOD(N,3)
  IF (M.EQ.0) GO TO 40
  DO I = 1,M
     STEMP = SX(I)
     SX(I) = SY(I)
     SY(I) = STEMP
  enddo
  IF (N.LT.3) RETURN
40 MP1 = M + 1
  DO I = MP1,N,3
     STEMP = SX(I)
     SX(I) = SY(I)
     SY(I) = STEMP
     STEMP = SX(I+1)
     SX(I+1) = SY(I+1)
     SY(I+1) = STEMP
     STEMP = SX(I+2)
     SX(I+2) = SY(I+2)
     SY(I+2) = STEMP
  enddo
  RETURN
END SUBROUTINE SSWAP


SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
  use precision
  
  implicit none

!     .. Scalar Arguments ..
  REAL(kind=rk), intent(in)    :: ALPHA
  INTEGER(kind=ik), intent(in) :: LDA, LDB, M, N
  CHARACTER, intent(in)        :: DIAG, SIDE, TRANSA, UPLO
!     ..
!     .. Array Arguments ..
  REAL(kind=rk), intent(in)    :: A(LDA,*)
  REAL(kind=rk), intent(inout) :: B(LDB,*)
!     ..
!
!  Purpose
!  =======
!
!  STRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Arguments
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - REAL             array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!  =====================================================================
!
  !     .. External Functions ..
  !      LOGICAL LSAME
  !      EXTERNAL LSAME
!     ..
  !     .. External Subroutines ..
  ! EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!   INTRINSIC MAX
!     ..
!     .. Local Scalars ..
  REAL(kind=rk)    :: TEMP
  INTEGER(kind=ik) :: I, INFO, J, K, NROWA
  LOGICAL          :: LSIDE, NOUNIT, UPPER
!     ..
!     .. Parameters ..
  REAL(kind=rk), parameter :: ONE=1.0_rk, ZERO=0.0_rk
!     ..
!
!     Test the input parameters.
!
  LSIDE = LSAME(SIDE,'L')
  IF (LSIDE) THEN
     NROWA = M
  ELSE
     NROWA = N
  END IF
  NOUNIT = LSAME(DIAG,'N')
  UPPER = LSAME(UPLO,'U')
  !
  INFO = 0
  IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
     INFO = 1
  ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
     INFO = 2
  ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
           (.NOT.LSAME(TRANSA,'T')) .AND. &
           (.NOT.LSAME(TRANSA,'C'))) THEN
     INFO = 3
  ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
     INFO = 4
  ELSE IF (M.LT.0) THEN
     INFO = 5
  ELSE IF (N.LT.0) THEN
     INFO = 6
  ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
     INFO = 9
  ELSE IF (LDB.LT.MAX(1,M)) THEN
     INFO = 11
  END IF
  IF (INFO.NE.0) THEN
     CALL XERBLA('STRSM ',INFO)
     RETURN
  END IF
!
!     Quick return if possible.
!
  IF (M.EQ.0 .OR. N.EQ.0) RETURN
!
!     And when  alpha.eq.zero.
!
  IF (ALPHA.EQ.ZERO) THEN
     DO J = 1,N
        DO I = 1,M
           B(I,J) = ZERO
        enddo
     enddo
     RETURN
  END IF
!
!     Start the operations.
!
  IF (LSIDE) THEN
     IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*inv( A )*B.
!
        IF (UPPER) THEN
           DO J = 1,N
              IF (ALPHA.NE.ONE) THEN
                 DO I = 1,M
                    B(I,J) = ALPHA*B(I,J)
                 enddo
              END IF
              DO K = M,1,-1
                 IF (B(K,J).NE.ZERO) THEN
                    IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                    DO I = 1,K - 1
                       B(I,J) = B(I,J) - B(K,J)*A(I,K)
                    enddo
                 END IF
              enddo
           enddo
        ELSE
           DO J = 1,N
              IF (ALPHA.NE.ONE) THEN
                 DO I = 1,M
                    B(I,J) = ALPHA*B(I,J)
                 enddo
              END IF
              DO K = 1,M
                 IF (B(K,J).NE.ZERO) THEN
                    IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                    DO I = K + 1,M
                       B(I,J) = B(I,J) - B(K,J)*A(I,K)
                    enddo
                 END IF
              enddo
           enddo
        END IF
     ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
        IF (UPPER) THEN
           DO J = 1,N
              DO I = 1,M
                 TEMP = ALPHA*B(I,J)
                 DO K = 1,I - 1
                    TEMP = TEMP - A(K,I)*B(K,J)
                 enddo

                 IF (NOUNIT) TEMP = TEMP/A(I,I)
                 B(I,J) = TEMP
              enddo
           enddo
        ELSE
           DO J = 1,N
              DO I = M,1,-1
                 TEMP = ALPHA*B(I,J)
                 DO K = I + 1,M
                    TEMP = TEMP - A(K,I)*B(K,J)
                 enddo
                 IF (NOUNIT) TEMP = TEMP/A(I,I)
                 B(I,J) = TEMP
              enddo
           enddo
        END IF
     END IF
  ELSE
     IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*B*inv( A ).
!
        IF (UPPER) THEN
           DO J = 1,N
              IF (ALPHA.NE.ONE) THEN
                 DO I = 1,M
                    B(I,J) = ALPHA*B(I,J)
                 ENDDO
              END IF
              DO K = 1,J - 1
                 IF (A(K,J).NE.ZERO) THEN
                    DO I = 1,M
                       B(I,J) = B(I,J) - A(K,J)*B(I,K)
                    ENDDO
                 END IF
              ENDDO
              IF (NOUNIT) THEN
                 TEMP = ONE/A(J,J)
                 DO I = 1,M
                    B(I,J) = TEMP*B(I,J)
                 ENDDO
              END IF
           ENDDO
        ELSE
           DO J = N,1,-1
              IF (ALPHA.NE.ONE) THEN
                 DO I = 1,M
                    B(I,J) = ALPHA*B(I,J)
                 ENDDO
              END IF
              DO K = J + 1,N
                 IF (A(K,J).NE.ZERO) THEN
                    DO I = 1,M
                       B(I,J) = B(I,J) - A(K,J)*B(I,K)
                    ENDDO
                 END IF
              ENDDO
              IF (NOUNIT) THEN
                 TEMP = ONE/A(J,J)
                 DO I = 1,M
                    B(I,J) = TEMP*B(I,J)
                 ENDDO
              END IF
           ENDDO
        END IF
     ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
        IF (UPPER) THEN
           DO K = N,1,-1
              IF (NOUNIT) THEN
                 TEMP = ONE/A(K,K)
                 DO I = 1,M
                    B(I,K) = TEMP*B(I,K)
                 ENDDO
              END IF
              DO J = 1,K - 1
                 IF (A(J,K).NE.ZERO) THEN
                    TEMP = A(J,K)
                    DO I = 1,M
                       B(I,J) = B(I,J) - TEMP*B(I,K)
                    ENDDO
                 END IF
              ENDDO
              IF (ALPHA.NE.ONE) THEN
                 DO I = 1,M
                    B(I,K) = ALPHA*B(I,K)
                 enddo
              END IF
           enddo
        ELSE
           DO K = 1,N
              IF (NOUNIT) THEN
                 TEMP = ONE/A(K,K)
                 DO I = 1,M
                    B(I,K) = TEMP*B(I,K)
                 enddo
              END IF
              DO J = K + 1,N
                 IF (A(J,K).NE.ZERO) THEN
                    TEMP = A(J,K)
                    DO I = 1,M
                       B(I,J) = B(I,J) - TEMP*B(I,K)
                    enddo
                 END IF
              enddo
              IF (ALPHA.NE.ONE) THEN
                 DO I = 1,M
                    B(I,K) = ALPHA*B(I,K)
                 enddo
              END IF
           enddo
        END IF
     END IF
  END IF
  !
  RETURN
  !
  !     End of STRSM .
  !
END SUBROUTINE STRSM


FUNCTION ISAMAX(N,SX,INCX) result(ISAMAX_OUT)

  use precision
  
  implicit none
!     .. Scalar Arguments ..
  INTEGER(kind=ik), intent(in)  :: INCX, N
!     ..
!     .. Array Arguments ..
  REAL(kind=rk), intent(in)     :: SX(*)
  INTEGER(kind=ik)              :: isamax_out  
!     ..
!
!  Purpose
!  =======
!
!     ISAMAX finds the index of element having max. absolute value.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(:)
!
!  =====================================================================
!
!     .. Local Scalars ..
  REAL(kind=rk)    :: SMAX
  INTEGER(kind=ik) :: I, IX
!     ..
!     .. Intrinsic Functions ..
!  INTRINSIC ABS
!     ..
  ISAMAX_OUT = 0
  IF (N.LT.1 .OR. INCX.LE.0) RETURN
  ISAMAX_OUT = 1
  IF (N.EQ.1) RETURN
  IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
  IX = 1
  SMAX = ABS(SX(1))
  IX = IX + INCX
  DO I = 2,N
     IF (ABS(SX(IX)).LE.SMAX) GO TO 5
     ISAMAX_OUT = I
     SMAX = ABS(SX(IX))
5    IX = IX + INCX
  enddo
  RETURN
!
!        code for increment equal to 1
!
20 SMAX = ABS(SX(1))
  DO I = 2,N
     IF (ABS(SX(I)).LE.SMAX) GO TO 30
     ISAMAX_OUT = I
     SMAX = ABS(SX(I))
  enddo
30 continue
  RETURN
END FUNCTION ISAMAX

SUBROUTINE DEGTRT( M, N, A, LDA, IPIV, INFO )

  use precision

  implicit none
!
!  -- LAPACK routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
  INTEGER(kind=ik), intent(out) :: INFO
  INTEGER(kind=ik), intent(in)  ::  LDA, M, N
!     ..
!     .. Array Arguments ..
  INTEGER(kind=ik), intent(out) :: IPIV(*)
  REAL(kind=rk), intent(inout)  :: A(LDA,*)
!     ..
!
!  Purpose
!  =======
!
!  DEGTRT computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
  REAL(kind=rk), parameter :: ONE=1.0_rk
!     ..
!     .. Local Scalars ..
  INTEGER(kind=ik)         :: I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMM, DGETF2, SLASWP, STRSM, XERBLA
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DEGTRT', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
!     Determine the block size for this environment.
!
  NB = ILAENV( 1, 'DEGTRT', ' ', M, N, -1, -1 )
  IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
     CALL DGETF2( M, N, A, LDA, IPIV, INFO )
  ELSE
!
!        Use blocked code.
!
     DO J = 1, MIN( M, N ), NB
        JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
        CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
        IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + J - 1
        DO I = J, MIN( M, J+JB-1 )
           IPIV( I ) = J - 1 + IPIV( I )
        enddo
!
!           Apply interchanges to columns 1:J-1.
!
        CALL SLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
        IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
           CALL SLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, IPIV, 1 )
!
!              Compute block row of U.
!
           CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                        N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),&
                        LDA )
           IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
              CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,   &
                           N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,     &
                           A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),   &
                           LDA )
           END IF
        END IF
     enddo
  END IF
  RETURN
!
!     End of DEGTRT
!
END SUBROUTINE DEGTRT


SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
  use precision
  
  implicit none

!
!  -- LAPACK routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
  CHARACTER, intent(in)         :: TRANS
  INTEGER(kind=ik), intent(out) :: INFO
  INTEGER(kind=ik), intent(in)  :: LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
  INTEGER(kind=ik), intent(in)  :: IPIV(*)
  REAL(kind=rk), intent(in)     :: A(LDA,*)
  REAL(kind=rk), intent(inout)  :: B(LDB,*)
!     ..
!
!  Purpose
!  =======
!
!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by DEGTRT.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) REAL array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DEGTRT.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DEGTRT; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
  REAL(kind=rk), parameter :: ONE=1.0_rk
!     ..
!     .. Local Scalars ..
  LOGICAL ::           NOTRAN
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           SLASWP, STRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  NOTRAN = LSAME( TRANS, 'N' )
  IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
         LSAME( TRANS, 'C' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( NRHS.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
     INFO = -5
  ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
     INFO = -8
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DGETRS', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
        RETURN
!
  IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
     CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
         !
!        Solve L*X = B, overwriting B with X.
!
     CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
                 ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
     CALL STRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
                   NRHS, ONE, A, LDA, B, LDB )
  ELSE
!
!        Solve A' * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
     CALL STRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
                    ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
     CALL STRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, &
          A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
     CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
  END IF
  !
  RETURN
!
!     End of DGETRS
!
END SUBROUTINE DGETRS

   
SUBROUTINE SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
  USE PRECISION
  
  IMPLICIT NONE

!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  INTEGER(kind=ik), intent(in) :: INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
  INTEGER(kind=ik), intent(in) :: IPIV(*)
  REAL(kind=rk), intent(inout) :: A(LDA,*)
!     ..
!
!  Purpose
!  =======
!
!  SLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
  INTEGER(kind=ik) :: I, I1, I2, INC, IP, IX, IX0, J, K, N32
  REAL(kind=rk)    :: TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
  IF( INCX.GT.0 ) THEN
     IX0 = K1
     I1 = K1
     I2 = K2
     INC = 1
  ELSE IF( INCX.LT.0 ) THEN
     IX0 = 1 + ( 1-K2 )*INCX
     I1 = K2
     I2 = K1
     INC = -1
  ELSE
     RETURN
  END IF
!
  N32 = ( N / 32 )*32
  IF( N32.NE.0 ) THEN
     DO J = 1, N32, 32
        IX = IX0
        DO I = I1, I2, INC
           IP = IPIV( IX )
           IF( IP.NE.I ) THEN
              DO K = J, J + 31
                 TEMP = A( I, K )
                 A( I, K ) = A( IP, K )
                 A( IP, K ) = TEMP
              enddo
           END IF
           IX = IX + INCX
        enddo
     enddo
  END IF
  IF( N32.NE.N ) THEN
     N32 = N32 + 1
     IX = IX0
     DO I = I1, I2, INC
        IP = IPIV( IX )
        IF( IP.NE.I ) THEN
           DO K = N32, N
              TEMP = A( I, K )
              A( I, K ) = A( IP, K )
              A( IP, K ) = TEMP
           enddo
        END IF
        IX = IX + INCX
     enddo
  END IF
!
  RETURN
!
!     End of SLASWP
!
END SUBROUTINE SLASWP


SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
  
  use precision
  
  implicit none
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
  INTEGER(kind=ik), intent(out) :: INFO
  INTEGER(kind=ik), intent(in)  :: LDA, M, N
!     ..
!     .. Array Arguments ..
  INTEGER(kind=ik), intent(out) :: IPIV(*)
  REAL(kind=rk), intent(inout)  :: A(LDA,*)
!     ..
!
!  Purpose
!  =======
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
  REAL(kind=rk), parameter :: ONE=1.0_rk, ZERO=0.0_rk
  !     ..
!     .. Local Scalars ..
  INTEGER(kind=ik)         :: J, JP
!     ..
!     .. External Functions ..
!      INTEGER            ISAMAX
!      EXTERNAL           ISAMAX
!     ..
!     .. External Subroutines ..
!      EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DGETF2', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( M.EQ.0 .OR. N.EQ.0 ) &
        RETURN
!
  DO J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
     JP = J - 1 + ISAMAX( M-J+1, A( J, J ), 1 )
     IPIV( J ) = JP
     IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
        IF( JP.NE.J ) &
             CALL SSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
        IF( J.LT.M ) &
              CALL SSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
     ELSE IF( INFO.EQ.0 ) THEN
!
        INFO = J
     END IF
     !
     IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
        CALL SGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, &
                   A( J+1, J+1 ), LDA )
     END IF
  enddo
  RETURN
!
!     End of DGETF2
!
END SUBROUTINE DGETF2

SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
  use precision
  
  implicit none

!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993 
!
!     .. Scalar Arguments ..
  INTEGER(kind=ik) :: INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
  INTEGER(kind=ik) :: IPIV( * )
  REAL(kind=rk) :: A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!     .. Parameters ..
  REAL(kind=rk), parameter :: ONE=1.0_rk
!     ..
!     .. Local Scalars ..
  INTEGER(kind=ik)         :: I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
!  EXTERNAL           SGEMM, SGETF2, SLASWP, STRSM, XERBLA
!     ..
!     .. External Functions ..
!  INTEGER            ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SGETRF', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
  !
  IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
!     Determine the block size for this environment.
!
  NB = ILAENV( 1, 'SGETRF', ' ', M, N, -1, -1 )
  IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
     CALL DGETF2( M, N, A, LDA, IPIV, INFO )
  ELSE
!
!        Use blocked code.
!
     DO J = 1, MIN( M, N ), NB
        JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
        CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
        IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + J - 1
        DO I = J, MIN( M, J+JB-1 )
           IPIV( I ) = J - 1 + IPIV( I )
        enddo
!
!           Apply interchanges to columns 1:J-1.
!
        CALL SLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
        IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
           CALL SLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, &
                        IPIV, 1 )
!
!              Compute block row of U.
!
           CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                        N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),&
                        LDA )
           IF( J+JB.LE.M ) THEN
!
              !                 Update trailing submatrix.
!
              CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                           N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,   &
                           A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), &
                          LDA )
           END IF
        END IF
     enddo
  END IF
  RETURN
!
!     End of SGETRF
!
END SUBROUTINE DGETRF


SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
  use precision

  implicit none

!     .. Scalar Arguments ..
  REAL(kind=rk), intent(in)    :: ALPHA, BETA
  INTEGER(kind=ik), intent(in) :: INCX, INCY, LDA, M, N
  CHARACTER, intent(in)        :: TRANS
!     ..
!     .. Array Arguments ..
  REAL(kind=rk), intent(in)    :: A(LDA,*), X(*)
  REAL(kind=rk), intent(inout) :: Y(*)
!     ..
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - REAL             array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - REAL            .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - REAL             array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!  =====================================================================
!
!     .. Parameters ..
  REAL(kind=rk), parameter :: ONE=1.0_rk, ZERO=0.0_rk
!     ..
!     .. Local Scalars ..
  REAL(kind=rk)            :: TEMP
  INTEGER(kind=ik)         :: I, INFO, IX, IY, J, JX, JY, KX, &
                              KY, LENX, LENY
!     ..
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
  INFO = 0
  IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
       .NOT.LSAME(TRANS,'C')) THEN
     INFO = 1
  ELSE IF (M.LT.0) THEN
     INFO = 2
  ELSE IF (N.LT.0) THEN
     INFO = 3
  ELSE IF (LDA.LT.MAX(1,M)) THEN
     INFO = 6
  ELSE IF (INCX.EQ.0) THEN
     INFO = 8
  ELSE IF (INCY.EQ.0) THEN
     INFO = 11
  END IF
  IF (INFO.NE.0) THEN
     CALL XERBLA('DGEMV ',INFO)
     RETURN
  END IF
!
!     Quick return if possible.
!
  IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
         ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
  IF (LSAME(TRANS,'N')) THEN
     LENX = N
     LENY = M
  ELSE
     LENX = M
     LENY = N
  END IF
  IF (INCX.GT.0) THEN
     KX = 1
  ELSE
     KX = 1 - (LENX-1)*INCX
  END IF
  IF (INCY.GT.0) THEN
     KY = 1
  ELSE
     KY = 1 - (LENY-1)*INCY
  END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
  IF (BETA.NE.ONE) THEN
     IF (INCY.EQ.1) THEN
        IF (BETA.EQ.ZERO) THEN
           DO I = 1,LENY
              Y(I) = ZERO
           enddo
        ELSE
           DO I = 1,LENY
              Y(I) = BETA*Y(I)
           enddo
        END IF
     ELSE
        IY = KY
        IF (BETA.EQ.ZERO) THEN
           DO I = 1,LENY
              Y(IY) = ZERO
              IY = IY + INCY
           enddo
        ELSE
           DO I = 1,LENY
              Y(IY) = BETA*Y(IY)
              IY = IY + INCY
           enddo
        END IF
     END IF
  END IF
  IF (ALPHA.EQ.ZERO) RETURN
  IF (LSAME(TRANS,'N')) THEN
     !
     !        Form  y := alpha*A*x + y.
!
     JX = KX
     IF (INCY.EQ.1) THEN
        DO J = 1,N
           IF (X(JX).NE.ZERO) THEN
              TEMP = ALPHA*X(JX)
              DO I = 1,M
                 Y(I) = Y(I) + TEMP*A(I,J)
              enddo
           END IF
           JX = JX + INCX
        enddo
     ELSE
        DO J = 1,N
           IF (X(JX).NE.ZERO) THEN
              TEMP = ALPHA*X(JX)
              IY = KY
              DO I = 1,M
                 Y(IY) = Y(IY) + TEMP*A(I,J)
                 IY = IY + INCY
              enddo
           END IF
           JX = JX + INCX
        enddo
     END IF
  ELSE
     !
!        Form  y := alpha*A'*x + y.
!
     JY = KY
     IF (INCX.EQ.1) THEN
        DO J = 1,N
           TEMP = ZERO
           DO I = 1,M
              TEMP = TEMP + A(I,J)*X(I)
           enddo
           Y(JY) = Y(JY) + ALPHA*TEMP
           JY = JY + INCY
        enddo
     ELSE
        DO J = 1,N
           TEMP = ZERO
           IX = KX
           DO I = 1,M
              TEMP = TEMP + A(I,J)*X(IX)
              IX = IX + INCX
           enddo
           Y(JY) = Y(JY) + ALPHA*TEMP
           JY = JY + INCY
        enddo
     END IF
  END IF
  !
  RETURN
!
  !     End of DGEMV .
  !
END SUBROUTINE DGEMV


end module lapack
