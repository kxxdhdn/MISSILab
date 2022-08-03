!******************************************************************************
!*
!*                    VARIOUS MATHEMATICAL DISTRIBUTIONS
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 09/2007
  !    - Updated 11/2012
  !    - 04/2016: remove ERF and ERFC (now intrinsic in F08).
  ! 
  ! 3) DESCRIPTION: implements various commonly used mathematical functions,
  !                 distributions, and operations to manipulate them.
  !==========================================================================

MODULE special_functions

  USE utilities, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: factorial_small, factorial, lngamma, igamma, igammac, ibeta
  PUBLIC :: expm1

  INTERFACE expm1
    MODULE PROCEDURE expm1_scl, expm1_1D
  END INTERFACE expm1

  INTERFACE lngamma
    MODULE PROCEDURE lngamma_scl, lngamma_1D
  END INTERFACE lngamma

  INTERFACE igamma
    MODULE PROCEDURE igamma_scl, igamma_1D
  END INTERFACE igamma

  INTERFACE igammac
    MODULE PROCEDURE igammac_scl, igammac_1D
  END INTERFACE igammac

  INTERFACE factorial
    MODULE PROCEDURE factorial_scl, factorial_1D
  END INTERFACE factorial

  INTERFACE betacf
    MODULE PROCEDURE betacf_scl, betacf_1D
  END INTERFACE betacf

  INTERFACE ibeta
    MODULE PROCEDURE ibeta_scl, ibeta_1D
  END INTERFACE ibeta


CONTAINS


  !==========================================================================
  ! EXP(X) - 1 = EXPM1(x)
  !
  !   Numerically accurate development of exp(x)-1 around 0, using the 
  ! expansion of Abramowitz & Stegun (4.2.45). It is suposed to be accurate
  ! at 2.E-10.
  !==========================================================================

  PURE FUNCTION expm1_1D (x)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x)) :: expm1_1D
 
    INTEGER, PARAMETER :: Ncoef = 7
    REAL(DP), PARAMETER :: xmax = 0.6931472_DP ! ln(2)
    REAL(DP), DIMENSION(Ncoef), PARAMETER :: &
      coef = [ 0.9999999995_DP, 0.4999999206_DP, 0.1666653019_DP, &
               0.0416573475_DP, 0.0083013598_DP, 0.0013298820_DP, &
               0.0001413161_DP ]

    !------------------------------------------------------------------------

    WHERE (ABS(x(:)) > xmax) 
      expm1_1D(:) = EXP(x(:)) - 1._DP
    ELSEWHERE
      expm1_1D(:) = coef(7)*x(:)**7 
      expm1_1D(:) = coef(6)*x(:)**6 + expm1_1D(:)
      expm1_1D(:) = coef(5)*x(:)**5 + expm1_1D(:)
      expm1_1D(:) = coef(4)*x(:)**4 + expm1_1D(:)
      expm1_1D(:) = coef(3)*x(:)**3 + expm1_1D(:)
      expm1_1D(:) = coef(2)*x(:)**2 + expm1_1D(:)
      expm1_1D(:) = coef(1)*x(:)**1 + expm1_1D(:)
    END WHERE

    !------------------------------------------------------------------------

  END FUNCTION expm1_1D

  PURE FUNCTION expm1_scl (x)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: expm1_scl
 
    !------------------------------------------------------------------------

    expm1_scl = MAXVAL(EXPM1_1D([x]))

    !------------------------------------------------------------------------

  END FUNCTION expm1_scl


  !==========================================================================
  ! n! = FACTORIAL_SMALL(n)
  !
  !   Determines the factorial of a small integer as an integer.
  !==========================================================================

  RECURSIVE FUNCTION factorial_small (n) RESULT(fact)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: fact

    !------------------------------------------------------------------------

    def: IF (n <= 1) THEN
      fact = 1
    ELSE
      fact = n*FACTORIAL_SMALL(n-1)
    END IF def

    !------------------------------------------------------------------------

  END FUNCTION factorial_small


  !==========================================================================
  ! n! = FACTORIAL(n)
  !
  !   Determines the factorial of an integer.
  !==========================================================================

  FUNCTION factorial_scl (n)

    USE utilities, ONLY: DP, cumprod, arth
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(DP) :: factorial_scl

    INTEGER, PARAMETER :: nmax = 32

    INTEGER, SAVE :: ntop = 0
    REAL(DP), DIMENSION(nmax), SAVE :: a

    !------------------------------------------------------------------------

    IF (n < ntop) THEN
      factorial_scl = a(n+1)
    ELSE IF (n < nmax) THEN
      ntop = nmax
      a(1) = 1._DP
      a(2:nmax) = CUMPROD(ARTH(1._DP,1._DP,nmax-1))
      factorial_scl = a(n+1)
    ELSE 
      factorial_scl = EXP(LNGAMMA(n+1._DP))
    END IF

    !------------------------------------------------------------------------

  END FUNCTION factorial_scl

    !------------------------------------------------------------------------

  FUNCTION factorial_1D (n)

    USE utilities, ONLY: DP, cumprod, arth
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN) :: n
    REAL(DP), DIMENSION(SIZE(n)) :: factorial_1D

    INTEGER, PARAMETER :: nmax = 32

    LOGICAL, DIMENSION(SIZE(n)) :: mask
    INTEGER, SAVE :: ntop = 0
    REAL(DP), DIMENSION(nmax), SAVE :: a

    !------------------------------------------------------------------------

    IF (ntop == 0) THEN
      ntop = nmax
      a(1) = 1._DP
      a(2:nmax) = CUMPROD(ARTH(1._DP,1._DP,nmax-1))
    END IF
    mask = (n >= nmax)
    factorial_1D = UNPACK(EXP(LNGAMMA(PACK(n,mask)+1.0_DP)),mask,0._DP)
    WHERE (.NOT. mask) factorial_1D = a(n+1)

    !------------------------------------------------------------------------

  END FUNCTION factorial_1D


  !==========================================================================
  ! y = LN(GAMMA(x))
  !
  !   Log of gamma function for x > 0.
  !==========================================================================

  ELEMENTAL FUNCTION lngamma_scl (xx)

    USE utilities, ONLY: DP, arth
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: xx
    REAL(DP) :: lngamma_scl

    REAL(DP), PARAMETER :: stp = 2.5066282746310005_DP
    REAL(DP), DIMENSION(6), PARAMETER :: &
      coef = (/ 76.18009172947146_DP,     -86.50532032941677_DP, &
                24.01409824083091_DP,     -1.231739572450155_DP, &
                0.1208650973866179E-2_DP, -0.5395239384953E-5_DP /)

    REAL(DP) :: tmp, x

    !------------------------------------------------------------------------

    x = xx
    tmp = x + 5.5_DP
    tmp = (x+0.5_DP)*LOG(tmp) - tmp
    lngamma_scl = tmp + LOG(stp*(1.000000000190015_DP &
                              +SUM(coef(:)/ARTH(x+1.0_DP,1.0_DP,SIZE(coef))))/x)

    !------------------------------------------------------------------------

  END FUNCTION lngamma_scl

    !------------------------------------------------------------------------

  PURE FUNCTION lngamma_1D (xx)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: xx
    REAL(DP), DIMENSION(SIZE(xx)) :: lngamma_1D

    REAL(DP), PARAMETER :: stp = 2.5066282746310005_DP
    REAL(DP), DIMENSION(6), PARAMETER :: &
      coef = (/ 76.18009172947146_DP,     -86.50532032941677_DP, &
                24.01409824083091_DP,     -1.231739572450155_DP, &
                0.1208650973866179E-2_DP, -0.5395239384953E-5_DP /)

    INTEGER :: i
    REAL(DP), DIMENSION(SIZE(xx)) :: ser, tmp, x, y

    !------------------------------------------------------------------------

    x = xx
    tmp = x + 5.5_DP
    tmp = (x+0.5_DP)*LOG(tmp) - tmp
    ser = 1.000000000190015_DP
    y = x
    DO i=1,SIZE(coef)
      y = y + 1.0_DP
      ser = ser + coef(i)/y
    END DO
    lngamma_1D = tmp + LOG(stp*ser/x)    

    !------------------------------------------------------------------------

  END FUNCTION lngamma_1D


  !==========================================================================
  ! IGAMMA(x) = 1/Gamma(a) * int(EXP(-t)*t^(a-1),{t=[0,x]}), a>0
  !
  !   Incomplete gamma function. We must have x >= 0 and a > 0. The preliminary
  ! functions GSER et GCF returns the two forms of the incomplete gamma function
  ! P(a,x) and Q(a,x) = 1 - P(a,x), with different numerical methods, depending
  ! on the parameter values.
  !==========================================================================

  ELEMENTAL FUNCTION gser_scl (a,x)

    USE utilities, ONLY: DP, epsDP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: a, x
    REAL(DP) :: gser_scl

    INTEGER, PARAMETER :: itmax = 100

    INTEGER :: n
    REAL(DP) :: ap, del, summ

    !------------------------------------------------------------------------

    IF (x == 0._DP) THEN
      gser_scl = 0._DP
      RETURN
    END IF
    ap = a
    summ = 1._DP / a
    del = summ
    DO n=1,itmax
      ap = ap + 1._DP
      del = del * x / ap
      summ = summ + del
      IF (ABS(del) < ABS(summ)*epsDP) EXIT
    END DO
    gser_scl = summ * EXP( -x + a*LOG(x) - LNGAMMA(a) )

    !------------------------------------------------------------------------

  END FUNCTION gser_scl

    !------------------------------------------------------------------------

  PURE FUNCTION gser_1D(a,x)
   
    USE utilities, ONLY: DP, epsDP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: a, x
    REAL(DP), DIMENSION(SIZE(a)) :: gser_1D

    INTEGER, PARAMETER :: itmax = 100

    INTEGER :: i, N
    REAL(DP), DIMENSION(SIZE(a)) :: ap, del, summ
    LOGICAL, DIMENSION(SIZE(a)) :: converged, zero

    !------------------------------------------------------------------------

    N = SIZE(x(:))
    zero(:) = ( x(:) == 0._DP )
    WHERE (zero) gser_1D(:) = 0._DP
    ap(:) = a(:)
    summ(:) = 1._DP / a(:)
    del(:) = summ(:)
    converged(:) = zero(:)
    DO i=1,itmax
      WHERE (.NOT. converged(:))
        ap(:) = ap(:) + 1._DP
        del(:) = del(:) * x(:) / ap(:)
        summ(:) = summ(:) + del(:)
        converged(:) = (ABS(del(:)) < ABS(summ(:))*epsDP)
      END WHERE
      IF (ALL(converged(:))) EXIT
    END DO
    WHERE (.NOT. zero(:)) &
      gser_1D(:) = summ(:) * EXP( -x(:) + a(:)*LOG(x(:)) - LNGAMMA(a(:)) )

    !------------------------------------------------------------------------

  END FUNCTION gser_1D

    !------------------------------------------------------------------------

  ELEMENTAL FUNCTION gcf_scl(a,x)

    USE utilities, ONLY: DP, epsDP, tinyDP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: a, x
    REAL(DP) :: gcf_scl

    INTEGER, PARAMETER :: itmax = 100
    REAL(DP), PARAMETER :: fpmin = tinyDP / epsDP

    INTEGER :: i
    REAL(DP) :: an, b, c, d, del, h

    !------------------------------------------------------------------------

    IF (x == 0._DP) THEN
      gcf_scl = 1._DP
      RETURN
    ENDIF
    b = x + 1._DP - a 
    c = 1._DP / fpmin
    d = 1._DP / b
    h = d
    DO i=1,itmax
      an = -i * (i-a)
      b = b + 2._DP
      d = an*d + b
      IF (ABS(d) < fpmin) d = fpmin
      c = b + an/c
      IF (ABS(c) < fpmin) c = fpmin
      d = 1._DP / d
      del = d * c
      h = h * del
      IF (ABS(del-1._DP) <= epsDP) EXIT
    END DO
    gcf_scl = EXP( -x + a*LOG(x) - LNGAMMA(a) ) * h

    !------------------------------------------------------------------------

  END FUNCTION gcf_scl

    !------------------------------------------------------------------------

  PURE FUNCTION gcf_1D(a,x)

    USE utilities, ONLY: DP, epsDP, tinyDP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: a, x
    REAL(DP), DIMENSION(SIZE(a)) :: gcf_1D

    INTEGER, PARAMETER :: itmax = 100
    REAL(DP), PARAMETER :: fpmin = tinyDP / epsDP

    INTEGER :: i, N
    REAL(DP), DIMENSION(SIZE(a)) :: an, b, c, d, del, h
    LOGICAL, DIMENSION(SIZE(a)) :: converged, zero

    !------------------------------------------------------------------------

    N = SIZE(x(:))
    zero(:) = ( x(:) == 0._DP )
    WHERE (zero(:))
      gcf_1D(:) = 1._DP
    ELSEWHERE
      b(:) = x(:) + 1._DP - a(:)
      c(:) = 1._DP / fpmin
      d(:) = 1._DP / b(:)
      h(:) = d(:)
    END WHERE
    converged(:) = zero(:)
    DO i=1,itmax
      WHERE (.NOT. converged(:))
        an(:) = -i * (i-a(:))
        b(:) = b(:) + 2._DP
        d(:) = an(:)*d(:) + b(:)
        d(:) = MERGE( fpmin, d(:), ABS(d(:))<fpmin )
        c(:) = b(:) + an(:)/c(:)
        c(:) = MERGE( fpmin, c(:), ABS(c(:))<fpmin )
        d(:) = 1._DP / d(:)
        del(:) = d(:) * c(:)
        h(:) = h(:) * del(:)
        converged(:) = (ABS(del(:)-1._DP)<=epsDP)
      END WHERE
      IF (ALL(converged(:))) EXIT
    END DO
    WHERE (.NOT. zero(:)) &
      gcf_1D(:) = EXP( -x(:) + a(:)*LOG(x(:)) - LNGAMMA(a(:)) ) * h(:)

    !------------------------------------------------------------------------

  END FUNCTION gcf_1D

    !------------------------------------------------------------------------

  ELEMENTAL FUNCTION igamma_scl (a,x)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: a, x
    REAL(DP) :: igamma_scl

    !------------------------------------------------------------------------

    IF ( x < a+1._DP ) THEN
      igamma_scl = GSER_SCL(a,x)
    ELSE
      igamma_scl = 1._DP - GCF_SCL(a,x)
    END IF

    !------------------------------------------------------------------------

  END FUNCTION igamma_scl

    !------------------------------------------------------------------------

  PURE FUNCTION igamma_1D (a,x)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: a, x
    REAL(DP), DIMENSION(SIZE(x)) :: igamma_1D

    INTEGER :: N
    LOGICAL, DIMENSION(SIZE(x)) :: mask

    !------------------------------------------------------------------------

    N = SIZE(x(:))
    mask(:) = ( x(:) < a(:)+1._DP )
    igamma_1D(:) = MERGE( GSER_1D(a(:),MERGE(x(:),0._DP,mask)), &
                          1._DP-GCF_1D(a(:),MERGE(x(:),0._DP,.NOT. mask(:))), &
                          mask(:) )

    !------------------------------------------------------------------------

  END FUNCTION igamma_1D


  !==========================================================================
  ! IGAMMAC(x) = 1/Gamma(a) * int(EXP(-t)*t^(a-1),{t=[x,infty]}), a>0
  !
  !   Complementary incomplete gamma function, Q(a,x) = 1 - P(a,x). We must 
  ! have x >= 0 and a > 0. 
  !==========================================================================

  ELEMENTAL FUNCTION igammac_scl (a,x)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: a, x
    REAL(DP) :: igammac_scl

    !------------------------------------------------------------------------

    IF ( x < a+1._DP) THEN
      igammac_scl = 1._DP - GSER_SCL(a,x) 
    ELSE
      igammac_scl = GCF_SCL(a,x)
    END IF

    !------------------------------------------------------------------------

  END FUNCTION igammac_scl

    !------------------------------------------------------------------------

  PURE FUNCTION igammac_1D(a,x)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: a, x
    REAL(DP), DIMENSION(SIZE(a)) :: igammac_1D
    LOGICAL, DIMENSION(SIZE(x)) :: mask
    INTEGER :: N

    !------------------------------------------------------------------------

    N = SIZE(x(:))
    mask(:) = ( x(:) < a(:)+1._DP )
    igammac_1D = MERGE( 1._DP - GSER_1D(a(:),MERGE(x(:),0._DP,mask(:))), &
                        GCF_1D(a(:),MERGE(x(:),0._DP,.NOT. mask(:))), &
                        mask(:) )

    !------------------------------------------------------------------------

  END FUNCTION igammac_1D


  !==========================================================================
  ! y[N] = BETACF(a[N],b[N],x[N])
  !
  !   Evaluates continued fraction for incomplete beta function by modified 
  ! Lenz's method.
  !==========================================================================

  ELEMENTAL FUNCTION betacf_scl (a,b,x)

    USE utilities, ONLY: DP, tinyDP, epsDP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: a, b, x
    REAL(DP) :: betacf_scl
  
    INTEGER, PARAMETER :: maxit = 100
    REAL(DP), PARAMETER :: fpmin = tinyDP/epsDP

    INTEGER :: m, m2
    REAL(DP) :: aa, c, d, del, h, qab, qam, qap
        
    !------------------------------------------------------------------------

    qab = a + b
    qap = a + 1._DP
    qam = a - 1._DP
    c = 1._DP
    d = 1._DP - qab*x/qap
    IF (ABS(d) < fpmin) d = fpmin
    d = 1._DP/d
    h = d
    DO m=1,maxit
      m2 = 2*m
      aa = m * (b-m) * x / ((qam+m2)*(a+m2))
      d = 1._DP + aa*d
      IF (ABS(d) < fpmin) d = fpmin
      c = 1._DP + aa/c
      IF (ABS(c) < fpmin) c = fpmin
      d = 1._DP/d
      h = h*d*c
      aa = - (a+m) * (qab+m) * x / ((a+m2)*(qap+m2))
      d = 1._DP + aa*d
      IF (ABS(d) < fpmin) d = fpmin
      c = 1._DP + aa/c
      IF (ABS(c) < fpmin) c = fpmin
      d = 1._DP / d
      del = d*c
      h = h*del
      IF (ABS(del-1._DP) <= epsDP) EXIT
    END DO

    betacf_scl = h

    !------------------------------------------------------------------------

  END FUNCTION betacf_scl

    !------------------------------------------------------------------------

  PURE FUNCTION betacf_1D (a,b,x)

    USE utilities, ONLY: DP, tinyDP, epsDP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: a, b, x
    REAL(DP), DIMENSION(SIZE(x)) :: betacf_1D

    INTEGER, PARAMETER :: maxit = 100
    REAL(DP), PARAMETER :: fpmin = tinyDP/epsDP

    INTEGER :: m
    INTEGER, DIMENSION(SIZE(x)) :: m2
    REAL(DP), DIMENSION(SIZE(x)) :: aa, c, d, del, h, qab, qam, qap
    LOGICAL, DIMENSION(SIZE(x)) :: converged

    !------------------------------------------------------------------------

    m = SIZE(x)
    qab = a + b
    qap = a + 1._DP
    qam = a - 1._DP
    c = 1._DP
    d = 1._DP - qab*x/qap
    WHERE (ABS(d) < fpmin) d = fpmin
    d = 1._DP/d
    h = d
    converged = .FALSE.
    DO m=1,maxit
      WHERE (.NOT. converged)
        m2 = 2*m
        aa = m * (b-m) * x / ((qam+m2)*(a+m2))
        d = 1._DP + aa*d
        d = MERGE(fpmin,d,ABS(d)<fpmin)
        c = 1._DP + aa/c
        c = MERGE(fpmin,c,ABS(c)<fpmin)
        d = 1._DP/d
        h = h*d*c
        aa = - (a+m) * (qab+m) * x / ((a+m2)*(qap+m2))
        d = 1._DP + aa*d
        d = MERGE(fpmin,d,ABS(d)<fpmin)
        c = 1._DP + aa/c
        c = MERGE(fpmin,c,ABS(c)<fpmin)
        d = 1._DP/d
        del = d*c
        h = h*del
        converged = (ABS(del-1._DP) <= epsDP)
      END WHERE
      IF (ALL(converged)) EXIT
    END DO

    betacf_1D = h

    !------------------------------------------------------------------------

  END FUNCTION betacf_1D


  !==========================================================================
  ! y[N] = IBETA(a[N],b[N],x[N])
  !
  !   Incomplete beta function.
  !==========================================================================

  ELEMENTAL FUNCTION ibeta_scl (a,b,x)

    USE utilities, ONLY: DP
    IMPLICIT NONE
  
    REAL(DP), INTENT(IN) :: a, b, x
    REAL(DP) :: ibeta_scl
    REAL(DP) :: bt

    !------------------------------------------------------------------------

    IF (x == 0._DP .OR. x == 1._DP) THEN
      bt = 0._DP
    ELSE
      bt = EXP( LNGAMMA(a+b) - LNGAMMA(a) - LNGAMMA(b) &
              + a*LOG(x) + b*LOG(1._DP-x) )
    END IF
    IF (x < (a+1._DP)/(a+b+2._DP)) THEN
      ibeta_scl = bt * BETACF(a,b,x) / a
    ELSE
      ibeta_scl = 1._DP - bt * BETACF(b,a,1._DP-x) / b
    END IF

    !------------------------------------------------------------------------

  END FUNCTION ibeta_scl

    !------------------------------------------------------------------------

  PURE FUNCTION ibeta_1D (a,b,x)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: a, b, x
    REAL(DP), DIMENSION(SIZE(a)) :: ibeta_1D
  
    INTEGER :: ndum
    REAL(DP), DIMENSION(SIZE(a)) :: bt
    LOGICAL, DIMENSION(SIZE(a)) :: mask

    !------------------------------------------------------------------------

    ndum = SIZE(x)
    WHERE (x == 0._DP .OR. x == 1._DP)
      bt = 0._DP
    ELSEWHERE
      bt = EXP( LNGAMMA(a+b) - LNGAMMA(a) - LNGAMMA(b) &
              + a*LOG(x) + b*LOG(1._DP-x) )
    END WHERE
    mask = (x < (a+1._DP)/(a+b+2._DP))
    ibeta_1D = bt * BETACF(MERGE(a,b,mask),MERGE(b,a,mask), &
                            MERGE(x,1._DP-x,mask)) / MERGE(a,b,mask)
    WHERE (.NOT. mask) ibeta_1D = 1._DP - ibeta_1D

    !------------------------------------------------------------------------

  END FUNCTION ibeta_1D


END MODULE special_functions
