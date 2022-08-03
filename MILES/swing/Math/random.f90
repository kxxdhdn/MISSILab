!******************************************************************************
!*
!*                 VARIOUS DISTRIBUTIONS OF RANDOM VARIABLES
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 10/2007 
  !    - 03/2015: add multivariate distributions
  !    - 03/2015: add GENERATE_NEWSEED
  !    - 04/2016: add the skew-normal generator
  ! 
  ! 3) DESCRIPTION: provides random generators of various well-known
  !                 probability distributions.
  !==========================================================================

MODULE random

  USE utilities, ONLY: DP
  USE arrays, ONLY:
  USE inout, ONLY:
  USE linear_system, ONLY:
  USE interpolation, ONLY:
  USE adaptative_grid, ONLY:
  PRIVATE

  PUBLIC :: rand_exp, rand_norm, rand_skewnorm, rand_splitnorm, rand_student
  PUBLIC :: rand_poisson, rand_general, rand_multinorm, rand_multistudent
  PUBLIC :: generate_newseed

  INTERFACE rand_exp
    MODULE PROCEDURE rand_exp_v, rand_exp_s
  END INTERFACE rand_exp

  INTERFACE rand_norm
    MODULE PROCEDURE rand_norm_m, rand_norm_v, rand_norm_s
  END INTERFACE rand_norm

  INTERFACE rand_skewnorm
    MODULE PROCEDURE rand_skewnorm_m, rand_skewnorm_v, rand_skewnorm_s
  END INTERFACE rand_skewnorm

  INTERFACE rand_splitnorm
    MODULE PROCEDURE rand_splitnorm_m, rand_splitnorm_v, rand_splitnorm_s
  END INTERFACE rand_splitnorm

  INTERFACE rand_student
    MODULE PROCEDURE rand_student_m, rand_student_v, rand_student_s
  END INTERFACE rand_student

  INTERFACE rand_multinorm
    MODULE PROCEDURE rand_multinorm_v, rand_multinorm_s
  END INTERFACE rand_multinorm

  INTERFACE rand_multistudent
    MODULE PROCEDURE rand_multistudent_v, rand_multistudent_s
  END INTERFACE rand_multistudent

  INTERFACE rand_general
    MODULE PROCEDURE rand_general_vec, rand_general_scl
  END INTERFACE rand_general
  

CONTAINS


  !==========================================================================
  ! CALL GENERATE_NEWSEED(fileIN,fileOUT)
  ! 
  !   IF FILEIN and/or FILEOUT are set, then the sequence of integers 
  ! constituing the seed is read/written. To ensure the seed is different at
  ! each run, we initialize it with date and time.
  !==========================================================================

  SUBROUTINE generate_newseed (fileIN,fileOUT)

    USE utilities, ONLY: verbatim, trimlr
    USE inout, ONLY: lenpath, ascext
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN), OPTIONAL :: fileIN
    CHARACTER(lenpath), INTENT(OUT), OPTIONAL :: fileOUT

    INTEGER, PARAMETER :: unitIN = 1, unitOUT = 2

    INTEGER :: i, Nseed
    INTEGER, DIMENSION(8) :: values
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    !-----------------------------------------------------------------------

    ! Get the new seed
    genseed: IF (PRESENT(fileIN)) THEN 
      OPEN(unitIN,FILE=fileIN,ACTION="READ")
        READ(unitIN,"(I20)") Nseed
        ALLOCATE (seed(Nseed))
        DO i=1,Nseed
          READ(unitIN,"(I20)") seed(i)
        END DO
      CLOSE(unitIN)
      IF (verbatim) PRINT*, " - File "//TRIMLR(fileIN)//" has been written."
    ELSE 
      CALL RANDOM_SEED(SIZE=Nseed)
      ALLOCATE (seed(Nseed))
      CALL RANDOM_SEED(GET=seed)
      CALL DATE_AND_TIME(VALUES=values)
      seed(:) = seed(:)*SUM(values)
    END IF genseed

    ! Writes the seed to a file
    IF (PRESENT(fileOUT)) THEN
      OPEN(unitOUT,FILE=fileOUT,STATUS="REPLACE",ACTION="WRITE")
        WRITE(unitOUT,"(I20)") Nseed
        DO i=1,Nseed
           WRITE(unitOUT,"(I20)") seed(i)
        END DO
      CLOSE(unitOUT)
    END IF

    ! Initialize the random sequence
    CALL RANDOM_SEED(PUT=seed)
    DEALLOCATE (seed)

    !-----------------------------------------------------------------------

  END SUBROUTINE generate_newseed


  !==========================================================================
  ! x = RAND_EXP (N,tau)
  !
  !   Returns an array[N] of random variables following an exponential law 
  ! with a timescale of TAU. Taken from http://users.bigpond.net.au/amiller/ .
  !==========================================================================

  FUNCTION rand_exp_v (N,tau)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N
    REAL(DP), INTENT(IN), OPTIONAL :: tau
    REAL(DP), DIMENSION(N) :: rand_exp_v

    INTEGER :: i
    REAL(DP), DIMENSION(N) :: r

    !-----------------------------------------------------------------------

    varN: DO i=1,N
      randgen: DO 
        CALL RANDOM_NUMBER (r(i))
        IF (r(i) > 0._DP) EXIT
      END DO randgen
    END DO varN
      
    ! p(t) = EXP(-t/tau)
    rand_exp_v(:) = - LOG(r(:)) 
    IF (PRESENT(tau)) rand_exp_v(:) = rand_exp_v(:) * tau

    !-----------------------------------------------------------------------

  END FUNCTION rand_exp_v

  ! Interface shell
  FUNCTION rand_exp_s(tau)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN), OPTIONAL :: tau
    REAL(DP) :: rand_exp_s
    rand_exp_s = MAXVAL(RAND_EXP_V(1,tau))
  END FUNCTION rand_exp_s


  !==========================================================================
  ! x = RAND_NORM (N,M,center,sigma)
  !
  !   Returns an array[N] of random variables following a normal distribution
  ! of mean CENTER and standard deviation SIGMA. Taken from 
  ! http://people.sc.fsu.edu/~jburkardt/f_src/normal/normal.f90, tested and 
  ! vectorized.
  !==========================================================================

  FUNCTION rand_norm_m (N,M,center,sigma)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    USE arrays, ONLY: reallocate
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N, M
    REAL(DP), INTENT(IN), OPTIONAL :: center, sigma
    REAL(DP), DIMENSION(N,M) :: rand_norm_m

    REAL(DP), DIMENSION(N,M) :: r1, r2, x

    !-----------------------------------------------------------------------

    CALL RANDOM_NUMBER(r1(:,:))
    CALL RANDOM_NUMBER(r2(:,:))
    x(:,:) = SQRT( - 2._DP * LOG(r1(:,:)) ) * COS( 2._DP * pi * r2(:,:) )
    rand_norm_m(:,:) = x(:,:)
    IF (PRESENT(sigma)) rand_norm_m(:,:) = rand_norm_m(:,:) * sigma
    IF (PRESENT(center)) rand_norm_m(:,:) = rand_norm_m(:,:) + center 

    !-----------------------------------------------------------------------

  END FUNCTION rand_norm_m

  ! Interface shells
  FUNCTION rand_norm_v (N,center,sigma)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL(DP), INTENT(IN), OPTIONAL :: center, sigma
    REAL(DP), DIMENSION(N) :: rand_norm_v
    rand_norm_v(:) = RESHAPE(RAND_NORM_M(N,1,CENTER=center,SIGMA=sigma),[N])
  END FUNCTION rand_norm_v
  FUNCTION rand_norm_s (center,sigma)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN), OPTIONAL :: center, sigma
    REAL(DP) :: rand_norm_s
    rand_norm_s = MAXVAL(RAND_NORM_M(1,1,CENTER=center,SIGMA=sigma))
  END FUNCTION rand_norm_s


  !==========================================================================
  ! x = RAND_SKEWNORM (N,M,ksi,omega,alpha)
  !
  !   Returns an array[N] of random variables following a skew-normal
  ! distribution of position parameter KSI, scale parameter OMEGA, and shape
  ! parameter ALPHA. Extracted from the MATLAB library.
  !==========================================================================

  FUNCTION rand_skewnorm_m (N,M,ksi,omega,alpha)

    USE utilities, ONLY: DP
    USE arrays, ONLY: reallocate
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N, M
    REAL(DP), INTENT(IN), OPTIONAL :: ksi, omega, alpha
    REAL(DP), DIMENSION(N,M) :: rand_skewnorm_m

    REAL(DP), DIMENSION(N,M) :: r1, r2, x
    REAL(DP) :: delta, ksi0, omega0, alpha0

    !-----------------------------------------------------------------------

    IF (PRESENT(ksi)) THEN ; ksi0 = ksi ; ELSE ; ksi0 = 0._DP ; END IF
    IF (PRESENT(omega)) THEN ; omega0 = omega ; ELSE ; omega0 = 1._DP ; END IF
    IF (PRESENT(alpha)) THEN ; alpha0 = alpha ; ELSE ; alpha0 = 0._DP ; END IF
    r1(:,:) = RAND_NORM_M(N,M)
    r2(:,:) = RAND_NORM_M(N,M)
    delta = alpha0 / SQRT(1+alpha0**2)
    x(:,:) = delta * r1(:,:) + SQRT(1-delta**2)*r2(:,:)
    rand_skewnorm_m(:,:) = SIGN(1._DP,r1(:,:))*x(:,:)*omega0 + ksi0
    
    !-----------------------------------------------------------------------

  END FUNCTION rand_skewnorm_m

  ! Interface shells
  FUNCTION rand_skewnorm_v (N,ksi,omega,alpha)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL(DP), INTENT(IN), OPTIONAL :: ksi, omega, alpha
    REAL(DP), DIMENSION(N) :: rand_skewnorm_v
    rand_skewnorm_v(:) &
      = RESHAPE(RAND_SKEWNORM_M(N,1,KSI=ksi,OMEGA=omega,ALPHA=alpha),[N])
  END FUNCTION rand_skewnorm_v
  FUNCTION rand_skewnorm_s (ksi,omega,alpha)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN), OPTIONAL :: ksi, omega, alpha
    REAL(DP) :: rand_skewnorm_s
    rand_skewnorm_s &
      = MAXVAL(RAND_SKEWNORM_M(1,1,KSI=ksi,OMEGA=omega,ALPHA=alpha))
  END FUNCTION rand_skewnorm_s


  !==========================================================================
  ! x = RAND_SPLITNORM (N,M,mu,lambda,tau)
  !
  !   Returns an array[N] of random variables following a split-normal
  ! distribution of position parameter MU, scale parameter LAMBDA, and shape
  ! parameter TAU. 
  !==========================================================================

  FUNCTION rand_splitnorm_m (N,M,mu,lambda,tau)

    USE utilities, ONLY: DP
    USE arrays, ONLY: reallocate
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N, M
    REAL(DP), INTENT(IN), OPTIONAL :: mu, lambda, tau
    REAL(DP), DIMENSION(N,M) :: rand_splitnorm_m

    REAL(DP), DIMENSION(N,M) :: r1, r2
    REAL(DP) :: mu0, lambda0, tau0

    !-----------------------------------------------------------------------

    IF (PRESENT(mu)) THEN ; mu0 = mu ; ELSE ; mu0 = 0._DP ; END IF
      IF (PRESENT(lambda)) THEN ; lambda0 = lambda
                           ELSE ; lambda0 = 1._DP ; END IF
    IF (PRESENT(tau)) THEN ; tau0 = tau ; ELSE ; tau0 = 0._DP ; END IF
    r1(:,:) = RAND_NORM_M(N,M)
    CALL RANDOM_NUMBER(r2(:,:))
    WHERE (r2(:,:) < 1._DP/(1._DP+tau))
      rand_splitnorm_m(:,:) = - ABS(r1(:,:)) * lambda + mu
    ELSEWHERE
      rand_splitnorm_m(:,:) = ABS(r1(:,:)) * lambda * tau + mu
    END WHERE
    
    !-----------------------------------------------------------------------

  END FUNCTION rand_splitnorm_m

  ! Interface shells
  FUNCTION rand_splitnorm_v (N,mu,lambda,tau)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL(DP), INTENT(IN), OPTIONAL :: mu, lambda, tau
    REAL(DP), DIMENSION(N) :: rand_splitnorm_v
    rand_splitnorm_v(:) &
      = RESHAPE(RAND_SPLITNORM_M(N,1,MU=mu,LAMBDA=lambda,TAU=tau),[N])
  END FUNCTION rand_splitnorm_v
  FUNCTION rand_splitnorm_s (mu,lambda,tau)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN), OPTIONAL :: mu, lambda, tau
    REAL(DP) :: rand_splitnorm_s
    rand_splitnorm_s = MAXVAL(RAND_SPLITNORM_M(1,1,MU=mu,LAMBDA=lambda,TAU=tau))
  END FUNCTION rand_splitnorm_s


  !==========================================================================
  ! x = RAND_MULTINORM(N,covar[M,M],mu)
  !
  !   Returns an array[M,N] of random variables following a multivariate normal
  ! distribution with covariance matrix COVAR. 
  !==========================================================================

  FUNCTION rand_multinorm_v (N,covar,mu) 

    USE utilities, ONLY: DP, strike
    USE linear_system, ONLY: cholesky_decomp
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: covar
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: mu
    REAL(DP), DIMENSION(SIZE(covar,1),N) :: rand_multinorm_v

    INTEGER :: i, M
    REAL(DP), DIMENSION(SIZE(covar,1),N) :: theta
    REAL(DP), DIMENSION(SIZE(covar,1),SIZE(covar,2)) :: L

    !-----------------------------------------------------------------------
 
    ! Covariance size
    M = SIZE(covar(:,:),1)
    IF (SIZE(covar(:,:),2) /= M) &
      CALL STRIKE("RAN_MULTINORM","wrong size for COVAR")
    
    ! Draw independent normal variables
    theta(:,:) = RAND_NORM(M,N)

    ! Perform Cholesky decompsition
    L(:,:) = CHOLESKY_DECOMP(covar(:,:)) 
    rand_multinorm_v(:,:) = MATMUL(L(:,:),theta(:,:))
    IF (PRESENT(mu)) THEN
      FORALL (i=1:N) rand_multinorm_v(:,i) = rand_multinorm_v(:,i) + mu(:)
    END IF

    !-----------------------------------------------------------------------
 
  END FUNCTION rand_multinorm_v

  ! Interface shell
  FUNCTION rand_multinorm_s (covar,mu) 
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: covar
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: mu
    REAL(DP), DIMENSION(SIZE(covar,1)) :: rand_multinorm_s
    rand_multinorm_s(:) = RESHAPE(RAND_MULTINORM_V(1,covar,mu),[SIZE(covar,1)])
  END FUNCTION rand_multinorm_s


  !==========================================================================
  ! x = RAND_STUDENT (N,M,Df,center,sigma)
  !
  !   Returns an array[N] of random variables following a Student's t 
  ! distribution of mean CENTER and standard deviation SIGMA, with Df degrees
  ! of freedom. Taken from http://www.netlib.org/random/.
  !==========================================================================

  ! Scalar function for a reduced random variable
  !----------------------------------------------
  FUNCTION rand_student_m(N,M,Df,center,sigma)

    USE utilities, ONLY: DP, strike, NaN
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N, M, Df
    REAL(DP), INTENT(IN), OPTIONAL :: center, sigma
    REAL(DP), DIMENSION(N,M) :: rand_student_m

    REAL(DP), PARAMETER :: zero = 0.0_DP, quart = 0.25_DP, half = 0.5_DP
    REAL(DP), PARAMETER :: one = 1.0_DP, two = 2.0_DP, three = 3.0_DP 
    REAL(DP), PARAMETER :: four = 4.0_DP, five = 5.0_DP, sixteen = 16.0_DP

    INTEGER :: mm, i, j
    REAL(DP), SAVE :: s, c, a, f, g
    REAL(DP) :: r, x, v

    !-----------------------------------------------------------------------

    IF (Df < 1) THEN
      rand_student_m(:,:) = NaN()
      RETURN
    END IF
  
    mm = 0
    DO i=1,N
      DO j=1,M
        IF (Df /= mm) THEN ! Initialization, if necessary
          s = Df
          c = - quart * (s + one)
          a = four / (one + one/s)**c
          f = sixteen / a
          IF (Df > 1) THEN
            g = s - one
            g = ((s + one)/g)**c * SQRT((s+s)/g)
          ELSE
            g = one
          END IF
          mm = Df
        END IF

        DO
          CALL RANDOM_NUMBER(r)
          IF (r <= zero) CYCLE
          CALL RANDOM_NUMBER(v)
          x = (two*v - one)*g/r
          v = x*x
          IF (v > five - a*r) THEN
            IF (Df >= 1 .AND. r*(v + three) > f) CYCLE
            IF (r > (one + v/s)**c) CYCLE
          END IF
          EXIT
        END DO
        rand_student_m(i,j) = x
      END DO
    END DO

    IF (PRESENT(sigma)) rand_student_m(:,:) = rand_student_m(:,:) * sigma
    IF (PRESENT(center)) rand_student_m(:,:) = rand_student_m(:,:) + center

    !-----------------------------------------------------------------------

  END FUNCTION rand_student_m

  ! Interface shell
  FUNCTION rand_student_v(N,Df,center,sigma)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, Df
    REAL(DP), INTENT(IN), OPTIONAL :: center, sigma
    REAL(DP), DIMENSION(N) :: rand_student_v
    rand_student_v(:) &
      = RESHAPE(RAND_STUDENT_M(N,1,Df,CENTER=center,SIGMA=sigma),[N])
  END FUNCTION rand_student_v
  FUNCTION rand_student_s(Df,center,sigma)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Df
    REAL(DP), INTENT(IN), OPTIONAL :: center, sigma
    REAL(DP) :: rand_student_s
    rand_student_s = MAXVAL(RAND_STUDENT_M(1,1,Df,CENTER=center,SIGMA=sigma))
  END FUNCTION rand_student_s


  !==========================================================================
  ! x = RAND_MULTISTUDENT (N,Df,covar[M,M],mu)
  !
  !   Returns an array[M,N] of random variables following a multivariate t
  ! distribution with covariance matrix COVAR. 
  !==========================================================================

  FUNCTION rand_multistudent_v (N,Df,covar,mu) 

    USE utilities, ONLY: DP, strike
    USE linear_system, ONLY: cholesky_decomp
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N, Df
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: covar
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: mu
    REAL(DP), DIMENSION(SIZE(covar,1),N) :: rand_multistudent_v

    INTEGER :: i, M
    REAL(DP), DIMENSION(SIZE(covar,1),N) :: theta
    REAL(DP), DIMENSION(SIZE(covar,1),SIZE(covar,2)) :: L

    !-----------------------------------------------------------------------
 
    ! Covariance size
    M = SIZE(covar(:,:),1)
    IF (SIZE(covar(:,:),2) /= M) &
      CALL STRIKE("RAN_MULTISTUDENT","wrong size for COVAR")
    
    ! Draw independent normal variables
    theta(:,:) = RAND_STUDENT(M,N,Df)

    ! Perform Cholesky decompsition
    L(:,:) = CHOLESKY_DECOMP(covar(:,:)) 
    rand_multistudent_v(:,:) = MATMUL(L(:,:),theta(:,:))
    IF (PRESENT(mu)) THEN
      FORALL (i=1:N) rand_multistudent_v(:,i) = rand_multistudent_v(:,i) + mu(:)
    END IF

    !-----------------------------------------------------------------------
 
  END FUNCTION rand_multistudent_v

  ! Interface shell
  FUNCTION rand_multistudent_s (Df,covar,mu) 
    USE utilities, ONLY: DP
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Df
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: covar
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: mu
    REAL(DP), DIMENSION(SIZE(covar,1)) :: rand_multistudent_s
    rand_multistudent_s(:) = RESHAPE(RAND_MULTISTUDENT_V(1,Df,covar,mu), &
                                     [SIZE(covar,1)])
  END FUNCTION rand_multistudent_s


  !==========================================================================
  ! x = GEN_POISSON (mu,first)
  !
  !   Returns a random poissonian deviate of average MU. FIRST for initiating
  ! the serie. Taken from http://users.bigpond.net.au/amiller/.
  !==========================================================================

  FUNCTION gen_poisson (mu,first) RESULT(ival)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: mu
    LOGICAL, INTENT(IN) :: first
    INTEGER :: ival

    REAL(DP), PARAMETER :: a0 = -0.5_DP, a1 = 0.3333333_DP, a2 = -0.2500068_DP
    REAL(DP), PARAMETER :: a3 = 0.2000118_DP, a4 = -0.1661269_DP
    REAL(DP), PARAMETER :: a5 = 0.1421878_DP, a6 = -0.1384794_DP
    REAL(DP), PARAMETER :: a7 = 0.1250060_DP
    REAL(DP), DIMENSION(10), PARAMETER :: &
      fact(10) = [ 1._DP, 1._DP, 2._DP, 6._DP, 24._DP, 120._DP, 720._DP, &
                   5040._DP, 40320._DP, 362880._DP ]

    INTEGER :: j, k, kflag
    INTEGER, SAVE :: l, m
    REAL(DP) :: b1, b2, c, c0, c1, c2, c3, del, difmuk, e, fk, fx, fy, g
    REAL(DP) :: omega, px, py, t, u, v, x, xx
    REAL(DP), SAVE :: s, d, p, q, p0
    LOGICAL, SAVE :: full_init
    REAL(DP), SAVE :: pp(35)

    !-----------------------------------------------------------------------

    fk = 0._DP
    difmuk = 0._DP
    ival = 0

    valofmu: IF (mu > 10.0) THEN
    ! CASE  A. (Recalculation OF S, D, L if MU has changed)
    !--------

      IF (first) THEN
        s = SQRT(mu)
        d = 6._DP*mu*mu
        ! The Poisson probabilities Pk exceed the discrete normal
        ! probabilities Fk whenever K >= M(MU). L=IFIX(MU-1.1484)
        ! is an upper bound to M(MU) for all MU >= 10 .
        l = INT(mu - 1.1484_DP)
        full_init = .False.
      END IF

      ! STEP N. Normal sample - random_normal() for standard normal deviate
      g = mu + s*MAXVAL(rand_norm(1))
      IF (g > 0.0) THEN
        ival = INT(g)
        ! STEP I. Immediate acceptance if IVAL is large enough
        IF (ival>=l) RETURN
        ! STEP S. Squeeze acceptance - Sample U
        fk = ival
        difmuk = mu - fk
        CALL RANDOM_NUMBER(u)
        IF (d*u >= difmuk*difmuk*difmuk) RETURN
      END IF

      ! STEP P. Preparations for steps Q and H.
      !         (Recpitulations of parameters if necessary)
      !         .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
      !         The quantities B1, B2, C3, C2, C1, C0 are for the Hermite
      !         Aproximations to the discrete normal probabilities Fk.
      !         C=.1069/MU guarantees majorization by the 'Hat'-function.
      IF (.NOT. full_init) THEN
        omega = 0.3989423_DP / s
        b1 = 0.4166667E-1/mu
        b2 = 0.3*b1*b1
        c3 = 0.1428571*b1*b2
        c2 = b2 - 15.*c3
        c1 = b1 - 6.*b2 + 45.*c3
        c0 = 1. - b1 + 3.*b2 - 15.*c3
        c = 0.1069/mu
        full_init = .True.
      ELSE 
        omega = 0._DP
        c3 = 0._DP
        c2 = 0._DP
        c1 = 0._DP
        c0 = 0._DP
        c = 0._DP
      END IF

      IF (g < 0.0) GOTO 50

      ! 'Subroutine' F is called (KFLAG=0 for correct return)
      kflag = 0
      GOTO 70

      ! STEP Q. Quotient acceptance (rare case)
      40 IF (fy-u*fy <= py*EXP(px-fx)) RETURN

      ! STEP E. Exponential sample - rand_exp() for standard exponential
      !         Deviate E and sample T from the Laplace 'Hat'
      !         (If T <= -.6744 then Pk < Fk for all MU >= 10.)
      50 e = MAXVAL(rand_exp(1))
      CALL RANDOM_NUMBER(u)
      u = u + u - 1.
      t = 1.8 + SIGN(e, u)
      IF (t <= (-.6744)) GOTO 50
      ival = INT(mu + s*t)
      fk = ival
      difmuk = mu - fk

      ! 'Subroutine' F is called (KFLAG=1 for correct return)
      kflag = 1
      GOTO 70

      ! STEP H. Hat acceptance (E is repeated on rejection)
      60 IF (c*ABS(u) > py*EXP(px+e) - fy*EXP(fx+e)) GOTO 50
      RETURN

      ! STEP F. 'Subroutine' F. Calculation of Px, Py, Fx, Fy.
      !         Case IVAL < 10 uses factorials from table fact
      70 IF (ival>=10) GOTO 80
      px = -mu
      py = mu**ival/fact(ival+1)
      GOTO 110

      ! Case IVAL >= 10 uses polynomial aproximation
      ! A0-A7 for accuracy when advisable
      ! .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
      80 del = .8333333E-1/fk
      del = del - 4.8*del*del*del
      v = difmuk/fk
      IF (ABS(v)>0.25) THEN
        px = fk*LOG(1. + v) - difmuk - del
      ELSE
        px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) - del
      END IF
      py = .3989423/SQRT(fk)
      110 x = (0.5 - difmuk)/s
      xx = x*x
      fx = -0.5*xx
      fy = omega * (((c3*xx + c2)*xx + c1)*xx + c0)
      IF (kflag <= 0) GOTO 40
      GOTO 60

    ELSE
    ! CASE B. mu < 10
    !-------
    ! Start new table and calculate P0 if necessary

      IF (first) THEN
        m = MAX(1,INT(mu))
        l = 0
        p = EXP(-mu)
        q = p
        p0 = p
      END IF

      ! STEP U. Uniform sample for inversion method
      DO
        CALL RANDOM_NUMBER(u)
        ival = 0
        IF (u <= p0) RETURN

        ! STEP T. Table comparison until the end PP(L) of the
        !         PP-Table of cumulative Poisson probabilities
        !         (0.458=PP(9) for MU=10)
        IF (l == 0) GOTO 150
        j = 1
        IF (u > 0.458) j = MIN(l, m)
        DO k=j,l
          IF (u <= pp(k)) GOTO 180
        END DO
        IF (l == 35) CYCLE

        ! STEP C. Creation of new Poisson probabilities P
        !         and their cumulatives Q=PP(K)
        150 l = l + 1
        DO k=l,35
          p = p*mu / k
          q = q + p
          pp(k) = q
          IF (u <= q) GOTO 170
        END DO
        l = 35
      END DO

      170 l = k
      180 ival = k
      RETURN

    END IF valofmu

    RETURN

    !-----------------------------------------------------------------------

  END FUNCTION gen_poisson


  !==========================================================================
  ! x = RAND_POISSON (N,tau,first)
  !
  !   Returns an array[N] of random poissonian variables with a timescale of 
  ! TAU. Taken from http://users.bigpond.net.au/amiller/ .
  !==========================================================================

  FUNCTION rand_poisson (N,tau,first)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N
    REAL(DP), INTENT(IN) :: tau
    INTEGER, DIMENSION(N) :: rand_poisson
    LOGICAL, INTENT(IN), OPTIONAL :: first

    INTEGER :: i
    LOGICAL :: prem

    !-----------------------------------------------------------------------

    IF (PRESENT(first)) THEN ; prem = first ; ELSE ; prem = .True. ; END IF

    rand_poisson(1) = GEN_POISSON(1._DP/tau,FIRST=prem)
    rand_poisson(2:N) = (/(GEN_POISSON(1._DP/tau,FIRST=.False.),i=2,N)/)

    !-----------------------------------------------------------------------

  END FUNCTION rand_poisson


  !==========================================================================
  !  x[N] = RAND_GENERAL(N,func,limits,xlog,ylog,limited,accmax)
  !
  !  Draw random sample from a user defined distribution.
  !==========================================================================

  FUNCTION rand_general_vec (N,func,limits,xlog,ylog,lnfunc,verbose,limited, &
                             accuracy,Nmax,Nsamp,xgrid,pgrid,Fgrid,proper)

    USE utilities, ONLY: DP    
    USE arrays, ONLY: ramp, uniq_sorted
    USE interpolation, ONLY: interp_lin_sorted
    USE adaptative_grid, ONLY: gridadapt1D
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N
    INTERFACE
      FUNCTION func (x)
        USE utilities, ONLY: DP
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP), DIMENSION(SIZE(x)) :: func
      END FUNCTION func
    END INTERFACE
    REAL(DP), DIMENSION(2), INTENT(IN) :: limits
    LOGICAL, INTENT(IN), OPTIONAL :: xlog, ylog, lnfunc, verbose
    LOGICAL, DIMENSION(2), INTENT(IN), OPTIONAL :: limited
    REAL(DP), INTENT(IN), OPTIONAL :: accuracy
    INTEGER, INTENT(IN), OPTIONAL :: Nmax
    INTEGER, INTENT(OUT), OPTIONAL :: Nsamp
    REAL(DP), INTENT(OUT), DIMENSION(:), ALLOCATABLE, OPTIONAL :: xgrid, pgrid
    REAL(DP), INTENT(OUT), DIMENSION(:), ALLOCATABLE, OPTIONAL :: Fgrid
    LOGICAL, INTENT(OUT), DIMENSION(N), OPTIONAL :: proper
    REAL(DP), DIMENSION(N) :: rand_general_vec

    INTEGER, PARAMETER :: Ncoarse = 10
    REAL(DP), PARAMETER :: accmax0 = 1.E-3_DP

    INTEGER :: Nfine, Ncdf
    REAL(DP) :: acc, xinf0, xsup1, accmax
    REAL(DP), DIMENSION(Ncoarse) :: xcoarse
    REAL(DP), DIMENSION(N) :: rand
    REAL(DP), DIMENSION(2) :: lim
    REAL(DP), DIMENSION(:), ALLOCATABLE :: xfine, yfine, prim, xuniq, Funiq
    LOGICAL :: xl, yl, lnfn, verb
    LOGICAL, DIMENSION(2) :: adjxlim
    LOGICAL, DIMENSION(N) :: prop
    LOGICAL, DIMENSION(:), ALLOCATABLE :: keepinf, keepsup, keep

    !-----------------------------------------------------------------------

    ! Log axes
    xl = .False.
    IF (PRESENT(xlog)) xl = xlog
    yl = .False.
    IF (PRESENT(ylog)) yl = ylog

    ! User defined function
    lnfn = .False.
    IF (PRESENT(lnfunc)) lnfn = lnfunc

    ! Limited or not
    adjxlim(:) = .False.
    IF (PRESENT(limited)) adjxlim(:) = ( .NOT. limited(:) )

    ! Random drawing
    IF (PRESENT(accuracy)) THEN ; accmax = accuracy 
                           ELSE ; accmax = accmax0 ; END IF
    CALL RANDOM_NUMBER(rand(:))
    acc = MIN( accmax*ABS(MINVAL(rand(:))), accmax*ABS(1._DP-MAXVAL(rand(:))) )
    
    ! Coarse grid
    xcoarse(:) = RAMP(Ncoarse,limits(1),limits(2),XLOG=xl)

    ! Adpatative grid
    IF (PRESENT(verbose)) THEN ; verb = verbose ; ELSE ; verb = .False. ; END IF
    CALL GRIDADAPT1D(xcoarse,xfine,yfine,func,acc,NMAX=Nmax,INTEG=.True., &
                     PRIMITIVE=prim,XLOG=xl,YLOG=yl,VERBOSE=verb,LNFUNC=lnfn, &
                     RESCALE=.True.,SLIM=.False.,ADJUSTXLIM=adjxlim(:))
    Nfine = SIZE(xfine(:))
    IF (PRESENT(Nsamp)) Nsamp = Nfine
    prim(:) = prim(:) / prim(Nfine)
    
    ! Ensure the most compact support: if there are several CDF=0, take the 
    ! largest x, if there are several CDF=1 take the smaller x.
    ALLOCATE (keepinf(Nfine),keepsup(Nfine),keep(Nfine))
    keepinf(:) = UNIQ_SORTED(prim(:),LAST=.True.)
    keepsup(:) = UNIQ_SORTED(prim(:),FIRST=.True.)
    keep(:) = MERGE(keepinf(:),keepsup(:),prim(:) < 0.5_DP)
    
    ! Remove values outside of support
    xinf0 = MAXVAL(PACK(xfine(:),prim(:) == 0._DP))
    xsup1 = MINVAL(PACK(xfine(:),prim(:) == 1._DP))
    keep(:) = MERGE(keep(:),.False.,xfine(:) >= xinf0)
    keep(:) = MERGE(keep(:),.False.,xfine(:) <= xsup1)
    
    ! Pack the CDF
    Ncdf = COUNT(keep)
    ALLOCATE (Funiq(Ncdf),xuniq(Ncdf))
    Funiq(:) = PACK(prim(:),keep(:))
    xuniq(:) = PACK(xfine(:),keep(:))
    
    ! Interpolate the values (y->x)
    rand_general_vec(:) = INTERP_LIN_SORTED(xuniq(:),Funiq(:),rand(:),YLOG=xl)
    
    ! Make sure the value does not get out of the range
    lim(1) = MAX( xuniq(1), INTERP_LIN_SORTED(xuniq(:),Funiq(:),0._DP) )
    lim(2) = MIN( xuniq(Ncdf), INTERP_LIN_SORTED(xuniq(:),Funiq(:),1._DP) )
    prop(:) = .True.
    WHERE (rand_general_vec(:) < lim(1))
      rand_general_vec(:) = lim(1)
      prop(:) = .False.
    END WHERE
    WHERE (rand_general_vec(:) > lim(2))
      rand_general_vec(:) = lim(2)
      prop(:) = .False.
    END WHERE
    IF (PRESENT(proper)) proper(:) = prop(:)
    
    ! For debugging
    IF (PRESENT(xgrid)) THEN
      ALLOCATE (xgrid(Ncdf))
      xgrid(:) = xuniq(:)
    END IF
    IF (PRESENT(pgrid)) THEN
      ALLOCATE (pgrid(Ncdf))
      pgrid(:) = PACK(yfine(:),keep(:))
    END IF
    IF (PRESENT(Fgrid)) THEN
      ALLOCATE (Fgrid(Ncdf))
      Fgrid(:) = Funiq(:)
    END IF

    ! Clean memory space
    DEALLOCATE (xfine,yfine,prim,xuniq,Funiq,keepinf,keepsup,keep)
    
    !-----------------------------------------------------------------------

  END FUNCTION rand_general_vec

  ! Interface shell
  FUNCTION rand_general_scl (func,limits,xlog,ylog,lnfunc,verbose,limited, &
                             accuracy,Nmax,Nsamp,xgrid,pgrid,Fgrid,proper)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    INTERFACE
      FUNCTION func (x)
        USE utilities, ONLY: DP
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP), DIMENSION(SIZE(x)) :: func
      END FUNCTION func
    END INTERFACE
    REAL(DP), DIMENSION(2), INTENT(IN) :: limits
    LOGICAL, INTENT(IN), OPTIONAL :: xlog, ylog, lnfunc, verbose
    LOGICAL, DIMENSION(2), INTENT(IN), OPTIONAL :: limited
    REAL(DP), INTENT(IN), OPTIONAL :: accuracy
    INTEGER, INTENT(IN), OPTIONAL :: Nmax
    INTEGER, INTENT(OUT), OPTIONAL :: Nsamp
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: xgrid, pgrid
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Fgrid
    LOGICAL, INTENT(OUT), OPTIONAL :: proper
    REAL(DP) :: rand_general_scl
    LOGICAL, DIMENSION(1) :: prop
    rand_general_scl = MINVAL(RAND_GENERAL_VEC(1,func,limits,XLOG=xlog, &
                                               YLOG=ylog,VERBOSE=verbose, &
                                               LNFUNC=lnfunc,LIMITED=limited, &
                                               ACCURACY=accuracy,NMAX=Nmax, &
                                               NSAMP=Nsamp,XGRID=xgrid, & 
                                               PGRID=pgrid,FGRID=Fgrid, &
                                               PROPER=prop))
    IF (PRESENT(proper)) proper = prop(1)
  END FUNCTION rand_general_scl

  !-------------------------------------------------------------------------

END MODULE random
