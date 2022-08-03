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
  !    - 05/2015: update with Student's t and multivariate forms.
  !    - 04/2016: add the skew-normal distribution.
  !    - 05/2016: add the split-normal distribution.
  ! 
  ! 3) DESCRIPTION: implements various commonly used mathematical 
  !                 distributions, and operations to manipulate them.
  !==========================================================================


! TODO: dist_splitnorm and param_splitnorm in 2D

MODULE distributions

  USE utilities, ONLY:
  USE constants, ONLY:
  USE arrays, ONLY:
  USE special_functions, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: dist_power, dist_gauss, dist_lorentz, dist_lognormal, &
            dist_skewnorm, param_skewnorm, dist_splitnorm, param_splitnorm, &
            dist_splitlorentz, &
            histogram1D, histogram1Dmulti, histogram2D

  INTERFACE dist_gauss
    MODULE PROCEDURE dist_gauss_1D_scl, dist_gauss_1D_vect
    MODULE PROCEDURE dist_gauss_ND_scl, dist_gauss_ND_vect
  END INTERFACE dist_gauss
  
  INTERFACE dist_skewnorm
    MODULE PROCEDURE dist_skewnorm_1D_scl, dist_skewnorm_1D_vect
  END INTERFACE dist_skewnorm

  INTERFACE param_skewnorm
    MODULE PROCEDURE param_skewnorm_0D, param_skewnorm_1D
    MODULE PROCEDURE param_skewnorm_2D, param_skewnorm_3D
  END INTERFACE param_skewnorm
  
  INTERFACE dist_splitnorm
    MODULE PROCEDURE dist_splitnorm_1D_scl, dist_splitnorm_1D_vect
  END INTERFACE dist_splitnorm

  INTERFACE param_splitnorm
    MODULE PROCEDURE param_splitnorm_0D, param_splitnorm_1D
    MODULE PROCEDURE param_splitnorm_2D, param_splitnorm_3D
  END INTERFACE param_splitnorm
  
  INTERFACE dist_lognormal
    MODULE PROCEDURE dist_lognormal_scl, dist_lognormal_1D
  END INTERFACE dist_lognormal

  INTERFACE dist_lorentz
    MODULE PROCEDURE dist_lorentz_scl, dist_lorentz_1D
  END INTERFACE dist_lorentz

  INTERFACE dist_splitlorentz
    MODULE PROCEDURE dist_splitlorentz_scl, dist_splitlorentz_1D
  END INTERFACE dist_splitlorentz

CONTAINS


  !==========================================================================
  ! f[N] = DIST_POWER(x[N],alpha,xinf,xsup)
  !
  !   Returns a normalised power-law distribution of index alpha, between
  ! xinf and xsup, for the input array x. It is singular in -1.
  !==========================================================================

  FUNCTION dist_power (x,alpha,xinf,xsup)

    USE utilities, ONLY: DP, strike
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), INTENT(IN) :: alpha, xinf, xsup
    REAL(DP), DIMENSION(SIZE(x)) :: dist_power

    REAL(DP) :: norm, ap1

    !-----------------------------------------------------------------------

    ! Check arguments
    IF (alpha == -1._DP) &
      CALL STRIKE ("DIST_POWER","singular distribution (alpha=-1).")

    ! Definition
    ap1 = alpha+1._DP
    norm = ap1 / ( xsup**ap1 - xinf**ap1 )
    WHERE (x >= xinf .AND. x <= xsup)
      dist_power = norm * x**alpha
    ELSEWHERE
      dist_power = 0._DP
    ENDWHERE

    !-----------------------------------------------------------------------

  END FUNCTION dist_power


  !==========================================================================
  ! f[N] = DIST_GAUSS(x[N],xmean,sigma)
  ! f[N] = DIST_GAUSS(x[N,M],xmean[M],covariance[M,M])
  !
  !   Returns a normalised Gauss distribution centered n xmean, with a 
  ! standard deviation sigma.
  !==========================================================================

  ! 1 variable, scalar
  PURE FUNCTION dist_gauss_1D_scl (x,xmean,sigma)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi, oneoversqrt2pi
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: xmean, sigma
    REAL(DP) :: dist_gauss_1D_scl

    REAL(DP) :: norm, twice_sig2

    !-----------------------------------------------------------------------

    norm = oneoversqrt2pi / sigma
    twice_sig2 = 2._DP*sigma**2
    dist_gauss_1D_scl = norm * EXP(-(x-xmean)**2/twice_sig2)

    !-----------------------------------------------------------------------

  END FUNCTION dist_gauss_1D_scl

  ! 1 variable, vector
  PURE FUNCTION dist_gauss_1D_vect (x,xmean,sigma)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi, oneoversqrt2pi
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), INTENT(IN) :: xmean, sigma
    REAL(DP), DIMENSION(SIZE(x)) :: dist_gauss_1D_vect

    REAL(DP) :: norm, twice_sig2

    !-----------------------------------------------------------------------

    norm = oneoversqrt2pi / sigma
    twice_sig2 = 2._DP*sigma**2
    dist_gauss_1D_vect(:) = norm * EXP(-(x(:)-xmean)**2/twice_sig2)

    !-----------------------------------------------------------------------

  END FUNCTION dist_gauss_1D_vect

  ! Multivariate, scalar 
  PURE FUNCTION dist_gauss_ND_scl (x,xmean,invcov,detcov)

    USE utilities, ONLY: DP
    USE constants, ONLY: twopi
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), DIMENSION(SIZE(x)), INTENT(IN) :: xmean
    REAL(DP), DIMENSION(SIZE(x),SIZE(x)), INTENT(IN) :: invcov
    REAL(DP), INTENT(IN) :: detcov
    REAL(DP) :: dist_gauss_ND_scl

    INTEGER :: N
    REAL(DP) :: norm, quadform
    REAL(DP), DIMENSION(SIZE(x)) :: x0

    !-----------------------------------------------------------------------

    N = SIZE(x(:))
    norm = 1._DP / SQRT( twopi**N * detcov )
    x0(:) = x(:) - xmean(:)
    quadform = DOT_PRODUCT( x0(:), MATMUL(invcov(:,:),x0(:)) )
    dist_gauss_ND_scl = norm * EXP(- 0.5_DP * quadform )    

    !-----------------------------------------------------------------------

  END FUNCTION dist_gauss_ND_scl


  ! Multivariate, scalar
  PURE FUNCTION dist_gauss_ND_vect (x,xmean,invcov,detcov)

    USE utilities, ONLY: DP
    USE constants, ONLY: twopi
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:,:) :: x
    REAL(DP), DIMENSION(SIZE(x,2)), INTENT(IN) :: xmean
    REAL(DP), DIMENSION(SIZE(x,2),SIZE(x,2)), INTENT(IN) :: invcov
    REAL(DP), INTENT(IN) :: detcov
    REAL(DP), DIMENSION(SIZE(x,1)) :: dist_gauss_ND_vect

    INTEGER :: i, N
    REAL(DP) :: norm
    REAL(DP), DIMENSION(SIZE(x,1)) :: quadform
    REAL(DP), DIMENSION(SIZE(x,1),SIZE(x,2)) :: x0

    !-----------------------------------------------------------------------

    N = SIZE(x(:,:),1)
    norm = 1._DP / SQRT( twopi**N * detcov )
    FORALL (i=1:N) x0(i,:) = x(i,:) - xmean(:)
    FORALL (i=1:N) &
      quadform(i) = DOT_PRODUCT( x0(i,:), MATMUL(invcov(:,:),x0(i,:)) )
    dist_gauss_ND_vect(:) = norm * EXP(- 0.5_DP * quadform(:) )

    !-----------------------------------------------------------------------

  END FUNCTION dist_gauss_ND_vect


  !==========================================================================
  ! f[N] = DIST_SKEWNORM(x[N],ksi,omega,alpha)
  !
  !   Returns a normalised skew-normal distribution with position parameter
  !  ksi, scale parameter omega and shape parameter alpha.
  !==========================================================================

  ! 1D vector
  PURE FUNCTION dist_skewnorm_1D_vect (x,ksi,omega,alpha)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), INTENT(IN) :: ksi, omega, alpha
    REAL(DP), DIMENSION(SIZE(x)) :: dist_skewnorm_1D_vect

    REAL(DP), DIMENSION(SIZE(x)) :: xred, gauss, repart

    !-----------------------------------------------------------------------

    xred(:) = ( x(:) - ksi ) / omega
    gauss(:) = EXP( -0.5_DP * xred(:)**2 ) / SQRT(2*pi)
    repart(:) = 0.5_DP * ( 1 + ERF(alpha*xred(:)/SQRT(2._DP)) )
    dist_skewnorm_1D_vect(:) = 2._DP/omega * gauss(:) * repart(:)

    !-----------------------------------------------------------------------

  END FUNCTION dist_skewnorm_1D_vect

  ! 1D scalar
  PURE FUNCTION dist_skewnorm_1D_scl (x,ksi,omega,alpha)
    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: ksi, omega, alpha
    REAL(DP) :: dist_skewnorm_1D_scl
    dist_skewnorm_1D_scl = MAXVAL(DIST_SKEWNORM_1D_VECT([x],ksi,omega,alpha))
  END FUNCTION dist_skewnorm_1D_scl


  !==========================================================================
  ! CALL PARAM_SKEWNORM(ksi,omega,alpha,mean,sigma,skewness,param2moment, &
  !                     moment2param)
  !
  !   Conversion between parameters of the skew-normal distribution and its
  ! moments.
  !==========================================================================

  PURE SUBROUTINE param_skewnorm_3D (ksi,omega,alpha,mean,sigma,skewness, &
                                     param2moment,moment2param)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: ksi, omega, alpha
    REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: mean, sigma, skewness
    LOGICAL, INTENT(IN), OPTIONAL :: param2moment, moment2param

    REAL(DP), PARAMETER :: eps = 1.E-4_DP
    REAL(DP), DIMENSION(SIZE(ksi,1),SIZE(ksi,2),SIZE(ksi,3)) :: delta
    LOGICAL :: param
    
    !-----------------------------------------------------------------------

    ! Determine which way to go: param -> moments or moments -> param
    IF (PRESENT(param2moment)) THEN
      param = ( .NOT. param2moment )
    ELSE IF (PRESENT(moment2param)) THEN
      param = moment2param
    ELSE
      param = .False.
    END IF

    ! Compute the parameters
    IF (param) THEN
      delta(:,:,:) = SIGN(1._DP,skewness(:,:,:)) &
                   * ABS(skewness(:,:,:))**(1._DP/3._DP) &
                   / SQRT( ((4._DP-pi)/2._DP)**(2._DP/3._DP)*2/pi &
                         + 2*ABS(skewness(:,:,:))**(2._DP/3._DP)/pi )
      WHERE (ABS(delta(:,:,:)) > 1._DP-eps) &
        delta(:,:,:) = SIGN(1._DP,delta(:,:,:)) * (1-eps)
      ksi(:,:,:) = mean(:,:,:) &
                 - delta(:,:,:) * sigma(:,:,:) &
                   * SQRT(2._DP / (pi-2*delta(:,:,:)**2))
      omega(:,:,:) = sigma(:,:,:) / SQRT( 1 - 2*delta(:,:,:)**2/pi )
      alpha(:,:,:) = delta(:,:,:) / SQRT( 1 - delta(:,:,:)**2 )
    ELSE
      delta(:,:,:) = alpha(:,:,:) / SQRT( 1 + alpha(:,:,:)**2 )
      mean(:,:,:) = ksi(:,:,:) + omega(:,:,:) * delta(:,:,:) * SQRT(2._DP/pi)
      sigma(:,:,:) = omega(:,:,:) * SQRT( 1._DP - 2*delta(:,:,:)**2/pi )
      skewness(:,:,:) = (4._DP-pi)/2._DP * (delta(:,:,:)*SQRT(2._DP/pi))**3 &
                      / (1._DP-2*delta(:,:,:)**2/pi)**1.5
    END IF
    
    !-----------------------------------------------------------------------
    
  END SUBROUTINE param_skewnorm_3D

  ! 2D
  PURE SUBROUTINE param_skewnorm_2D (ksi,omega,alpha,mean,sigma,skewness, &
                                     param2moment,moment2param)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: ksi, omega, alpha
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mean, sigma, skewness
    LOGICAL, INTENT(IN), OPTIONAL :: param2moment, moment2param

    REAL(DP), PARAMETER :: eps = 1.E-4_DP
    REAL(DP), DIMENSION(SIZE(ksi,1),SIZE(ksi,2)) :: delta
    LOGICAL :: param
    
    !-----------------------------------------------------------------------

    ! Determine which way to go: param -> moments or moments -> param
    IF (PRESENT(param2moment)) THEN
      param = ( .NOT. param2moment )
    ELSE IF (PRESENT(moment2param)) THEN
      param = moment2param
    ELSE
      param = .False.
    END IF

    ! Compute the parameters
    IF (param) THEN
      delta(:,:) = SIGN(1._DP,skewness(:,:)) &
                 * ABS(skewness(:,:))**(1._DP/3._DP) &
                 / SQRT( ((4._DP-pi)/2._DP)**(2._DP/3._DP)*2/pi &
                       + 2*ABS(skewness(:,:))**(2._DP/3._DP)/pi )
      WHERE (ABS(delta(:,:)) > 1._DP-eps) &
        delta(:,:) = SIGN(1._DP,delta(:,:)) * (1-eps)
      ksi(:,:) = mean(:,:) &
               - delta(:,:) * sigma(:,:) * SQRT(2._DP / (pi-2*delta(:,:)**2))
      omega(:,:) = sigma(:,:) / SQRT( 1 - 2*delta(:,:)**2/pi )
      alpha(:,:) = delta(:,:) / SQRT( 1 - delta(:,:)**2 )
    ELSE
      delta(:,:) = alpha(:,:) / SQRT( 1 + alpha(:,:)**2 )
      mean(:,:) = ksi(:,:) + omega(:,:) * delta(:,:) * SQRT(2._DP/pi)
      sigma(:,:) = omega(:,:) * SQRT( 1._DP - 2*delta(:,:)**2/pi )
      skewness(:,:) = (4._DP-pi)/2._DP * (delta(:,:)*SQRT(2._DP/pi))**3 &
                  / (1._DP-2*delta(:,:)**2/pi)**1.5
    END IF
    
    !-----------------------------------------------------------------------
    
  END SUBROUTINE param_skewnorm_2D

  ! 1D
  PURE SUBROUTINE param_skewnorm_1D (ksi,omega,alpha,mean,sigma,skewness, &
                                     param2moment,moment2param)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(INOUT) :: ksi, omega, alpha
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: mean, sigma, skewness
    LOGICAL, INTENT(IN), OPTIONAL :: param2moment, moment2param

    REAL(DP), PARAMETER :: eps = 1.E-4_DP
    REAL(DP), DIMENSION(SIZE(ksi)) :: delta
    LOGICAL :: param
    
    !-----------------------------------------------------------------------

    ! Determine which way to go: param -> moments or moments -> param
    IF (PRESENT(param2moment)) THEN
      param = ( .NOT. param2moment )
    ELSE IF (PRESENT(moment2param)) THEN
      param = moment2param
    ELSE
      param = .False.
    END IF

    ! Compute the parameters
    IF (param) THEN
      delta(:) = SIGN(1._DP,skewness(:)) * ABS(skewness(:))**(1._DP/3._DP) &
               / SQRT( ((4._DP-pi)/2._DP)**(2._DP/3._DP)*2/pi &
                     + 2*ABS(skewness(:))**(2._DP/3._DP)/pi )
      WHERE (ABS(delta(:)) > 1._DP-eps) &
        delta(:) = SIGN(1._DP,delta(:)) * (1-eps)
      ksi(:) = mean(:) - delta(:) * sigma(:) * SQRT(2._DP / (pi-2*delta(:)**2))
      omega(:) = sigma(:) / SQRT( 1 - 2*delta(:)**2/pi )
      alpha(:) = delta(:) / SQRT( 1 - delta(:)**2 )
    ELSE
      delta(:) = alpha(:) / SQRT( 1 + alpha(:)**2 )
      mean(:) = ksi(:) + omega(:) * delta(:) * SQRT(2._DP/pi)
      sigma(:) = omega(:) * SQRT( 1._DP - 2*delta(:)**2/pi )
      skewness(:) = (4._DP-pi)/2._DP * (delta(:)*SQRT(2._DP/pi))**3 &
                  / (1._DP-2*delta(:)**2/pi)**1.5
    END IF
    
    !-----------------------------------------------------------------------
    
  END SUBROUTINE param_skewnorm_1D

  ! Scalar interface
  ELEMENTAL SUBROUTINE param_skewnorm_0D (ksi,omega,alpha, &
                                          mean,sigma,skewness, &
                                          param2moment,moment2param)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    IMPLICIT NONE

    REAL(DP), INTENT(INOUT) :: ksi, omega, alpha, mean, sigma, skewness
    LOGICAL, INTENT(IN), OPTIONAL :: param2moment, moment2param
    
    REAL(DP), PARAMETER :: eps = 1.E-4_DP
    REAL(DP) :: delta
    LOGICAL :: param
    
    !-----------------------------------------------------------------------

    ! Determine which way to go: param -> moments or moments -> param
    IF (PRESENT(param2moment)) THEN
      param = ( .NOT. param2moment )
    ELSE IF (PRESENT(moment2param)) THEN
      param = moment2param
    ELSE
      param = .False.
    END IF

    ! Compute the parameters
    IF (param) THEN
      delta = SIGN(1._DP,skewness) * ABS(skewness)**(1._DP/3._DP) &
            / SQRT( ((4._DP-pi)/2._DP)**(2._DP/3._DP)*2/pi &
                  + 2*ABS(skewness)**(2._DP/3._DP)/pi )
      IF (ABS(delta) > 1._DP-eps) delta = SIGN(1._DP,delta) * (1-eps)
      ksi = mean - delta * sigma * SQRT(2._DP / (pi-2*delta**2))
      omega = sigma / SQRT( 1 - 2*delta**2/pi )
      alpha = delta / SQRT( 1 - delta**2 )
    ELSE
      delta = alpha / SQRT( 1 + alpha**2 )
      mean = ksi + omega * delta * SQRT(2._DP/pi)
      sigma = omega * SQRT( 1._DP - 2*delta**2/pi )
      skewness = (4._DP-pi)/2._DP * (delta*SQRT(2._DP/pi))**3 &
               / (1._DP-2*delta**2/pi)**1.5
    END IF
    
    !-----------------------------------------------------------------------
    
  END SUBROUTINE param_skewnorm_0D

  
  !==========================================================================
  ! f[N] = DIST_SPLITNORM(x[N],ksi,omega,alpha)
  !
  !   Returns a normalised split-normal distribution with position parameter
  !  mu, scale parameter lambda and shape parameter tau.
  !==========================================================================

  ! 1D vector
  PURE FUNCTION dist_splitnorm_1D_vect (x,mu,lambda,tau)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), INTENT(IN) :: mu, lambda, tau
    REAL(DP), DIMENSION(SIZE(x)) :: dist_splitnorm_1D_vect

    REAL(DP) :: A, tau2
    REAL(DP), DIMENSION(SIZE(x)) :: xred2
    
    !-----------------------------------------------------------------------

    A = SQRT(2/pi) / lambda / (1+tau)
    tau2 = tau**2
    xred2(:) = - 0.5_DP * ( x(:) - mu )**2 / lambda**2
    WHERE (x(:) < mu)
      dist_splitnorm_1D_vect(:) = A * EXP(xred2(:))
    ELSEWHERE
      dist_splitnorm_1D_vect(:) = A * EXP(xred2(:)/tau2)
    END WHERE
    
    !-----------------------------------------------------------------------

  END FUNCTION dist_splitnorm_1D_vect

  ! 1D scalar
  PURE FUNCTION dist_splitnorm_1D_scl (x,mu,lambda,tau)
    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: mu, lambda, tau
    REAL(DP) :: dist_splitnorm_1D_scl
    dist_splitnorm_1D_scl = MAXVAL(DIST_SPLITNORM_1D_VECT([x],mu,lambda,tau))
  END FUNCTION dist_splitnorm_1D_scl


  !==========================================================================
  ! CALL PARAM_SPLITNORM(mu,lambda,tau,mean,sigma,skewness,param2moment, &
  !                     moment2param)
  !
  !   Conversion between parameters of the split-normal distribution and its
  ! moments.
  !==========================================================================

  PURE SUBROUTINE param_splitnorm_3D (mu,lambda,tau,mean,sigma,skewness, &
                                      param2moment,moment2param)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    USE arrays, ONLY: ramp
    USE interpolation, ONLY: locate_sorted
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: mu, lambda, tau
    REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: mean, sigma, skewness
    LOGICAL, INTENT(IN), OPTIONAL :: param2moment, moment2param

    INTEGER, PARAMETER :: N = 10000
    INTEGER :: x, y, z, Nx, Ny, Nz
    INTEGER, DIMENSION(SIZE(mu,1),SIZE(mu,2),SIZE(mu,3)) :: i1, i2
    REAL(DP) :: sqrt2ovPi
    REAL(DP), DIMENSION(N) :: taugrid, bgrid, skewgrid, lntaugrid
    REAL(DP), DIMENSION(SIZE(mu,1),SIZE(mu,2),SIZE(mu,3)) :: b
    LOGICAL :: param
    
    !-----------------------------------------------------------------------

    ! Determine which way to go: param -> moments or moments -> param
    IF (PRESENT(param2moment)) THEN
      param = ( .NOT. param2moment )
    ELSE IF (PRESENT(moment2param)) THEN
      param = moment2param
    ELSE
      param = .False.
    END IF

    ! Compute the parameters (compute tau numerically)
    sqrt2ovPi = SQRT(2._DP/pi)
    IF (param) THEN
      Nx = SIZE(mu(:,:,:),1)
      Ny = SIZE(mu(:,:,:),2)
      Nz = SIZE(mu(:,:,:),3)
      taugrid(:) = RAMP(N,1.E-2_DP,1.E2_DP,XLOG=.True.) ! tau=1.E2 => skew=0.99
      lntaugrid(:) = LOG(taugrid(:))
      bgrid(:) = (pi-2._DP)/pi * (taugrid(:)-1._DP)**2 + taugrid(:)
      skewgrid(:) = bgrid(:)**(-1.5_DP) * sqrt2ovPi * (taugrid(:)-1._DP) &
                  * ( (4/pi-1._DP)*(taugrid(:)-1._DP)**2 + taugrid(:) )
      FORALL (x=1:Nx,y=1:Ny,z=1:Nz) &
        i1(x,y,z) = LOCATE_SORTED(skewgrid(:),skewness(x,y,z))
      i2(:,:,:) = i1(:,:,:) + 1
      FORALL (x=1:Nx,y=1:Ny,z=1:Nz) &
        tau(x,y,z) = EXP( lntaugrid(i1(x,y,z)) &
                        + (lntaugrid(i2(x,y,z))-lntaugrid(i1(x,y,z))) &
                        / (skewgrid(i2(x,y,z))-skewgrid(i1(x,y,z))) &
                        * (skewness(x,y,z)-skewgrid(i1(x,y,z))) )
      b(:,:,:) = (pi-2._DP)/pi * (tau(:,:,:)-1._DP)**2 + tau(:,:,:)
      lambda(:,:,:) = sigma(:,:,:) / SQRT(b(:,:,:))
      mu(:,:,:) = mean(:,:,:) - sqrt2ovPi * lambda(:,:,:) * (tau(:,:,:)-1._DP)
    ELSE
      b(:,:,:) = (pi-2._DP)/pi * (tau(:,:,:)-1._DP)**2 + tau(:,:,:)
      mean(:,:,:) = mu(:,:,:) + sqrt2ovPi * lambda(:,:,:) * (tau(:,:,:)-1._DP)
      sigma(:,:,:) = SQRT(b(:,:,:)) * lambda(:,:,:)
      skewness(:,:,:) = b(:,:,:)**(-1.5_DP) * sqrt2ovPi * (tau(:,:,:)-1._DP) &
                      * ( (4/pi-1._DP)*(tau(:,:,:)-1._DP)**2 + tau(:,:,:) )
    END IF
    
    !-----------------------------------------------------------------------
    
  END SUBROUTINE param_splitnorm_3D

  ! 2D
  PURE SUBROUTINE param_splitnorm_2D (mu,lambda,tau,mean,sigma,skewness, &
                                      param2moment,moment2param)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    USE arrays, ONLY: ramp
    USE interpolation, ONLY: locate_sorted
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mu, lambda, tau
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mean, sigma, skewness
    LOGICAL, INTENT(IN), OPTIONAL :: param2moment, moment2param

    INTEGER, PARAMETER :: N = 10000
    INTEGER :: x, y, Nx, Ny
    INTEGER, DIMENSION(SIZE(mu,1),SIZE(mu,2)) :: i1, i2
    REAL(DP) :: sqrt2ovPi
    REAL(DP), DIMENSION(N) :: taugrid, bgrid, skewgrid, lntaugrid
    REAL(DP), DIMENSION(SIZE(mu,1),SIZE(mu,2)) :: b
    LOGICAL :: param
    
    !-----------------------------------------------------------------------

    ! Determine which way to go: param -> moments or moments -> param
    IF (PRESENT(param2moment)) THEN
      param = ( .NOT. param2moment )
    ELSE IF (PRESENT(moment2param)) THEN
      param = moment2param
    ELSE
      param = .False.
    END IF

    ! Compute the parameters (compute tau numerically)
    sqrt2ovPi = SQRT(2._DP/pi)
    IF (param) THEN
      Nx = SIZE(mu(:,:),1)
      Ny = SIZE(mu(:,:),2)
      taugrid(:) = RAMP(N,1.E-2_DP,1.E2_DP,XLOG=.True.) ! tau=1.E2 => skew=0.99
      lntaugrid(:) = LOG(taugrid(:))
      bgrid(:) = (pi-2._DP)/pi * (taugrid(:)-1._DP)**2 + taugrid(:)
      skewgrid(:) = bgrid(:)**(-1.5_DP) * sqrt2ovPi * (taugrid(:)-1._DP) &
                  * ( (4/pi-1._DP)*(taugrid(:)-1._DP)**2 + taugrid(:) )
      FORALL (x=1:Nx,y=1:Ny) &
        i1(x,y) = LOCATE_SORTED(skewgrid(:),skewness(x,y))
      i2(:,:) = i1(:,:) + 1
      FORALL (x=1:Nx,y=1:Ny) &
        tau(x,y) = EXP( lntaugrid(i1(x,y)) &
                      + (lntaugrid(i2(x,y))-lntaugrid(i1(x,y))) &
                      / (skewgrid(i2(x,y))-skewgrid(i1(x,y))) &
                      * (skewness(x,y)-skewgrid(i1(x,y))) )
      b(:,:) = (pi-2._DP)/pi * (tau(:,:)-1._DP)**2 + tau(:,:)
      lambda(:,:) = sigma(:,:) / SQRT(b(:,:))
      mu(:,:) = mean(:,:) - sqrt2ovPi * lambda(:,:) * (tau(:,:)-1._DP)
    ELSE
      b(:,:) = (pi-2._DP)/pi * (tau(:,:)-1._DP)**2 + tau(:,:)
      mean(:,:) = mu(:,:) + sqrt2ovPi * lambda(:,:) * (tau(:,:)-1._DP)
      sigma(:,:) = SQRT(b(:,:)) * lambda(:,:)
      skewness(:,:) = b(:,:)**(-1.5_DP) * sqrt2ovPi * (tau(:,:)-1._DP) &
                    * ( (4/pi-1._DP)*(tau(:,:)-1._DP)**2 + tau(:,:) )
    END IF
    
    !-----------------------------------------------------------------------
    
  END SUBROUTINE param_splitnorm_2D

  ! 1D
  PURE SUBROUTINE param_splitnorm_1D (mu,lambda,tau,mean,sigma,skewness, &
                                      param2moment,moment2param)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    USE arrays, ONLY: ramp
    USE interpolation, ONLY: locate_sorted
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(INOUT) :: mu, lambda, tau
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: mean, sigma, skewness
    LOGICAL, INTENT(IN), OPTIONAL :: param2moment, moment2param

    INTEGER, PARAMETER :: N = 10000
    INTEGER :: x, Nx
    INTEGER, DIMENSION(SIZE(mu)) :: i1, i2
    REAL(DP) :: sqrt2ovPi
    REAL(DP), DIMENSION(N) :: taugrid, bgrid, skewgrid, lntaugrid
    REAL(DP), DIMENSION(SIZE(mu)) :: b
    LOGICAL :: param
    
    !-----------------------------------------------------------------------

    ! Determine which way to go: param -> moments or moments -> param
    IF (PRESENT(param2moment)) THEN
      param = ( .NOT. param2moment )
    ELSE IF (PRESENT(moment2param)) THEN
      param = moment2param
    ELSE
      param = .False.
    END IF

    ! Compute the parameters (compute tau numerically)
    sqrt2ovPi = SQRT(2._DP/pi)
    IF (param) THEN
      Nx = SIZE(mu(:))
      taugrid(:) = RAMP(N,1.E-2_DP,1.E2_DP,XLOG=.True.) ! tau=1.E2 => skew=0.99
      lntaugrid(:) = LOG(taugrid(:))
      bgrid(:) = (pi-2._DP)/pi * (taugrid(:)-1._DP)**2 + taugrid(:)
      skewgrid(:) = bgrid(:)**(-1.5_DP) * sqrt2ovPi * (taugrid(:)-1._DP) &
                  * ( (4/pi-1._DP)*(taugrid(:)-1._DP)**2 + taugrid(:) )
      FORALL (x=1:Nx) &
        i1(x) = LOCATE_SORTED(skewgrid(:),skewness(x))
      i2(:) = i1(:) + 1
      FORALL (x=1:Nx) &
        tau(x) = EXP( lntaugrid(i1(x)) &
                    + (lntaugrid(i2(x))-lntaugrid(i1(x))) &
                    / (skewgrid(i2(x))-skewgrid(i1(x))) &
                    * (skewness(x)-skewgrid(i1(x))) )
      b(:) = (pi-2._DP)/pi * (tau(:)-1._DP)**2 + tau(:)
      lambda(:) = sigma(:) / SQRT(b(:))
      mu(:) = mean(:) - sqrt2ovPi * lambda(:) * (tau(:)-1._DP)
    ELSE
      b(:) = (pi-2._DP)/pi * (tau(:)-1._DP)**2 + tau(:)
      mean(:) = mu(:) + sqrt2ovPi * lambda(:) * (tau(:)-1._DP)
      sigma(:) = SQRT(b(:)) * lambda(:)
      skewness(:) = b(:)**(-1.5_DP) * sqrt2ovPi * (tau(:)-1._DP) &
                    * ( (4/pi-1._DP)*(tau(:)-1._DP)**2 + tau(:) )
    END IF
    
    !-----------------------------------------------------------------------
    
  END SUBROUTINE param_splitnorm_1D

  ! 1D
  PURE SUBROUTINE param_splitnorm_0D (mu,lambda,tau,mean,sigma,skewness, &
                                      param2moment,moment2param)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    USE arrays, ONLY: ramp
    USE interpolation, ONLY: locate_sorted
    IMPLICIT NONE

    REAL(DP), INTENT(INOUT) :: mu, lambda, tau
    REAL(DP), INTENT(INOUT) :: mean, sigma, skewness
    LOGICAL, INTENT(IN), OPTIONAL :: param2moment, moment2param

    INTEGER, PARAMETER :: N = 10000
    INTEGER :: i1, i2
    REAL(DP) :: sqrt2ovPi
    REAL(DP), DIMENSION(N) :: taugrid, bgrid, skewgrid, lntaugrid
    REAL(DP) :: b
    LOGICAL :: param
    
    !-----------------------------------------------------------------------

    ! Determine which way to go: param -> moments or moments -> param
    IF (PRESENT(param2moment)) THEN
      param = ( .NOT. param2moment )
    ELSE IF (PRESENT(moment2param)) THEN
      param = moment2param
    ELSE
      param = .False.
    END IF

    ! Compute the parameters (compute tau numerically)
    sqrt2ovPi = SQRT(2._DP/pi)
    IF (param) THEN
      taugrid(:) = RAMP(N,1.E-2_DP,1.E2_DP,XLOG=.True.) ! tau=1.E2 => skew=0.99
      lntaugrid(:) = LOG(taugrid(:))
      bgrid(:) = (pi-2._DP)/pi * (taugrid(:)-1._DP)**2 + taugrid(:)
      skewgrid(:) = bgrid(:)**(-1.5_DP) * sqrt2ovPi * (taugrid(:)-1._DP) &
                  * ( (4/pi-1._DP)*(taugrid(:)-1._DP)**2 + taugrid(:) )
      i1 = LOCATE_SORTED(skewgrid(:),skewness)
      i2 = i1 + 1
      tau = EXP( lntaugrid(i1) &
               + (lntaugrid(i2)-lntaugrid(i1)) / (skewgrid(i2)-skewgrid(i1)) &
               * (skewness-skewgrid(i1)) )
      b = (pi-2._DP)/pi * (tau-1._DP)**2 + tau
      lambda = sigma / SQRT(b)
      mu = mean - sqrt2ovPi * lambda * (tau-1._DP)
    ELSE
      b = (pi-2._DP)/pi * (tau-1._DP)**2 + tau
      mean = mu + sqrt2ovPi * lambda * (tau-1._DP)
      sigma = SQRT(b) * lambda
      skewness = b**(-1.5_DP) * sqrt2ovPi * (tau-1._DP) &
               * ( (4/pi-1._DP)*(tau-1._DP)**2 + tau )
    END IF
    
    !-----------------------------------------------------------------------
    
  END SUBROUTINE param_splitnorm_0D


  !==========================================================================
  ! f[N] = DIST_LOGNORMAL(x[N],lnxmean,sigma)
  !
  !   Returns a normalised log-normal distribution centered in xmean, with a 
  ! standard deviation sigma.
  !==========================================================================

  PURE FUNCTION dist_lognormal_1D (x,lnxmean,sigma)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), INTENT(IN) :: lnxmean, sigma
    REAL(DP), DIMENSION(SIZE(x)) :: dist_lognormal_1D

    REAL(DP) :: norm, twice_sig2

    !-----------------------------------------------------------------------

    norm = 1._DP / (sigma*SQRT(2._DP*pi))
    twice_sig2 = 2._DP*sigma**2
    dist_lognormal_1D(:) = norm * EXP(-(LOG(x(:))-lnxmean)**2/twice_sig2) / x(:)

    !-----------------------------------------------------------------------

  END FUNCTION dist_lognormal_1D

  !----------------------------------------------------------------------------

  PURE FUNCTION dist_lognormal_scl (x,lnxmean,sigma)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: lnxmean, sigma
    REAL(DP) :: dist_lognormal_scl
    dist_lognormal_scl = MAXVAL(DIST_LOGNORMAL_1D([x],lnxmean,sigma))
  END FUNCTION dist_lognormal_scl


  !==========================================================================
  ! f[N] = DIST_LORENTZ(x[N],xmean,sigma)
  !
  !   Returns a normalised Lorentz distribution centered in xmean, with a 
  ! width WIDTH.
  !==========================================================================

  PURE FUNCTION dist_lorentz_1D (x,xmean,width)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), INTENT(IN) :: xmean, width
    REAL(DP), DIMENSION(SIZE(x)) :: dist_lorentz_1D

    REAL(DP) :: norm, quarter_wid2

    !-----------------------------------------------------------------------

    norm = width / (2._DP*pi)
    quarter_wid2 = width**2/4._DP
    dist_lorentz_1D = norm / ((x-xmean)**2 + quarter_wid2)

    !-----------------------------------------------------------------------

  END FUNCTION dist_lorentz_1D

  !----------------------------------------------------------------------------

  PURE FUNCTION dist_lorentz_scl (x,xmean,width)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: xmean, width
    REAL(DP) :: dist_lorentz_scl
    dist_lorentz_scl = MAXVAL(DIST_LORENTZ_1D([x],xmean,width))
  END FUNCTION dist_lorentz_scl


  !==========================================================================
  ! f[N] = DIST_SPLITLORENTZ(x[N],mu,lambda,tau)
  !
  !   Returns a normalised Split-Lorentz distribution centered in xmean, with a 
  ! width WIDTH and a shape parameter TAU.
  !==========================================================================

  PURE FUNCTION dist_splitlorentz_1D (x,mu,lambda,tau)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), INTENT(IN) :: mu, lambda, tau
    REAL(DP), DIMENSION(SIZE(x)) :: dist_splitlorentz_1D

    REAL(DP) :: norm

    !-----------------------------------------------------------------------

    norm = lambda / pi / ( 1._DP + tau )
    WHERE (x(:) < mu)
      dist_splitlorentz_1D(:) = norm / ((x(:)-mu)**2 + lambda**2/4._DP)
    ELSEWHERE
      dist_splitlorentz_1D(:) = norm*tau**2 / ((x(:)-mu)**2 &
                                              + lambda**2*tau**2/4._DP)
    END WHERE

    !-----------------------------------------------------------------------

  END FUNCTION dist_splitlorentz_1D

  !----------------------------------------------------------------------------

  PURE FUNCTION dist_splitlorentz_scl (x,mu,lambda,tau)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: mu, lambda, tau
    REAL(DP) :: dist_splitlorentz_scl
    dist_splitlorentz_scl = MAXVAL(DIST_SPLITLORENTZ_1D([x],mu,lambda,tau))
  END FUNCTION dist_splitlorentz_scl


  !==========================================================================
  ! CALL HISTOGRAM1D (x[N](IN),xbin[Nb](OUT),pbin[Nb](OUT),xlimits,
  !                   Nperbinmax(IN)])
  !
  !   Make an histogram of the array X. XBIN contains the center of each bin
  ! and PBIN the probability distribution (number-per-bin/N). The grid is
  ! regular and depends on the parameter Nperbinmax (30 by default). 
  !==========================================================================

  SUBROUTINE histogram1D (x,xbin,pbin,xlimits,Nperbinmax,Nbinmax,norm)

    USE utilities, ONLY: DP
    USE arrays, ONLY: sort, reallocate
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: xbin, pbin
    REAL(DP), DIMENSION(2), INTENT(IN), OPTIONAL :: xlimits
    INTEGER, INTENT(IN), OPTIONAL :: Nperbinmax, Nbinmax
    LOGICAL, INTENT(IN), OPTIONAL :: norm

    INTEGER :: i
    INTEGER :: N, Npbmax, Nmax, Nref, Nbin0, Nbin1
    INTEGER, DIMENSION(:), ALLOCATABLE :: Ndist0, Ndist1
    REAL(DP), DIMENSION(:), ALLOCATABLE :: xinf0, xinf1, xsup0, xsup1
    REAL(DP), DIMENSION(2) :: xrange
    LOGICAL, DIMENSION(:), ALLOCATABLE :: refine0, refine1, inbin
    LOGICAL :: normalize

    !-----------------------------------------------------------------------

    ! 1) Preparation
    !---------------
    N = SIZE(x)
    IF (PRESENT(xlimits)) THEN
      xrange(:) = xlimits(:)
      IF (xlimits(1) < MINVAL(x(:))) xrange(1) = MINVAL(x(:))
      IF (xlimits(2) > MAXVAL(x(:))) xrange(2) = MAXVAL(x(:))
    ELSE
      xrange(:) = [MINVAL(x(:)),MAXVAL(x(:))]
    END IF
    IF (PRESENT(norm)) THEN
      normalize = norm
    ELSE
      normalize = .False.
    END IF

    ! Maximum number of elements per bin
    IF (PRESENT(Nbinmax)) THEN 
      Nmax = Nbinmax 
    ELSE
      Nmax = 30000
    END IF
    IF (PRESENT(Nperbinmax)) THEN 
      Npbmax = Nperbinmax 
    ELSE
      Npbmax = 30
    END IF
    
    
    ! 2) Iteratively adpated bin grid
    !--------------------------------
    ! a. Initialization
    Nbin0 = 1
    ALLOCATE (refine0(Nbin0),Ndist0(Nbin0),xinf0(Nbin0),xsup0(Nbin0))
    refine0(:) = .True.
    xinf0(:) = xrange(1)
    xsup0(:) = xrange(2)
    Ndist0(:) = N
    Nref = 1
    
    ! b. Iteration
    ALLOCATE (inbin(N))
    adaptive: DO 
      ! Initialize
      Nbin1 = 2*Nbin0
      CALL REALLOCATE(refine1,Nbin1)
      CALL REALLOCATE(xinf1,Nbin1)
      CALL REALLOCATE(xsup1,Nbin1)
      CALL REALLOCATE(Ndist1,Nbin1)
      ! Refine
      xinf1(1:Nbin1-1:2) = xinf0(1:Nbin0)
      xinf1(2:Nbin1:2) = ( xsup0(1:Nbin0) + xinf0(1:Nbin0) ) / 2._DP
      xsup1(1:Nbin1-1) = xinf1(2:Nbin1)
      xsup1(Nbin1) = xsup0(Nbin0)
      DO i=1,Nbin1
        inbin(:) = ( x(:) >= xinf1(i) .AND. x(:) < xsup1(i) )
        Ndist1(i) = COUNT(inbin(:))
      END DO
      refine1(:) = ( Ndist1(:) >= Npbmax )
      ! Exit?
      Nref = COUNT(refine1(:))
      IF (Nref == 0 .OR. Nbin1 >= Nmax) THEN
        EXIT
      ELSE
        Nbin0 = Nbin1
        CALL REALLOCATE(refine0,Nbin0)
        CALL REALLOCATE(xinf0,Nbin0)
        CALL REALLOCATE(xsup0,Nbin0)
        CALL REALLOCATE(Ndist0,Nbin0)
        refine0(:) = refine1(:)
        xinf0(:) = xinf1(:)
        xsup0(:) = xsup1(:)
        Ndist0(:) = Ndist1(:)
      END IF
    END DO adaptive

    ! 3) Output
    !----------
    ALLOCATE (xbin(Nbin1))
    ALLOCATE (pbin(Nbin1))
    xbin(:) = ( xsup1(:) + xinf1(:) ) / 2._DP
    pbin(:) = REAL(Ndist1(:),DP) 
    IF (normalize) pbin(:) = pbin(:) / REAL(N,DP)

    ! Free memory space
    DEALLOCATE (Ndist0,Ndist1,xinf0,xinf1,xsup0,xsup1,refine0,refine1,inbin)

    !-----------------------------------------------------------------------

  END SUBROUTINE histogram1D


  !==========================================================================
  ! CALL HISTOGRAM1DMULTI (x[N,M](IN),xbin[Nb](OUT),pbin[Nb,M](OUT),xlimits,
  !                        Nperbinmax(IN)])
  !
  !   Make an histogram of the array X. Bin along N only, and keep the M
  ! dimension unchanged.
  !==========================================================================

  SUBROUTINE histogram1Dmulti (x,xbin,pbin,xlimits,Nperbinmax,Nbinmax,norm)

    USE utilities, ONLY: DP
    USE arrays, ONLY: sort, reallocate
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: xbin
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: pbin
    REAL(DP), DIMENSION(2), INTENT(IN), OPTIONAL :: xlimits
    INTEGER, INTENT(IN), OPTIONAL :: Nperbinmax, Nbinmax
    LOGICAL, INTENT(IN), OPTIONAL :: norm

    INTEGER :: i, j
    INTEGER :: Nx, Ny, Nmax, Npbmax, Nref, Nbin0, Nbin1
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: Ndist0, Ndist1
    REAL(DP), DIMENSION(:), ALLOCATABLE :: xinf0, xinf1, xsup0, xsup1
    REAL(DP), DIMENSION(2) :: xrange
    LOGICAL, DIMENSION(:), ALLOCATABLE :: refine0, refine1, inbin
    LOGICAL :: normalize

    !-----------------------------------------------------------------------

    ! 1) Preparation
    !---------------
    Nx = SIZE(x,1)
    Ny = SIZE(x,2)
    IF (PRESENT(xlimits)) THEN
      xrange(:) = xlimits(:)
      IF (xlimits(1) < MINVAL(x(:,:))) xrange(1) = MINVAL(x(:,:))
      IF (xlimits(2) > MAXVAL(x(:,:))) xrange(2) = MAXVAL(x(:,:))
    ELSE
      xrange(:) = [ MINVAL(x(:,:)), MAXVAL(x(:,:)) ]
    END IF
    IF (PRESENT(norm)) THEN
      normalize = norm
    ELSE
      normalize = .False.
    END IF

    ! Maximum number of elements per bin
    IF (PRESENT(Nbinmax)) THEN 
      Nmax = Nbinmax 
    ELSE
      Nmax = 30000
    END IF
    IF (PRESENT(Nperbinmax)) THEN 
      Npbmax = Nperbinmax 
    ELSE
      Npbmax = 30
    END IF
    
    ! 2) Iteratively adpated bin grid
    !--------------------------------
    ! a. Initialization
    Nbin0 = 1
    ALLOCATE (refine0(Nbin0),Ndist0(Nbin0,Ny),xinf0(Nbin0),xsup0(Nbin0))
    refine0(:) = .True.
    xinf0(:) = xrange(1)
    xsup0(:) = xrange(2)
    Ndist0(:,:) = Nx
    Nref = 1
    
    ! b. Iteration
    ALLOCATE (inbin(Nx))
    adaptive: DO 
      ! Initialize
      Nbin1 = 2*Nbin0
      CALL REALLOCATE(refine1,Nbin1)
      CALL REALLOCATE(xinf1,Nbin1)
      CALL REALLOCATE(xsup1,Nbin1)
      CALL REALLOCATE(Ndist1,Nbin1,Ny)
      ! Refine
      xinf1(1:Nbin1-1:2) = xinf0(1:Nbin0)
      xinf1(2:Nbin1:2) = ( xsup0(1:Nbin0) + xinf0(1:Nbin0) ) / 2._DP
      xsup1(1:Nbin1-1) = xinf1(2:Nbin1)
      xsup1(Nbin1) = xsup0(Nbin0)
      DO i=1,Nbin1
        DO j=1,Ny
          inbin(:) = ( x(:,j) >= xinf1(i) .AND. x(:,j) < xsup1(i) )
          Ndist1(i,j) = COUNT(inbin(:))
        END DO
        refine1(i) = ANY(Ndist1(i,:) >= Npbmax)
      END DO
      ! Exit?
      Nref = COUNT(refine1(:))
      IF (Nref == 0 .OR. Nbin1 >= Nmax) THEN
        EXIT
      ELSE
        Nbin0 = Nbin1
        CALL REALLOCATE(refine0,Nbin0)
        CALL REALLOCATE(xinf0,Nbin0)
        CALL REALLOCATE(xsup0,Nbin0)
        CALL REALLOCATE(Ndist0,Nbin0,Ny)
        refine0(:) = refine1(:)
        xinf0(:) = xinf1(:)
        xsup0(:) = xsup1(:)
        Ndist0(:,:) = Ndist1(:,:)
      END IF
    END DO adaptive

    ! 3) Output
    !----------
    ALLOCATE (xbin(Nbin1))
    ALLOCATE (pbin(Nbin1,Ny))
    xbin(:) = ( xsup1(:) + xinf1(:) ) / 2._DP
    pbin(:,:) = REAL(Ndist1(:,:),DP)
    IF (normalize) pbin(:,:) = pbin(:,:) / REAL(Nx,DP)

    ! Free memory space
    DEALLOCATE (Ndist0,Ndist1,xinf0,xinf1,xsup0,xsup1,refine0,refine1,inbin)

    !-----------------------------------------------------------------------

  END SUBROUTINE histogram1Dmulti


  !==========================================================================
  ! CALL HISTOGRAM2D (x[N](IN),y[N](IN),xbin[Nb](OUT),ybin[Nb](OUT), 
  !                   pbin[Nb,Nb](OUT),Nperbinmax(IN),xlimits(IN),ylimits(IN)])
  !
  !   Makes a 2D histogram of the arrays X and Y. It is similar to HISTOGRAM1D.
  !==========================================================================

  SUBROUTINE histogram2D (x,y,xbin,ybin,pbin,xlimits,ylimits,Nperbinmax, &
                          Nbinmax,norm)

    USE utilities, ONLY: DP, strike
    USE arrays, ONLY: reallocate
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: xbin, ybin
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: pbin
    REAL(DP), DIMENSION(2), INTENT(IN), OPTIONAL :: xlimits, ylimits
    INTEGER, INTENT(IN), OPTIONAL :: Nperbinmax, Nbinmax
    LOGICAL, INTENT(IN), OPTIONAL :: norm

    INTEGER :: ix, iy
    INTEGER :: N, Nmax, Npbmax, Nref, Nbin0, Nbin1
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: Ndist0, Ndist1
    REAL(DP), DIMENSION(:), ALLOCATABLE :: xinf0, xinf1, xsup0, xsup1
    REAL(DP), DIMENSION(:), ALLOCATABLE :: yinf0, yinf1, ysup0, ysup1
    REAL(DP), DIMENSION(2) :: xrange, yrange
    LOGICAL, DIMENSION(:), ALLOCATABLE :: inbin
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: refine0, refine1
    LOGICAL :: normalize

    !-----------------------------------------------------------------------

    ! 1) Preparation
    !---------------
    N = SIZE(x)
    IF (SIZE(y) /= N) CALL STRIKE ("HISTOGRAM2D","wrong dimensions of X and Y")
    IF (PRESENT(xlimits)) THEN
      xrange(:) = xlimits(:)
      IF (xlimits(1) < MINVAL(x(:))) xrange(1) = MINVAL(x(:))
      IF (xlimits(2) > MAXVAL(x(:))) xrange(2) = MAXVAL(x(:))
    ELSE
      xrange(:) = [MINVAL(x(:)),MAXVAL(x(:))]
    END IF
    IF (PRESENT(ylimits)) THEN
      yrange(:) = ylimits(:)
      IF (ylimits(1) < MINVAL(y(:))) yrange(1) = MINVAL(y(:))
      IF (ylimits(2) > MAXVAL(y(:))) yrange(2) = MAXVAL(y(:))
    ELSE
      yrange(:) = [MINVAL(y(:)),MAXVAL(y(:))]
    END IF
    IF (PRESENT(norm)) THEN
      normalize = norm
    ELSE
      normalize = .False.
    END IF

    ! Maximum number of elements per bin
    IF (PRESENT(Nbinmax)) THEN 
      Nmax = Nbinmax 
    ELSE
      Nmax = 30000
    END IF
    IF (PRESENT(Nperbinmax)) THEN 
      Npbmax = Nperbinmax 
    ELSE
      Npbmax = 30
    END IF
    
    ! 2) Iteratively adpated bin grid
    !--------------------------------
    ! a. Initialization
    Nbin0 = 1
    ALLOCATE ( refine0(Nbin0,Nbin0), Ndist0(Nbin0,Nbin0), &
               xinf0(Nbin0), xsup0(Nbin0), &
               yinf0(Nbin0), ysup0(Nbin0) )
    refine0(:,:) = .True.
    xinf0(:) = xrange(1)
    xsup0(:) = xrange(2)
    yinf0(:) = yrange(1)
    ysup0(:) = yrange(2)
    Ndist0(:,:) = N
    Nref = 1
    
    ! b. Iteration
    ALLOCATE (inbin(N))
    adaptive: DO 
      ! Initialize
      Nbin1 = 2*Nbin0
      CALL REALLOCATE(refine1,Nbin1,Nbin1)
      CALL REALLOCATE(xinf1,Nbin1)
      CALL REALLOCATE(xsup1,Nbin1)
      CALL REALLOCATE(yinf1,Nbin1)
      CALL REALLOCATE(ysup1,Nbin1)
      CALL REALLOCATE(Ndist1,Nbin1,Nbin1)
      ! Refine
      xinf1(1:Nbin1-1:2) = xinf0(1:Nbin0)
      xinf1(2:Nbin1:2) = ( xsup0(1:Nbin0) + xinf0(1:Nbin0) ) / 2._DP
      xsup1(1:Nbin1-1) = xinf1(2:Nbin1)
      xsup1(Nbin1) = xsup0(Nbin0)
      yinf1(1:Nbin1-1:2) = yinf0(1:Nbin0)
      yinf1(2:Nbin1:2) = ( ysup0(1:Nbin0) + yinf0(1:Nbin0) ) / 2._DP
      ysup1(1:Nbin1-1) = yinf1(2:Nbin1)
      ysup1(Nbin1) = ysup0(Nbin0)
      DO ix=1,Nbin1
        DO iy=1,Nbin1
          inbin(:) = ( x(:) >= xinf1(ix) .AND. x(:) < xsup1(ix) &
                     .AND. y(:) >= yinf1(iy) .AND. y(:) < ysup1(iy) )
          Ndist1(ix,iy) = COUNT(inbin(:))
        END DO
      END DO
      refine1(:,:) = ( Ndist1(:,:) >= Npbmax )
      ! Exit?
      Nref = COUNT(refine1(:,:))
      IF (Nref == 0 .OR. Nbin1 >= Nmax) THEN
        EXIT
      ELSE
        Nbin0 = Nbin1
        CALL REALLOCATE(refine0,Nbin0,Nbin0)
        CALL REALLOCATE(xinf0,Nbin0)
        CALL REALLOCATE(xsup0,Nbin0)
        CALL REALLOCATE(yinf0,Nbin0)
        CALL REALLOCATE(ysup0,Nbin0)
        CALL REALLOCATE(Ndist0,Nbin0,Nbin0)
        refine0(:,:) = refine1(:,:)
        xinf0(:) = xinf1(:)
        xsup0(:) = xsup1(:)
        yinf0(:) = yinf1(:)
        ysup0(:) = ysup1(:)
        Ndist0(:,:) = Ndist1(:,:)
      END IF
    END DO adaptive

    ! 3) Output
    !----------
    ALLOCATE (xbin(Nbin1))
    ALLOCATE (ybin(Nbin1))
    ALLOCATE (pbin(Nbin1,Nbin1))
    xbin(:) = ( xsup1(:) + xinf1(:) ) / 2._DP
    ybin(:) = ( ysup1(:) + yinf1(:) ) / 2._DP
    pbin(:,:) = REAL(Ndist1(:,:),DP)
    IF (normalize) pbin(:,:) = pbin(:,:) / REAL(N,DP)

    ! Free memory space
    DEALLOCATE (Ndist0,Ndist1,xinf0,xinf1,xsup0,xsup1,yinf0,yinf1,ysup0,ysup1, &
                refine0,refine1,inbin)

    !-----------------------------------------------------------------------

  END SUBROUTINE histogram2D


END MODULE distributions
