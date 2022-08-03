!******************************************************************************
!*
!*                     ROUTINES FOR LEAST-SQUARE METHODS
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 06/2014.
  ! 
  ! 3) DESCRIPTION: contains various methods performing chi2 minimization
  !                 for general multi-dimensional non-linear functions.
  !==========================================================================

MODULE chi2_minimization

  USE utilities, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: chi2min_LM


CONTAINS


  !==========================================================================
  ! CALL CHI2MIN_LM (funcresidual,Nobs,par,tol,resid,status,verbose, &
  !                  limited,limits,fixed,itied,chi2,chi2red,parerr,covar,Niter,
  !                  step,relstep,twoside)
  !
  !   Modified Levenberg-Marquardt chi2 minimization routine. Greatly based 
  ! on MINPACK's function LMDIF1 and dependencies, painfully converted from F77,
  ! and subsequently, but still painfully, partly vectorized. The Jacobian is 
  ! computed using finite difference approximation. Finally, I implemented 
  ! several functionalities of C. Markwardt's IDL MPFIT (limiting, fixing 
  ! parameters).
  !   Npar is the number of parameters par. Nobs is the number of observations. 
  ! PAR[Npar] is the parameter array. RESID[Nobs] is the result of the function.
  ! TOL is the termination tolerance. FUNCRESIDUAL is the function to be called.
  !   FUNCRESIDUAL is a function returning the weighted residuals between the
  ! model and the observations. It should have the following interface:
  !    INTERFACE
  !      FUNCTION funcresidual (par,Nobs)
  !        USE utilities, ONLY: DP
  !        IMPLICIT NONE
  !        INTEGER, INTENT(IN) :: Nobs
  !        REAL(DP), DIMENSION(:), INTENT(IN) :: par
  !        REAL(DP), DIMENSION(Nobs) :: funcresidual
  !      END FUNCTION funcresidual
  !    END INTERFACE
  ! STATUS is an integer output variable. If the user has terminated execution, 
  ! STATUS is set to the (negative) value of iflag. Otherwise, status is set as 
  ! follows:
  !         status = 0: improper input parameters.
  !         status = 1: algorithm estimates that the relative error
  !                     in the sum of squares is at most tol.
  !         status = 2: algorithm estimates that the relative error
  !                     between par and the solution is at most tol.
  !         status = 3: conditions for status = 1 and status = 2 both hold.
  !         status = 4: resid is orthogonal to the columns of the
  !                     jacobian to machine precision.
  !         status = 5: number of calls to funcresidual has reached or
  !                     exceeded 200*(n+1).
  !         status = 6: tol is too small. No further reduction in
  !                     the sum of squares is possible.
  !         status = 7: tol is too small. No further improvement in
  !                     the approximate solution par is possible.
  !==========================================================================

  SUBROUTINE chi2min_LM (funcresidual,Nobs,par,tol,resid,status,verbose, &
                         limited,limits,fixed,itied,parname,chi2,chi2red, &
                         parerr,covar,Niter,step,relstep,twoside,Nitermax)

    USE utilities, ONLY: DP, trimlr, pring, strike, verbatim
    USE inout, ONLY: lenpar
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nobs
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: par 
    INTEGER, INTENT(OUT), OPTIONAL :: status
    INTEGER, INTENT(OUT), OPTIONAL :: Niter
    INTEGER, DIMENSION(SIZE(par)), INTENT(IN), OPTIONAL :: itied
    REAL(DP), INTENT(IN), OPTIONAL :: tol
    REAL(DP), INTENT(OUT), OPTIONAL :: chi2, chi2red
    REAL(DP), DIMENSION(Nobs), INTENT(OUT), OPTIONAL :: resid 
    REAL(DP), DIMENSION(SIZE(par),2), INTENT(IN), OPTIONAL :: limits
    REAL(DP), DIMENSION(SIZE(par)), INTENT(OUT), OPTIONAL :: parerr
    REAL(DP), DIMENSION(SIZE(par),SIZE(par)), INTENT(OUT), OPTIONAL :: covar
    LOGICAL, INTENT(IN), OPTIONAL :: verbose
    LOGICAL, DIMENSION(SIZE(par),2), INTENT(IN), OPTIONAL :: limited
    LOGICAL, DIMENSION(SIZE(par)), INTENT(IN), OPTIONAL :: fixed
    CHARACTER(*), DIMENSION(SIZE(par)), INTENT(IN), OPTIONAL :: parname
    REAL(DP), DIMENSION(SIZE(par)), INTENT(IN), OPTIONAL :: step
    LOGICAL, DIMENSION(SIZE(par)), INTENT(IN), OPTIONAL :: relstep, twoside
    INTEGER, INTENT(IN), OPTIONAL :: Nitermax
    INTERFACE
      FUNCTION funcresidual (par,Nobs)
        USE utilities, ONLY: DP
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Nobs
        REAL(DP), DIMENSION(:), INTENT(IN) :: par
        REAL(DP), DIMENSION(Nobs) :: funcresidual
      END FUNCTION funcresidual
    END INTERFACE

    REAL(DP), PARAMETER :: factor = 1.E2_DP
    REAL(DP), PARAMETER :: tol0 = 1.E-10_DP
    REAL(DP), PARAMETER :: tolR = 1.E-14_DP
    LOGICAL, PARAMETER :: fastnorm = .False.

    INTEGER :: i, maxfev, mode, Nfev, Npar, Nparfree, info, iter
    INTEGER, DIMENSION(:), ALLOCATABLE :: ifree, ipvt
    INTEGER, DIMENSION(SIZE(par)) :: jtied
    REAL(DP) :: tol1, restol, gtol, partol
    REAL(DP), DIMENSION(Nobs) :: resid0
    REAL(DP), DIMENSION(SIZE(par),2) :: limpar
    REAL(DP), DIMENSION(:), ALLOCATABLE :: parfree, diag, diffstep
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: limitsfree, fjac, covarfree
    LOGICAL, DIMENSION(SIZE(par),2) :: limbool
    LOGICAL, DIMENSION(SIZE(par)) :: fix
    LOGICAL, DIMENSION(:), ALLOCATABLE :: reldiffstep, difftwoside
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: limitedfree
    LOGICAL :: printing
    CHARACTER(lenpar), DIMENSION(SIZE(par)) :: parnam


    !------------------------------------------------------------------------
    !                             Preparation
    !------------------------------------------------------------------------

    ! General parameters
    !-------------------
    Npar = SIZE(par(:))   ! Number of model parameters

    ! Fitter settings
    IF (.NOT. PRESENT(Nitermax)) THEN
      maxfev = 200 * (Npar + 1)     ! Maximum number of function evaluation
    ELSE
      maxfev = Nitermax * (Npar + 1)
    END IF
    tol1 = tol0
    IF (PRESENT(tol)) tol1 = tol  ! General tolerance
    restol = tol1     ! Tolerance on the residual
    partol = tol1     ! Tolerance on the parameter
    gtol = 0._DP      ! Tolerance on cos between residuals and jacobian
    mode = 1          ! 1 => internal scaling of parameters

    ! Options
    printing = verbatim            ! Compute function for printing
    IF (PRESENT(verbose)) printing = verbose
  
    ! Parameter name
    IF (PRESENT(parname)) THEN
      parnam(:) = parname(:)
    ELSE
      DO i=1,Npar 
        parnam(i) = "Par "//TRIMLR(PRING(i))
      END DO 
    END IF


    ! Parameter constraints
    !----------------------
    ! Tied parameters
    IF (PRESENT(itied)) THEN 
      IF (SIZE(itied) /= Npar) CALL STRIKE("CHI2MIN_LM","Wrong size for ITIED")
      jtied(:) = itied(:)

      ! Make sure that the parameter has been computed before tying its value:
      ! the index of the tied parameter should be larger than the index of the
      ! parameter it is tied to.
      DO i=1,Npar
        IF (jtied(i) > 0) THEN 
          IF (jtied(i) >= i) THEN
            jtied(jtied(i)) = i
            jtied(i) = 0
          ENDIF
        END IF
      END DO
    ELSE
      jtied(:) = 0
    END IF

    ! Fixed parameters
    fix(:) = .False.
    IF (PRESENT(fixed)) fix(:) = fixed(:)
    Nparfree = COUNT((.NOT. fix(:)) .AND. (jtied(:) <= 0))
    IF (Nobs < Nparfree) &
      CALL STRIKE ("CHI2MIN_LM","Not enough free parameters")
    ALLOCATE (parfree(Nparfree),ifree(Nparfree))
    parfree(:) = PACK(par(:),(.NOT. fix(:)) .AND. (jtied(:) <= 0))
    ifree(:) = PACK([(i,i=1,Npar)],(.NOT. fix(:)) .AND. (jtied(:) <= 0))

    ! Parameter limits
    limbool(:,:) = .False.
    IF (PRESENT(limited)) limbool(:,:) = limited(:,:)
    limpar(:,:) = 0._DP
    IF (PRESENT(limits)) limpar(:,:) = limits(:,:)

    ! Check consistency of limits and parameters
    IF (ANY(limbool(:,:))) THEN
      DO i=1,Npar
        IF (limbool(i,1) .AND. limbool(i,2)) THEN
          IF (limpar(i,1) > limpar(i,2)) &
            CALL STRIKE("CHI2MIN_LM","Limits of parameter "//TRIMLR(PRING(i)) &
                                   //" are inconsistent")
        END IF
        IF (limbool(i,1)) THEN
          IF (par(i) < limpar(i,1)) &
            CALL STRIKE ("CHI2MIN_LM","Parameter "//TRIMLR(PRING(i)) &
                                    //" is lower than its limit")
        END IF
        IF (limbool(i,2)) THEN
          IF (par(i) > limpar(i,2)) &
            CALL STRIKE ("CHI2MIN_LM","Parameter "//TRIMLR(PRING(i)) &
                                    //" is larger than its limit")
        END IF
      END DO
    END IF

    ! Effective limits on free parameters
    ALLOCATE (limitedfree(Nparfree,2),limitsfree(Nparfree,2))
    DO i=1,2
      limitedfree(:,i) = PACK(limbool(:,i),(.NOT. fix(:)) .AND. (jtied(:) <= 0))
      limitsfree(:,i) = PACK(limpar(:,i),(.NOT. fix(:)) .AND. (jtied(:) <= 0))
    END DO

    ! Optional scaling
    ALLOCATE(diag(Nparfree))
    diag(:) = 0._DP

    ! Optional stepsize for the Jacobian
    ALLOCATE (diffstep(Nparfree),reldiffstep(Nparfree))
    IF (PRESENT(step) .AND. PRESENT(relstep)) THEN
      diffstep(:) = PACK(step(:),(.NOT. fix(:)) .AND. (jtied(:) <= 0))
      reldiffstep(:) = PACK(relstep(:),(.NOT. fix(:)) .AND. (jtied(:) <= 0))
    ELSE
      diffstep(:) = 0._DP
      reldiffstep(:) = .True.
    END IF
    ALLOCATE (difftwoside(Nparfree))
    IF (PRESENT(twoside)) &
      difftwoside(:) = PACK(twoside(:),(.NOT. fix(:)) .AND. (jtied(:) <= 0))
       
    
    !------------------------------------------------------------------------
    !                        Actual System Solution
    !------------------------------------------------------------------------

    ! LMDIF
    ALLOCATE (fjac(Nobs,Nparfree),ipvt(Nparfree))
    CALL LMDIF(funcresidual,Nobs,Npar,Nparfree,par,parfree,ifree,jtied,resid0, &
               restol,partol,gtol,maxfev,diffstep,reldiffstep,difftwoside, &
               diag,mode,factor,printing,info,Nfev,fastnorm,limitedfree, &
               limitsfree,parnam,ipvt,fjac,Iter)

    ! Set (very carefully) the covariance matrix
    IF (PRESENT(covar)) covar(:,:) = 0._DP
    IF (PRESENT(parerr)) parerr(:) = 0._DP
    IF (info > 0 .AND. (PRESENT(covar) .OR. PRESENT(parerr))) THEN

      ! Covariance matrix for free parameters 
      ALLOCATE (covarfree(Nparfree,Nparfree))
      covarfree(:,:) = COVARMAT(fjac(:,:),Nobs,Nparfree,ipvt(:),tolR)
   
      ! Fill in actual covariance matrix, accounting for fixed parameters
      IF (PRESENT(covar)) THEN
        FORALL (i=1:Nparfree) covar(ifree(:),ifree(i)) = covarfree(:,i)
        DO i=1,Npar
          IF (jtied(i) > 0) THEN
            covar(i,:) = covar(jtied(i),:)
            covar(:,i) = covar(:,jtied(i))
          END IF
        END DO
      END IF
          
      ! Comput errors in parameters
      IF (PRESENT(parerr)) THEN
        FORALL (i=1:Nparfree) parerr(ifree(i)) = SQRT(covarfree(i,i))
        DO i=1,Npar
          IF (jtied(i) > 0) parerr(i) = parerr(jtied(i))
        END DO
      END IF

      ! Free memory space
      DEALLOCATE (covarfree)

    END IF

    ! Output status
    IF (info == 8) info = 4
    IF (printing) THEN
      PRINT*, "Status = "//TRIMLR(PRING(info))
      SELECT CASE (info)
        CASE (0)
          PRINT*, "Termination of CHI2MIN_LM: improper input parameters!"
        CASE (1)
          PRINT*, "Termination of CHI2MIN_LM: " &
                //"the relative error in the residuals is at most TOL."
        CASE (2)
          PRINT*, "Termination of CHI2MIN_LM: " &
                //"the relative error between the parameters and the " &
                //"solution is at most TOL."
        CASE (3)
          PRINT*, "Termination of CHI2MIN_LM: " &
                //"the relative errors in both the parameters and the " &
                //"residuals are at most TOL."
        CASE (4)
          PRINT*, "Termination of CHI2MIN_LM: " &
                //"the residuals are orthogonal to the columns of the " &
                //"jacobian to machine precision."
        CASE (5)
          PRINT*, "Termination of CHI2MIN_LM: " &
                //"number of calls to the function has reached or " &
                //"exceeded the maximum number."
        CASE (6)
          PRINT*, "Termination of CHI2MIN_LM: " &
                //"TOL is too small. No further reduction in the sum of " &
                //"squares is possible."
        CASE (7)
          PRINT*, "Termination of CHI2MIN_LM: " &
                //"TOL is too small. No further improvement in the " &
                //"approximate solution par is possible."
        CASE DEFAULT
          PRINT*, "Termination of CHI2MIN_LM: improper termination!"
      END SELECT
    END IF

    ! Optional ouput
    IF (PRESENT(status)) status = info
    IF (PRESENT(resid)) resid = resid0
    IF (PRESENT(chi2)) chi2 = SUM(resid0(:)**2)
    IF (PRESENT(chi2red)) chi2red = SUM(resid0(:)**2) / (Nobs-Nparfree)
    IF (PRESENT(Niter)) Niter = Iter

    ! Free memory space
    DEALLOCATE (ifree,ipvt,parfree,diag,limitsfree,fjac,reldiffstep, &
                difftwoside,limitedfree)

    !------------------------------------------------------------------------

  END SUBROUTINE chi2min_LM


  !--------------------------------------------------------------------------
  !   General Levenberg-Marquardt subroutine called by CHI2MIN_LM.
  !--------------------------------------------------------------------------

  SUBROUTINE lmdif (funcresidual,Nobs,Npar,Nparfree,par,parfree,ifree,itied, &
                    resid,restol,partol,gtol,maxfev,step,relstep,twoside,diag, &
                    mode,factor,printing,info,Nfev,fastnorm,limited,limits, &
                    parname,ipvt,fjac,Iter)
      
    USE utilities, ONLY: DP, epsDP, trimlr, pring
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nobs, Npar, Nparfree, maxfev
    INTEGER, INTENT(INOUT) :: mode
    INTEGER, INTENT(OUT) :: info, Nfev, Iter
    INTEGER, DIMENSION(Nparfree), INTENT(IN) :: ifree
    INTEGER, DIMENSION(Npar), INTENT(IN) :: itied
    INTEGER, DIMENSION(Nparfree), INTENT(OUT) :: ipvt
    REAL(DP), INTENT(IN) :: restol, partol, gtol, factor
    REAL(DP), DIMENSION(Nparfree), INTENT(IN) :: step
    REAL(DP), DIMENSION(Npar), INTENT(INOUT) :: par
    REAL(DP), DIMENSION(Nparfree), INTENT(INOUT) :: parfree, diag
    REAL(DP), DIMENSION(Nobs), INTENT(INOUT) :: resid
    REAL(DP), DIMENSION(Nparfree,2), INTENT(IN) :: limits
    REAL(DP), DIMENSION(Nobs,Nparfree), INTENT(OUT) :: fjac
    LOGICAL, DIMENSION(Nparfree,2), INTENT(IN) :: limited
    LOGICAL, DIMENSION(Nparfree), INTENT(IN) :: relstep, twoside
    LOGICAL, INTENT(IN) :: printing, fastnorm
    CHARACTER(*), DIMENSION(Npar), INTENT(IN) :: parname
    INTERFACE
      FUNCTION funcresidual (par,Nobs)
        USE utilities, ONLY: DP
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Nobs
        REAL(DP), DIMENSION(:), INTENT(IN) :: par
        REAL(DP), DIMENSION(Nobs) :: funcresidual
      END FUNCTION funcresidual
    END INTERFACE

    INTEGER :: i, iflag, j, l
    REAL(DP) :: actred, delta, dirder, fnorm, fnorm1, gnorm, alpha
    REAL(DP) :: parLM, pnorm, prered, ratio, temp, temp1, temp2, parnorm
    REAL(DP), DIMENSION(Npar) :: wa2par
    REAL(DP), DIMENSION(Nparfree) :: ulim1, llim1, qtf, wa1, wa2, wa3
    REAL(DP), DIMENSION(Nobs) :: wa4
    LOGICAL :: getout
    LOGICAL, DIMENSION(Nparfree) :: dwa1, limsup, liminf

    !------------------------------------------------------------------------

    info = 0
    iflag = 0
    Nfev = 0
    delta = 0._DP
    getout = .False.

    ! Check the input parameters for errors
    IF (Nparfree > 0 .AND. Nobs >= Nparfree &
        .AND. restol >= 0._DP .AND. partol >= 0._DP .AND. gtol >= 0._DP &
        .AND. maxfev > 0 .AND. factor > 0._DP) THEN

      IF (mode == 2) THEN
        IF (ANY(diag(:) <= 0._DP)) getout = .True.
      END IF

      ! Evaluate the function at the starting point and calculate its norm
      iflag = 1
      WHERE (itied(:) > 0) par(:) = par(itied(:))
      resid(:) = FUNCRESIDUAL(par(:),Nobs)
      Nfev = 1
      IF (iflag < 0) getout = .True.
      fnorm = ENORM(resid(:),fastnorm)

      ! Initialize levenberg-marquardt parameter and iteration counter
      parLM = 0._DP
      iter = 1
      parnorm = 1._DP

      ! Beginning of the outer loop
      outerloop: DO 
        IF (getout) EXIT

        ! Calculate the jacobian matrix
        iflag = 2
        CALL FDJAC2(funcresidual,Nobs,Npar,Nparfree,par,parfree,ifree,itied, &
                    resid,fjac,iflag,step,relstep,twoside,wa4,limited(1,:), &
                    limited(:,2),limits(1,:),limits(:,2))
        
        Nfev = Nfev + Nparfree
        IF (iflag < 0) EXIT

        ! Limits
        liminf(:) = ( limited(:,1) .AND. parfree(:) == limits(:,1) )
        limsup(:) = ( limited(:,2) .AND. parfree(:) == limits(:,2) )

        ! Total derivative of sum wrt lower/upper pegged parameters
        ! Note: SUM(resid*fjac(:,i)) is d(CHI^2)/dpar(i)
        IF (ANY(limited(:,1))) THEN
          DO i=1,Nparfree
            IF (liminf(i)) THEN
              IF (SUM(resid(:) * fjac(:,i)) > 0._DP) fjac(:,i) = 0._DP
            END IF
          END DO
        END IF
        IF (ANY(limited(:,2))) THEN
          DO i=1,Nparfree
            IF (limsup(i)) THEN
              IF (SUM(resid(:) * fjac(:,i)) < 0._DP) fjac(:,i) = 0._DP
            END IF
          END DO
        END IF

        ! Print iterates
        IF (printing) CALL PRINTPAR(iter,Nobs,Npar,Nparfree,par,resid, &
                                    limited,limits,parname,ifree,itied)

        ! Compute the QR factorization of the jacobian
        CALL QRFAC(Nobs,Nparfree,fjac,.True.,ipvt,Nparfree,wa1,wa2,wa3,fastnorm)
        
        ! On the first iteration and if mode is 1, scale according
        ! to the norms of the columns of the initial jacobian
        IF (iter == 1) THEN
          IF (mode /= 2) diag(:) = MERGE(wa2(:),1._DP,wa2(:) /= 0._DP)

          ! On the first iteration, calculate the norm of the scaled par
          ! and initialize the step bound delta
          wa3(:) = diag(:) * parfree(:)
          parnorm = ENORM(wa3(:),fastnorm)
          delta = factor * parnorm
          IF (delta == 0._DP) delta = factor
        END IF

        ! Form (Q transpose)*resid and store the first n components in qtf
        wa4(:) = resid(:)
        DO j=1,Nparfree
          IF (fjac(j,j) /= 0._DP) THEN
            temp = - SUM(fjac(j:Nobs,j)*wa4(j:Nobs)) / fjac(j,j)
            wa4(j:Nobs) = wa4(j:Nobs) + fjac(j:Nobs,j)*temp
          END IF
          fjac(j,j) = wa1(j)
          qtf(j) = wa4(j)
        END DO

        ! Compute the norm of the scaled gradient
        gnorm = 0._DP
        IF (fnorm /= 0._DP) THEN
          DO j=1,Nparfree
            l = ipvt(j)
            IF (wa2(l) /= 0._DP) &
              gnorm = MAX(gnorm,ABS(SUM(fjac(1:j,j)*qtf(1:j)/fnorm)/wa2(l)))
          END DO
        END IF

        ! Test for convergence of the gradient norm
        IF (gnorm <= gtol) info = 4
        IF (info /= 0) EXIT

        ! Rescale if necessary
        IF (mode /= 2) THEN
          DO j=1,Nparfree
            diag(j) = MAX(diag(j),wa2(j))
          END DO
        END IF

        ! Beginning of the inner loop
        innerloop: DO

          ! Determine the Levenberg-Marquardt parameter
          CALL LMPAR(Nparfree,fjac,Nobs,ipvt,diag,qtf,delta,parLM, &
                     wa1,wa2,wa3,wa4,fastnorm)

          ! Store the direction p and par + p. Calculate the norm of p
          wa1(:) = - wa1(:)
          IF (ANY(limited(:,:))) THEN

            ! Respect the limits. If a step were to go out of bounds, then
            ! we should take a step in the same direction but shorter distance.
            ! The step should take us right to the limit in that case.
            alpha = 1._DP

            ! Do not allow any steps out of bounds
            WHERE (liminf(:)) wa1(:) = MERGE(wa1(:),0._DP,wa1(:) > 0._DP)
            WHERE (limsup(:)) wa1(:) = MERGE(wa1(:),0._DP,wa1(:) < 0._DP)
            dwa1(:) = ( ABS(wa1(:)) > epsDP )
            alpha = MIN( alpha, &
                         MINVAL( (limits(:,1)-parfree(:))/wa1(:), &
                                 MASK=(dwa1(:) .AND. limited(:,1) &
                                      .AND. parfree(:)+wa1(:) < limits(:,1)) ) )
            alpha = MIN( alpha, &
                         MINVAL( (limits(:,2)-parfree(:))/wa1(:), &
                                 MASK=(dwa1(:) .AND. limited(:,2) &
                                      .AND. parfree(:)+wa1(:) > limits(:,2)) ) )

            ! Scale the resulting vector
            wa1(:) = wa1(:) * alpha
            wa2(:) = parfree(:) + wa1(:)

            ! Adjust the final output values. If the step put us exactly
            ! on a boundary, make sure we peg it there.
            ! Handles case of: nonzero*LIM +/- zero*LIM
            ulim1(:) = limits(:,2) * ( 1._DP - SIGN(1._DP,limits(:,2))*epsDP ) &
                     - MERGE(1._DP,0._DP,limits(:,2) == 0._DP) * epsDP
            llim1(:) = limits(:,1) * ( 1._DP + SIGN(1._DP,limits(:,1))*epsDP ) &
                     + MERGE(1._DP,0._DP,limits(:,1) == 0._DP) * epsDP
            wa2(:) = MERGE( limits(:,2), wa2(:), &
                            limited(:,2) .AND. wa2(:) >= ulim1(:) )
            wa2(:) = MERGE( limits(:,1), wa2(:), &
                            limited(:,1) .AND. wa2(:) <= llim1(:) )

          ELSE

            ! No parameter limits, so just move to new position WA2
            alpha = 1._DP
            wa2(:) = parfree(:) + wa1(:)

          END IF
          wa3(:) = diag(:) * wa1(:)
          pnorm = ENORM(wa3(:),fastnorm)

          ! On the first iteration, adjust the initial step bound
          IF (iter == 1) delta = MIN(delta,pnorm)

          ! Evaluate the function at par + p and calculate its norm
          iflag = 1
          wa2par(:) = par(:)
          wa2par(ifree(:)) = wa2(:)
          WHERE (itied(:) > 0) wa2par(:) = wa2par(itied(:))
          wa4(:) = FUNCRESIDUAL(wa2par(:),Nobs)
          Nfev = Nfev + 1
          IF (iflag < 0) THEN
            getout = .True.
            EXIT
          END IF
          fnorm1 = ENORM(wa4(:),fastnorm)
          
          ! Compute the scaled actual reduction
          actred = - 1._DP
          IF (0.1_DP*fnorm1 < fnorm) actred = 1._DP - (fnorm1/fnorm)**2
          
          ! Compute the scaled predicted reduction and the scaled directional 
          ! derivative
          DO j=1,Nparfree
            wa3(j) = 0._DP
            wa3(1:j) = wa3(1:j) + fjac(1:j,j)*wa1(ipvt(j))
          END DO
      
          ! Remember, alpha is the fraction of the full LM step actually taken
          temp1 = ENORM(alpha*wa3(:),fastnorm) / fnorm
          temp2 = SQRT(alpha*parLM) * pnorm / fnorm
          prered = temp1**2 + temp2**2/0.5_DP
          dirder = - (temp1**2 + temp2**2)

          ! Compute the ratio of the actual to the predicted reduction
          ratio = 0._DP
          IF (prered /= 0._DP) ratio = actred / prered

          ! Update the step bound
          IF (ratio <= 0.25_DP) THEN
            IF (actred >= 0._DP) temp = 0.5_DP
            IF (actred < 0._DP) &
              temp = 0.5_DP * dirder / (dirder + 0.5_DP*actred)
            IF (0.1_DP*fnorm1 >= fnorm .OR. temp < 0.1_DP) temp = 0.1_DP
            delta = temp * MIN(delta,pnorm/0.1_DP)
            parLM = parLM / temp
          ELSE
            IF (parLM == 0._DP .OR. ratio >= 0.75_DP) THEN
              delta = pnorm / 0.5_DP
              parLM = 0.5_DP * parLM
            END IF
          END IF

          ! Test for successful iteration
          IF (ratio >= 1.E-4_DP) THEN

            ! Successful iteration. Update par, resid, and their norms
            parfree(:) = wa2(:)
            par(ifree(:)) = parfree(:)
            WHERE (itied(:) > 0) par(:) = par(itied(:))
            wa2(:) = diag(:) * parfree(:)
            resid(:) = wa4(:)
            parnorm = ENORM(wa2(:),fastnorm)
            fnorm = fnorm1
            iter = iter + 1

          END IF

          ! Tests for convergence
          IF (ABS(actred) <= restol .AND. prered <= restol &
              .AND. 0.5_DP*ratio <= 1._DP) info = 1
          IF (delta <= partol*parnorm) info = 2
          IF (ABS(actred) <= restol .AND. prered <= restol &
              .AND. 0.5_DP*ratio <= 1._DP .AND. info == 2) info = 3
          IF (info /= 0) THEN
            getout = .True.
            EXIT
          END IF

          ! Tests for termination and stringent tolerances
          IF (Nfev >= maxfev) info = 5
          IF (ABS(actred) <= epsDP .AND. prered <= epsDP &
              .AND. 0.5_DP*ratio <= 1._DP) info = 6
          IF (delta <= epsDP*parnorm) info = 7
          IF (gnorm <= epsDP) info = 8
          IF (info /= 0) THEN
            getout = .True.
            EXIT
          END IF

          ! End of the inner loop. Repeat if iteration unsuccessful
          IF (ratio >= 1.E-4_DP) EXIT

       END DO innerloop
        
      END DO outerloop

    END IF

    ! Termination, either normal or user imposed
    IF (iflag < 0) info = iflag
    WHERE (itied(:) > 0) par(:) = par(itied(:))
    resid(:) = FUNCRESIDUAL(par(:),Nobs)
    IF (printing) CALL PRINTPAR(iter,Nobs,Npar,Nparfree,par,resid, &
                                limited,limits,parname,ifree,itied)
    RETURN

    !------------------------------------------------------------------------

  END SUBROUTINE lmdif 


  !--------------------------------------------------------------------------

  SUBROUTINE printpar (iter,Nobs,Npar,Nparfree,par,resid,limited,limits, &
                       parname,ifree,itied)

    USE utilities, ONLY: DP, trimlr, pring
    USE inout, ONLY: lenpar
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iter, Nobs, Npar, Nparfree
    INTEGER, DIMENSION(Nparfree), INTENT(IN) :: ifree
    INTEGER, DIMENSION(Npar), INTENT(IN) :: itied
    REAL(DP), DIMENSION(Npar), INTENT(IN) :: par
    REAL(DP), DIMENSION(Nobs), INTENT(IN) :: resid
    REAL(DP), DIMENSION(Nparfree,2), INTENT(IN) :: limits
    LOGICAL, DIMENSION(Nparfree,2), INTENT(IN) :: limited
    CHARACTER(*), DIMENSION(Npar), INTENT(IN) :: parname

    INTEGER, PARAMETER :: Ndec = 6

    INTEGER :: i, j
    CHARACTER(lenpar) :: parstr
    CHARACTER(lenpar*2) :: infostr

    !------------------------------------------------------------------------

    PRINT*, "  Iter "//TRIMLR(PRING(Iter))//" of CHI2MIN_LM" &
          //"       N(d.o.f.) = "//TRIMLR(PRING(Nobs-Nparfree)) &
          //"       chi2/N(d.o.f.) = " &
          //TRIMLR(PRING(SUM(resid(:)**2)/(Nobs-Nparfree),NDEC=Ndec))
    j = 1
    DO i=1,Npar
      IF (ALL(ifree(:) /= i)) THEN
        IF (itied(i) <= 0) THEN
          infostr = "fixed"
        ELSE
          infostr = "tied to "//TRIMLR(parname(itied(i)))
        END IF
      ELSE 
        IF (.NOT. limited(j,1) .AND. .NOT. limited(j,2)) THEN
          infostr = "free"
        ELSE IF (limited(j,1) .AND. .NOT. limited(j,2)) THEN
          infostr = "> "//TRIMLR(PRING(limits(j,1),NDEC=Ndec))
        ELSE IF (.NOT. limited(j,1) .AND. limited(j,2)) THEN
          infostr = "< "//TRIMLR(PRING(limits(j,2),NDEC=Ndec))
        ELSE IF (limited(j,1) .AND. limited(j,2)) THEN
          infostr = "> "//TRIMLR(PRING(limits(j,1),NDEC=Ndec))//" & < " &
                  //TRIMLR(PRING(limits(j,2),NDEC=Ndec))
        END IF
        j = j + 1
      END IF
      parstr = TRIMLR(PRING(par(i),NDEC=Ndec))
      PRINT*, ADJUSTR(parname(i))//"  =  "//ADJUSTL(parstr)//ADJUSTL(infostr)
    END DO                

    !------------------------------------------------------------------------

  END SUBROUTINE printpar


  !--------------------------------------------------------------------------
  !   Function calculating the Euclidian norm of a vector. It is computed by 
  ! accumulating the sum of squares in three different sums. The sums of 
  ! squares for the small and large components are scaled so that no overflows
  ! occur. Non-destructive underflows are permitted. Underflows and overflows 
  ! do not occur in the computation of the unscaled sum of squares for the 
  ! intermediate components. The definitions of small, intermediate and large 
  ! components depend on two constants, rdwarf and rgiant. The main
  ! restrictions on these constants are that rdwarf**2 not underflow and 
  ! rgiant**2 not overflow. The constants given here are suitable for every
  ! known computer.
  !--------------------------------------------------------------------------

  FUNCTION enorm (parfree,fastnorm)

    USE utilities, ONLY: DP, tinyDP, hugeDP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: parfree
    LOGICAL, INTENT(IN) :: fastnorm
    REAL(DP) :: enorm

    INTEGER :: i, Nparfree
    REAL(DP) :: agiant, adwarf, s1, s2, s3, parabs, par1max, par3max
    REAL(DP) :: rdwarf, rgiant

    !------------------------------------------------------------------------

    fastorslow: IF (fastnorm) THEN
  
      ! Fast but risky normalization
      !-----------------------------
      enorm = SQRT(SUM(parfree(:)**2))

    ELSE

      ! Numerically robust, but slower normalization
      !---------------------------------------------
      Nparfree = SIZE(parfree(:))
      rdwarf = SQRT(tinyDP) * 10._DP
      rgiant = SQRT(hugeDP) * 0.1_DP
      s1 = 0._DP
      s2 = 0._DP
      s3 = 0._DP
      par1max = 0._DP
      par3max = 0._DP
      agiant = rgiant / Nparfree
      adwarf = rdwarf * Nparfree
      DO i=1,Nparfree
        parabs = ABS(parfree(i))
        IF (parabs <= adwarf .OR. parabs >= agiant) THEN
          IF (parabs > adwarf) THEN

            ! Sum for large components
            IF (parabs > par1max) THEN
              s1 = 1._DP + s1*(par1max/parabs)**2
              par1max = parabs
            ELSE
              s1 = s1 + (parabs/par1max)**2
            END IF

          ELSE

            ! Sum for small components
            IF (parabs > par3max) THEN
              s3 = 1._DP + s3*(par3max/parabs)**2
              par3max = parabs
            ELSE
              IF (parabs /= 0._DP) s3 = s3 + (parabs/par3max)**2
            END IF
  
          END IF

        ELSE

          ! Sum for intermediate components
          s2 = s2 + parabs**2

        END IF
      END DO

      ! Calculation of norm
      IF (s1 /= 0._DP) THEN
        enorm = par1max * SQRT( s1 + (s2/par1max) / par1max )
      ELSE
        IF (s2 /= 0._DP) THEN
          IF (s2 >= par3max) &
            enorm = SQRT( s2 * (1._DP+(par3max/s2)*(par3max*s3)) )
          IF (s2 < par3max) &
            enorm = SQRT( par3max * ((s2/par3max)+(par3max*s3)) )
        ELSE
          enorm = par3max * SQRT(s3)
        END IF
      END IF
  
    END IF fastorslow
    
    RETURN

    !------------------------------------------------------------------------

  END FUNCTION enorm


  !--------------------------------------------------------------------------
  !   Computes a forward-difference approximation to the m by n jacobian 
  ! matrix associated with a specified problem of m functions in n variables.
  !--------------------------------------------------------------------------

  SUBROUTINE fdjac2 (funcresidual,Nobs,Npar,Nparfree,par,parfree,ifree,itied, &
                     resid,fjac,iflag,step,relstep,twoside,wa,downlimited, &
                     uplimited,downlimit,uplimit)

    USE utilities, ONLY: DP, epsDP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nobs, Nparfree, Npar
    INTEGER, INTENT(INOUT) :: iflag
    INTEGER, DIMENSION(Nparfree), INTENT(IN) :: ifree
    INTEGER, DIMENSION(Npar), INTENT(IN) :: itied
    REAL(DP), DIMENSION(Nparfree), INTENT(IN) :: step
    REAL(DP), DIMENSION(Npar), INTENT(INOUT) :: par
    REAL(DP), DIMENSION(Nparfree), INTENT(INOUT) :: parfree
    REAL(DP), DIMENSION(Nobs), INTENT(IN) :: resid
    REAL(DP), DIMENSION(Nobs,Nparfree), INTENT(OUT) :: fjac
    REAL(DP), DIMENSION(Nobs), INTENT(INOUT) :: wa
    REAL(DP), DIMENSION(Nparfree), INTENT(IN) :: downlimit, uplimit
    LOGICAL, DIMENSION(Nparfree), INTENT(IN) :: downlimited, uplimited
    LOGICAL, DIMENSION(Nparfree), INTENT(IN) :: relstep, twoside
    INTERFACE
      FUNCTION funcresidual (par,Nobs)
        USE utilities, ONLY: DP
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Nobs
        REAL(DP), DIMENSION(:), INTENT(IN) :: par
        REAL(DP), DIMENSION(Nobs) :: funcresidual
      END FUNCTION funcresidual
    END INTERFACE

    INTEGER :: j
    REAL(DP) :: eps, h, temp
    REAL(DP), DIMENSION(Nobs) :: wainf
    LOGICAL :: twosided
    
    !------------------------------------------------------------------------

    DO j=1,Nparfree
      twosided = twoside(j)
       
      ! Step size
      temp = parfree(j)
      eps = MAX(step(j),SQRT(epsDP))
      IF (relstep(j)) THEN 
        h = eps * ABS(temp)
        IF (h == 0._DP) h = eps
      ELSE
        h = eps
      END IF

      ! Adjust the limits
      IF (uplimited(j)) THEN
        IF (parfree(j) + h > uplimit(j)) THEN
          h = - h
          twosided = .False.
        END IF
      END IF
      IF (twosided) THEN
        IF (downlimited(j)) THEN
          IF (parfree(j) - h < downlimit(j)) twosided = .False.
        END IF
      END IF

      ! Compute the finite difference
      parfree(j) = temp + h
      par(ifree(j)) = parfree(j)
      WHERE (itied(:) > 0) par(:) = par(itied(:))
      wa(:) = FUNCRESIDUAL(par(:),Nobs)
      IF (iflag >= 0) THEN
        IF (twosided) THEN
          parfree(j) = temp - h
          par(ifree(j)) = parfree(j)
          WHERE (itied(:) > 0) par(:) = par(itied(:))
          wainf(:) = FUNCRESIDUAL(par(:),Nobs)          
          fjac(:,j) = (wa(:) - wainf(:)) / (2._DP*h)
        ELSE
          fjac(:,j) = (wa(:) - resid(:)) / h          
        END IF
        parfree(j) = temp
        par(ifree(j)) = parfree(j)      
        WHERE (itied(:) > 0) par(:) = par(itied(:))
      ELSE 
        EXIT
      END IF
    END DO

    RETURN

    !------------------------------------------------------------------------

  END SUBROUTINE fdjac2


  !--------------------------------------------------------------------------
  !   This subroutine uses householder transformations with column pivoting 
  ! (optional) to compute a QR factorization of the m by n matrix A. That is, 
  ! qrfac determines an orthogonal matrix Q, a permutation matrix P, and an 
  ! upper trapezoidal matrix R with diagonal elements of nonincreasing 
  ! magnitude, such that A*P = Q*R. The householder transformation for
  ! column k, k = 1,2,...,min(m,n), is of the form:
  !                         t
  !         i - (1/u(k))*u*u
  ! where u has zeros in the first k-1 positions. The form of this 
  ! transformation and the method of pivoting first appeared in the 
  ! corresponding linpack subroutine.
  !--------------------------------------------------------------------------

  SUBROUTINE qrfac (Nobs,Nparfree,A,pivot,ipvt,lipvt,rdiag,acnorm,wa,fastnorm)

    USE utilities, ONLY: DP, epsDP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nobs, Nparfree, lipvt
    INTEGER, DIMENSION(lipvt), INTENT(OUT) :: ipvt
    LOGICAL, INTENT(IN) :: pivot, fastnorm
    REAL(DP), DIMENSION(Nobs,Nparfree), INTENT(INOUT) :: A
    REAL(DP), DIMENSION(Nparfree), INTENT(OUT) :: rdiag, acnorm
    REAL(DP), DIMENSION(Nparfree), INTENT(INOUT) :: wa

    INTEGER :: i, j, jp1, k, kmax, minmn
    REAL(DP) :: ajnorm, temp

    !------------------------------------------------------------------------

    ! Compute the initial column norms and initialize several arrays
    DO j=1,Nparfree
      acnorm(j) = ENORM(A(:,j),fastnorm)
    END DO
    rdiag(:) = acnorm(:)
    wa(:) = rdiag(:)
    IF (pivot) ipvt(:) = [(j,j=1,Nparfree)]

    ! Reduce A to R with householder transformations
    minmn = MIN(Nobs,Nparfree)
    DO j=1,minmn

      IF (pivot) THEN
        ! Bring the column of largest norm into the pivot position
        kmax = MAXLOC(rdiag(j:Nparfree),DIM=1) + j - 1
        IF (kmax /= j) THEN
          DO i=1,Nobs
            temp = A(i,j)
            A(i,j) = A(i,kmax)
            A(i,kmax) = temp
          END DO
          rdiag(kmax) = rdiag(j)
          wa(kmax) = wa(j)
          k = ipvt(j)
          ipvt(j) = ipvt(kmax)
          ipvt(kmax) = k
        END IF
      END IF

      ! Compute the householder transformation to reduce the j-th column of a 
      ! to a multiple of the j-th unit vector
      ajnorm = ENORM(A(j:Nobs,j),fastnorm)
      IF (ajnorm /= 0._DP) THEN
        IF (A(j,j) < 0._DP) ajnorm = - ajnorm
        A(j:Nobs,j) = A(j:Nobs,j) / ajnorm
        A(j,j) = A(j,j) + 1._DP

        ! Apply the transformation to the remaining columns and update the norms
        jp1 = j + 1
        IF (Nparfree >= jp1) THEN
          DO k=jp1,Nparfree
            temp = SUM(A(j:Nobs,j)*A(j:Nobs,k)) / A(j,j)
            A(j:Nobs,k) = A(j:Nobs,k) - temp * A(j:Nobs,j)
            IF (pivot .AND. rdiag(k) /= 0._DP) THEN
              temp = A(j,k) / rdiag(k)
              rdiag(k) = rdiag(k) * SQRT(MAX(0._DP,1._DP-temp**2))
              IF (5.E-2_DP*(rdiag(k)/wa(k))**2 <= epsDP) THEN
                rdiag(k) = ENORM(A(jp1:Nobs,k),fastnorm)
                wa(k) = rdiag(k)
              END IF
            END IF
          END DO
        END IF
      END IF
      rdiag(j) = - ajnorm

    END DO

    RETURN

    !------------------------------------------------------------------------

  END SUBROUTINE qrfac


  !--------------------------------------------------------------------------
  !   Given an m by n matrix A, an n by n nonsingular diagonal matrix d, an 
  ! m-vector b, and A positive number delta, the problem is to determine a 
  ! value for the parameter parLM such that if par solves the system:
  !         A*par = b ,     sqrt(parLM)*d*par = 0 ,
  ! in the least squares sense, and dparnorm is the euclidean norm of d*par, 
  ! then either parLM is zero and:
  !         (dparnorm-delta) <= 0.1*delta ,
  ! or parLM is positive and:
  !         abs(dparnorm-delta) <= 0.1*delta .
  ! This subroutine completes the solution of the problem if it is provided 
  ! with the necessary information from the QR factorization, with column 
  ! pivoting, of A. That is, if A*P = Q*R, where P is a permutation matrix, Q 
  ! has orthogonal columns, and R is an upper triangular matrix with diagonal
  ! elements of nonincreasing magnitude, then lmpar expects the full upper 
  ! triangle of R, the permutation matrix P, and the first n components of 
  ! (Q transpose)*b. On output lmpar also provides an upper triangular matrix 
  ! S such that:
  !          t   t                   t
  !         P *(A *A + parLM*d*d)*P = S *S .
  ! S is employed within lmpar and may be of separate interest.
  ! Only a few iterations are generally needed for convergence of the 
  ! algorithm. If, however, the limit of 10 iterations is reached, then the 
  ! output parLM will contain the best value obtained so far.
  !--------------------------------------------------------------------------

  SUBROUTINE lmpar (Nparfree,R,Nobs,ipvt,diag,qtb,delta,parLM, &
                    parfree,sdiag,wa1,wa2,fastnorm)

    USE utilities, ONLY: DP, tinyDP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nparfree, Nobs
    INTEGER, DIMENSION(Nparfree), INTENT(IN) :: ipvt
    REAL(DP), INTENT(IN) :: delta
    REAL(DP), INTENT(INOUT) :: parLM
    REAL(DP), DIMENSION(Nobs,Nparfree), INTENT(INOUT) :: R
    REAL(DP), DIMENSION(Nparfree), INTENT(IN) :: diag, qtb
    REAL(DP), DIMENSION(Nparfree), INTENT(INOUT) :: parfree, sdiag
    REAL(DP), DIMENSION(Nparfree), INTENT(INOUT) :: wa1, wa2
    LOGICAL, INTENT(IN) :: fastnorm

    INTEGER :: iter, j, jm1, jp1, k, Nsing
    REAL(DP) :: dparnorm, fp, gnorm, parc, parl, paru, temp

    !------------------------------------------------------------------------

    ! Compute and store in par the gauss-newton direction. If the
    ! jacobian is rank-deficient, obtain a least squares solution.
    wa1(:) = qtb(:)
    Nsing = MINVAL(PACK([(j-1,j=1,Nparfree)],[(R(j,j),j=1,Nparfree)] == 0._DP, &
                        SPREAD(Nparfree,1,Nparfree)))
    IF (Nsing < Nparfree) wa1(Nsing+1:Nparfree) = 0._DP
    IF (Nsing >= 1) THEN
      DO k=1,Nsing
        j = Nsing - k + 1
        wa1(j) = wa1(j) / R(j,j)
        temp = wa1(j)
        jm1 = j - 1
        IF (jm1 >= 1) wa1(1:jm1) = wa1(1:jm1) - R(1:jm1,j) * temp
      END DO
    END IF
    parfree(ipvt(:)) = wa1(:)

    ! Initialize the iteration counter. Evaluate the function at the origin, 
    ! and test for acceptance of the gauss-newton direction.
    iter = 0
    wa2(:) = diag(:) * parfree(:)
    dparnorm = ENORM(wa2(:),fastnorm)
    fp = dparnorm - delta
    IF (fp > 0.1_DP*delta) THEN

      ! If the jacobian is not rank deficient, the Newton step provides a 
      ! lower bound, parl, for the zero of the function. Otherwise set this 
      ! bound to zero.
      parl = 0._DP
      IF (Nsing >= Nparfree) THEN
        wa1(:) = diag(ipvt(:)) * wa2(ipvt(:)) / dparnorm
        wa1(1) = wa1(1) / R(1,1)
        DO j=2,Nparfree
          jm1 = j - 1
          wa1(j) = (wa1(j) - SUM(R(1:jm1,j)*wa1(1:jm1)) ) / R(j,j)
        END DO
        temp = ENORM(wa1(:),fastnorm)
        parl = ( (fp/delta) / temp ) / temp
      END IF

      ! Calculate an upper bound, paru, for the zero of the function.
      FORALL (j=1:Nparfree) wa1(j) = SUM(R(1:j,j)*qtb(1:j)) / diag(ipvt(j))
      gnorm = ENORM(wa1(:),fastnorm)
      paru = gnorm / delta
      IF (paru == 0._DP) paru = tinyDP / MIN(delta,0.1_DP)

      ! If the input parLM lies outside of the interval (parl,paru),
      ! set parLM to the closer endpoint.
      parLM = MAX(parLM,parl)
      parLM = MIN(parLM,paru)
      IF (parLM == 0._DP) parLM = gnorm / dparnorm

      ! Beginning of an iteration
      iteration: DO 
         iter = iter + 1

         ! Evaluate the function at the current value of parLM
         IF (parLM == 0._DP) parLM = MAX(tinyDP,1.E-3_DP*paru)
         wa1(:) = SQRT(parLM) * diag(:)
         CALL QRSOLV(Nparfree,r,Nobs,ipvt,wa1,qtb,parfree,sdiag,wa2)
         wa2(:) = diag(:) * parfree(:)
         dparnorm = ENORM(wa2(:),fastnorm)
         temp = fp
         fp = dparnorm - delta

         ! If the function is small enough, accept the current value of parLM. 
         ! Also test for the exceptional cases where parl is zero or the 
         ! number of iterations has reached 10.
         IF (ABS(fp) <= 0.1_DP*delta .OR. parl == 0._DP .AND. fp <= temp &
             .AND. temp < 0._DP .OR. iter == 10) EXIT

         ! Compute the newton correction
         wa1(:) = diag(ipvt(:)) * wa2(ipvt(:)) / dparnorm
         DO j=1,Nparfree-1
           wa1(j) = wa1(j) / sdiag(j)
           jp1 = j + 1
           wa1(jp1:Nparfree) = wa1(jp1:Nparfree) - R(jp1:Nparfree,j) * wa1(j)
         END DO
         wa1(Nparfree) = wa1(Nparfree) / sdiag(Nparfree)
         temp = ENORM(wa1(:),fastnorm)
         parc = ( (fp/delta) / temp ) / temp

         ! Depending on the sign of the function, update parl or paru
         IF (fp > 0._DP) parl = MAX(parl,parLM)
         IF (fp < 0._DP) paru = MIN(paru,parLM)

         ! Compute an improved estimate for 
         parLM = MAX(parl,parLM+parc)

      END DO iteration

    END IF

    ! Termination
    IF (iter == 0) parLM = 0._DP
      
    RETURN

    !------------------------------------------------------------------------

  END SUBROUTINE lmpar


  !--------------------------------------------------------------------------
  !   Given an m by n matrix A, an n by n diagonal matrix D, and an m-vector b,
  ! the problem is to determine an par which solves the system:
  !         A*par = b ,     D*par = 0 ,
  ! in the least squares sense. This subroutine completes the solution of the 
  ! problem if it is provided with the necessary information from the QR 
  ! factorization, with column pivoting, of A. That is, if A*P = Q*R, where P 
  ! is a permutation matrix, Q has orthogonal columns, and R is an upper 
  ! triangular matrix with diagonal elements of nonincreasing magnitude, then 
  ! qrsolv expects the full upper triangle of R, the permutation matrix P, and 
  ! the first n components of (Q transpose)*b. the system
  !   A*par = b, D*par = 0, is then equivalent to:
  !                t       t
  !         R*z = Q *b ,  P *D*P*z = 0 ,
  ! where par = P*z. If this system does not have full rank, then a least 
  ! squares solution is obtained. On output qrsolv also provides an upper 
  ! triangular matrix s such that:
  !          t   t               t
  !         P *(A *A + D*D)*P = S *S .
  ! S is computed within qrsolv and may be of separate interest.
  !--------------------------------------------------------------------------

  SUBROUTINE qrsolv (Nparfree,r,Nobs,ipvt,diag,qtb,parfree, &
                     sdiag,wa)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nparfree, Nobs
    INTEGER, DIMENSION(Nparfree), INTENT(IN) :: ipvt
    REAL(DP), DIMENSION(Nobs,Nparfree), INTENT(INOUT) :: R
    REAL(DP), DIMENSION(Nparfree), INTENT(IN) :: diag, qtb
    REAL(DP), DIMENSION(Nparfree), INTENT(INOUT) :: parfree, sdiag
    REAL(DP), DIMENSION(Nparfree), INTENT(INOUT) :: wa
 
    INTEGER :: j, jp1, k, kp1, l, Nsing
    REAL(DP) :: cosine, cotan, qtbpj, sine, tang, temp
    REAL(DP), DIMENSION(Nparfree) :: tempvec

    !------------------------------------------------------------------------

    ! Copy R and (Q transpose)*b to preserve input and initialize S. In 
    ! particular, save the diagonal elements of R in par.
    DO j=1,Nparfree
      R(j:Nparfree,j) = R(j,j:Nparfree)
      parfree(j) = R(j,j)
    END DO
    wa(:) = qtb(:)

    ! Eliminate the diagonal matrix D using a givens rotation
    DO j=1,Nparfree

      ! Prepare the row of D to be eliminated, locating the diagonal element 
      ! using P from the QR factorization
      l = ipvt(j)
      IF (diag(l) /= 0._DP) THEN

        sdiag(j:Nparfree) = 0._DP
        sdiag(j) = diag(l)

        ! The transformations to eliminate the row of D modify only a single 
        ! element of (Q transpose)*b beyond the first n, which is initially 
        ! zero
        qtbpj = 0._DP
        DO k=j,Nparfree

          ! Determine a givens rotation which eliminates the appropriate 
          ! element in the current row of D
          IF (sdiag(k) /= 0._DP) THEN

            IF (ABS(R(k,k)) < ABS(sdiag(k))) THEN
              cotan = R(k,k) / sdiag(k)
              sine = 0.5_DP / SQRT(0.25_DP+0.25_DP*cotan**2)
              cosine = sine * cotan
            ELSE
              tang = sdiag(k) / R(k,k)
              cosine = 0.5_DP / SQRT(0.25_DP+0.25_DP*tang**2)
              sine = cosine * tang
            END IF

            ! Compute the modified diagonal element of R and the modified 
            ! element of ((Q transpose)*b,0)
            R(k,k) = cosine * R(k,k) + sine * sdiag(k)
            temp = cosine * wa(k) + sine * qtbpj
            qtbpj = - sine * wa(k) + cosine * qtbpj
            wa(k) = temp

            ! Accumulate the tranformation in the row of S
            kp1 = k + 1
            IF (Nparfree >= kp1) THEN
              tempvec(kp1:Nparfree) = cosine*R(kp1:Nparfree,k) &
                                    + sine*sdiag(kp1:Nparfree)
              sdiag(kp1:Nparfree) = - sine*R(kp1:Nparfree,k) &
                                    + cosine*sdiag(kp1:Nparfree)
              R(kp1:Nparfree,k) = tempvec(kp1:Nparfree)
            END IF
          END IF

        END DO

      END IF   

      ! Store the diagonal element of S and restore the corresponding diagonal
      ! element of R
      sdiag(j) = R(j,j)
      R(j,j) = parfree(j)

    END DO

    ! Solve the triangular system for z. If the system is singular, then obtain
    ! a least squares solution
    Nsing = MINVAL(PACK([(j-1,j=1,Nparfree)],sdiag(:) == 0._DP, &
                   SPREAD(Nparfree,1,Nparfree)))
    IF (Nsing < Nparfree) wa(Nsing+1:Nparfree) = 0._DP
    IF (Nsing >= 1) THEN
      DO k=1,Nsing
        j = Nsing - k + 1
        jp1 = j + 1
        wa(j) = (wa(j) &
                - MERGE(SUM(R(jp1:Nsing,j)*wa(jp1:Nsing)),0._DP,Nsing >= jp1)) &
              / sdiag(j)
      END DO
    END IF

    ! Permute the components of z back to components of par
    parfree(ipvt(:)) = wa(:)

    RETURN

    !------------------------------------------------------------------------

  END SUBROUTINE qrsolv


  !--------------------------------------------------------------------------
  !   Given an m by n matrix A, the problem is to determine the covariance 
  ! matrix corresponding to A, defined as:
  !          t
  ! inverse(A *A) .
  ! This subroutine completes the solution of the problem if it is provided 
  ! with the necessary information from the QR factorization, with column 
  ! pivoting, of A. that is, if
  ! A*P = Q*R, where P is a permutation matrix, Q has orthogonal columns, and 
  ! R is an upper triangular matrix with diagonal elements of nonincreasing 
  ! magnitude, then covar expects the full upper triangle of R and the 
  ! permutation matrix P. The covariance matrix is then computed as:
  !  t          t
  ! P *inverse(R *R)*P .
  ! If A is nearly rank deficient, it may be desirable to compute the 
  ! covariance matrix corresponding to the linearly independent columns of A. 
  ! To define the numerical rank of A, covar uses the tolerance tol. If l is 
  ! the largest integer such that:
  ! abs(R(l,l)) .gt. tol*abs(R(1,1)) ,
  ! then covar computes the covariance matrix corresponding to the first l 
  ! columns of R. For k greater than l, column c and row ipvt(k) of the 
  ! covariance matrix are set to zero.
  !--------------------------------------------------------------------------

  FUNCTION covarmat (R,Nobs,Nparfree,ipvt,tol)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nobs, Nparfree
    INTEGER, DIMENSION(Nparfree), INTENT(IN) :: ipvt
    REAL(DP), DIMENSION(Nobs,Nparfree), INTENT(IN) :: R
    REAL(DP), INTENT(IN) :: tol
    REAL(DP), DIMENSION(Nparfree,Nparfree) :: covarmat

    INTEGER :: i, ii, j, jj, k, km1, l
    LOGICAL :: sing
    REAL(DP) :: temp, tolr
    REAL(DP), DIMENSION(Nparfree) :: wa
    REAL(DP), DIMENSION(Nobs,Nparfree) :: cov

    !------------------------------------------------------------------------

    cov(:,:) = R(:,:)

    ! Form the inverse of r in the full upper triangle of r
    tolr = tol * ABS(cov(1,1))
    l = 0
    DO k=1,Nparfree
      IF (ABS(cov(k,k)) <= tolr) EXIT
      cov(k,k) = 1._DP / cov(k,k)
      km1 = k - 1
      IF (km1 >= 1) THEN
        DO j=1,km1
          temp = cov(k,k) * cov(j,k)
          cov(j,k) = 0._DP
          cov(1:j,k) = cov(1:j,k) - temp * cov(1:j,j)
        END DO
      END IF
      l = k
    END DO

    ! Form the full upper triangle of the inverse of (R transpose)*R in the 
    ! full upper triangle of R
    IF (l >= 1) THEN
      DO k=1,l
        km1 = k - 1
        IF (km1 >= 1) THEN
          DO j=1,km1
            temp = cov(j,k)
            cov(1:j,j) = cov(1:j,j) + temp * cov(1:j,k)
          END DO
        END IF
        temp = cov(k,k)
        cov(1:k,k) = temp * cov(1:k,k)
      END DO
    END IF

    ! Form the full lower triangle of the covariance matrix in the strict 
    ! lower triangle of R and in wa
    DO j=1,Nparfree
      jj = ipvt(j)
      sing = (j > l)
      DO i=1,j
        IF (sing) cov(i,j) = 0._DP
        ii = ipvt(i)
        IF (ii > jj) cov(ii,jj) = cov(i,j)
        IF (ii < jj) cov(jj,ii) = cov(i,j)
      END DO
      wa(jj) = cov(j,j)
    END DO

    ! Symmetrize the covariance matrix in R
    DO j=1,Nparfree
      cov(1:j,j) = cov(j,1:j)
      cov(j,j) = wa(j)
    END DO
    covarmat(:,:) = cov(1:Nparfree,1:Nparfree)

    RETURN

    !------------------------------------------------------------------------

  END FUNCTION covarmat


END MODULE chi2_minimization
