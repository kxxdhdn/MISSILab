!******************************************************************************
!*
!*                              Chi2 Fitting Module
!*
!******************************************************************************


MODULE chi2

  USE utilities, ONLY: DP
  USE core, ONLY: parinfo_type, indpar_type, Qabs_type
  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PUBLIC :: Nx, Ny, NwOBS, xOBS, yOBS, jw
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: iwfree
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wOBS, nuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Fnu_mod, resid
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC :: extinct
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC :: FnuOBS, dFnuOBS
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE, PUBLIC :: invLcovarOBS
  LOGICAL, SAVE, PUBLIC :: iid ! indpdt identically dist. (wvl)
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC :: mask

  !! Init param
  !!------------
  TYPE(indpar_type), SAVE, PUBLIC :: ind
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parinfo
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Qabs
  
  PUBLIC :: residuals

CONTAINS

  !! Main function for L-M
  !!-----------------------
  FUNCTION residuals(par, NwOBS0)
    
    USE utilities, ONLY: DP
    USE core, ONLY: specModel
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NwOBS0
    REAL(DP), DIMENSION(:), INTENT(IN) :: par
    INTEGER :: i
    
    REAL(DP), DIMENSION(NwOBS0) :: residuals

    Fnu_mod(:) = specModel(wOBS(:), INDPAR=ind, PARVAL=par(:), &
                           QABS=Qabs(:), EXTINCT=extinct(:,:))

    !! Unweighted residuals    
    WHERE (mask(xOBS,yOBS,:))
      resid(:) = FnuOBS(xOBS,yOBS,:) - Fnu_mod(:)
    ELSEWHERE
      resid(:) = 0._DP

    END WHERE
    
    !! Weighted residuals
    ! residuals(:) = MATMUL(invLcovarOBS(xOBS,yOBS,:,:), resid(:))
    
    !! With keyword iid
    IF (iid) THEN
      FORALL (i=1:NwOBS0) &
        residuals(i) = invLcovarOBS(xOBS,yOBS,i,i) * resid(i)
    ELSE
      residuals(:) = MATMUL(invLcovarOBS(xOBS,yOBS,:,:), resid(:))

    END IF

  END FUNCTION residuals

END MODULE chi2
