!******************************************************************************
!*
!*                             BASIC STATISTICS
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 08/2007.
  !    - 09/2013: add autocorrel.
  !    - 08/2014: transfer Rmat2corr.
  !    - 05/2015: transfer corr2Rmat.
  !    - 03/2016: add MASK for SIGMA and MEDIAN.
  !    - 07/2016: add the computation of the integrated autocorrelation time
  !               and of the effective sample size of an MCMC.
  ! 
  ! 3) DESCRIPTION: computes the basic statistical quantities for a given
  !                 data set.
  !==========================================================================

MODULE statistics

  USE utilities, ONLY:
  USE arrays, ONLY:
  USE special_functions, ONLY:
  USE FFT_specials, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: moment_data, median_data, median, median_conf, info_data, fwhm_data
  PUBLIC :: mean, sigma, Rmat2corr, corr2Rmat, Vmat2corr, N_corr, correlate
  PUBLIC :: autocorrel, correl_parlist, correl_index, intautocorrtime

  INTERFACE mean
    MODULE PROCEDURE mean_1D, mean_2D, mean_3D, mean_4D
    MODULE PROCEDURE mean_4D_dim, mean_3D_dim, mean_2D_dim
  END INTERFACE mean

  INTERFACE median ! alias
    MODULE PROCEDURE median_data_1D, median_data_2D, median_data_3D
    MODULE PROCEDURE median_data_4D
    MODULE PROCEDURE median_data_2D_dim, median_data_3D_dim, median_data_4D_dim
    MODULE PROCEDURE median_data_wei
  END INTERFACE median

  INTERFACE median_data
    MODULE PROCEDURE median_data_1D, median_data_2D, median_data_3D
    MODULE PROCEDURE median_data_4D
    MODULE PROCEDURE median_data_2D_dim, median_data_3D_dim, median_data_4D_dim
    MODULE PROCEDURE median_data_wei
  END INTERFACE median_data

  INTERFACE median_conf
    MODULE PROCEDURE median_conf_1D, median_conf_2D, median_conf_3D
    MODULE PROCEDURE median_conf_4D
  END INTERFACE median_conf

  INTERFACE sigma
    MODULE PROCEDURE sigma_1D, sigma_2D, sigma_3D, sigma_4D
    MODULE PROCEDURE sigma_2D_dim, sigma_3D_dim, sigma_4D_dim
  END INTERFACE sigma

  INTERFACE correlate
    MODULE PROCEDURE correlate_1D, correlate_2D, correlate_3D, correlate_4D
  END INTERFACE correlate

  INTERFACE info_data
    MODULE PROCEDURE info_data_1D, info_data_2D
  END INTERFACE info_data

  INTERFACE Rmat2corr
    MODULE PROCEDURE Rmat2corr_DP, Rmat2corr_log, Rmat2corr_char
  END INTERFACE Rmat2corr

  INTERFACE corr2Rmat
    MODULE PROCEDURE corr2Rmat_DP
  END INTERFACE corr2Rmat

  INTERFACE correl_parlist
    MODULE PROCEDURE correl_parlist_DP, correl_parlist_char
  END INTERFACE correl_parlist

  LOGICAL, PARAMETER :: upper = .True.


CONTAINS


  !==========================================================================
  ! m(X) = MEAN(X[],weight[],mask[],dim)
  !
  !   Computes the mean of the variable X, with weights w. If w is not 
  ! supplied then the distribution is assumed to be uniform.
  !==========================================================================

  PURE FUNCTION mean_4D (x,weight,mask)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)), INTENT(IN), &
      OPTIONAL :: weight
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)), INTENT(IN), &
      OPTIONAL :: mask
    REAL(DP) :: mean_4D

    !-----------------------------------------------------------------------

    IF (PRESENT(weight)) THEN 
      IF (PRESENT(mask)) THEN
        mean_4D = SUM(x(:,:,:,:)*weight(:,:,:,:),MASK=mask(:,:,:,:)) &
                / SUM(weight(:,:,:,:),MASK=mask(:,:,:,:))
      ELSE
        mean_4D = SUM(x(:,:,:,:)*weight(:,:,:,:)) / SUM(weight(:,:,:,:))
      END IF
    ELSE 
      IF (PRESENT(mask)) THEN
        mean_4D = SUM(x(:,:,:,:),MASK=mask(:,:,:,:)) / COUNT(mask(:,:,:,:))
      ELSE
        mean_4D = SUM(x(:,:,:,:)) / SIZE(x(:,:,:,:))
      END IF
    END IF

    !-----------------------------------------------------------------------

  END FUNCTION mean_4D


  ! Shell Interfaces
  !-----------------
  ! 3D
  PURE FUNCTION mean_3D (x,weight,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)), INTENT(IN), &
      OPTIONAL :: weight
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)), INTENT(IN), &
      OPTIONAL :: mask
    REAL(DP) :: mean_3D
    INTEGER, DIMENSION(4) :: N 
    N = [SIZE(x,1),SIZE(x,2),SIZE(x,3),1]
    IF (.NOT. PRESENT(weight)) THEN
      IF (.NOT. PRESENT(mask)) THEN
        mean_3D = MEAN_4D(RESHAPE(x(:,:,:),N(:)))
      ELSE
        mean_3D = MEAN_4D(RESHAPE(x(:,:,:),N(:)),MASK=RESHAPE(mask(:,:,:),N(:)))
      END IF
    ELSE
      IF (.NOT. PRESENT(mask)) THEN
        mean_3D = MEAN_4D(RESHAPE(x(:,:,:),N(:)), &
                          WEIGHT=RESHAPE(weight(:,:,:),N(:)))
      ELSE
        mean_3D = MEAN_4D(RESHAPE(x(:,:,:),N(:)), &
                          WEIGHT=RESHAPE(weight(:,:,:),N(:)), &
                          MASK=RESHAPE(mask(:,:,:),N(:)))
      END IF
    END IF
  END FUNCTION mean_3D

  ! 2D
  PURE FUNCTION mean_2D (x,weight,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x,1),SIZE(x,2)), INTENT(IN), OPTIONAL :: weight
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2)), INTENT(IN), OPTIONAL :: mask
    REAL(DP) :: mean_2D
    INTEGER, DIMENSION(4) :: N 
    N = [SIZE(x,1),SIZE(x,2),1,1]
    IF (.NOT. PRESENT(weight)) THEN
      IF (.NOT. PRESENT(mask)) THEN
        mean_2D = MEAN_4D(RESHAPE(x(:,:),N(:)))
      ELSE
        mean_2D = MEAN_4D(RESHAPE(x(:,:),N(:)),MASK=RESHAPE(mask(:,:),N(:)))
      END IF
    ELSE
      IF (.NOT. PRESENT(mask)) THEN
        mean_2D = MEAN_4D(RESHAPE(x(:,:),N(:)),WEIGHT=RESHAPE(weight(:,:),N(:)))
      ELSE
        mean_2D = MEAN_4D(RESHAPE(x(:,:),N(:)), &
                          WEIGHT=RESHAPE(weight(:,:),N(:)), &
                          MASK=RESHAPE(mask(:,:),N(:)))
      END IF
    END IF
  END FUNCTION mean_2D

  ! 1D
  PURE FUNCTION mean_1D (x,weight,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x,1)), INTENT(IN), OPTIONAL :: weight
    LOGICAL, DIMENSION(SIZE(x,1)), INTENT(IN), OPTIONAL :: mask
    REAL(DP) :: mean_1D
    INTEGER, DIMENSION(4) :: N 
    N = [SIZE(x,1),1,1,1]
    IF (.NOT. PRESENT(weight)) THEN
      IF (.NOT. PRESENT(mask)) THEN
        mean_1D = MEAN_4D(RESHAPE(x(:),N(:)))
      ELSE
        mean_1D = MEAN_4D(RESHAPE(x(:),N(:)),MASK=RESHAPE(mask(:),N(:)))
      END IF
    ELSE
      IF (.NOT. PRESENT(mask)) THEN
        mean_1D = MEAN_4D(RESHAPE(x(:),N(:)),WEIGHT=RESHAPE(weight(:),N(:)))
      ELSE
        mean_1D = MEAN_4D(RESHAPE(x(:),N(:)),WEIGHT=RESHAPE(weight(:),N(:)), &
                          MASK=RESHAPE(mask(:),N(:)))
      END IF
    END IF
  END FUNCTION mean_1D

  !-------------------------------------------------------------------------

  PURE FUNCTION mean_4D_dim (x,weight,mask,dim)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)), INTENT(IN), &
      OPTIONAL :: weight
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)), INTENT(IN), &
      OPTIONAL :: mask
    INTEGER, INTENT(IN) :: dim
    REAL(DP), DIMENSION(SIZE(x,MERGE(2,1,dim == 1)), &
                        SIZE(x,MERGE(3,2,dim == 1 .OR. dim == 2)), &
                        SIZE(x,MERGE(3,4,dim == 4))) :: mean_4D_dim

    INTEGER :: i1, i2, i3, i4, N1, N2, N3, N4

    !-----------------------------------------------------------------------

    N1 = SIZE(x(:,:,:,:),1)
    N2 = SIZE(x(:,:,:,:),2)
    N3 = SIZE(x(:,:,:,:),3)
    N4 = SIZE(x(:,:,:,:),4)
    SELECT CASE (dim)

      CASE (1)
        IF (.NOT. PRESENT(weight)) THEN 
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i2=1:N2,i3=1:N3,i4=1:N4) &
              mean_4D_dim(i2,i3,i4) = MEAN_1D(x(:,i2,i3,i4))
          ELSE
            FORALL (i2=1:N2,i3=1:N3,i4=1:N4) &
              mean_4D_dim(i2,i3,i4) = MEAN_1D(x(:,i2,i3,i4), &
                                              MASK=mask(:,i2,i3,i4))
          END IF
        ELSE
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i2=1:N2,i3=1:N3,i4=1:N4) &
              mean_4D_dim(i2,i3,i4) = MEAN_1D(x(:,i2,i3,i4), &
                                              WEIGHT=weight(:,i2,i3,i4))
          ELSE
            FORALL (i2=1:N2,i3=1:N3,i4=1:N4) &
              mean_4D_dim(i2,i3,i4) = MEAN_1D(x(:,i2,i3,i4), &
                                              MASK=mask(:,i2,i3,i4), &
                                              WEIGHT=weight(:,i2,i3,i4))
          END IF
        END IF

      CASE (2)
        IF (.NOT. PRESENT(weight)) THEN 
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1,i3=1:N3,i4=1:N4) &
              mean_4D_dim(i1,i3,i4) = MEAN_1D(x(i1,:,i3,i4))
          ELSE
            FORALL (i1=1:N1,i3=1:N3,i4=1:N4) &
              mean_4D_dim(i1,i3,i4) = MEAN_1D(x(i1,:,i3,i4), &
                                              MASK=mask(i1,:,i3,i4))
          END IF
        ELSE
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1,i3=1:N3,i4=1:N4) &
              mean_4D_dim(i1,i3,i4) = MEAN_1D(x(i1,:,i3,i4), &
                                              WEIGHT=weight(i1,:,i3,i4))
          ELSE
            FORALL (i1=1:N1,i3=1:N3,i4=1:N4) &
              mean_4D_dim(i1,i3,i4) = MEAN_1D(x(i1,:,i3,i4), &
                                              MASK=mask(i1,:,i3,i4), &
                                              WEIGHT=weight(i1,:,i3,i4))
          END IF
        END IF

      CASE (3)
        IF (.NOT. PRESENT(weight)) THEN 
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1,i2=1:N2,i4=1:N4) &
              mean_4D_dim(i1,i2,i4) = MEAN_1D(x(i1,i2,:,i4))
          ELSE
            FORALL (i1=1:N1,i2=1:N2,i4=1:N4) &
              mean_4D_dim(i1,i2,i4) = MEAN_1D(x(i1,i2,:,i4), &
                                              MASK=mask(i1,i2,:,i4))
          END IF
        ELSE
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1,i2=1:N2,i4=1:N4) &
              mean_4D_dim(i1,i2,i4) = MEAN_1D(x(i1,i2,:,i4), &
                                              WEIGHT=weight(i1,i2,:,i4))
          ELSE
            FORALL (i1=1:N1,i2=1:N2,i4=1:N4) &
              mean_4D_dim(i1,i2,i4) = MEAN_1D(x(i1,i2,:,i4), &
                                              MASK=mask(i1,i2,:,i4), &
                                              WEIGHT=weight(i1,i2,:,i4))
          END IF
        END IF

      CASE (4)
        IF (.NOT. PRESENT(weight)) THEN 
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1,i2=1:N2,i3=1:N3) &
              mean_4D_dim(i1,i2,i3) = MEAN_1D(x(i1,i2,i3,:))
          ELSE
            FORALL (i1=1:N1,i2=1:N2,i3=1:N3) &
              mean_4D_dim(i1,i2,i3) = MEAN_1D(x(i1,i2,i3,:), &
                                              MASK=mask(i1,i2,i3,:))
          END IF
        ELSE
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1,i2=1:N2,i3=1:N3) &
              mean_4D_dim(i1,i2,i3) = MEAN_1D(x(i1,i2,i3,:), &
                                              WEIGHT=weight(i1,i2,i3,:))
          ELSE
            FORALL (i1=1:N1,i2=1:N2,i3=1:N3) &
              mean_4D_dim(i1,i2,i3) = MEAN_1D(x(i1,i2,i3,:), &
                                              MASK=mask(i1,i2,i3,:), &
                                              WEIGHT=weight(i1,i2,i3,:))
          END IF
        END IF

    END SELECT

    !-----------------------------------------------------------------------

  END FUNCTION mean_4D_dim

  !-------------------------------------------------------------------------

  PURE FUNCTION mean_3D_dim (x,weight,mask,dim)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)), INTENT(IN), &
      OPTIONAL :: weight
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)), INTENT(IN), &
      OPTIONAL :: mask
    INTEGER, INTENT(IN) :: dim
    REAL(DP), DIMENSION(SIZE(x,MERGE(2,1,dim == 1)), &
                        SIZE(x,MERGE(2,3,dim == 3))) :: mean_3D_dim

    INTEGER :: i1, i2, i3, N1, N2, N3

    !-----------------------------------------------------------------------

    N1 = SIZE(x(:,:,:),1)
    N2 = SIZE(x(:,:,:),2)
    N3 = SIZE(x(:,:,:),3)
    SELECT CASE (dim)

      CASE (1)
        IF (.NOT. PRESENT(weight)) THEN 
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i2=1:N2,i3=1:N3) &
              mean_3D_dim(i2,i3) = MEAN_1D(x(:,i2,i3))
          ELSE
            FORALL (i2=1:N2,i3=1:N3) &
              mean_3D_dim(i2,i3) = MEAN_1D(x(:,i2,i3),MASK=mask(:,i2,i3))
          END IF
        ELSE
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i2=1:N2,i3=1:N3) &
              mean_3D_dim(i2,i3) = MEAN_1D(x(:,i2,i3),WEIGHT=weight(:,i2,i3))
          ELSE
            FORALL (i2=1:N2,i3=1:N3) &
              mean_3D_dim(i2,i3) = MEAN_1D(x(:,i2,i3),MASK=mask(:,i2,i3), &
                                           WEIGHT=weight(:,i2,i3))
          END IF
        END IF

      CASE (2)
        IF (.NOT. PRESENT(weight)) THEN 
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1,i3=1:N3) &
              mean_3D_dim(i1,i3) = MEAN_1D(x(i1,:,i3))
          ELSE
            FORALL (i1=1:N1,i3=1:N3) &
              mean_3D_dim(i1,i3) = MEAN_1D(x(i1,:,i3),MASK=mask(i1,:,i3))
          END IF
        ELSE
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1,i3=1:N3) &
              mean_3D_dim(i1,i3) = MEAN_1D(x(i1,:,i3),WEIGHT=weight(i1,:,i3))
          ELSE
            FORALL (i1=1:N1,i3=1:N3) &
              mean_3D_dim(i1,i3) = MEAN_1D(x(i1,:,i3),MASK=mask(i1,:,i3), &
                                           WEIGHT=weight(i1,:,i3))
          END IF
        END IF

      CASE (3)
        IF (.NOT. PRESENT(weight)) THEN 
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1,i2=1:N2) &
              mean_3D_dim(i1,i2) = MEAN_1D(x(i1,i2,:))
          ELSE
            FORALL (i1=1:N1,i2=1:N2) &
              mean_3D_dim(i1,i2) = MEAN_1D(x(i1,i2,:),MASK=mask(i1,i2,:))
          END IF
        ELSE
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1,i2=1:N2) &
              mean_3D_dim(i1,i2) = MEAN_1D(x(i1,i2,:),WEIGHT=weight(i1,i2,:))
          ELSE
            FORALL (i1=1:N1,i2=1:N2) &
              mean_3D_dim(i1,i2) = MEAN_1D(x(i1,i2,:),MASK=mask(i1,i2,:), &
                                           WEIGHT=weight(i1,i2,:))
          END IF
        END IF

    END SELECT

    !-----------------------------------------------------------------------

  END FUNCTION mean_3D_dim

  !-------------------------------------------------------------------------

  PURE FUNCTION mean_2D_dim (x,weight,mask,dim)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x,1),SIZE(x,2)), INTENT(IN), OPTIONAL :: weight
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2)), INTENT(IN), OPTIONAL :: mask
    INTEGER, INTENT(IN) :: dim
    REAL(DP), DIMENSION(SIZE(x,MERGE(2,1,dim == 1))) :: mean_2D_dim

    INTEGER :: i1, i2, N1, N2

    !-----------------------------------------------------------------------

    N1 = SIZE(x(:,:),1)
    N2 = SIZE(x(:,:),2)
    SELECT CASE (dim)

      CASE (1)
        IF (.NOT. PRESENT(weight)) THEN 
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i2=1:N2) mean_2D_dim(i2) = MEAN_1D(x(:,i2))
          ELSE
            FORALL (i2=1:N2) &
              mean_2D_dim(i2) = MEAN_1D(x(:,i2),MASK=mask(:,i2))
          END IF
        ELSE
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i2=1:N2) &
              mean_2D_dim(i2) = MEAN_1D(x(:,i2),WEIGHT=weight(:,i2))
          ELSE
            FORALL (i2=1:N2) &
              mean_2D_dim(i2) = MEAN_1D(x(:,i2),MASK=mask(:,i2), &
                                        WEIGHT=weight(:,i2))
          END IF
        END IF

      CASE (2)
        IF (.NOT. PRESENT(weight)) THEN 
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1) mean_2D_dim(i1) = MEAN_1D(x(i1,:))
          ELSE
            FORALL (i1=1:N1) &
              mean_2D_dim(i1) = MEAN_1D(x(i1,:),MASK=mask(i1,:))
          END IF
        ELSE
          IF (.NOT. PRESENT(mask)) THEN
            FORALL (i1=1:N1) &
              mean_2D_dim(i1) = MEAN_1D(x(i1,:),WEIGHT=weight(i1,:))
          ELSE
            FORALL (i1=1:N1) &
              mean_2D_dim(i1) = MEAN_1D(x(i1,:),MASK=mask(i1,:), &
                                        WEIGHT=weight(i1,:))
          END IF
        END IF

    END SELECT

    !-----------------------------------------------------------------------

  END FUNCTION mean_2D_dim


  !==========================================================================
  ! s(X) = SIGMA(X,dim,mask)
  !        SIGMA(X[N],weight[N])
  !
  !   Computes the standard deviation of the variable X, with weights w. If w 
  ! is not supplied then the distribution is assumed to be uniform.
  !==========================================================================

  PURE FUNCTION sigma_4D (x,weight,mask)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)), INTENT(IN), &
      OPTIONAL :: weight
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)), INTENT(IN), &
      OPTIONAL :: mask
    REAL(DP) :: sigma_4D

    INTEGER :: N
    REAL(DP) :: norm, av
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)) :: mask4D
    
    !-----------------------------------------------------------------------

    IF (PRESENT(mask)) THEN ; mask4D(:,:,:,:) = mask(:,:,:,:)
                       ELSE ; mask4D(:,:,:,:) = .True. ; END IF
    N = COUNT(mask4D(:,:,:,:))
    IF (.NOT. PRESENT(weight)) THEN 
      av = SUM(x(:,:,:,:),MASK=mask4D(:,:,:,:)) / N
      sigma_4D = SQRT( SUM((x(:,:,:,:)-av)**2,MASK=mask4D(:,:,:,:))/(N-1._DP) )
    ELSE
      norm = SUM(weight(:,:,:,:),MASK=mask4D(:,:,:,:))
      av = SUM(x(:,:,:,:)*weight(:,:,:,:),MASK=mask4D(:,:,:,:)) / norm
      sigma_4D = SQRT( SUM((x(:,:,:,:)-av)**2*weight(:,:,:,:), &
                           MASK=mask4D(:,:,:,:)) &
                     / norm * N/(N-1._DP) )
    END IF

    !-----------------------------------------------------------------------

  END FUNCTION sigma_4D

  ! Interface shells
  !-----------------
  ! 3D
  PURE FUNCTION sigma_3D (x,weight,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)), INTENT(IN), &
      OPTIONAL :: weight
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)), INTENT(IN), &
      OPTIONAL :: mask
    REAL(DP) :: sigma_3D
    INTEGER, DIMENSION(4) :: N
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),1) :: mask4D
    N = [SIZE(x,1),SIZE(x,2),SIZE(x,3),1]
    IF (PRESENT(mask)) THEN ; mask4D(:,:,:,1) = mask(:,:,:)
                       ELSE ; mask4D(:,:,:,:) = .True. ; END IF
    IF (.NOT. PRESENT(weight)) THEN
      sigma_3D = SIGMA_4D(RESHAPE(x(:,:,:),N(:)),MASK=mask4D(:,:,:,:))
    ELSE
      sigma_3D = SIGMA_4D(RESHAPE(x(:,:,:),N(:)), &
                          WEIGHT=RESHAPE(weight(:,:,:),N(:)), &
                          MASK=mask4D(:,:,:,:))
    END IF
  END FUNCTION sigma_3D

  ! 2D
  PURE FUNCTION sigma_2D (x,weight,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x,1),SIZE(x,2)), INTENT(IN), OPTIONAL :: weight
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2)), INTENT(IN), OPTIONAL :: mask
    REAL(DP) :: sigma_2D
    INTEGER, DIMENSION(4) :: N
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),1,1) :: mask4D
    N = [SIZE(x,1),SIZE(x,2),1,1]
    IF (PRESENT(mask)) THEN ; mask4D(:,:,1,1) = mask(:,:)
                       ELSE ; mask4D(:,:,:,:) = .True. ; END IF
    IF (.NOT. PRESENT(weight)) THEN
      sigma_2D = SIGMA_4D(RESHAPE(x(:,:),N(:)),MASK=mask4D(:,:,:,:))
    ELSE
      sigma_2D = SIGMA_4D(RESHAPE(x(:,:),N(:)), &
                          WEIGHT=RESHAPE(weight(:,:),N(:)),MASK=mask4D(:,:,:,:))
    END IF
  END FUNCTION sigma_2D

  ! 1D
  PURE FUNCTION sigma_1D (x,weight,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(SIZE(x,1)), INTENT(IN), OPTIONAL :: weight
    LOGICAL, DIMENSION(SIZE(x,1)), INTENT(IN), OPTIONAL :: mask
    REAL(DP) :: sigma_1D
    INTEGER, DIMENSION(4) :: N
    LOGICAL, DIMENSION(SIZE(x,1),1,1,1) :: mask4D
    N = [SIZE(x,1),1,1,1]
    IF (PRESENT(mask)) THEN ; mask4D(:,1,1,1) = mask(:)
                       ELSE ; mask4D(:,:,:,:) = .True. ; END IF
    IF (.NOT. PRESENT(weight)) THEN
      sigma_1D = SIGMA_4D(RESHAPE(x(:),N(:)),MASK=mask4D(:,:,:,:))
    ELSE
      sigma_1D = SIGMA_4D(RESHAPE(x(:),N(:)),WEIGHT=RESHAPE(weight(:),N(:)), &
                          MASK=mask4D(:,:,:,:))
    END IF
  END FUNCTION sigma_1D

  !-------------------------------------------------------------------------

  ! 4D dim
  PURE FUNCTION sigma_4D_dim (x,dim,mask)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: dim
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)), INTENT(IN), &
      OPTIONAL:: mask
    REAL(DP), DIMENSION(SIZE(x,MERGE(2,1,dim == 1)), &
                        SIZE(x,MERGE(3,2,dim == 1 .OR. dim == 2)), &
                        SIZE(x,MERGE(3,4,dim == 4))) :: sigma_4D_dim

    INTEGER :: i1, i2, i3, i4, N1, N2, N3, N4
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)) :: mask4D

    !-----------------------------------------------------------------------

    N1 = SIZE(x(:,:,:,:),1)
    N2 = SIZE(x(:,:,:,:),2)
    N3 = SIZE(x(:,:,:,:),3)
    N4 = SIZE(x(:,:,:,:),4)
    IF (PRESENT(mask)) THEN ; mask4D(:,:,:,:) = mask(:,:,:,:)
                       ELSE ; mask4D(:,:,:,:) = .True. ; END IF
    SELECT CASE (dim)
      CASE (1)
        FORALL (i2=1:N2,i3=1:N3,i4=1:N4) &
          sigma_4D_dim(i2,i3,i4) = SIGMA_1D(x(:,i2,i3,i4), &
                                            MASK=mask4D(:,i2,i3,i4))
      CASE (2)
        FORALL (i1=1:N1,i3=1:N3,i4=1:N4) &
          sigma_4D_dim(i1,i3,i4) = SIGMA_1D(x(i1,:,i3,i4), &
                                            MASK=mask4D(i1,:,i3,i4))
      CASE (3)
        FORALL (i1=1:N1,i2=1:N2,i4=1:N4) &
          sigma_4D_dim(i1,i2,i4) = SIGMA_1D(x(i1,i2,:,i4), &
                                            MASK=mask4D(i1,i2,:,i4))
      CASE (4)
        FORALL (i1=1:N1,i2=1:N2,i3=1:N3) &
          sigma_4D_dim(i1,i2,i3) = SIGMA_1D(x(i1,i2,i3,:), &
                                            MASK=mask4D(i1,i2,i3,:))
    END SELECT

    !-----------------------------------------------------------------------

  END FUNCTION sigma_4D_dim

  !-------------------------------------------------------------------------

  ! 3D dim
  PURE FUNCTION sigma_3D_dim (x,dim,mask)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: dim
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)), INTENT(IN), &
      OPTIONAL :: mask
    REAL(DP), DIMENSION(SIZE(x,MERGE(2,1,dim == 1)), &
                        SIZE(x,MERGE(2,3,dim == 3))) :: sigma_3D_dim

    INTEGER :: i1, i2, i3, N1, N2, N3
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)) :: mask3D

    !-----------------------------------------------------------------------

    N1 = SIZE(x(:,:,:),1)
    N2 = SIZE(x(:,:,:),2)
    N3 = SIZE(x(:,:,:),3)
    IF (PRESENT(mask)) THEN ; mask3D(:,:,:) = mask(:,:,:)
                       ELSE ; mask3D(:,:,:) = .True. ; END IF
    SELECT CASE (dim)
      CASE (1)
        FORALL (i2=1:N2,i3=1:N3) &
          sigma_3D_dim(i2,i3) = SIGMA_1D(x(:,i2,i3),MASK=mask3D(:,i2,i3))
      CASE (2)
        FORALL (i1=1:N1,i3=1:N3) &
          sigma_3D_dim(i1,i3) = SIGMA_1D(x(i1,:,i3),MASK=mask3D(i1,:,i3))
      CASE (3)
        FORALL (i1=1:N1,i2=1:N2) &
          sigma_3D_dim(i1,i2) = SIGMA_1D(x(i1,i2,:),MASK=mask3D(i1,i2,:))
    END SELECT

    !-----------------------------------------------------------------------

  END FUNCTION sigma_3D_dim

  !-------------------------------------------------------------------------

  ! 2D dim
  PURE FUNCTION sigma_2D_dim (x,dim,mask)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: dim
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2)), INTENT(IN), OPTIONAL :: mask
    REAL(DP), DIMENSION(SIZE(x,MERGE(2,1,dim == 1))) :: sigma_2D_dim

    INTEGER :: i1, i2, N1, N2
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2)) :: mask2D

    !-----------------------------------------------------------------------

    N1 = SIZE(x(:,:),1)
    N2 = SIZE(x(:,:),2)
    IF (PRESENT(mask)) THEN ; mask2D(:,:) = mask(:,:)
                       ELSE ; mask2D(:,:) = .True. ; END IF
    SELECT CASE (dim)
      CASE (1)
        FORALL (i2=1:N2) &
          sigma_2D_dim(i2) = SIGMA_1D(x(:,i2),MASK=mask2D(:,i2))
      CASE (2)
        FORALL (i1=1:N1) &
          sigma_2D_dim(i1) = SIGMA_1D(x(i1,:),MASK=mask2D(i1,:))
    END SELECT

    !-----------------------------------------------------------------------

  END FUNCTION sigma_2D_dim


  !==========================================================================
  ! m(X) = MOMENT_DATA(X[N],order,weight[N])
  !
  !   Computes the orderieth moment of the variable X, with weights w. If
  ! w is not supplied then the distribution is assumed to be uniform.
  !==========================================================================

  FUNCTION moment_data (x,order,weight)

    USE utilities, ONLY: DP, strike, flEQ, warning
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    INTEGER, INTENT(IN) :: order
    REAL(DP), INTENT(IN), DIMENSION(:), OPTIONAL :: weight
    REAL(DP) :: moment_data

    INTEGER :: N
    REAL(DP) :: norm, mean
    REAL(DP), DIMENSION(SIZE(x)) :: w

    !-----------------------------------------------------------------------

    ! Check input
    N = SIZE(x)
    IF (PRESENT(weight)) THEN
      w = weight 
    ELSE 
      w = SPREAD(1._DP/N,1,N)
    END IF
    IF (SIZE(w) /= N) CALL STRIKE ("MOMENT_DATA","bad input.")

    ! Check normalisation if order>0, otherwise the function returns the
    ! normalisation factor.
    normalisation: IF (order > 0) THEN
      norm = SUM(w)
      IF (.NOT. flEQ(norm,1._DP,TOL=1.E-2_DP)) &
        CALL WARNING ("MOMENT_DATA","the distribution was not normalised.")
      w = w/norm
    ELSE IF (order < 0) THEN
      CALL STRIKE ("MOMENT_DATA","order must be positive or zero.")
    END IF normalisation

    ! Compute the moments
    moment: SELECT CASE (order-1)
      CASE (:-1)
        moment_data = SUM(w)
      CASE (0)
        moment_data = SUM(x*w)
      CASE (+1:)
        mean = SUM(x*w)
        moment_data = SUM((x-mean)**order*w)
    END SELECT moment

    !-----------------------------------------------------------------------

  END FUNCTION moment_data


  !==========================================================================
  ! m = MEDIAN_DATA(X[1D to 4D],dim=,mask)
  !     MEDIAN_DATA(X[N],weight[N],mask)
  !
  !   Computes the median of a variable X, using part of the Hoare sorting
  ! method. If the weight W is not supplied then the distribution is assumed
  ! to be uniform.
  !==========================================================================

  FUNCTION median_data_4D (x,mask)

    USE utilities, ONLY: DP
    USE arrays, ONLY: sort
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: x
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)), INTENT(IN), &
      OPTIONAL :: mask
    REAL(DP) :: median_data_4D

    INTEGER :: N
    REAL(DP), DIMENSION(SIZE(x)) :: xs
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)) :: mask4D
    
    !-----------------------------------------------------------------------

    ! Mask
    IF (PRESENT(mask)) THEN ; mask4D(:,:,:,:) = mask(:,:,:,:)
                       ELSE ; mask4D(:,:,:,:) = .True. ; END IF
      
    ! Sort the arrays
    N = COUNT(mask4D(:,:,:,:))
    xs(:) = SORT(PACK(x(:,:,:,:),MASK=mask4D(:,:,:,:)))

    ! Even / odd cases
    median_data_4D = MERGE( (xs(N/2) + xs(N/2+1)) / 2._DP, &
                            xs(N/2+1), ( MOD(N,2) == 0 ) )

    !-----------------------------------------------------------------------

  END FUNCTION median_data_4D

  ! Interface
  FUNCTION median_data_3D (x,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: x
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)), INTENT(IN), &
      OPTIONAL :: mask
    REAL(DP) :: median_data_3D
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),1) :: mask4D
    IF (PRESENT(mask)) THEN ; mask4D(:,:,:,1) = mask(:,:,:)
                       ELSE ; mask4D(:,:,:,:) = .True. ; END IF
    median_data_3D = MEDIAN_DATA_4D(RESHAPE(x(:,:,:), &
                                            [SIZE(x,1),SIZE(x,2),SIZE(x,3),1]),&
                                    MASK=mask4D(:,:,:,:))
  END FUNCTION median_data_3D
  FUNCTION median_data_2D (x,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: x
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2)), INTENT(IN), OPTIONAL :: mask
    REAL(DP) :: median_data_2D
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),1,1) :: mask4D
    IF (PRESENT(mask)) THEN ; mask4D(:,:,1,1) = mask(:,:)
                       ELSE ; mask4D(:,:,:,:) = .True. ; END IF
    median_data_2D = MEDIAN_DATA_4D(RESHAPE(x(:,:),[SIZE(x,1),SIZE(x,2),1,1]), &
                                    MASK=mask4D(:,:,:,:))
  END FUNCTION median_data_2D
  FUNCTION median_data_1D (x,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    LOGICAL, DIMENSION(SIZE(x)), INTENT(IN), OPTIONAL :: mask
    REAL(DP) :: median_data_1D
    LOGICAL, DIMENSION(SIZE(x),1,1,1) :: mask4D
    IF (PRESENT(mask)) THEN ; mask4D(:,1,1,1) = mask(:)
                       ELSE ; mask4D(:,:,:,:) = .True. ; END IF
    median_data_1D = MEDIAN_DATA_4D(RESHAPE(x(:),[SIZE(x,1),1,1,1]), &
                                    MASK=mask4D(:,:,:,:))
  END FUNCTION median_data_1D

  !-------------------------------------------------------------------------

  FUNCTION median_data_4D_dim (x,dim,mask)

    USE utilities, ONLY: DP, strike, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: dim
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)), INTENT(IN), &
      OPTIONAL :: mask
    REAL(DP), DIMENSION(SIZE(x,MERGE(2,1,dim == 1)), &
                        SIZE(x,MERGE(3,2,dim == 1 .OR. dim == 2)), &
                        SIZE(x,MERGE(3,4,dim == 4))) :: median_data_4D_dim

    INTEGER :: i1, i2, i3, i4
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)) :: mask4D

    !-----------------------------------------------------------------------

    IF (PRESENT(mask)) THEN ; mask4D(:,:,:,:) = mask(:,:,:,:)
                       ELSE ; mask4D(:,:,:,:) = .True. ; END IF
    SELECT CASE (dim)
      CASE (1)
        DO i2=1,SIZE(x(:,:,:,:),2)
          DO i3=1,SIZE(x(:,:,:,:),3)
            DO i4=1,SIZE(x(:,:,:,:),4)
              median_data_4D_dim(i2,i3,i4) = MEDIAN_DATA_1D(x(:,i2,i3,i4), &
                                                        MASK=mask4D(:,i2,i3,i4))
            END DO
          END DO
        END DO
      CASE (2)
        DO i1=1,SIZE(x(:,:,:,:),1)
          DO i3=1,SIZE(x(:,:,:,:),3)
            DO i4=1,SIZE(x(:,:,:,:),4)
              median_data_4D_dim(i1,i3,i4) = MEDIAN_DATA_1D(x(i1,:,i3,i4), &
                                                        MASK=mask4D(i1,:,i3,i4))
            END DO
          END DO
        END DO
      CASE (3)
        DO i1=1,SIZE(x(:,:,:,:),1)
          DO i2=1,SIZE(x(:,:,:,:),2)
            DO i4=1,SIZE(x(:,:,:,:),4)
              median_data_4D_dim(i1,i2,i4) = MEDIAN_DATA_1D(x(i1,i2,:,i4), &
                                                        MASK=mask4D(i1,i2,:,i4))
            END DO
          END DO
        END DO
      CASE (4)
        DO i1=1,SIZE(x(:,:,:,:),1)
          DO i2=1,SIZE(x(:,:,:,:),2)
            DO i3=1,SIZE(x(:,:,:,:),3)
              median_data_4D_dim(i1,i2,i3) = MEDIAN_DATA_1D(x(i1,i2,i3,:), &
                                                        MASK=mask4D(i1,i2,i3,:))
            END DO
          END DO
        END DO
      CASE DEFAULT
        median_data_4D_dim(:,:,:) = NaN()
        CALL STRIKE("MEDIAN_DATA","Wrong value for keyword DIM")
    END SELECT

    !-----------------------------------------------------------------------

  END FUNCTION median_data_4D_dim

  !-------------------------------------------------------------------------

  FUNCTION median_data_3D_dim (x,dim,mask)

    USE utilities, ONLY: DP, strike, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: dim
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)), INTENT(IN), &
      OPTIONAL :: mask
    REAL(DP), DIMENSION(SIZE(x,MERGE(2,1,dim == 1)), &
                        SIZE(x,MERGE(2,3,dim == 3))) :: median_data_3D_dim

    INTEGER :: i1, i2, i3
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)) :: mask3D

    !-----------------------------------------------------------------------

    IF (PRESENT(mask)) THEN ; mask3D(:,:,:) = mask(:,:,:)
                       ELSE ; mask3D(:,:,:) = .True. ; END IF
    SELECT CASE (dim)
      CASE (1)
        DO i2=1,SIZE(x(:,:,:),2)
          DO i3=1,SIZE(x(:,:,:),3)
            median_data_3D_dim(i2,i3) = MEDIAN_DATA_1D(x(:,i2,i3), &
                                                       MASK=mask3D(:,i2,i3))
          END DO
        END DO
      CASE (2)
        DO i1=1,SIZE(x(:,:,:),1)
          DO i3=1,SIZE(x(:,:,:),3)
            median_data_3D_dim(i1,i3) = MEDIAN_DATA_1D(x(i1,:,i3), &
                                                       MASK=mask3D(i1,:,i3))
          END DO
        END DO
      CASE (3)
        DO i1=1,SIZE(x(:,:,:),1)
          DO i2=1,SIZE(x(:,:,:),2)
            median_data_3D_dim(i1,i2) = MEDIAN_DATA_1D(x(i1,i2,:), &
                                                       MASK=mask3D(i1,i2,:))
          END DO
        END DO
      CASE DEFAULT
        median_data_3D_dim(:,:) = NaN()
        CALL STRIKE("MEDIAN_DATA","Wrong value for keyword DIM")
    END SELECT

    !-----------------------------------------------------------------------

  END FUNCTION median_data_3D_dim

  !-------------------------------------------------------------------------

  FUNCTION median_data_2D_dim (x,dim,mask)

    USE utilities, ONLY: DP, strike, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: dim
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2)), INTENT(IN), OPTIONAL :: mask
    REAL(DP), DIMENSION(SIZE(x,MERGE(2,1,dim == 1))) :: median_data_2D_dim

    INTEGER :: i1, i2
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2)) :: mask2D

    !-----------------------------------------------------------------------

    IF (PRESENT(mask)) THEN ; mask2D(:,:) = mask(:,:)
                       ELSE ; mask2D(:,:) = .True. ; END IF
    SELECT CASE (dim)
      CASE (1)
        DO i2=1,SIZE(x(:,:),2)
          median_data_2D_dim(i2) = MEDIAN_DATA_1D(x(:,i2),MASK=mask2D(:,i2))
        END DO
      CASE (2)
        DO i1=1,SIZE(x(:,:),1)
          median_data_2D_dim(i1) = MEDIAN_DATA_1D(x(i1,:),MASK=mask2D(i1,:))
        END DO
      CASE DEFAULT
        median_data_2D_dim(:) = NaN()
        CALL STRIKE("MEDIAN_DATA","Wrong value for keyword DIM")
    END SELECT

    !-----------------------------------------------------------------------

  END FUNCTION median_data_2D_dim

  !-------------------------------------------------------------------------

  FUNCTION median_data_wei (x,weight)

    USE utilities, ONLY: DP, strike, flEQ, warning
    USE arrays, ONLY: sort, closest
    USe interpolation, ONLY: interp_lin_sorted
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), INTENT(IN), DIMENSION(:) :: weight
    REAL(DP) :: median_data_wei

    INTEGER :: N, i, Nun
    REAL(DP) :: norm
    INTEGER, DIMENSION(SIZE(x)) :: ind
    REAL(DP), DIMENSION(SIZE(x)) :: xs, ws, F
    LOGICAL, DIMENSION(SIZE(x)) :: uniq

    !-----------------------------------------------------------------------

    ! Check input
    N = SIZE(x)
    ws = weight 
    IF (SIZE(ws) /= N) CALL STRIKE ("MEDIAN_DATA","bad input.")

    ! Check normalisation
    norm = SUM(ws)
    IF (.NOT. flEQ(norm,1._DP,TOL=1.E-2_DP)) &
      CALL WARNING ("MEDIAN_DATA","the distribution was not normalised.")
    ws = ws / norm

    ! Sort the arrays
    xs = SORT(x,IND=ind)
    ws = ws(ind)

    ! Repartition function (cumulate the duplicate values in one bin)
    F(1) = 0._DP
    uniq = .True.
    DO i=2,N
      F(i) = F(i-1) + ws(i)
      uniq(i-1) = MERGE(.True.,.False.,xs(i) > xs(i-1))
    END DO

    ! Interpolation to get a good estimate of the median
    Nun = COUNT(uniq)
    IF (Nun == 1) THEN
      median_data_wei = x(1)
    ELSE
      median_data_wei = INTERP_LIN_SORTED(xs(:),F(:),0.5_DP)
    END IF

    !-----------------------------------------------------------------------

  END FUNCTION median_data_wei


  !==========================================================================
  ! CALL MEDIAN_CONF(X[1-4D],CONF=0.99,x0,dxinf,dxsup)
  !
  !   Computes the median of a variable X, and the lower and upper limits
  ! corresponding to a confidence inter val of CONF.
  !==========================================================================

  SUBROUTINE median_conf_4D (x,conf,x0,dxinf,dxsup)

    USE utilities, ONLY: DP
    USE arrays, ONLY: sort
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:,:,:,:) :: x
    REAL(DP), INTENT(IN) :: conf
    REAL(DP), INTENT(OUT) :: x0, dxinf, dxsup

    INTEGER :: N, N0, Ninf, Nsup
    INTEGER, DIMENSION(SIZE(x)) :: ind
    REAL(DP), DIMENSION(SIZE(x)) :: xs

    !-----------------------------------------------------------------------

    ! Sort the arrays
    N = SIZE(x(:,:,:,:))
    xs(:) = SORT(RESHAPE(x(:,:,:,:),[N]),IND=ind)

    ! Even / odd cases
    N0 = N / 2
    x0 = MERGE( (xs(N0) + xs(N0+1)) / 2._DP, xs(N0+1), ( MOD(N,2) == 0 ) )
    Ninf = NINT( (1._DP-conf)/2._DP * N )
    Nsup = NINT( (1._DP+conf)/2._DP * N )
    dxinf = x0 - xs(Ninf)
    dxsup = xs(Nsup) - x0

    !-----------------------------------------------------------------------

  END SUBROUTINE median_conf_4D

  ! 3D
  SUBROUTINE median_conf_3D (x,conf,x0,dxinf,dxsup)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(:,:,:) :: x
    REAL(DP), INTENT(IN) :: conf
    REAL(DP), INTENT(OUT) :: x0, dxinf, dxsup
    INTEGER :: N1, N2, N3
    N1 = SIZE(x(:,:,:),1)
    N2 = SIZE(x(:,:,:),2)
    N3 = SIZE(x(:,:,:),3)
    CALL MEDIAN_CONF_4D(RESHAPE(x(:,:,:),[N1,N2,N3,1]),conf,x0,dxinf,dxsup)
  END SUBROUTINE median_conf_3D

  ! 2D
  SUBROUTINE median_conf_2D (x,conf,x0,dxinf,dxsup)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: x
    REAL(DP), INTENT(IN) :: conf
    REAL(DP), INTENT(OUT) :: x0, dxinf, dxsup
    INTEGER :: N1, N2
    N1 = SIZE(x(:,:),1)
    N2 = SIZE(x(:,:),2)
    CALL MEDIAN_CONF_4D(RESHAPE(x(:,:),[N1,N2,1,1]),conf,x0,dxinf,dxsup)
  END SUBROUTINE median_conf_2D

  ! 1D
  SUBROUTINE median_conf_1D (x,conf,x0,dxinf,dxsup)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), INTENT(IN) :: conf
    REAL(DP), INTENT(OUT) :: x0, dxinf, dxsup
    INTEGER :: N1
    N1 = SIZE(x(:),1)
    CALL MEDIAN_CONF_4D(RESHAPE(x(:),[N1,1,1,1]),conf,x0,dxinf,dxsup)
  END SUBROUTINE median_conf_1D


  !==========================================================================
  ! CALL INFO_DATA (X[N],WEIGHT[N],name)
  !
  !   Computes the basic statistical quantities of a set of data.
  !==========================================================================

  SUBROUTINE info_data_1D (x,weight,name)

    USE utilities, ONLY: DP, trimLR
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), INTENT(IN), DIMENSION(:), OPTIONAL :: weight
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    INTEGER, PARAMETER :: l=20

    INTEGER :: N
    CHARACTER(l) :: title
    REAL(DP) :: xinf, xsup, xmean, sigx
    REAL(DP), DIMENSION(SIZE(x)) :: w

    !-----------------------------------------------------------------------

    ! Computes the quantities    
    N = SIZE(x)
    xinf = MINVAL(x) 
    xsup = MAXVAL(x)
    IF (PRESENT(weight)) THEN
      w = weight 
    ELSE 
      w = SPREAD(1._DP/N,1,N)
    END IF
    xmean = MOMENT_DATA(x,ORDER=1,WEIGHT=w)
    sigx = SQRT(MOMENT_DATA(x,ORDER=2,WEIGHT=w))

    ! Print the summary
    IF (PRESENT(name)) THEN 
      title = name 
    ELSE
      title = REPEAT(" ",l)
    END IF
    PRINT*, " => INFO_DATA: "//TRIMLR(title)//" = ",xmean,"+-",sigx
    PRINT*, "               "//TRIMLR(title)//" in [",xinf,",",xsup,"]"
    PRINT*, "               "//TRIMLR(title)//" has ",N," elements."

    !-----------------------------------------------------------------------

  END SUBROUTINE info_data_1D

  !==========================================================================

  SUBROUTINE info_data_2D (x,name)

    USE utilities, ONLY: DP, trimLR
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:,:) :: x
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    INTEGER, PARAMETER :: l=20

    INTEGER :: N, N1, N2
    CHARACTER(l) :: title
    REAL(DP) :: xinf, xsup, xmean, sigx

    !-----------------------------------------------------------------------

    ! Computes the quantities    
    N = SIZE(x)
    N1 = SIZE(x,1)
    N2 = SIZE(x,2)
    xinf = MINVAL(x) 
    xsup = MAXVAL(x)
    xmean = MOMENT_DATA(RESHAPE(x,SHAPE=(/N/)),ORDER=1)
    sigx = SQRT(MOMENT_DATA(RESHAPE(x,SHAPE=(/N/)),ORDER=2))

    ! Print the summary
    IF (PRESENT(name)) THEN 
      title = name 
    ELSE
      title = REPEAT(" ",l)
    END IF
    PRINT*, " => INFO_DATA: "//TRIMLR(title)//" = ",xmean,"+-",sigx
    PRINT*, "               "//TRIMLR(title)//" in [",xinf,",",xsup,"]"
    PRINT*, "               "//TRIMLR(title)//" has the profile [",N1,N2,"]"

    !-----------------------------------------------------------------------

  END SUBROUTINE info_data_2D


  !==========================================================================
  ! fwhm = FWHM_DATA(x[N],F[X])
  !
  !   Computes the full width at half maximum of any distribution. The
  ! distribution is supposed to be monomodal and well sampled.
  !==========================================================================

  FUNCTION fwhm_data (x,F)
 
    USE utilities, ONLY: DP
    IMPLICIT NONE
 
    REAL(DP), DIMENSION(:) :: x, F
    REAL(DP) :: fwhm_data

    INTEGER :: Nhm
    REAL(DP) :: hm
    REAL(DP), DIMENSION(:), ALLOCATABLE :: xhm

    !-----------------------------------------------------------------------

    hm = MAXVAL(F)/2._DP
    Nhm = COUNT(F>=hm)
    ALLOCATE(xhm(Nhm))
    xhm = PACK(x,F>=hm)

    fwhm_data = MAXVAL(xhm) - MINVAL(xhm)
    DEALLOCATE (xhm)

    !-----------------------------------------------------------------------

  END FUNCTION fwhm_data


  !==========================================================================
  ! corr[Ncorr] = RMAT2CORR(Rmat[Npar,Npar],Npar,Ncorr)
  !
  !   Extract the correlation coefficients from the correlation matrix.
  !==========================================================================

  PURE FUNCTION Rmat2corr_DP (Rmat,Npar,Ncorr)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(Npar,Npar), INTENT(IN) :: Rmat
    INTEGER, INTENT(IN) :: Npar, Ncorr
    REAL(DP), DIMENSION(Ncorr) :: Rmat2corr_DP

    INTEGER :: i, j, k

    !------------------------------------------------------------------------

    k = 1
    DO i=1,Npar
      DO j=1,Npar
        IF (MERGE(j > i,i > j,upper)) THEN
          Rmat2corr_DP(k) = Rmat(i,j)
          k = k + 1
        END IF
      END DO
    END DO

    !------------------------------------------------------------------------

  END FUNCTION Rmat2corr_DP

  !--------------------------------------------------------------------------

  PURE FUNCTION Rmat2corr_log (Rmat,Npar,Ncorr)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    LOGICAL, DIMENSION(Npar,Npar), INTENT(IN) :: Rmat
    INTEGER, INTENT(IN) :: Npar, Ncorr
    LOGICAL, DIMENSION(Ncorr) :: Rmat2corr_log

    INTEGER :: i, j, k

    !------------------------------------------------------------------------

    k = 1
    DO i=1,Npar
      DO j=1,Npar
        IF (MERGE(j > i,i > j,upper)) THEN
          Rmat2corr_log(k) = Rmat(i,j)
          k = k + 1
        END IF
      END DO
    END DO

    !------------------------------------------------------------------------

  END FUNCTION Rmat2corr_log

  !--------------------------------------------------------------------------

  PURE FUNCTION Rmat2corr_char (Rmat,Npar,Ncorr)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    CHARACTER(*), DIMENSION(Npar,Npar), INTENT(IN) :: Rmat
    INTEGER, INTENT(IN) :: Npar, Ncorr
    CHARACTER(LEN(Rmat)), DIMENSION(Ncorr) :: Rmat2corr_char

    INTEGER :: i, j, k

    !------------------------------------------------------------------------

    k = 1
    DO i=1,Npar
      DO j=1,Npar
        IF (MERGE(j > i,i > j,upper)) THEN
          Rmat2corr_char(k) = Rmat(i,j)
          k = k + 1
        END IF
      END DO
    END DO
    
    !------------------------------------------------------------------------

  END FUNCTION Rmat2corr_char


  !==========================================================================
  ! Rmat[Npar,Npar] = CORR2RMAT(corr[Ncorr],Npar)
  !
  !   Build a correlation matrix, from the correlation coefficients.
  !==========================================================================

  PURE FUNCTION corr2Rmat_DP (corr,Npar)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: corr
    INTEGER, INTENT(IN) :: Npar
    REAL(DP), DIMENSION(Npar,Npar) :: corr2Rmat_DP

    INTEGER :: i, j, k

    !------------------------------------------------------------------------
    
    ! Upper and lower triangles
    k = 1
    DO i=1,Npar
      DO j=1,Npar
        IF (j > i) THEN
          corr2Rmat_DP(i,j) = corr(k)
          corr2Rmat_DP(j,i) = corr(k)
          k = k + 1
        END IF
      END DO
    END DO

    ! Diagonal
    FORALL (i=1:Npar) corr2Rmat_DP(i,i) = 1._DP

    !------------------------------------------------------------------------

  END FUNCTION corr2Rmat_DP


  !==========================================================================
  ! corr[Ncorr] = VMAT2CORR(Vmat[Npar,Npar],Npar,Ncorr)
  !
  !   Extract the correlation coefficients from the covariance matrix.
  !==========================================================================

  PURE FUNCTION Vmat2corr (Vmat,Npar,Ncorr)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(Npar,Npar), INTENT(IN) :: Vmat
    INTEGER, INTENT(IN) :: Npar, Ncorr
    REAL(DP), DIMENSION(Ncorr) :: Vmat2corr

    INTEGER :: i, j, k
    REAL(DP), DIMENSION(Npar) :: sig

    !------------------------------------------------------------------------

    FORALL (i=1:Npar) sig(i) = SQRT(Vmat(i,i))
    k = 1
    DO i=1,Npar
      DO j=1,Npar
        IF (MERGE(j > i,i > j,upper)) THEN
          Vmat2corr(k) = Vmat(i,j) / sig(i) / sig(j)
          k = k + 1
        END IF
      END DO
    END DO

    !------------------------------------------------------------------------

  END FUNCTION Vmat2corr


  !==========================================================================
  ! CALL CORREL_INDEX(Npar,Ncorr,icorr2ij)
  !
  !   Returns the pair of parameter inidices corresponding to unique correlation
  ! coefficients.
  !==========================================================================

  SUBROUTINE correl_index (Npar,Ncorr,icorr2ij)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: Npar, Ncorr
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: icorr2ij

    INTEGER :: i, j, k

    !-----------------------------------------------------------------------

    ALLOCATE (icorr2ij(Ncorr,2))
    k = 1
    DO i=1,Npar
      DO j=1,Npar
        IF (MERGE(j > i,i > j,upper)) THEN
          icorr2ij(k,:) = [ i, j ]
          k = k + 1
        END IF
      END DO
    END DO

    !-----------------------------------------------------------------------

  END SUBROUTINE correl_index


  !==========================================================================
  ! Ncorr = N_CORR(Npar)
  !
  !   Number of unique correlations for Npar parameters.
  !==========================================================================

  FUNCTION N_corr (Npar)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: Npar
    INTEGER :: N_corr

    !------------------------------------------------------------------------

    N_corr = NINT( REAL(Npar*(Npar-1))/2._DP )

    !------------------------------------------------------------------------

  END FUNCTION N_corr


  !==========================================================================
  ! r = CORRELATE(x[N],y[N])
  !
  !   Computes the Pearson correlation coefficient of two arrays. 
  !==========================================================================

  PURE FUNCTION correlate_4D (x,y,mask)

    USE utilities, ONLY: DP, tinyDP
    USE special_functions, ONLY: ibeta
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: x, y
    LOGICAL, DIMENSION(:,:,:,:), INTENT(IN), OPTIONAL :: mask
    REAL(DP) :: correlate_4D

    REAL(DP) :: N, ax, ay, sxx, sxy, syy
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3),SIZE(x,4)) :: mask0
    
    !-----------------------------------------------------------------------

    mask0(:,:,:,:) = .True.
    IF (PRESENT(mask)) mask0(:,:,:,:) = mask(:,:,:,:)
    N = COUNT(mask0(:,:,:,:))
    ax = SUM(x(:,:,:,:),MASK=mask0(:,:,:,:)) / N
    ay = SUM(y(:,:,:,:),MASK=mask0(:,:,:,:)) / N
    sxx = SUM((x(:,:,:,:)-ax)**2,MASK=mask0(:,:,:,:)) / (N-1)
    syy = SUM((y(:,:,:,:)-ay)**2,MASK=mask0(:,:,:,:)) / (N-1)
    sxy = SUM((x(:,:,:,:)-ax)*(y(:,:,:,:)-ay),MASK=mask0(:,:,:,:)) / (N-1)
    correlate_4D = sxy / (SQRT(sxx*syy) + tinyDP)

    !-----------------------------------------------------------------------

  END FUNCTION correlate_4D

  ! 3D interface
  PURE FUNCTION correlate_3D(x,y,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: x, y
    LOGICAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: mask
    REAL(DP) :: correlate_3D
    INTEGER :: N1, N2, N3
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2),SIZE(x,3)) :: mask0
    N1 = SIZE(x(:,:,:),1)
    N2 = SIZE(x(:,:,:),2)
    N3 = SIZE(x(:,:,:),3)
    mask0(:,:,:) = .True.
    IF (PRESENT(mask)) mask0(:,:,:) = mask(:,:,:)
    correlate_3D = CORRELATE_4D(RESHAPE(x(:,:,:),[N1,N2,N3,1]), &
                                RESHAPE(y(:,:,:),[N1,N2,N3,1]), &
                                MASK=RESHAPE(mask0(:,:,:),[N1,N2,N3,1]))
  END FUNCTION correlate_3D

  ! 2D interface
  PURE FUNCTION correlate_2D(x,y,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: x, y
    LOGICAL, DIMENSION(:,:), INTENT(IN), OPTIONAL :: mask
    REAL(DP) :: correlate_2D
    INTEGER :: N1, N2
    LOGICAL, DIMENSION(SIZE(x,1),SIZE(x,2)) :: mask0
    N1 = SIZE(x(:,:),1)
    N2 = SIZE(x(:,:),2)
    mask0(:,:) = .True.
    IF (PRESENT(mask)) mask0(:,:) = mask(:,:)
    correlate_2D = CORRELATE_4D(RESHAPE(x(:,:),[N1,N2,1,1]), &
                                RESHAPE(y(:,:),[N1,N2,1,1]), &
                                MASK=RESHAPE(mask0(:,:),[N1,N2,1,1]))
  END FUNCTION correlate_2D

  ! 1D interface
  PURE FUNCTION correlate_1D(x,y,mask)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
    LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: mask
    REAL(DP) :: correlate_1D
    INTEGER :: N1
    LOGICAL, DIMENSION(SIZE(x)) :: mask0
    N1 = SIZE(x(:),1)
    mask0(:) = .True.
    IF (PRESENT(mask)) mask0(:) = mask(:)
    correlate_1D = CORRELATE_4D(RESHAPE(x(:),[N1,1,1,1]), &
                                RESHAPE(y(:),[N1,1,1,1]), &
                                MASK=RESHAPE(mask0(:),[N1,1,1,1]))
  END FUNCTION correlate_1D


  !==========================================================================
  ! CALL CORREL_PARLIST(x[N,Npar]IN,rho[N_CORR(Npar)]OUT)
  !
  !   Computes all the non-trivial correlation coefficients among a list of 
  ! parameters. If the parmaters are CHARACTER, then the correlation name is
  ! built. 
  !==========================================================================

  SUBROUTINE correl_parlist_DP (par,corr)
  
    USE utilities, ONLY: DP
    IMPLICIT NONE
   
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: par
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: corr

    INTEGER :: i, j, N, Npar, Ncorr
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Rmat

    !-----------------------------------------------------------------------

    N = SIZE(par(:,:),1)
    Npar = SIZE(par(:,:),2)
    Ncorr = N_CORR(Npar)
    ALLOCATE (Rmat(Npar,Npar))
    FORALL (i=1:Npar,j=1:Npar,MERGE(j > i,i > j,upper)) &
      Rmat(i,j) = CORRELATE(par(:,i),par(:,j))
    ALLOCATE (corr(Ncorr))
    corr(:) = RMAT2CORR(Rmat(:,:),Npar,Ncorr)
    DEALLOCATE (Rmat)

    !-----------------------------------------------------------------------

  END SUBROUTINE correl_parlist_DP

  !-------------------------------------------------------------------------

  SUBROUTINE correl_parlist_char (par,corr)
  
    USE utilities, ONLY: DP, trimlr, lenmax
    IMPLICIT NONE
   
    CHARACTER(*), DIMENSION(:), INTENT(IN) :: par
    CHARACTER(*), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: corr

    INTEGER :: i, j, Npar, Ncorr
    CHARACTER(lenmax), DIMENSION(:,:), ALLOCATABLE :: Rmat

    !-----------------------------------------------------------------------

    Npar = SIZE(par(:))
    Ncorr = N_CORR(Npar)
    ALLOCATE (Rmat(Npar,Npar))
    FORALL (i=1:Npar,j=1:Npar,MERGE(j > i,i > j,upper)) &
      Rmat(i,j) = TRIMLR(par(i))//"-"//TRIMLR(par(j))    
    ALLOCATE (corr(Ncorr))
    corr(:) = RMAT2CORR(Rmat(:,:),Npar,Ncorr)
    DEALLOCATE (Rmat)    

    !-----------------------------------------------------------------------

  END SUBROUTINE correl_parlist_char


  !==========================================================================
  ! CALL AUTOCORREL(data[N],rho[Nlag],Nlag)
  !
  !   Computes the reduced autocorrelation function of a time serie, as a
  ! a function of positive lag [1,Nlag].
  !==========================================================================

  SUBROUTINE autocorrel (data,rho,Nlag)
   
    USE utilities, ONLY: DP
    USE FFT_specials, ONLY: correl
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: data
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: rho
    INTEGER, INTENT(OUT), OPTIONAL :: Nlag

    INTEGER :: Ndata, Nacf, Nl
    REAL(DP) :: meandata, vardata
    REAL(DP), DIMENSION(:), ALLOCATABLE :: acf, dat

    !-------------------------------------------------------------------------

    ! Properties of the time serie
    Ndata = SIZE(data)
    meandata = MEAN(data)
    vardata = (SIGMA(data))**2

    ! Raw autocorrelation
    Nacf = 2**(FLOOR(LOG(REAL(Ndata,DP))/LOG(2._DP)))
    ALLOCATE (acf(Nacf),dat(Nacf))
    dat(:) = data(1:Nacf) - meandata
    acf(:) = CORREL(dat(:),dat(:))
    
    ! Normalization
    Nl = Nacf / 2
    ALLOCATE (rho(Nl))
!    rho(:) = ( acf(1:Nl) - meandata**2 ) / vardata
    rho(:) = acf(1:Nl) / acf(1)
    IF (PRESENT(Nlag)) Nlag = Nl
    DEALLOCATE(acf,dat)

    !-------------------------------------------------------------------------

  END SUBROUTINE autocorrel


  !==========================================================================
  ! t_int = INTAUTOCORRTIME(acf[Nlag],Nmcmc,Neff(OUT,OPT),low(OPT),high(OPT), &
  !                         step(OPT),c(OPT),FORCE=T/F)
  !
  !   Computes the integrated autocorrelation time, following the method by
  ! Sokal (http://www.stat.unc.edu/faculty/cji/Sokal.pdf), as implemented in
  ! Emcee (https://github.com/dfm/emcee/blob/master/emcee/autocorr.py).
  !   RHO[Nlag] is the precomputed ACF and Nmcmc the length of the chain, both
  ! after burn-in.
  !   LOW is the minimum window size to test (default = 10).
  !   HIGH is the maximum window size to test (default = Nmcmc/(2c)).
  !   STEP is the step size for the window search (default = 1).
  !   C is the minimum number of autocorrelation times needed to trust the
  ! estimate (default = 10).
  !==========================================================================

  FUNCTION intautocorrtime (acf,Nmcmc,Neff,low,high,step,c,warn,force)

    USE utilities, ONLY: DP, warning, warnings, hugeDP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: acf
    INTEGER, INTENT(IN) :: Nmcmc
    INTEGER, INTENT(OUT), OPTIONAL :: Neff
    INTEGER, INTENT(IN), OPTIONAL :: low, high, step
    REAL(DP), INTENT(IN), OPTIONAL :: c
    LOGICAL, INTENT(IN), OPTIONAL :: warn, force
    REAL(DP) :: intautocorrtime

    INTEGER :: siz, low0, high0, step0, m
    REAL(DP) :: c0, tau
    LOGICAL :: warn0, force0
    
    !-------------------------------------------------------------------------
    
    ! Default values
    siz = NINT(Nmcmc/2._DP)
    IF (PRESENT(warn)) THEN ; warn0 = warn ; ELSE ; warn0 = warnings ; END IF
    IF (PRESENT(force)) THEN ; force0 = force ; ELSE ; force0 = .False. ; END IF
    IF (PRESENT(low)) THEN ; low0 = low ; ELSE ; low0 = 10 ; END IF
    IF (PRESENT(c)) THEN ; c0 = c ; ELSE ; c0 = 10._DP ; END IF
    IF (PRESENT(high)) THEN ; high0 = high
                       ELSE ; high0 = NINT(Nmcmc/(2*c0)) ; END IF
    IF (PRESENT(step)) THEN ; step0 = step ; ELSE ; step0 = 1 ; END IF
    IF (NINT(c0*low0) >= siz) THEN
      IF (warn0) CALL WARNING("INTAUTOCORRTIME","The chain is too short")
      intautocorrtime = hugeDP
      IF (PRESENT(Neff)) Neff = 0
      RETURN
    END IF
    
    ! Loop over proposed window sizes until convergence is reached
    DO m=low0,high0,step0
      ! Compute the autocorrelation time with the given window
      tau = 1 + 2 * SUM(acf(1:m)) / acf(1)
      ! Accept the window size if it satisfies the convergence criterion
      IF (tau > 1._DP .AND. m > c0*tau) THEN
        intautocorrtime = tau
        IF (PRESENT(Neff)) Neff = NINT(Nmcmc/tau)
        RETURN
      ELSE IF (c0*tau >= siz) THEN 
        IF (warn0) CALL WARNING("INTAUTOCORRTIME", &
                            "The chain is too short to reliably estimate the " &
                            //"autocorrelation time")
        IF (force0) THEN ! make t_int negative to flag it as unreliable
          intautocorrtime = - ( 1 + 2 * SUM(acf(:)) / acf(1) ) 
          IF (PRESENT(Neff)) Neff = NINT(Nmcmc/intautocorrtime)
        ELSE 
          intautocorrtime = hugeDP
          IF (PRESENT(Neff)) Neff = 0
        END IF
        RETURN
      END IF
    END DO

    ! To quiet the compiler
    IF (force0) THEN ! make t_int negative to flag it as unreliable
      intautocorrtime = - ( 1 + 2 * SUM(acf(:)) / acf(1) ) 
      IF (PRESENT(Neff)) Neff = NINT(Nmcmc/intautocorrtime)
    ELSE 
      intautocorrtime = hugeDP
      IF (PRESENT(Neff)) Neff = 0
    END IF
    RETURN
    
    !-------------------------------------------------------------------------

  END FUNCTION intautocorrtime

  
END MODULE statistics
