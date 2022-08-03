!******************************************************************************
!*
!*                        OPERATIONS ON MATRICES
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 11/2012
  ! 
  ! 3) DESCRIPTION: implements various operations on matrices.
  !==========================================================================

MODULE matrices

  USE utilities, ONLY:
  USE linear_system, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: determinant_matrix, invert_matrix
  PUBLIC :: invert_cholesky, determinant_cholesky


CONTAINS


  !==========================================================================
  ! detA = DETERMINANT_LU(A[N,N][,NODECOMPOSITION])
  !
  !   Compute the determinant DETA of a square matrix A, using the LU 
  ! decomposition method. The keyword NODECOMPOSITION has to be used
  ! if the input matrix has already been decomposed.
  !==========================================================================

  FUNCTION determinant_LU (A,d,nodecomposition)

    USE utilities, ONLY: DP
    USE linear_system, ONLY: LU_decomp
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:,:) :: A
    INTEGER, INTENT(IN), OPTIONAL :: d
    LOGICAL, INTENT(IN), OPTIONAL :: nodecomposition
    REAL(DP) :: determinant_LU

    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: LU
    REAL(DP), DIMENSION(SIZE(A,1)) :: diag
    LOGICAL :: decomp
    INTEGER, DIMENSION(SIZE(A,1)) :: indx
    INTEGER :: N, i, dd

    !-----------------------------------------------------------------------

    N = SIZE(A,1)
    decomp = .True.
    IF (PRESENT(nodecomposition)) decomp = MERGE(.False.,.True.,nodecomposition)
    IF (decomp) THEN 
      LU(:,:) = LU_DECOMP(A,indx,dd)
    ELSE 
     LU(:,:) = A 
     dd = d
    END IF
    FORALL (i=1:N) diag(i) = LU(i,i)
    determinant_LU = REAL(dd,DP) * PRODUCT(diag)

    !-----------------------------------------------------------------------

  END FUNCTION determinant_LU


  !==========================================================================
  ! detA = DETERMINANT_CHOLESKY(A[N,N][,NODECOMPOSITION,NOPOSDEF])
  !
  !   Compute the determinant DETA of a positive definite matrix A, using the 
  ! Cholesky decomposition method. The keyword NODECOMPOSITION has to be used
  ! if the input matrix has already been decomposed.
  !==========================================================================

  FUNCTION determinant_cholesky (A,nodecomposition,noposdef)

    USE utilities, ONLY: DP, NaN
    USE linear_system, ONLY: cholesky_decomp
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: A
    LOGICAL, INTENT(IN), OPTIONAL :: nodecomposition
    LOGICAL, INTENT(OUT), OPTIONAL :: noposdef
    REAL(DP) :: determinant_cholesky

    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Chdcp
    REAL(DP), DIMENSION(SIZE(A,1)) :: diag
    LOGICAL :: decomp
    INTEGER :: i, N

    !-----------------------------------------------------------------------

    ! Cholesky decomposition
    N = SIZE(A,1)
    decomp = .True.
    IF (PRESENT(noposdef)) noposdef = .False.
    IF (PRESENT(nodecomposition)) decomp = MERGE(.False.,.True.,nodecomposition)
    IF (decomp) THEN 
      IF (PRESENT(noposdef)) THEN 
        Chdcp = CHOLESKY_DECOMP(A,noposdef) 
        IF (noposdef) THEN
          determinant_cholesky = NaN()
          RETURN
        END IF
      ELSE
        Chdcp = CHOLESKY_DECOMP(A) 
      END IF
    ELSE 
      Chdcp = A
    END IF

    ! Compute the determinant
    FORALL (i=1:N) diag(i) = Chdcp(i,i)
    determinant_cholesky = EXP(2._DP*SUM(LOG(diag)))
    
    !-----------------------------------------------------------------------

  END FUNCTION determinant_cholesky


  !==========================================================================
  ! detA = DETERMINANT_MATRIX(A[N,N][,LU,Cholesky])
  !
  !   Overload of the computation of a determinant matrix A[N,N].
  !==========================================================================

  FUNCTION determinant_matrix (A,lu,Cholesky)
 
    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:,:) :: A
    LOGICAL, INTENT(IN), OPTIONAL :: lu, Cholesky
    REAL(DP) :: determinant_matrix

    CHARACTER(8) :: method

    !-----------------------------------------------------------------------

    method = "LU      "
    IF (PRESENT(LU)) method = MERGE("LU      ","Cholesky",LU)
    IF (PRESENT(Cholesky)) method = MERGE("Cholesky",method,Cholesky)

    SELECT CASE (method)
      CASE ("LU      ")
        determinant_matrix = DETERMINANT_LU(A)
      CASE ("Cholesky")
        determinant_matrix = DETERMINANT_CHOLESKY(A)
      CASE DEFAULT
        determinant_matrix = 0._DP
    END SELECT

    !-----------------------------------------------------------------------

  END FUNCTION determinant_matrix


  !==========================================================================
  ! invA[N,N] = INVERT_GAUSS(A[N,N])
  !
  !   Compute the inverse INVA of a square matrix A, with Gauss-Jordan 
  ! elimination and full pivot.
  !==========================================================================

  FUNCTION invert_gauss (A)

    USE utilities, ONLY: DP
    USE linear_system, ONLY: linsyst_gauss
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:,:) :: A
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: invert_gauss

    REAL(DP), DIMENSION(SIZE(A,1),1) :: B

    !-----------------------------------------------------------------------
   
    B(:,:) = 0._DP
    B = LINSYST_GAUSS(A,B,invert_gauss,ONLY_INVERSE=.True.)

    !-----------------------------------------------------------------------

  END FUNCTION invert_gauss


  !==========================================================================
  ! invA[N,N] = INVERT_LU(A[N,N][,determinant])
  !
  !   Compute the inverse INVA of a square matrix A, with Gauss-Jordan 
  ! elimination and full pivot.
  !==========================================================================

  FUNCTION invert_LU (A,determinant)

    USE utilities, ONLY: DP
    USE linear_system, ONLY: LU_decomp, linsyst_LU
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:,:) :: A
    REAL(DP), INTENT(OUT), OPTIONAL :: determinant
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: invert_LU

    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: B, LU
    INTEGER, DIMENSION(SIZE(A,1)) :: indx
    INTEGER :: i, N, d

    !-----------------------------------------------------------------------
   
    ! Solve for the identity matrix
    N = SIZE(A,1)
    B(:,:) = 0._DP
    FORALL (i=1:N) B(i,i) = 1._DP
    LU(:,:) = LU_DECOMP(A,indx,d)

    ! Loop on the number of systems
    DO i=1,N
      invert_LU(:,i) = LINSYST_LU(LU,B(:,i),indx)
    END DO

    ! Determinant
    IF (PRESENT(determinant)) &
      determinant = DETERMINANT_LU(LU,NODECOMPOSITION=.True.)

    !-----------------------------------------------------------------------

  END FUNCTION invert_LU


  !==========================================================================
  ! invA[N,N] = INVERT_CHOLESKY(A[N,N][,determinant,noposdef])
  !
  !   Compute the inverse INVA of a positive definite matrix A, with Cholesky
  ! decomposition.
  !==========================================================================

  FUNCTION invert_cholesky (A,determinant,noposdef)

    USE utilities, ONLY: DP
    USE linear_system, ONLY: cholesky_decomp
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(IN) :: A
    REAL(DP), INTENT(OUT), OPTIONAL :: determinant
    LOGICAL, INTENT(OUT), OPTIONAL :: noposdef
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: invert_cholesky

    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Chdcp
    INTEGER :: i, j, N

    !-----------------------------------------------------------------------

    ! Perform the Cholesky decomposition
    N = SIZE(A,1)
    IF (PRESENT(noposdef)) THEN
      Chdcp = CHOLESKY_DECOMP(A,noposdef)
      IF (noposdef) RETURN
    ELSE
      Chdcp = CHOLESKY_DECOMP(A)
    END IF  
  
    ! Compute the inverse matrix
    row_low: DO i=1,N 
      FORALL (j=1:i) &
        invert_cholesky(j,i) = (MERGE(1._DP,0._DP,i==j) &
                       - DOT_PRODUCT(Chdcp(i,j:i-1),invert_cholesky(j,j:i-1))) &
                       / Chdcp(i,i)
    END DO row_low
    row_up: DO i=N,1,-1
      FORALL (j=1:i)
        invert_cholesky(j,i) = (MERGE(0._DP,invert_cholesky(j,i),i<j) &
                       - DOT_PRODUCT(Chdcp(i+1:N,i),invert_cholesky(j,i+1:N))) &
                      / Chdcp(i,i)
        invert_cholesky(i,j) = invert_cholesky(j,i)
      END FORALL
    END DO row_up

    ! Determinant
    IF (PRESENT(determinant)) &
      determinant = DETERMINANT_CHOLESKY(Chdcp,NODECOMPOSITION=.True.)

    !-----------------------------------------------------------------------
    
  END FUNCTION invert_cholesky


  !==========================================================================
  ! invA[N,N] = INVERT_MATRIX(A[N,N][,gauss,lu,cholesky,determinant])
  !
  !   Overload the inverse of a matrix.
  !==========================================================================

  FUNCTION invert_matrix (A,gauss,LU,Cholesky,determinant)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:,:) :: A
    REAL(DP), INTENT(OUT), OPTIONAL :: determinant
    LOGICAL, INTENT(IN), OPTIONAL :: gauss, LU, Cholesky
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: invert_matrix

    CHARACTER(8) :: method

    !-----------------------------------------------------------------------
 
    ! Choice of the method   
    method = "Gauss"
    IF (PRESENT(LU)) method = MERGE("LU      ",method,LU)
    IF (PRESENT(gauss)) method = MERGE("Gauss   ",method,gauss)
    IF (PRESENT(Cholesky)) method = MERGE("Cholesky",method,Cholesky)

    ! Compute the inverse
    SELECT CASE (method)
      CASE ("Gauss   ")
        invert_matrix = INVERT_GAUSS(A)
      CASE ("LU      ")
        IF (PRESENT(determinant)) THEN 
          invert_matrix = INVERT_LU(A,determinant)
        ELSE 
          invert_matrix = INVERT_LU(A)
        END IF
      CASE ("Cholesky")
        IF (PRESENT(determinant)) THEN 
          invert_matrix = INVERT_CHOLESKY(A,determinant)
        ELSE 
          invert_matrix = INVERT_CHOLESKY(A)
        END IF
    END SELECT

    !-----------------------------------------------------------------------

  END FUNCTION invert_matrix


END MODULE matrices
