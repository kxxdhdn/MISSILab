PROGRAM test_invert

  USE auxil, ONLY: check_SM, invert_SM, invert_mSM
  USE utilities, ONLY: DP
  USE matrices, ONLY: invert_cholesky
  IMPLICIT NONE

  INTEGER, DIMENSION(2) :: pos
  REAL(DP) :: delta
  REAL(DP), DIMENSION(4,4) :: C0, invC0, C, invC
  REAL(DP), DIMENSION(3,3) :: A0, invA0, A1, A, invA, invA_ch

  !! Sherman-Morrison (ex via Sherman&Morrison50)
  !!------------------
  C0(1,:) = [2.384_DP, 1.238_DP, 0.861_DP, 2.413_DP]
  C0(2,:) = [0.648_DP, 1.113_DP, 0.761_DP, 0.137_DP]
  C0(3,:) = [1.119_DP, 0.643_DP, 3.172_DP, 1.139_DP]
  C0(4,:) = [0.745_DP, 2.137_DP, 1.268_DP, 0.542_DP]

  C(:,:) = C0(:,:)
  C(2,4) = C0(2,4) + 0.4

  invC0(1,:) = [0.2220_DP, 2.5275_DP, -0.1012_DP, -1.4145_DP]
  invC0(2,:) = [-0.04806_DP, -0.2918_DP, -0.1999_DP, 0.7079_DP]
  invC0(3,:) = [-0.1692_DP, 0.01195_DP, 0.3656_DP, -0.01824_DP]
  invC0(4,:) = [0.2801_DP, -2.3517_DP, 0.07209_DP, 1.0409_DP]

  CALL check_SM(C, C0, pos, delta)
  invC(:,:) = invert_SM(C, invC0, pos, delta)
  PRINT*, '- invC via Sherman-Morrison: '
  PRINT*, invC
  PRINT*

  !! Modified Sherman-Morrison
  !!---------------------------
  A0(1,:) = [2._DP, -1._DP, 0._DP]
  A0(2,:) = [-1._DP, 2._DP, -1._DP]
  A0(3,:) = [0._DP, -1._DP, 2._DP]

  A1(:,:) = A0(:,:)
  A1(2,3) = A0(2,3) + 0.4
  PRINT*, '- A1: '
  PRINT*, A1
  PRINT*

  A(:,:) = A1(:,:)
  A(3,2) = A1(3,2) + 0.4
  PRINT*, '- A: '
  PRINT*, A
  PRINT*
  
  invA0(:,:) = invert_cholesky(A0)
  PRINT*, '- invA0 via Cholesky: '
  PRINT*, invA0
  PRINT*

  CALL check_SM(A, A1, pos, delta)
  PRINT*, 'pos = ', pos
  PRINT*, 'delta = ', delta
  PRINT*
  invA(:,:) = invert_mSM(A, invA0, pos, delta)
  invA_ch(:,:) = invert_cholesky(A)
  PRINT*, '- invA via modified Sherman-Morrison: '
  PRINT*, invA
  PRINT*
  PRINT*, '- discrepancy between invA and invA_ch: '
  PRINT*, invA - invA_ch
  PRINT*
  
END PROGRAM test_invert
