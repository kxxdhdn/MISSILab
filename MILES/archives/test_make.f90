PROGRAM test_make

  USE auxiltest
  USE utilities, ONLY: DP
  IMPLICIT NONE

  REAL(DP) :: a=22._DP
  CHARACTER(15) b
  
  CALL coucou(a, b)
  PRINT*, 'OUT', a, b
  
END PROGRAM test_make
