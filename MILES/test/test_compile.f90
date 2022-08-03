PROGRAM test_compile
  
  USE utilities, ONLY: DP, pring, trimLR
  USE arrays, ONLY: reallocate
  USE constants, ONLY: pi
  IMPLICIT NONE
  ! PRIVATE

  ! PUBLIC
  INTEGER i, N
  REAL(DP) :: a, b
  CHARACTER(30), DIMENSION(50) :: cIN
  CHARACTER(30), DIMENSION(:), ALLOCATABLE :: c

  a = 3.14_DP
  b = pi
  PRINT*, "a = "//TRIMLR(PRING(a))
  PRINT*, "b = "//TRIMLR(PRING(b))

  !! Append inputs
  PRINT*, "Maximum 50 inputs; 'exit' to quit "
  ALLOCATE(c(0))
  DO i=1,50
    READ(*,*) cIN(i)
    IF (cIN(i)=='exit') THEN
      !! Finish input
      IF (i.EQ.1) THEN
        PRINT*, "No input. "

      END IF
      EXIT !! Only way to exit (except ERRORs)

    END IF
  END DO
  CALL reallocate(c, i)
  c = cIN(1:i-1)
  N=SIZE(c)
  PRINT*, "The size of c is ", N
  PRINT*, "c = ", c

  PRINT*, "Bravo! "

END PROGRAM test_compile
