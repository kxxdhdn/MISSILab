PROGRAM test_modifBB

  USE auxil
  USE utilities, ONLY: DP
  USE inout, ONLY: write_hdf5, h5ext, write_ascii, ascext
  IMPLICIT NONE

  INTEGER Nspec
  ! CHARACTER(30), DIMENSION(50)              :: keyIN
  CHARACTER(30), DIMENSION(:), ALLOCATABLE  :: label
  TYPE(Qabs_STR), DIMENSION(:), ALLOCATABLE :: nQ_str
  REAL(DP) :: distance = 1._DP ! kpc
  ! REAL(DP), DIMENSION(2,50)                 :: key2IN
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: M, T
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: Fnu

  !! QUICK INPUTS
  INTEGER :: n = 1
  ALLOCATE(label(n), nQ_str(n), M(n), T(n))
  label = (/'Sil_D03 '/)
  M = 3._DP
  T = 100._DP

  !! Inputs -> To create config file with Python interface
  ! DO i=1,50
  !   PRINT*, "Enter a dust component: "
  !   READ(*,*) keyIN(i)
  !   IF (keyIN(i)=='exit') THEN
  !     !! Finish input
  !     IF (i.EQ.1) THEN
  !       PRINT*, "No input. "
  !       STOP

  !     END IF
  !     EXIT ! Only way to exit (except ERRORs)

  !   END IF
  !   PRINT*, "Enter the mass (Msun): "
  !   READ(*,*) key2IN(1,i)
  !   PRINT*, "Enter the temperature (K): "
  !   READ(*,*) key2IN(2,i)
    
  ! END DO
  ! ALLOCATE(label(i-1), M(i-1), T(i-1)) ! exclude 'exit'
  ! label = keyIN(i-1)
  ! M = key2IN(1,i-1)
  ! T = key2IN(2,i-1)

  Nspec=SIZE(label)
  PRINT*, "Number of components: ", Nspec

  !! Get param
  CALL optics_INIT(label, nQ_str)
  !! Calculate Fnu
  Fnu = M(1)/distance**2 * modifBB(T(1), nQ_str(1))
  
  CALL write_hdf5(nQ_str(1)%wave, FILE="output/test_modifBB"//h5ext, &
                  COMPRESS=.False., NAME="Wavelength (micron)")
  CALL write_hdf5(nQ_str(1)%nu, FILE="output/test_modifBB"//h5ext, &
                  COMPRESS=.False., NAME="Wavelength (Hz)", &
                  APPEND=.True.)
  CALL write_hdf5(Fnu, FILE="output/test_modifBB"//h5ext, &
                  COMPRESS=.False., NAME="Fnu (MJy/sr)", &
                  APPEND=.True.)
  CALL write_ascii(VEC1=nQ_str(1)%wave, VEC2=Fnu, &
                   FILE="output/test_modifBB"//ascext)

  !! Free memory space
  DEALLOCATE(label, M, T, Fnu)

END PROGRAM test_modifBB
