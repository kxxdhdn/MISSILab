PROGRAM test_profiles

  USE auxil
  USE arrays, ONLY: ramp
  USE utilities, ONLY: DP
  ! USE constants, ONLY: MKS
  USE integration, ONLY: integ_tabulated
  USE inout, ONLY: write_hdf5, h5ext, write_ascii, ascext
  IMPLICIT NONE

  ! INTEGER i
  CHARACTER(30), DIMENSION(:), ALLOCATABLE  :: label
  TYPE(Qabs_STR), DIMENSION(:), ALLOCATABLE :: nQ_str
  REAL(DP) :: distance, sigma, sigmaS, sigmaL, ref1, ref2
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: massBB, tempBB
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: FnuLINE, FnuBAND
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: FnuBB, Fnu
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: FnuTEST

  !! QUICK INPUTS
  !!--------------
  INTEGER :: Nbb = 2
  ! INTEGER :: Nline = 1
  ! INTEGER :: Nband = 1
  ! PRINT*, "Number of components (dust continuum): ", Nbb
  ! PRINT*, "Number of lines: ", Nline
  ! PRINT*, "Number of bands: ", Nband

  !! Integration test (homogeneous x)
  !!----------------------------------
  INTEGER, PARAMETER :: Nt = 1000
  REAL(DP), DIMENSION(:), ALLOCATABLE :: xt
  ALLOCATE(xt(Nt))
  xt(:) = ramp(Nt, 1._DP, 40._DP)
  ! FnuTEST = gaussLine(xt, 15.555_DP, 15.555_DP*0.0055722841_DP)
  FnuTEST = lorentzBand(xt, 6.2_DP, 0.060000000_DP, 0.031300317_DP)
  print*, 'Should be 1 (homo x): ', integ_tabulated(xt, FnuTEST)
  !!----------------------------------------------------------------

  ALLOCATE(label(Nbb), nQ_str(Nbb), massBB(Nbb), tempBB(Nbb))
  
  label = (/'Sil_D03      ', 'PAHn_DL07_G11'/)
  massBB = 3._DP
  tempBB = 200._DP
  distance = 1._DP ! kpc
  ref1 = 15.555_DP ! [NeIII] 1 line
  ref2 = 6.2_DP ! 6.2 micron band
  sigma = ref1 * 0.0055722841_DP
  sigmaL = 0.060000000_DP
  sigmaS = 0.031300317_DP

  !! GET PARAM
  !!-----------
  CALL optics_INIT(label, nQ_str)

  !! Calculate Fnu
  !!---------------
  FnuBB = massBB(1)/distance**2 * modifBB(tempBB(1), nQ_str(1))
  FnuLINE = gaussLine(nQ_str(1)%wave, ref1, sigma)
  FnuBAND = lorentzBand(nQ_str(1)%wave, ref2, sigmaS, sigmaL)
  PRINT*, 'Should be 1 (gauss): ', integ_tabulated(nQ_str(1)%nu, FnuLINE, &
          XLOG=.True., YLOG=.False.)
  PRINT*, 'Should be 1 (lorentz): ', integ_tabulated(nQ_str(1)%nu, FnuBAND, &
          XLOG=.True., YLOG=.False.)
  Fnu = FnuBB + FnuLINE + FnuBAND
  
  !! OUTPUTS
  !!---------
  CALL write_hdf5(nQ_str(1)%wave, NAME="Wavelength (micron)", &
                  COMPRESS=.False., APPEND=.False., &
                  FILE="output/test_profiles"//h5ext)
  CALL write_hdf5(nQ_str(1)%nu, NAME="Wavelength (Hz)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/test_profiles"//h5ext)
  CALL write_hdf5(FnuBB, NAME="FnuBB (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/test_profiles"//h5ext)
  CALL write_hdf5(FnuLINE, NAME="FnuLINE (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/test_profiles"//h5ext)
  CALL write_hdf5(FnuBAND, NAME="FnuBAND (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/test_profiles"//h5ext)
  CALL write_hdf5(Fnu, NAME="Fnu (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/test_profiles"//h5ext)
  CALL write_hdf5(FnuTEST, NAME="FnuTEST (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/test_profiles"//h5ext)
  CALL write_hdf5(xt, NAME="x", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/test_profiles"//h5ext)

  !! Free memory space
  DEALLOCATE(label, massBB, tempBB, &
             FnuLINE, FnuBAND, FnuBB, Fnu, FnuTEST)

END PROGRAM test_profiles
