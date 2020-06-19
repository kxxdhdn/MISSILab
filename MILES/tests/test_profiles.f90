PROGRAM test_profiles

  USE auxil
  USE arrays, ONLY: ramp, reverse
  USE utilities, ONLY: DP, TRIMLR
  USE constants, ONLY: MKS
  USE integration, ONLY: integ_tabulated
  USE interpolation, ONLY: interp_lin_sorted
  USE inout, ONLY: write_hdf5, h5ext, write_ascii, ascext
  IMPLICIT NONE

  ! INTEGER i
  CHARACTER(30)                             :: path
  CHARACTER(30), DIMENSION(:), ALLOCATABLE  :: label
  TYPE(Qabs_STR), DIMENSION(:), ALLOCATABLE :: nQabs
  REAL(DP) :: sigma, sigmaS, sigmaL, ref1, ref2, Iline, Iband
  INTEGER, PARAMETER                        :: Nw = 1000
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: wave, nu
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: Pabs
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: massBB, tempBB
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: FnuLINE, FnuBAND
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: FnuBB, Fnu
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: FnuTEST

  !! QUICK INPUTS
  !!--------------
  INTEGER :: Nbb = 2
  ! INTEGER :: Nline = 1
  ! INTEGER :: Nband = 1

  !! Integration test (homogeneous x)
  !!----------------------------------
  INTEGER, PARAMETER                        :: NwTEST = 582 !+ 1000
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: wTEST, nuTEST
  ALLOCATE(wTEST(NwTEST), nuTEST(NwTEST))
  
  !! nu + log
  !!----------
  ! nuTEST(:) = ramp(NwTEST, 1.E12_DP, 1.E14_DP, .true.) ! this sampling may miss the pic
  ! FnuTEST = gaussLine(nuTEST, 1.9E13_DP, 1.E11_DP) ! /Hz
  ! FnuTEST = lorentzBand(nuTEST, 1.9E13_DP, 1.E11_DP, 5.E10_DP)
  ! PRINT*, 'I(test) = ', integ_tabulated(nuTEST, FnuTEST, xlog=.true.)

  !! w
  !!---
  ! wTEST(:) = ramp(NwTEST, 1._DP, 40._DP)
  ! wTEST(:) = ramp(NwTEST, 1.E-2_DP, 3.E4_DP) ! undersampling WARNING
  ! nuTEST = MKS%clight/MKS%micron / wTEST
  ! FnuTEST = gaussLine(wTEST, 15.555_DP, 15.555_DP*0.0055722841_DP) ! /um
  ! FnuTEST = lorentzBand(wTEST, 6.2_DP, 0.060000000_DP, 0.031300317_DP) ! /um
  ! PRINT*, 'I(test) = ', integ_tabulated(wTEST, FnuTEST)

  !! w + log
  !!---------
  ! wTEST(:) = ramp(NwTEST, 1.E-2_DP, 3.E4_DP, .true.) ! try larger NwTEST
  ! FnuTEST = gaussLine(wTEST, 15.555_DP, 15.555_DP*0.0055722841_DP) ! /um
  ! PRINT*, 'I(test) = ', integ_tabulated(wTEST, FnuTEST, xlog=.true.)

  !! w -> nu
  !!---------
  ! wTEST(:) = ramp(NwTEST, 1._DP, 40._DP)
  wTEST(:) = ramp(NwTEST, 1.E-2_DP, 3.E4_DP, .true.) ! try larger NwTEST
  nuTEST = MKS%clight/MKS%micron / wTEST ! as reciprocal of wTEST, nuTEST is not in log
  FnuTEST = gaussLine(wTEST, 15.555_DP, 15.555_DP*0.0055722841_DP, .TRUE.) ! /Hz
  ! FnuTEST = lorentzBand(wTEST, 6.2_DP, 0.060000000_DP, 0.031300317_DP, .TRUE.) ! /Hz
  PRINT*, 'I(test no reverse) = ', integ_tabulated(nuTEST, FnuTEST)!, xlog=.true.)
  PRINT*, 'I(test) = ', integ_tabulated(reverse(nuTEST), reverse(FnuTEST))!, xlog=.true.)

  PRINT*, 'Integration test [Done]'
  PRINT*, '------------------------------'
  !!----------------------------------------------------------------

  PRINT*, "Number of components (dust continuum): ", Nbb
  ! PRINT*, "Number of lines: ", Nline
  ! PRINT*, "Number of bands: ", Nband
  
  ALLOCATE(label(Nbb), nQabs(Nbb), massBB(Nbb), tempBB(Nbb))
  
  label = (/'Sil_D03 ', 'ACAR_Z96'/)
  massBB = 3.E-2_DP ! Msun/pc2
  tempBB = 200._DP
  Iline = 1.E-7 ! W/m2
  ref1 = 15.555_DP ! [NeIII] 1 line
  sigma = 0.0055722841_DP * ref1
  ref2 = 6.2_DP ! 6.2 micron band
  Iband = 1.E-7 ! W/m2
  sigmaL = 0.060000000_DP * ref2
  sigmaS = 0.031300317_DP * ref2

  !! GET PARAM
  !!-----------
  CALL make_Qabs(label, NQABS=nQabs)

  !! Calculate Fnu
  !!---------------
  wave = ramp(Nw, 1._DP, 40._DP)
  nu = MKS%clight/MKS%micron / wave
  !! extinction
  Pabs = EXP(-.1_DP/1.086_DP * extCurve(wave))
  !! mbb
  FnuBB = interp_lin_sorted(modifBB(tempBB(1), nQabs(1)), &
          nQabs(1)%nu, nu, xlog=.True., ylog=.True.) * massBB(1) * MKS%Msun/MKS%pc**2
  FnuBB = FnuBB + interp_lin_sorted(modifBB(tempBB(2), nQabs(2)), &
          nQabs(2)%nu, nu, xlog=.True., ylog=.True.) * massBB(2) * MKS%Msun/MKS%pc**2
  !! gauss line
  FnuLINE = gaussLine(wave, ref1, sigma, .TRUE.)
  PRINT*, 'I(gauss) = ', integ_tabulated(reverse(nu), reverse(FnuLINE))
  FnuLINE = FnuLINE * Iline /MKS%Jy/1.E6 ! divided by MJy to convert W/m2/Hz to MJy
  !! lorentz band
  FnuBAND = lorentzBand(wave, ref2, sigmaS, sigmaL, .TRUE.)
  PRINT*, 'I(lorentz) = ', integ_tabulated(reverse(nu), reverse(FnuBAND))
  FnuBAND = FnuBAND * Iband /MKS%Jy/1.E6
  !! total
  Fnu = (FnuBB + FnuBAND) * Pabs + FnuLINE
  PRINT*, 'MAX(BB) = ', MAXVAL(FnuBB)
  PRINT*, 'MAX(gauss) = ', MAXVAL(FnuLINE)
  PRINT*, 'MAX(lorentz) = ', MAXVAL(FnuBAND)

  !! OUTPUTS
  !!---------
  path = 'out/test_profiles'
  CALL write_hdf5(wave, NAME="Wavelength (microns)", &
                  COMPRESS=.False., APPEND=.False., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5(nu, NAME="Wavelength (Hz)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5(FnuBB, NAME="FnuBB (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5(FnuLINE, NAME="FnuLINE (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5(FnuBAND, NAME="FnuBAND (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5(Fnu, NAME="Fnu (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5(FnuTEST, NAME="FnuTEST (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5(wTEST, NAME="x", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)

  !! Free memory space
  DEALLOCATE(wave, nu, Pabs, label, massBB, tempBB, &
             FnuLINE, FnuBAND, FnuBB, Fnu, FnuTEST)

END PROGRAM test_profiles
