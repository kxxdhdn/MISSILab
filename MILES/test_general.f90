PROGRAM test_general

  USE auxil
  USE factable
  USE arrays, ONLY: ramp
  USE utilities, ONLY: DP, PRING, TRIMLR
  USE constants, ONLY: MKS
  USE inout, ONLY: write_hdf5, h5ext
  IMPLICIT NONE
  
  ! INTEGER i
  ! INTEGER Nparm
  CHARACTER(30) :: ins, path
  REAL(DP)      :: a
  INTEGER       :: i, Nt=1000
  REAL(DP), DIMENSION(1000) :: Ftest, ext
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wave, nu
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: tau_tab, Fnu_cont_tab
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Fnu_line_tab, Fnu_band_tab
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Pabs, Fnu_cont, Fnu_star
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Fnu_line, Fnu_band
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Fnu
  ! REAL(DP), DIMENSION(1000) :: Fnu

  INTEGER :: Nbb, Nline, Nband, NAv, Nstar
  CHARACTER(30), DIMENSION(:), ALLOCATABLE :: label
  REAL(DP), DIMENSION(:), ALLOCATABLE :: parm
  TYPE(Qabs_STR), DIMENSION(:), ALLOCATABLE :: nQ_str

  wave = ramp(Nt, 1._DP, 40._DP)
  nu = MKS%clight/MKS%micron / wave
  print*, 'nu_GEN done'

  !! extCurve
  !!----------
  ext = extCurve(nu)
  print*, 'extCurve test (should be 1000): ', SHAPE(ext)
  print*, ext

  !! spec_MODEL
  !!------------
  ! Ftest = gaussLine(wave, 0._DP, 1._DP)
  Ftest = lorentzBand(wave, 5._DP, .5_DP, 2._DP)
  ! print*, Ftest
  print*, 'id func test done'

  CALL chi2_INIT(label, parm, NAv, Nbb, Nline, Nband, Nstar)
  print*, 'chi2_INIT done'
  ALLOCATE(nQ_str(Nbb))
  CALL optics_INIT(label, nQ_str)
  print*, 'optics_INIT done'
  
  ALLOCATE(Pabs(Nt), Fnu_cont(Nt), Fnu_line(Nt), &
           Fnu_band(Nt), Fnu_star(Nt))
  ALLOCATE(tau_tab(NAv, Nt), Fnu_cont_tab(Nbb,Nt), &
           Fnu_line_tab(Nline,Nt), Fnu_band_tab(Nband,Nt))
  ! ALLOCATE(Fnu(Nt))
  Fnu = specModel(wave, parm, nQ_str, &
                  NAv, Nbb, Nline, Nband, Nstar, &
                  tau_tab, Fnu_cont_tab, Fnu_line_tab, Fnu_band_tab, &
                  PABS=Pabs, FNU_CONT=Fnu_cont, &
                  FNU_LINE=Fnu_line, FNU_BAND=Fnu_band, &
                  FNU_STAR=Fnu_star)
  ! print*, Fnu_cont
  ! print*, shape(Fnu_cont_tab)
  ! Fnu = specModel(wave, parm, nQ_str, &
  !                 NAv, Nbb, Nline, Nband, Nstar)
  print*, 'specModel done'

  ! !! subroutine option (in obsolete script)
  ! CALL spec_MODEL(label, parm, &
  !                 NAv, Nbb, Nline, Nband, Nstar, &
  !                 Pabs, Fnu_cont, Fnu_line, Fnu_band, &
  !                 Fnu_star, Fnu, wave)
  path = 'output/test_total'
  CALL write_hdf5(wave, NAME="Wavelength (micron)", &
                  COMPRESS=.False., APPEND=.False., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5(Fnu_cont*Pabs, NAME="FnuBB (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5(Fnu_line+(Fnu_cont+Fnu_star)*Pabs, NAME="FnuLINE (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5((Fnu_band+Fnu_cont+Fnu_star)*Pabs, NAME="FnuBAND (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5(Fnu_star*Pabs, NAME="FnuSTAR (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)
  CALL write_hdf5(Fnu, NAME="Fnu (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(path)//h5ext)
  !! Write components
  DO i=1,Nbb
    CALL write_hdf5(Fnu_cont_tab(i,:)*Pabs(:), &
                    NAME='FnuBB_'//TRIMLR(PRING(i)), &
                    COMPRESS=.False., APPEND=.True., &
                    FILE=TRIMLR(path)//h5ext)

  END DO
  print*, 'write_RESULTS done'

  print*, '----------------------------------'

  !! BIN
  !!-----
  ! DO i=1,33
  !   print*, BIN(i)%sigmaS
  ! END DO

  !! LIN
  !!-----
  ! DO i=1,46
  !   print*, LIN(i)%name
  ! END DO
  
  !! RES
  !!-----
  ins = 'SL-LL'
  a = degradeRes(6._DP, .001_DP, ins)
  print*, 'a = ', a
  print*, 'AKARI_Ns RES = ', RES%AKARI_Ns

END PROGRAM test_general
