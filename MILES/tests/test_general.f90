MODULE gen_external

  USE utilities, ONLY: DP
  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PUBLIC                             :: Npar, NwOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wOBS, nuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: FnuOBS, dFnuOBS, FnuMOD
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Pabs
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Fnu_cont
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Fnu_line
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Fnu_band
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Fnu_star
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: specModel
    
  PUBLIC :: funcResid

CONTAINS

  !! Main function for L-M
  !!-----------------------
  FUNCTION funcResid(par, NwOBS0)
    
    USE auxil, ONLY: Qabs_str, specModel
    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN)                      :: NwOBS0
    REAL(DP), DIMENSION(:), INTENT(IN)       :: par
    TYPE(Qabs_str), DIMENSION(:), INTENT(IN) :: nQabs
    REAL(DP), DIMENSION(NwOBS0)              :: funcResid

    FnuMOD(:) = specModel(wOBS(:), par, nQabs, NAv, Nbb, Nline, Nband, Nstar, &
                          tau_tab,  Fnu_cont_tab, Fnu_line_tab, Fnu_band_tab, &
                          Pabs, Fnu_cont, Fnu_line, Fnu_band, Fnu_star))
    
    funcResid = FnuOBS(:) - FnuMOD(:)

  END FUNCTION funcResid
  
END MODULE gen_external

!!==========================================================================

PROGRAM test_general

  USE auxil, ONLY: 
  USE datable
  USE arrays, ONLY: ramp
  USE random, ONLY: rand_norm
  USE utilities, ONLY: DP, PRING, TRIMLR
  USE constants, ONLY: MKS
  USE inout, ONLY: write_hdf5, h5ext
  USE chi2_minimization, ONLY: chi2min_LM
  USE gen_external, ONLY: funcResid, Npar, NwOBS, wOBS, nuOBS, &
                          FnuOBS, dFnuOBS, FnuMOD
  IMPLICIT NONE
  
  LOGICAL, DIMENSION(Npar)                  :: fixed
  LOGICAL, DIMENSION(Npar,2)                :: limited
  CHARACTER(20), DIMENSION(Npar)            :: parname
  CHARACTER(30)                             :: instr, path
  INTEGER                                   :: i, status, Nparfree
  INTEGER, DIMENSION(Npar)                  :: itied
  INTEGER, PARAMETER                        :: Nt=1000
  REAL(DP)                                  :: chi2red
  REAL(DP), DIMENSION(Npar)                 :: parerr
  REAL(DP), DIMENSION(Npar,Npar)            :: covar
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: wave, nu
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: tau_tab, Fnu_cont_tab
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: Fnu_line_tab, Fnu_band_tab
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: Pabs, Fnu_cont, Fnu_star
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: Fnu_line, Fnu_band
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: Fnu

  INTEGER                                   :: Nbb, Nline, Nband, NAv, Nstar
  CHARACTER(30), DIMENSION(:), ALLOCATABLE  :: label
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: parr
  TYPE(Qabs_STR), DIMENSION(:), ALLOCATABLE :: nQabs

  wave = ramp(Nt, 1._DP, 40._DP)
  nu = MKS%clight/MKS%micron / wave

  !! spec_MODEL
  !!------------
  CALL chi2_INIT(label, parr, NAv, Nbb, Nline, Nband, Nstar)
  print*, 'chi2_INIT [done]'
  ALLOCATE(nQabs(Nbb))
  CALL make_Qabs(label, nQabs)
  print*, 'make_Qabs [done]'
  
  ALLOCATE(Pabs(Nt), Fnu_cont(Nt), Fnu_line(Nt), &
           Fnu_band(Nt), Fnu_star(Nt))
  ALLOCATE(tau_tab(NAv, Nt), Fnu_cont_tab(Nbb,Nt), &
           Fnu_line_tab(Nline,Nt), Fnu_band_tab(Nband,Nt))
  ! ALLOCATE(Fnu(Nt))
  Fnu = specModel(wave, parr, nQabs, &
                  NAv, Nbb, Nline, Nband, Nstar, &
                  tau_tab, Fnu_cont_tab, Fnu_line_tab, Fnu_band_tab, &
                  PABS=Pabs, FNU_CONT=Fnu_cont, &
                  FNU_LINE=Fnu_line, FNU_BAND=Fnu_band, &
                  FNU_STAR=Fnu_star)
  print*, 'MAX(cont) = ', MAXVAL(Fnu_cont*Pabs)
  print*, 'MAX(line) = ', MAXVAL(Fnu_line)
  print*, 'MAX(band) = ', MAXVAL(Fnu_band*Pabs)
  print*, 'MAX(star) = ', MAXVAL(Fnu_star*Pabs)
  print*, 'MIN(Pabs) = ', MINVAL(Pabs)
  print*, 'specModel [done]'

  !! OUTPUTS
  !!---------
  path = 'out/test_total'
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
  print*, 'generate synthetic spectrum [done]'

  print*, '=================================================='

  
  print*, 'chi2 fit of syn spec [done]'

  print*, '=================================================='

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
  
  !! res
  !!-----
  instr = 'SL-LL'
  a = degradeRes(6._DP, .001_DP, instr)
  print*, 'a = ', a
  print*, 'AKARI_Ns line width = ', res%dw_w_AKARI_Ns
  print*, 'test degradeRes [done]'

END PROGRAM test_general
