MODULE gen_external

  USE utilities, ONLY: DP
  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PUBLIC                             :: NwOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wOBS, nuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: FnuOBS, dFnuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: FnuMOD, resid
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Pabs
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Fnu_cont
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Fnu_line
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Fnu_band
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Fnu_star

  !!------------
  !! Init param
  !!------------
  INTEGER, PARAMETER, PUBLIC :: Nbb=3, Nline=3, Nband=2, NAv=1, Nstar=1
  !! Npar = 2*Nbb + 3*Nline + 4*Nband + NAv + Nstar
  INTEGER, PARAMETER, PUBLIC :: Npar=25
  CHARACTER(*), PARAMETER, DIMENSION(Nbb), PUBLIC :: &
    compdust = ['ACH2_Z96             ', &
                'BE_Z96               ', &
                'Sil_D03              ']
  
  PUBLIC :: funcResid

CONTAINS

  !! Main function for L-M
  !!-----------------------
  FUNCTION funcResid(parr, NwOBS0)
    
    USE auxil, ONLY: Qabs_str, make_Qabs, specModel
    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN)                        :: NwOBS0
    REAL(DP), DIMENSION(:), INTENT(IN)         :: parr
    TYPE(Qabs_str), DIMENSION(Nbb)             :: nQabs
    REAL(DP), DIMENSION(NwOBS0)                :: funcResid

    CALL make_Qabs(compdust, nQabs)
    FnuMOD = specModel(wOBS, parr, nQabs, Nbb, Nline, Nband, NAv, Nstar)
    
    funcResid = FnuOBS(:) - FnuMOD(:)

  END FUNCTION funcResid
  
END MODULE gen_external

!!==========================================================================

PROGRAM test_general

  USE auxil, ONLY: par_str, Qabs_str, make_Qabs, degradeRes, specModel
  USE datable, ONLY: res, LIN, BIN
  USE arrays, ONLY: ramp
  USE random, ONLY: rand_norm
  USE utilities, ONLY: DP, PRING, TRIMLR
  USE constants, ONLY: MKS
  USE inout, ONLY: write_hdf5, h5ext, write_ascii
  USE chi2_minimization, ONLY: chi2min_LM
  USE gen_external, ONLY: funcResid, NwOBS, wOBS, nuOBS, &
                          FnuOBS, dFnuOBS, FnuMOD, resid, &
                          Pabs, Fnu_cont, Fnu_line, Fnu_band, Fnu_star, &
                          Nbb, Nline, Nband, NAv, Nstar, Npar, compdust
  IMPLICIT NONE

  CHARACTER(30)                             :: filOUT
  INTEGER                                   :: i
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: tau_tab, Fnu_cont_tab
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: Fnu_line_tab, Fnu_band_tab
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: Fnu

  INTEGER, PARAMETER                        :: Nt=1000
  REAL(DP), DIMENSION(Npar)                 :: parr, parr0
  TYPE(par_str)                             :: par
  TYPE(Qabs_str), DIMENSION(Nbb)            :: nQabs
  REAL(DP)                                  :: a

  !! chi2min_LM parinfo
  !!--------------------
  LOGICAL, DIMENSION(Npar)                  :: fixed
  LOGICAL, DIMENSION(Npar,2)                :: limited
  CHARACTER(20), DIMENSION(Npar)            :: parname
  INTEGER                                   :: status, Nparfree
  INTEGER, DIMENSION(Npar)                  :: itied
  REAL(DP), PARAMETER                       :: tol = 1.E-10_DP
  REAL(DP)                                  :: chi2red
  REAL(DP), DIMENSION(Npar)                 :: parerr
  REAL(DP), DIMENSION(Npar,2)               :: limits
  REAL(DP), DIMENSION(Npar,Npar)            :: covar

  wOBS = ramp(Nt, 1._DP, 40._DP, xlog=.TRUE.)
  NwOBS = SIZE(wOBS)
  nuOBS = MKS%clight/MKS%micron / wOBS

  ALLOCATE(Fnu(NwOBS), FnuOBS(NwOBS), dFnuOBS(NwOBS), FnuMOD(NwOBS))
  ALLOCATE(Pabs(NwOBS), Fnu_cont(NwOBS), Fnu_line(NwOBS), &
           Fnu_band(NwOBS), Fnu_star(NwOBS))
  ALLOCATE(tau_tab(NAv, NwOBS), Fnu_cont_tab(Nbb,NwOBS), &
           Fnu_line_tab(Nline,NwOBS), Fnu_band_tab(Nband,NwOBS))
  
  !! parr = [massBB [Msun/pc2], tempBB [T], &
  parr = [2.E-2_DP, 100._DP, 5.E-2_DP, 100._DP, 3.E-2_DP, 100._DP, & 
  !!      Iline, Cline, Wline, &
          2.E-9_DP, LIN(5)%wave, degradeRes(LIN(5)%wave, .01_DP, 'SL-LL'), & 
          2.E-9_DP, LIN(10)%wave, degradeRes(LIN(10)%wave, .01_DP, 'SL-LL'), &
          1.5E-9_DP, LIN(25)%wave, degradeRes(LIN(25)%wave, .01_DP, 'SL-LL'), &
  !!      Iband, Cband, WSband, WLband, &
          9.E-9_DP, BIN(10)%wave, BIN(10)%sigmaS, BIN(10)%sigmaL, & 
          7.E-9_DP, BIN(15)%wave, BIN(15)%sigmaS, BIN(15)%sigmaL, &
  !!      Av [mag], &
          0._DP, & 
  !!      Fstar [Lsun/pc2]]
          1.E4_DP]

  CALL make_Qabs(compdust, nQabs)
  print*, 'coucou1'
  FnuOBS(:) = specModel(wOBS, parr, nQabs, Nbb, Nline, Nband, NAv, Nstar)!, &
                     ! tau_tab, Fnu_cont_tab, Fnu_line_tab, Fnu_band_tab, &
                     ! Pabs, Fnu_cont, Fnu_line, Fnu_band, Fnu_star)
  dFnuOBS = rand_norm(NwOBS) * 0.01_DP*MAXVAL(FnuOBS)
  print*, 'coucou2'
  FnuOBS = FnuOBS + dFnuOBS
  ! print*, 'MAX(cont) = ', MAXVAL(Fnu_cont*Pabs)
  ! print*, 'MAX(line) = ', MAXVAL(Fnu_line)
  ! print*, 'MAX(band) = ', MAXVAL(Fnu_band*Pabs)
  ! print*, 'MAX(star) = ', MAXVAL(Fnu_star*Pabs)
  ! print*, 'MIN(Pabs) = ', MINVAL(Pabs)

  filOUT = 'out/test_total'
  CALL write_hdf5(wOBS, NAME="Wavelength (micron)", &
                  COMPRESS=.False., APPEND=.False., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(dFnuOBS, NAME="FnuUNC (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(FnuOBS, NAME="FnuOBS (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  !! Write components
  ! DO i=1,Nbb
  !   CALL write_hdf5(Fnu_cont_tab(i,:)*Pabs(:), &
  !                   NAME='FnuBB_'//TRIMLR(PRING(i)), &
  !                   COMPRESS=.False., APPEND=.True., &
  !                   FILE=TRIMLR(filOUT)//h5ext)

  ! END DO
  print*, 'Gen synthetic spectrum [done]'

  print*, '=================================================='
  
  !! Parameter constraints
  !! parr = [massBB [Msun/pc2], tempBB [T], &
  fixed(:) = [0,1, 0,1, 0,1, & 
  !!          Iline, Cline, Wline, &
              0,1,1, 0,1,1, 0,1,1, &
  !!          Iband, Cband, WSband, WLband, &
              0,1,1,1, 0,1,1,1, &
  !!          Av [mag], &
              1, & 
  !!          Fstar [Lsun/pc2]]
              0]
  ! fixed(:) = 0
  ! fixed(24) = 1 ! Av
  itied(:) = 0
  limited(:,:) = 0
  limits(:,:) = 0._DP
  Nparfree = COUNT((.NOT. fixed(:)) .AND. (itied(:) <= 0))
  parname = ['massBB1', 'tempBB1', 'massBB2', 'tempBB2', 'massBB3', 'tempBB3', &
             'Iline1 ', 'Cline1 ', 'Wline1 ', &
             'Iline2 ', 'Cline2 ', 'Wline2 ', &
             'Iline3 ', 'Cline3 ', 'Wline3 ', &
             'Iband1 ', 'Cband1 ', 'WSband1', 'WLband1', &
             'Iband2 ', 'Cband2 ', 'WSband2', 'WLband2', &
             'Av     ', &
             'Fstar  ']

  !! parr0 = [massBB [Msun/pc2], tempBB [T], &
  parr0 = [1._DP, 100._DP, 1._DP, 100._DP, 1._DP, 100._DP, & 
  !!       Iline, Cline, Wline, &
           1._DP, LIN(5)%wave, degradeRes(LIN(5)%wave, .01_DP, 'SL-LL'), & 
           1._DP, LIN(10)%wave, degradeRes(LIN(10)%wave, .01_DP, 'SL-LL'), &
           1._DP, LIN(25)%wave, degradeRes(LIN(25)%wave, .01_DP, 'SL-LL'), &
  !!       Iband, Cband, WSband, WLband, &
           1._DP, BIN(10)%wave, BIN(10)%sigmaS, BIN(10)%sigmaL, & 
           1._DP, BIN(15)%wave, BIN(15)%sigmaS, BIN(15)%sigmaL, &
  !!       Av [mag], &
           0._DP, & 
  !!       Fstar [Lsun/pc2]]
           1._DP]

  !! Run the fitter
  CALL chi2min_LM (funcResid, NwOBS, parr0, tol, resid, status, VERBOSE=.True., &
                   LIMITED=limited, LIMITS=limits, FIXED=fixed, ITIED=itied, &
                   PARNAME=parname, CHI2RED=chi2red, PARERR=parerr, COVAR=covar)

  !! Compute model result
  PRINT*
  PRINT*, "Checking output:"
  PRINT*, "  status = "//TRIMLR(PRING(status))
  ! chi2red = SUM(deviate(:)**2) / (Nobs-Nparfree)
  PRINT*, "  chi2 = "//TRIMLR(PRING(chi2red,NDEC=10))
  CALL write_ascii(VEC1=parr0,VEC2=parerr,FILE="out/chi2min.txt")
  PRINT*
  PRINT*, "Covariance matrix:"
  DO i=1,Npar
    PRINT*, REAL(covar(:,i), KIND(0.))
  END DO
  PRINT*

  FnuMOD = specModel(wOBS, parr0, nQabs, Nbb, Nline, Nband, NAv, Nstar, &
                     tau_tab, Fnu_cont_tab, Fnu_line_tab, Fnu_band_tab, &
                     PABS=Pabs, FNU_CONT=Fnu_cont, &
                     FNU_LINE=Fnu_line, FNU_BAND=Fnu_band, &
                     FNU_STAR=Fnu_star)
  CALL write_hdf5(Fnu_cont*Pabs, NAME="FnuBB (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(Fnu_line+(Fnu_cont+Fnu_star)*Pabs, NAME="FnuLINE (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5((Fnu_band+Fnu_cont+Fnu_star)*Pabs, NAME="FnuBAND (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(Fnu_star*Pabs, NAME="FnuSTAR (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(FnuMOD, NAME="FnuMOD (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  
  print*, 'Chi2 fit synt spec [done]'

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
  a = degradeRes(6._DP, .001_DP, 'SL-LL')
  print*, 'a = ', a
  print*, 'AKARI_Ns line width = ', res%dw_w_AKARI_Ns
  print*, 'Test degradeRes [done]'

END PROGRAM test_general
