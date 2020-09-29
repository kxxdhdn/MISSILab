MODULE fitChi2_external

  USE auxil, ONLY: parinfo_type, Qabs_type
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

  !! Init param
  !!------------
  INTEGER, SAVE, PUBLIC :: Nbb, Nline, Nband, NAv, Nstar, Npar
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parinfo
  CHARACTER(30), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC      :: compdust
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC    :: nQabs
  
  PUBLIC :: funcResid

CONTAINS

  !! Main function for L-M
  !!-----------------------
  FUNCTION funcResid(pargen, NwOBS0)
    
    USE auxil, ONLY: specModel
    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN)                :: NwOBS0
    REAL(DP), DIMENSION(:), INTENT(IN) :: pargen
    REAL(DP), DIMENSION(NwOBS0)        :: funcResid

    FnuMOD = specModel(wOBS, pargen, nQabs, Nbb, Nline, Nband, NAv, Nstar)
    
    funcResid = FnuOBS(:) - FnuMOD(:)

  END FUNCTION funcResid

END MODULE fitChi2_external

!!==========================================================================

PROGRAM test_fitChi2

  USE auxil, ONLY: make_par, make_Qabs, specModel, degradeRes
  USE datable, ONLY: res, LIN, BIN
  USE arrays, ONLY: ramp
  USE random, ONLY: rand_norm
  USE utilities, ONLY: DP, PRING, TRIMLR
  USE constants, ONLY: MKS
  USE inout, ONLY: write_hdf5, h5ext, write_ascii
  USE chi2_minimization, ONLY: chi2min_LM
  USE fitChi2_external, ONLY: funcResid, NwOBS, wOBS, nuOBS, &
                          FnuOBS, dFnuOBS, FnuMOD, resid, &
                          Pabs, Fnu_cont, Fnu_line, Fnu_band, Fnu_star, &
                          Nbb, Nline, Nband, NAv, Nstar, Npar, parinfo, compdust, nQabs
  IMPLICIT NONE

  CHARACTER(30)                            :: filOUT
  INTEGER                                  :: i

  INTEGER, PARAMETER                       :: Ngrid=1000
  REAL(DP), DIMENSION(:), ALLOCATABLE      :: pargen, parr0
  REAL(DP)                                 :: a

  !! chi2min_LM parinfo
  !!--------------------
  ! LOGICAL, DIMENSION(:), ALLOCATABLE       :: fixed
  LOGICAL, DIMENSION(:,:), ALLOCATABLE     :: limited
  ! CHARACTER(20), DIMENSION(:), ALLOCATABLE :: parname
  INTEGER                                  :: status, Nparfree
  ! INTEGER, DIMENSION(:), ALLOCATABLE       :: itied
  REAL(DP), PARAMETER                      :: tol = 1.E-10_DP
  REAL(DP)                                 :: chi2red
  REAL(DP), DIMENSION(:), ALLOCATABLE      :: parerr
  REAL(DP), DIMENSION(:,:), ALLOCATABLE    :: limits
  REAL(DP), DIMENSION(:,:), ALLOCATABLE    :: covar

  !!----------------
  !! Generate param
  !!----------------
  Nbb=3
  Nline=12
  Nband=14
  NAv=1
  Nstar=1

  !! pargen = [massBB [Msun/pc2], tempBB [K], &
  pargen = [2.E-2_DP, 100._DP, 5.E-2_DP, 100._DP, 3.E-2_DP, 100._DP, & 
  !!      Iline, Cline, Wline, &
          1.E-9_DP, LIN(3)%wave, degradeRes(LIN(3)%wave, .01_DP, 'SL-LL'), & 
          2.5E-9_DP, LIN(9)%wave, degradeRes(LIN(9)%wave, .01_DP, 'SL-LL'), & 
          3.E-9_DP, LIN(11)%wave, degradeRes(LIN(11)%wave, .01_DP, 'SL-LL'), &
          1.5E-9_DP, LIN(20)%wave, degradeRes(LIN(20)%wave, .01_DP, 'SL-LL'), &
          1.5E-9_DP, LIN(23)%wave, degradeRes(LIN(23)%wave, .01_DP, 'SL-LL'), &
          2.5E-9_DP, LIN(24)%wave, degradeRes(LIN(24)%wave, .01_DP, 'SL-LL'), &
          1.5E-9_DP, LIN(26)%wave, degradeRes(LIN(26)%wave, .01_DP, 'SL-LL'), &
          3.E-9_DP, LIN(27)%wave, degradeRes(LIN(27)%wave, .01_DP, 'SL-LL'), &
          3.5E-9_DP, LIN(29)%wave, degradeRes(LIN(29)%wave, .01_DP, 'SL-LL'), &
          4.5E-9_DP, LIN(33)%wave, degradeRes(LIN(33)%wave, .01_DP, 'SL-LL'), &
          1.5E-9_DP, LIN(34)%wave, degradeRes(LIN(34)%wave, .01_DP, 'SL-LL'), &
          1.E-9_DP, LIN(36)%wave, degradeRes(LIN(36)%wave, .01_DP, 'SL-LL'), &
  !!      Iband, Cband, WSband, WLband, &
          1.E-9_DP, BIN(1)%wave, BIN(1)%sigmaS, BIN(1)%sigmaL, & 
          .8E-9_DP, BIN(2)%wave, BIN(2)%sigmaS, BIN(2)%sigmaL, & 
          .5E-9_DP, BIN(7)%wave, BIN(7)%sigmaS, BIN(7)%sigmaL, & 
          .7E-9_DP, BIN(8)%wave, BIN(8)%sigmaS, BIN(8)%sigmaL, & 
          .3E-9_DP, BIN(12)%wave, BIN(12)%sigmaS, BIN(12)%sigmaL, & 
          .5E-9_DP, BIN(13)%wave, BIN(13)%sigmaS, BIN(13)%sigmaL, & 
          .4E-9_DP, BIN(14)%wave, BIN(14)%sigmaS, BIN(14)%sigmaL, & 
          1.E-9_DP, BIN(16)%wave, BIN(16)%sigmaS, BIN(16)%sigmaL, & 
          .3E-9_DP, BIN(19)%wave, BIN(19)%sigmaS, BIN(19)%sigmaL, & 
          .7E-9_DP, BIN(20)%wave, BIN(20)%sigmaS, BIN(20)%sigmaL, &
          .6E-9_DP, BIN(21)%wave, BIN(21)%sigmaS, BIN(21)%sigmaL, & 
          .5E-9_DP, BIN(24)%wave, BIN(24)%sigmaS, BIN(24)%sigmaL, & 
          .5E-9_DP, BIN(25)%wave, BIN(25)%sigmaS, BIN(25)%sigmaL, & 
          .9E-9_DP, BIN(30)%wave, BIN(30)%sigmaS, BIN(30)%sigmaL, & 
  !!      Av [mag], &
          0._DP, & 
  !!      Fstar [Lsun/pc2]]
          1.E4_DP]
  
  !!-------------
  !! Create grid
  !!-------------
  wOBS = ramp(Ngrid, 1._DP, 40._DP, xlog=.TRUE.)
  NwOBS = SIZE(wOBS)
  nuOBS = MKS%clight/MKS%micron / wOBS

  !! Init optical properties
  ALLOCATE(compdust(Nbb), nQabs(Nbb))
  
  compdust = ['ACH2_Z96             ', &
              'BE_Z96               ', &
              'Sil_D03              ']
  CALL make_Qabs(compdust, nQabs, WAVEALL=wOBS)

  !! Build synthetic spectrum
  ALLOCATE(FnuOBS(NwOBS), dFnuOBS(NwOBS))
  ! ALLOCATE(Pabs(NwOBS),Fnu_cont(NwOBS),Fnu_line(NwOBS),Fnu_band(NwOBS),Fnu_star(NwOBS))

  FnuOBS(:) = specModel(wOBS, pargen, nQabs, Nbb, Nline, Nband, NAv, Nstar)!, &
                     ! Pabs, Fnu_cont, Fnu_line, Fnu_band, Fnu_star)
  dFnuOBS = rand_norm(NwOBS) * 0.01_DP*MAXVAL(FnuOBS)
  FnuOBS = FnuOBS + dFnuOBS
  ! print*, 'MAX(cont) = ', MAXVAL(Fnu_cont*Pabs)
  ! print*, 'MAX(line) = ', MAXVAL(Fnu_line)
  ! print*, 'MAX(band) = ', MAXVAL(Fnu_band*Pabs)
  ! print*, 'MAX(star) = ', MAXVAL(Fnu_star*Pabs)
  ! print*, 'MIN(Pabs) = ', MINVAL(Pabs)

  filOUT = 'out/test_fitChi2syn'
  CALL write_hdf5(wOBS, NAME="Wavelength (microns)", &
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
  
  !! parr0 = [massBB [Msun/pc2], tempBB [K], &
  parr0 = [1._DP, 50._DP, 1._DP, 50._DP, 1._DP, 50._DP, & 
  !!      Iline, Cline, Wline, &
          1._DP, LIN(3)%wave, degradeRes(LIN(3)%wave, .01_DP, 'SL-LL'), & 
          1._DP, LIN(9)%wave, degradeRes(LIN(9)%wave, .01_DP, 'SL-LL'), & 
          1._DP, LIN(11)%wave, degradeRes(LIN(11)%wave, .01_DP, 'SL-LL'), &
          1._DP, LIN(20)%wave, degradeRes(LIN(20)%wave, .01_DP, 'SL-LL'), &
          1._DP, LIN(23)%wave, degradeRes(LIN(23)%wave, .01_DP, 'SL-LL'), &
          1._DP, LIN(24)%wave, degradeRes(LIN(24)%wave, .01_DP, 'SL-LL'), &
          1._DP, LIN(26)%wave, degradeRes(LIN(26)%wave, .01_DP, 'SL-LL'), &
          1._DP, LIN(27)%wave, degradeRes(LIN(27)%wave, .01_DP, 'SL-LL'), &
          1._DP, LIN(29)%wave, degradeRes(LIN(29)%wave, .01_DP, 'SL-LL'), &
          1._DP, LIN(33)%wave, degradeRes(LIN(33)%wave, .01_DP, 'SL-LL'), &
          1._DP, LIN(34)%wave, degradeRes(LIN(34)%wave, .01_DP, 'SL-LL'), &
          1._DP, LIN(36)%wave, degradeRes(LIN(36)%wave, .01_DP, 'SL-LL'), &
  !!      Iband, Cband, WSband, WLband, &
          1._DP, BIN(1)%wave, BIN(1)%sigmaS, BIN(1)%sigmaL, & 
          1._DP, BIN(2)%wave, BIN(2)%sigmaS, BIN(2)%sigmaL, & 
          1._DP, BIN(7)%wave, BIN(7)%sigmaS, BIN(7)%sigmaL, & 
          1._DP, BIN(8)%wave, BIN(8)%sigmaS, BIN(8)%sigmaL, & 
          1._DP, BIN(12)%wave, BIN(12)%sigmaS, BIN(12)%sigmaL, & 
          1._DP, BIN(13)%wave, BIN(13)%sigmaS, BIN(13)%sigmaL, & 
          1._DP, BIN(14)%wave, BIN(14)%sigmaS, BIN(14)%sigmaL, & 
          1._DP, BIN(16)%wave, BIN(16)%sigmaS, BIN(16)%sigmaL, & 
          1._DP, BIN(19)%wave, BIN(19)%sigmaS, BIN(19)%sigmaL, & 
          1._DP, BIN(20)%wave, BIN(20)%sigmaS, BIN(20)%sigmaL, &
          1._DP, BIN(21)%wave, BIN(21)%sigmaS, BIN(21)%sigmaL, & 
          1._DP, BIN(24)%wave, BIN(24)%sigmaS, BIN(24)%sigmaL, & 
          1._DP, BIN(25)%wave, BIN(25)%sigmaS, BIN(25)%sigmaL, & 
          1._DP, BIN(30)%wave, BIN(30)%sigmaS, BIN(30)%sigmaL, & 
  !!      Av [mag], &
          0._DP, & 
  !!      Fstar [Lsun/pc2]]
          1._DP]

  CALL make_par(parr0, Nbb, Nline, Nband, NAv, Nstar, NPAR=Npar, PARINFO=parinfo)
  Nparfree = COUNT((.NOT. parinfo(:)%fixed) .AND. (parinfo(:)%itied <= 0))
  
  ALLOCATE(limits(Npar,2), limited(Npar,2), parerr(Npar), covar(Npar,Npar))
  
  FORALL (i=1:Npar)
    limits(i,:) = parinfo(i)%limits(:)
    limited(i,:) = parinfo(i)%limited(:)
  END FORALL
  
  !! Run the fitter
  CALL chi2min_LM (funcResid, NwOBS, parr0, tol, resid, status, VERBOSE=.True., &
                   LIMITED=limited, LIMITS=limits, &
                   FIXED=parinfo(:)%fixed, ITIED=parinfo(:)%itied, &
                   PARNAME=parinfo(:)%name, CHI2RED=chi2red, PARERR=parerr, COVAR=covar)

  !! Compute model result
  PRINT*
  PRINT*, "Checking output:"
  PRINT*, "  status = "//TRIMLR(PRING(status))
  ! chi2red = SUM(deviate(:)**2) / (Nobs-Nparfree)
  PRINT*, "  chi2 = "//TRIMLR(PRING(chi2red,NDEC=10))
  CALL write_ascii(VEC1=parr0,VEC2=parerr,FILE="out/chi2min.txt")
  PRINT*
  PRINT*, "Covariance matrix:"
  ! DO i=1,Npar
  !   PRINT*, REAL(covar(:,i), KIND(0.))
  ! END DO
  PRINT*

  ALLOCATE(Pabs(NwOBS),Fnu_cont(NwOBS),Fnu_line(NwOBS),Fnu_band(NwOBS),Fnu_star(NwOBS))
  
  FnuMOD = specModel(wOBS, parr0, nQabs, Nbb, Nline, Nband, NAv, Nstar, &
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

  DEALLOCATE(compdust, nQabs, wOBS, nuOBS, FnuOBS, dFnuOBS)
  DEALLOCATE(pargen, parr0, limits, limited, parerr, covar)
  DEALLOCATE(Pabs,Fnu_cont,Fnu_line,Fnu_band,Fnu_star)

END PROGRAM test_fitChi2
