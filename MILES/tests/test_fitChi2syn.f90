MODULE fitChi2syn_external

  USE auxil, ONLY: parinfo_type, indpar_type, Qabs_type
  USE utilities, ONLY: DP
  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PUBLIC                             :: NwOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wOBS, nuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: FnuOBS, dFnuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: Fnu_mod, resid

  !! Init param
  !!------------
  TYPE(indpar_type), SAVE, PUBLIC                             :: ind
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parinfo
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC    :: Qabs
  
  PUBLIC :: residuals

CONTAINS

  !! Main function for L-M
  !!-----------------------
  FUNCTION residuals(par, NwOBS0)
    
    USE auxil, ONLY: specModel
    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN)                :: NwOBS0
    REAL(DP), DIMENSION(:), INTENT(IN) :: par
    REAL(DP), DIMENSION(NwOBS0)        :: residuals

    Fnu_mod(:) = specModel(wOBS(:), INDPAR=ind, PARVAL=par(:), QABS=Qabs(:))
    
    residuals(:) = (FnuOBS(:) - Fnu_mod(:)) / dFnuOBS(:)

  END FUNCTION residuals

END MODULE fitChi2syn_external


!!==========================================================================
!!                           Main Program
!!==========================================================================


PROGRAM test_fitChi2syn

  USE auxil, ONLY: read_master, degradeRes, set_indpar, specModel
  USE datable, ONLY: TABLine, TABand
  USE utilities, ONLY: DP, pring, trimLR, trimeq, timinfo, &
                       banner_program, ustd, initiate_clock, time_type
  USE arrays, ONLY: ramp, iwhere
  USE constants, ONLY: MKS
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, write_ascii, ascext, &
                   lenpar, lenpath
  USE grain_optics, ONLY: lendustQ
  USE random, ONLY: generate_newseed, rand_norm
  USE chi2_minimization, ONLY: chi2min_LM
  USE fitChi2syn_external, ONLY: NwOBS, wOBS, nuOBS, FnuOBS, dFnuOBS, resid, &
                                 residuals, Fnu_mod, ind, parinfo, Qabs
  IMPLICIT NONE

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  REAL(DP), PARAMETER :: tol = 1.E-10_DP
  CHARACTER(*), PARAMETER :: dirOUT = "./out/"
  CHARACTER(*), PARAMETER :: filOUT = dirOUT//'test_fitChi2syn'//h5ext
  CHARACTER(*), PARAMETER :: fiLOG = "./tmp/log_fitChi2syn"//ascext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Gen spec
  INTEGER, PARAMETER :: Nw=1000
  REAL(DP), DIMENSION(:), ALLOCATABLE :: pargen

  !! Input variables
  INTEGER :: i0, i, Npar, NiniMC, Niter
  INTEGER :: Ncont, Nband, Nline
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: parini
  ! CHARACTER(lenpar) :: spec_unit
  CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB, labL
  LOGICAL :: verbose, newseed!, calib, newinit, dostop
  !! chi2min_LM parameters
  ! LOGICAL, DIMENSION(:), ALLOCATABLE :: fixed
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: limited
  ! CHARACTER(20), DIMENSION(:), ALLOCATABLE :: parname
  INTEGER :: status, Nparfree
  INTEGER, DIMENSION(:), ALLOCATABLE :: itied
  REAL(DP) :: chi2red
  REAL(DP), DIMENSION(:), ALLOCATABLE :: par, parerr
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: limits
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: covar
  !! Output variables
  INTEGER, DIMENSION(2) :: unitlog
  TYPE(time_type) :: timestr

  !! Analysis
  REAL(DP), DIMENSION(:), ALLOCATABLE :: FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  
  !! Output settings
  !!-----------------
  CALL INITIATE_CLOCK(timestr)
  OPEN (ulog, FILE=fiLOG, STATUS="REPLACE", ACTION="WRITE")
  unitlog(:) = [ ulog, ustd ]

  !!-----------------
  !! Read the inputs
  !!-----------------
  labQ = (/'ACH2_Z96             ', &
           'BE_Z96               ', &
           'Sil_D03              '/)
  labB = (/'Main 3.3     ', & ! 1
             'Main 3.4     ', & ! 2
             'Main 6.2 (1) ', & ! 7
             'Main 6.2 (2) ', & ! 8
             'Plateau 7.7  ', & ! 12
             'Main 7.7 (1) ', & ! 13
             'Main 7.7 (2) ', & ! 14
             'Main 8.6     ', & ! 16
             'Small 11.0   ', & ! 19
             'Main 11.2    ', & ! 20
             'Plateau 11.3 ', & !21
             'Main 12.7 (1)', & ! 24
             'Main 12.7 (2)', & ! 25
             'Plateau 17.0 '/) ! 30
  labL = (/'H2S7  ', & ! 3
             'H2S5  ', & ! 9
             'ArII  ', & ! 11
             'ArIII1', & ! 20
             'H2S3  ', & ! 23
             'HeII2 ', & ! 24
             'SIV   ', & ! 26
             'H2S2  ', & ! 27
             'NeII  ', & ! 29
             'NeIII1', & ! 33
             'H2S1  ', & ! 34
             'SIII1 '/) ! 36

  !! No READ_MASTER buffer
  verbose = .TRUE.
  newseed = .TRUE.
  NiniMC = 0

  IF (newseed) CALL GENERATE_NEWSEED()
  
  !!---------------------
  !! Generate a spectrum
  !!---------------------

  !! pargen = [lnMovd2 LOG[Msun/pc2], lnT LOG[K], &
  pargen = [46._DP, 6._DP, 48._DP, 5._DP, 50._DP, 4._DP, & 
  !!        lnIline, Cline, Wline, &
            4._DP, TABLine(3)%wave, degradeRes(TABLine(3)%wave,.01_DP,'SL-LL'), & 
            3._DP, TABLine(9)%wave, degradeRes(TABLine(9)%wave,.01_DP,'SL-LL'), & 
            2._DP, TABLine(11)%wave, degradeRes(TABLine(11)%wave,.01_DP,'SL-LL'), &
            1._DP, TABLine(20)%wave, degradeRes(TABLine(20)%wave,.01_DP,'SL-LL'), &
            1._DP, TABLine(23)%wave, degradeRes(TABLine(23)%wave,.01_DP,'SL-LL'), &
            2._DP, TABLine(24)%wave, degradeRes(TABLine(24)%wave,.01_DP,'SL-LL'), &
            3._DP, TABLine(26)%wave, degradeRes(TABLine(26)%wave,.01_DP,'SL-LL'), &
            4._DP, TABLine(27)%wave, degradeRes(TABLine(27)%wave,.01_DP,'SL-LL'), &
            4._DP, TABLine(29)%wave, degradeRes(TABLine(29)%wave,.01_DP,'SL-LL'), &
            3._DP, TABLine(33)%wave, degradeRes(TABLine(33)%wave,.01_DP,'SL-LL'), &
            2._DP, TABLine(34)%wave, degradeRes(TABLine(34)%wave,.01_DP,'SL-LL'), &
            1._DP, TABLine(36)%wave, degradeRes(TABLine(36)%wave,.01_DP,'SL-LL'), &
  !!        lnIband, Cband, WSband, WLband, &
            1.9_DP, TABand(1)%wave, TABand(1)%sigmaS, TABand(1)%sigmaL, &
            2.8_DP, TABand(2)%wave, TABand(2)%sigmaS, TABand(2)%sigmaL, &
            3.7_DP, TABand(7)%wave, TABand(7)%sigmaS, TABand(7)%sigmaL, &
            3.6_DP, TABand(8)%wave, TABand(8)%sigmaS, TABand(8)%sigmaL, &
            2.5_DP, TABand(12)%wave, TABand(12)%sigmaS, TABand(12)%sigmaL, &
            1.5_DP, TABand(13)%wave, TABand(13)%sigmaS, TABand(13)%sigmaL, &
            1.6_DP, TABand(14)%wave, TABand(14)%sigmaS, TABand(14)%sigmaL, &
            2.7_DP, TABand(16)%wave, TABand(16)%sigmaS, TABand(16)%sigmaL, &
            3.8_DP, TABand(19)%wave, TABand(19)%sigmaS, TABand(19)%sigmaL, &
            3.9_DP, TABand(20)%wave, TABand(20)%sigmaS, TABand(20)%sigmaL, &
            2.9_DP, TABand(21)%wave, TABand(21)%sigmaS, TABand(21)%sigmaL, &
            1.8_DP, TABand(24)%wave, TABand(24)%sigmaS, TABand(24)%sigmaL, &
            1.7_DP, TABand(25)%wave, TABand(25)%sigmaS, TABand(25)%sigmaL, &
            2.6_DP, TABand(30)%wave, TABand(30)%sigmaS, TABand(30)%sigmaL, &
  !!        lnAv LOG[mag], &
            0._DP, &
  !!        lnFstar LOG[Lsun/pc2]]
            45._DP]

  Ncont = SIZE(labQ(:))
  Nband = SIZE(labB(:))
  Nline = SIZE(labL(:))
  Npar = 2*Ncont + 4*Nband + 3*Nline + 2
  
  ALLOCATE(parinfo(Npar))
  i0 = 0
  CONTinfo: DO i=1,Ncont
    parinfo(i0+2*i-1)%name = "lnMovd2"//TRIMLR(PRING(i))
    parinfo(i0+2*i-1)%comp = "CONT"
    ! parinfo(i0+2*i-1)%limits = [10._DP, 20._DP]
    ! parinfo(i0+2*i-1)%limited = [.TRUE., .TRUE.]
    parinfo(i0+2*i-1)%fixed = .FALSE.
    parinfo(i0+2*i)%name = "lnT"//TRIMLR(PRING(i))
    parinfo(i0+2*i)%comp = "CONT"
    parinfo(i0+2*i)%limits = [LOG(50._DP), LOG(500._DP)]
    parinfo(i0+2*i)%limited = [.TRUE., .TRUE.]
    parinfo(i0+2*i)%fixed = .FALSE.
  END DO CONTinfo
  i0 = i0 + 2*Ncont
  LINEinfo: DO i=1,Nline
    parinfo(i0+3*i-2)%name = "lnIline"//TRIMLR(PRING(i))
    parinfo(i0+3*i-2)%comp = "LINE"
    ! parinfo(i0+3*i-2)%limits = [0._DP, 5._DP]
    ! parinfo(i0+3*i-2)%limited = [.TRUE., .TRUE.]
    parinfo(i0+3*i-2)%fixed = .FALSE.
    parinfo(i0+3*i-1)%name = "Cline"//TRIMLR(PRING(i))
    parinfo(i0+3*i-1)%comp = "LINE"
    parinfo(i0+3*i-1)%fixed = .TRUE.
    parinfo(i0+3*i)%name = "Wline"//TRIMLR(PRING(i))
    parinfo(i0+3*i)%comp = "LINE"
    parinfo(i0+3*i)%fixed = .TRUE.
  END DO LINEinfo
  i0 = i0 + 3*Nline
  BANDinfo: DO i=1,Nband
    parinfo(i0+4*i-3)%name = "lnIband"//TRIMLR(PRING(i))
    parinfo(i0+4*i-3)%comp = "BAND"
    ! parinfo(i0+4*i-3)%limits = [0._DP,5._DP]
    ! parinfo(i0+4*i-3)%limited = [.TRUE.,.TRUE.]
    parinfo(i0+4*i-3)%fixed = .FALSE.
    parinfo(i0+4*i-2)%name = "Cband"//TRIMLR(PRING(i))
    parinfo(i0+4*i-2)%comp = "BAND"
    parinfo(i0+4*i-2)%fixed = .TRUE.
    parinfo(i0+4*i-1)%name = "WSband"//TRIMLR(PRING(i))
    parinfo(i0+4*i-1)%comp = "BAND"
    parinfo(i0+4*i-1)%fixed = .TRUE.
    parinfo(i0+4*i)%name = "WLband"//TRIMLR(PRING(i))
    parinfo(i0+4*i)%comp = "BAND"
    parinfo(i0+4*i)%fixed = .TRUE.
  END DO BANDinfo
  i0 = i0 + 4*Nband
  parinfo(i0+1)%name = "lnAv1"
  parinfo(i0+1)%comp = "PABS"
  parinfo(i0+1)%fixed = .TRUE.
  i0 = i0 + 1
  parinfo(i0+1)%name = "lnFstar1"
  parinfo(i0+1)%comp = "STAR"
  ! parinfo(i0+1)%limits = [10._DP,20._DP]
  ! parinfo(i0+1)%limited = [.TRUE.,.TRUE.]
  parinfo(i0+1)%fixed = .FALSE.

  parinfo(:)%ind = [(i,i=1,Npar)]
  CALL set_indpar(ind, parinfo(:))
  
  !! Create grid
  !!-------------
  wOBS = RAMP(Nw, 1._DP, 40._DP, XLOG=.TRUE.)
  NwOBS = SIZE(wOBS(:))
  nuOBS = MKS%clight/MKS%micron / wOBS(:)

  !! Make Qabs
  CALL READ_MASTER(WAVALL=wOBS(:), &
                   ! VERBOSE=verbose, NiniMC=NiniMC, &
                   ! CALIB=calib, NEWSEED=newseed, NEWINIT=newinit, &
                   LABQ=labQ, QABS=Qabs)!, LABL=labL, LABB=labB, &
                   ! NCONT=Ncont, NBAND=Nband, NLINE=Nline, DOSTOP=dostop, &
                   ! PARINFO=parinfo, INDPAR=ind, NPAR=Npar)
  
  !! Build synthetic spectrum
  !!--------------------------
  ALLOCATE(FnuOBS(NwOBS), dFnuOBS(NwOBS))

  FnuOBS(:) = specModel(wOBS(:), INDPAR=ind, PARVAL=pargen(:), QABS=Qabs(:), verbose=.TRUE.)
  dFnuOBS(:) = 0.2_DP*ABS(FnuOBS(:))
  FnuOBS(:) = FnuOBS(:) + dFnuOBS(:) * RAND_NORM(NwOBS)

  CALL WRITE_HDF5(wOBS, NAME="Wavelength (microns)", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
  CALL WRITE_HDF5(dFnuOBS, NAME="FnuUNC (MKS)", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(FnuOBS, NAME="FnuOBS (MKS)", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  
  print*, 'Gen synthetic spectrum [done]'

  print*, '=================================================='

  !!----------
  !! Chi2 fit
  !!----------

  IF (verbose) PRINT*
  DO i=1,MERGE(2,1,verbose)
    CALL BANNER_PROGRAM("LE MIROIR: LEast-squares fitting of Mid-IR emission " &
                        //"OptImized Routine", UNIT=unitlog(i), SWING=.True.)

  END DO
  
  !! Initialize the parameters
  !!---------------------------
  ALLOCATE(parini(Npar,1))
  parini(:,:) = 0._DP

  !! parini(:,1) = [ lnMovd2 LOG[Msun/pc2], lnT LOG[K], &
  parini(:,1) = [47._DP, 4.5_DP, 46._DP, 5.2_DP, 45._DP, 5.9_DP, &
  !!             lnIline, Cline, Wline, &
                 1._DP, TABLine(3)%wave, degradeRes(TABLine(3)%wave, .01_DP, 'SL-LL'), & 
                 1._DP, TABLine(9)%wave, degradeRes(TABLine(9)%wave, .01_DP, 'SL-LL'), & 
                 1._DP, TABLine(11)%wave, degradeRes(TABLine(11)%wave, .01_DP, 'SL-LL'), &
                 1._DP, TABLine(20)%wave, degradeRes(TABLine(20)%wave, .01_DP, 'SL-LL'), &
                 1._DP, TABLine(23)%wave, degradeRes(TABLine(23)%wave, .01_DP, 'SL-LL'), &
                 1._DP, TABLine(24)%wave, degradeRes(TABLine(24)%wave, .01_DP, 'SL-LL'), &
                 1._DP, TABLine(26)%wave, degradeRes(TABLine(26)%wave, .01_DP, 'SL-LL'), &
                 1._DP, TABLine(27)%wave, degradeRes(TABLine(27)%wave, .01_DP, 'SL-LL'), &
                 1._DP, TABLine(29)%wave, degradeRes(TABLine(29)%wave, .01_DP, 'SL-LL'), &
                 1._DP, TABLine(33)%wave, degradeRes(TABLine(33)%wave, .01_DP, 'SL-LL'), &
                 1._DP, TABLine(34)%wave, degradeRes(TABLine(34)%wave, .01_DP, 'SL-LL'), &
                 1._DP, TABLine(36)%wave, degradeRes(TABLine(36)%wave, .01_DP, 'SL-LL'), &
  !!             lnIband, Cband, WSband, WLband, &
                 1._DP, TABand(1)%wave, TABand(1)%sigmaS, TABand(1)%sigmaL, & 
                 1._DP, TABand(2)%wave, TABand(2)%sigmaS, TABand(2)%sigmaL, & 
                 1._DP, TABand(7)%wave, TABand(7)%sigmaS, TABand(7)%sigmaL, & 
                 1._DP, TABand(8)%wave, TABand(8)%sigmaS, TABand(8)%sigmaL, & 
                 1._DP, TABand(12)%wave, TABand(12)%sigmaS, TABand(12)%sigmaL, & 
                 1._DP, TABand(13)%wave, TABand(13)%sigmaS, TABand(13)%sigmaL, & 
                 1._DP, TABand(14)%wave, TABand(14)%sigmaS, TABand(14)%sigmaL, & 
                 1._DP, TABand(16)%wave, TABand(16)%sigmaS, TABand(16)%sigmaL, & 
                 1._DP, TABand(19)%wave, TABand(19)%sigmaS, TABand(19)%sigmaL, & 
                 1._DP, TABand(20)%wave, TABand(20)%sigmaS, TABand(20)%sigmaL, &
                 1._DP, TABand(21)%wave, TABand(21)%sigmaS, TABand(21)%sigmaL, & 
                 1._DP, TABand(24)%wave, TABand(24)%sigmaS, TABand(24)%sigmaL, & 
                 1._DP, TABand(25)%wave, TABand(25)%sigmaS, TABand(25)%sigmaL, & 
                 1._DP, TABand(30)%wave, TABand(30)%sigmaS, TABand(30)%sigmaL, & 
  !!             lnAv LOG[mag], &
                 0._DP, & 
  !!             lnFstar LOG[Lsun/pc2]]
                 44._DP]

  ALLOCATE (itied(Npar))
  itied(:) = 0
  DO i=1,Npar
    IF (.NOT. TRIMEQ(parinfo(i)%tied,"")) &
      CALL IWHERE(TRIMEQ(parinfo(:)%name,parinfo(i)%tied),itied(i))
  END DO
  Nparfree = COUNT((.NOT. parinfo(:)%fixed) .AND. (itied(:) <= 0))
  
  !! Run the fitter
  ALLOCATE(par(Npar), Fnu_mod(NwOBS), resid(NwOBS), &
           limits(Npar,2), limited(Npar,2), parerr(Npar), covar(Npar,Npar))
  
  FORALL (i=1:Npar)
    limits(i,:) = parinfo(i)%limits(:)
    limited(i,:) = parinfo(i)%limited(:)
  END FORALL

  par(:) = parini(:,1)
  CALL chi2min_LM (residuals, NwOBS, PAR=par(:), &
                   TOL=tol, RESID=resid(:), STATUS=status, VERBOSE=debug, &
                   LIMITED=limited(:,:), LIMITS=limits(:,:), &
                   FIXED=parinfo(:)%fixed, ITIED=itied(:), &
                   PARNAME=parinfo(:)%name, CHI2RED=chi2red, NITER=Niter, &
                   PARERR=parerr(:), COVAR=covar(:,:))

  print*, 'Chi2 fit synt spec [done]'

  print*, '=================================================='
  
  !!----------
  !! Analysis
  !!----------
  
  !! Compute model result
  !!----------------------
  PRINT*
  PRINT*, "Checking output:"
  PRINT*, "  status = "//TRIMLR(PRING(status))
  ! chi2red = SUM(deviate(:)**2) / (Nobs-Nparfree)
  PRINT*, "  chi2red = "//TRIMLR(PRING(chi2red,NDEC=10))
  PRINT*, "  Niter = "//TRIMLR(PRING(Niter))
  CALL WRITE_ASCII(VEC1=par(:),VEC2=parerr(:),FILE="out/chi2min.txt")
  PRINT*
  PRINT*, "Covariance matrix:"
  ! DO i=1,Npar
  !   PRINT*, REAL(covar(:,i), KIND(0.))
  ! END DO
  PRINT*

  CALL WRITE_HDF5(DBLARR1D=par(:), NAME='Best fitted parameter value', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

  !! Calculate model
  !!-----------------
  ALLOCATE(FnuMOD(NwOBS))
  
  FnuMOD(:) = specModel(wOBS(:), INDPAR=ind, PARVAL=par(:), QABS=Qabs(:), &
                        FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
                        PABS=Pabs, FNULINE=FnuLINE)

  CALL WRITE_HDF5(DBLARR1D=FnuCONT*Pabs, NAME="FnuCONT (MKS)", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=FnuLINE+(FnuCONT+FnuSTAR)*Pabs, NAME="FnuLINE (MKS)", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=(FnuBAND+FnuCONT+FnuSTAR)*Pabs, NAME="FnuBAND (MKS)", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=FnuSTAR*Pabs, NAME="FnuSTAR (MKS)", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=FnuMOD, NAME="FnuMOD (MKS)", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  
  print*, 'fitChi2syn analysis [done]'

  print*, '=================================================='
  
  !! Final
  !!-------
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) "PROGRAM EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
    WRITE(unitlog(i),*)
  END DO
  IF (verbose) THEN
    PRINT*, " - File "//fiLOG//" has been written."
    PRINT*
  END IF

  ! DEALLOCATE(parinfo, FnuOBS, dFnuOBS, wOBS, Fnu_mod, resid, parini, itied, par)
  
END PROGRAM test_fitChi2syn
