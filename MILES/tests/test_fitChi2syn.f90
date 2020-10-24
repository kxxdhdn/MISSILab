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
    
    residuals(:) = FnuOBS(:) - Fnu_mod(:)

  END FUNCTION residuals

END MODULE fitChi2syn_external


!!==========================================================================
!!                           Main Program
!!==========================================================================


PROGRAM test_fitChi2syn

  USE auxil, ONLY: read_master, degradeRes, set_indpar, specModel
  USE datable, ONLY: BIN, LIN
  USE utilities, ONLY: DP, pring, trimLR, timinfo, ustd, initiate_clock, time_type
  USE arrays, ONLY: ramp
  USE constants, ONLY: MKS
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, write_ascii, ascext, lenpar, lenpath
  USE grain_optics, ONLY: lendustQ
  USE random, ONLY: generate_newseed, rand_norm
  USE statistics, ONLY: mean, sigma
  USE chi2_minimization, ONLY: chi2min_LM
  USE fitChi2syn_external, ONLY: NwOBS, wOBS, nuOBS, FnuOBS, dFnuOBS, resid, &
                                 residuals, Fnu_mod, ind, parinfo, Qabs
  IMPLICIT NONE

  INTEGER :: Ncont, Nband, Nline

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  REAL(DP), PARAMETER :: tol = 1.E-10_DP
  CHARACTER(*), PARAMETER :: dirOUT = "./out/Chi2syn/"
  CHARACTER(*), PARAMETER :: filog = dirOUT//"log_fitChi2syn"//ascext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Gen spec
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labBIN, labLIN
  INTEGER, PARAMETER :: Nw=1000
  REAL(DP), DIMENSION(:), ALLOCATABLE :: pargen

  !! Input variables
  INTEGER :: i0, i, Npar, NiniMC
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: parini
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ
  LOGICAL :: verbose, calib, newseed, newinit, dostop
  !! chi2min_LM parameters
  ! LOGICAL, DIMENSION(:), ALLOCATABLE :: fixed
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: limited
  ! CHARACTER(20), DIMENSION(:), ALLOCATABLE :: parname
  INTEGER :: status, Nparfree
  ! INTEGER, DIMENSION(:), ALLOCATABLE :: itied
  REAL(DP) :: chi2red
  REAL(DP), DIMENSION(:), ALLOCATABLE :: par, parerr
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: limits
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: covar
  !! Output variables
  CHARACTER(lenpar) :: filOUT
  INTEGER, DIMENSION(2) :: unitlog
  TYPE(time_type) :: timestr
  
  !! Analysis
  REAL(DP), DIMENSION(:), ALLOCATABLE :: FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  CHARACTER(lenpar) :: filANAL
  
  !! Output settings
  !!-----------------
  CALL INITIATE_CLOCK(timestr)
  OPEN (ulog, FILE=filog, STATUS="REPLACE", ACTION="WRITE")
  unitlog(:) = [ ulog, ustd ]

  !!-----------------
  !! Read the inputs
  !!-----------------

  !! Create grid
  !!-------------
  wOBS = RAMP(Nw, 3._DP, 40._DP, XLOG=.TRUE.)
  NwOBS = SIZE(wOBS(:))
  nuOBS = MKS%clight/MKS%micron / wOBS(:)

  CALL READ_MASTER(WAVALL=wOBS(:), &
                   VERBOSE=verbose, NiniMC=NiniMC, &
                   LABQ=labQ, LABL=labL, LABB=labB, QABS=Qabs, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, DOSTOP=dostop, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar)

  !!---------------------
  !! Generate a spectrum
  !!---------------------

  !! pargen = [lnMovd2 LOG[Msun/pc2], lnT LOG[K], &
  pargen = [LOG(2.E-2_DP), LOG(100._DP), LOG(5.E-2_DP), LOG(100._DP), LOG(3.E-2_DP), LOG(100._DP), & 
  !!        lnIline, Cline, Wline, &
            LOG(1.E-9_DP), LIN(3)%wave, degradeRes(LIN(3)%wave,.01_DP,'SL-LL'), & 
            LOG(2.5E-9_DP), LIN(9)%wave, degradeRes(LIN(9)%wave,.01_DP,'SL-LL'), & 
            LOG(3.E-9_DP), LIN(11)%wave, degradeRes(LIN(11)%wave,.01_DP,'SL-LL'), &
            LOG(1.5E-9_DP), LIN(20)%wave, degradeRes(LIN(20)%wave,.01_DP,'SL-LL'), &
            LOG(1.5E-9_DP), LIN(23)%wave, degradeRes(LIN(23)%wave,.01_DP,'SL-LL'), &
            LOG(2.5E-9_DP), LIN(24)%wave, degradeRes(LIN(24)%wave,.01_DP,'SL-LL'), &
            LOG(1.5E-9_DP), LIN(26)%wave, degradeRes(LIN(26)%wave,.01_DP,'SL-LL'), &
            LOG(3.E-9_DP), LIN(27)%wave, degradeRes(LIN(27)%wave,.01_DP,'SL-LL'), &
            LOG(3.5E-9_DP), LIN(29)%wave, degradeRes(LIN(29)%wave,.01_DP,'SL-LL'), &
            LOG(4.5E-9_DP), LIN(33)%wave, degradeRes(LIN(33)%wave,.01_DP,'SL-LL'), &
            LOG(1.5E-9_DP), LIN(34)%wave, degradeRes(LIN(34)%wave,.01_DP,'SL-LL'), &
            LOG(1.E-9_DP), LIN(36)%wave, degradeRes(LIN(36)%wave,.01_DP,'SL-LL'), &
  !!        lnIband, Cband, WSband, WLband, &
            LOG(1.E-9_DP), BIN(1)%wave, BIN(1)%sigmaS, BIN(1)%sigmaL, &
            LOG(.8E-9_DP), BIN(2)%wave, BIN(2)%sigmaS, BIN(2)%sigmaL, &
            LOG(.5E-9_DP), BIN(7)%wave, BIN(7)%sigmaS, BIN(7)%sigmaL, &
            LOG(.7E-9_DP), BIN(8)%wave, BIN(8)%sigmaS, BIN(8)%sigmaL, &
            LOG(.3E-9_DP), BIN(12)%wave, BIN(12)%sigmaS, BIN(12)%sigmaL, &
            LOG(.5E-9_DP), BIN(13)%wave, BIN(13)%sigmaS, BIN(13)%sigmaL, &
            LOG(.4E-9_DP), BIN(14)%wave, BIN(14)%sigmaS, BIN(14)%sigmaL, &
            LOG(1.E-9_DP), BIN(16)%wave, BIN(16)%sigmaS, BIN(16)%sigmaL, &
            LOG(.3E-9_DP), BIN(19)%wave, BIN(19)%sigmaS, BIN(19)%sigmaL, &
            LOG(.7E-9_DP), BIN(20)%wave, BIN(20)%sigmaS, BIN(20)%sigmaL, &
            LOG(.6E-9_DP), BIN(21)%wave, BIN(21)%sigmaS, BIN(21)%sigmaL, &
            LOG(.5E-9_DP), BIN(24)%wave, BIN(24)%sigmaS, BIN(24)%sigmaL, &
            LOG(.5E-9_DP), BIN(25)%wave, BIN(25)%sigmaS, BIN(25)%sigmaL, &
            LOG(.9E-9_DP), BIN(30)%wave, BIN(30)%sigmaS, BIN(30)%sigmaL, &
  !!        lnAv LOG[mag], &
            LOG(1._DP), &
  !!        lnFstar LOG[Lsun/pc2]]
            LOG(2.E-3_DP)]

  !! Build synthetic spectrum
  !!--------------------------
  ALLOCATE(FnuOBS(NwOBS), dFnuOBS(NwOBS))

  FnuOBS(:) = specModel(wOBS(:), INDPAR=ind, PARVAL=pargen(:), QABS=Qabs(:))
  dFnuOBS(:) = RAND_NORM(NwOBS) * 0.01_DP*MAXVAL(FnuOBS(:))
  FnuOBS(:) = FnuOBS(:) + dFnuOBS(:)

  filOUT = 'out/test_fitChi2syn'
  CALL WRITE_HDF5(wOBS, FILE=TRIMLR(filOUT)//h5ext, &
                  NAME="Wavelength (microns)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.FALSE.)
  CALL WRITE_HDF5(dFnuOBS, FILE=TRIMLR(filOUT)//h5ext, &
                  NAME="FnuUNC (MJyovsr)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(FnuOBS, FILE=TRIMLR(filOUT)//h5ext, &
                  NAME="FnuOBS (MJyovsr)", &
                  COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  
  print*, 'Gen synthetic spectrum [done]'

  print*, '=================================================='

  !!----------
  !! Chi2 fit
  !!----------
  CALL READ_MASTER(WAVALL=wave, &
                   VERBOSE=verbose, NiniMC=NiniMC, &
                   CALIB=calib, NEWSEED=newseed, NEWINIT=newinit, &
                   LABQ=labQ, LABL=labL, LABB=labB, QABS=Qabs, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, DOSTOP=dostop, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar)

  IF (newseed) CALL GENERATE_NEWSEED()
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
  parini(:,1) = [LOG(1.E-2_DP), LOG(50._DP), LOG(1.E-2_DP), LOG(50._DP), LOG(1.E-2_DP), LOG(50._DP), & 
  !!             lnIline, Cline, Wline, &
                 LOG(1.E-9_DP), LIN(3)%wave, degradeRes(LIN(3)%wave, .01_DP, 'SL-LL'), & 
                 LOG(1.E-9_DP), LIN(9)%wave, degradeRes(LIN(9)%wave, .01_DP, 'SL-LL'), & 
                 LOG(1.E-9_DP), LIN(11)%wave, degradeRes(LIN(11)%wave, .01_DP, 'SL-LL'), &
                 LOG(1.E-9_DP), LIN(20)%wave, degradeRes(LIN(20)%wave, .01_DP, 'SL-LL'), &
                 LOG(1.E-9_DP), LIN(23)%wave, degradeRes(LIN(23)%wave, .01_DP, 'SL-LL'), &
                 LOG(1.E-9_DP), LIN(24)%wave, degradeRes(LIN(24)%wave, .01_DP, 'SL-LL'), &
                 LOG(1.E-9_DP), LIN(26)%wave, degradeRes(LIN(26)%wave, .01_DP, 'SL-LL'), &
                 LOG(1.E-9_DP), LIN(27)%wave, degradeRes(LIN(27)%wave, .01_DP, 'SL-LL'), &
                 LOG(1.E-9_DP), LIN(29)%wave, degradeRes(LIN(29)%wave, .01_DP, 'SL-LL'), &
                 LOG(1.E-9_DP), LIN(33)%wave, degradeRes(LIN(33)%wave, .01_DP, 'SL-LL'), &
                 LOG(1.E-9_DP), LIN(34)%wave, degradeRes(LIN(34)%wave, .01_DP, 'SL-LL'), &
                 LOG(1.E-9_DP), LIN(36)%wave, degradeRes(LIN(36)%wave, .01_DP, 'SL-LL'), &
  !!             lnIband, Cband, WSband, WLband, &
                 LOG(1.E-9_DP), BIN(1)%wave, BIN(1)%sigmaS, BIN(1)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(2)%wave, BIN(2)%sigmaS, BIN(2)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(7)%wave, BIN(7)%sigmaS, BIN(7)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(8)%wave, BIN(8)%sigmaS, BIN(8)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(12)%wave, BIN(12)%sigmaS, BIN(12)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(13)%wave, BIN(13)%sigmaS, BIN(13)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(14)%wave, BIN(14)%sigmaS, BIN(14)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(16)%wave, BIN(16)%sigmaS, BIN(16)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(19)%wave, BIN(19)%sigmaS, BIN(19)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(20)%wave, BIN(20)%sigmaS, BIN(20)%sigmaL, &
                 LOG(1.E-9_DP), BIN(21)%wave, BIN(21)%sigmaS, BIN(21)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(24)%wave, BIN(24)%sigmaS, BIN(24)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(25)%wave, BIN(25)%sigmaS, BIN(25)%sigmaL, & 
                 LOG(1.E-9_DP), BIN(30)%wave, BIN(30)%sigmaS, BIN(30)%sigmaL, & 
  !!             lnAv LOG[mag], &
                 LOG(1._DP), & 
  !!             lnFstar LOG[Lsun/pc2]]
                 LOG(1.E-3_DP)]

  Nparfree = COUNT((.NOT. parinfo(:)%fixed) .AND. (parinfo(:)%itied <= 0))
  
  !! Run the fitter
  ALLOCATE(par(Npar), Fnu_mod(NwOBS), resid(NwOBS), &
           limits(Npar,2), limited(Npar,2), parerr(Npar), covar(Npar,Npar))
  
  FORALL (i=1:Npar)
    limits(i,:) = parinfo(i)%limits(:)
    limited(i,:) = parinfo(i)%limited(:)
  END FORALL

  par(:) = parini(:,1)
  CALL chi2min_LM (residuals, NwOBS, PAR=par(:), &
                   TOL=tol, RESID=resid(:), STATUS=status, VERBOSE=.True., &
                   LIMITED=limited(:,:), LIMITS=limits(:,:), &
                   FIXED=parinfo(:)%fixed, ITIED=parinfo(:)%itied, &
                   PARNAME=parinfo(:)%name, CHI2RED=chi2red, &
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
  PRINT*, "  chi2 = "//TRIMLR(PRING(chi2red,NDEC=10))
  CALL write_ascii(VEC1=par(:),VEC2=parerr(:),FILE="out/chi2min.txt")
  PRINT*
  PRINT*, "Covariance matrix:"
  ! DO i=1,Npar
  !   PRINT*, REAL(covar(:,i), KIND(0.))
  ! END DO
  PRINT*

  !! Calculate model
  !!-----------------
  ALLOCATE(FnuMOD(NwOBS))
  
  FnuMOD(:) = specModel(wOBS(:), INDPAR=ind, PARVAL=par(:), QABS=Qabs(:), &
                        FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
                        PABS=Pabs, FNULINE=FnuLINE)

  filANAL = filOUT
  CALL WRITE_HDF5(DBLARR1D=FnuCONT*Pabs, FILE=TRIMLR(filANAL)//h5ext, &
                  NAME="FnuCONT (MJyovsr)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=FnuLINE+(FnuCONT+FnuSTAR)*Pabs, FILE=TRIMLR(filANAL)//h5ext, &
                  NAME="FnuLINE (MJyovsr)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=(FnuBAND+FnuCONT+FnuSTAR)*Pabs, FILE=TRIMLR(filANAL)//h5ext, &
                  NAME="FnuBAND (MJyovsr)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=FnuSTAR*Pabs, FILE=TRIMLR(filANAL)//h5ext, &
                  NAME="FnuSTAR (MJyovsr)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=FnuMOD, FILE=TRIMLR(filANAL)//h5ext, &
                  NAME="FnuMOD (MJyovsr)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  
  print*, 'fitChi2syn analysis [done]'

  print*, '=================================================='
  
  ! !! Final
  ! !!-------
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) "PROGRAM EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
    WRITE(unitlog(i),*)
  END DO
  IF (verbose) THEN
    PRINT*, " - File "//filog//" has been written."
    PRINT*
  END IF
  
END PROGRAM test_fitChi2syn
