MODULE gibbs_external

  USE core, ONLY: parinfo_type, indpar_type, Qabs_type
  USE utilities, ONLY: DP
  USE inout, ONLY: lenpar
  USE grain_optics, ONLY: lendustQ
  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PUBLIC :: NwOBS, ipar
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wOBS, nuOBS
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC :: extinct
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parcurr
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: FnuOBS, dFnuOBS

  !! Init param
  !!------------
  TYPE(indpar_type), SAVE, PUBLIC                             :: ind
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parinfo
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC    :: Qabs
  
  PUBLIC :: lnpost_par, lnlhobs_par

CONTAINS

  !!-------------------------------------------------------
  !!
  !!            Likelihoods of the observations
  !!
  !!-------------------------------------------------------
  
  !! For sampling physical parameters
  !!----------------------------------
  FUNCTION lnLHobs_par(pargrid)

    USE core, ONLY: specModel
    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: pargrid
    REAL(DP), DIMENSION(SIZE(pargrid)) :: lnLHobs_par

    INTEGER :: iw
    REAL(DP), DIMENSION(SIZE(pargrid),NwOBS) :: Fnu_mod, varred, lnpind

    !! Model
    Fnu_mod(:,:) = specModel(wOBS(:), PARVEC=pargrid(:), PARNAME=parinfo(ipar)%name, &
                             PARINFO=parinfo(:), INDPAR=ind, PARVAL=parcurr(:), &
                             QABS=Qabs, EXTINCT=extinct)

    !! Likelihoods
    varred(:,:) = 0._DP
    FORALL (iw=1:NwOBS) &
      varred(:,iw) = ( FnuOBS(iw) - Fnu_mod(:,iw) ) / dFnuOBS(iw)

    ! DO iw=1,NwOBS
      ! print*, "obs, min, max", iw, REAL(FnuOBS(iw)), REAL(MINVAL(Fnu_mod(:,iw))), REAL(MAXVAL(Fnu_mod(:,iw)))
    ! END DO

    lnpind(:,:) = - 0.5_DP * varred(:,:)**2

    !! Global likelihood, accounting for the excesses
    lnLHobs_par(:) = SUM(lnpind(:,:), DIM=2)

  END FUNCTION lnLHobs_par

  !!-------------------------------------------------------
  !!
  !!                Posterior distributions
  !!
  !!-------------------------------------------------------

  !! For sampling physical parameters
  !!----------------------------------
  FUNCTION lnpost_par(pargrid)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: pargrid
    REAL(DP), DIMENSION(SIZE(pargrid)) :: lnpost_par

    lnpost_par(:) = LNLHOBS_PAR(pargrid(:))

  END FUNCTION lnpost_par

END MODULE gibbs_external


!!==========================================================================
!!                           Main Program
!!==========================================================================


PROGRAM test_gibbs

  USE core, ONLY: read_master, degradeRes, set_indpar, specModel, make_Qabs, extCurve
  USE auxil, ONLY: TABLine, TABand
  USE utilities, ONLY: DP, pring, trimLR, swap, timinfo, strike, &
                       ustd, isNaN, warning, &
                       initiate_clock, time_type
  USE arrays, ONLY: ramp, closest
  USE constants, ONLY: MKS
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, write_ascii, ascext, &
                   lenpar, lenpath
  USE grain_optics, ONLY: lendustQ
  USE random, ONLY: generate_newseed, rand_general, rand_norm
  USE statistics, ONLY: mean, sigma
  USE gibbs_external, ONLY: NwOBS, wOBS, nuOBS, FnuOBS, dFnuOBS, ind, &
                            ipar, parcurr, lnpost_par, lnlhobs_par, parinfo, &
                            Qabs, extinct
  IMPLICIT NONE


  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: Ngibbsmax = 3000
  REAL(DP), PARAMETER :: accrand = 1.E-3_DP
  CHARACTER(*), PARAMETER :: dirOUT = "./out/"
  CHARACTER(*), PARAMETER :: filOUT = dirOUT//'test_gibbs'//h5ext
  CHARACTER(*), PARAMETER :: nampar = "Parameter values"
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Gen spec
  INTEGER, PARAMETER :: Nw=500
  REAL(DP), DIMENSION(:), ALLOCATABLE :: pargen
  
  !! Input variables
  INTEGER :: i0, i, Npar, iw0, Nmcmc!, NiniMC
  INTEGER :: Ncont, Nband, Nline, Nextc
  REAL(DP) :: refw
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: parini, parmcmc
  LOGICAL :: verbose, newseed!, calib, newinit, dostop
  CHARACTER(lendustQ) :: refB
  CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ, labB, labL, labE
  !! MCMC parameters
  INTEGER :: icurr, iprev, counter
  REAL(DP), DIMENSION(2) :: lim
  CHARACTER(lenpar) :: filMCMC = 'out/parlog_gibbs'//h5ext
  !! Output variables
  INTEGER, DIMENSION(2) :: unitlog
  TYPE(time_type) :: timestr

  !! Analysis
  INTEGER :: t_burnin, t_end, Nmcmc_eff
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meanpar, stdevpar
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: parpix
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: par2D
  REAL(DP), DIMENSION(:), ALLOCATABLE :: FnuMOD
  !! Tests
  ! INTEGER :: Ntestgrid, itest
  ! REAL(DP) :: pargen0
  ! REAL(DP), DIMENSION(:), ALLOCATABLE :: testgrid, lnLHobs
  ! REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Fnu1par

  !! Output settings
  CALL INITIATE_CLOCK(timestr)
  unitlog(:) = [ ulog, ustd ]

  !!-----------------
  !! Read the inputs
  !!-----------------
  labQ = (/'ACH2_Z96             ', &
           'Sil_D03              '/)
  labB = (/'Main 3.3     ', &
           'Main 6.2 (1) ', &
           'Main 6.2 (2) ', &
           'Plateau 7.7  ', &
           'Main 7.7 (1) ', &
           'Main 7.7 (2) ', &
           'Main 8.6     ', &
           'Small 11.0   ', &
           'Main 11.2    '/)
  labL = (/'H2S7  '/)
  labE = (/'D03'/)
  Ncont = SIZE(labQ(:))
  Nband = SIZE(labB(:))
  Nline = SIZE(labL(:))
  Nextc = SIZE(labE(:))

  !! Create grid
  !!-------------
  wOBS = RAMP(Nw, 2.5_DP, 40._DP, XLOG=.TRUE.)
  NwOBS = SIZE(wOBS(:))
  nuOBS = MKS%clight/MKS%micron / wOBS(:)

  !! No READ_MASTER buffer
  !!-----------------------
  !! Make Qabs
  ALLOCATE(Qabs(Ncont), extinct(Nextc,NwOBS))
  CALL MAKE_QABS(LABEL=labQ(:), QABS=Qabs(:), WAVALL=wOBS(:))
  !! Extinction curve
  extinct(1,:) = extCurve(labE(1),wOBS(:),.TRUE.)
  !! Set ref band
  refB = 'Main 11.2    '
  !! Set ref wavelength
  refw = 15.0_DP
  iw0 = CLOSEST(wOBS(:), refw)
  !! Set master
  verbose = .TRUE.
  newseed = .FALSE.
  Nmcmc = 100 ! Length of Markov chain Monte Carlo

  IF (newseed) CALL GENERATE_NEWSEED()

  !!----------------
  !! Generate param
  !!----------------

  !! pargen = [lnFcont LOG[Msun/pc2], lnT LOG[K], &
  pargen = [4._DP, 6._DP, 3._DP, 4.2_DP, &
  !!        lnRline, Cline, Wline, &
            -0.5_DP, TABLine(3)%wave, degradeRes(TABLine(3)%wave,.01_DP,'SL-LL'), &
  !!        lnRband, Cband, WSband, WLband, &
            -0.8_DP, TABand(1)%wave, TABand(1)%sigmaS, TABand(1)%sigmaL, &
            -0.2_DP, TABand(8)%wave, TABand(8)%sigmaS, TABand(8)%sigmaL, &
            -0.8_DP, TABand(9)%wave, TABand(9)%sigmaS, TABand(9)%sigmaL, &
            -0.2_DP, TABand(13)%wave, TABand(13)%sigmaS, TABand(13)%sigmaL, &
            -0.2_DP, TABand(14)%wave, TABand(14)%sigmaS, TABand(14)%sigmaL, &
            -0.4_DP, TABand(15)%wave, TABand(15)%sigmaS, TABand(15)%sigmaL, &
            -0.8_DP, TABand(17)%wave, TABand(17)%sigmaS, TABand(17)%sigmaL, &
            -0.5_DP, TABand(20)%wave, TABand(20)%sigmaS, TABand(20)%sigmaL, &
            34.0_DP, TABand(21)%wave, TABand(21)%sigmaS, TABand(21)%sigmaL, &
  !!        lnAv LOG[mag], &
            0._DP, &
  !!        lnFstar LOG[Lsun/pc2]]
            42._DP]

  Npar = 2*Ncont + 4*Nband + 3*Nline + Nextc + 1
  
  ALLOCATE(parinfo(Npar))
  i0 = 0
  CONTinfo: DO i=1,Ncont
    parinfo(i0+2*i-1)%name = "lnFcont"//TRIMLR(PRING(i))
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
    parinfo(i0+3*i-2)%name = "lnRline"//TRIMLR(PRING(i))
    parinfo(i0+3*i-2)%comp = "LINE"
    ! parinfo(i0+3*i-2)%limits = [-25._DP, 10._DP]
    ! parinfo(i0+3*i-2)%limited = [.TRUE., .TRUE.]
    parinfo(i0+3*i-2)%fixed = .FALSE.
    parinfo(i0+3*i-1)%name = "Cline"//TRIMLR(PRING(i))
    parinfo(i0+3*i-1)%comp = "LINE"
    parinfo(i0+3*i)%name = "Wline"//TRIMLR(PRING(i))
    parinfo(i0+3*i)%comp = "LINE"
    parinfo(i0+3*i)%fixed = .TRUE.
  END DO LINEinfo
  i0 = i0 + 3*Nline
  BANDinfo: DO i=1,Nband
    parinfo(i0+4*i-3)%name = "lnRband"//TRIMLR(PRING(i))
    parinfo(i0+4*i-3)%comp = "BAND"
    ! parinfo(i0+4*i-3)%limits = [-25._DP,10._DP]
    ! parinfo(i0+4*i-3)%limited = [.TRUE.,.TRUE.]
    parinfo(i0+4*i-3)%fixed = .FALSE.
    parinfo(i0+4*i-2)%name = "Cband"//TRIMLR(PRING(i))
    parinfo(i0+4*i-2)%comp = "BAND"
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
  i0 = i0 + Nextc
  parinfo(i0+1)%name = "lnFstar1"
  parinfo(i0+1)%comp = "STAR"
  ! parinfo(i0+1)%limits = [-25._DP,10._DP]
  ! parinfo(i0+1)%limited = [.TRUE.,.TRUE.]
  parinfo(i0+1)%fixed = .FALSE.

  parinfo(:)%ind = [(i,i=1,Npar)]
  CALL set_indpar(INDPAR=ind, PARINFO=parinfo(:), &
                  REFB=refB, LABB=labB, REFW=iw0, LABQ=labQ)
  
  !! Build synthetic spectrum
  !!--------------------------
  ALLOCATE(FnuOBS(NwOBS), dFnuOBS(NwOBS))

  FnuOBS(:) = specModel(wOBS(:), INDPAR=ind, PARVAL=pargen(:), &
                        QABS=Qabs(:), EXTINCT=extinct(:,:), verbose=.TRUE.)
  dFnuOBS(:) = 0.2_DP*ABS(FnuOBS(:))
  FnuOBS(:) = FnuOBS(:) + dFnuOBS(:) * RAND_NORM(NwOBS)

  CALL WRITE_HDF5(wOBS, FILE=filOUT, NAME="Wavelength (microns)", &
                  COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
  CALL WRITE_HDF5(dFnuOBS, FILE=filOUT, NAME="FnuUNC (MJyovsr)", &
                  COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(FnuOBS, FILE=filOUT, NAME="FnuOBS (MJyovsr)", &
                  COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(pargen, FILE=filOUT, NAME="True parameter value", &
                  COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

  print*, 'Gen synthetic spectrum [done]'

  print*, '=================================================='

  ! CALL WRITE_HDF5(INITDBLARR=[NwOBS,5], FILE=filOUT, &
  !                 NAME="var1par", COMPRESS=compress, verbose=debug, &
  !                 APPEND=.TRUE.)
  ! ALLOCATE(Fnu1par(NwOBS,5))

  ! itest = 1
  ! pargen0 = pargen(itest)
  ! DO i=1,5
  !   ! pargen(itest) = pargen0 * i!**.25_DP
  !   pargen(itest) = pargen0 + LOG(i**1._DP) ! LOG cases
  !   Fnu1par(:,i) = specModel(wOBS(:), INDPAR=ind, PARVAL=pargen(:), &
                     ! QABS=Qabs(:), EXTINCT=extinct(:,:))
  !   CALL WRITE_HDF5(DBLARR2D=Fnu1par(:,i:i), FILE=filOUT, &
  !                   NAME="var1par", COMPRESS=compress, verbose=debug, &
  !                   IND2=[i,i])

  ! END DO

  ! print*, 'PARVEC (', TRIMLR(parinfo(itest)%name), ') varying test [done]'

  ! print*, '=================================================='


  
  !!------------------
  !! Prepare the MCMC
  !!------------------

  !! Initialize the parameters
  !!---------------------------
  ALLOCATE(parini(Npar,1))
  
  parini(:,:) = 0._DP
  !! pargen = [lnFcont LOG[Msun/pc2], lnT LOG[K], &
  parini(:,1) = [5._DP, 5._DP, 5._DP, 5._DP, &
  !!        lnRline, Cline, Wline, &
            0._DP, TABLine(3)%wave, degradeRes(TABLine(3)%wave,.01_DP,'SL-LL'), &
  !!        lnRband, Cband, WSband, WLband, &
            0._DP, TABand(1)%wave, TABand(1)%sigmaS, TABand(1)%sigmaL, &
            0._DP, TABand(8)%wave, TABand(8)%sigmaS, TABand(8)%sigmaL, &
            0._DP, TABand(9)%wave, TABand(9)%sigmaS, TABand(9)%sigmaL, &
            0._DP, TABand(13)%wave, TABand(13)%sigmaS, TABand(13)%sigmaL, &
            0._DP, TABand(14)%wave, TABand(14)%sigmaS, TABand(14)%sigmaL, &
            0._DP, TABand(15)%wave, TABand(15)%sigmaS, TABand(15)%sigmaL, &
            0._DP, TABand(17)%wave, TABand(17)%sigmaS, TABand(17)%sigmaL, &
            0._DP, TABand(20)%wave, TABand(20)%sigmaS, TABand(20)%sigmaL, &
            30._DP, TABand(21)%wave, TABand(21)%sigmaS, TABand(21)%sigmaL, &
  !!        lnAv LOG[mag], &
            0._DP, &
  !!        lnFstar LOG[Lsun/pc2]]
            40._DP]
! parini(:,1) = pargen
  !! Range of parameters for Gibbs sampling
  dolimit: DO i=1,Npar
    IF (.NOT. parinfo(i)%limited(1)) &
      parinfo(i)%limits(1) = parini(i,1) - 5._DP
    IF (.NOT. parinfo(i)%limited(2)) &
      parinfo(i)%limits(2) = parini(i,1) + 5._DP
    
  END DO dolimit

  !! Initialize the MCMC
  !!---------------------

  !! Initialize the parameters of the chain
  ALLOCATE(parmcmc(Npar,2))
  DO i=1,2!!
    parmcmc(:,i) = parini(:,1)
  END DO
  
  !! Outputs, counter and allocation
  !!---------------------------------

  !! Initiate counters
  icurr = 2
  iprev = 1
  counter = 1

  !! Allocate arrays
  ALLOCATE(parcurr(Npar))

  !! Initialize the arrays
  CALL WRITE_HDF5(INTARR1D=[Nmcmc], FILE=filMCMC, &
                  NAME="Length of MCMC", COMPRESS=compress, verbose=debug, &
                  APPEND=.FALSE.)
  CALL WRITE_HDF5(INITDBLARR=[Npar,Nmcmc], FILE=filMCMC, &
                  NAME=nampar, COMPRESS=compress, verbose=debug, &
                  APPEND=.TRUE.)

  !!--------------
  !! Run the MCMC
  !!--------------

  !! Banner
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) "RUNNING THE MARKOV CHAIN MONTE-CARLO (" &
                        //TRIMLR(TIMINFO(timestr))//")"
    WRITE(unitlog(i),*)
  END DO

  !! MCMC
  !!======
  big_loop: DO
    ! IF ((MODULO(counter,printevery) == 0)) THEN
      DO i=1,MERGE(2,1,verbose)
        WRITE(unitlog(i),*) " Markov Chain Monte Carlo iteration " &
                            //TRIMLR(PRING(counter)) //" (" &
                            //TRIMLR(TIMINFO(timestr))//")"
      END DO
    ! END IF

    !! Physical parameters
    !!---------------------
    parmcmc(:,icurr) = parmcmc(:,iprev)
    ! parmcmc(:,:,:,icurr) = parmcmc(:,:,:,iprev)
    param: DO ipar=1,Npar
      IF (.NOT. parinfo(ipar)%fixed) THEN
        lim(:) = parinfo(ipar)%limits(:)
! print*, "limit", ipar, lim
        parcurr(:) = parmcmc(:,icurr)
! print*, counter, ipar, parinfo(ipar)%name
! print*, "lnFcont1 = (-3.9) ", parcurr(1)!, parcurr(3), parcurr(5)
! print*, "lnT1 = (4.6) ", parcurr(2)!, parcurr(4), parcurr(6)
! print*, "lnRline1 = (-20) ", parcurr(7)
! print*, "lnRband1 = (-20) ", parcurr(10)
! print*, "lnFstar = (9.9) ", parcurr(15)
! print*, "max val = ", MAXVAL(specModel(wOBS(:), &

! IF (ipar==1 .AND. counter==1) THEN
  ! Ntestgrid = 1
  ! testgrid = [parini(ipar,1)]
  ! testgrid = RAMP(Ntestgrid, -25._DP, 10._DP, XLOG=.FALSE.)
  ! ALLOCATE(lnLHobs(Ntestgrid))
  ! lnLHobs(:) = lnpost_par(testgrid(:))

  ! CALL WRITE_HDF5(DBLARR1D=testgrid(:), FILE=filOUT, &
  !                 NAME="testgrid", COMPRESS=compress, verbose=debug, &
  !                 APPEND=.TRUE.)
  ! CALL WRITE_HDF5(DBLARR1D=lnLHobs(:), FILE=filOUT, &
  !                 NAME="lnLHobs", COMPRESS=compress, verbose=debug, &
  !                 APPEND=.TRUE.)
  ! STOP
! END IF
! print*, 'LH dist test [done]'
! print*, '=================================================='

        parmcmc(ipar,icurr) &
          = RAND_GENERAL(lnpost_par, lim(:), YLOG=.TRUE., VERBOSE=debug, &
                         LNFUNC=.TRUE., ACCURACY=accrand, NMAX=Ngibbsmax)
        IF (debug) PRINT*, "par(",TRIMLR(PRING(ipar)),") = ", parmcmc(ipar,icurr)

      END IF
    END DO param

    !! Outputs
    !!---------
    CALL WRITE_HDF5(DBLARR2D=parmcmc(:,icurr:icurr), FILE=filMCMC, &
                    NAME=nampar, COMPRESS=compress, verbose=debug, &
                    IND2=[counter,counter])

    !! Check output for debugging
    IF (debug) THEN
      PRINT*, "Counter = ", counter
      PRINT*
    END IF
    
    !! Normal termination condition
    IF (counter == Nmcmc) EXIT

    !! Update the indices
    CALL SWAP(icurr,iprev)
    counter = counter + 1
  
  END DO big_loop

  print*, 'gibbs mcmc [done]'

  print*, '=================================================='
    
  !!----------
  !! Analysis
  !!----------
  
  !! Compute the statistics
  !!------------------------
  t_burnin = 2000000
  t_end = Nmcmc ! end of Markov chain Monte Carlo
  IF (t_burnin <= 0 .OR. t_burnin >= t_end) t_burnin = t_end/5
  Nmcmc_eff = t_end - t_burnin + 1

  ALLOCATE(parpix(Npar,t_end), meanpar(Npar), stdevpar(Npar))

  CALL READ_HDF5(DBLARR2D=par2D, FILE=filMCMC, &
                 NAME=nampar, IND2=[1,t_end])
  parpix(:,:) = par2D(:,:)
  
  meanpar(:) = MEAN(parpix(:,t_burnin:t_end), DIM=2)
  stdevpar(:) = SIGMA(parpix(:,t_burnin:t_end), DIM=2)

  CALL WRITE_HDF5(DBLARR1D=meanpar(:), FILE=filOUT, NAME="Mean of parameter value", &
                  COMPRESS=compress, verbose=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=stdevpar(:), FILE=filOUT, NAME="Sigma of parameter value", &
                  COMPRESS=compress, verbose=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(INTARR1D=[t_burnin], FILE=filOUT, NAME="t_burnin", &
                  COMPRESS=compress, verbose=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(INTARR1D=[t_end], FILE=filOUT, NAME="t_end", &
                  COMPRESS=compress, verbose=debug, APPEND=.TRUE.)

  !! Calculate model
  !!-----------------
  ALLOCATE(FnuMOD(NwOBS))
  
  FnuMOD(:) = specModel(wOBS(:), INDPAR=ind, PARVAL=meanpar(:), &
                        QABS=Qabs(:), EXTINCT=extinct(:,:), verbose=.TRUE.)
  
  CALL WRITE_HDF5(DBLARR1D=FnuMOD, FILE=filOUT, &
                  NAME="FnuMOD (MJyovsr)", COMPRESS=compress, verbose=debug, &
                  APPEND=.TRUE.)

  print*, 'gibbs analysis [done]'

  print*, '=================================================='
    
END PROGRAM test_gibbs
