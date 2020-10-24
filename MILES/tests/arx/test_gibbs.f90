MODULE gibbs_external

  USE auxil, ONLY: parinfo_type, indpar_type, Qabs_type
  USE utilities, ONLY: DP
  USE inout, ONLY: lenpar
  USE grain_optics, ONLY: lendustQ
  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PUBLIC :: NwOBS, ipar
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wOBS, nuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parcurr
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: FnuOBS, dFnuOBS

  !! Init param
  !!------------
  TYPE(indpar_type), SAVE, PUBLIC                             :: indpar
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parinfo
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC    :: Qabs
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC            :: indBIN, indLIN
  
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

    USE auxil, ONLY: specModel
    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: pargrid
    REAL(DP), DIMENSION(SIZE(pargrid)) :: lnLHobs_par

    INTEGER :: iw
    REAL(DP), DIMENSION(SIZE(pargrid),NwOBS) :: Fnu_mod, varred, lnpind

    !! Model
    Fnu_mod(:,:) = specModel(wOBS(:), PARVEC=pargrid(:), PARNAME=parinfo(ipar)%name, &
                             PARINFO=parinfo(:), INDPAR=indpar, PARVAL=parcurr(:), &
                             QABS=Qabs, INDBIN=indBIN, INDLIN=indLIN)

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

  USE auxil, ONLY: read_master, degradeRes, set_indpar, specModel
  USE datable, ONLY: BIN, LIN
  USE utilities, ONLY: DP, pring, trimLR, swap, timinfo, strike, &
                       ustd, isNaN, warning, &
                       initiate_clock, time_type
  USE arrays, ONLY: ramp
  USE constants, ONLY: MKS
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, write_ascii, ascext, &
                   lenpar, lenpath, write_hdf5_frac, read_hdf5_frac
  USE grain_optics, ONLY: lendustQ
  USE random, ONLY: generate_newseed, rand_general, rand_norm
  USE statistics, ONLY: mean, sigma
  USE gibbs_external, ONLY: NwOBS, wOBS, nuOBS, &
                            FnuOBS, dFnuOBS, indpar, &
                            ipar, parcurr, &
                            lnpost_par, lnlhobs_par, &
                            parinfo, Qabs, indBIN, indLIN
  IMPLICIT NONE

  INTEGER :: Ncont, Nband, Nline

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: Ngibbsmax = 3000
  REAL(DP), PARAMETER :: accrand = 1.E-3_DP
  CHARACTER(*), PARAMETER :: nampar = "Parameter values"
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Gen spec
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labBIN, labLIN
  INTEGER, PARAMETER :: Nw=500
  REAL(DP), DIMENSION(:), ALLOCATABLE :: pargen
  
  !! Input variables
  INTEGER :: i0, i
  INTEGER :: Npar, Nmcmc
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: parini, parmcmc
  LOGICAL :: verbose, newseed
  CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ
  !! MCMC parameters
  INTEGER :: icurr, iprev, counter
  REAL(DP), DIMENSION(2) :: lim
  CHARACTER(lenpar) :: filMCMC
  !! Output variables
  CHARACTER(lenpar) :: filOUT
  INTEGER, DIMENSION(2) :: unitlog
  TYPE(time_type) :: timestr

  !! Analysis
  INTEGER :: t_burnin, t_end, Nmcmc_eff
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meanpar, stdevpar
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: parpix
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: par2D
  REAL(DP), DIMENSION(:), ALLOCATABLE :: FnuMOD
  CHARACTER(lenpar) :: filANAL
  !! Tests
  ! INTEGER :: Ntestgrid, itest
  ! REAL(DP) :: pargen0
  ! REAL(DP), DIMENSION(:), ALLOCATABLE :: testgrid, lnLHobs
  ! REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Fnu1par

  !! Output settings
  CALL INITIATE_CLOCK(timestr)
  unitlog(:) = [ ulog, ustd ]

  !! Read the inputs
  !!-----------------
  labQ = (/'ACH2_Z96             ', &
           'BE_Z96               ', &
           'Sil_D03              '/)
  labBIN = (/'Main 3.3     '/)
  labLIN = (/'H2S7  '/)

  !! No READ_MASTER buffer
  verbose = .TRUE.
  newseed = .TRUE.
  Nmcmc = 200 ! Length of Markov chain Monte Carlo

  IF (newseed) CALL GENERATE_NEWSEED()

  !! Generate param
  !!----------------

  !! pargen = [lnMovd2 LOG[Msun/pc2], lnT LOG[K], &
  !! -3.9, 4.6, -3, 4.6, -3.5, 4.6
  pargen = [LOG(2.E-2_DP), LOG(100._DP), LOG(5.E-2_DP), LOG(100._DP), LOG(3.E-2_DP), LOG(100._DP), &
  !!        lnIline, Cline, Wline, &
  !! -19.6, 
            LOG(3.E-9_DP), LIN(3)%wave, degradeRes(LIN(3)%wave,.01_DP,'SL-LL'), &
  !!        lnIband, Cband, WSband, WLband, &
  !! -20
            LOG(2.E-9_DP), BIN(1)%wave, BIN(1)%sigmaS, BIN(1)%sigmaL, &
  !!        lnAv LOG[mag], &
            LOG(1._DP), &
  !!        lnFstar LOG[MJy/sr]]
  !! -5.5
            LOG(4.E-3_DP)]

  Ncont = SIZE(labQ(:))
  Nband = SIZE(labBIN(:))
  Nline = SIZE(labLIN(:))
  Npar = 2*Ncont + 4*Nband + 3*Nline + 2
  
  ALLOCATE(parinfo(Npar))
  i0 = 0
  CONTinfo: DO i=1,Ncont
    parinfo(i0+2*i-1)%name = "lnMovd2"//TRIMLR(PRING(i))
    parinfo(i0+2*i-1)%comp = "CONT"
    parinfo(i0+2*i-1)%limits = [-25._DP, 10._DP]
    parinfo(i0+2*i-1)%limited = [.TRUE., .TRUE.]
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
    parinfo(i0+3*i-2)%limits = [-25._DP, 10._DP]
    parinfo(i0+3*i-2)%limited = [.TRUE., .TRUE.]
    parinfo(i0+3*i-2)%fixed = .FALSE.
    parinfo(i0+3*i-1)%name = "Cline"//TRIMLR(PRING(i))
    parinfo(i0+3*i-1)%comp = "LINE"
    parinfo(i0+3*i)%name = "Wline"//TRIMLR(PRING(i))
    parinfo(i0+3*i)%comp = "LINE"
    parinfo(i0+3*i)%fixed = .TRUE.
  END DO LINEinfo
  i0 = i0 + 3*Nline
  BANDinfo: DO i=1,Nband
    parinfo(i0+4*i-3)%name = "lnIband"//TRIMLR(PRING(i))
    parinfo(i0+4*i-3)%comp = "BAND"
    parinfo(i0+4*i-3)%limits = [-25._DP,10._DP]
    parinfo(i0+4*i-3)%limited = [.TRUE.,.TRUE.]
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
  parinfo(i0+1)%name = "lnAv"
  parinfo(i0+1)%comp = "PABS"
  parinfo(i0+1)%fixed = .TRUE.
  i0 = i0 + 1
  parinfo(i0+1)%name = "lnFstar"
  parinfo(i0+1)%comp = "STAR"
  parinfo(i0+1)%limits = [-25._DP,10._DP]
  parinfo(i0+1)%limited = [.TRUE.,.TRUE.]
  parinfo(i0+1)%fixed = .FALSE.

  parinfo(:)%ind = [(i,i=1,Npar)]
  CALL set_indpar(indpar, parinfo(:))

  !! Create grid
  !!-------------
  wOBS = RAMP(Nw, 3._DP, 40._DP, XLOG=.TRUE.)
  NwOBS = SIZE(wOBS(:))
  nuOBS = MKS%clight/MKS%micron / wOBS(:)

  CALL READ_MASTER(LABQ=labQ(:), &
                   LABBIN=labBIN(:), LABLIN=labLIN(:), &
                   WAVEALL=wOBS(:), &
                   QABS=Qabs, INDBIN=indBIN, INDLIN=indLIN)

  !! Build synthetic spectrum
  ALLOCATE(FnuOBS(NwOBS), dFnuOBS(NwOBS))

  FnuOBS(:) = specModel(wOBS(:), INDPAR=indpar, PARVAL=pargen(:), QABS=Qabs(:))
  dFnuOBS(:) = RAND_NORM(NwOBS) * 0.01_DP*MAXVAL(FnuOBS(:))
  FnuOBS(:) = FnuOBS(:) + dFnuOBS(:)

  filOUT = 'out/test_gibbs'
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

  ! CALL WRITE_HDF5(INITDBLARR=[NwOBS,5], FILE=TRIMLR(filOUT)//h5ext, &
  !                 NAME="var1par", COMPRESS=compress, verbose=debug, &
  !                 APPEND=.TRUE.)
  ! ALLOCATE(Fnu1par(NwOBS,5))

  ! itest = 1
  ! pargen0 = pargen(itest)
  ! DO i=1,5
  !   ! pargen(itest) = pargen0 * i!**.25_DP
  !   pargen(itest) = pargen0 + LOG(i**1._DP) ! LOG cases
  !   Fnu1par(:,i) = specModel(wOBS(:), INDPAR=indpar, PARVAL=pargen(:), QABS=Qabs(:))
  !   CALL WRITE_HDF5(DBLARR2D=Fnu1par(:,i:i), FILE=TRIMLR(filOUT)//h5ext, &
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
  !! pargen = [lnMovd2 LOG[Msun/pc2], lnT LOG[K], &
  parini(:,1) = [LOG(1.E-2_DP), LOG(50._DP), LOG(1.E-2_DP), LOG(50._DP), LOG(1.E-2_DP), LOG(50._DP), &
  !!        lnIline, Cline, Wline, &
            LOG(1.E-9_DP), LIN(3)%wave, degradeRes(LIN(3)%wave,.01_DP,'SL-LL'), &
  !!        lnIband, Cband, WSband, WLband, &
            LOG(1.E-9_DP), BIN(1)%wave, BIN(1)%sigmaS, BIN(1)%sigmaL, &
  !!        lnAv LOG[mag], &
            LOG(1._DP), &
  !!        lnFstar LOG[MJy/sr]]
            LOG(1.E-3_DP)]
! parini(:,1) = pargen
  !! Range of parameters for Gibbs sampling
  dolimit: DO i=1,Npar
    IF (.NOT. parinfo(i)%limited(1)) &
      parinfo(i)%limits(1) = parini(i,1) - 10._DP
    IF (.NOT. parinfo(i)%limited(2)) &
      parinfo(i)%limits(2) = parini(i,1) + 10._DP
    
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
  filMCMC = 'out/parlog_gibbs'
  CALL WRITE_HDF5(INTARR1D=[Nmcmc], FILE=TRIMLR(filMCMC)//h5ext, &
                  NAME="Length of MCMC", COMPRESS=compress, verbose=debug, &
                  APPEND=.FALSE.)
  CALL WRITE_HDF5(INITDBLARR=[Npar,Nmcmc], FILE=TRIMLR(filMCMC)//h5ext, &
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
! print*, "lnMovd21 = (-3.9) ", parcurr(1)!, parcurr(3), parcurr(5)
! print*, "lnT1 = (4.6) ", parcurr(2)!, parcurr(4), parcurr(6)
! print*, "lnIline1 = (-20) ", parcurr(7)
! print*, "lnIband1 = (-20) ", parcurr(10)
! print*, "lnFstar = (9.9) ", parcurr(15)
! print*, "max val = ", MAXVAL(specModel(wOBS(:), &

! IF (ipar==1 .AND. counter==1) THEN
  ! Ntestgrid = 1
  ! testgrid = [parini(ipar,1)]
  ! testgrid = RAMP(Ntestgrid, -25._DP, 10._DP, XLOG=.FALSE.)
  ! ALLOCATE(lnLHobs(Ntestgrid))
  ! lnLHobs(:) = lnpost_par(testgrid(:))

  ! CALL WRITE_HDF5(DBLARR1D=testgrid(:), FILE=TRIMLR(filOUT)//h5ext, &
  !                 NAME="testgrid", COMPRESS=compress, verbose=debug, &
  !                 APPEND=.TRUE.)
  ! CALL WRITE_HDF5(DBLARR1D=lnLHobs(:), FILE=TRIMLR(filOUT)//h5ext, &
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
    CALL WRITE_HDF5(DBLARR2D=parmcmc(:,icurr:icurr), FILE=TRIMLR(filMCMC)//h5ext, &
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
  IF (t_burnin <= 0 .OR. t_burnin >= t_end) t_burnin = t_end/10
  Nmcmc_eff = t_end - t_burnin + 1

  ALLOCATE(parpix(Npar,t_end), meanpar(Npar), stdevpar(Npar))

  CALL READ_HDF5(DBLARR2D=par2D, FILE=TRIMLR(filMCMC)//h5ext, &
                 NAME=nampar, IND2=[1,t_end])
  parpix(:,:) = par2D(:,:)
  
  meanpar(:) = MEAN(parpix(:,t_burnin:t_end), DIM=2)
  stdevpar(:) = SIGMA(parpix(:,t_burnin:t_end), DIM=2)

  filANAL = filOUT
  CALL WRITE_HDF5(DBLARR1D=meanpar(:), FILE=TRIMLR(filANAL)//h5ext, &
                  NAME="Mean of parameter value", COMPRESS=compress, verbose=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=stdevpar(:), FILE=TRIMLR(filANAL)//h5ext, &
                  NAME="Sigma of parameter value", COMPRESS=compress, verbose=debug, &
                  APPEND=.TRUE.)

  !! Calculate model
  !!-----------------
  ALLOCATE(FnuMOD(NwOBS))
  
  FnuMOD(:) = specModel(wOBS(:), INDPAR=indpar, PARVAL=meanpar(:), QABS=Qabs(:))
  
  CALL WRITE_HDF5(DBLARR1D=FnuMOD, FILE=TRIMLR(filOUT)//h5ext, &
                  NAME="FnuMOD (MJyovsr)", COMPRESS=compress, verbose=debug, &
                  APPEND=.TRUE.)

  print*, 'gibbs analysis [done]'

  print*, '=================================================='
    
END PROGRAM test_gibbs
