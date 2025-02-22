!******************************************************************************
!*
!*                         FIT MIR SPECTRA (HB)
!*
!******************************************************************************


PROGRAM fit_hb

  USE utilities, ONLY: DP, isNaN, pring, trimLR, trimeq, swap, tinyDP, &
                       banner_program, ustd, strike, warning, &
                       timinfo, initiate_clock, time_type, today
  USE arrays, ONLY: iwhere, closest, reallocate
  USE constants, ONLY: MKS
  USE matrices, ONLY: invert_cholesky
  USE statistics, ONLY: mean, sigma, correlate, correl_index, corr2Rmat
  USE inout, ONLY: read_hdf5, write_hdf5, check_hdf5, h5ext, ascext, lenpar, lenpath
  USE grain_optics, ONLY: lendustQ
  USE random, ONLY: generate_newseed, rand_general
  USE core, ONLY: read_master, initparam, lencorr, specModel, set_indcal
  USE ext_hb, ONLY: Nx, Ny, NwOBS, wOBS, nuOBS, xOBS, yOBS, &
                    NwFIT, iwFIT, wFIT, nuFIT, &
                    maskwall, maskparall, maskhypall, &
                    FnuOBS, dFnuOBS, cenlnS0, siglnS0, &
                    i2ih, ipar, ihpar, icorr, icorr2ij, &
                    parcurr, parhypcurr, allparhypcurr, &
                    cov_prev, invcov_prev, invcorr_prev, &
                    mucurr, sigcurr, corrcurr, Nparhyp, Ncorrhyp, &
                    lnpost_par, Fnu_model, &
                    detcov_prev, detcorr, noposdef_prev, &
                    lnpost_mu, lnpost_sig, lnpost_corr, covariance, &
                    mask, maskpar, maskhyp, maskhypcurr, maskS, &!maskextra, &
                    ind, parinfo, parhypinfo, Qabs, extinct, labS, &
                    robust_RMS, skew_RMS, robust_cal, calib, &
                    jcal, iw2ical, Ncalib, sigcal, calibool, NwCAL, &
                    delp1, lnpost_del, specOBS ! Calibration errors
  IMPLICIT NONE

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: textwid = 60
  INTEGER, PARAMETER :: Ngibbsmax = 3000
  REAL(DP), PARAMETER :: Nsigcalran = 3._DP
  REAL(DP), PARAMETER :: accrand = 1.E-3_DP
  CHARACTER(*), PARAMETER :: nampar = "Parameter values"
  CHARACTER(*), PARAMETER :: namln1pd = "ln(1+delta)"
  CHARACTER(*), PARAMETER :: nammu = "Mean of hyperdistribution"
  CHARACTER(*), PARAMETER :: namsig = "Sigma of hyperdistribution"
  CHARACTER(*), PARAMETER :: namcorr = "Correlation of hyperdistribution"
  CHARACTER(*), PARAMETER :: namind = "Last index"
  LOGICAL, DIMENSION(2), PARAMETER :: limT = [ .TRUE., .TRUE. ]
  ! LOGICAL, DIMENSION(2), PARAMETER :: limF = [ .FALSE., .FALSE. ]
  LOGICAL, PARAMETER :: compress = .TRUE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Input variables
  INTEGER :: x, y, i, ih, ical, indhyp, indres0
  INTEGER :: Npar, Nmcmc, NiniMC, Nsou, Nparfree
  INTEGER :: Ncont, Nband, Nline, Nextc, Nstar, Nextra
  INTEGER, DIMENSION(:), ALLOCATABLE :: indresume, itied
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: maskint ! convert mask=0 to mask=T
  LOGICAL :: verbose, newseed, newinit, dostop, resume, nohi
  CHARACTER(lenpath) :: dirIN, dirOUT, filOBS, filOUT, filMCMC, fiLOG
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lencorr), DIMENSION(:), ALLOCATABLE :: corrhypname
  CHARACTER(lenpath), DIMENSION(:), ALLOCATABLE :: path1d
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB, labL

  !! MCMC parameters
  INTEGER :: icurr, iprev, counter
  REAL(DP), DIMENSION(2) :: lim, limlnS0, limR
  REAL(DP), DIMENSION(:), ALLOCATABLE :: mu0, sig0, corr0
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mumcmc, sigmcmc, corrmcmc
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: limln1pd
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: ln1pd0
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parini, parmcmc, ln1pdmcmc
  LOGICAL :: fileOK
  LOGICAL, DIMENSION(:), ALLOCATABLE :: fixed, hyper

  !! Output variables
  INTEGER :: Ncorr
  INTEGER, DIMENSION(2) :: unitlog
  TYPE(time_type) :: timestr

  
  !! Fixed seed
  !!------------
  ! INTEGER :: k, n
  ! INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  ! REAL(DP) :: omega
  
  ! CALL RANDOM_SEED(SIZE = n)
  ! ALLOCATE(seed(n))
  ! seed = (/ (3*k+6, k=1,n) /)
  ! CALL RANDOM_SEED(PUT = seed)
  ! CALL RANDOM_SEED(GET = seed)
  
  ! CALL RANDOM_NUMBER(omega)
  ! PRINT*, omega
  ! STOP

  
  !!------------------------------------------------------------------------
  !!                            I. Read the inputs
  !!------------------------------------------------------------------------

  CALL READ_HDF5(STRARR1D=path1d, FILE='../out/set_input'//h5ext, NAME='input dir')
  dirIN = path1d(1)
  filOBS = TRIMLR(dirIN)//'observation_MIR'//h5ext
  
  !! 1) Main settings
  !!------------------
  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='wavelength (microns)', N1=NwOBS)
  CALL READ_MASTER(WAVALL=wOBS(:), DIRIN=dirIN, DIROUT=dirOUT, &
                   VERBOSE=verbose, NMCMC=Nmcmc, NINIMC=NiniMC, &
                   ROBUST_RMS=robust_RMS, SKEW_RMS=skew_RMS, &
                   CALIB=calib, ROBUST_CAL=robust_cal, &
                   NEWSEED=newseed, NEWINIT=newinit, &
                   LABL=labL, LABB=labB, QABS=Qabs, EXTINCT=extinct, LABS=labS, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, &
                   NEXTC=Nextc, NSTAR=Nstar, NEXTRA=Nextra, NCORR=Ncorr, &
                   DOSTOP=dostop, RESUME=resume, INDRESUME=indres0, NOHI=nohi, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, &
                   CORRHYPNAME=corrhypname, SPEC_UNIT=spec_unit, &
                   PARHYPINFO=parhypinfo, NPARHYP=Nparhyp, NCORRHYP=Ncorrhyp)

  !! Output settings
  IF (nohi) THEN
    filOUT = TRIMLR(dirOUT)//'fit_bb'//h5ext
    filMCMC = TRIMLR(dirOUT)//'parlog_fit_bb'//h5ext
    fiLOG = TRIMLR(dirOUT)//'log_fit_bb'//ascext
  ELSE
    filOUT = TRIMLR(dirOUT)//'fit_hb'//h5ext
    filMCMC = TRIMLR(dirOUT)//'parlog_fit_hb'//h5ext
    fiLOG = TRIMLR(dirOUT)//'log_fit_hb'//ascext
  END IF
  
  CALL INITIATE_CLOCK(timestr)
  OPEN (ulog,FILE=fiLOG,STATUS=MERGE("OLD    ","REPLACE",resume), &
        ACTION="WRITE",POSITION=MERGE("APPEND","REWIND",resume))
  unitlog(:) = [ ulog, ustd ]

  IF (newseed) CALL GENERATE_NEWSEED()
  IF (verbose) PRINT*
  DO i=1,MERGE(2,1,verbose)
    CALL BANNER_PROGRAM("HISTOIRE (HIerarchical bayeSian fitting Tool Of mid-IR Emission)" &
                        //NEW_LINE('')// "written by D. Hu"//NEW_LINE(''), &
                        UNIT=unitlog(i), SWING=.TRUE.)
    IF (resume) THEN
      WRITE(unitlog(i),*) 
      WRITE(unitlog(i),*) REPEAT("*",80)
      WRITE(unitlog(i),*) 
      WRITE(unitlog(i),*) "HISTOIRE: ...RESUMING WHERE THE RUN STOPPED!!!"
      WRITE(unitlog(i),*) 
      WRITE(unitlog(i),"(A80)") ADJUSTR("("//TRIMLR(TODAY())//")")
      WRITE(unitlog(i),*) 
    END IF
  END DO
  
  !! 2) Observational sample
  !!-------------------------
  CALL READ_HDF5(STRARR1D=specOBS, FILE=filOBS, &
                 NAME='spectroscopic module labels')
  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, &
                 NAME='FnuOBS ('//TRIMLR(spec_unit)//')', N1=Nx, N2=Ny)
  !! On input, the convention is mask=1 => blocks; mask=0 => passes.
  !! However, within Fortran, the convention is the Fortran mask=T => passes.
  CALL READ_HDF5(INTARR3D=maskint, FILE=filOBS, NAME='NaN mask')

  !! Masks
  !! - MASK(Nx,Ny,NwOBS) is the 3D mask of the observed spectra. It is False if
  !!   one pixel at one wavelength is missing or flagged.
  !! - MASKEXTRA(Nx,Ny,Nextra) is the 3D mask of the observed extra parameters.
  !!   It is False if one pixel has a missing or flagged extra parameter.
  !! - MASKHYP(Nx,Ny,Nparhyp) selects the parameters of pixels which having to
  !!   be treated hierarchically. For extra parameters, it is MASKEXTRA.
  !! - MASKPAR(Nx,Ny,Npar) is the same as MASKHYP but for the parameter list.
  !! - MASKWALL(NwOBS) selects wavelengths presented in at least one pixel
  !! - MASKHYPALL(Nparhyp) selects hyper parameter not masked in all pixels
  ALLOCATE(mask(Nx,Ny,NwOBS))
  mask(:,:,:) = ( maskint(:,:,:) == 0 )
  WHERE (isNaN(FnuOBS(:,:,:)) .OR. .NOT. mask(:,:,:)) FnuOBS(:,:,:) = 0._DP

  !! Noise
  CALL READ_HDF5(DBLARR3D=dFnuOBS, FILE=filOBS, &
                 NAME='dFnuOBS ('//TRIMLR(spec_unit)//')')
  WHERE ( isNaN(dFnuOBS(:,:,:)) .OR. dFnuOBS<0._DP )
    dFnuOBS(:,:,:) = 0._DP
    mask(:,:,:) = .FALSE.
  END WHERE

  !! Not all-NaN wavelength grid wFIT
  ALLOCATE(maskwall(NwOBS))
  maskwall(:) = ANY(ANY(mask(:,:,:),DIM=1),DIM=1)
  NwFIT = COUNT(maskwall(:))
  CALL IWHERE(maskwall(:),iwFIT) ! i=1,NwFIT => iwFIT(i)=1,NwOBS
  ALLOCATE(wFIT(NwFIT))
  wFIT(:) = wOBS(iwFIT(:)) ! new wavelength grid excluding NaN wvl

  ALLOCATE(nuOBS(NwOBS),nuFIT(NwFIT))
  nuOBS(:) = MKS%clight/MKS%micron / wOBS(:)
  nuFIT(:) = MKS%clight/MKS%micron / wFIT(:)

  !! Read the instrumental covariance matrix
  calibread: IF (calib) THEN
    Ncalib = SIZE(specOBS)
    ALLOCATE (iw2ical(Ncalib)) ! iw=1,NwOBS => iw2ical(ical)%indw=1,Nwcal
    ALLOCATE(calibool(Ncalib,NwOBS),sigcal(Ncalib))
    DO ical=1,Ncalib
      CALL set_indcal(iw2ical(ical)%arr1d, wOBS(:), INSTR=specOBS(ical), &
                      CALIBOOL=calibool(ical,:), CALIBERR=sigcal(ical))
    END DO
  ELSE
    Ncalib = 1
    ALLOCATE(iw2ical(Ncalib))
    ALLOCATE(iw2ical(1)%arr1d(NwOBS))
    ALLOCATE(calibool(Ncalib,NwOBS),sigcal(Ncalib))
    iw2ical(1)%arr1d(:) = [(i,i=1,NwOBS)]
    calibool(1,:) = .TRUE.
    sigcal(1) = 0._DP
  END IF calibread
  
  !! Indices of tied parameters
  ALLOCATE(itied(Npar))
  itied(:) = 0
  DO i=1,Npar
    IF (.NOT. TRIMEQ(parinfo(i)%tied,"")) &
      CALL IWHERE(TRIMEQ(parinfo(:)%name,parinfo(i)%tied),itied(i))
  END DO

  !! Number of free parameters
  Nparfree = COUNT( (.NOT. parinfo(:)%fixed) .AND. (itied(:) == 0) )

  !! Read the extra parameters (TBD)

  !! 3) Summary of the sample properties
  !!-------------------------------------
  Nsou = COUNT(ANY(mask(:,:,:),DIM=3))
  summary: IF (.NOT. resume) THEN
    DO i=1,MERGE(2,1,verbose)
      WRITE(unitlog(i),*)
      WRITE(unitlog(i),*) " - Number of spectra/pixels to be fitted = " &
                          //TRIMLR(PRING(Nsou))
      WRITE(unitlog(i),*) " - Size of spectral sampling = " &
                          //TRIMLR(PRING(NwFIT))
      IF (calib) THEN
        WRITE(unitlog(i),*) " - Calibration errors are taken into " &
                            //"account"
      ELSE
        WRITE(unitlog(i),*) " - Calibration errors are not taken " &
                            //"into account"
      END IF
      WRITE(unitlog(i),*) " - Nparfree = "//TRIMLR(PRING(Nparfree))
      IF (.NOT. nohi) THEN
        WRITE(unitlog(i),*) " - Nparhyp = "//TRIMLR(PRING(Nparhyp))
        WRITE(unitlog(i),*) " - Ncorrhyp = "//TRIMLR(PRING(Ncorrhyp))
      END IF
      WRITE(unitlog(i),*)
    END DO
  END IF summary

  !! 4) aux
  !!--------
  !! Conversion between indices of par and parhyp
  !! - parhyp(:) = par(i2ih(:)) ; SIZE(i2ih)=Nparhyp ; MAXVAL(i2ih)<=Npar
  ALLOCATE (i2ih(Nparhyp))
  DO ih=1,Nparhyp 
    CALL IWHERE(TRIMEQ(parinfo(:)%name,parhypinfo(ih)%name),indhyp)
    i2ih(ih) = indhyp
  END DO

  !! 3D hyper mask
  ALLOCATE (maskhyp(Nx,Ny,Nparhyp))
  DO i=1,Nparhyp
    IF (parinfo(i2ih(i))%model) THEN
      maskhyp(:,:,i) = ANY(mask(:,:,:),DIM=3)
    ELSE
      maskhyp(:,:,i) = .FALSE.
    END IF
  END DO

  !! 3D parameter mask
  ALLOCATE (maskpar(Nx,Ny,Npar))
  DO i=1,Npar
    IF (parinfo(i)%model) THEN
      maskpar(:,:,i) = ANY(mask(:,:,:),DIM=3)
    ELSE
      maskpar(:,:,i) = .FALSE.
    END IF
  END DO
  
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'HISTOIRE: Read the inputs [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO
  
  !!------------------------------------------------------------------------
  !!                          II. Prepare the MCMC
  !!------------------------------------------------------------------------
  
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*) "ESTIMATING THE INITIAL VALUES OF THE PARAMETERS (" &
                        //TRIMLR(TIMINFO(timestr))//")"//NEW_LINE('')
  END DO

  !! a. Rough spectrum estimators
  !!------------------------------

  !! b. Automatic initial parameter estimates
  !!------------------------------------------
  ALLOCATE(parini(Nx,Ny,Npar,1))
  parini(:,:,:,:) = 0._DP
  
  CALL INITPARAM(NiniMC, IND=ind, PAR=parini(:,:,:,:), PARINFO=parinfo(:), &
                 ITIED=itied(:), MASK=maskpar(:,:,:), &
                 NEWINIT=newinit, FILOBS=filOBS, &
                 LABB=labB(:), LABL=labL(:), QABS=Qabs(:), LABS=labS(:))

  !! Update 3D hyper mask
  FORALL (x=1:Nx,y=1:Ny,i=1:Nparhyp,.NOT.maskpar(x,y,i2ih(i))) &
    maskhyp(x,y,i) = .FALSE.

  ALLOCATE(maskparall(Npar),maskhypall(Nparhyp))
  maskparall(:) = ANY(ANY(maskpar(:,:,:),DIM=1),DIM=1)
  maskhypall(:) = ANY(ANY(maskhyp(:,:,:),DIM=1),DIM=1)
  
  !! c. Range of parameters for Gibbs sampling
  !!-------------------------------------------
  !! Calibration errors
  ALLOCATE (limln1pd(Ncalib,2))
  IF (calib) THEN
    limln1pd(:,1) = - Nsigcalran * sigcal(:)
    limln1pd(:,2) = + Nsigcalran * sigcal(:)
  END IF
  
  !! Parameters
  dolimit: DO i=1,Npar
    IF (.NOT. parinfo(i)%limited(1)) THEN
        parinfo(i)%limits(1) &
          = MERGE( MINVAL(PACK(parini(:,:,i,1),maskpar(:,:,i))) &
                     - 3._DP*SIGMA(PACK(parini(:,:,i,1),maskpar(:,:,i))), &
                   parini(1,1,i,1) - 5._DP, (Nsou > 1) )
    END IF
    IF (.NOT. parinfo(i)%limited(2)) THEN
        parinfo(i)%limits(2) &
          = MERGE( MAXVAL(PACK(parini(:,:,i,1),maskpar(:,:,i))) &
                     + 3._DP*SIGMA(PACK(parini(:,:,i,1),maskpar(:,:,i))), &
                   parini(1,1,i,1) + 5._DP, (Nsou > 1) )
    END IF
  END DO dolimit

  !! d. Initialize the MCMC
  !!------------------------
  IF (nohi) parinfo(:)%hyper = .FALSE.

  !! Initial hyper parameters
  IF (ANY(parinfo(:)%hyper)) THEN
    ALLOCATE (mu0(Nparhyp),sig0(Nparhyp),corr0(Ncorrhyp))
    mu0(:) = 0._DP
    sig0(:) = 1._DP
    FORALL (i=1:Nparhyp,maskhypall(i))
      mu0(i) = MEAN(parini(:,:,i2ih(i),1),MASK=maskhyp(:,:,i))
      sig0(i) = SIGMA(parini(:,:,i2ih(i),1),MASK=maskhyp(:,:,i))
    END FORALL
    WHERE (ABS(sig0(:)) < 1.E-5_DP*ABS(mu0(:)/sig0(:))) &
      sig0(:) = MERGE( ABS(mu0(:)), 1._DP, (mu0(:) > 0._DP) )
    WHERE (sig0(:) <= 0._DP) &
      sig0(:) = 0.5_DP * ( parinfo(i2ih(:))%limits(2) &
                           - parinfo(i2ih(:))%limits(1) )
    corr0(:) = 0._DP
  ELSE
    ALLOCATE (mu0(Npar),sig0(Npar),corr0(Ncorr))
    mu0(:) = 0._DP
    sig0(:) = 1._DP
    FORALL (i=1:Nparhyp,maskhypall(i))
      mu0(i) = MEAN(parini(:,:,i2ih(i),1),MASK=maskpar(:,:,i))
      sig0(i) = SIGMA(parini(:,:,i2ih(i),1),MASK=maskpar(:,:,i))
    END FORALL
    WHERE (ABS(sig0(:)) < tinyDP) &
      sig0(:) = MERGE( ABS(mu0(:)), 1._DP, (mu0(:) > 0._DP) )
    WHERE (sig0(:) <= 0._DP) &
      sig0(:) = 0.5_DP * ( parinfo(:)%limits(2) - parinfo(:)%limits(1) )
    corr0(:) = 0._DP
  END IF
  ALLOCATE (ln1pd0(Nx,Ny,Ncalib))
  ln1pd0(:,:,:) = 0._DP

  !! Print the initial values
  IF (.NOT. nohi) THEN
    DO i=1,MERGE(2,1,debug)
      WRITE(unitlog(i),*) "mu0 = ", mu0
      WRITE(unitlog(i),*) "sig0 = ", sig0
      ! WRITE(unitlog(i),*) "corr0 = ", corr0
      ! WRITE(unitlog(i),*) "ln1pd = ", ln1pd0
    END DO
  END IF
  
  !! Initialize the parameters of the chain
  ALLOCATE (parmcmc(Nx,Ny,Npar,2),ln1pdmcmc(Nx,Ny,Ncalib,2))
  DO i=1,2 
    FORALL (x=1:Nx,y=1:Ny,ANY(maskpar(x,y,:)))
      parmcmc(x,y,:,i) = parini(x,y,:,1)
      ln1pdmcmc(x,y,:,i) = ln1pd0(x,y,:)
    END FORALL
  END DO

  !! Initialize the hyperparameters. If the run is non hierarchical, we
  !! replace the hyperparameters, by the measured average and covariance of the 
  !! parameters, for visualisation. However, in this case, there is no parameter
  !! prior.
  ALLOCATE (mumcmc(Nparhyp,2),sigmcmc(Nparhyp,2),corrmcmc(Ncorrhyp,2))
  IF (ANY(parinfo(:)%hyper)) THEN
    DO i=1,2
      mumcmc(:,i) = mu0(:)
      sigmcmc(:,i) = sig0(:)
      corrmcmc(:,i) = corr0(:)
    END DO
  END IF
  CALL CORREL_INDEX(Nparhyp,Ncorrhyp,icorr2ij)

  !! Covariance matrix
  IF (ANY(parinfo(:)%hyper)) THEN
    limlnS0(:) = [ -1._DP, 1._DP ] * siglnS0
    ALLOCATE (cenlnS0(Nparhyp))
    cenlnS0(:) = LOG(sig0(:))
    limR(:) = [ -0.999_DP, 0.999_DP ]
  END IF

  !! Mask for the standard deviation matrix
  IF (ANY(parinfo(:)%hyper)) THEN
    ALLOCATE (maskS(Nparhyp,Nparhyp))
    maskS(:,:) = .FALSE.
    FORALL (ihpar=1:Nparhyp) maskS(ihpar,ihpar) = .TRUE.
  END IF
  
  !! e. Output files, counter and allocation
  !!-----------------------------------------

  !! Initiate counters
  counter = 1
  icurr = 2
  iprev = 1
  
  !! Allocate arrays
  ALLOCATE (Fnu_model(Nx,Ny,NwOBS)) ! spectral model
  ! ALLOCATE (ln1pd(Nx,Ny,Ncalib)) ! LN(1+delta) for the current iteration
  ALLOCATE (delp1(Nx,Ny,NwOBS)) ! delta+1 for the current iteration
  ALLOCATE (parcurr(Npar)) ! current parameters
  ALLOCATE (fixed(Npar), hyper(Npar))
  IF (ANY(parinfo(:)%hyper)) THEN
    ALLOCATE (allparhypcurr(Nx,Ny,Nparhyp))  ! parameter grid
    ALLOCATE (cov_prev(Nparhyp,Nparhyp))     ! correlation matrix
    ALLOCATE (invcov_prev(Nparhyp,Nparhyp))  ! inverse correlation matrix
    ALLOCATE (invcorr_prev(Nparhyp,Nparhyp)) ! inverse correlation matrix
    ALLOCATE (mucurr(Nparhyp))               ! current average
    ALLOCATE (sigcurr(Nparhyp))              ! current average
    ALLOCATE (corrcurr(Ncorrhyp))            ! current correlation coefficient
    ALLOCATE (parhypcurr(Nparhyp))           ! current parameters
    ALLOCATE (maskhypcurr(Nparhyp))          ! current parameters' mask
  END IF

  !! Open the output file
  fileresume: IF (.NOT. resume) THEN

    !! Initialize the files
    CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME="Parameter label", &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
    CALL WRITE_HDF5(STRARR1D=parhypinfo(:)%name,NAME="Hyper parameter label",&
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(STRARR1D=corrhypname,NAME="Correlation label",&
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(INTARR1D=[Nmcmc], NAME="Length of MCMC", &
                    FILE=filMCMC, COMPRESS=compress, verbose=debug, APPEND=.TRUE.)
    !! - parameters
    CALL WRITE_HDF5(INITDBLARR=[Nx,Ny,Npar,Nmcmc], NAME=nampar, &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    !! - calibration errors
    CALL WRITE_HDF5(INITDBLARR=[Nx,Ny,Ncalib,Nmcmc], NAME=namln1pd, &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    !! - hyperparameters or moments of the parameters (/NOHYPER)
    CALL WRITE_HDF5(INITDBLARR=[Nparhyp,Nmcmc], NAME=nammu, &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(INITDBLARR=[Nparhyp,Nmcmc], NAME=namsig, &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(INITDBLARR=[Ncorrhyp,Nmcmc], NAME=namcorr, &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    !! - Current index
    CALL WRITE_HDF5(INITINTARR=[1], NAME=namind, &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  ELSE

    !! Resume the run from one index before the stop (Last index)
    fileOK = CHECK_HDF5(FILE=filMCMC)
    IF (fileOK) THEN
      !! Auto-resume
      CALL READ_HDF5(INTARR1D=indresume,FILE=filMCMC,NAME=namind)
      
      IF (indres0 > 0) THEN
        IF (indres0 > indresume(1)) THEN
          CALL WARNING( 'HISTOIRE', &
            'Resuming from an index after the stop! ' &
            //'Returned to the last index = '//TRIMLR(PRING(indresume(1))) )
        ELSE
          !! Input-resume
          CALL REALLOCATE(indresume,1)
          indresume(1) = indres0
          
        END IF
      END IF
      
      IF (counter < indresume(1)) counter = indresume(1)

      IF (counter > 1) THEN
        !! Read the lastest values
        !! - Parameters
        CALL READ_HDF5(DBLARR4D=parmcmc,FILE=filMCMC,NAME=nampar, &
                       IND4=[counter-1,counter])
        !! - calibration errors
        CALL READ_HDF5(DBLARR4D=ln1pdmcmc,FILE=filMCMC,NAME=namln1pd, &
                       IND4=[counter-1,counter])
        !! - Hyperparameters or moments of parameters (BB)
        CALL READ_HDF5(DBLARR2D=mumcmc,FILE=filMCMC,NAME=nammu, &
                       IND2=[counter-1,counter])
        CALL READ_HDF5(DBLARR2D=sigmcmc,FILE=filMCMC,NAME=namsig, &
                       IND2=[counter-1,counter])
        CALL READ_HDF5(DBLARR2D=corrmcmc,FILE=filMCMC,NAME=namcorr, &
                       IND2=[counter-1,counter])

      END IF
    END IF

  END IF fileresume
  
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'Prepare the MCMC [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO

  !!------------------------------------------------------------------------
  !!                          III. Run the MCMC
  !!------------------------------------------------------------------------

  !! Banner
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*) "RUNNING THE MARKOV CHAIN MONTE-CARLO (" &
                        //TRIMLR(TIMINFO(timestr))//")"//NEW_LINE('')
  END DO

  !! Backup initial parinfo(:)%fixed and %hyper
  fixed(:) = parinfo(:)%fixed
  hyper(:) = parinfo(:)%hyper

  !! MCMC
  !!======
  big_loop: DO
    DO i=1,MERGE(2,1,verbose)
      WRITE(unitlog(i),*) " Markov Chain Monte Carlo iteration " &
                          //TRIMLR(PRING(counter)) //" (" &
                          //TRIMLR(TIMINFO(timestr))//")"
    END DO

    !! 1) Calibration errors
    !!-----------------------
    delp1(:,:,:) = 1._DP
    calibration: IF (calib) THEN
      ln1pdmcmc(:,:,:,icurr) = ln1pdmcmc(:,:,:,iprev)

      !! a. Compute the model for the current parameters
      Fnu_model(:,:,:) &
        = specModel( wOBS(:), INDPAR=ind, PARVAL=parmcmc(:,:,:,iprev), &
                     MASK=mask(:,:,:), &
                     QABS=Qabs(:), EXTINCT=extinct(:,:), LABS=labS(:) )

      !! b. Draw the calibration errors
      DO jcal=1,Ncalib
        NwCAL = COUNT(calibool(jcal,:))
        
        DO xOBS=1,Nx
          DO yOBS=1,Ny
            notmasked: IF (ANY( mask(xOBS,yOBS,:) )) THEN
              !! Take the first as ref
              ref: IF (jcal/=1) THEN
                IF (debug) PRINT*, " - ln(1+d("//TRIMLR(specOBS(jcal))//"))"
                ! ln1pd(xOBS,yOBS,jcal) = ln1pdmcmc(xOBS,yOBS,jcal,icurr)
                
                ln1pdmcmc(xOBS,yOBS,jcal,icurr) &
                  = RAND_GENERAL(lnpost_del, limln1pd(jcal,:), &
                                 YLOG=.TRUE., VERBOSE=debug, &
                                 LNFUNC=.TRUE., ACCURACY=accrand, &
                                 LIMITED=limT(:), NMAX=Ngibbsmax)
                IF (debug) PRINT*, "ln1pd(",xOBS,yOBS,jcal,") = ", &
                                   ln1pdmcmc(xOBS,yOBS,jcal,icurr)
                
              END IF ref
              WHERE (calibool(jcal,:)) &
                delp1(xOBS,yOBS,:) = EXP(ln1pdmcmc(xOBS,yOBS,jcal,icurr))
            END IF notmasked
          END DO
        END DO
      END DO
      
    END IF calibration
    ! DO i=1,MERGE(2,1,verbose)
    !   WRITE(unitlog(i),*) 'Calibration error sampling [done] - '// &
    !                         TRIMLR(TIMINFO(timestr))
    ! END DO
    
    !! 2) Physical parameters
    !!------------------------

    !! Covariance matrix of the hyperparameters
    IF (ANY(hyper(:))) THEN
      cov_prev(:,:) = COVARIANCE(sigmcmc(:,iprev),corrmcmc(:,iprev))
      invcov_prev(:,:) = INVERT_CHOLESKY(cov_prev(:,:))
      invcorr_prev(:,:) = INVERT_CHOLESKY( &
                            corr2Rmat(corrmcmc(:,iprev),Nparhyp), &
                            DETERMINANT=detcorr, &
                            NOPOSDEF=noposdef_prev ) ! new
      mucurr(:) = mumcmc(:,iprev)
    END IF

    !! Draw the model parameters
    parmcmc(:,:,:,icurr) = parmcmc(:,:,:,iprev)
    param: DO ipar=1,Npar
      IF (.NOT. fixed(ipar)) THEN
        IF (debug) PRINT*, " - "//TRIMLR(parinfo(ipar)%name)

        !! Even the parinfo-not-limited par need this limits
        !! See c. Range of parameters for Gibbs sampling
        lim(:) = parinfo(ipar)%limits(:)

        !! Update parinfo(:)%fixed with maskpar
        IF (maskpar(xOBS,yOBS,ipar)) THEN
          parinfo(ipar)%fixed = fixed(ipar)
          parinfo(ipar)%hyper = hyper(ipar)
        ELSE
          parinfo(ipar)%fixed = .TRUE.
          parinfo(ipar)%hyper = .FALSE.
        END IF
        xsource: DO xOBS=1,Nx
          ysource: DO yOBS=1,Ny
            IF (maskpar(xOBS,yOBS,ipar)) THEN
              parcurr(:) = parmcmc(xOBS,yOBS,:,icurr)
              IF (hyper(ipar)) THEN
                maskhypcurr(:) = maskhyp(xOBS,yOBS,:)
                WHERE (maskhyp(xOBS,yOBS,:))
                  parhypcurr(:) = parmcmc(xOBS,yOBS,i2ih(:),icurr)
                ELSEWHERE
                  parhypcurr(:) = mucurr(:)
                END WHERE
              ENDIF
              
              parmcmc(xOBS,yOBS,ipar,icurr) &
                = RAND_GENERAL(lnpost_par, lim(:), YLOG=.TRUE., VERBOSE=debug, &
                               LNFUNC=.TRUE., ACCURACY=accrand, NMAX=Ngibbsmax)

              IF (debug) PRINT*, "par(",xOBS,yOBS,ipar,") = ", &
                                 parmcmc(xOBS,yOBS,ipar,icurr)
            END IF
          END DO ysource
        END DO xsource
      END IF
    END DO param
    ! PRINT*, 'parmcmc(1,1,1,icurr)', parmcmc(1,1,1,icurr) ! with fixed seed
    ! DO i=1,MERGE(2,1,verbose)
    !   WRITE(unitlog(i),*) 'Parameter sampling [done] - '// &
    !                         TRIMLR(TIMINFO(timestr))
    ! END DO

    !! 3) Hyperparameters
    !!--------------------
    hierarchy: IF (ANY(hyper(:))) THEN
      FORALL (i=1:Nparhyp,maskhypall(i)) allparhypcurr(:,:,i) = mucurr(i)
      FORALL (xOBS=1:Nx,yOBS=1:Ny,i=1:Nparhyp,maskhyp(xOBS,yOBS,i)) &
        allparhypcurr(xOBS,yOBS,i) = parmcmc(xOBS,yOBS,i2ih(i),icurr)

      !! a. Average
      mumcmc(:,icurr) = mumcmc(:,iprev)      
      average: DO ihpar=1,Nparhyp
        IF (maskhypall(ihpar)) THEN
          IF (debug) PRINT*, " - mu("//TRIMLR(parinfo(i2ih(ihpar))%name)//")"
          lim(:) = parinfo(i2ih(ihpar))%limits(:)
          mucurr(:) = mumcmc(:,icurr)
          mumcmc(ihpar,icurr) = RAND_GENERAL(lnpost_mu, lim(:), YLOG=.TRUE., &
                                             VERBOSE=debug, LNFUNC=.TRUE., &
                                             ACCURACY=accrand, NMAX=Ngibbsmax)
          ! IF (ihpar==1) PRINT*, 'mumcmc(1,icurr)', mumcmc(ihpar,icurr) ! with fixed seed

        END IF
      END DO average
      ! DO i=1,MERGE(2,1,verbose)
      !   WRITE(unitlog(i),*) 'Average sampling [done] - '// &
      !                         TRIMLR(TIMINFO(timestr))
      ! END DO

      !! b. Variance
      mucurr(:) = mumcmc(:,icurr)
      sigmcmc(:,icurr) = sigmcmc(:,iprev)
      corrcurr(:) = corrmcmc(:,iprev)
      variance: DO ihpar=1,Nparhyp
        IF (maskhypall(ihpar)) THEN
          IF (debug) PRINT*, " - sig("//TRIMLR(parinfo(i2ih(ihpar))%name)//")"
          lim(:) = cenlnS0(ihpar) + limlnS0(:)
          sigcurr(:) = sigmcmc(:,icurr)
  
          sigmcmc(ihpar,icurr) = EXP(RAND_GENERAL(lnpost_sig, lim(:), YLOG=.TRUE., &
                                                  VERBOSE=debug, LNFUNC=.TRUE., &
                                                  ACCURACY=accrand, NMAX=Ngibbsmax))
          ! IF (ihpar==1) PRINT*, 'sigmcmc(1,icurr)', sigmcmc(ihpar,icurr) ! with fixed seed

        END IF
      END DO variance
      ! DO i=1,MERGE(2,1,verbose)
      !   WRITE(unitlog(i),*) 'Variance sampling [done] - '// &
      !                         TRIMLR(TIMINFO(timestr))
      ! END DO

      !! c. Correlation
      !! Correlation coefficients are the tricky part. Indeed, the covariance
      !! matrix has to be positive definite. The Wishart prior is here to ensure
      !! that. However, if due to numerical inaccuracies, a value of rho is drawn
      !! slightly off the range that would allow a positive definite matrix, we
      !! simply ignore this drawing (i.e. we keep
      !! the previous MCMC parameter value). Of course, it is only at step
      !! iMCMC+1 that we realize that rho(iMCMC) was out, so we ignore also step
      !! iMCMC. This way we do not interfere with the MCMC statistics and avoid
      !! having a NaN introduced in the chain.
      corrmcmc(:,icurr) = corrmcmc(:,iprev) ! à tester
      IF (MODULO(counter,10)==1 .OR. counter==Nmcmc) THEN
        sigcurr(:) = sigmcmc(:,icurr)
        correlation: DO icorr=1,Ncorrhyp
          IF (maskhypall(icorr2ij(icorr,1)).AND.maskhypall(icorr2ij(icorr,2))) THEN
            IF (debug) PRINT*, " - corr("//TRIMLR(corrhypname(icorr))//")"
            corrcurr(:) = corrmcmc(:,icurr)
            !! Used by S-M, cannot wait the next big loop to update
            cov_prev(:,:) = COVARIANCE(sigmcmc(:,icurr),corrmcmc(:,icurr))
            invcov_prev(:,:) = INVERT_CHOLESKY( cov_prev(:,:), &
                                                DETERMINANT=detcov_prev, &
                                                NOPOSDEF=noposdef_prev )
            corrmcmc(icorr,icurr) = RAND_GENERAL(lnpost_corr, limR(:), YLOG=.TRUE., &
                                                 VERBOSE=debug, LNFUNC=.TRUE.,  &
                                                 ! ACCURACY=accrand, NMAX=Ngibbsmax)
                                                 ACCURACY=1.E-2_DP, NMAX=Ngibbsmax)
            ! IF (icorr==1) PRINT*, 'corrmcmc(1,icurr)', corrmcmc(icorr,icurr) ! with fixed seed
            
            IF (isNaN(corrmcmc(icorr,icurr))) THEN
              IF (icorr > 1) THEN
                corrmcmc(icorr-1:icorr,icurr) = corrmcmc(icorr-1:icorr,iprev)
              ELSE 
                corrmcmc(icorr,icurr) = corrmcmc(icorr,iprev)
                CALL WARNING("HISTOIRE", &
                             "first element in the correlation matrix was a NaN")
              END IF
            END IF
          
          END IF
        END DO correlation
        DO i=1,MERGE(2,1,verbose)
          WRITE(unitlog(i),*) 'Correlation sampling [done] - '// &
                                TRIMLR(TIMINFO(timestr))
        END DO
      END IF

    ELSE
      PRINT*, "  HISTOIRE: Non-hierarchical run."
      !! If the run is non hierarchical, we save in place of the hyperparameters
      !! the moments of the parameters
      FORALL (ipar=1:Nparhyp,maskhypall(ipar))
        mumcmc(ipar,icurr) = MEAN(parmcmc(:,:,i2ih(ipar),icurr), &
                                  MASK=maskhyp(:,:,ipar))
        sigmcmc(ipar,icurr) = SIGMA(parmcmc(:,:,i2ih(ipar),icurr), &
                                    MASK=maskhyp(:,:,ipar))
      END FORALL

      FORALL (icorr=1:Ncorrhyp, &
              maskhypall(icorr2ij(icorr,1)).AND.maskhypall(icorr2ij(icorr,2))) &
        corrmcmc(icorr,icurr) &
          = CORRELATE(parmcmc(:,:,i2ih(icorr2ij(icorr,1)),icurr),&
                      parmcmc(:,:,i2ih(icorr2ij(icorr,2)),icurr), &
                      MASK=(maskhyp(:,:,icorr2ij(icorr,1)) &
                            .AND. maskhyp(:,:,icorr2ij(icorr,2))))

    END IF hierarchy
    
    !! 4) Outputs
    !!------------
    CALL WRITE_HDF5(DBLARR4D=parmcmc(:,:,:,icurr:icurr), FILE=filMCMC, &
                    NAME=nampar, COMPRESS=compress, VERBOSE=debug, &
                    IND4=[counter,counter])
    CALL WRITE_HDF5(DBLARR4D=ln1pdmcmc(:,:,:,icurr:icurr), FILE=filMCMC, &
                    NAME=namln1pd, COMPRESS=compress, VERBOSE=debug, &
                    IND4=[counter,counter])
    CALL WRITE_HDF5(DBLARR2D=mumcmc(:,icurr:icurr), FILE=filMCMC, &
                    NAME=nammu, COMPRESS=compress, VERBOSE=debug, &
                    IND2=[counter,counter])
    CALL WRITE_HDF5(DBLARR2D=sigmcmc(:,icurr:icurr), FILE=filMCMC, &
                    NAME=namsig, COMPRESS=compress, VERBOSE=debug, &
                    IND2=[counter,counter])
    CALL WRITE_HDF5(DBLARR2D=corrmcmc(:,icurr:icurr),FILE=filMCMC(:), &
                    NAME=namcorr,COMPRESS=compress,VERBOSE=debug, &
                    IND2=[counter,counter])
    CALL WRITE_HDF5(INTARR1D=[counter],FILE=filMCMC, &
                    NAME=namind,COMPRESS=compress,VERBOSE=debug, &
                    IND1=[1,1])
    
    !! Check output for debugging
    IF (debug) THEN
      PRINT*, 'Counter = ', counter
      PRINT*, "  mu = ", mumcmc(:,icurr)
      IF (Nx*Ny > 1) THEN  
        PRINT*, "  sigma = ", sigmcmc(:,icurr)
        PRINT*, "  corr = ", corrmcmc(:,icurr)
      END IF
      IF (calib) PRINT*, "  ln1pd = ",ln1pdmcmc(:,:,:,icurr)
      WRITE(unitlog(i),*)
    END IF
    
    !! Check for NaNs and eventually stop the code
    IF (ANY(isNaN(mumcmc(:,icurr)))) &
      CALL STRIKE("HISTOIRE", &
                  "a NaN appeared in the MCMC at iteration "//TRIMLR(PRING(counter)))

    !! Normal termination condition
    IF (counter == Nmcmc) EXIT

    !! Update the indices
    CALL SWAP(icurr,iprev)
    counter = counter + 1
    
  END DO big_loop

  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'HISTOIRE: HB fit [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO
  
  !! Free memory space
  DEALLOCATE(wOBS, FnuOBS, dFnuOBS, maskint, mask, maskpar, maskhyp, &
             maskwall, maskparall, maskhypall)

END PROGRAM fit_hb
