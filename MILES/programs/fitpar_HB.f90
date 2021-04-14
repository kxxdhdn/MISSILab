!******************************************************************************
!*
!*                      Fit Simulation Parameters (HB)
!*
!******************************************************************************


  !==========================================================================
  ! 1) AUTHOR: D. HU
  ! 
  ! 2) DESCRIPTION: NONE
  !
  ! 3) HISTORY: 
  !    - 20210120: Created.
  !==========================================================================

PROGRAM fitpar_HB

  USE utilities, ONLY: DP, pring, trimLR, trimEQ, swap, tinyDP, &
                       banner_program, ustd, isNaN, strike, warning, &
                       timinfo, initiate_clock, time_type, today
  USE arrays, ONLY: iwhere, closest, reallocate
  USE constants, ONLY: MKS
  USE matrices, ONLY: invert_cholesky
  USE statistics, ONLY: mean, sigma, correlate, correl_index, corr2Rmat
  USE inout, ONLY: read_hdf5, write_hdf5, check_hdf5, h5ext, ascext, lenpar, lenpath
  USE grain_optics, ONLY: lendustQ
  USE random, ONLY: generate_newseed, rand_general
  USE auxil, ONLY: read_master, specModel, initparam
  USE HB_kit, ONLY: Nx, Ny, NwOBS, wOBS, nuOBS, xOBS, yOBS, &
                    FnuOBS, dFnuOBS, cenlnS0, siglnS0, &
                    i2ih, ipar, ihpar, icorr, icorr2ij, &
                    parcurr, parhypcurr, allparhypcurr, &
                    cov_prev, invcov_prev, invcorr_prev, &
                    mucurr, sigcurr, corrcurr, Nparhyp, Ncorrhyp, &
                    lnpost_par, lnlhobs_par, detcorr, noposdef_prev, &
                    lnpost_mu, lnpost_sig, lnpost_corr, covariance, &
                    mask, maskpar, maskhyp, maskhypcurr, maskS, &!maskextra, &
                    ind, parinfo, parhypinfo, Qabs, extinct
  IMPLICIT NONE

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: textwid = 60
  INTEGER, PARAMETER :: Ngibbsmax = 3000 ! determines Ngrid
  REAL(DP), PARAMETER :: accrand = 1.E-3_DP
  CHARACTER(*), PARAMETER :: nampar = "Parameter values"
  ! CHARACTER(*), PARAMETER :: namln1pd = "ln(1+delta)"
  CHARACTER(*), PARAMETER :: nammu = "Mean of hyperdistribution"
  CHARACTER(*), PARAMETER :: namsig = "Sigma of hyperdistribution"
  CHARACTER(*), PARAMETER :: namcorr = "Correlation of hyperdistribution"
  CHARACTER(*), PARAMETER :: dirIN = '../out1/'
  CHARACTER(*), PARAMETER :: filOBS = dirIN//'galspec'//h5ext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Input variables
  INTEGER :: x, y, i, ih, indhyp, indres0!, iw, j
  INTEGER :: Npar, Nmcmc, NiniMC, Nsou, Nparfree!, Nwfree
  INTEGER :: Ncont, Nband, Nline, Nextc, Nstar, Nextra
  INTEGER, DIMENSION(:), ALLOCATABLE :: indresume, itied
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: maskint ! convert mask=0 to mask=T
  LOGICAL :: verbose, calib, newseed, newinit, dostop, resume, nohi
  CHARACTER(lenpath) :: dirOUT, filOUT, filMCMC, fiLOG
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB, labL

  !! MCMC parameters
  INTEGER :: icurr, iprev, counter
  REAL(DP), DIMENSION(2) :: lim, limlnS0, limR
  REAL(DP), DIMENSION(:), ALLOCATABLE :: mu0, sig0, corr0
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meanmu, stdevmu, meansig, stdevsig
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meancorr, stdevcorr
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dblarr2d
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mumcmc, sigmcmc, corrmcmc
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parini, parmcmc
  LOGICAL :: fileOK

  !! Output variables
  INTEGER :: Ncorr
  INTEGER, DIMENSION(2) :: unitlog
  TYPE(time_type) :: timestr

  !! Analysis variables
  INTEGER :: t_burnin, t_end, Nmcmc_eff
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: meanpar, stdevpar
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT, &
    FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FnuCONT_tab, &
    FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parpix, par4D

  
  !!------------------------------------------------------------------------
  !!                            I. Read the inputs
  !!------------------------------------------------------------------------

  !! 1) Main settings
  !!------------------
  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='Wavelength (microns)', N1=NwOBS)
  CALL READ_MASTER(WAVALL=wOBS(:), DIRIN=dirIN, DIROUT=dirOUT, &
                   VERBOSE=verbose, NMCMC=Nmcmc, NINIMC=NiniMC, &
                   CALIB=calib, NEWSEED=newseed, NEWINIT=newinit, &
                   LABL=labL, LABB=labB, QABS=Qabs, EXTINCT=extinct, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, &
                   NEXTC=Nextc, NSTAR=Nstar, NEXTRA=Nextra, NCORR=Ncorr, &
                   DOSTOP=dostop, RESUME=resume, INDRESUME=indres0, NOHI=nohi, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit, &
                   PARHYPINFO=parhypinfo, NPARHYP=Nparhyp, NCORRHYP=Ncorrhyp)

  !! Output settings
  !!-----------------
  IF (nohi) THEN
    filOUT = TRIMLR(dirOUT)//'fitpar_BB'//h5ext
    filMCMC = TRIMLR(dirOUT)//'parlog_fitpar_BB'//h5ext
    fiLOG = TRIMLR(dirOUT)//'log_fitpar_BB'//ascext
  ELSE
    filOUT = TRIMLR(dirOUT)//'fitpar_HB'//h5ext
    filMCMC = TRIMLR(dirOUT)//'parlog_fitpar_HB'//h5ext
    fiLOG = TRIMLR(dirOUT)//'log_fitpar_HB'//ascext
  END IF
  
  CALL INITIATE_CLOCK(timestr)
  OPEN (ulog,FILE=fiLOG,STATUS=MERGE("OLD    ","REPLACE",resume), &
        ACTION="WRITE",POSITION=MERGE("APPEND","REWIND",resume))
  unitlog(:) = [ ulog, ustd ]
  
  IF (newseed) CALL GENERATE_NEWSEED()
  IF (verbose) PRINT*
  DO i=1,MERGE(2,1,verbose)
    IF (.NOT. resume) THEN
      CALL BANNER_PROGRAM("HIBARI: HIerarchical BAyesian fitting Routine of mid-IR emission", &
                          UNIT=unitlog(i), SWING=.True.)
    ELSE
      WRITE(unitlog(i),*) 
      WRITE(unitlog(i),*) REPEAT("*",80)
      WRITE(unitlog(i),*) 
      WRITE(unitlog(i),*) "...RESUMING WHERE THE RUN STOPPED!!!"
      WRITE(unitlog(i),*) 
      WRITE(unitlog(i),"(A80)") ADJUSTR("("//TRIMLR(TODAY())//")")
      WRITE(unitlog(i),*) 
    END IF
  END DO
  
  !! 2) Observation sample
  !!-----------------------
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
  ALLOCATE(mask(Nx,Ny,NwOBS))
  mask(:,:,:) = ( maskint(:,:,:) == 0 )
  WHERE (isNaN(FnuOBS(:,:,:)) .OR. .NOT. mask(:,:,:)) FnuOBS(:,:,:) = 0._DP

  !! Noise
  CALL READ_HDF5(DBLARR3D=dFnuOBS, FILE=filOBS, &
                 NAME='dFnuOBS ('//TRIMLR(spec_unit)//')')
  WHERE (isNaN(dFnuOBS(:,:,:)))
    dFnuOBS(:,:,:) = 0._DP
    mask(:,:,:) = .FALSE.
  END WHERE
  
  ALLOCATE(nuOBS(NwOBS))
  nuOBS(:) = MKS%clight/MKS%micron / wOBS(:)

  !! Read the instrumental covariance matrix (TBD)


  !! Indices of tied parameters
  ALLOCATE(itied(Npar))
  itied(:) = 0
  DO i=1,Npar
    IF (.NOT. TRIMEQ(parinfo(i)%tied,"")) &
      CALL IWHERE(TRIMEQ(parinfo(:)%name,parinfo(i)%tied),itied(i))
  END DO

  !! Number of free parameters
  Nparfree = COUNT( (.NOT. parinfo(:)%fixed) .AND. (itied(:) == 0) )
  ! CALL IWHERE(( .NOT. parinfo(:)%fixed ) .AND. ( itied(:) == 0 ),ifree)
  ! Ncons = NwOBS - Nfreefilt
  ! Ndof = Ncons - Nparfree

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
                          //TRIMLR(PRING(NwOBS))
      WRITE(unitlog(i),*) " - Nparhyp = "//TRIMLR(PRING(Nparhyp))
      WRITE(unitlog(i),*) " - Ncorrhyp = "//TRIMLR(PRING(Ncorrhyp))
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
    ! IF (parinfo(i2ih(i))%model) THEN
      maskhyp(:,:,i) = ANY(mask(:,:,:),DIM=3)
    ! ELSE
      ! maskhyp(:,:,i) = maskextra(:,:,i2ih(i)-iext0+1)
    ! END IF
  END DO

  !! 3D parameter mask
  ALLOCATE (maskpar(Nx,Ny,Npar))
  DO i=1,Npar
    ! IF (parinfo(i)%model) THEN
      maskpar(:,:,i) = ANY(mask(:,:,:),DIM=3)
    ! ELSE
    !   maskpar(:,:,i) = maskextra(:,:,i-iext0+1)
    ! END IF
  END DO  
  
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'Read the inputs [done]'
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
                 LABB=labB(:), LABL=labL(:), QABS=Qabs(:))

  !! c. Range of parameters for Gibbs sampling
  !!-------------------------------------------
  !! Calibration errors
  
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
    DO i=1,Nparhyp
      mu0(i) = MEAN(PACK(parini(:,:,i2ih(i),1),maskhyp(:,:,i)))
      sig0(i) = SIGMA(PACK(parini(:,:,i2ih(i),1),maskhyp(:,:,i)))
    END DO
    WHERE (ABS(sig0(:)) < 1.E-5_DP*ABS(mu0(:)/sig0(:))) &
      sig0(:) = MERGE( ABS(mu0(:)), 1._DP, (mu0(:) > 0._DP) )
    WHERE (sig0(:) <= 0._DP) &
      sig0(:) = 0.5_DP * ( parinfo(i2ih(:))%limits(2) &
                         - parinfo(i2ih(:))%limits(1) )
    corr0(:) = 0._DP
  ELSE
    ALLOCATE (mu0(Npar),sig0(Npar),corr0(Ncorr))
    DO i=1,Npar
      mu0(i) = MEAN(PACK(parini(:,:,i,1),maskpar(:,:,i)))
      sig0(i) = SIGMA(PACK(parini(:,:,i,1),maskpar(:,:,i)))
    END DO
    WHERE (ABS(sig0(:)) < tinyDP) &
      sig0(:) = MERGE( ABS(mu0(:)), 1._DP, (mu0(:) > 0._DP) )
    WHERE (sig0(:) <= 0._DP) &
      sig0(:) = 0.5_DP * ( parinfo(:)%limits(2) - parinfo(:)%limits(1) )
    corr0(:) = 0._DP
  END IF

  !! Print the initial values
  DO i=1,MERGE(2,1,debug)
    WRITE(unitlog(i),*) "mu0 = ", mu0
    WRITE(unitlog(i),*) "sig0 = ", sig0
    WRITE(unitlog(i),*) "corr0 = ", corr0
    ! WRITE(unitlog(i),*) "ln1pd = ", ln1pd0
  END DO
  
  !! Initialize the parameters of the chain
  ALLOCATE (parmcmc(Nx,Ny,Npar,2))!,ln1pdmcmc(Ncalib,2))
  DO i=1,2 
    FORALL (x=1:Nx,y=1:Ny,ANY(maskpar(x,y,:))) parmcmc(x,y,:,i) = parini(x,y,:,1)
    ! ln1pdmcmc(:,i) = ln1pd0(:)
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
  ! ALLOCATE (ln1pd(Ncalib)) ! LN(1+delta) for the current iteration
  ! ALLOCATE (delp1(NwOBS))  ! delta+1 for the current iteration
  ALLOCATE (parcurr(Npar)) ! current parameters
  IF (ANY(parinfo(:)%hyper)) THEN
    ALLOCATE (allparhypcurr(Nx,Ny,Nparhyp)) ! parameter grid
    ALLOCATE (cov_prev(Nparhyp,Nparhyp))    ! correlation matrix
    ALLOCATE (invcov_prev(Nparhyp,Nparhyp)) ! inverse covariance matrix
    ALLOCATE (invcorr_prev(Nparhyp,Nparhyp)) ! inverse correlation matrix
    ALLOCATE (mucurr(Nparhyp))              ! current average
    ALLOCATE (sigcurr(Nparhyp))             ! current average
    ALLOCATE (corrcurr(Ncorrhyp))           ! current correlation coefficient
    ALLOCATE (parhypcurr(Nparhyp))          ! current parameters
    ALLOCATE (maskhypcurr(Nparhyp))         ! current parameters' mask
  END IF

  !! Open the output file
  fileresume: IF (.NOT. resume) THEN

    !! Initialize the files
    CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME="Parameter label", &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
    CALL WRITE_HDF5(STRARR1D=parhypinfo(:)%name,NAME="Hyper parameter label",&
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(INTARR1D=[Nmcmc], NAME="Length of MCMC", &
                    FILE=filMCMC, COMPRESS=compress, verbose=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(INITDBLARR=[Nx,Ny,Npar,Nmcmc], NAME=nampar, &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(INITDBLARR=[Nparhyp,Nmcmc], NAME=nammu, &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(INITDBLARR=[Nparhyp,Nmcmc], NAME=namsig, &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(INITDBLARR=[Ncorrhyp,Nmcmc], NAME=namcorr, &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(INITINTARR=[1], NAME="Last index", &
                    FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  ELSE

    !! Resume the run from one index before the stop (Last index)
    fileOK = CHECK_HDF5(FILE=filMCMC)
    IF (fileOK) THEN
      !! Auto-resume
      CALL READ_HDF5(INTARR1D=indresume,FILE=filMCMC,NAME="Last index")
      
      IF (indres0 > 0) THEN
        IF (indres0 > indresume(1)) THEN
          CALL WARNING( 'FITPAR_HB', &
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
    
    !! 2) Physical parameters
    !!------------------------

    !! Covariance matrix of the hyperparameters
    IF (ANY(parinfo(:)%hyper)) THEN
      cov_prev(:,:) = COVARIANCE(sigmcmc(:,iprev),corrmcmc(:,iprev))
      invcov_prev(:,:) = INVERT_CHOLESKY(cov_prev(:,:))
      invcorr_prev(:,:) = INVERT_CHOLESKY( &
                            corr2Rmat(corrmcmc(:,iprev),Nparhyp), &
                            DETERMINANT=detcorr, NOPOSDEF=noposdef_prev ) ! new
      mucurr(:) = mumcmc(:,iprev)
    END IF

    !! Draw the model parameters
    parmcmc(:,:,:,icurr) = parmcmc(:,:,:,iprev)
    param: DO ipar=1,Npar
      IF (.NOT. parinfo(ipar)%fixed) THEN
        IF (debug) PRINT*, " - "//TRIMLR(parinfo(ipar)%name) 

        !! Even the parinfo-not-limited par need this limits
        !! See c. Range of parameters for Gibbs sampling
        lim(:) = parinfo(ipar)%limits(:)
        xsource: DO xOBS=1,Nx
          ysource: DO yOBS=1,Ny
            IF (maskpar(xOBS,yOBS,ipar)) THEN
              parcurr(:) = parmcmc(xOBS,yOBS,:,icurr)
              IF (parinfo(ipar)%hyper) THEN
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
! DO i=1,MERGE(2,1,verbose)
!   WRITE(unitlog(i),*) 'Parameter sampling [done] - '// &
!                         TRIMLR(TIMINFO(timestr))
! END DO

    !! 3) Hyperparameters
    !!--------------------
    hierarchy: IF (ANY(parinfo(:)%hyper)) THEN

      !! a. Average
      mumcmc(:,icurr) = mumcmc(:,iprev)      
      average: DO ihpar=1,Nparhyp    
        IF (debug) PRINT*, " - mu("//TRIMLR(parinfo(i2ih(ihpar))%name)//")"
        lim(:) = parinfo(i2ih(ihpar))%limits(:)
        mucurr(:) = mumcmc(:,icurr)
        mumcmc(ihpar,icurr) = RAND_GENERAL(lnpost_mu,lim(:),YLOG=.True., &
                                           VERBOSE=debug,LNFUNC=.True., &
                                           ACCURACY=accrand,NMAX=Ngibbsmax)
      END DO average
DO i=1,MERGE(2,1,verbose)
  WRITE(unitlog(i),*) 'Average sampling [done] - '// &
                        TRIMLR(TIMINFO(timestr))
END DO

      !! b. Variance
      mucurr(:) = mumcmc(:,icurr)
      sigmcmc(:,icurr) = sigmcmc(:,iprev)
      corrcurr(:) = corrmcmc(:,iprev)
      variance: DO ihpar=1,Nparhyp
        IF (debug) PRINT*, " - sig("//TRIMLR(parinfo(i2ih(ihpar))%name)//")"
        lim(:) = cenlnS0(ihpar) + limlnS0(:)
        sigcurr(:) = sigmcmc(:,icurr)
        sigmcmc(ihpar,icurr) = EXP(RAND_GENERAL(lnpost_sig,lim(:),YLOG=.True., &
                                                VERBOSE=debug,LNFUNC=.True., &
                                                ACCURACY=accrand,NMAX=Ngibbsmax))
      END DO variance
DO i=1,MERGE(2,1,verbose)
  WRITE(unitlog(i),*) 'Variance sampling [done] - '// &
                        TRIMLR(TIMINFO(timestr))
END DO

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
      sigcurr(:) = sigmcmc(:,icurr)
      corrmcmc(:,icurr) = corrmcmc(:,iprev)
      correlation: DO icorr=1,Ncorrhyp
        ! IF (debug) PRINT*, " - corr("//TRIMLR(corrhypname(icorr))//")"
        corrcurr(:) = corrmcmc(:,icurr)
        corrmcmc(icorr,icurr) = RAND_GENERAL(lnpost_corr,limR(:),YLOG=.True., &
                                             VERBOSE=debug,LNFUNC=.True.,  &
                                             ! ACCURACY=accrand,NMAX=Ngibbsmax)
                                             ACCURACY=1.E-2_DP,NMAX=Ngibbsmax)
        IF (isNaN(corrmcmc(icorr,icurr))) THEN
          IF (icorr > 1) THEN
            corrmcmc(icorr-1:icorr,icurr) = corrmcmc(icorr-1:icorr,iprev)
          ELSE 
            corrmcmc(icorr,icurr) = corrmcmc(icorr,iprev)
            CALL WARNING("FITPAR_HB", &
                         "first element in the correlation matrix was a NaN")
          END IF
        END IF
      END DO correlation
DO i=1,MERGE(2,1,verbose)
  WRITE(unitlog(i),*) 'Correlation sampling [done] - '// &
                        TRIMLR(TIMINFO(timestr))
END DO

    ELSE
      PRINT*, "  HIBARI: Non-hierarchical run."
      !! If the run is non hierarchical, we save in place of the hyperparameters
      !! the moments of the parameters
      FORALL (ipar=1:Nparhyp)
        mumcmc(ipar,icurr) = MEAN(parmcmc(:,:,i2ih(ipar),icurr), &
                                  MASK=maskhyp(:,:,ipar))
        sigmcmc(ipar,icurr) = SIGMA(parmcmc(:,:,i2ih(ipar),icurr), &
                                    MASK=maskhyp(:,:,ipar))
      END FORALL
      FORALL (icorr=1:Ncorrhyp) &
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
                    NAME="Last index",COMPRESS=compress,VERBOSE=debug, &
                    IND1=[1,1])
    
    !! Check output for debugging
    IF (debug) THEN
      PRINT*, 'Counter = ', counter
      PRINT*, "  mu = ", mumcmc(:,icurr)
      IF (Nx*Ny > 1) THEN  
        PRINT*, "  sigma = ", sigmcmc(:,icurr)
        PRINT*, "  corr = ", corrmcmc(:,icurr)
      END IF
      ! IF (calib) PRINT*, "  ln1pd = ",ln1pdmcmc(:,icurr)
      WRITE(unitlog(i),*)
    END IF
    
    !! Check for NaNs and eventually stop the code
    IF (ANY(isNaN(mumcmc(:,icurr)))) &
      CALL STRIKE("FITPAR_HB", &
                  "a NaN appeared in the MCMC at iteration "//TRIMLR(PRING(counter)))

    !! Normal termination condition
    IF (counter == Nmcmc) EXIT

    !! Update the indices
    CALL SWAP(icurr,iprev)
    counter = counter + 1
    
  END DO big_loop

  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'HB fit galspec (simulated spectra) [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO
  
  !!------------------------------------------------------------------------
  !!                                 Analysis
  !!------------------------------------------------------------------------

  !! Compute the statistics
  !!------------------------
  t_burnin = 2000000
  t_end = Nmcmc ! end of Markov chain Monte Carlo
  IF (t_burnin <= 0 .OR. t_burnin >= t_end) t_burnin = t_end/10
  Nmcmc_eff = t_end - t_burnin + 1

  ALLOCATE(parpix(Nx,Ny,Npar,t_end), meanpar(Nx,Ny,Npar), stdevpar(Nx,Ny,Npar))

  CALL READ_HDF5(DBLARR4D=par4D, FILE=filMCMC, &
                 NAME=nampar, IND4=[1,t_end])
  parpix(:,:,:,:) = par4D(:,:,:,:)
  meanpar(:,:,:) = MEAN(parpix(:,:,:,t_burnin:t_end), DIM=4)
  stdevpar(:,:,:) = SIGMA(parpix(:,:,:,t_burnin:t_end), DIM=4)

  CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME="Parameter label", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
  CALL WRITE_HDF5(DBLARR3D=meanpar(:,:,:), NAME="Mean of parameter value", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=stdevpar(:,:,:), NAME="Sigma of parameter value", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

  !! Compute the average, standard deviations and correlations of each quantity
  !!----------------------------------------------------------------------------
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*) " - Compute the statistics ("//TRIMLR(TIMINFO(timestr)) &
                        //")"
    WRITE(unitlog(i),*)
  END DO

  ALLOCATE (meanmu(Nparhyp),stdevmu(Nparhyp),meansig(Nparhyp),stdevsig(Nparhyp))

  CALL READ_HDF5(DBLARR2D=dblarr2d, FILE=filMCMC, &
                 NAME=nammu, IND2=[1,t_end])
  meanmu(:) = MEAN(dblarr2d(:,t_burnin:t_end), DIM=2)
  stdevmu(:) = SIGMA(dblarr2d(:,t_burnin:t_end), DIM=2)
  CALL READ_HDF5(DBLARR2D=dblarr2d, FILE=filMCMC, &
                 NAME=namsig, IND2=[1,t_end])
  meansig(:) = MEAN(dblarr2d(:,t_burnin:t_end), DIM=2)
  stdevsig(:) = SIGMA(dblarr2d(:,t_burnin:t_end), DIM=2)
  
  ALLOCATE (meancorr(Ncorrhyp),stdevcorr(Ncorrhyp))

  CALL READ_HDF5(DBLARR2D=dblarr2d, FILE=filMCMC, &
                 NAME=namcorr, IND2=[1,t_end])
  meancorr(:) = MEAN(dblarr2d(:,t_burnin:t_end), DIM=2)
  stdevcorr(:) = SIGMA(dblarr2d(:,t_burnin:t_end), DIM=2)
  
  CALL WRITE_HDF5(DBLARR1D=meanmu(:),NAME="Mean of parameter average", &
                  FILE=filOUT,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)
  CALL WRITE_HDF5(DBLARR1D=stdevmu(:),NAME="Sigma of parameter average", &
                  FILE=filOUT,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)
  CALL WRITE_HDF5(DBLARR1D=meansig(:), &
                  NAME="Mean of parameter standard deviation", &
                  FILE=filOUT,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)
  CALL WRITE_HDF5(DBLARR1D=stdevsig(:), &
                  NAME="Sigma of parameter standard deviation", &
                  FILE=filOUT,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)
  CALL WRITE_HDF5(DBLARR1D=meancorr(:),NAME="Mean of parameter correlation", &
                  FILE=filOUT,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)
  CALL WRITE_HDF5(DBLARR1D=stdevcorr(:),NAME="Sigma of parameter correlation", &
                  FILE=filOUT,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)

  !! Calculate model
  !!-----------------
  ALLOCATE(FnuMOD(Nx,Ny,NwOBS))

  FnuMOD(:,:,:) = specModel( wOBS(:), INDPAR=ind, PARVAL=meanpar(:,:,:), &
                             QABS=Qabs(:), EXTINCT=extinct(:,:), &
                             FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
                             PABS=Pabs, FNULINE=FnuLINE, &
                             FNUCONT_TAB=FnuCONT_tab, FNUBAND_TAB=FnuBAND_tab, &
                             FNUSTAR_TAB=FnuSTAR_tab, PABS_TAB=Pabs_tab, &
                             FNULINE_TAB=FnuLINE_tab )

  CALL WRITE_HDF5(DBLARR4D=FnuCONT_tab, NAME='FnuCONT ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR4D=FnuLINE_tab, NAME='FnuLINE ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR4D=FnuBAND_tab, NAME='FnuBAND ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR4D=FnuSTAR_tab, NAME='FnuSTAR ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR4D=Pabs_tab, NAME='PABS', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=FnuMOD, NAME='FnuMOD ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'fitpar_HB analysis [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO

  !! Final
  !!-------
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) "PROGRAM EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
    WRITE(unitlog(i),*)
  END DO
  
  !! Free memory space
  DEALLOCATE(wOBS, FnuOBS, dFnuOBS, FnuMOD, maskint, mask)

END PROGRAM fitpar_HB
