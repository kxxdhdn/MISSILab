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

  USE utilities, ONLY: DP, pring, trimLR, trimeq, swap, timinfo, tinyDP, &
                       banner_program, ustd, isNaN, initiate_clock, time_type
  USE arrays, ONLY: iwhere, closest
  USE constants, ONLY: MKS
  USE matrices, ONLY: invert_cholesky
  USE statistics, ONLY: mean, sigma
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, ascext, lenpar, lenpath
  USE grain_optics, ONLY: lendustQ
  USE random, ONLY: generate_newseed, rand_general
  USE auxil, ONLY: read_master, specModel, initparam
  USE HB_kit, ONLY: Nx, Ny, NwOBS, wOBS, nuOBS, xOBS, yOBS, &
                    FnuOBS, dFnuOBS, cenlnS0, siglnS0, &
                    i2ih, ipar, ihpar, &!icorr, icorr2ij, &
                    parcurr, parhypcurr, allparhypcurr, &
                    cov_prev, invcov_prev, &
                    mucurr, sigcurr, corrcurr, Nparhyp, Ncorrhyp, &
                    lnpost_par, lnlhobs_par, &
                    lnpost_mu, lnpost_sig, lnpost_corr, covariance, &
                    mask, maskpar, maskhyp, maskhypcurr, maskS, &!maskextra, &
                    ind, parinfo, parhypinfo, Qabs, extinct
  IMPLICIT NONE

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: textwid = 60
  INTEGER, PARAMETER :: Ngibbsmax = 3000
  REAL(DP), PARAMETER :: accrand = 1.E-3_DP
  CHARACTER(*), PARAMETER :: nampar = "Parameter values"
  CHARACTER(*), PARAMETER :: namln1pd = "ln(1+delta)"
  CHARACTER(*), PARAMETER :: nammu = "Mean of hyperdistribution"
  CHARACTER(*), PARAMETER :: namsig = "Sigma of hyperdistribution"
  CHARACTER(*), PARAMETER :: namcorr = "Correlation of hyperdistribution"
  CHARACTER(*), PARAMETER :: namrat = "Band ratio values"
  CHARACTER(*), PARAMETER :: dirIN = '../cache/'
  CHARACTER(*), PARAMETER :: dirOUT = '../out1/'
  CHARACTER(*), PARAMETER :: filOBS = dirIN//'galspec'//h5ext
  CHARACTER(*), PARAMETER :: filMCMC = dirIN//'parlog_fitpar_HB'//h5ext
  CHARACTER(*), PARAMETER :: filOUT = dirOUT//'fitpar_HB'//h5ext
  CHARACTER(*), PARAMETER :: fiLOG = dirOUT//'log_fitpar_HB'//ascext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Input variables
  INTEGER :: x, y, i, ih, indhyp!, iw, j
  INTEGER :: Npar, Nmcmc, NiniMC, Nsou, Nrat, Nparfree!, Nwfree
  INTEGER :: Ncont, Nband, Nline, Npabs, Nstar, Nextra
  INTEGER, DIMENSION(:), ALLOCATABLE :: itied, indBcorr
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: maskint ! convert mask=0 to mask=T
  LOGICAL :: verbose, calib, newseed, newinit, dostop
  CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB, labL, ratname, namBcorr
  CHARACTER(lenpar) :: spec_unit

  !! MCMC parameters
  INTEGER :: icurr, iprev, counter
  REAL(DP), DIMENSION(2) :: lim, limlnS0, limR
  REAL(DP), DIMENSION(:), ALLOCATABLE :: mu0, sig0, corr0
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mumcmc, sigmcmc, corrmcmc
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parini, parmcmc, ratmcmc

  !! Output variables
  INTEGER :: Ncorr
  INTEGER, DIMENSION(2) :: unitlog
  TYPE(time_type) :: timestr

  !! Analysis variables
  INTEGER :: t_burnin, t_end, Nmcmc_eff
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: meanpar, stdevpar
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: meanrat, stdevrat
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT, &
    FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: par4D, parpix, ratpix
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FnuCONT_tab, &
    FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
  
  ! Output settings
  !----------------
  CALL INITIATE_CLOCK(timestr)
  OPEN (ulog,FILE=fiLOG,STATUS="REPLACE",ACTION="WRITE")
  unitlog(:) = [ ulog, ustd ]
  
  !!------------------------------------------------------------------------
  !!                            I. Read the inputs
  !!------------------------------------------------------------------------

  !! 1) Main settings
  !!------------------
  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='Wavelength (microns)', N1=NwOBS)
  CALL READ_MASTER(WAVALL=wOBS(:), DIRIN=dirIN, &
                   VERBOSE=verbose, NMCMC=Nmcmc, NINIMC=NiniMC, &
                   CALIB=calib, NEWSEED=newseed, NEWINIT=newinit, &
                   LABQ=labQ, LABL=labL, LABB=labB, QABS=Qabs, EXTINCT=extinct, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, &
                   NPABS=Npabs, NSTAR=Nstar, NEXTRA=Nextra, &
                   DOSTOP=dostop, NCORR=Ncorr, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit, &
                   PARHYPINFO=parhypinfo, NPARHYP=Nparhyp, NCORRHYP=Ncorrhyp)

  !! Band ratio correlations
  CALL READ_HDF5(STRARR1D=ratname, FILE=filOBS, NAME='Band ratio name', N1=Nrat)
  CALL READ_HDF5(STRARR1D=namBcorr, FILE=filOBS, NAME='Correlated band name')
  CALL READ_HDF5(INTARR1D=indBcorr, FILE=filOBS, NAME='Correlated band indpar')
  
  IF (newseed) CALL GENERATE_NEWSEED()
  IF (verbose) PRINT*
  DO i=1,MERGE(2,1,verbose)
    CALL BANNER_PROGRAM("HIBARI: HIerarchical BAyesian fitting Routine of mid-IR emission", &
                        UNIT=unitlog(i), SWING=.True.)
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
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) " - Number of spectra/pixels to be fitted = " &
                        //TRIMLR(PRING(Nsou))
    WRITE(unitlog(i),*) " - Size of spectral sampling = " &
                        //TRIMLR(PRING(NwOBS))
    WRITE(unitlog(i),*)
  END DO

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
print*, 'cband init= ', parini(:,:,indBcorr(1)+1,1)
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
  ALLOCATE (parmcmc(Nx,Ny,Npar,2), ratmcmc(Nx,Ny,Nrat,2))!,ln1pdmcmc(Ncalib,2))
  DO i=1,2 
    FORALL (x=1:Nx,y=1:Ny,ANY(maskpar(x,y,:)))
      parmcmc(x,y,:,i) = parini(x,y,:,1)
      ! ln1pdmcmc(:,i) = ln1pd0(:)
    END FORALL
    ratmcmc(:,:,:,i) = 0._DP
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
    ALLOCATE (invcov_prev(Nparhyp,Nparhyp)) ! inverse correlation matrix
    ALLOCATE (mucurr(Nparhyp))              ! current average
    ALLOCATE (sigcurr(Nparhyp))             ! current average
    ALLOCATE (corrcurr(Ncorrhyp))           ! current correlation coefficient
    ALLOCATE (parhypcurr(Nparhyp))          ! current parameters
    ALLOCATE (maskhypcurr(Nparhyp))         ! current parameters' mask
  END IF

  !! Initialize the files
  CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME="Parameter label", &
                  FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
  CALL WRITE_HDF5(STRARR1D=ratname(:), NAME="Band ratio label", &
                  FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(STRARR1D=namBcorr(:), NAME="Correlated band label", &
                  FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(INTARR1D=[Nmcmc], NAME="Length of MCMC", &
                  FILE=filMCMC, COMPRESS=compress, verbose=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(INITDBLARR=[Nx,Ny,Npar,Nmcmc], NAME=nampar, &
                  FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(INITDBLARR=[Nx,Ny,Nrat,Nmcmc], NAME=namrat, &
                  FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(STRARR1D=parhypinfo(:)%name,NAME="Hyper parameter label",&
                  FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  ! CALL WRITE_HDF5(INITDBLARR=[Nparhyp,Nmcmc], NAME=nammu, &
  !                 FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  ! CALL WRITE_HDF5(INITDBLARR=[Nparhyp,Nmcmc], NAME=namsig, &
  !                 FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'Initial parameters [done]'
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
    ! IF (ANY(parinfo(:)%hyper)) THEN
    !   cov_prev(:,:) = COVARIANCE(sigmcmc(:,iprev),corrmcmc(:,iprev))
    !   invcov_prev(:,:) = INVERT_CHOLESKY(cov_prev(:,:))
    !   mucurr(:) = mumcmc(:,iprev)
    ! END IF

    !! Draw the model parameters
    parmcmc(:,:,:,icurr) = parmcmc(:,:,:,iprev)
    param: DO ipar=1,Npar
      IF (.NOT. parinfo(ipar)%fixed) THEN
        ! DO i=1,MERGE(2,1,debug)
        !   WRITE(unitlog(i),*) " - "//TRIMLR(parinfo(ipar)%name)
        ! END DO
        !! Even the parinfo-not-limited par need this limits
        !! See c. Range of parameters for Gibbs sampling
        lim(:) = parinfo(ipar)%limits(:)
        xsource: DO xOBS=1,Nx
          ysource: DO yOBS=1,Ny
            IF (maskpar(xOBS,yOBS,ipar)) THEN
              parcurr(:) = parmcmc(xOBS,yOBS,:,icurr)
              ! IF (parinfo(ipar)%hyper) THEN
              !   maskhypcurr(:) = maskhyp(xOBS,yOBS,:)
              !   WHERE (maskhyp(xOBS,yOBS,:))
              !     parhypcurr(:) = parmcmc(xOBS,yOBS,i2ih(:),icurr)
              !   ELSEWHERE
              !     parhypcurr(:) = mucurr(:)
              !   END WHERE
              ! ENDIF
              parmcmc(xOBS,yOBS,ipar,icurr) &
                = RAND_GENERAL(lnpost_par, lim(:), YLOG=.TRUE., VERBOSE=debug, &
                               LNFUNC=.TRUE., ACCURACY=accrand, NMAX=Ngibbsmax)

              ! DO i=1,MERGE(2,1,debug)
              !   WRITE(unitlog(i),*) "par(",xOBS,yOBS,ipar,") = ", &
              !                       parmcmc(xOBS,yOBS,ipar,icurr)
              ! END DO

            END IF
          END DO ysource
        END DO xsource
      END IF
    END DO param
print*, 'cband mcmc icurr = ', parmcmc(:,:,indBcorr(1)+1,icurr)

    ratmcmc(:,:,1,icurr) = EXP(parmcmc(:,:,indBcorr(1),icurr)) / EXP(parmcmc(:,:,indBcorr(4),icurr))
    ratmcmc(:,:,2,icurr) = EXP(parmcmc(:,:,indBcorr(2),icurr)) / EXP(parmcmc(:,:,indBcorr(4),icurr))
    ratmcmc(:,:,3,icurr) = EXP(parmcmc(:,:,indBcorr(3),icurr)) / EXP(parmcmc(:,:,indBcorr(4),icurr))
    ratmcmc(:,:,4,icurr) = EXP(parmcmc(:,:,indBcorr(5),icurr)) / EXP(parmcmc(:,:,indBcorr(4),icurr))

    !! 3) Hyperparameters
    !!--------------------
    ! hierarchy: IF (ANY(parinfo(:)%hyper)) THEN

    !   !! a. Average
    !   !! b. Variance
    !   !! c. Correlation
      
    ! ELSE
      
    !   !! If the run is non hierarchical, we save in place of the hyperparameters
    !   !! the moments of the parameters
    !   FORALL (ipar=1:Nparhyp)
    !     mumcmc(ipar,icurr) = MEAN(parmcmc(:,:,i2ih(ipar),icurr), &
    !                               MASK=maskhyp(:,:,ipar))
    !     sigmcmc(ipar,icurr) = SIGMA(parmcmc(:,:,i2ih(ipar),icurr), &
    !                                 MASK=maskhyp(:,:,ipar))
    !   END FORALL
    !   FORALL (icorr=1:Ncorrhyp) &
    !     corrmcmc(icorr,icurr) &
    !       = CORRELATE(parmcmc(:,:,i2ih(icorr2ij(icorr,1)),icurr),&
    !                   parmcmc(:,:,i2ih(icorr2ij(icorr,2)),icurr), &
    !                   MASK=(maskhyp(:,:,icorr2ij(icorr,1)) &
    !                         .AND. maskhyp(:,:,icorr2ij(icorr,2))))

    ! END IF hierarchy
    
    !! 4) Outputs
    !!------------
    CALL WRITE_HDF5(DBLARR4D=parmcmc(:,:,:,icurr:icurr), FILE=filMCMC, &
                    NAME=nampar, COMPRESS=compress, VERBOSE=debug, &
                    IND4=[counter,counter])
    CALL WRITE_HDF5(DBLARR4D=ratmcmc(:,:,:,icurr:icurr), FILE=filMCMC, &
                    NAME=namrat, COMPRESS=compress, VERBOSE=debug, &
                    IND4=[counter,counter])
    ! CALL WRITE_HDF5(DBLARR2D=mumcmc(:,icurr:icurr), FILE=filMCMC, &
    !                 NAME=nammu, COMPRESS=compress, VERBOSE=debug, &
    !                 IND2=[counter,counter])
    ! CALL WRITE_HDF5(DBLARR2D=sigmcmc(:,icurr:icurr), FILE=filMCMC, &
    !                 NAME=namsig, COMPRESS=compress, VERBOSE=debug, &
    !                 IND2=[counter,counter])

    !! Check output for debugging
    DO i=1,MERGE(2,1,debug)
      WRITE(unitlog(i),*) 'Counter = ', counter
      ! PRINT*, "  mu = ", mumcmc(:,icurr)
      ! IF (Nx*Ny > 1) THEN  
        ! PRINT*, "  sigma = ", sigmcmc(:,icurr)
        ! PRINT*, "  corr = ", corrmcmc(:,icurr)
      ! END IF
      ! IF (calib) PRINT*, "  ln1pd = ",ln1pdmcmc(:,icurr)
      WRITE(unitlog(i),*)
    END DO
    
    !! Check for NaNs and eventually stop the code
    ! IF (ANY(isNaN(mumcmc(:,icurr)))) &
    !   CALL STRIKE("FITHB", &
    !               "a NaN appeared in the MCMC at iteration "//TRIMLR(PRING(counter)))

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
  ALLOCATE(ratpix(Nx,Ny,Nrat,t_end), meanrat(Nx,Ny,Nrat), stdevrat(Nx,Ny,Nrat))

  CALL READ_HDF5(DBLARR4D=par4D, FILE=filMCMC, &
                 NAME=nampar, IND4=[1,t_end])
  parpix(:,:,:,:) = par4D(:,:,:,:)
  meanpar(:,:,:) = MEAN(parpix(:,:,:,t_burnin:t_end), DIM=4)
  stdevpar(:,:,:) = SIGMA(parpix(:,:,:,t_burnin:t_end), DIM=4)

  CALL READ_HDF5(DBLARR4D=par4D, FILE=filMCMC, &
                 NAME=namrat, IND4=[1,t_end])
  ratpix(:,:,:,:) = par4D(:,:,:,:)
  meanrat(:,:,:) = MEAN(ratpix(:,:,:,t_burnin:t_end), DIM=4)
  stdevrat(:,:,:) = SIGMA(ratpix(:,:,:,t_burnin:t_end), DIM=4)
  
  CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME="Parameter label", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
  CALL WRITE_HDF5(DBLARR3D=meanpar(:,:,:), NAME="Mean of parameter value", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=stdevpar(:,:,:), NAME="Sigma of parameter value", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(STRARR1D=namBcorr(:), NAME="Correlated band label", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(STRARR1D=ratname(:), NAME="Band ratio label", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=meanrat(:,:,:), NAME="Mean of band ratio value", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=stdevrat(:,:,:), NAME="Sigma of band ratio value", &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  
  !! Calculate model
  !!-----------------
  ALLOCATE(FnuMOD(Nx,Ny,NwOBS))

  FnuMOD(:,:,:) = specModel( wOBS(:), INDPAR=ind, PARVAL=meanpar(:,:,:), &
                             QABS=Qabs(:), EXTINCT=extinct(:), &
                             FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
                             PABS=Pabs, FNULINE=FnuLINE, &
                             FNUCONT_TAB=FnuCONT_tab, FNUBAND_TAB=FnuBAND_tab, &
                             FNUSTAR_TAB=FnuSTAR_tab, PABS_TAB=Pabs_tab, &
                             FNULINE_TAB=FnuLINE_tab )

  DO i=1,Npabs
    FnuCONT_tab(:,:,:,i) = FnuCONT_tab(:,:,:,i) * Pabs(:,:,:)
  END DO 
  CALL WRITE_HDF5(DBLARR4D=FnuCONT_tab, NAME='FnuCONT ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  DO i=1,Nline
    FnuLINE_tab(:,:,:,i) = FnuLINE_tab(:,:,:,i) + (FnuCONT(:,:,:)+FnuSTAR(:,:,:))*Pabs(:,:,:)
  END DO
  CALL WRITE_HDF5(DBLARR4D=FnuLINE_tab, NAME='FnuLINE ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  DO i=1,Nband
    FnuBAND_tab(:,:,:,i) = FnuBAND_tab(:,:,:,i) + (FnuCONT(:,:,:)+FnuSTAR(:,:,:))*Pabs(:,:,:)
  END DO
  CALL WRITE_HDF5(DBLARR4D=FnuBAND_tab, NAME='FnuBAND ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  DO i=1,Nstar
    FnuSTAR_tab(:,:,:,i) = FnuSTAR_tab(:,:,:,i) * Pabs(:,:,:)
  END DO
  CALL WRITE_HDF5(DBLARR4D=FnuSTAR_tab, NAME='FnuSTAR ('//TRIMLR(spec_unit)//')', &
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
