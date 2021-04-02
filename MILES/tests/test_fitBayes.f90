MODULE fitBayes_external

  USE auxil, ONLY: parinfo_type, indpar_type, Qabs_type
  USE utilities, ONLY: DP, isNaN
  USE inout, ONLY: lenpar
  IMPLICIT NONE
  PRIVATE

  ! INTEGER, PARAMETER, PUBLIC :: Df_prior = 8

  INTEGER, SAVE, PUBLIC :: Nx, Ny, NwOBS
  INTEGER, SAVE, PUBLIC :: xOBS, yOBS, jw, ipar, ihpar
  INTEGER, SAVE, PUBLIC :: Nparhyp, Ncorrhyp
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wOBS, nuOBS
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC :: extinct
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parcurr, parhypcurr
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC :: FnuOBS, dFnuOBS
  ! LOGICAL, DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC :: maskS
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC :: mask, maskextra
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC :: maskhyp, maskpar

  !! Init param
  !!------------
  TYPE(indpar_type), SAVE, PUBLIC                             :: ind
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parinfo
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC    :: Qabs
  
  PUBLIC :: lnpost_par, lnlhobs_par!, lnpost_mu, lnpost_sig
  ! PUBLIC :: covariance

CONTAINS

  !! Build the covariance matrix
  !!-----------------------------
  ! PURE FUNCTION covariance(sig, corr)

  !   USE utilities, ONLY: DP
  !   USE statistics, ONLY: corr2Rmat
  !   IMPLICIT NONE

  !   REAL(DP), INTENT(IN), DIMENSION(Nparhyp)  :: sig
  !   REAL(DP), INTENT(IN), DIMENSION(Ncorrhyp) :: corr
  !   REAL(DP), DIMENSION(Nparhyp,Nparhyp)      :: covariance

  !   REAL(DP), DIMENSION(Nparhyp,Nparhyp)      :: Smat, Rmat

  !   Smat(:,:) = UNPACK(sig(:), maskS(:,:), FIELD=0._DP)
  !   Rmat(:,:) = CORR2RMAT(corr(:), Nparhyp) ! correlation matrix
  !   covariance = MATMUL( MATMUL(Smat(:,:), Rmat(:,:)), Smat(:,:) )

  ! END FUNCTION covariance

  !!-------------------------------------------------------
  !!
  !!            Likelihoods of the observations
  !!
  !!-------------------------------------------------------

  !! For sampling delta
  !!--------------------
  ! FUNCTION lnLHobs_del(ln1pgrid)

  !   USE utilities, ONLY: DP
  !   IMPLICIT NONE

  !   REAL(DP), INTENT(IN), DIMENSION(:)   :: ln1pdgrid
  !   REAL(DP), DIMENSION(SIZE(ln1pdgrid)) :: lnLHobs_del
    
  ! END FUNCTION lnLHobs_del
  
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
                             PARINFO=parinfo(:), INDPAR=ind, PARVAL=parcurr(:), &
                             QABS=Qabs(:), EXTINCT=extinct(:,:))!, debug=.TRUE.)

    !! Likelihoods
    varred(:,:) = 0._DP
    FORALL (iw=1:NwOBS) &
      varred(:,iw) = ( FnuOBS(xOBS,yOBS,iw) - Fnu_mod(:,iw) ) / dFnuOBS(xOBS,yOBS,iw)

    lnpind(:,:) = - 0.5_DP * varred(:,:)**2

    !! Global likelihood, accounting for the excesses
    lnLHobs_par(:) = SUM(lnpind(:,:), DIM=2)!, MASK=mask2D(:,:))
    ! lnLHobs_par(:) = SUM(Fnu_mod(:,:), DIM=2)!, MASK=mask2D(:,:))

    !! Extra parameters
    !!-----------------
    !! Likelihood

  END FUNCTION lnLHobs_par

  !!-------------------------------------------------------
  !!
  !!                  Prior distributions
  !!
  !!-------------------------------------------------------

  !! For calibration errors
  !!------------------------
  ! FUNCTION lnprior_del(ln1pdgrid)
    
  !   USE utilities, ONLY: DP
  !   IMPLICIT NONE

  !   REAL(DP), DIMENSION(:), INTENT(IN)   :: ln1pdgrid
  !   REAL(DP), DIMENSION(SIZE(ln1pdgrid)) :: lnprior_del
    
  ! END FUNCTION lnprior_del
  
  !! For sampling physical parameters
  !!----------------------------------
  ! FUNCTION lnprior_par(pargrid)

  !   USE utilities, ONLY: DP
  !   USE arrays, ONLY: iwhere
  !   IMPLICIT NONE

  !   REAL(DP), DIMENSION(:), INTENT(IN) :: pargrid
  !   REAL(DP), DIMENSION(SIZE(pargrid)) :: lnprior_par

  !   INTEGER :: ip, i, ih, Ngrid

  !   !! Cube of parameters
  !   Ngrid = SIZE(pargrid(:))
  !   CALL IWHERE(i2ih(:) == ipar,ip)
  !   allpargrid(:,:) = 0._DP
  !   FORALL (ih=1:Nparhyp,maskhypcurr(ih)) &
  !     allpargrid(ih,:) = MERGE( SPREAD(parhypcurr(ih),DIM=1,NCOPIES=Ngrid), &
  !                               pargrid(:), (ip /= ih) ) - mucurr(ih)
    
  !   !! Hyperdistribution
  !   FORALL (i=1:Ngird) &
  !     lnprior_par(i) = LOG( DOT_PRODUCT(allpargrid(:,i), &
  !                             MATMUL(invcov_prev(:,:), allpargrid(:,i))) &
  !                           / Df_prior + 1._DP ) * (-Df_prior+Nparhyp)/2._DP)

  ! END FUNCTION lnprior_par

  !! For the variance
  !!------------------
  ! FUNCTION lnprior_sig(lnSgrid)

  !   USE utilities, ONLY: DP
  !   IMPLICIT NONE

  !   REAL(DP), DIMENSION(:), INTENT(IN) :: lnSgrid
  !   REAL(DP), DIMENSION(SIZE(lnSgrid)) :: lnprior_sig
    
  !   lnprior_sig(:) = - 0.5_DP * ( (lnSgrid(:)-cenlnS0(ihpar)) / siglnS0 )**2

  ! END FUNCTION lnprior_sig

  !! For the correlation coefficient
  !!---------------------------------
  ! FUNCTION lnprior_corr(corrgrid)

  !   USE utilities, ONLY: DP
  !   IMPLICIT NONE

  !   REAL(DP), DIMENSION(:), INTENT(IN) :: corrgrid
  !   REAL(DP), DIMENSION(SIZE(corrgrid)) :: lnprior_corr

  ! END FUNCTION lnprior_corr

  !!-------------------------------------------------------
  !!
  !!             Hyperparameter distributions
  !!
  !!-------------------------------------------------------

  !! For sampling the variance
  !!---------------------------
  ! FUNCTION lnhyper_sig(lnSgrid)

  !   USE utilities, ONLY: DP
  !   USE matrices, ONLY: invert_cholesky
  !   IMPLICIT NONE

  !   REAL(DP), DIMENSION(:), INTENT(IN) :: lnSgrid
  !   REAL(DP), DIMENSION(SIZE(lnSgrid)) :: lnhyper_sig

  ! END FUNCTION lnhyper_sig

  !! For sampling the correlation
  !!------------------------------
  ! FUNCTION lnhyper_corr(corrgrid)

  !   USE utilities, ONLY: DP
  !   USE matrices, ONLY: invert_cholesky
  !   IMPLICIT NONE

  !   REAL(DP), DIMENSION(:), INTENT(IN)   :: corrgrid
  !   REAL(DP), DIMENSION(SIZE(corrgrid))  :: lnhyper_corr

  ! END FUNCTION lnhyper_corr

  !!-------------------------------------------------------
  !!
  !!                Posterior distributions
  !!
  !!-------------------------------------------------------

  !! For sampling delta
  !!--------------------
  ! FUNCTION lnpost_del

  !   USE utilities, ONLY: DP
  !   IMPLICIT NONE

  !   REAL(DP), INTENT(IN), DIMENSION(:)   :: ln1pdgrid
  !   REAL(DP), DIMENSION(SIZE(ln1pdgrid)) :: lnpost_del
    
  !   IF (usefilt(jw)) THEN
  !     lnpost_del(:) = LNLHOBS_DEL(ln1pdgrid(:)) + LNPRIOR_DEL(ln1pdgrid(:))
  !   ELSE 
  !     lnpost_del(:) = LNPRIOR_DEL(ln1pdgrid(:))
  !   END IF

  ! END FUNCTION lnpost_del

  !! For sampling physical parameters
  !!----------------------------------
  FUNCTION lnpost_par(pargrid)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: pargrid
    REAL(DP), DIMENSION(SIZE(pargrid)) :: lnpost_par

    ! IF (parinfo(ipar)%hyper) THEN
    !   lnpost_par(:) = LNLHOBS_PAR(pargrid(:)) + LNPRIOR_PAR(pargrid(:))
    ! ELSE
      lnpost_par(:) = LNLHOBS_PAR(pargrid(:))
    ! END IF

  END FUNCTION lnpost_par

  !! For sampling the mean
  !!-----------------------
  ! FUNCTION lnpost_mu(mugrid)

  !   USE utilities, ONLY: DP
  !   IMPLICIT NONE

  !   REAL(DP), DIMENSION(:), INTENT(IN) :: mugrid
  !   REAL(DP), DIMENSION(SIZE(mugrid))  :: lnpost_mu

  ! END FUNCTION lnpost_mu

  !! For sampling the variance
  !!---------------------------
  ! FUNCTION lnpost_sig(lnSgrid)

  !   USE utilities, ONLY: DP
  !   IMPLICIT NONE

  !   REAL(DP), DIMENSION(:), INTENT(IN) :: lnSgrid
  !   REAL(DP), DIMENSION(SIZE(lnSgrid)) :: lnpost_sig

  !   lnpost_sig(:) = LNHYPER_SIG(lnSgrid(:)) + LNPRIOR_SIG(lnSgrid(:))

  ! END FUNCTION lnpost_sig

  !! For sampling the correlation
  !!------------------------------
!   FUNCTION lnpost_corr(corrgrid)

!     USE utilities, ONLY: DP
!     IMPLICIT NONE

!     REAL(DP), DIMENSION(:), INTENT(IN)      :: corrgrid
!     REAL(DP), DIMENSION(SIZE(corrgrid(:)))  :: lnpost_corr

!     lnpost_corr(:) = LNHYPER_CORR(corrgrid(:)) + LNPRIOR_CORR(corrgrid(:))

!   END FUNCTION lnpost_corr

END MODULE fitBayes_external


!!==========================================================================
!!                           Main Program
!!==========================================================================


PROGRAM test_fitBayes

  USE datable, ONLY: cband_sig
  USE auxil, ONLY: read_master, set_indpar, specModel, initparam, degradeRes
  USE utilities, ONLY: DP, pring, trimLR, trimeq, swap, timinfo, &
                       banner_program, ustd, isNaN, initiate_clock, time_type
  USE arrays, ONLY: iwhere, closest
  USE constants, ONLY: MKS
  USE statistics, ONLY: mean, sigma
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, lenpar, lenpath
  USE grain_optics, ONLY: lendustQ
  USE random, ONLY: generate_newseed, rand_general
  USE fitBayes_external, ONLY: Nx, Ny, NwOBS, wOBS, nuOBS, xOBS, yOBS, &
                            FnuOBS, dFnuOBS, &
                            ipar, &!ihpar, &
                            parcurr, &!parhypcurr, maskS, &
                            lnpost_par, lnlhobs_par, &!lnpost_mu, lnpost_sig, lnpost_corr, &
                            ! covariance,
                            mask, maskpar, &!maskextra, maskhyp, &
                            ind, parinfo, Qabs, extinct
  IMPLICIT NONE

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: Ngibbsmax = 3000
  REAL(DP), PARAMETER :: accrand = 1.E-3_DP
  CHARACTER(*), PARAMETER :: nampar = "Parameter values"
  CHARACTER(*), PARAMETER :: namln1pd = "ln(1+delta)"
  CHARACTER(*), PARAMETER :: nammu = "Mean of hyperdistribution"
  CHARACTER(*), PARAMETER :: namsig = "Sigma of hyperdistribution"
  CHARACTER(*), PARAMETER :: namcorr = "Correlation of hyperdistribution"
  CHARACTER(*), PARAMETER :: filOBS = './dat/observations_fitMIR'//h5ext
  CHARACTER(*), PARAMETER :: dirOUT = './out/'
  CHARACTER(*), PARAMETER :: filOUT = dirOUT//'test_fitBayes'//h5ext
  CHARACTER(*), PARAMETER :: filMCMC = dirOUT//'parlog_fitBayes'//h5ext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Input variables
  INTEGER :: i, iw, x, y!, j, ih
  INTEGER :: Npar, Nmcmc, NiniMC, Nsou!, Nparfree, Nwfree
  INTEGER :: Ncont, Nband, Nline, Nextc, Nstar, Nextra
  INTEGER, DIMENSION(:), ALLOCATABLE :: itied
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: maskint ! convert mask=0 to mask=T
  REAL(DP) :: val, sig
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parini, parmcmc
  LOGICAL :: verbose, calib, newseed, newinit, dostop
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: maskxy3D
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB, labL

  !! MCMC parameters
  INTEGER :: icurr, iprev, counter
  REAL(DP), DIMENSION(2) :: lim
  
  !! Output variables
  INTEGER, DIMENSION(2) :: unitlog
  TYPE(time_type) :: timestr

  !! Analysis variables
  INTEGER :: t_burnin, t_end, Nmcmc_eff
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: meanpar, stdevpar
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT, &
    FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parpix, par4D
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FnuCONT_tab, &
    FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab

  ! Output settings
  !----------------
  CALL INITIATE_CLOCK(timestr)
  unitlog(:) = [ ulog, ustd ]
  
  !!------------------------------------------------------------------------
  !!                            I. Read the inputs
  !!------------------------------------------------------------------------

  !! a. Main settings
  !!------------------
  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='Wavelength (microns)', N1=NwOBS)
  CALL READ_MASTER(WAVALL=wOBS(:), &
                   VERBOSE=verbose, NMCMC=Nmcmc, NINIMC=NiniMC, &
                   CALIB=calib, NEWSEED=newseed, NEWINIT=newinit, &
                   LABQ=labQ, LABL=labL, LABB=labB, QABS=Qabs, EXTINCT=extinct, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, &
                   NEXTC=Nextc, NSTAR=Nstar, NEXTRA=Nextra, DOSTOP=dostop, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit)

  IF (newseed) CALL GENERATE_NEWSEED()
  IF (verbose) PRINT*
  DO i=1,MERGE(2,1,verbose)
    CALL BANNER_PROGRAM("HIBARI: HIerarchical BAyesian fitting Routine of mid-IR emission", &
                        UNIT=unitlog(i), SWING=.True.)
  END DO
  
  !! b. Observation sample
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
  !! - MASKXY3D(Nx,Ny,NwOBS) is a 2D mask (Nx,Ny) replicated NwOBS times
  !!   identifying pixels which have at least one defined wavelength. Thus, it
  !!   is False if all wavelengths of the pixel are missing. [ANY(wOBS)]
  !! - MASKHYP(Nx,Ny,Nparhyp) selects the parameters of pixels which having to
  !!   be treated hierarchically. For dust model parameters, this mask is a
  !!   replicate of MASKXY3D(:,:,1). For extra parameters, it is MASKEXTRA.
  !! - MASKPAR(Nx,Ny,Npar) is the same as MASKHYP but for the parameter list.
  ALLOCATE(mask(Nx,Ny,NwOBS), maskxy3D(Nx,Ny,NwOBS))
  mask(:,:,:) = ( maskint(:,:,:) == 0 )
  WHERE (isNaN(FnuOBS(:,:,:)) .OR. .NOT. mask(:,:,:)) FnuOBS(:,:,:) = 0._DP
  DO iw=1,NwOBS
    maskxy3D(:,:,iw) = ANY(mask(:,:,:),DIM=3)
  END DO

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
  ! CALL IWHERE(( .NOT. parinfo(:)%fixed ) .AND. ( itied(:) == 0 ),ifree)
  ! Ncons = NwOBS - Nfreefilt
  ! Nfreepar = COUNT( (.NOT. parinfo(:)%fixed ) .AND. ( itied(:) == 0 ) )
  ! Ndof = Ncons - Nfreepar

  !! Read the extra parameters (TBD)

  !! Units
  ! unitmks = TRIMEQ(SED_unit,"W.m-2.Hz-1")

  Nsou = COUNT(ANY(mask(:,:,:),DIM=3))

  !! 3D hyper mask
  ! ALLOCATE (maskhyp(Nx,Ny,Nparhyp))
  ! DO i=1,Nparhyp
  !   IF (parinfo(i2ih(i))%model) THEN
  !     maskhyp(:,:,i) = ANY(mask(:,:,:),DIM=3)
  !   ELSE
  !     maskhyp(:,:,i) = maskextra(:,:,i2ih(i)-iext0+1)
  !   END IF
  ! END DO

  !! 3D parameter mask
  ALLOCATE (maskpar(Nx,Ny,Npar))
  DO i=1,Npar
    ! IF (parinfo(i)%model) THEN
      maskpar(:,:,i) = ALL(mask(:,:,:),DIM=3)
    ! ELSE
    !   maskpar(:,:,i) = maskextra(:,:,i-iext0+1)
    ! END IF
  END DO  
  
  print*, 'Read the inputs [done]'//NEW_LINE('')

  print*, '=================================================='
  
  !!------------------------------------------------------------------------
  !!                          II. Prepare the MCMC
  !!------------------------------------------------------------------------
  
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) "ESTIMATING THE INITIAL VALUES OF THE PARAMETERS (" &
                        //TRIMLR(TIMINFO(timestr))//")"
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
  
  ! Initialize the parameters of the chain
  ALLOCATE (parmcmc(Nx,Ny,Npar,2))!,ln1pdmcmc(Ncalib,2))
  DO i=1,2 
    FORALL (x=1:Nx,y=1:Ny,ANY(maskpar(x,y,:))) parmcmc(x,y,:,i) = parini(x,y,:,1)
    ! ln1pdmcmc(:,i) = ln1pd0(:)
  END DO
  
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
  ! IF (ANY(parinfo(:)%hyper)) THEN
  !   ALLOCATE (allparhypcurr(Nx,Ny,Nparhyp)) ! parameter grid
  !   ALLOCATE (cov_prev(Nparhyp,Nparhyp))    ! correlation matrix
  !   ALLOCATE (invcov_prev(Nparhyp,Nparhyp)) ! inverse correlation matrix
  !   ALLOCATE (mucurr(Nparhyp))              ! current average
  !   ALLOCATE (sigcurr(Nparhyp))             ! current average
  !   ALLOCATE (corrcurr(Ncorrhyp))           ! current correlation coefficient
  !   ALLOCATE (parhypcurr(Nparhyp))          ! current parameters
  !   ALLOCATE (maskhypcurr(Nparhyp))         ! current parameters' mask
  ! END IF

  !! Initialize the files
  CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME="Parameter label", &
                  FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
  CALL WRITE_HDF5(INTARR1D=[Nmcmc], NAME="Length of MCMC", &
                  FILE=filMCMC, COMPRESS=compress, verbose=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(INITDBLARR=[Nx,Ny,Npar,Nmcmc], NAME=nampar, &
                  FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  ! CALL WRITE_HDF5(INITDBLARR=[Nparhyp,Nmcmc], NAME=nammu, &
  !                 FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  ! CALL WRITE_HDF5(INITDBLARR=[Nparhyp,Nmcmc], NAME=namsig, &
  !                 FILE=filMCMC, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

  print*, 'Initial parameters [done]'//NEW_LINE('')

  print*, '=================================================='

  !!------------------------------------------------------------------------
  !!                          III. Run the MCMC
  !!------------------------------------------------------------------------

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

    !! Draw the model parameters
    parmcmc(:,:,:,icurr) = parmcmc(:,:,:,iprev)
    param: DO ipar=1,Npar
      IF (.NOT. parinfo(ipar)%fixed) THEN
        IF (debug) PRINT*, " - "//TRIMLR(parinfo(ipar)%name) 
        lim(:) = parinfo(ipar)%limits(:)
        xsource: DO xOBS=1,1
          ysource: DO yOBS=1,1
        ! xsource: DO xOBS=1,Nx
        !   ysource: DO yOBS=1,Ny
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

              IF (debug) PRINT*, "par(",xOBS,yOBS,ipar,") = ", &
                                 parmcmc(xOBS,yOBS,ipar,icurr)
              
            END IF
          END DO ysource
        END DO xsource
      END IF
    END DO param

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
    !                    parmcmc(:,:,i2ih(icorr2ij(icorr,2)),icurr), &
    !                    MASK=(maskhyp(:,:,icorr2ij(icorr,1)) &
    !                          .AND. maskhyp(:,:,icorr2ij(icorr,2))))

    ! END IF hierarchy
    
    !! 4) Outputs
    !!------------
    CALL WRITE_HDF5(DBLARR4D=parmcmc(:,:,:,icurr:icurr), FILE=filMCMC, &
                    NAME=nampar, COMPRESS=compress, VERBOSE=debug, &
                    IND4=[counter,counter])
    ! CALL WRITE_HDF5(DBLARR2D=mumcmc(:,icurr:icurr), FILE=filMCMC, &
    !                 NAME=nammu, COMPRESS=compress, VERBOSE=debug, &
    !                 IND2=[counter,counter])
    ! CALL WRITE_HDF5(DBLARR2D=sigmcmc(:,icurr:icurr), FILE=filMCMC, &
    !                 NAME=namsig, COMPRESS=compress, VERBOSE=debug, &
    !                 IND2=[counter,counter])

    !! Check output for debugging
    IF (debug) THEN
      PRINT*, "Counter = ", counter
      ! PRINT*, "  mu = ", mumcmc(:,icurr)
      ! IF (Nx*Ny > 1) THEN  
        ! PRINT*, "  sigma = ", sigmcmc(:,icurr)
        ! PRINT*, "  corr = ", corrmcmc(:,icurr)
      ! END IF
      ! IF (calib) PRINT*, "  ln1pd = ",ln1pdmcmc(:,icurr)
      PRINT*
    END IF
    
    !! Check for NaNs and eventually stop the code
    ! IF (ANY(isNaN(mumcmc(:,icurr)))) &
    !   CALL STRIKE("FITBayes", &
    !               "a NaN appeared in the MCMC at iteration "//TRIMLR(PRING(counter)))

    !! Normal termination condition
    IF (counter == Nmcmc) EXIT

    !! Update the indices
    CALL SWAP(icurr,iprev)
    counter = counter + 1
    
  END DO big_loop

  
  print*, 'Bayes fit M83 IRS obs [done]'//NEW_LINE('')

  print*, '=================================================='
  
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

  DO i=1,Nextc
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
  
  print*, 'fitBayes analysis [done]'//NEW_LINE('')

  print*, '=================================================='

  !! Final
  !!-------
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) "PROGRAM EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
    WRITE(unitlog(i),*)
  END DO
  
  !! Free memory space
  DEALLOCATE(wOBS, FnuOBS, dFnuOBS, FnuMOD, maskint, mask)

END PROGRAM test_fitBayes
