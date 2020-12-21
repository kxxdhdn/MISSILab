MODULE fitHB_external

  USE auxil, ONLY: parinfo_type, indpar_type, Qabs_type
  USE utilities, ONLY: DP
  USE inout, ONLY: lenpar
  USE grain_optics, ONLY: lendustQ
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: Df_prior = 8

  INTEGER, SAVE, PUBLIC :: Nx, Ny, NwOBS
  INTEGER, SAVE, PUBLIC :: xOBS, yOBS, jw, ipar, ihpar
  INTEGER, SAVE, PUBLIC :: Nparhyp, Ncorrhyp
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wOBS, nuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parcurr, parhypcurr
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: FnuOBS, dFnuOBS
  ! LOGICAL, DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC    :: maskS
  ! LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC  :: maskhyp, maskpar

  !! Init param
  !!------------
  TYPE(indpar_type), SAVE, PUBLIC                             :: ind
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parinfo
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC    :: Qabs
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC            :: indBIN, indLIN
  
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
     ! print*, SIZE(specModel(wOBS(:), PARVEC=pargrid(:), PARNAME=parinfo(ipar)%name, &
     !                         PARINFO=parinfo(:), INDPAR=ind, PARVAL=parcurr(:), &
     !                         QABS=Qabs, INDBIN=indBIN, INDLIN=indLIN),1)
     ! print*, SIZE(specModel(wOBS(:), PARVEC=pargrid(:), PARNAME=parinfo(ipar)%name, &
     !                         PARINFO=parinfo(:), INDPAR=ind, PARVAL=parcurr(:), &
     !                         QABS=Qabs, INDBIN=indBIN, INDLIN=indLIN),2)
     ! print*, shape(fnu_mod)
    Fnu_mod(:,:) = specModel(wOBS(:), PARVEC=pargrid(:), PARNAME=parinfo(ipar)%name, &
                             PARINFO=parinfo(:), INDPAR=ind, PARVAL=parcurr(:), &
                             QABS=Qabs(:), INDBIN=indBIN, INDLIN=indLIN)

    !! Likelihoods
    varred(:,:) = 0._DP
    FORALL (iw=1:NwOBS) &
      varred(:,iw) = ( FnuOBS(iw) - Fnu_mod(:,iw) ) / dFnuOBS(iw)

    lnpind(:,:) = - 0.5_DP * varred(:,:)**2

    !! Global likelihood, accounting for the excesses
    lnLHobs_par(:) = SUM(lnpind(:,:), DIM=2)!, MASK=mask2D(:,:))

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

END MODULE fitHB_external


!!==========================================================================
!!                           Main Program
!!==========================================================================


PROGRAM test_fitHB

  USE auxil, ONLY: read_master, set_indpar, specModel, initparam
  USE utilities, ONLY: DP, pring, trimLR, trimeq, timinfo, &
                       banner_program, ustd, isNaN, initiate_clock, time_type
  USE arrays, ONLY: iwhere, closest
  USE constants, ONLY: MKS
  USE statistics, ONLY: MEDIAN
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, write_ascii, ascext, lenpar, lenpath
  USE grain_optics, ONLY: lendustQ
  USE random, ONLY: generate_newseed
  USE fitHB_external, ONLY: Nx, Ny, NwOBS, wOBS, nuOBS, FnuOBS, dFnuOBS, &
                            residuals, Fnu_mod, resid, &
                            ipar, &!ihpar, &
                            parcurr, &!parhypcurr, maskS, &
                            lnpost_par, lnlhobs_par, &!lnpost_mu, lnpost_sig, lnpost_corr, &
                            ! covariance, maskhyp, maskpar, &
                            ind, parinfo, Qabs, xOBS, yOBS
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
  CHARACTER(*), PARAMETER :: filOUT = dirOUT//'test_fitHB'//h5ext
  CHARACTER(*), PARAMETER :: filMCMC = dirOUT//'parlog_fitHB'//h5ext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Input variables
  INTEGER :: i, j, x, y, Npar, Nparfree, Nwfree, Nmcmc
  INTEGER :: Ncont, Nband, Nline, Nextra
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: maskint ! convert mask=0 to mask=T
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parini
  LOGICAL :: verbose, calib, newseed, newinit, dostop
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
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD

  ! Output settings
  !----------------
  CALL INITIATE_CLOCK(timestr)
  unitlog(:) = [ ulog, ustd ]
  
  !! Wavelengths -> indpdt identically dist. (diagonal covar matrix)
  iid = .TRUE.
  
  !!------------------------------------------------------------------------
  !!                            I. Read the inputs
  !!------------------------------------------------------------------------

  !! a. Main settings
  !!------------------
  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='Wavelength (microns)', N1=NwOBS)
  CALL READ_MASTER(WAVALL=wOBS(:), &
                   VERBOSE=verbose, NMCMC=Nmcmc, &
                   CALIB=calib, NEWSEED=newseed, NEWINIT=newinit, &
                   LABQ=labQ, LABL=labL, LABB=labB, QABS=Qabs, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, NEXTRA=Nextra, DOSTOP=dostop, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit)

  IF (newseed) CALL GENERATE_NEWSEED()
  IF (verbose) PRINT*
  DO i=1,MERGE(2,1,verbose)
    CALL BANNER_PROGRAM("LE MIROIR (LEast-squares fitting of Mid-IR emission) " &
                        //"OptImized Routine", UNIT=unitlog(i), SWING=.True.)
  END DO
  
  !! b. Observation sample
  !!-----------------------
  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, NAME='FnuOBS (MKS)', N1=Nx, N2=Ny)
  CALL READ_HDF5(DBLARR3D=dFnuOBS, FILE=filOBS, NAME='dFnuOBS (MKS)')
  !! On input, the convention is mask=1 => blocks; mask=0 => passes.
  !! However, within Fortran, the convention is the Fortran mask=T => passes.
  CALL READ_HDF5(INTARR3D=maskint, FILE=filOBS, NAME='NaN mask')
  print*, size(maskint, dim=3)
  ALLOCATE(mask(Nx,Ny,NwOBS))
  mask(:,:,:) = ( maskint(:,:,:) == 0 )
  
  ALLOCATE(nuOBS(NwOBS))
  nuOBS(:) = MKS%clight/MKS%micron / wOBS(:)

  !! dfnu positive -> 0 ???
  ! WHERE (isNaN(FnuOBS(:,:,:))) FnuOBS(:,:,:) = 0._DP
  ! WHERE (isNaN(dFnuOBS(:,:,:))) dFnuOBS(:,:,:) = 0._DP
  
  
  Nwfree = COUNT( ALL(ALL(.NOT. mask(:,:,:) &
                            .AND. dFnuOBS(:,:,:) > 0._DP,DIM=1),DIM=1) )
  IF (Nwfree > 0) &
    CALL IWHERE(ALL(ALL(.NOT. mask(:,:,:) &
                        .AND. dFnuOBS(:,:,:) > 0._DP,DIM=1),DIM=1), iwfree)

  !! Median signal-to-noise ratio
  ! ALLOCATE(medSovN(Nx,Ny))
  ! medSovN(:,:) = 0._DP
  ! DO x=1,Nx
  !   DO y=1,Ny
  !     IF (ANY(mask(x,y,:))) THEN
  !       medSovN(x,y) = MEDIAN(FnuOBS(x,y,:)/dFnuOBS(x,y,:),MASK=mask(x,y,:))
  !       IF (medSovN(x,y) < 0._DP) medSovN(x,y) = 0._DP

  !     END IF
  !   END DO
  ! END DO

  !! Read the instrumental covariance matrix (TBD)


  !! Constraints
  ALLOCATE(itied(Npar))
  itied(:) = 0
  DO i=1,Npar
    IF (.NOT. TRIMEQ(parinfo(i)%tied,"")) &
      CALL IWHERE(TRIMEQ(parinfo(:)%name,parinfo(i)%tied),itied(i))
  END DO
  ! CALL IWHERE(( .NOT. parinfo(:)%fixed ) .AND. ( itied(:) == 0 ),ifree)
  ! Ncons = NwOBS - Nfreefilt
  ! Nfreepar = COUNT( (.NOT. parinfo(:)%fixed ) .AND. ( itied(:) == 0 ) )
  ! Ndof = Ncons - Nfreepar

  !! Units
  ! unitmks = TRIMEQ(SED_unit,"W.m-2.Hz-1")
  ! power_unit = MERGE("W.m-2","Lsun ",unitmks)
  
  !! c. Compute the covariance matrix of the uncertainties
  !!-------------------------------------------------------
  !! Mask for the RMS covariance (diagonal)
  ALLOCATE(maskS(NwOBS,NwOBS))
  maskS(:,:) = .FALSE.
  FORALL (i=1:NwOBS) maskS(i,i) = .TRUE.

  !! 4D mask for the covariance matrix
  ! ALLOCATE(mask4D(Nx,Ny,NwOBS,NwOBS))
  ! mask4D(:,:,:,:) = .TRUE.
  ! FORALL (x=1:Nx,y=1:Ny,i=1:NwOBS, .NOT. mask(x,y,i))
  !   mask4D(x,y,i,:) = .FALSE.
  !   mask4D(x,y,:,i) = .FALSE.
  ! END FORALL

  !! Covariance matrix for each pixel and its inverse
  ALLOCATE(covarOBS(Nx,Ny,NwOBS,NwOBS), invLcovarOBS(Nx,Ny,NwOBS,NwOBS))

  covarOBS(:,:,:,:) = 0._DP
  !! invLcovarOBS is a lower triangle matrix L (Likelihood) results from 
  !! the Cholesky decomposition of the covariance matrix (Appendix C2, Galliano18a)
  invLcovarOBS(:,:,:,:) = 0._DP
  FORALL (x=1:Nx,y=1:Ny,ALL(mask(x,y,:)))
    covarOBS(x,y,:,:) = UNPACK(dFnuOBS(x,y,:)**2, maskS, FIELD=0._DP)
    invLcovarOBS(x,y,:,:) = UNPACK(1._DP/dFnuOBS(x,y,:), maskS, FIELD=0._DP)

  END FORALL

  !! No wave calib (between diff instr/module)
  ! DO x=1,Nx
  !   DO y=1,Ny
  !     DO i=1,NwOBS
  !       IF (.NOT. mask(x,y,i)) THEN
  !         invLcovarOBS(x,y,i,:) = 0._DP
  !         invLcovarOBS(x,y,:,i) = 0._DP

  !       END IF
  !     END DO
  !   END DO
  ! END DO

  ! DEALLOCATE(covarOBS)

  print*, 'Read the inputs [done]'//NEW_LINE('')

  print*, '=================================================='
  
  !!------------------------------------------------------------------------
  !!                          II. Initial parameters
  !!------------------------------------------------------------------------
  
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) "ESTIMATING THE INITIAL VALUES OF THE PARAMETERS (" &
                        //TRIMLR(TIMINFO(timestr))//")"
  END DO

  !! a. Rough spectrum estimators
  !!------------------------------

  !! b. Automatic init param estimates
  !!-----------------------------------
  ALLOCATE(parini(Nx, Ny, Npar, MAX(NiniMC,1)))
  parini(:,:,:,:) = 0._DP
  CALL INITPARAM(NiniMC, IND=ind, PAR=parini(:,:,:,:), PARINFO=parinfo(:), &
                 ITIED=itied(:), MASK=mask(:,:,:), &
                 NEWINIT=newinit, FILOBS=filOBS, &
                 LABB=labB(:), LABL=labL(:), QABS=Qabs(:))
  
  Nparfree = COUNT((.NOT. parinfo(:)%fixed) .AND. (itied(:) <= 0))
  PRINT*, 'Number of free param: '//TRIMLR(PRING(Nparfree))

  print*, 'Initial parameters [done]'//NEW_LINE('')

  print*, '=================================================='

  !!------------------------------------------------------------------------
  !!                            Run the fitter
  !!------------------------------------------------------------------------
  
  ALLOCATE(par(Nx,Ny,Npar), status(Nx,Ny), resid(NwOBS), chi2red(Nx,Ny), &
           parerr(Nx,Ny,Npar), covpar(Nx,Ny,Npar,Npar), Niter(Nx,Ny), &
           Fnu_mod(NwOBS), limits(Npar,2), limited(Npar,2))

  FORALL (i=1:Npar)
    limits(i,:) = parinfo(i)%limits(:)
    limited(i,:) = parinfo(i)%limited(:)
    
  END FORALL
  
  ALLOCATE (statusMC(MAX(NiniMC,1)),NiterMC(MAX(NiniMC,1)), &
            parerrMC(Npar,MAX(NiniMC,1)),chi2redMC(MAX(NiniMC,1)), &
            covarMC(Npar,Npar,MAX(NiniMC,1)))

  par(:,:,:) = parini(:,:,:,1)
  CALL WRITE_HDF5(DBLARR3D=par(:,:,:), FILE=filOUT, &
                  NAME='Best fitted parameter value', COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.FALSE.)

  chi2fitx: DO xOBS=6,6
    chi2fity: DO yOBS=7,7
  ! chi2fitx: DO xOBS=1,Nx
  !   chi2fity: DO yOBS=1,Ny
      notmasked: IF (ANY( mask(xOBS,yOBS,:) )) THEN
        PRINT*
        PRINT*, 'pos: ('//TRIMLR(PRING(xOBS))//', '//TRIMLR(PRING(yOBS))//'): '

        !! a. Levenberg-Marquardt method
        !!-------------------------------
        pariniMC: IF (NiniMC == 0) THEN

          !! With single init param
          ! par(xOBS,yOBS,:) = parini(xOBS,yOBS,:,1)
          CALL CHI2MIN_LM (residuals, NwOBS, PAR=par(xOBS,yOBS,:), VERBOSE=.FALSE., &
                           STATUS=status(xOBS,yOBS), PARNAME=parinfo(:)%name, &
                           LIMITED=limited(:,:), LIMITS=limits(:,:), &
                           FIXED=parinfo(:)%fixed, ITIED=itied(:), &
                           CHI2RED=chi2red(xOBS,yOBS), NITER=Niter(xOBS,yOBS), &
                           PARERR=parerr(xOBS,yOBS,:), COVAR=covpar(xOBS,yOBS,:,:))
        ELSE

          !! Vary the init param and keep the best fit
          DO j=1,NiniMC
            CALL CHI2MIN_LM (residuals, NwOBS, PAR=parini(xOBS,yOBS,:,j), VERBOSE=.FALSE., &
                             STATUS=statusMC(j), PARNAME=parinfo(:)%name, &
                             LIMITED=limited(:,:), LIMITS=limits(:,:), &
                             FIXED=parinfo(:)%fixed, ITIED=itied(:), &
                             CHI2RED=chi2redMC(j), NITER=NiterMC(j), &
                             PARERR=parerrMC(:,j), COVAR=covarMC(:,:,j))
            
          END DO
          ibest = MAXVAL(MINLOC(chi2redMC(:)))
          par(xOBS,yOBS,:) = parini(xOBS,yOBS,:,ibest)
          parerr(xOBS,yOBS,:) = parerrMC(:,ibest)
          covpar(xOBS,yOBS,:,:) = covarMC(:,:,ibest)
          status(xOBS,yOBS) = statusMC(ibest)
          Niter(xOBS,yOBS) = NiterMC(ibest)
          chi2red(xOBS,yOBS) = chi2redMC(ibest)
          DO i=1,MERGE(2,1,verbose)
            WRITE(unitlog(i),*) "  MC chi2red is between " &
              //TRIMLR(PRING(MINVAL(chi2redMC(:)),NDEC=6))//" and " &
              //TRIMLR(PRING(MAXVAL(chi2redMC(:)),NDEC=6))
            
          END DO
        
        END IF pariniMC

        !! b. Derived quantities
        !!-----------------------

        !! c. Correlation between parameters
        !!-----------------------------------

        !! d. Summary
        !!------------
        DO i=1,MERGE(2,1,verbose)
          PRINT*, '>>>'
          PRINT*, "Checking output ("//TRIMLR(PRING(xOBS))//', '//TRIMLR(PRING(yOBS))//'): '
          PRINT*, "  status = "//TRIMLR(PRING(status(xOBS,yOBS)))
          PRINT*, "  chi2 = "//TRIMLR(PRING(chi2red(xOBS,yOBS),NDEC=10))
          PRINT*, "  Niter = "//TRIMLR(PRING(Niter(xOBS,yOBS)))
          ! CALL WRITE_ASCII(VEC1=par(xOBS,yOBS,:), VEC2=parerr(xOBS,yOBS,:), &
          !                  FILE="out/chi2min_fitChi2.txt")
          CALL WRITE_HDF5(DBLARR3D=par(:,:,:), FILE=filOUT, &
                          NAME='Best fitted parameter value', COMPRESS=compress, VERBOSE=debug, &
                          APPEND=.TRUE.)
          !! Renew init param for HB inputs
          ! CALL WRITE_HDF5(DBLARR3D=par(:,:,:), FILE=filOBS, &
          !                 NAME='Initial parameter value', COMPRESS=compress, VERBOSE=debug, &
          !                 APPEND=.TRUE.)
          PRINT*
          PRINT*, "Covariance matrix:"
          ! DO i=1,Npar
          !   PRINT*, REAL(covpar(xOBS,yOBS,:,i), KIND(0.))
          ! END DO
          PRINT*, '<<<'
          PRINT*

        END DO
       
      END IF notmasked
    END DO chi2fity
  END DO chi2fitx

  !! Renew init param for HB inputs
  CALL WRITE_HDF5(DBLARR3D=par(:,:,:), FILE=filOBS, &
                  NAME='Initial parameter value', COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  
  print*, 'Chi2 fit M83 IRS obs (15 pix * 15 pix) [done]'//NEW_LINE('')

  print*, '=================================================='
  
  !!------------------------------------------------------------------------
  !!                                 Analysis
  !!------------------------------------------------------------------------

  !! Calculate model
  !!-----------------
  ALLOCATE(FnuMOD(Nx,Ny,NwOBS))

  FnuMOD(:,:,:) = specModel(wOBS(:), INDPAR=ind, PARVAL=par(:,:,:), QABS=Qabs(:), &
                            FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
                            PABS=Pabs, FNULINE=FnuLINE)
  
  CALL WRITE_HDF5(DBLARR3D=FnuCONT*Pabs, FILE=filOUT, &
                  NAME="FnuCONT (MKS)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=FnuLINE+(FnuCONT+FnuSTAR)*Pabs, FILE=filOUT, &
                  NAME="FnuLINE (MKS)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=(FnuBAND+FnuCONT+FnuSTAR)*Pabs, FILE=filOUT, &
                  NAME="FnuBAND (MKS)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=FnuSTAR*Pabs, FILE=filOUT, &
                  NAME="FnuSTAR (MKS)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=FnuMOD, FILE=filOUT, &
                  NAME="FnuMOD (MKS)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  
  print*, 'fitChi2 analysis [done]'//NEW_LINE('')

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

END PROGRAM test_fitHB
