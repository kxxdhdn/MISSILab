!******************************************************************************
!*
!*                               HB Fitting Module
!*
!******************************************************************************


MODULE hb

  USE utilities, ONLY: DP
  USE inout, ONLY: lenpar
  USE core, ONLY: parinfo_type, indpar_type, Qabs_type
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: Df_prior = 8
  REAL(DP), PARAMETER, PUBLIC :: siglnS0 = 10._DP

  INTEGER, SAVE, PUBLIC :: Nx, Ny, NwOBS, Nparhyp, Ncorrhyp
  INTEGER, SAVE, PUBLIC :: xOBS, yOBS, ipar, ihpar, icorr
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: i2ih
  INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC :: icorr2ij
  REAL(DP), SAVE, PUBLIC :: detcov_prev, detcorr
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: cenlnS0
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wOBS, nuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parcurr, parhypcurr
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: mucurr, corrcurr, sigcurr
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC :: extinct
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC :: cov_prev, invcov_prev ! S-1
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC :: invcorr_prev ! R^-1
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC :: FnuOBS, dFnuOBS
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC :: allparhypcurr
  LOGICAL, SAVE, PUBLIC :: noposdef_prev
  LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: maskhypcurr
  LOGICAL, DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC :: maskS
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC :: mask, maskextra
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC :: maskhyp, maskpar

  !! Init param
  !!------------
  TYPE(indpar_type), SAVE, PUBLIC                             :: ind
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parinfo, parhypinfo
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC    :: Qabs
  
  PUBLIC :: lnpost_par, lnlhobs_par, lnpost_mu, lnpost_sig
  PUBLIC :: lnpost_corr, covariance

CONTAINS

  !! Build the covariance matrix
  !!-----------------------------
  PURE FUNCTION covariance(sig, corr)

    USE utilities, ONLY: DP
    USE statistics, ONLY: corr2Rmat
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(Nparhyp)  :: sig
    REAL(DP), INTENT(IN), DIMENSION(Ncorrhyp) :: corr
    REAL(DP), DIMENSION(Nparhyp,Nparhyp)      :: covariance

    REAL(DP), DIMENSION(Nparhyp,Nparhyp)      :: Smat, Rmat

    Smat(:,:) = UNPACK(sig(:), maskS(:,:), FIELD=0._DP)
    Rmat(:,:) = CORR2RMAT(corr(:), Nparhyp) ! correlation matrix
    covariance = MATMUL( MATMUL(Smat(:,:), Rmat(:,:)), Smat(:,:) )

  END FUNCTION covariance

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

    USE utilities, ONLY: DP
    USE core, ONLY: specModel
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: pargrid
    REAL(DP), DIMENSION(SIZE(pargrid)) :: lnLHobs_par

    INTEGER :: iw
    REAL(DP), DIMENSION(SIZE(pargrid),NwOBS) :: Fnu_mod, varred, lnpind

    !! Model
    Fnu_mod(:,:) = specModel(wOBS(:), PARVEC=pargrid(:), PARNAME=parinfo(ipar)%name, &
                             PARINFO=parinfo(:), INDPAR=ind, PARVAL=parcurr(:), &
                             QABS=Qabs(:), EXTINCT=extinct(:,:))

    !! Likelihoods
    varred(:,:) = 0._DP
    FORALL (iw=1:NwOBS) &
      varred(:,iw) = ( FnuOBS(xOBS,yOBS,iw) - Fnu_mod(:,iw) ) / dFnuOBS(xOBS,yOBS,iw)

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
  FUNCTION lnprior_par(pargrid)

    USE utilities, ONLY: DP
    USE arrays, ONLY: iwhere
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: pargrid
    REAL(DP), DIMENSION(SIZE(pargrid)) :: lnprior_par

    INTEGER :: ip, i, ih, Ngrid
    REAL(DP), DIMENSION(Nparhyp,SIZE(pargrid)) :: allpargrid

    !! Cube of parameters
    Ngrid = SIZE(pargrid(:))
    CALL IWHERE(i2ih(:) == ipar,ip)
    allpargrid(:,:) = 0._DP
    FORALL (ih=1:Nparhyp,maskhypcurr(ih)) &
      allpargrid(ih,:) = MERGE( SPREAD(parhypcurr(ih),DIM=1,NCOPIES=Ngrid), &
                                pargrid(:), (ip /= ih) ) - mucurr(ih)
    
    !! Hyperdistribution
    FORALL (i=1:Ngrid) &
      lnprior_par(i) = LOG( DOT_PRODUCT(allpargrid(:,i), &
                              MATMUL(invcov_prev(:,:), allpargrid(:,i))) &
                            / Df_prior + 1._DP ) * (-(Df_prior+Nparhyp)/2._DP)

  END FUNCTION lnprior_par

  !! For the variance
  !!------------------
  FUNCTION lnprior_sig(lnSgrid)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: lnSgrid
    REAL(DP), DIMENSION(SIZE(lnSgrid)) :: lnprior_sig
    
    lnprior_sig(:) = - 0.5_DP * ( (lnSgrid(:)-cenlnS0(ihpar)) / siglnS0 )**2

  END FUNCTION lnprior_sig

  !! For the correlation coefficient
  !!---------------------------------
  FUNCTION lnprior_corr (corrgrid)

    USE utilities, ONLY: DP, hugeDP
    USE matrices, ONLY: determinant_cholesky
    USE statistics, ONLY: corr2Rmat
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: corrgrid
    REAL(DP), DIMENSION(SIZE(corrgrid)) :: lnprior_corr

    INTEGER :: i, ip, Df_R, Ngrid
    REAL(DP), DIMENSION(Ncorrhyp) :: allcorr
    REAL(DP), DIMENSION(Nparhyp,SIZE(corrgrid)) :: lndetsubmat
    REAL(DP), DIMENSION(Nparhyp,Nparhyp) :: Rmatgrid
    LOGICAL, DIMENSION(Nparhyp,SIZE(corrgrid)) :: noposdef

    ! Grid of correlation matrices
    Df_R = Nparhyp + 1
    lndetsubmat(1,:) = 0._DP
    Ngrid = SIZE(corrgrid(:))

    DO i=1,Ngrid
      allcorr(:) = corrcurr(:)
      allcorr(icorr) = corrgrid(i)
      Rmatgrid(:,:) = CORR2RMAT(allcorr(:),Nparhyp)
      DO ip=2,Nparhyp
        lndetsubmat(ip,i) = LOG(DETERMINANT_CHOLESKY( Rmatgrid(1:ip,1:ip), &
                                                      NOPOSDEF=noposdef(ip,i) ))
      END DO
    END DO

    ! Distribution
    WHERE (ALL(.NOT. noposdef(2:Nparhyp,:),DIM=1))
      lnprior_corr(:) = ( 0.5_DP*(Df_R-1)*(Nparhyp-1) - 1._DP ) &
                      * lndetsubmat(Nparhyp,:) &
                        - 0.5_DP*Df_R * SUM(lndetsubmat(:,:),DIM=1)
    ELSEWHERE
      lnprior_corr(:) = - hugeDP
    END WHERE

  END FUNCTION lnprior_corr

  !!-------------------------------------------------------
  !!
  !!             Hyperparameter distributions
  !!
  !!-------------------------------------------------------

  !! For sampling the variance
  !!---------------------------
  FUNCTION lnhyper_sig (lnSgrid)

    USE utilities, ONLY: DP, hugeDP
    ! USE matrices, ONLY: invert_cholesky
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: lnSgrid
    REAL(DP), DIMENSION(SIZE(lnSgrid)) :: lnhyper_sig

    INTEGER :: ip, ix, iy, i, Ngrid
    REAL(DP), DIMENSION(Nparhyp,SIZE(lnSgrid)) :: allinvSgrid
    ! REAL(DP), DIMENSION(Nparhyp,SIZE(lnSgrid)) :: allSgrid
    ! REAL(DP), DIMENSION(Nparhyp,Nparhyp,SIZE(lnSgrid)) :: covar
    REAL(DP), DIMENSION(Nx,Ny,SIZE(lnSgrid)) :: lnpind
    REAL(DP), DIMENSION(Nx,Ny,Nparhyp) :: allpar0
    REAL(DP), DIMENSION(Nparhyp,Nparhyp,SIZE(lnSgrid)) :: invcov, allinvSmat
    REAL(DP), DIMENSION(SIZE(lnSgrid)) :: detcov
    ! LOGICAL, DIMENSION(SIZE(lnSgrid)) :: noposdef
    
    ! Covariance matrices and inverses for the grid
    Ngrid = SIZE(lnSgrid(:))

    ! Calculate invcov and detcov (opt.1: Cholesky)
    ! FORALL (ip=1:Nparhyp) &
    !   allSgrid(ip,:) = MERGE( SPREAD(sigcurr(ip),DIM=1,NCOPIES=Ngrid), &
    !                           EXP(lnSgrid(:)), (ip /= ihpar) )
    ! FORALL (i=1:Ngrid) &
    !   covar(:,:,i) = COVARIANCE( allSgrid(:,i), corrcurr(:) )
    ! DO i=1,Ngrid
    !   invcov(:,:,i) = INVERT_CHOLESKY( covar(:,:,i), DETERMINANT=detcov(i), &
    !                                    NOPOSDEF=noposdef(i) )
    ! END DO
    
    !! Calculate invcov and detcov (opt.2: invSmat+invcorr)
    FORALL (ip=1:Nparhyp) &
      allinvSgrid(ip,:) = MERGE( SPREAD(1._DP/sigcurr(ip),DIM=1,NCOPIES=Ngrid), &
                                 1._DP/EXP(lnSgrid(:)), (ip /= ihpar) )
    allinvSmat(:,:,:) = 0._DP
    FORALL (ip=1:Nparhyp) &
      allinvSmat(ip,ip,:) = allinvSgrid(ip,:)

    FORALL (i=1:Ngrid)
      invcov(:,:,i) = MATMUL( MATMUL( allinvSmat(:,:,i),invcorr_prev(:,:) ), &
                              allinvSmat(:,:,i) )
      detcov(i) = detcorr * 1._DP/PRODUCT( allinvSgrid(:,i) )**2
    END FORALL

    ! Individual distributions
    allpar0(:,:,:) = 0._DP
    FORALL (ix=1:Nx,iy=1:Ny,ip=1:Nparhyp,maskhyp(ix,iy,ip)) &
      allpar0(ix,iy,ip) = allparhypcurr(ix,iy,ip) - mucurr(ip)
    FORALL (ix=1:Nx,iy=1:Ny,i=1:Ngrid)
      lnpind(ix,iy,i) = MERGE( &
        LOG( DOT_PRODUCT(allpar0(ix,iy,:), &
                         MATMUL(invcov(:,:,i),allpar0(ix,iy,:))) &
             / Df_prior + 1._DP ) * (-(Df_prior+Nparhyp)/2._DP) &
        - 0.5_DP * LOG(detcov(i)), - hugeDP, (.NOT. noposdef_prev) )
    END FORALL

    ! Total distribution
    lnhyper_sig(:) = 0._DP
    FORALL (i=1:Ngrid) &
      lnhyper_sig(i) = SUM(lnpind(:,:,i),MASK=maskhyp(:,:,ihpar))

  END FUNCTION lnhyper_sig

  !! For sampling the correlation
  !!------------------------------
  FUNCTION lnhyper_corr (corrgrid)

    USE utilities, ONLY: DP, hugeDP!, &
                         ! pring, trimLR, timinfo, initiate_clock, time_type
    USE matrices, ONLY: invert_cholesky, determinant_cholesky
    USE core, ONLY: invert_mSM
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: corrgrid
    REAL(DP), DIMENSION(SIZE(corrgrid)) :: lnhyper_corr

    INTEGER :: ix, iy, ipar, jpar, i, Ngrid!, ic
    REAL(DP), DIMENSION(SIZE(corrgrid)) :: delta
    ! REAL(DP), DIMENSION(Ncorrhyp,SIZE(corrgrid)) :: allcorrgrid
    REAL(DP), DIMENSION(Nx,Ny,SIZE(corrgrid)) :: lnpind
    REAL(DP), DIMENSION(Nx,Ny,Nparhyp) :: allpar0
    REAL(DP), DIMENSION(Nparhyp,Nparhyp,SIZE(corrgrid)) :: invcov!, covar
    REAL(DP), DIMENSION(SIZE(corrgrid)) :: detcov
    LOGICAL, DIMENSION(SIZE(corrgrid)) :: noposdef

    ! TYPE(time_type) :: timestr
    ! CALL INITIATE_CLOCK(timestr)
    Ngrid = SIZE(corrgrid(:))

    !! Calculate invcov and detcov (opt.1: Cholesky)
    !!-----------------------------------------------
    !! Covariance matrices and inverses for the grid
    ! FORALL (ic=1:Ncorrhyp) &
    !   allcorrgrid(ic,:) = MERGE( SPREAD(corrcurr(ic),DIM=1,NCOPIES=Ngrid), &
    !                              corrgrid(:), (ic /= icorr) )
    ! FORALL (i=1:Ngrid) &
    !   covar(:,:,i) = COVARIANCE( sigcurr(:), allcorrgrid(:,i) )

    ! DO i=1,Ngrid
    !   invcov(:,:,i) = INVERT_CHOLESKY( covar(:,:,i), DETERMINANT=detcov(i), &
    !                                    NOPOSDEF=noposdef(i) )
    ! END DO

    !! Calculate invcov and detcov (opt.2: modified Sherman-Morrison)
    !!----------------------------------------------------------------
    ipar = icorr2ij(icorr,1)
    jpar = icorr2ij(icorr,2)
    delta(:) = corrgrid(:)*sigcurr(ipar)*sigcurr(jpar) - cov_prev(ipar,jpar)
    DO i=1,Ngrid
      invcov(:,:,i) = INVERT_MSM( cov_prev(:,:), invcov_prev(:,:), &
                                  icorr2ij(icorr,:), delta(i), Nparhyp, &
                                  detcov_prev, detcov(i) )!, noposdef(i) )
    END DO

    FORALL (i=1:Ngrid) noposdef(i) = ANY( isNaN(invcov(:,:,i)) )
    ! PRINT*, 'icorr='//TRIMLR(pring(icorr))//' end: '//TRIMLR(TIMINFO(timestr))

    !! Individual distributions
    allpar0(:,:,:) = 0._DP
    FORALL (ix=1:Nx,iy=1:Ny,i=1:Nparhyp,maskhyp(ix,iy,i)) &
      allpar0(ix,iy,i) = allparhypcurr(ix,iy,i) - mucurr(i)
    FORALL (ix=1:Nx,iy=1:Ny,i=1:Ngrid)
      lnpind(ix,iy,i) = MERGE( &
        LOG( DOT_PRODUCT(allpar0(ix,iy,:), &
                         MATMUL(invcov(:,:,i),allpar0(ix,iy,:))) &
             / Df_prior + 1._DP ) * (-(Df_prior+Nparhyp)/2._DP) &
        - 0.5_DP * LOG(detcov(i)), - hugeDP, (.NOT. noposdef(i)) )
    END FORALL
    
    !! Total distribution
    lnhyper_corr(:) = 0._DP
    FORALL (i=1:Ngrid) &
      lnhyper_corr(i) = SUM(lnpind(:,:,i),MASK=(maskhyp(:,:,icorr2ij(icorr,1)) &
                                          .AND. maskhyp(:,:,icorr2ij(icorr,2))))
    
  END FUNCTION lnhyper_corr

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

    IF (parinfo(ipar)%hyper) THEN
      lnpost_par(:) = LNLHOBS_PAR(pargrid(:)) + LNPRIOR_PAR(pargrid(:))
    ELSE
      lnpost_par(:) = LNLHOBS_PAR(pargrid(:))
    END IF

  END FUNCTION lnpost_par

  !! For sampling the mean
  !!-----------------------
  FUNCTION lnpost_mu (mugrid)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: mugrid
    REAL(DP), DIMENSION(SIZE(mugrid)) :: lnpost_mu

    INTEGER :: ip, ix, iy, i, ih, Ngrid
    REAL(DP), DIMENSION(Nparhyp,SIZE(mugrid)) :: allmugrid
    REAL(DP), DIMENSION(Nx,Ny,Nparhyp,SIZE(mugrid)) :: allmugrid0
    REAL(DP), DIMENSION(Nx,Ny,SIZE(mugrid)) :: lnpind

    ! Grid of parameters
    Ngrid = SIZE(mugrid(:))
    FORALL (ih=1:Nparhyp) &
      allmugrid(ih,:) = MERGE( SPREAD(mucurr(ih),DIM=1,NCOPIES=Ngrid), &
                               mugrid(:), (ih /= ihpar) )

    ! Individual distributions
    allmugrid0(:,:,:,:) = 0._DP
    FORALL (ix=1:Nx,iy=1:Ny,ip=1:Nparhyp,maskhyp(ix,iy,ip)) &
      allmugrid0(ix,iy,ip,:) = allparhypcurr(ix,iy,ip) - allmugrid(ip,:)
    FORALL (ix=1:Nx,iy=1:Ny,i=1:Ngrid)
      lnpind(ix,iy,i) = LOG( DOT_PRODUCT(allmugrid0(ix,iy,:,i), &
                             MATMUL(invcov_prev(:,:),allmugrid0(ix,iy,:,i))) &
                           / Df_prior + 1._DP ) * (-(Df_prior+Nparhyp)/2._DP)
    END FORALL

    ! Total distribution
    lnpost_mu(:) = 0._DP
    FORALL (i=1:Ngrid) &
      lnpost_mu(i) = SUM(lnpind(:,:,i),MASK=maskhyp(:,:,ihpar))

  END FUNCTION lnpost_mu

  !! For sampling the variance
  !!---------------------------
  FUNCTION lnpost_sig(lnSgrid)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: lnSgrid
    REAL(DP), DIMENSION(SIZE(lnSgrid)) :: lnpost_sig

    lnpost_sig(:) = LNHYPER_SIG(lnSgrid(:)) + LNPRIOR_SIG(lnSgrid(:))

  END FUNCTION lnpost_sig

  !! For sampling the correlation
  !!------------------------------
  FUNCTION lnpost_corr(corrgrid)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN)      :: corrgrid
    REAL(DP), DIMENSION(SIZE(corrgrid(:)))  :: lnpost_corr

    lnpost_corr(:) = LNHYPER_CORR(corrgrid(:)) + LNPRIOR_CORR(corrgrid(:))

  END FUNCTION lnpost_corr

END MODULE hb
