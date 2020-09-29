MODULE fitHB_external

  USE auxil, ONLY: parinfo_type, indpar_type, Qabs_type
  USE utilities, ONLY: DP
  USE inout, ONLY: lenpar
  USE grain_optics, ONLY: lendustQ
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: Df_prior = 8

  INTEGER, SAVE, PUBLIC :: Nx, Ny, NwOBS, Nw
  INTEGER, SAVE, PUBLIC :: xOBS, yOBS, jw, ipar, ihpar
  INTEGER, SAVE, PUBLIC :: Nparhyp, Ncorrhyp
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wOBS, nuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parcurr, parhypcurr
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: FnuOBS, dFnuOBS
  ! LOGICAL, DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC    :: maskS
  ! LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC  :: maskhyp, maskpar

  !! Init param
  !!------------
  TYPE(indpar_type), SAVE, PUBLIC                             :: indpar
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
     !                         PARINFO=parinfo(:), INDPAR=indpar, PARVAL=parcurr(:), &
     !                         QABS=Qabs, INDBIN=indBIN, INDLIN=indLIN),1)
     ! print*, SIZE(specModel(wOBS(:), PARVEC=pargrid(:), PARNAME=parinfo(ipar)%name, &
     !                         PARINFO=parinfo(:), INDPAR=indpar, PARVAL=parcurr(:), &
     !                         QABS=Qabs, INDBIN=indBIN, INDLIN=indLIN),2)
     ! print*, shape(fnu_mod)
    Fnu_mod(:,:) = specModel(wOBS(:), PARVEC=pargrid(:), PARNAME=parinfo(ipar)%name, &
                             PARINFO=parinfo(:), INDPAR=indpar, PARVAL=parcurr(:), &
                             QABS=Qabs, INDBIN=indBIN, INDLIN=indLIN)

    !! Likelihoods
    varred(:,:) = 0._DP
    FORALL (iw=1:Nw) &
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
  USE fitHB_external, ONLY: NwOBS, wOBS, nuOBS, &!Nx, Ny, xOBS, yOBS, &
                            FnuOBS, dFnuOBS, indpar, &!Nw, &
                            ! Nparhyp, &
                            ipar, &!ihpar, &
                            parcurr, &!parhypcurr, maskS, &
                            lnpost_par, lnlhobs_par, &!lnpost_mu, lnpost_sig, lnpost_corr, &
                            ! covariance, maskhyp, maskpar, &
                            parinfo, Qabs, indBIN, indLIN
  IMPLICIT NONE

  INTEGER :: Ncont, Nband, Nline

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: Ngibbsmax = 3000
  REAL(DP), PARAMETER :: accrand = 1.E-3_DP
  CHARACTER(*), PARAMETER :: nampar = "Parameter values"
  CHARACTER(*), PARAMETER :: namln1pd = "ln(1+delta)"
  CHARACTER(*), PARAMETER :: nammu = "Mean of hyperdistribution"
  CHARACTER(*), PARAMETER :: namsig = "Sigma of hyperdistribution"
  CHARACTER(*), PARAMETER :: namcorr = "Correlation of hyperdistribution"
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Gen spec
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labBIN, labLIN
  INTEGER, PARAMETER :: Ngrid=1000
  REAL(DP), DIMENSION(:), ALLOCATABLE :: pargen
  
  !! Input variables
  INTEGER :: i0, i!, x, y
  INTEGER :: Npar, Nmcmc!, Nextra, NiniMC, Nsou
  ! INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: maskint
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: parini, parmcmc
  ! REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parini, parmcmc
  ! REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mumcmc, sigmcmc
  LOGICAL :: verbose, newseed
  ! LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: maskxy3D
  CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ
  !! MCMC parameters
  INTEGER :: icurr, iprev, counter!, printevery
  ! INTEGER, DIMENSION(:), ALLOCATABLE :: ifirst, ilast
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
  REAL(DP), DIMENSION(:), ALLOCATABLE :: FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  CHARACTER(lenpar) :: filANAL

  !! Output settings
  CALL INITIATE_CLOCK(timestr)
  unitlog(:) = [ ulog, ustd ]

  !!-----------------
  !! Read the inputs
  !!-----------------
  labQ = (/'ACH2_Z96             ', &
           'BE_Z96               ', &
           'Sil_D03              '/)
  labBIN = (/'Main 3.3     ', & ! 1
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
  labLIN = (/'H2S7  ', & ! 3
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
  Nmcmc = 10 ! Length of Markov chain Monte Carlo

  IF (newseed) CALL GENERATE_NEWSEED()

  !!----------------
  !! Generate param
  !!----------------

  !! pargen = [Mcont [Msun/pc2], Tcont [K], &
  pargen = [2.E-2_DP, 100._DP, 5.E-2_DP, 100._DP, 3.E-2_DP, 100._DP, & 
  !!        Iline, Cline, Wline, &
            1.E-9_DP, LIN(3)%wave, degradeRes(LIN(3)%wave,.01_DP,'SL-LL'), & 
            2.5E-9_DP, LIN(9)%wave, degradeRes(LIN(9)%wave,.01_DP,'SL-LL'), & 
            3.E-9_DP, LIN(11)%wave, degradeRes(LIN(11)%wave,.01_DP,'SL-LL'), &
            1.5E-9_DP, LIN(20)%wave, degradeRes(LIN(20)%wave,.01_DP,'SL-LL'), &
            1.5E-9_DP, LIN(23)%wave, degradeRes(LIN(23)%wave,.01_DP,'SL-LL'), &
            2.5E-9_DP, LIN(24)%wave, degradeRes(LIN(24)%wave,.01_DP,'SL-LL'), &
            1.5E-9_DP, LIN(26)%wave, degradeRes(LIN(26)%wave,.01_DP,'SL-LL'), &
            3.E-9_DP, LIN(27)%wave, degradeRes(LIN(27)%wave,.01_DP,'SL-LL'), &
            3.5E-9_DP, LIN(29)%wave, degradeRes(LIN(29)%wave,.01_DP,'SL-LL'), &
            4.5E-9_DP, LIN(33)%wave, degradeRes(LIN(33)%wave,.01_DP,'SL-LL'), &
            1.5E-9_DP, LIN(34)%wave, degradeRes(LIN(34)%wave,.01_DP,'SL-LL'), &
            1.E-9_DP, LIN(36)%wave, degradeRes(LIN(36)%wave,.01_DP,'SL-LL'), &
  !!        Iband, Cband, WSband, WLband, &
            1.E-9_DP, BIN(1)%wave, BIN(1)%sigmaS, BIN(1)%sigmaL, &
            .8E-9_DP, BIN(2)%wave, BIN(2)%sigmaS, BIN(2)%sigmaL, &
            .5E-9_DP, BIN(7)%wave, BIN(7)%sigmaS, BIN(7)%sigmaL, &
            .7E-9_DP, BIN(8)%wave, BIN(8)%sigmaS, BIN(8)%sigmaL, &
            .3E-9_DP, BIN(12)%wave, BIN(12)%sigmaS, BIN(12)%sigmaL, &
            .5E-9_DP, BIN(13)%wave, BIN(13)%sigmaS, BIN(13)%sigmaL, &
            .4E-9_DP, BIN(14)%wave, BIN(14)%sigmaS, BIN(14)%sigmaL, &
            1.E-9_DP, BIN(16)%wave, BIN(16)%sigmaS, BIN(16)%sigmaL, &
            .3E-9_DP, BIN(19)%wave, BIN(19)%sigmaS, BIN(19)%sigmaL, &
            .7E-9_DP, BIN(20)%wave, BIN(20)%sigmaS, BIN(20)%sigmaL, &
            .6E-9_DP, BIN(21)%wave, BIN(21)%sigmaS, BIN(21)%sigmaL, &
            .5E-9_DP, BIN(24)%wave, BIN(24)%sigmaS, BIN(24)%sigmaL, &
            .5E-9_DP, BIN(25)%wave, BIN(25)%sigmaS, BIN(25)%sigmaL, &
            .9E-9_DP, BIN(30)%wave, BIN(30)%sigmaS, BIN(30)%sigmaL, &
  !!        Av [mag], &
            0._DP, &
  !!        Fstar [Lsun/pc2]]
            2.E4_DP]

  !!-------------
  !! Create grid
  !!-------------
  wOBS = RAMP(Ngrid, 1._DP, 40._DP, XLOG=.TRUE.)
  NwOBS = SIZE(wOBS(:))
  nuOBS = MKS%clight/MKS%micron / wOBS(:)

  Ncont = SIZE(labQ(:))
  Nband = SIZE(labBIN(:))
  Nline = SIZE(labLIN(:))
  Npar = 2*Ncont + 4*Nband + 3*Nline + 2
  
  ALLOCATE(parinfo(Npar))
  i0 = 0
  CONTinfo: DO i=1,Ncont
    parinfo(i0+2*i-1)%name = "Mcont"//TRIMLR(PRING(i))
    parinfo(i0+2*i-1)%comp = "CONT"
    parinfo(i0+2*i-1)%limits = [0._DP, 1._DP]
    parinfo(i0+2*i-1)%limited = [.TRUE., .FALSE.]
    parinfo(i0+2*i-1)%fixed = .FALSE.
    parinfo(i0+2*i)%name = "Tcont"//TRIMLR(PRING(i))
    parinfo(i0+2*i)%comp = "CONT"
    parinfo(i0+2*i)%limits = [50._DP, 500._DP]
    parinfo(i0+2*i)%limited = [.TRUE., .TRUE.]
    parinfo(i0+2*i)%fixed = .FALSE.
  END DO CONTinfo
  i0 = i0 + 2*Ncont
  LINEinfo: DO i=1,Nline
    parinfo(i0+3*i-2)%name = "Iline"//TRIMLR(PRING(i))
    parinfo(i0+3*i-2)%comp = "LINE"
    parinfo(i0+3*i-2)%limits = [0._DP, 1._DP]
    parinfo(i0+3*i-2)%limited = [.TRUE., .FALSE.]
    parinfo(i0+3*i-1)%name = "Cline"//TRIMLR(PRING(i))
    parinfo(i0+3*i-1)%comp = "LINE"
    parinfo(i0+3*i)%name = "Wline"//TRIMLR(PRING(i))
    parinfo(i0+3*i)%comp = "LINE"
    parinfo(i0+3*i)%fixed = .TRUE.
  END DO LINEinfo
  i0 = i0 + 3*Nline
  BANDinfo: DO i=1,Nband
    parinfo(i0+4*i-3)%name = "Iband"//TRIMLR(PRING(i))
    parinfo(i0+4*i-3)%comp = "BAND"
    parinfo(i0+4*i-3)%limits = [0._DP,1._DP]
    parinfo(i0+4*i-3)%limited = [.TRUE.,.FALSE.]
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
  parinfo(i0+1)%name = "Av"
  parinfo(i0+1)%comp = "PABS"
  parinfo(i0+1)%fixed = .TRUE.
  i0 = i0 + 1
  parinfo(i0+1)%name = "Fstar"
  parinfo(i0+1)%comp = "STAR"
  parinfo(i0+1)%limits = [0._DP,1._DP]
  parinfo(i0+1)%limited = [.TRUE.,.FALSE.]

  parinfo(:)%ind = [(i,i=1,Npar)]
  CALL set_indpar(indpar, parinfo(:))

  CALL READ_MASTER(LABQ=labQ(:), &
                   LABBIN=labBIN(:), LABLIN=labLIN(:), &
                   WAVEALL=wOBS(:), &
                   QABS=Qabs, INDBIN=indBIN, INDLIN=indLIN)

  !! Build synthetic spectrum
  ALLOCATE(FnuOBS(NwOBS), dFnuOBS(NwOBS))

  FnuOBS(:) = specModel(wOBS(:), INDPAR=indpar, PARVAL=pargen(:), QABS=Qabs(:))
  dFnuOBS(:) = RAND_NORM(NwOBS) * 0.01_DP*MAXVAL(FnuOBS(:))
  FnuOBS(:) = FnuOBS(:) + dFnuOBS(:)

  filOUT = 'out/test_fitHBsyn'
  CALL write_hdf5(wOBS, NAME="Wavelength (microns)", &
                  COMPRESS=.False., APPEND=.False., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(dFnuOBS, NAME="FnuUNC (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(FnuOBS, NAME="FnuOBS (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)

  print*, 'Gen synthetic spectrum [done]'

  print*, '=================================================='
  
  !! Observations
  ! filOBS = "dat/M83"//h5ext
  ! CALL READ_HDF5(wOBS, FILE=filOBS, NAME='Wavelength (microns)')
  ! CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, NAME='FnuOBS (MJyovsr)')
  ! CALL READ_HDF5(DBLARR3D=dFnuOBS, FILE=filOBS, NAME='dFnuOBS (MJyovsr)')
  ! CALL READ_HDF5(INTARR3D=maskint, FILE=filOBS, NAME='Mask')

  !! 3D parameter mask
  ! ALLOCATE (maskpar(Nx,Ny,Npar))
  ! maskpar(:,:,:) = .TRUE.

  !!------------------
  !! Prepare the MCMC
  !!------------------

  !! Initialize the parameters
  !!---------------------------
  ALLOCATE(parini(Npar,1))
  parini(:,:) = 0._DP
  ! ALLOCATE(parini(Nx,Ny,Npar,1))
  ! parini(:,:,:,:) = 0._DP
  
  !! Automatic initial parameter estimates (via Chi2 result)
  ! CALL INITPARAM()

  !! Alternative manual INITPARAM
  !! parini(:,1) = [ Mcont[Msun/pc2], Tcont [K], &
  parini(:,1) = [1.E-2_DP, 50._DP, 1.E-2_DP, 50._DP, 1.E-2_DP, 50._DP, & 
  !!             Iline, Cline, Wline, &
                 1.E-9_DP, LIN(3)%wave, degradeRes(LIN(3)%wave, .01_DP, 'SL-LL'), & 
                 1.E-9_DP, LIN(9)%wave, degradeRes(LIN(9)%wave, .01_DP, 'SL-LL'), & 
                 1.E-9_DP, LIN(11)%wave, degradeRes(LIN(11)%wave, .01_DP, 'SL-LL'), &
                 1.E-9_DP, LIN(20)%wave, degradeRes(LIN(20)%wave, .01_DP, 'SL-LL'), &
                 1.E-9_DP, LIN(23)%wave, degradeRes(LIN(23)%wave, .01_DP, 'SL-LL'), &
                 1.E-9_DP, LIN(24)%wave, degradeRes(LIN(24)%wave, .01_DP, 'SL-LL'), &
                 1.E-9_DP, LIN(26)%wave, degradeRes(LIN(26)%wave, .01_DP, 'SL-LL'), &
                 1.E-9_DP, LIN(27)%wave, degradeRes(LIN(27)%wave, .01_DP, 'SL-LL'), &
                 1.E-9_DP, LIN(29)%wave, degradeRes(LIN(29)%wave, .01_DP, 'SL-LL'), &
                 1.E-9_DP, LIN(33)%wave, degradeRes(LIN(33)%wave, .01_DP, 'SL-LL'), &
                 1.E-9_DP, LIN(34)%wave, degradeRes(LIN(34)%wave, .01_DP, 'SL-LL'), &
                 1.E-9_DP, LIN(36)%wave, degradeRes(LIN(36)%wave, .01_DP, 'SL-LL'), &
  !!             Iband, Cband, WSband, WLband, &
                 1.E-9_DP, BIN(1)%wave, BIN(1)%sigmaS, BIN(1)%sigmaL, & 
                 1.E-9_DP, BIN(2)%wave, BIN(2)%sigmaS, BIN(2)%sigmaL, & 
                 1.E-9_DP, BIN(7)%wave, BIN(7)%sigmaS, BIN(7)%sigmaL, & 
                 1.E-9_DP, BIN(8)%wave, BIN(8)%sigmaS, BIN(8)%sigmaL, & 
                 1.E-9_DP, BIN(12)%wave, BIN(12)%sigmaS, BIN(12)%sigmaL, & 
                 1.E-9_DP, BIN(13)%wave, BIN(13)%sigmaS, BIN(13)%sigmaL, & 
                 1.E-9_DP, BIN(14)%wave, BIN(14)%sigmaS, BIN(14)%sigmaL, & 
                 1.E-9_DP, BIN(16)%wave, BIN(16)%sigmaS, BIN(16)%sigmaL, & 
                 1.E-9_DP, BIN(19)%wave, BIN(19)%sigmaS, BIN(19)%sigmaL, & 
                 1.E-9_DP, BIN(20)%wave, BIN(20)%sigmaS, BIN(20)%sigmaL, &
                 1.E-9_DP, BIN(21)%wave, BIN(21)%sigmaS, BIN(21)%sigmaL, & 
                 1.E-9_DP, BIN(24)%wave, BIN(24)%sigmaS, BIN(24)%sigmaL, & 
                 1.E-9_DP, BIN(25)%wave, BIN(25)%sigmaS, BIN(25)%sigmaL, & 
                 1.E-9_DP, BIN(30)%wave, BIN(30)%sigmaS, BIN(30)%sigmaL, & 
  !!             Av [mag], &
                 0._DP, & 
  !!             Fstar [Lsun/pc2]]
                 1.E4_DP]

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
  ! ALLOCATE(parmcmc(Nx,Ny,Npar,2))
  DO i=1,2!!
    parmcmc(:,i) = parini(:,1)
    ! FORALL(x=1:Nx,y=1:Ny) parmcmc(x,y,:,1) = parini(x,y,:,1)
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
  filMCMC = 'out/parlog_fitHBsyn'
  CALL WRITE_HDF5(INTARR1D=[Nmcmc], FILE=TRIMLR(filMCMC)//h5ext, &
                  NAME="Length of MCMC", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.FALSE.)
  CALL WRITE_HDF5(INITDBLARR=[Npar,Nmcmc], FILE=TRIMLR(filMCMC)//h5ext, &
                  NAME=nampar, COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  ! CALL WRITE_HDF5(INITDBLARR=[Nparhyp,Nmcmc], FILE=TRIMLR(filMCMC)//h5ext, &
  !                 NAME=nammu, COMPRESS=compress, VERBOSE=debug, &
  !                 APPEND=.TRUE.)
  ! CALL WRITE_HDF5(INITDBLARR=[Nparhyp,Nmcmc], FILE=TRIMLR(filMCMC)//h5ext, &
  !                 NAME=namsig, COMPRESS=compress, VERBOSE=debug, &
  !                 APPEND=.TRUE.)
  
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

  !! Frequency of printing
  ! printevery = MAX( NINT(1000._DP/(REAL(Nsou,DP)*Npar))*10, 1 )
  ! IF (debug) printevery = 1

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
        ! xsource: DO xobs=1,Nx
          ! ysource: DO yobs=1,Ny
            ! IF (maskpar(xobs,yobs,ipar)) THEN
              parcurr(:) = parmcmc(:,icurr)
  ! print*, counter, ipar
  ! print*, MAXVAL(specModel(wOBS(:), INDPAR=indpar, PARVAL=parcurr(:), QABS=Qabs(:)))
              ! parcurr(:) = parmcmc(xobs,yobs,:,icurr)
              parmcmc(ipar,icurr) &
              ! parmcmc(xobs,yobs,ipar,icurr) &
                = RAND_GENERAL(lnpost_par, lim(:), YLOG=.TRUE., VERBOSE=debug, &
                               LNFUNC=.TRUE., ACCURACY=accrand, NMAX=Ngibbsmax)
              IF (debug) PRINT*, "par(",TRIMLR(PRING(ipar)),") = ", parmcmc(ipar,icurr)
              ! IF (debug) PRINT*, "par(",xobs,yobs,ipar,") = ", &
                ! parmcmc(xobs,yobs,ipar,icurr)
              
            ! END IF
          ! END DO ysource
        ! END DO xsource
      END IF
    END DO param

    !! Hyperparameters
    !!-----------------

      !! a. Average
      !! b. Variance
      !! c. Correlation

    !! Outputs
    !!---------
    CALL WRITE_HDF5(DBLARR2D=parmcmc(:,icurr:icurr), FILE=TRIMLR(filMCMC)//h5ext, &
                    NAME=nampar, COMPRESS=compress, VERBOSE=debug, &
                    IND2=[counter,counter])
    ! CALL WRITE_HDF5(DBLARR2D=mumcmc(:,icurr:icurr), FILE=TRIMLR(filMCMC)//h5ext, &
    !                 NAME=nampar, COMPRESS=compress, VERBOSE=debug, &
    !                 IND2=[counter,counter])
    ! CALL WRITE_HDF5(DBLARR2D=sigmcmc(:,icurr:icurr), FILE=TRIMLR(filMCMC)//h5ext, &
    !                 NAME=nampar, COMPRESS=compress, VERBOSE=debug, &
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
    !   CALL STRIKE("FITHB", &
    !               "a NaN appeared in the MCMC at iteration "//TRIMLR(PRING(counter)))

    !! Normal termination condition
    IF (counter == Nmcmc) EXIT

    !! Update the indices
    CALL SWAP(icurr,iprev)
    counter = counter + 1
  
  END DO big_loop

  print*, 'HB fit synt spec [done]'

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
                  NAME="Mean of parameter value", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=stdevpar(:), FILE=TRIMLR(filANAL)//h5ext, &
                  NAME="Sigma of parameter value", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)

  !! Calculate model
  !!-----------------
  ALLOCATE(FnuMOD(NwOBS))
  
  FnuMOD(:) = specModel(wOBS(:), INDPAR=indpar, PARVAL=meanpar(:), QABS=Qabs(:), &
                        FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
                        PABS=Pabs, FNULINE=FnuLINE)
  
  CALL WRITE_HDF5(DBLARR1D=FnuCONT*Pabs, FILE=TRIMLR(filOUT)//h5ext, &
                  NAME="FnuCONT (MJyovsr)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=FnuLINE+(FnuCONT+FnuSTAR)*Pabs, FILE=TRIMLR(filOUT)//h5ext, &
                  NAME="FnuLINE (MJyovsr)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=(FnuBAND+FnuCONT+FnuSTAR)*Pabs, FILE=TRIMLR(filOUT)//h5ext, &
                  NAME="FnuBAND (MJyovsr)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=FnuSTAR*Pabs, FILE=TRIMLR(filOUT)//h5ext, &
                  NAME="FnuSTAR (MJyovsr)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=FnuMOD, FILE=TRIMLR(filOUT)//h5ext, &
                  NAME="FnuMOD (MJyovsr)", COMPRESS=compress, VERBOSE=debug, &
                  APPEND=.TRUE.)
  
  ! !! Final
  ! !!-------
  ! DO i=1,MERGE(2,1,verbose)
  !   WRITE(unitlog(i),*)
  !   WRITE(unitlog(i),*) "PROGRAM EXECUTED IN " &
  !                       //TRIMLR(TIMINFO(timestr))//"."
  !   WRITE(unitlog(i),*)
  ! END DO
  ! IF (verbose) THEN
  !   PRINT*, " - File "//filog//" has been written."
  !   PRINT*
  ! END IF

  print*, 'HB fit analysis [done]'

  print*, '=================================================='
    
END PROGRAM test_fitHB
