!******************************************************************************
!*
!*         SIMULATING MIDINFRARED SPECTRA WITH NOISE AND CORRELATIONS
!*
!******************************************************************************


  !==========================================================================
  ! 1) AUTHOR: D. HU
  ! 
  ! 2) DESCRIPTION: Generating synthetic spectra for testing MILES fitting
  !
  ! 3) HISTORY: 
  !    - 20210120: Created.
  !==========================================================================

  
PROGRAM simulate_MIR

  USE auxil, ONLY: parinfo_type, indpar_type, Qabs_type, read_master, specModel
  USE utilities, ONLY: DP, trimeq, trimlr, pring
  USE arrays, ONLY: reallocate, closest, iwhere, incrarr
  USE statistics, ONLY: median_data, mean, sigma, N_corr, &
                        corr2Rmat, correlate
  USE random, ONLY: rand_norm, rand_multinorm
  USE interpolation, ONLY: interp_lin_sorted
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, lenpar
  USE grain_optics, ONLY: lendustQ
  IMPLICIT NONE

  INTEGER, PARAMETER :: textwid = 60
  CHARACTER(*), PARAMETER :: dirIN = '../out1/'
  CHARACTER(*), PARAMETER :: filOBS = dirIN//'galgen'//h5ext
  CHARACTER(*), PARAMETER :: filGEN = dirIN//'pargen'//h5ext
  CHARACTER(*), PARAMETER :: filOUT = dirIN//'galspec'//h5ext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .TRUE.
  
  INTEGER :: i, p, q, iw!, icorr, ipar
  INTEGER :: Np, Nq, Nw, NwOBS, Npar, Ncorr, Nrat, Nbandr, Ncont
  INTEGER :: counter
  INTEGER, DIMENSION(:), ALLOCATABLE :: indB ! 1D (Nbandr)
  REAL(DP) :: wvl0!, lnIref, rescaling
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dblarr1d
  REAL(DP), DIMENSION(:), ALLOCATABLE :: mulnR, siglnR ! 1D (Nrat)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: corr ! 1D (Ncorr)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wOBS ! 1D (NwOBS)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wvl, normSovN ! 1D (Nw)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: multiCONT, SovN ! 1D (Nq)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dblarr2d
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: extinct ! 2D (Nextc, Nw)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: covar, Smat, Rmat ! 2D (Nrat,Nrat)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: lnR ! 3D (Nrat,Np,Nq)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: pargen ! 3D (1,1,Npar)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: par ! 3D (Np,Nq,Npar)
  ! REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Rband ! 3D (Np,Nq,Nbandr)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuSIM, dFnuSIM ! 3D (Np,Nq,Nw) [U_in]
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuOBS, dFnuOBS ! 3D (1,1,NwOBS) [U_in]
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: fnuband_tab
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: ratname ! 1D (Nrat)
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: corrname ! 1D (Ncorr)
  LOGICAL :: varCONT, varSovN, flatSovN, flatUNC
  TYPE(indpar_type) :: ind
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE :: parinfo
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE :: Qabs ! 1D (Ncont)
  
  !!---------------
  !! Preliminaries
  !!---------------

  !! Banner
  PRINT*,REPEAT("=",textwid)
  PRINT*, REPEAT(" ",textwid/3)//'Galaxy Spectrum Gallery'
  PRINT*,REPEAT("=",textwid)
  PRINT*

  !! Sample size
  Np = 6 ! sampling along correlation axis
  Nq = 9 ! sampling along option (FnuCONT or SovN) axis

  !! Options
  varCONT = .TRUE.
  varSovN = .TRUE.
  
  !! Wavelength reference
  wvl0 = 15._DP
  
  !! Opt 1. Continuum level
  !!------------------------
  ALLOCATE(multiCONT(Nq))
  
  IF (varCONT) THEN
    multiCONT(:) = [ .1_DP, 1._DP, 10._DP, &
                     .1_DP, 1._DP, 10._DP, &
                     .1_DP, 1._DP, 10._DP ]
  ELSE
    multiCONT(:) = 1._DP
  END IF

  !! Opt 2. S/N ratio
  !!------------------
  ALLOCATE(SovN(Nq))

  IF (varSovN) THEN
    SovN(:) = [ 5._DP, 5._DP, 5._DP, &
                25._DP, 25._DP, 25._DP, &
                125._DP, 125._DP, 125._DP ]
  ELSE
    SovN(:) = 10._DP
  END IF

  !! Read the inputs
  !!--------------------
  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='Wavelength (microns)', N1=NwOBS)
  !! Interpolate wavelength grid
  wvl = wOBS(:)
  Nw = SIZE(wvl)

  CALL READ_MASTER(WAVALL=wvl(:), DIRIN=dirIN, &
                   LABB=labB, &
                   NCONT=Ncont, &
                   QABS=Qabs, EXTINCT=extinct, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit)

  !! Spitzer/IRS obs of M83 (gen SovN)
  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, &
                 NAME='FnuOBS ('//TRIMLR(spec_unit)//')')
  CALL READ_HDF5(DBLARR3D=dFnuOBS, FILE=filOBS, &
                 NAME='dFnuOBS ('//TRIMLR(spec_unit)//')')

  !! Chi2 fitted par
  CALL READ_HDF5(DBLARR3D=pargen, FILE=filGEN, &
                 NAME='Best fitted parameter value')

  !! Init par
  ALLOCATE(par(Np,Nq,Npar))

  FORALL (p=1:Np,q=1:Nq) par(p,q,:) = pargen(1,1,:)

  !!----------------------
  !! Correlation settings
  !!----------------------

  !! Number of bands used to calculate intensity ratios
  Nbandr = 5

  ALLOCATE(indB(Nbandr))

  !! Find band index
  !!-----------------
  CALL IWHERE( TRIMEQ(labB(:),'Main 6.2 (1)'), indB(1) )
  CALL IWHERE( TRIMEQ(labB(:),'Main 7.7 (1)'), indB(2) )
  CALL IWHERE( TRIMEQ(labB(:),'Main 8.6'), indB(3) )
  CALL IWHERE( TRIMEQ(labB(:),'Main 11.2'), indB(4) )
  CALL IWHERE( TRIMEQ(labB(:),'Main 12.7 (1)'), indB(5) )
  
  !! Number of indpdt intensity ratios should be no larger
  Nrat = Nbandr - 1 ! fix I11.2

  ALLOCATE(ratname(Nrat), mulnR(Nrat), siglnR(Nrat))
  
  !! Denote ratio name
  !!-------------------
  ratname(1) = 'I6.2/I11.2'
  ratname(2) = 'I7.7/I11.2'
  ratname(3) = 'I8.6/I11.2'
  ratname(4) = 'I12.7/I11.2'

  mulnR(1) = LOG(.9_DP)
  mulnR(2) = LOG(1.8_DP)
  mulnR(3) = LOG(1.4_DP)
  mulnR(4) = LOG(.8_DP)

  siglnR(1) = LOG(1.2_DP) - mulnR(1)
  siglnR(2) = LOG(2.1_DP) - mulnR(2)
  siglnR(3) = LOG(1.5_DP) - mulnR(3)
  siglnR(4) = LOG(1._DP) - mulnR(4)

  !! Number of correlations (incl. non-correlated ones)
  Ncorr = N_CORR(Nrat)

  ALLOCATE(corrname(Ncorr), corr(Ncorr), lnR(Nrat,Np,Nq), &
           covar(Nrat,Nrat), Smat(Nrat,Nrat), Rmat(Nrat,Nrat))

  !! Denote correlation name (X - Y)
  !!---------------------------------
  corrname(1) = 'I6.2/I11.2 - I7.7/I11.2' ! r12
  corrname(2) = 'I6.2/I11.2 - I8.6/I11.2' ! r13
  corrname(3) = 'I6.2/I11.2 - I12.7/I11.2' ! r14
  corrname(4) = 'I7.7/I11.2 - I8.6/I11.2' ! r23
  corrname(5) = 'I7.7/I11.2 - I12.7/I11.2' ! r24
  corrname(6) = 'I8.6/I11.2 - I12.7/I11.2' ! r34

  !! Define correlation coefficients
  !!---------------------------------
  corr(1) = 0.95_DP ! r12
  corr(2) = 0.9_DP ! r13
  corr(3) = 0._DP ! r14
  corr(4) = 0.9_DP ! r23
  corr(5) = 0._DP ! r24
  corr(6) = 0._DP ! r34

  !! MC sampling of ratios
  !!-----------------------
  Smat(:,:) = 0._DP ! Std var matrix
  FORALL (i=1:Nrat) Smat(i,i) = siglnR(i)
  Rmat(:,:) = CORR2RMAT(corr(:), Nrat) ! Correlation matrix
  covar(:,:) = MATMUL( Smat(:,:),MATMUL(Rmat(:,:),Smat(:,:)) )

  !! Generate lnR(Nrat,Np,Nq)
  ! lnR(:,:,1) = RAND_MULTINORM(Np,COVAR=covar(:,:),MU=mulnR(:))
  DO q=1,Nq
    ! lnR(:,:,q) = lnR(:,:,1)
    lnR(:,:,q) = RAND_MULTINORM(Np,COVAR=covar(:,:),MU=mulnR(:))

  END DO

  !! Get band intensity ratios (lnRband in par)
  !!---------------------------------------
  
  !! Set ref band intensity (see also cont level)
  ! lnIref = pargen(1,1,ind%lnRband(ind%refB)) ! lnI11.2

  par(:,:,ind%lnRband(indB(1))) = lnR(1,:,:) ! lnI6.2 (1)
  ! par(:,:,ind%lnRband(indB(1)+1)) = LOG(.2_DP) ! lnI6.2 (2)
  par(:,:,ind%lnRband(indB(2))) = lnR(2,:,:) ! lnI7.7
  par(:,:,ind%lnRband(indB(3))) = lnR(3,:,:) ! lnI8.6
  par(:,:,ind%lnRband(indB(5))) = lnR(4,:,:) ! lnI12.7
  
  CALL WRITE_HDF5(STRARR1D=[spec_unit], NAME='Spectral unit', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
  CALL WRITE_HDF5(STRARR1D=labB(indB(:)), NAME='Correlated band name', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(INTARR1D=ind%lnRband(indB(:)), NAME='Correlated band indpar', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(STRARR1D=ratname(:), NAME='Band ratio name', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(STRARR1D=corrname(:), NAME='Correlation name', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  
  !!---------------------------------
  !! Generate spectra with the model
  !!--------------------------------
  FORALL (p=1:Np,q=1:Nq) &
    par(p,q,ind%lnMovd2(:Ncont)) = par(p,q,ind%lnMovd2(:Ncont)) + LOG(multiCONT(q))

  par(:,:,ind%lnAv(:)) = 0.5_DP
  
  ALLOCATE(FnuSIM(Np,Nq,Nw), dFnuSIM(Np,Nq,Nw))

  FnuSIM(:,:,:) = specModel(wvl(:), INDPAR=ind, PARVAL=par(:,:,:), &
                            FNUBAND_TAB=fnuband_tab, &
                            QABS=Qabs(:), EXTINCT=extinct(:,:))

  !! Statistical noises
  !!--------------------
  
  ! iw = MINLOC(dblarr1d(:),DIM=1) ! 1. min Fnu wvl index
  ! CALL IWHERE( dblarr1d(:)==MEDIAN_DATA(dblarr1d(:)), iw ) ! 2. median Fnu wvl index
  iw = CLOSEST(wOBS(:), wvl0) ! 3. ref wvl index
  print*, 'wvl (S/N) = ', wOBS(iw)
  
  !! Normalized S/N profile
  flatSovN = .FALSE.
  ! rescaling = 3._DP ! reduce uncertainties if > 1
  ALLOCATE(dblarr1d(NwOBS), normSovN(Nw))
  dblarr1d(:) = FnuOBS(1,1,:)/dFnuOBS(1,1,:)

  IF (flatSovN) THEN
    normSovN(:) = dblarr1d(iw) ! * rescaling
  ELSE
    normSovN(:) = INTERP_LIN_SORTED( dblarr1d(:)/dblarr1d(iw), &
                                     wOBS(:), wvl(:), XLOG=.FALSE., YLOG=.FALSE.)  
  END IF
  
  flatUNC = .FALSE.

  IF (flatUNC) THEN
    FORALL (q=1:Nq) dFnuSIM(:,q,:) = FnuSIM(1,q,iw) / SovN(q)
  ELSE
    FORALL (q=1:Nq,iw=1:Nw) dFnuSIM(:,q,iw) = FnuSIM(:,q,iw) / normSovN(iw) / SovN(q)

  END IF

  ALLOCATE(dblarr2d(Np,Nw))
  
  counter = 0
  CALL REALLOCATE(dblarr1d,Nw)
  ! dblarr1d(:) = RAND_NORM(Nw)
  spec_noising: DO
    q = counter + 1
    ! dblarr1d(:) = RAND_NORM(Nw)
    ! FORALL (p=1:Np) &
    !   dblarr2d(p,:) = dFnuSIM(p,q,:) * dblarr1d(:) + FnuSIM(p,q,:)
    DO p=1,Np
      dblarr2d(p,:) = dFnuSIM(p,q,:) * RAND_NORM(Nw) + FnuSIM(p,q,:)
      
    END DO
    
    ! IF (ALL(dblarr2d .GE. 0._DP)) THEN
      counter = counter + 1
      FnuSIM(:,q,:) = dblarr2d(:,:)

      IF (debug) THEN
        PRINT*
        PRINT*, 'spec_noising counter (q) = '//TRIMLR(PRING(counter))//' [done]'
      
      END IF
    ! END IF
    
    IF (counter .EQ. Nq) EXIT

  END DO spec_noising

  !! Write the output spectra
  !!--------------------------
  CALL WRITE_HDF5(DBLARR1D=wvl(:), NAME='Wavelength (microns)', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=FnuSIM(:,:,:), NAME='FnuOBS ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=dFnuSIM(:,:,:), NAME='dFnuOBS ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR4D=fnuband_tab(:,:,:,:), NAME='FnuBAND ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=par(:,:,:), NAME='Simulated parameter value', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

  PRINT*, 'Simulate spectra [done]'//NEW_LINE('')
  PRINT*,REPEAT("=",textwid)
  
END PROGRAM simulate_MIR

