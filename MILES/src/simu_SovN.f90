!******************************************************************************
!*
!*         SIMULATING MID-INFRARED GALAXY SPECTRA (varying S/N ratio)
!*
!******************************************************************************

  
PROGRAM simu_SovN

  USE utilities, ONLY: DP, trimeq, trimlr, pring
  USE arrays, ONLY: reallocate, closest, iwhere, incrarr
  USE statistics, ONLY: median_data, mean, sigma, N_corr, &
                        corr2Rmat, correlate
  USE random, ONLY: generate_newseed, rand_norm, rand_multinorm
  USE interpolation, ONLY: interp_lin_sorted
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, lenpar
  USE grain_optics, ONLY: lendustQ
  USE core, ONLY: parinfo_type, indpar_type, Qabs_type, &
                  read_master, specModel, set_indcal
  USE ext_chi2, ONLY: mask
  USE ext_hb, ONLY: specOBS
  IMPLICIT NONE

  INTEGER, PARAMETER :: textwid = 60
  CHARACTER(*), PARAMETER :: dirIN = '../out/'
  CHARACTER(*), PARAMETER :: filOBS = dirIN//'observation_MIR'//h5ext
  CHARACTER(*), PARAMETER :: filGEN = dirIN//'presimu_bb'//h5ext
  CHARACTER(*), PARAMETER :: filOUT = dirIN//'simu_SovN'//h5ext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .TRUE.
  
  INTEGER :: i, p, q, iw, iw0!, icorr, ipar
  INTEGER :: Np, Nq, Nw, NwOBS, Npar, Ncorr, Nrat, NindB, Ncont, Ncalib
  INTEGER :: counter
  INTEGER, DIMENSION(:), ALLOCATABLE :: indB, indBp ! 1D (NindB)
  INTEGER, DIMENSION(:), ALLOCATABLE :: indcal ! (NwCAL)
  REAL(DP) :: dblval
  REAL(DP) :: wvl0!, lnIref, rescaling
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dblarr1d
  REAL(DP), DIMENSION(:), ALLOCATABLE :: pargen ! 1D (Npar)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: ln1pd ! 1D (Ncalib)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: muR, sigR ! 1D (Nrat)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: corr ! 1D (Ncorr)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wOBS ! 1D (NwOBS)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wvl, normSovN ! 1D (Nw)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: multiCONT, SovN ! 1D (Nq)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dblarr2d
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: extinct ! 2D (Nextc, Nw)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: covar, Smat, Rmat ! 2D (Nrat,Nrat)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: rat ! 3D (Nrat,Np,Nq)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: par ! 3D (Np,Nq,Npar)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuSIM, dFnuSIM ! 3D (Np,Nq,Nw) [U_in]
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuOBS, dFnuOBS ! 3D (1,1,NwOBS) [U_in]
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: fnuband_tab
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: qpar, qcal
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB, labS
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: nameR ! 1D (Nrat)
  LOGICAL :: varCONT, varSovN, flatSovN, flatUNC
  TYPE(indpar_type) :: ind
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE :: parinfo
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE :: Qabs ! 1D (Ncont)
  
  !!---------------
  !! Preliminaries
  !!---------------

  CALL GENERATE_NEWSEED()

  !! Banner
  PRINT*,REPEAT("=",textwid)
  PRINT*, REPEAT(" ",textwid/3)//'Synthetic Spectra of Galaxies'
  PRINT*,REPEAT("=",textwid)
  PRINT*

  !! Sample size
  Np = 10 ! sampling along correlation axis
  Nq = 4 ! sampling along option (FnuCONT or SovN) axis

  !! Options
  varCONT = .TRUE.
  varSovN = .TRUE.
  
  !! Wavelength reference
  wvl0 = 15._DP
  
  !! Opt 1. Continuum level
  !!------------------------
  ALLOCATE(multiCONT(Nq))
  
  IF (varCONT) THEN
    multiCONT(:) = [ 1._DP, 10._DP, &
                     1._DP, 10._DP ]
  ELSE
    multiCONT(:) = 1._DP
  END IF

  !! Opt 2. S/N ratio
  !!------------------
  ALLOCATE(SovN(Nq))

  IF (varSovN) THEN
    SovN(:) = [ 10._DP, 10._DP, &
                100._DP, 100._DP ]
  ELSE
    SovN(:) = 20._DP
  END IF

  !! Read the inputs
  !!--------------------
  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='wavelength (microns)', N1=NwOBS)
  !! Interpolate wavelength grid
  wvl = wOBS(:)
  Nw = SIZE(wvl)

  CALL READ_MASTER(WAVALL=wvl(:), DIRIN=dirIN, &
                   LABB=labB, &
                   NCONT=Ncont, &
                   QABS=Qabs, EXTINCT=extinct, LABS=labS, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit)

  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, &
                 NAME='FnuOBS ('//TRIMLR(spec_unit)//')')
  CALL READ_HDF5(DBLARR3D=dFnuOBS, FILE=filOBS, &
                 NAME='dFnuOBS ('//TRIMLR(spec_unit)//')')
  CALL READ_HDF5(STRARR1D=specOBS, FILE=filOBS, &
                 NAME='spectroscopic module labels')

  !! BB fitted par
  CALL READ_HDF5(DBLARR4D=qpar, FILE=filGEN, &
                 NAME='Quantiles of parameter value')
  CALL READ_HDF5(DBLARR4D=qcal, FILE=filGEN, &
                 NAME='Quantiles of ln(1+delta)')

  Ncalib = SIZE(specOBS)
  ALLOCATE(pargen(Npar),ln1pd(Ncalib))

  pargen(:) = qpar(1,1,:,2)
  ln1pd(:) = qcal(1,1,:,2)
  
  !! Init par
  ALLOCATE(par(Np,Nq,Npar))

  FORALL (p=1:Np,q=1:Nq) par(p,q,:) = pargen(:)

  !!----------------------
  !! Correlation settings
  !!----------------------

  !! Number of individual bands used to calculate intensity ratios of complexes
  NindB = 12

  ALLOCATE(indB(NindB),indBp(NindB))

  !! Find band index
  !!-----------------
  CALL IWHERE( TRIMEQ(labB(:),'Main 11.2    '), indB(1) ) ! 1
  CALL IWHERE( TRIMEQ(labB(:),'Plateau 11.3 '), indB(2) ) ! 1
  CALL IWHERE( TRIMEQ(labB(:),'Main 3.3     '), indB(3) ) ! 2
  CALL IWHERE( TRIMEQ(labB(:),'Main 3.4     '), indB(4) ) ! 3
  CALL IWHERE( TRIMEQ(labB(:),'Small 3.5    '), indB(5) ) ! 3
  CALL IWHERE( TRIMEQ(labB(:),'Main 6.2 (1) '), indB(6) ) ! 4
  CALL IWHERE( TRIMEQ(labB(:),'Main 6.2 (2) '), indB(7) ) ! 4
  CALL IWHERE( TRIMEQ(labB(:),'Plateau 7.7  '), indB(8) ) ! 5
  CALL IWHERE( TRIMEQ(labB(:),'Main 7.7 (1) '), indB(9) ) ! 5
  CALL IWHERE( TRIMEQ(labB(:),'Main 7.7 (2) '), indB(10) ) ! 5
  CALL IWHERE( TRIMEQ(labB(:),'Main 12.7 (1)'), indB(11) ) ! 6
  CALL IWHERE( TRIMEQ(labB(:),'Main 12.7 (2)'), indB(12) ) ! 6

  DO i=1,NindB
    CALL IWHERE( TRIMEQ(parinfo(:)%name,'lnRband'//TRIMLR(PRING(indB(i)))), indBp(i) )
  END DO
  
  !! Number of indpdt intensity ratios
  Nrat = 5

  ALLOCATE(nameR(Nrat), muR(Nrat), sigR(Nrat))
  
  !! Denote ratio name
  !!-------------------
  nameR(1) = '3.3/11.3'
  nameR(2) = '3.4/11.3'
  nameR(3) = '6.2/11.3'
  nameR(4) = '7.7/11.3'
  nameR(5) = '12.7/11.3'

  ! muR(:) = 0
  muR(1) = LOG( EXP(pargen(indBp(3))) &
                / (1+EXP(pargen(indBp(2)))) ) !- 1
  muR(2) = LOG( (EXP(pargen(indBp(4)))+EXP(pargen(indBp(5)))) &
                / (1+EXP(pargen(indBp(2)))) ) !- 1
  muR(3) = LOG( (EXP(pargen(indBp(6)))+EXP(pargen(indBp(7)))) &
                / (1+EXP(pargen(indBp(2)))) ) !- 1
  muR(4) = LOG( (EXP(pargen(indBp(8)))+EXP(pargen(indBp(9)))+ &
                 EXP(pargen(indBp(10)))) &
                / (1+EXP(pargen(indBp(2)))) ) !- 1
  muR(5) = LOG( (EXP(pargen(indBp(11)))+EXP(pargen(indBp(12)))) &
                / (1+EXP(pargen(indBp(2)))) ) !- 1
  IF (debug) THEN
    PRINT*
    PRINT*, 'muR (1-5) = ', muR
  END IF

  sigR(:) = 2._DP
  ! sigR(1) = 4._DP
  ! sigR(2) = 2._DP
  ! sigR(3) = 4._DP
  ! sigR(4) = 3._DP
  ! sigR(5) = 1._DP

  !! Number of correlations (incl. non-correlated ones)
  Ncorr = N_CORR(Nrat)

  ALLOCATE(corr(Ncorr), rat(Nrat,Np,Nq), &
           covar(Nrat,Nrat), Smat(Nrat,Nrat), Rmat(Nrat,Nrat))

  !! Define correlation coefficients
  !!---------------------------------
  corr(1) = -0.9_DP ! r12
  corr(2) =  0.75_DP ! r13
  corr(3) =  0.75_DP ! r14
  corr(4) =  0._DP ! r15
  corr(5) = -0.6_DP ! r23
  corr(6) = -0.6_DP ! r24
  corr(7) = -0._DP ! r25
  corr(8) =  0.9_DP ! r34
  corr(9) =  0._DP ! r35
  corr(10) = 0._DP ! r45
  
  !! MC sampling of ratios
  !!-----------------------
  Smat(:,:) = 0._DP ! Std var matrix
  FORALL (i=1:Nrat) Smat(i,i) = sigR(i)
  Rmat(:,:) = CORR2RMAT(corr(:), Nrat) ! Correlation matrix
  covar(:,:) = MATMUL( Smat(:,:),MATMUL(Rmat(:,:),Smat(:,:)) )

  !! Generate rat(Nrat,Np,Nq)
  ALLOCATE(dblarr2d(Nrat,Np))
  
  dblarr2d(:,:) = RAND_MULTINORM(Np,COVAR=covar(:,:),MU=muR(:))
  ! DO r=1,Nrat
  !   DO p=1,Np
  !     IF dblarr2d()
  !   END DO
  ! END DO
  
  DO q=1,Nq
    rat(:,:,q) = dblarr2d(:,:)
    ! PRINT*, 'r12', q, rat(1,:,q)/rat(2,:,q)
    ! PRINT*, 'r13', q, rat(1,:,q)/rat(3,:,q)
    ! PRINT*, 'r14', q, rat(1,:,q)/rat(4,:,q)
    ! PRINT*, 'r23', q, rat(2,:,q)/rat(3,:,q)
    ! PRINT*, 'r24', q, rat(2,:,q)/rat(4,:,q)
    ! PRINT*, 'r34', q, rat(3,:,q)/rat(4,:,q)
    ! PRINT*,

  END DO

  !! Get band intensity ratios (lnRband in par)
  !!---------------------------------------
  
  !! Set ref band intensity (see also cont level)
  ! lnIref = pargen(ind%lnRband(ind%refB)) ! lnI11.2
  
  par(:,:,ind%lnRband(indB(3)))  = &
       LOG( EXP(rat(1,:,:))*(1+EXP(pargen(indBp(2)))) ) ! lnR3.3
  CALL RANDOM_NUMBER( dblval )
  par(:,:,ind%lnRband(indB(4)))  = &
       LOG( EXP(rat(2,:,:))*(1+EXP(pargen(indBp(2))))*dblval ) ! lnR3.4
  par(:,:,ind%lnRband(indB(5)))  = &
       LOG( EXP(rat(2,:,:))*(1+EXP(pargen(indBp(2))))*(1-dblval) ) ! lnR3.5
  CALL RANDOM_NUMBER( dblval )
  par(:,:,ind%lnRband(indB(6)))  = &
       LOG( EXP(rat(3,:,:))*(1+EXP(pargen(indBp(2))))*dblval ) ! lnR6.2 (1)
  par(:,:,ind%lnRband(indB(7)))  = &
       LOG( EXP(rat(3,:,:))*(1+EXP(pargen(indBp(2))))*(1-dblval) ) ! lnR6.2 (2)
  CALL REALLOCATE( dblarr1d,2 )
  CALL RANDOM_NUMBER( dblarr1d(:) )
  par(:,:,ind%lnRband(indB(8)))  = &
       LOG( EXP(rat(4,:,:))*(1+EXP(pargen(indBp(2))))*dblarr1d(1) ) ! lnR7.7 plateau
  par(:,:,ind%lnRband(indB(9)))  = &
       LOG( EXP(rat(4,:,:))*(1+EXP(pargen(indBp(2))))*dblarr1d(2)*(1-dblarr1d(1)) ) ! lnR7.7 (1)
  par(:,:,ind%lnRband(indB(10))) = &
       LOG( EXP(rat(4,:,:))*(1+EXP(pargen(indBp(2)))) &
            *(1-dblarr1d(1)-dblarr1d(2)*(1-dblarr1d(1))) ) ! lnR7.7 (2)
  CALL RANDOM_NUMBER( dblval )
  par(:,:,ind%lnRband(indB(11))) = &
       LOG( EXP(rat(5,:,:))*(1+EXP(pargen(indBp(2))))*dblval ) ! lnR12.7 (1)
  par(:,:,ind%lnRband(indB(12))) = &
       LOG( EXP(rat(5,:,:))*(1+EXP(pargen(indBp(2))))*(1-dblval) ) ! lnR12.7 (2)
  ! IF (debug) THEN
  !   print*, "3.3 >>>>> ", pargen(ind%lnRband(indB(3))) - par(:,:,ind%lnRband(indB(3)))
  !   print*, "3.4 >>>>> ", pargen(ind%lnRband(indB(4))) - par(:,:,ind%lnRband(indB(4)))
  !   print*, "3.5 >>>>> ", pargen(ind%lnRband(indB(5))) - par(:,:,ind%lnRband(indB(5)))
  !   print*, "6.2 (1) >>>>> ", pargen(ind%lnRband(indB(6))) - par(:,:,ind%lnRband(indB(6)))
  !   print*, "6.2 (2) >>>>> ", pargen(ind%lnRband(indB(7))) - par(:,:,ind%lnRband(indB(7)))
  !   print*, "7.7 Plateau >>>>> ", pargen(ind%lnRband(indB(8))) - par(:,:,ind%lnRband(indB(8)))
  !   print*, "7.7 (1) >>>>> ", pargen(ind%lnRband(indB(9))) - par(:,:,ind%lnRband(indB(9)))
  !   print*, "7.7 (2) >>>>> ", pargen(ind%lnRband(indB(10))) - par(:,:,ind%lnRband(indB(10)))
  !   print*, "12.7 (1) >>>>> ", pargen(ind%lnRband(indB(11))) - par(:,:,ind%lnRband(indB(11)))
  !   print*, "12.7 (2) >>>>> ", pargen(ind%lnRband(indB(12))) - par(:,:,ind%lnRband(indB(12)))
  ! END IF
  CALL WRITE_HDF5(STRARR1D=[spec_unit], NAME='spectral unit', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
  CALL WRITE_HDF5(STRARR1D=labB(indB(:)), NAME='Correlated band name', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(INTARR1D=ind%lnRband(indB(:)), NAME='Correlated band indpar', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(STRARR1D=nameR(:), NAME='Band ratio name', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  
  !!---------------------------------
  !! Generate spectra with the model
  !!---------------------------------

  FORALL (p=1:Np,q=1:Nq) &
    par(p,q,ind%lnFcont(:Ncont)) = par(p,q,ind%lnFcont(:Ncont)) + LOG(multiCONT(q))

  ! par(:,:,ind%lnAv(:)) = 0.5_DP
  
  ALLOCATE(FnuSIM(Np,Nq,Nw), dFnuSIM(Np,Nq,Nw), mask(Np,Nq,Nw))
  mask(:,:,:) = .TRUE.

  FnuSIM(:,:,:) = specModel( wvl(:), INDPAR=ind, PARVAL=par(:,:,:), &
                             MASK=mask(:,:,:), FNUBAND_TAB=fnuband_tab, &
                             QABS=Qabs(:), EXTINCT=extinct(:,:), LABS=labS(:) )

  !! Apply calib factor
  DO i=1,Ncalib
    CALL set_indcal(indcal, wOBS(:), INSTR=specOBS(i))
    IF (debug) THEN
      PRINT*
      PRINT*, 'ln1pd ('//TRIMLR(PRING(i))//') = ', ln1pd(i)
    END IF
    DO iw=1,SIZE(indcal)
      FnuSIM(:,:,indcal(iw)) = FnuSIM(:,:,indcal(iw)) * EXP(ln1pd(i))
    END DO
    DEALLOCATE(indcal)
  END DO

  !! Statistical noises
  !!--------------------
  
  ! iw = MINLOC(dblarr1d(:),DIM=1) ! 1. min Fnu wvl index
  ! CALL IWHERE( dblarr1d(:)==MEDIAN_DATA(dblarr1d(:)), iw ) ! 2. median Fnu wvl index
  iw0 = CLOSEST(wOBS(:), wvl0) ! 3. ref wvl index
  ! print*, 'wvl (S/N) = ', wOBS(iw)
  
  !! Normalized S/N profile if flatSovN = .False.
  flatSovN = .False.
  
  ! rescaling = 3._DP ! reduce uncertainties if > 1
  ALLOCATE(normSovN(Nw))
  CALL REALLOCATE(dblarr1d,NwOBS)
  dblarr1d(:) = FnuOBS(1,1,:)/dFnuOBS(1,1,:)

  IF (flatSovN) THEN
    normSovN(:) = dblarr1d(iw0) ! * rescaling
  ELSE
    normSovN(:) = INTERP_LIN_SORTED( dblarr1d(:)/dblarr1d(iw0), &
                                     wOBS(:), wvl(:), XLOG=.FALSE., YLOG=.FALSE.)
    !! Clean AKARI band
    iw = CLOSEST(wvl(:),5._DP)
    normSovN(:iw) = normSovN(:iw)
  END IF
  
  flatUNC = .False.

  IF (flatUNC) THEN
    FORALL (q=1:Nq) dFnuSIM(:,q,:) = FnuSIM(1,q,iw0) / SovN(q)
  ELSE
    FORALL (q=1:Nq,iw=1:Nw) dFnuSIM(:,q,iw) = FnuSIM(:,q,iw) / normSovN(iw) / SovN(q)
    ! FORALL (q=1:Nq, p=1:Np) FnuSIM(p,q,:) = FnuSIM(p,q,:) * FnuOBS(p,q,iw0)/FnuSIM(p,q,iw0)
    ! FORALL (q=1:Nq, p=1:Np) dFnuSIM(p,q,:) = dFnuSIM(p,q,:) * FnuOBS(p,q,iw0)/FnuSIM(p,q,iw0)

  END IF

  CALL REALLOCATE(dblarr2d,Np,Nw)
  
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
  CALL WRITE_HDF5(DBLARR1D=wvl(:), NAME='wavelength (microns)', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=FnuSIM(:,:,:), NAME='FnuOBS ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=dFnuSIM(:,:,:), NAME='dFnuOBS ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR4D=fnuband_tab(:,:,:,:), NAME='FnuBAND ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=par(:,:,:), NAME='Simulated parameter value', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR1D=ln1pd(:), NAME='Simulated ln(1+delta)', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  
  PRINT*, 'Simulate spectra [done]'//NEW_LINE('')
  PRINT*,REPEAT("=",textwid)
  
END PROGRAM simu_SovN
