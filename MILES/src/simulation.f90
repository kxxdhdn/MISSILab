!******************************************************************************
!*
!*         SIMULATING MID-INFRARED SPECTRA WITH NOISES AND CORRELATIONS
!*
!******************************************************************************

  
PROGRAM simulation

  USE utilities, ONLY: DP, trimeq, trimlr, pring
  USE arrays, ONLY: reallocate, closest, iwhere, incrarr
  USE statistics, ONLY: median_data, mean, sigma, N_corr, &
                        corr2Rmat, correlate
  USE random, ONLY: rand_norm, rand_multinorm
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
  CHARACTER(*), PARAMETER :: filGEN = dirIN//'presimu_hb'//h5ext
  CHARACTER(*), PARAMETER :: filOUT = dirIN//'simulation_MIR'//h5ext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .TRUE.
  
  INTEGER :: i, p, q, iw!, icorr, ipar
  INTEGER :: Np, Nq, Nw, NwOBS, Npar, Ncorr, Nrat, Nbandr, Ncont, Ncalib
  INTEGER :: counter
  INTEGER, DIMENSION(:), ALLOCATABLE :: indB, indBp ! 1D (Nbandr)
  INTEGER, DIMENSION(:), ALLOCATABLE :: indcal ! (NwCAL)
  REAL(DP) :: wvl0!, lnIref, rescaling
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dblarr1d
  REAL(DP), DIMENSION(:), ALLOCATABLE :: ln1pd ! 1D (Ncalib)
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
  LOGICAL :: varCONT, varSovN, flatSovN, flatUNC
  TYPE(indpar_type) :: ind
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE :: parinfo
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE :: Qabs ! 1D (Ncont)
  
  !!---------------
  !! Preliminaries
  !!---------------

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
    multiCONT(:) = [ .4_DP, 20._DP, &
                     .4_DP, 20._DP ]
  ELSE
    multiCONT(:) = 1._DP
  END IF

  !! Opt 2. S/N ratio
  !!------------------
  ALLOCATE(SovN(Nq))

  IF (varSovN) THEN
    SovN(:) = [ 2._DP, 2._DP, &
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
                   QABS=Qabs, EXTINCT=extinct, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit)

  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, &
                 NAME='FnuOBS ('//TRIMLR(spec_unit)//')')
  CALL READ_HDF5(DBLARR3D=dFnuOBS, FILE=filOBS, &
                 NAME='dFnuOBS ('//TRIMLR(spec_unit)//')')
  CALL READ_HDF5(STRARR1D=specOBS, FILE=filOBS, &
                 NAME='spectroscopic module labels')

  !! BB fitted par
  CALL READ_HDF5(DBLARR3D=pargen, FILE=filGEN, &
                 NAME='Mean of parameter value')

  !! Init par
  ALLOCATE(par(Np,Nq,Npar))

  FORALL (p=1:Np,q=1:Nq) par(p,q,:) = pargen(1,1,:)

  !!----------------------
  !! Correlation settings
  !!----------------------

  !! Number of bands used to calculate intensity ratios
  Nbandr = 6

  ALLOCATE(indB(Nbandr),indBp(Nbandr))

  !! Find band index
  !!-----------------
  CALL IWHERE( TRIMEQ(labB(:),'Main 3.3     '), indB(1) )
  CALL IWHERE( TRIMEQ(labB(:),'Main 6.2 (1) '), indB(2) )
  CALL IWHERE( TRIMEQ(labB(:),'Main 7.7 (1) '), indB(3) )
  CALL IWHERE( TRIMEQ(labB(:),'Main 8.6     '), indB(4) )
  CALL IWHERE( TRIMEQ(labB(:),'Main 12.7 (1)'), indB(5) )
  CALL IWHERE( TRIMEQ(labB(:),'Main 11.2    '), indB(6) )

  DO i=1,Nbandr
    CALL IWHERE( TRIMEQ(parinfo(:)%name,'lnRband'//TRIMLR(PRING(indB(i)))), indBp(i) )
  END DO
  
  !! Number of indpdt intensity ratios should be no larger
  Nrat = Nbandr - 1 ! fix I11.2

  ALLOCATE(ratname(Nrat), mulnR(Nrat), siglnR(Nrat))
  
  !! Denote ratio name
  !!-------------------
  ratname(1) = '3.3/11.2'
  ratname(2) = '6.2/11.2'
  ratname(3) = '7.7/11.2'
  ratname(4) = '8.6/11.2'
  ratname(5) = '12.7/11.2'

  mulnR(1) = pargen(1,1,indBp(1))
  mulnR(2) = pargen(1,1,indBp(2))
  mulnR(3) = pargen(1,1,indBp(3))
  mulnR(4) = pargen(1,1,indBp(4))
  mulnR(5) = pargen(1,1,indBp(5))
  IF (debug) THEN
    PRINT*
    PRINT*, 'mulnR (1-5) = ', mulnR
  END IF

  ! siglnR(:) = LOG(2._DP)
  siglnR(1) = LOG(2._DP)
  siglnR(2) = LOG(2._DP)
  siglnR(3) = LOG(2._DP)
  siglnR(4) = LOG(2._DP)
  siglnR(5) = LOG(2._DP)

  !! Number of correlations (incl. non-correlated ones)
  Ncorr = N_CORR(Nrat)

  ALLOCATE(corr(Ncorr), lnR(Nrat,Np,Nq), &
           covar(Nrat,Nrat), Smat(Nrat,Nrat), Rmat(Nrat,Nrat))

  !! Define correlation coefficients
  !!---------------------------------
  corr(1) = -.95_DP ! r12
  corr(2) = -.6_DP ! r13
  corr(3) = -.6_DP ! r14
  corr(4) = 0._DP ! r15
  corr(5) = .6_DP ! r23
  corr(6) = .6_DP ! r24
  corr(7) = 0._DP ! r25
  corr(8) = .95_DP ! r34
  corr(9) = 0._DP ! r35
  corr(10) = 0._DP ! r45

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
    ! print*, 'lnR', q, lnR(:,:,q)
    print*,

  END DO

  !! Get band intensity ratios (lnRband in par)
  !!---------------------------------------
  
  !! Set ref band intensity (see also cont level)
  ! lnIref = pargen(1,1,ind%lnRband(ind%refB)) ! lnI11.2

  par(:,:,ind%lnRband(indB(1))) = lnR(1,:,:) ! lnI3.3
  par(:,:,ind%lnRband(indB(2))) = lnR(2,:,:) ! lnI6.2 (1)
  ! par(:,:,ind%lnRband(indB(1)+1)) = LOG(.2_DP) ! lnI6.2 (2)
  par(:,:,ind%lnRband(indB(3))) = lnR(3,:,:) ! lnI7.7
  par(:,:,ind%lnRband(indB(4))) = lnR(4,:,:) ! lnI8.6
  par(:,:,ind%lnRband(indB(5))) = lnR(5,:,:) ! lnI12.7
  IF (debug) THEN
    print*, "1>>>>>", pargen(1,1,ind%lnRband(indB(1))) - par(:,:,ind%lnRband(indB(1)))
    print*, "2>>>>>", pargen(1,1,ind%lnRband(indB(2))) - par(:,:,ind%lnRband(indB(2)))
    print*, "3>>>>>", pargen(1,1,ind%lnRband(indB(3))) - par(:,:,ind%lnRband(indB(3)))
    print*, "4>>>>>", pargen(1,1,ind%lnRband(indB(4))) - par(:,:,ind%lnRband(indB(4)))
    print*, "5>>>>>", pargen(1,1,ind%lnRband(indB(5))) - par(:,:,ind%lnRband(indB(5)))
  END IF
  CALL WRITE_HDF5(STRARR1D=[spec_unit], NAME='spectral unit', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
  CALL WRITE_HDF5(STRARR1D=labB(indB(:)), NAME='Correlated band name', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(INTARR1D=ind%lnRband(indB(:)), NAME='Correlated band indpar', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(STRARR1D=ratname(:), NAME='Band ratio name', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  
  !!---------------------------------
  !! Generate spectra with the model
  !!---------------------------------

  FORALL (p=1:Np,q=1:Nq) &
    par(p,q,ind%lnFcont(:Ncont)) = par(p,q,ind%lnFcont(:Ncont)) + LOG(multiCONT(q))

  par(:,:,ind%lnAv(:)) = 0.5_DP
  
  ALLOCATE(FnuSIM(Np,Nq,Nw), dFnuSIM(Np,Nq,Nw), mask(Np,Nq,Nw))
  mask(:,:,:) = .TRUE.

  FnuSIM(:,:,:) = specModel(wvl(:), INDPAR=ind, PARVAL=par(:,:,:), &
                            MASK=mask(:,:,:), FNUBAND_TAB=fnuband_tab, &
                            QABS=Qabs(:), EXTINCT=extinct(:,:))

  !! Calib errors
  Ncalib = SIZE(specOBS)
  ALLOCATE (ln1pd(Ncalib))

  DO i=1,Ncalib
    CALL set_indcal(indcal, wOBS(:), INSTR=specOBS(i))
    ln1pd(i) = LOG(1 - (i-1) * .1_DP)
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
  iw = CLOSEST(wOBS(:), wvl0) ! 3. ref wvl index
  ! print*, 'wvl (S/N) = ', wOBS(iw)
  
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

  PRINT*, 'Simulate spectra [done]'//NEW_LINE('')
  PRINT*,REPEAT("=",textwid)
  
END PROGRAM simulation
