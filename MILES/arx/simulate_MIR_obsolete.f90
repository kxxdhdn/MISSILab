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
                        corr2Rmat
  USE random, ONLY: rand_norm, rand_multinorm
  USE interpolation, ONLY: interp_lin_sorted
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, lenpar
  USE grain_optics, ONLY: lendustQ
  IMPLICIT NONE

  INTEGER, PARAMETER :: textwid = 60
  CHARACTER(*), PARAMETER :: dirOUT = '../cache/'
  CHARACTER(*), PARAMETER :: filOBS = dirOUT//'galgen'//h5ext
  CHARACTER(*), PARAMETER :: filGEN = dirOUT//'pargen'//h5ext
  CHARACTER(*), PARAMETER :: filOUT = dirOUT//'galspec'//h5ext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .TRUE.
  
  INTEGER :: i, p, q, iw, icorr!, ipar
  INTEGER :: Np, Nq, Nw, NwOBS, Npar, Ncorr, Nratio, Nbandr
  INTEGER :: counter
  INTEGER, DIMENSION(:), ALLOCATABLE :: indB ! 1D (Nbandr)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: indRxy ! 2D (Ncorr,2)
  REAL(DP) :: wvl0, lnIband0, rescaling, sig
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dblarr1d
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dlnR ! 1D (Np)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: axy, bxy, sigxy, rhoxy ! 1D (Ncorr)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wOBS ! 1D (NwOBS)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wvl, extinct, normSovN ! 1D (Nw)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: theta ! 1D (Nratio)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: SovN0 ! 1D (Nq)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dblarr2d
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: lnRlim ! (Np,2) -> ratio limits LOG(min,max)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Rxy ! (Ncorr,2) -> ratio pair (X,Y)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: covar, matS, matR
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: pargen ! 3D (1,1,Npar)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: par ! 3D (Np,Nq,Npar)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: lnRxy ! 3D (Np,Ncorr,2) -> ratio pairs LOG(X,Y)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Iband ! 3D (Np,Nq,Nbandr) [(U_in * Hz)]
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: SovN ! 3D (Np,Nq,Nw)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuSIM, dFnuSIM ! 3D (Np,Nq,Nw) [U_in]
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuOBS, dFnuOBS ! 3D (1,1,NwOBS) [U_in]
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: ALLabQ
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: ALLabB, ALLabL
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB ! 1D (Nbandr)
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: ratname ! 1D (Nratio)
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: corrname ! 1D (Ncorr)
  LOGICAL :: accepted, varCONT, varSovN
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
  Np = 16 ! sampling along correlation axis
  Nq = 3 ! sampling along option (FnuCONT or SovN) axis

  !! Options
  varCONT = .FALSE.
  varSovN = .TRUE.
  
  !! Wavelength reference
  wvl0 = 15._DP

  !! Read the inputs
  !!--------------------
  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='Wavelength (microns)', N1=NwOBS)
  !! Interpolate wavelength grid
  wvl = wOBS(:)
  Nw = SIZE(wvl)

  CALL READ_MASTER(WAVALL=wvl(:), DIRIN=dirOUT, &
                   LABQ=ALLabQ, LABL=ALLabL, LABB=ALLabB, &
                   QABS=Qabs, EXTINCT=extinct(:), &
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
  CALL IWHERE( TRIMEQ(ALLabB(:),'Main 6.2 (1)'), indB(1) )
  CALL IWHERE( TRIMEQ(ALLabB(:),'Main 7.7 (1)'), indB(2) )
  CALL IWHERE( TRIMEQ(ALLabB(:),'Main 8.6'), indB(3) )
  CALL IWHERE( TRIMEQ(ALLabB(:),'Main 11.2'), indB(4) )
  CALL IWHERE( TRIMEQ(ALLabB(:),'Main 12.7 (1)'), indB(5) )
  
  !! Number of indpdt intensity ratios should be no larger
  Nratio = 5

  ALLOCATE(ratname(Nratio))
  
  !! Note ratio name
  !!-----------------
  ratname(1) = 'I6.2/I11.3'
  ratname(2) = 'I7.7/I11.3'
  ratname(3) = 'I8.6/I11.3'
  ratname(4) = 'I12.7/I7.7'
  ratname(5) = 'I6.2/I7.7' ! R(5) = R(1)/R(2)
  
  !! Number of correlations (incl. non-correlated ones)
  Ncorr = 3 ! 1 correlation corresponds to 1 ratio pair (X,Y)

  ALLOCATE(axy(Ncorr),bxy(Ncorr),sigxy(Ncorr), rhoxy(Ncorr), &
           lnRlim(Ncorr,2), dlnR(Ncorr), corrname(Ncorr), &
           indRxy(Ncorr,2), Rxy(Ncorr,2), lnRxy(Np,Ncorr,2))

  !! Note correlation name
  !!-----------------------
  corrname(1) = 'I6.2/I11.3 - I7.7/I11.3' ! R(1)/R(2)
  corrname(2) = 'I8.6/I11.3 - I7.7/I11.3' ! R(3)/R(2)
  corrname(3) = 'I12.7/I7.7 - I6.2/I7.7' ! R(4)/R(5)

  !! Link correlations with band ratios
  !!------------------------------------
  indRxy(1,:) = [2,1]
  indRxy(2,:) = [2,3]
  indRxy(3,:) = [5,4]

  !! Define correlations
  !!---------------------
  axy(:) = 0._DP
  bxy(:) = 0._DP
  sigxy(:) = 0._DP ! dispersion of (Y - axy*X - bxy)
  rhoxy(:) = 0._DP ! corr coeff (non-correlated)
  !! I6.2/I11.3 - I7.7/I11.3 (M83, Lorentz fit, Galliano08)
  axy(1) = 2.19_DP
  bxy(1) = 0._DP
  sigxy(1) = 0.54_DP
  rhoxy(1) = 0.67_DP
  !! I8.6/I11.3 - I7.7/I11.3 (M83, Lorentz fit, Galliano08)
  axy(2) = 3.52_DP
  bxy(2) = 0._DP
  sigxy(2) = 0.70_DP
  rhoxy(2) = 0.66_DP
  !! I12.7/I7.7 - I6.2/I7.7 (supposed to be not correlated)

  !! Set ratio ranges
  !!------------------
  lnRlim(:,1) = LOG(0.1_DP) ! lnRmin
  lnRlim(:,2) = LOG(10._DP) ! lnRmax
  dlnR(:) = lnRlim(:,2) - lnRlim(:,1)
  
  !! MC sampling of ratios
  !!-----------------------
  ALLOCATE(theta(Nratio))
  
  counter = 0
  MC_sampling: DO
    CALL RANDOM_NUMBER(theta(:))
    accepted = .TRUE.

    !! Satisfy all correlations and then accept
    corr_checking: DO icorr=1,Ncorr
      !! Randomly generated ratios X & Y (in LOG)
      Rxy(icorr,:) = EXP(theta(indRxy(icorr,:)) * dlnR(icorr) + lnRlim(icorr,1))
      
      IF (rhoxy(icorr) .NE. 0._DP) THEN
        sig = ABS( Rxy(icorr,2) - Rxy(icorr,1) * axy(icorr) - bxy(icorr) )
        !! Accept only if all (X,Y) pairs accepted
        IF ( .NOT. (sig .LE. sigxy(icorr) .AND. accepted) ) THEN
          accepted = .FALSE.
          EXIT corr_checking

        END IF
      END IF
      
      IF (debug .AND. accepted) THEN
        PRINT*, 'Correlation No. '//TRIMLR(PRING(icorr))//' satisfied. '
        
      END IF
    END DO corr_checking

    !! Save accepted sampling elements
    IF (accepted) THEN
      counter = counter + 1
      lnRxy(counter,:,:) = LOG(Rxy(:,:))
    END IF
    
    IF (debug) THEN
      PRINT*, 'MC_sampling counter (p) = '//TRIMLR(PRING(counter))
      PRINT*
    END IF
    
    IF (counter .EQ. Np) EXIT
    
  END DO MC_sampling

  !! Get band intensities (lnIband in par)
  !!---------------------------------------
  
  !! Set Iband benchmark (see also cont level)
  lnIband0 = pargen(1,1,ind%lnIband(indB(4))) ! lnI11.3

  !! Modify and extract Iband
  DO p=1,Np
    par(p,:,ind%lnIband(indB(1))) = lnRxy(p,1,2) + lnIband0 ! lnI6.2
    par(p,:,ind%lnIband(indB(2))) = lnRxy(p,1,1) + lnIband0 ! lnI7.7
    par(p,:,ind%lnIband(indB(3))) = lnRxy(p,2,2) + lnIband0 ! lnI8.6
    par(p,:,ind%lnIband(indB(5))) = lnRxy(p,3,2) + lnRxy(p,1,1) + lnIband0 ! lnI12.7

  END DO

  ALLOCATE(labB(Nbandr), Iband(Np,Nq,Nbandr))
  
  DO i=1,Nbandr
    labB(i) = ALLabB(indB(i))
    !! Intensity unit in W/m2/sr
    SELECT CASE(spec_unit)
      CASE ('MKS')
        Iband(:,:,i) = EXP(par(:,:,ind%lnIband(indB(i))))
      CASE ('MJyovsr')
        Iband(:,:,i) = EXP(par(:,:,ind%lnIband(indB(i)))) * 1.E-20_DP
      CASE DEFAULT
        Iband(:,:,i) = EXP(par(:,:,ind%lnIband(indB(i))))

    END SELECT
  END DO

  CALL WRITE_HDF5(STRARR1D=[spec_unit], NAME='Spectral unit', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
  CALL WRITE_HDF5(STRARR1D=labB(:), NAME='Correlated band name', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(INTARR1D=ind%lnIband(indB(:)), NAME='Correlated band indpar', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=Iband(:,:,:), NAME='Correlated band intensity (Wovm2ovsr)', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(STRARR1D=corrname(:), NAME='Correlation name', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
print*, ind%lnIband(indB(:))
  !!------------------------
  !! Opt 1. Continuum level
  !!------------------------

  !!------------------
  !! Opt 2. S/N ratio
  !!------------------
  ALLOCATE(SovN(Np,Nq,Nw), SovN0(Nq), normSovN(Nw), dblarr1d(NwOBS))

  !! S/N profile
  !!-------------
  dblarr1d(:) = FnuOBS(1,1,:)/dFnuOBS(1,1,:)
  iw = MINLOC(dblarr1d(:),DIM=1)
  
  !! SovN normalized at its minimum value
  normSovN(:) = INTERP_LIN_SORTED( dblarr1d(:)/dblarr1d(iw), &
                                   wOBS(:), wvl(:), XLOG=.FALSE., YLOG=.FALSE.)  

  IF (varSovN) THEN
    SovN0(:) = [0.3_DP, 1._DP,3._DP]
  ELSE
    SovN0(:) = 100._DP

  END IF
  
  FORALL (q=1:Nq,iw=1:Nw) &
    SovN(:,q,iw) = normSovN(iw) * SovN0(q)

  !!---------------------------------
  !! Generate spectra with the model
  !!---------------------------------
  ALLOCATE(FnuSIM(Np,Nq,Nw), dFnuSIM(Np,Nq,Nw))

  FnuSIM(:,:,:) = specModel(wvl(:), INDPAR=ind, PARVAL=par(:,:,:), &
                            QABS=Qabs(:), EXTINCT=extinct(:))

  !! Statistical noises
  !!--------------------
  rescaling = 3._DP
  !! Vary FnuSIM by fixing dFnuSIM for the whole sample
  FORALL (iw=1:Nw) dFnuSIM(:,:,iw) = FnuSIM(:,:,iw) / normSovN(iw) / rescaling

  ALLOCATE(dblarr2d(Np,Nw))
  
  counter = 0
  CALL REALLOCATE(dblarr1d, Nw)
  spec_noising: DO
    q = counter + 1
    dblarr1d(:) = RAND_NORM(Nw)
    FORALL (p=1:Np) &
      dblarr2d(p,:) = dFnuSIM(p,q,:) * dblarr1d(:) &
                      + FnuSIM(p,q,:) / rescaling * SovN0(q)
    
    IF (ALL(dblarr2d .GE. 0._DP)) THEN
      counter = counter + 1
      FnuSIM(:,q,:) = dblarr2d(:,:)

      IF (debug) THEN
        PRINT*
        PRINT*, 'spec_noising counter (q) = '//TRIMLR(PRING(counter))//' [done]'
      
      END IF
    END IF
    
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

  PRINT*, 'Simulate spectra [done]'//NEW_LINE('')
  PRINT*,REPEAT("=",textwid)
  
END PROGRAM simulate_MIR

