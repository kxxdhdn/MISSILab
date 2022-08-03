!******************************************************************************
!*
!*                             MILES CORE FUNCTIONS
!*
!******************************************************************************


MODULE core

  USE auxil, ONLY: Ncont_max, Nline_max, Nband_max, Nextc_max, Nstar_max, Cband_sig
  USE utilities, ONLY: DP, CDP
  ! USE constants, ONLY: 
  USE inout, ONLY: lenpar
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: NparCONT=2, NparLINE=3, NparBAND=4, NparEXTC=1, NparSTAR=1
  INTEGER, PARAMETER, PUBLIC :: lencorr = 41
  
  PUBLIC :: set_indpar, set_indref, set_indcal
  PUBLIC :: make_Qabs, read_master, initparam, read_analysis
  PUBLIC :: check_SM, invert_SM, invert_mSM, NaN2zero, zero2NaN
  PUBLIC :: degradeRes, modifBB, gaussLine, lorentzBand, extCurve, specModel

  INTERFACE NaN2zero
    MODULE PROCEDURE NaN2zero_int1D, NaN2zero_int2D, NaN2zero_int3D, NaN2zero_int4D
    MODULE PROCEDURE NaN2zero_dbl1D, NaN2zero_dbl2D, NaN2zero_dbl3D, NaN2zero_dbl4D
  END INTERFACE NaN2zero
  
  INTERFACE zero2NaN
    MODULE PROCEDURE zero2NaN_int1D, zero2NaN_int2D, zero2NaN_int3D, zero2NaN_int4D
    MODULE PROCEDURE zero2NaN_dbl1D, zero2NaN_dbl2D, zero2NaN_dbl3D, zero2NaN_dbl4D
  END INTERFACE zero2NaN
  
  INTERFACE modifBB
    MODULE PROCEDURE modifBB_0, modifBB_1
  END INTERFACE modifBB

  INTERFACE gaussLine
    MODULE PROCEDURE gaussLine_0, gaussLine_1
  END INTERFACE gaussLine

  INTERFACE lorentzBand
    MODULE PROCEDURE lorentzBand_0, lorentzBand_1
  END INTERFACE lorentzBand
  
  INTERFACE specModel
    MODULE PROCEDURE specModel_3D, specModel_2D, specModel_1D
    MODULE PROCEDURE specModel_gen, specModel_scl
  END INTERFACE specModel
 
  TYPE, PUBLIC :: parinfo_type
    CHARACTER(lenpar) :: name = ''
    CHARACTER(lenpar) :: comp = ''
    CHARACTER(lenpar) :: tied = ''
    LOGICAL :: model = .TRUE.
    LOGICAL :: hyper = .FALSE.
    LOGICAL :: fixed = .TRUE.
    LOGICAL, DIMENSION(2) :: limited = [.FALSE., .FALSE.]
    REAL(DP), DIMENSION(2) :: limits = [0._DP, 0._DP]
    REAL(DP) :: value = 0._DP
    REAL(DP) :: mean = 0._DP
    REAL(DP) :: sigma = 0._DP
    INTEGER :: ind = -1
  END TYPE parinfo_type

  TYPE, PUBLIC :: indpar_type
    !! The actual unit depends on the input spectrum.
    !! Here as an example, suppose the input is in MKS.
    ! INTEGER, DIMENSION(Ncont_max) :: lnMovd2 = -1 ! dust MBB coeff LOG[Msun/pc2]
    INTEGER :: refw = 1 ! index of ref wavelength to define lnFcont
    INTEGER, DIMENSION(Ncont_max) :: lnFcont = -1 ! dust MBB flux at ref wavelength LOG[W/m2/sr]
    !! If there are identical cont composants, the second lnT becomes ln(deltaT), etc.
    INTEGER, DIMENSION(Ncont_max) :: lnT = -1 ! dust temperature (MBB) LOG[K]
    INTEGER, DIMENSION(Ncont_max) :: grpQ = 0 ! dust compo group nb (1,2,3,...), for unsorted labQ
    INTEGER, DIMENSION(Ncont_max) :: ordQ = 0 ! dust compo group order (0,1,2,...)
    INTEGER, DIMENSION(Nline_max) :: lnRline = -1 ! line intensity over ref band Ratio LOG[-]
    INTEGER, DIMENSION(Nline_max) :: Cline = -1 ! Center [um]
    INTEGER, DIMENSION(Nline_max) :: Wline = -1 ! Width [um]
    INTEGER :: refB = 1 ! index of ref band in labB to define lnRband and lnRline
    !! lnRband of ref band is LOG(intensity [W/m2/sr])
    INTEGER, DIMENSION(Nband_max) :: lnRband = -1 ! band intensity over ref band Ratio LOG[-]
    INTEGER, DIMENSION(Nband_max) :: Cband = -1 ! Center (peak) [um]
    INTEGER, DIMENSION(Nband_max) :: WSband = -1 ! Width Short nu side [um]
    INTEGER, DIMENSION(Nband_max) :: WLband = -1 ! Width Long nu side [um]
    INTEGER, DIMENSION(Nextc_max) :: lnAv = -1 ! Extinction LOG[mag]
    !! Total star surface brightness with dilution factor Omega = r/d LOG[W/m2/sr]
    INTEGER, DIMENSION(Nstar_max) :: lnFstar = -1
    INTEGER, DIMENSION(:), ALLOCATABLE :: extra
  END TYPE indpar_type

  TYPE, PUBLIC :: irrintarr_type
    !! Used to create irregular integer arrays
    INTEGER, DIMENSION(:), ALLOCATABLE :: arr1d
  END TYPE irrintarr_type

  TYPE, PUBLIC :: irrdblarr_type
    !! Used to create irregular real arrays
    REAL(DP), DIMENSION(:), ALLOCATABLE :: arr1d
  END TYPE irrdblarr_type
  
  TYPE, PUBLIC :: irrcmparr_type
    !! Used to create irregular complex arrays
    COMPLEX(CDP), DIMENSION(:), ALLOCATABLE :: arr1d
  END TYPE irrcmparr_type

  TYPE, PUBLIC :: irrlogarr_type
    !! Used to create irregular logical arrays
    LOGICAL, DIMENSION(:), ALLOCATABLE :: arr1d
  END TYPE irrlogarr_type

  TYPE, PUBLIC :: irrchararr_type
    !! Used to create irregular character arrays
    CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: arr1d
  END TYPE irrchararr_type
  
  TYPE, PUBLIC :: Qabs_type
    !! rho(kg/m3); wave(micron); nu(Hz); kappa(m2/kg)
    REAL(DP) :: rho
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave, nu
    REAL(DP), DIMENSION(:), ALLOCATABLE :: kappa
  END TYPE Qabs_type
 
CONTAINS

  !!-------------------------------------------------------
  !!
  !!                Read optical properties
  !!
  !!-------------------------------------------------------
  SUBROUTINE make_Qabs( label, Qabs, wavALL )
    !! Qabs % rho(kg/m3); wave(micron); nu(Hz); kappa(m2/kg)
    USE utilities, ONLY: DP
    USE constants, ONLY: MKS
    USE arrays, ONLY: closest
    USE interpolation, ONLY: interp_lin_sorted
    USE grain_optics, ONLY: rho_grain, read_optics
    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:), INTENT(IN) :: label
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: wavALL
    INTEGER :: i, Ngrain, Nr, Nw
    INTEGER, DIMENSION(SIZE(label)) :: indrad ! radius index
    REAL(DP), PARAMETER :: a0 = 1.E-2_DP ! grain radius
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave0 ! READ_OPTICS default grid
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: radius
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Qabsall ! (Ngrain, Nr, Nw)
    TYPE(Qabs_type), DIMENSION(SIZE(label)), INTENT(OUT) :: Qabs

    CALL READ_OPTICS(label, WAVE=wave0, RADIUSALL=radius, QABSALL=Qabsall)    

    wave = wave0(:)
    IF (PRESENT(wavAll)) wave = wavAll(:)
    Nw = SIZE(wave(:))

    Ngrain = SIZE(label)
    DO i=1,Ngrain
      ALLOCATE(Qabs(i)%wave(Nw),Qabs(i)%nu(Nw),Qabs(i)%kappa(Nw))
      
      Nr = SIZE(radius(i,:))
      Qabs(i)%rho = rho_grain(label(i))
      Qabs(i)%wave(:) = wave(:)
      Qabs(i)%nu(:) = MKS%clight/MKS%micron / wave(:)
      indrad(i) = CLOSEST(radius(i,:), a0)
      ! PRINT*, 'Radius of ', label(i), ': ', radius(i,indrad(i)), '->', indrad(i)
      IF (PRESENT(wavAll)) THEN
        Qabs(i)%kappa(:) = 3._DP/4._DP/rho_grain(label(i))/radius(i,indrad(i)) * &
                           INTERP_LIN_SORTED( Qabsall(i,indrad(i),:), &
                             wave0(:),wave(:),XLOG=.TRUE.,YLOG=.TRUE.,FORCE=.TRUE. )
      ELSE
        Qabs(i)%kappa(:) = 3._DP/4._DP/rho_grain(label(i))/radius(i,indrad(i)) * &
                           Qabsall(i,indrad(i),:)
    
      END IF
    END DO
    
    !! Free memory space
    DEALLOCATE(wave0, radius, Qabsall)
  
  END SUBROUTINE make_Qabs

  !!-------------------------------------------------------
  !!
  !! Fill the INDPAR_TYPE structure, from a PARINFO_TYPE structure
  !!
  !!-------------------------------------------------------
  SUBROUTINE set_indpar( indpar, parinfo, refB, labB, refw, labQ )
    
    USE utilities, ONLY: trimeq, trimlr, pring, strike
    USE arrays, ONLY: iwhere
    IMPLICIT NONE

    TYPE(indpar_type), INTENT(OUT) :: indpar
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo
    INTEGER, INTENT(IN), OPTIONAL :: refw
    CHARACTER(*), INTENT(IN), OPTIONAL :: refB
    CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: labB, labQ
    
    INTEGER :: i, Ngrp
    INTEGER :: Ncont, Nline, Nband, Nextc, Nstar, Nextra
    INTEGER, DIMENSION(:), ALLOCATABLE :: igrp
    LOGICAL :: CONTset, LINEset, BANDset, EXTCset, STARset, extraset

    CONTset = ( ANY(TRIMEQ(parinfo(:)%comp, 'CONT')) )
    LINEset = ( ANY(TRIMEQ(parinfo(:)%comp, 'LINE')) )
    BANDset = ( ANY(TRIMEQ(parinfo(:)%comp, 'BAND')) )
    EXTCset = ( ANY(TRIMEQ(parinfo(:)%comp, 'EXTC')) )
    STARset = ( ANY(TRIMEQ(parinfo(:)%comp, 'STAR')) )
    extraset = ( ANY(TRIMEQ(parinfo(:)%comp, 'extra')) )

    IF (CONTset) THEN
      Ncont = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'CONT')) ) / NparCONT
      Ngrp = 0
      DO i=1,Ncont
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'lnFcont'//TRIMLR(PRING(i))), indpar%lnFcont(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'lnT'//TRIMLR(PRING(i))), indpar%lnT(i) )

        !! Define dust compo groups indpar%grpQ
        IF (PRESENT(labQ)) THEN
          CALL IWHERE( TRIMEQ(labQ(1:i), labQ(i)), igrp )
          IF (SIZE(igrp)==1) THEN
            Ngrp = Ngrp + 1
            indpar%grpQ(i) = Ngrp
          ELSE
            indpar%grpQ(i) = indpar%grpQ(igrp(1))
            indpar%ordQ(i) = SIZE(igrp) - 1
          END IF
        END IF
        
      END DO
    END IF
    
    IF (LINEset) THEN
      Nline = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'LINE')) ) / NparLINE
      DO i=1,Nline
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'lnRline'//TRIMLR(PRING(i))), indpar%lnRline(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'Cline'//TRIMLR(PRING(i))), indpar%Cline(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'Wline'//TRIMLR(PRING(i))), indpar%Wline(i) )

      END DO
    END IF

    IF (BANDset) THEN
      Nband = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'BAND')) ) / NparBAND
      DO i=1,Nband
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'lnRband'//TRIMLR(PRING(i))), indpar%lnRband(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'Cband'//TRIMLR(PRING(i))), indpar%Cband(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'WSband'//TRIMLR(PRING(i))), indpar%WSband(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'WLband'//TRIMLR(PRING(i))), indpar%WLband(i) )
        
      END DO
    END IF

    IF (EXTCset) THEN
      Nextc = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'EXTC')) ) / NparEXTC
      DO i=1,Nextc
        CALL IWHERE( TRIMEQ(parinfo%name(:),'lnAv'//TRIMLR(PRING(i))), indpar%lnAv(i) )
        
      END DO
    END IF
    
    IF (STARset) THEN
      Nstar = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'STAR')) ) / NparSTAR
      DO i=1,Nstar
        CALL IWHERE( TRIMEQ(parinfo%name(:),'lnFstar'//TRIMLR(PRING(i))), indpar%lnFstar(i) )
        
      END DO
    END IF

    IF (extraset) THEN
      Nextra = COUNT(TRIMEQ(parinfo(:)%comp, 'extra'))
      IF (Nextra > 0) CALL IWHERE( TRIMEQ(parinfo(:)%comp,'extra'),indpar%extra )
      
    END IF
    
    IF (PRESENT(refB)) THEN
      IF (PRESENT(labB)) THEN
        CALL IWHERE( TRIMEQ(labB, refB), indpar%refB )
      ELSE
        CALL STRIKE('SET_INDPAR','Input band list (labB)')
      END IF
    END IF

    IF (PRESENT(refw)) indpar%refw = refw

  END SUBROUTINE set_indpar

  !!-------------------------------------------------------
  !!
  !! Find index of the reference band intensity to calculated ratios
  !!
  !!-------------------------------------------------------
  SUBROUTINE set_indref( indref, namref, labB )

    USE utilities, ONLY: trimeq
    USE arrays, ONLY: iwhere

    INTEGER, INTENT(OUT) :: indref
    CHARACTER(*), INTENT(IN) :: namref
    CHARACTER(*), DIMENSION(:), INTENT(IN) :: labB

    CALL IWHERE( TRIMEQ(labB, namref), indref ) ! indref = index of ref in labB

  END SUBROUTINE set_indref

  !!-------------------------------------------------------
  !!
  !!             Divide wavelengths by modules
  !!
  !!-------------------------------------------------------
  SUBROUTINE set_indcal ( indcal, x, instr, calibool, caliberr )
    !! 
    USE auxil, ONLY: ran, err
    USE utilities, ONLY: DP
    USE arrays, ONLY: iwhere

    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: indcal
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    CHARACTER(*), INTENT(IN) :: instr
    REAL(DP) :: caliberr0
    REAL(DP), DIMENSION(2) :: lim
    REAL(DP), INTENT(OUT), OPTIONAL :: caliberr
    LOGICAL, DIMENSION(SIZE(x)), INTENT(OUT), OPTIONAL :: calibool

    caliberr0 = 0._DP
    SELECT CASE(instr)
      CASE('IRC_NG')
        lim(:) = ran%wlim_AKARI_NG
        caliberr0 = err%caliberr_AKARI_NG
      CASE('IRS_SL2')
        lim(:) = ran%wlim_SL2
        caliberr0 = err%caliberr_SL
      CASE('IRS_SL1')
        lim(:) = ran%wlim_SL1
        caliberr0 = err%caliberr_SL
      CASE('IRS_LL2')
        lim(:) = ran%wlim_LL2
        caliberr0 = err%caliberr_LL
      CASE('IRS_LL1')
        lim(:) = ran%wlim_LL1
        caliberr0 = err%caliberr_LL
      CASE('IRS_SH')
        lim(:) = ran%wlim_SH
        caliberr0 = err%caliberr_SH
      CASE('IRS_LH')
        lim(:) = ran%wlim_LH
        caliberr0 = err%caliberr_LH
      CASE DEFAULT
        lim(:) = [2.50, 40.00]
        caliberr0 = 0._DP
      END SELECT

    IF (PRESENT(caliberr)) caliberr = caliberr0

    IF (PRESENT(calibool)) calibool(:) = .FALSE.
    IF (ANY(x(:)>lim(1) .AND. x(:)<lim(2))) THEN
      CALL IWHERE( x(:)>lim(1) .AND. x(:)<lim(2), indcal )
      IF (PRESENT(calibool)) calibool(indcal(:)) = .TRUE.
    ELSE
      ALLOCATE(indcal(0))
    END IF

  END SUBROUTINE set_indcal
  
  !!-------------------------------------------------------
  !!
  !! Read the input master files for the Chi2/Bayesian run
  !!
  !!-------------------------------------------------------
  SUBROUTINE read_master( wavAll, dirIN, dirOUT, spec_unit, &
                          Nmcmc, NiniMC, verbose, &
                          calib, robust_cal, robust_RMS, skew_RMS, &
                          resume, indresume, newseed, newinit, dostop, nohi, &
                          labL, labB, refB, labQ, Qabs, refw, labE, extinct, &
                          Ncont, Nband, Nline, Nextc, Nstar, Nextra, &
                          corrhypname, corrname, &
                          parinfo, parmodinfo, parhypinfo, parextinfo, &
                          indpar, Npar, Nparmod, Nparhyp, Ncorrhyp, Ncorr)

    USE utilities, ONLY: DP, strike, warning, trimeq, trimlr, pring!, verbatim
    USE inout, ONLY: lenpar, lenpath, read_input_line, read_hdf5, h5ext
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: reallocate, closest
    USE interpolation, ONLY: interp_lin_sorted
    USE statistics, ONLY: correl_parlist, N_corr
    USE grain_optics, ONLY: lendustQ, rho_grain, read_optics
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN), OPTIONAL :: dirIN ! default: ./
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: wavAll ! Interp Qabs if presents

    CHARACTER(lenpath) :: filmas, filmod, filext
    ! CHARACTER(lenpath) :: dirOUT0
    ! CHARACTER(lenpar) :: spec_unit0
    CHARACTER(lendustQ) :: refB0
    CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ0, labL0, labB0, labE0
    ! INTEGER :: Nmcmc0, NiniMC0
    INTEGER :: Ncont0, Nband0, Nline0, Nextc0, Nstar0, Nextra0, Npar0
    INTEGER :: Nparmod0, Nparhyp0, Ncorrhyp0, Ncorr0, Nparall
    ! INTEGER :: i, iostat, ipar
    REAL(DP) :: refw0
    TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE :: parinfall, parinfo0
    ! TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE :: parinfoextra
    ! LOGICAL :: verbose0!, robust_RMS0, robust_cal0, skew_RMS0
    ! LOGICAL :: calib0, newseed0, newinit0, dostop0
    LOGICAL, DIMENSION(:), ALLOCATABLE :: boolparmod, boolpar

    CHARACTER(lenpath), DIMENSION(:), ALLOCATABLE :: path1d
    CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: parhypname
    CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: strarr1d
    CHARACTER(lenpar), DIMENSION(:,:), ALLOCATABLE :: strarr2d
    INTEGER, DIMENSION(:), ALLOCATABLE :: intarr1d
    REAL(DP), DIMENSION(:), ALLOCATABLE :: dblarr1d
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dblarr2d
    INTEGER :: i, Nr, Nw
    INTEGER, DIMENSION(:), ALLOCATABLE :: indrad ! radius index
    REAL(DP), PARAMETER :: a0 = 1.E-2_DP ! grain radius
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave0 ! READ_OPTICS default grid
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: radius ! (Nr, Nw)
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Qabsall ! (Ncont, Nr, Nw)
    
    CHARACTER(lenpath), INTENT(OUT), OPTIONAL :: dirOUT
    CHARACTER(lenpar), INTENT(OUT), OPTIONAL :: spec_unit
    CHARACTER(lendustQ), INTENT(OUT), OPTIONAL :: refB
    CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      labQ, labB, labL, labE
    CHARACTER(lencorr), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      corrhypname, corrname
    INTEGER, INTENT(OUT), OPTIONAL :: Nmcmc, NiniMC, indresume, &
      Ncont, Nband, Nline, Nextc, Nstar, Nextra
    INTEGER, INTENT(OUT), OPTIONAL :: Npar, Nparmod, Nparhyp, Ncorrhyp, Ncorr
    REAL(DP), INTENT(OUT), OPTIONAL :: refw
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: extinct
    TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Qabs
    TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      parinfo, parmodinfo, parhypinfo, parextinfo
    TYPE(indpar_type), INTENT(OUT), OPTIONAL :: indpar
    LOGICAL, INTENT(OUT), OPTIONAL :: verbose, robust_RMS, robust_cal, skew_RMS
    LOGICAL, INTENT(OUT), OPTIONAL :: calib, newseed, newinit, dostop, resume, nohi

    !! Read the input master file
    !!----------------------------
    IF (PRESENT(dirIN)) THEN
      filmas = TRIMLR(dirIN)//'input_master'//h5ext
      filmod = TRIMLR(dirIN)//'input_model'//h5ext
      filext = TRIMLR(dirIN)//'input_extra'//h5ext
    ELSE
      CALL WARNING('READ_MASTER','Default input path')
      filmas = './input_master'//h5ext
      filmod = './input_model'//h5ext
      filext = './input_extra'//h5ext
    END IF

    !! Output path
    IF (PRESENT(dirOUT)) THEN
      CALL READ_HDF5(STRARR1D=path1d, FILE=filmas, NAME='output dir')
      dirOUT = path1d(1)
    END IF

    !! Stopping condition
    IF (PRESENT(dostop)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='dostop')
      dostop = TRIMEQ(strarr1d(1),'T')
    END IF
    
    !! Default values
    ! Nmcmc0 = 1000000
    ! verbose0 = verbatim
    ! robust_RMS0 = .FALSE.
    ! robust_cal0 = .FALSE.
    ! skew_RMS0 = .FALSE.
    ! NiniMC0 = 0
    ! calib0 = .FALSE.
    ! newseed0 = .FALSE.
    ! newinit0 = .FALSE.
    
    !! input_master
    IF (PRESENT(resume)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='resume')
      resume = TRIMEQ(strarr1d(1),'T')
    END IF
    IF (PRESENT(indresume)) THEN
      CALL READ_HDF5(INTARR1D=intarr1d, FILE=filmas, NAME='indresume')
      indresume = intarr1d(1)
    END IF
    IF (PRESENT(Nmcmc)) THEN
      CALL READ_HDF5(INTARR1D=intarr1d, FILE=filmas, NAME='Nmcmc')
      Nmcmc = intarr1d(1)
    END IF
    IF (PRESENT(verbose)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='verbose')
      verbose = TRIMEQ(strarr1d(1),'T')
    END IF
    IF (PRESENT(NiniMC)) THEN
      CALL READ_HDF5(INTARR1D=intarr1d, FILE=filmas, NAME='NiniMC')
      NiniMC = intarr1d(1)
    END IF
    IF (PRESENT(calib)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='calib')
      calib = TRIMEQ(strarr1d(1),'T')
    END IF
    IF (PRESENT(robust_RMS)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='robust_RMS')
      robust_RMS = TRIMEQ(strarr1d(1),'T')
    END IF
    IF (PRESENT(robust_cal)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='robust_cal')
      robust_cal = TRIMEQ(strarr1d(1),'T')
    END IF
    IF (PRESENT(skew_RMS)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='skew_RMS')
      skew_RMS = TRIMEQ(strarr1d(1),'T')
    END IF
    IF (PRESENT(newseed)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='newseed')
      newseed = TRIMEQ(strarr1d(1),'T')
    END IF
    IF (PRESENT(newinit)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='newinit')
      newinit = TRIMEQ(strarr1d(1),'T')
    END IF
    IF (PRESENT(nohi)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='nohi')
      nohi = TRIMEQ(strarr1d(1),'T')
    END IF
    
    !! input_model
    IF (PRESENT(spec_unit)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='spectral unit')
      spec_unit = strarr1d(1)
    END IF

    !! Get Ncont
    CALL READ_HDF5(STRARR1D=labQ0, FILE=filmod, NAME='label cont')
    Ncont0 = SIZE(labQ0(:))
    IF (PRESENT(labQ)) labQ = labQ0(:)

    !! Get Nline
    CALL READ_HDF5(STRARR1D=labL0, FILE=filmod, NAME='label line')
    Nline0 = SIZE(labL0(:))
    IF (PRESENT(labL)) labL = labL0(:)
  
    !! Get Nband
    CALL READ_HDF5(STRARR1D=labB0, FILE=filmod, NAME='label band')
    Nband0 = SIZE(labB0(:))
    IF (PRESENT(labB)) labB = labB0(:)
  
    !! Get Nextc
    CALL READ_HDF5(STRARR1D=labE0, FILE=filmod, NAME='label extc')
    Nextc0 = SIZE(labE0(:))
    IF (PRESENT(labE)) labE = labE0(:)
    
    !! Get Nstar
    Nstar0 = 1
    
    IF (PRESENT(Ncont)) Ncont = Ncont0
    IF (PRESENT(Nband)) Nband = Nband0
    IF (PRESENT(Nline)) Nline = Nline0
    IF (PRESENT(Nextc)) Nextc = Nextc0
    IF (PRESENT(Nstar)) Nstar = Nstar0

    !! Get ref band
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='ref band')
    refB0 = strarr1d(1)
    IF (PRESENT(refB)) refB = refB0

    !! Get ref wavelength
    CALL READ_HDF5(DBLARR1D=dblarr1d, FILE=filmod, NAME='ref wavelength')
    refw0 = dblarr1d(1)
    IF (PRESENT(refw)) refw = refw0
    
    !! Read optical properties (Qabs)
    !!--------------------------------
    IF (PRESENT(Qabs)) THEN
      !! Qabs % rho(kg/m3); wave(micron); nu(Hz); kappa(m2/kg)
      CALL READ_OPTICS(labQ0, WAVE=wave0, RADIUSALL=radius, QABSALL=Qabsall)    
      
      wave = wave0(:)
      IF (PRESENT(wavAll)) wave = wavAll(:)
      Nw = SIZE(wave(:))

      ALLOCATE(indrad(Ncont0), Qabs(Ncont0))
      
      DO i=1,Ncont0
        ALLOCATE(Qabs(i)%wave(Nw),Qabs(i)%nu(Nw),Qabs(i)%kappa(Nw))
        
        Nr = SIZE(radius(i,:))
        Qabs(i)%rho = rho_grain(labQ0(i))
        Qabs(i)%wave(:) = wave(:)
        Qabs(i)%nu(:) = MKS%clight/MKS%micron / wave(:)
        indrad(i) = CLOSEST(radius(i,:), a0)
        ! PRINT*, 'Radius of ', labQ0(i), ': ', radius(i,indrad(i)), '->', indrad(i)
        IF (PRESENT(wavAll)) THEN
          Qabs(i)%kappa(:) = 3._DP/4._DP/rho_grain(labQ0(i))/radius(i,indrad(i)) * &
                             INTERP_LIN_SORTED( Qabsall(i,indrad(i),:), &
                               wave0(:),wave(:),XLOG=.TRUE.,YLOG=.TRUE.,FORCE=.TRUE. )
        ELSE
          Qabs(i)%kappa(:) = 3._DP/4._DP/rho_grain(labQ0(i))/radius(i,indrad(i)) * &
                             Qabsall(i,indrad(i),:)

        END IF
      END DO
      !! Free memory space
      DEALLOCATE(wave0, radius, Qabsall)
      
    END IF

    !! Read extinction curves
    !!------------------------
    IF (PRESENT(extinct)) THEN
      IF (PRESENT(wavAll)) THEN
        wave = wavAll(:)
        Nw = SIZE(wave(:))
        
        ALLOCATE(extinct(Nextc0,Nw))
        
        DO i=1,Nextc0
          extinct(i,:) = extCurve(labE0(i),wave(:),.TRUE.)
        
        END DO
      ELSE
        CALL STRIKE('READ_MASTER','Input wavelength grid.')

      END IF
    END IF

    !! Extra parameters
    !!------------------
    Nextra0 = 0

    !! input_extra
    IF (PRESENT(Nextra)) THEN
    !   CALL READ_HDF5(INTARR1D=intarr1d, FILE=filext, NAME='Nextra')
    !   Nextra0 = intarr1d(1)
      Nextra = Nextra0
    END IF

    ! setEXTRA: IF (Nextra0 > 0) THEN
    !   ALLOCATE (parinfoextra(Nextra0))
    !   CALL READ_PARTUNING('extra',parinfoextra)
    ! ELSE
    !   ALLOCATE (parinfoextra(1))
    ! END IF setEXTRA
    
    ! IF (PRESENT(parextinfo) .AND. Nextra0 > 0) THEN
    !   ALLOCATE(parextinfo(Nextra0))
    !   parextinfo(:) = parinfoextra(:)
    ! END IF

    !! Rearrange the data (parinfo)
    !!------------------------------
    Nparall = 2*Ncont0 + 4*Nband0 + 3*Nline0 + Nextc0 + Nstar0 + Nextra0

    ALLOCATE(parinfall(Nparall), boolpar(Nparall), boolparmod(Nparall))

    !! Read input model parameters
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo name')
    parinfall(1:Nparall-Nextra0)%name = strarr1d(:)
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo comp')
    parinfall(1:Nparall-Nextra0)%comp = strarr1d(:)
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo fixed')
    parinfall(1:Nparall-Nextra0)%fixed = TRIMEQ(strarr1d(:),'T')
    CALL READ_HDF5(STRARR2D=strarr2d, FILE=filmod, NAME='parinfo limited')
    CALL READ_HDF5(DBLARR2D=dblarr2d, FILE=filmod, NAME='parinfo limits')
    DO i=1,2
      parinfall(1:Nparall-Nextra0)%limited(i) = TRIMEQ(strarr2d(i,:),'T')
      parinfall(1:Nparall-Nextra0)%limits(i) = dblarr2d(i,:)
    END DO
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo model')
    parinfall(1:Nparall-Nextra0)%model = TRIMEQ(strarr1d(:),'T')
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo hyper')
    parinfall(1:Nparall-Nextra0)%hyper = TRIMEQ(strarr1d(:),'T')
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo tied')
    parinfall(1:Nparall-Nextra0)%tied = strarr1d(:)
    CALL READ_HDF5(DBLARR1D=dblarr1d, FILE=filmod, NAME='parinfo value')
    parinfall(1:Nparall-Nextra0)%value = dblarr1d(:)

    !! Read input extra parameters
    ! CALL READ_HDF5(STRARR1D=strarr1d, FILE=filext, NAME='parextinfo name')
    ! parinfall(Nparall-Nextra0+1:)%name = strarr1d(:)

    boolparmod(:) = .FALSE.
    boolpar(:) = .TRUE.
    boolparmod(:) = parinfall(1:Nparall-Nextra0)%model
    boolpar(:) = parinfall(1:Nparall-Nextra0)%hyper

    !! Number of parameters for the modelling
    Nparmod0 = COUNT(boolparmod(:))
    IF (PRESENT(Nparmod)) Nparmod = Nparmod0
    
    !! Number of parameters for the modelling,
    !! including the fixed ones, plus the parameters to be correlated
    Npar0 = Nparmod0 + Nextra0
    IF (PRESENT(Npar)) Npar = Npar0
    !! Nparall is the number of param calculated by the full model
    !! Technically any certain param can be deactivated by parinfo%model=.FALSE.
    !! This happens when there are tied param
    !! Otherwise Npar0=Nparall (thus Nparall-Nextra0=Nparmod0) is valid

    !! General parameter structure
    ALLOCATE (parinfo0(Npar0))
    parinfo0(:) = PACK(parinfall(:), boolpar(:) .OR. boolparmod(:))
    parinfo0(:)%ind = [(i,i=1,Npar0)]
    IF (PRESENT(parinfo)) THEN
      ALLOCATE (parinfo(Npar0)) 
      parinfo(:) = parinfo0(:)
    END IF

    !! No correlation for tied parameters
    boolpar(1:Nparmod0) = ( boolpar(1:Nparmod0) &
                            .AND. TRIMEQ(parinfall(1:Nparmod0)%tied, "") )
    WHERE ( .NOT. TRIMEQ(parinfall(:)%tied, "") )
      parinfall(:)%hyper = .FALSE.
    END WHERE

    !! Number of parameters to be correlated (hyperstructure)
    Nparhyp0 = COUNT(boolpar(:))
    IF (Nparhyp0 == 0) &
      Nparhyp0 = COUNT(.NOT. parinfall(:)%fixed .AND. boolparmod(:))
    IF (PRESENT(Nparhyp)) Nparhyp = Nparhyp0

    !! Number of non-trivial correlations (C(N,2))
    Ncorrhyp0 = N_CORR(Nparhyp0)
    IF (PRESENT(Ncorrhyp)) Ncorrhyp = Ncorrhyp0
    Ncorr0 = N_CORR(Npar0)
    IF (PRESENT(Ncorr)) Ncorr = Ncorr0

    !! Parameter structures
    IF (PRESENT(parhypinfo)) THEN
      ALLOCATE (parhypinfo(Nparhyp0)) 
      IF (COUNT(boolpar(:)) > 0) THEN
        parhypinfo(:) = PACK(parinfall(:), boolpar(:))
      ELSE
        parhypinfo(:) = PACK(parinfall(:), &
                             .NOT. parinfall(:)%fixed .AND. boolparmod(:))
      END IF
    END IF
    IF (PRESENT(parmodinfo)) THEN
      ALLOCATE (parmodinfo(Nparmod0)) 
      parmodinfo(:) = PACK(parinfall(:), boolparmod(:))
    END IF

    !! Correlation names
    IF (PRESENT(corrhypname)) THEN
      ALLOCATE (parhypname(Nparhyp0))
      IF (COUNT(boolpar(:)) > 0) THEN
        parhypname(:) = PACK(parinfall(:)%name,boolpar(:))
        CALL CORREL_PARLIST(parhypname(:),corrhypname)
      ELSE
        parhypname(:) = PACK(parinfall(:)%name, &
                             .NOT. parinfall(:)%fixed .AND. boolparmod(:))
        CALL CORREL_PARLIST(parhypname(:),corrhypname)
      END IF
      DEALLOCATE (parhypname)
    END IF

    !! Correlation names for all parameters
    IF (PRESENT(corrname)) CALL CORREL_PARLIST(parinfall(:)%name,corrname)
    
    !! Indices
    IF (PRESENT(indpar)) THEN
      !! Convert refw from value to index
      i = 1
      IF (PRESENT(wavALL)) THEN
        wave = wavALL(:)
        i = CLOSEST(wave(:),refw0)
      END IF
      
      CALL SET_INDPAR(INDPAR=indpar, PARINFO=parinfo0(:), &
                      REFB=refB0, LABB=labB0, REFW=i, LABQ=labQ0)
    END IF
      
    !! Free memory space
    !!-------------------
    DEALLOCATE (labQ0, labL0, labB0, labE0, &
                boolpar, boolparmod, parinfo0, parinfall)
    
  END SUBROUTINE read_master

  !!-------------------------------------------------------
  !!
  !!     Automatic initialization of model parameters
  !!
  !!-------------------------------------------------------
  SUBROUTINE initparam( NiniMC, ind, par, parinfo, itied, mask, &
                        newinit, filobs, labB, labL, Qabs )
    
    USE auxil, ONLY: TABLine, TABand
    USE utilities, ONLY: DP, trimeq, trimlr, pring, isNaN
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: reallocate, iwhere, incrarr, closest, reverse
    USE inout, ONLY: read_hdf5, lenpar
    USE statistics, ONLY: mean
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NiniMC
    TYPE(indpar_type), INTENT(IN) :: ind
    REAL(DP), DIMENSION(:,:,:,:), INTENT(INOUT) :: par ! (Nx,Ny,Npar,NiniMC)
    TYPE(parinfo_type), DIMENSION(:), INTENT(INOUT) :: parinfo
    INTEGER, DIMENSION(:), INTENT(IN) :: itied
    LOGICAL, DIMENSION(:,:,:), INTENT(INOUT) :: mask ! maskpar(Nx,Ny,Npar)
    LOGICAL, INTENT(IN), OPTIONAL :: newinit
    CHARACTER(*), INTENT(IN), OPTIONAL :: filobs
    CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: labB, labL
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN), OPTIONAL :: Qabs
    
    INTEGER :: i, j, k, ipar, iw, x, y, Nx, Ny, Nw, Nmc, NwCONT, itab, indref
    INTEGER :: Npar, Ncont, Nline, Nband, Nextc, Nstar!, Nextra
    INTEGER, DIMENSION(:), ALLOCATABLE :: indw
    INTEGER, DIMENSION(SIZE(par,1),SIZE(par,2)) :: iw2
    REAL(DP) :: limi, lims, val, sig, Tstar
    REAL(DP) :: Cband, WSband, WLband, Cline, Wline
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wOBS, nuOBS
    REAL(DP), DIMENSION(MAX(NiniMC,1)) :: theta
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mu
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: iniparval, FnuOBS, val3
    REAL(DP), DIMENSION(SIZE(par,1),SIZE(par,2)) :: limi2, lims2, val2, lnT
    REAL(DP), DIMENSION(SIZE(par,1),SIZE(par,2),MAX(NiniMC,1)) :: theta3, lnIref
    LOGICAL :: newini, chi2ini
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: maskxy
    CHARACTER(lenpar) :: spec_unit
    CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: iniparname, strarr1d
    TYPE(irrintarr_type), DIMENSION(SIZE(par,1),SIZE(par,2)) :: &
      iwFIT, iwBL, iwCONT

    !!---------------
    !! Preliminaries
    !!---------------
    Nx = SIZE(par(:,:,:,:),1)
    Ny = SIZE(par(:,:,:,:),2)
    Npar = SIZE(par(:,:,:,:),3)
    ! Ncont = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'CONT')) ) / NparCONT
    ! Nline = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'LINE')) ) / NparLINE
    ! Nband = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'BAND')) ) / NparBAND
    ! Nextc = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'EXTC')) ) / NparEXTC
    ! Nstar = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'STAR')) ) / NparSTAR
    Ncont = SIZE(Qabs(:))
    Nband = COUNT(ind%lnRband(:) .NE. -1)
    Nline = COUNT(ind%lnRline(:) .NE. -1)
    Nextc = COUNT(ind%lnAv(:) .NE. -1)
    Nstar = COUNT(ind%lnFstar(:) .NE. -1)

    IF (PRESENT(newinit)) THEN
      newini = newinit
    ELSE
      newini = .FALSE.
    END IF

    Tstar = 5.E4_DP ! [K] high enough to stay in Rayleigh-Jeans limit (T>>5000K for 3um)

    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filobs, NAME='spectral unit')
    spec_unit = strarr1d(1)
    CALL READ_HDF5(DBLARR1D=wOBS, FILE=filobs, NAME='wavelength (microns)')
    CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filobs, NAME='FnuOBS ('//TRIMLR(spec_unit)//')')

    Nw = SIZE(wOBS(:))
    ALLOCATE(nuOBS(Nw), val3(Nx,Ny,Nw))
    nuOBS(:) = MKS%clight / wOBS(:) /MKS%micron

    indref = ind%refB
    lnIref(:,:,:) = 1._DP
    
    !! Mask
    ALLOCATE(maskxy(Nx,Ny))
    FORALL (x=1:Nx,y=1:Ny) maskxy(x,y) = ANY(mask(x,y,:))
    !! Find (indices of) all wavelengths where FnuOBS is not NaN
    DO x=1,Nx
      DO y=1,Ny
        IF (maskxy(x,y)) THEN
          !! iwFIT is NOT-NaN wvl index list
          CALL IWHERE( .NOT. isNaN(FnuOBS(x,y,:)), iwFIT(x,y)%arr1d )
        ELSE
          ALLOCATE(iwFIT(x,y)%arr1d(0))
        END IF
      END DO
    END DO
    
    newinipar: IF (.NOT. newini) THEN
      
      !!----------------------------------------------------
      !! I. Automatically generate initial parameter values
      !!----------------------------------------------------

      IF (NiniMC==0) THEN

        !! a. Simple initilization
        !!-------------------------

        DO i=1,Nband
          CALL IWHERE( TRIMEQ(TABand(:)%label, labB(i)), itab ) ! itab = ind of TABand
          Cband = TABand(itab)%wave
          WSband = TABand(itab)%sigmaS + degradeRes(Cband, .01_DP, 'NG-SL-LL')
          WLband = TABand(itab)%sigmaL + degradeRes(Cband, .01_DP, 'NG-SL-LL')
          
          !! Cband
          IF (parinfo(ind%Cband(i))%fixed) THEN
            par(:,:,ind%Cband(i),:) = parinfo(ind%Cband(i))%value
          ELSE
            par(:,:,ind%Cband(i),:) = Cband
            !! If Cband is not fixed, limits are indispensable
            sig = Cband_sig
            IF (.NOT. parinfo(ind%Cband(i))%limited(1)) &
              parinfo(ind%Cband(i))%limits(1) = Cband-sig
            IF (.NOT. parinfo(ind%Cband(i))%limited(2)) &
              parinfo(ind%Cband(i))%limits(2) = Cband+sig
            parinfo(ind%Cband(i))%limited = .TRUE.
          END IF
          
          !! WSband
          IF (parinfo(ind%WSband(i))%fixed) THEN
            par(:,:,ind%WSband(i),:) = parinfo(ind%WSband(i))%value
          ELSE
            par(:,:,ind%WSband(i),:) = WSband
            IF (.NOT. parinfo(ind%WSband(i))%limited(1)) &
              parinfo(ind%WSband(i))%limits(1) = WSband*0.5_DP
            IF (.NOT. parinfo(ind%WSband(i))%limited(2)) &
              parinfo(ind%WSband(i))%limits(2) = WSband*2._DP
            parinfo(ind%WSband(i))%limited = .TRUE.
          END IF
     
          !! WLband
          IF (parinfo(ind%WLband(i))%fixed) THEN
            par(:,:,ind%WLband(i),:) = parinfo(ind%WLband(i))%value
          ELSE
            par(:,:,ind%WLband(i),:) = WLband
            IF (.NOT. parinfo(ind%WLband(i))%limited(1)) &
              parinfo(ind%WLband(i))%limits(1) = WLband*0.5_DP
            IF (.NOT. parinfo(ind%WLband(i))%limited(2)) &
              parinfo(ind%WLband(i))%limits(2) = WLband*2._DP
            parinfo(ind%WLband(i))%limited = .TRUE.
          END IF
          
          !! Check wavelength mask to decide if mask par
          DO x=1,Nx
            DO y=1,Ny
              IF (maskxy(x,y)) THEN
                !! Update par
                Cband = par(x,y,ind%Cband(i),1)
                WSband = par(x,y,ind%WSband(i),1)
                WLband = par(x,y,ind%WLband(i),1)
                limi = Cband - WSband
                lims = Cband + WLband
                IF ( ALL( wOBS(iwFIT(x,y)%arr1d(:)).LE.limi &
                          .OR. wOBS(iwFIT(x,y)%arr1d(:)).GE.lims ) ) THEN
                  mask(x,y,ind%lnRband(i)) = .FALSE.
                  mask(x,y,ind%Cband(i)) = .FALSE.
                  mask(x,y,ind%WSband(i)) = .FALSE.
                  mask(x,y,ind%WLband(i)) = .FALSE.
                  !! Update maskxy
                  maskxy(x,y) = ANY(mask(x,y,:))
                ELSE
                  CALL IWHERE( wOBS(:)>limi .AND. wOBS(:)<lims, indw )
                  DO iw=1,SIZE(indw)
                    CALL INCRARR( iwBL(x,y)%arr1d, indw(iw) )
                  END DO
                  DEALLOCATE(indw)
                END IF
              END IF
            END DO
          END DO
          
          !! lnRband
          IF (parinfo(ind%lnRband(i))%fixed) THEN
            par(:,:,ind%lnRband(i),:) = parinfo(ind%lnRband(i))%value
          ELSE
            FORALL (x=1:Nx,y=1:Ny,mask(x,y,ind%lnRband(i)))
              iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cband(i),1) )
              !! Same estimation as lnRline
              par(x,y,ind%lnRband(i),:) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) / &
                MAXVAL(LORENTZBAND(wOBS(:),par(x,y,ind%Cband(i),1), &
                                   par(x,y,ind%WSband(i),1), &
                                   par(x,y,ind%WLband(i),1),.TRUE.)) )
            END FORALL

            IF (i == indref) lnIref(:,:,:) = par(:,:,ind%lnRband(i),:)
          END IF

        END DO

        DO i=1,Nband
          IF (i /= indref) &
            par(:,:,ind%lnRband(i),:) = par(:,:,ind%lnRband(i),:) - lnIref(:,:,:)
        END DO
        
        DO i=1,Nline
          CALL IWHERE( TRIMEQ(TABLine(:)%label, labL(i)), itab ) ! itab = ind of TABLine
          Cline = TABLine(itab)%wave
          Wline = degradeRes(Cline, .01_DP, 'NG-SL-LL')

          !! Cline
          IF (parinfo(ind%Cline(i))%fixed) THEN
            par(:,:,ind%Cline(i),:) = parinfo(ind%Cline(i))%value
          ELSE
            par(:,:,ind%Cline(i),:) = Cline
            !! If Cline is not fixed, limits are indispensable
            !! FWHM = 2*sqrt(2*log(2))*sigma ~ 2.355 (suppose PSF is gaussian)
            !! R = lambda / FWHM
            !! => sigma = lambda / R / 2.355
            sig = Wline / 2.355_DP
            IF (.NOT. parinfo(ind%Cline(i))%limited(1)) &
              parinfo(ind%Cline(i))%limits(1) = Cline-sig
            IF (.NOT. parinfo(ind%Cline(i))%limited(2)) &
              parinfo(ind%Cline(i))%limits(2) = Cline+sig
            parinfo(ind%Cline(i))%limited = .TRUE.
          END IF
          
          !! Wline
          IF (parinfo(ind%Wline(i))%fixed) THEN
            par(:,:,ind%Wline(i),:) = parinfo(ind%Wline(i))%value
          ELSE
            par(:,:,ind%Wline(i),:) = Wline
            IF (.NOT. parinfo(ind%Wline(i))%limited(1)) &
              parinfo(ind%Wline(i))%limits(1) = Wline*0.5_DP
            IF (.NOT. parinfo(ind%Wline(i))%limited(2)) &
              parinfo(ind%Wline(i))%limits(2) = Wline*2._DP
            parinfo(ind%Wline(i))%limited = .TRUE.
          END IF
          
          !! Check wavelength mask to decide if mask par
          DO x=1,Nx
            DO y=1,Ny
              IF (maskxy(x,y)) THEN
                !! Update param
                Cline = par(x,y,ind%Cline(i),1)
                Wline = par(x,y,ind%Wline(i),1)
                limi = Cline - Wline
                lims = Cline + Wline
                IF ( ALL( wOBS(iwFIT(x,y)%arr1d(:)).LE.limi &
                          .OR. wOBS(iwFIT(x,y)%arr1d(:)).GE.lims ) ) THEN
                  mask(x,y,ind%lnRline(i)) = .FALSE.
                  mask(x,y,ind%Cline(i)) = .FALSE.
                  mask(x,y,ind%Wline(i)) = .FALSE.
                  !! Update maskxy
                  maskxy(x,y) = ANY(mask(x,y,:))
                ELSE
                  CALL IWHERE( wOBS(:)>limi .AND. wOBS(:)<lims, indw )
                  DO iw=1,SIZE(indw)
                    CALL INCRARR( iwBL(x,y)%arr1d, indw(iw) )
                  END DO
                  DEALLOCATE(indw)
                END IF
              END IF
            END DO
          END DO
          
          !! lnRline
          IF (parinfo(ind%lnRline(i))%fixed) THEN
            par(:,:,ind%lnRline(i),:) = parinfo(ind%lnRline(i))%value
          ELSE
            FORALL (x=1:Nx,y=1:Ny,mask(x,y,ind%lnRline(i)))
              iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cline(i),1) )
              !! lnRline is auto determined by the spectrum to fit (idem. for lnRband)
              !! Can be overestimated (should subtract other compo from FnuOBS)
              !! or underestimated (due to extinction Pabs).
              !! Here we just need a reasonable order of magnitude.
              par(x,y,ind%lnRline(i),:) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) / &
                MAXVAL(GAUSSLINE(wOBS(:),par(x,y,ind%Cline(i),1), &
                                 par(x,y,ind%Wline(i),1),.TRUE.)) ) - &
                lnIref(x,y,:)
            END FORALL
          END IF
          
        END DO

        !! Create wvl grid with indices for band & line free spectra
        wcontx: DO x=1,Nx
          wconty: DO y=1,Ny
            IF (maskxy(x,y)) THEN
              NwCONT = Nw - SIZE(iwBL(x,y)%arr1d(:))
              DO iw=1,Nw
                IF ( ALL(iwBL(x,y)%arr1d(:).NE.iw) ) THEN
                  CALL INCRARR( iwCONT(x,y)%arr1d, iw )
                END IF
              END DO
              !! Larger T corresponds to smaller wvl
              iwCONT(x,y)%arr1d(:) = REVERSE(iwCONT(x,y)%arr1d(:))
        
              DO i=1,Ncont
                !! lnT LOG[K]
                IF (parinfo(ind%lnT(i))%fixed) THEN
                  par(x,y,ind%lnT(i),:) = parinfo(ind%lnT(i))%value
                ELSE
                  !! Create lnT grid: LOG(50,500) - (Wien Law: 58-5.8um)
                  IF (ind%ordQ(i)>0) THEN
                    !! lnT(i) becomes lndT(i) = LOG( EXP(lnT(i)) - EXP(lnT(i-1)) )
                    limi = MERGE(parinfo(ind%lnT(i))%limits(1), &
                                 0._DP, &
                                 parinfo(ind%lnT(i))%limited(1))
                    lims = MERGE(parinfo(ind%lnT(i))%limits(2), &
                                 LOG(500._DP-50._DP), &
                                 parinfo(ind%lnT(i))%limited(2))
                    par(x,y,ind%lnT(i),:) = LOG( (EXP(lims)-EXP(limi))/Ncont )
                  ELSE
                    limi = MERGE(parinfo(ind%lnT(i))%limits(1), &
                                 LOG(50._DP), &
                                 parinfo(ind%lnT(i))%limited(1))
                    lims = MERGE(parinfo(ind%lnT(i))%limits(2), &
                                 LOG(500._DP), &
                                 parinfo(ind%lnT(i))%limited(2))
                    par(x,y,ind%lnT(i),:) = LOG( (EXP(lims)-EXP(limi)) &
                                                 *REAL(i)/Ncont+EXP(limi) ) !!! i/Ncont
                  END IF
                  parinfo(ind%lnT(i))%limits(1) = limi
                  parinfo(ind%lnT(i))%limits(2) = lims
                  parinfo(ind%lnT(i))%limited = .TRUE.

                END IF
          
                !! lnFcont LOG[W/m2/sr]
                IF ( parinfo(ind%lnFcont(i))%fixed ) THEN
                  par(x,y,ind%lnFcont(i),:) = parinfo(ind%lnFcont(i))%value
                ELSE
                  !! lnFcont is auto determined by modifBB(lnT) & the spectrum to fit
                  IF (NwCONT/Ncont.GE.1) THEN
                    iw = INT(NwCONT/Ncont) * i ! iw is NOT wOBS's index, but iwCONT(x,y)%arr1d's
                  ELSE
                    !! Case where there are more CONT compo than its wvl grid (Ncont>NwCONT)
                    IF (i>NwCONT) THEN
                      iw = NwCONT
                    ELSE
                      iw = i
                    END IF
                  END IF
                  k = ind%ordQ(i)
                  lnT(x,y) = par(x,y,ind%lnT(i-k),1)
                  IF (k>0) &
                    lnT(x,y) = LOG( EXP(lnT(x,y)) &
                                    +SUM(EXP(par(x,y,ind%lnT(i-k+1:i),1))) )
                  val3(x,y,:) = MODIFBB(wOBS(:),EXP(lnT(x,y)),Qabs(i),.TRUE.,ind%refw)
                  par(x,y,ind%lnFcont(i),:) = LOG( ABS(FnuOBS(x,y,iwCONT(x,y)%arr1d(iw))) &
                                                   / val3(x,y,iwCONT(x,y)%arr1d(iw)) / Ncont )
                END IF
              END DO

              DO i=1,Nstar
                !! lnFstar LOG[W/m2/sr]
                IF (parinfo(ind%lnFstar(i))%fixed) THEN
                  par(x,y,ind%lnFstar(i),:) = parinfo(ind%lnFstar(i))%value
                ELSE
                  !! lnFstar is auto determined by BB & the spectrum to fit
                  iw = NwCONT - 10 ! avoid wvl edge
                  val3(x,y,:) = BLACKBODY(nuOBS(:),Tstar) / MKS%stefan/Tstar**4
                  val2(x,y) = LOG( ABS(FnuOBS(x,y,iwCONT(x,y)%arr1d(iw))) &
                                   / val3(x,y,iwCONT(x,y)%arr1d(iw)) /pi )
                  par(x,y,ind%lnFstar(i),:) = val2(x,y)
                END IF
              END DO

            END IF
          END DO wconty
        END DO wcontx
        
        DO i=1,Nextc
          !! lnAv LOG[mag]
          IF (parinfo(ind%lnAv(i))%fixed) THEN
            par(:,:,ind%lnAv(i),:) = parinfo(ind%lnAv(i))%value
          ELSE
            par(:,:,ind%lnAv(i),:) = 0._DP ! 1 [mag] = no extinction
            IF (.NOT. parinfo(ind%lnAv(i))%limited(1)) &
              parinfo(ind%lnAv(i))%limits(1) = 0._DP
            IF (.NOT. parinfo(ind%lnAv(i))%limited(2)) &
              parinfo(ind%lnAv(i))%limits(2) = 5._DP
            parinfo(ind%lnAv(i))%limited = .TRUE.
          END IF

        END DO
        
      ELSE

        !! b. MC initialization
        !!----------------------

        DO i=1,Nband
          CALL IWHERE( TRIMEQ(TABand(:)%label, labB(i)), itab ) ! itab = ind of TABand
          Cband = TABand(itab)%wave
          WSband = TABand(itab)%sigmaS + degradeRes(Cband, .01_DP, 'NG-SL-LL')
          WLband = TABand(itab)%sigmaL + degradeRes(Cband, .01_DP, 'NG-SL-LL')
          
          !! Cband
          IF (parinfo(ind%Cband(i))%fixed) THEN
            par(:,:,ind%Cband(i),:) = parinfo(ind%Cband(i))%value
          ELSE
            ! sig = degradeRes(Cband, 0.01_DP, 'NG-SL-LL') / 2.355_DP
            sig = Cband_sig
            limi = MERGE( parinfo(ind%Cband(i))%limits(1), &
                          Cband-sig, &
                          parinfo(ind%Cband(i))%limited(1) ) ! lim inf
            lims = MERGE( parinfo(ind%Cband(i))%limits(2), &
                          Cband+sig, &
                          parinfo(ind%Cband(i))%limited(2) ) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%Cband(i),j) = Cband
              ELSE
                !! Uniform distribution
                par(:,:,ind%Cband(i),j) = (lims - limi) * theta3(:,:,j) + limi
              END IF
            END DO
            IF (.NOT. parinfo(ind%Cband(i))%limited(1)) &
              parinfo(ind%Cband(i))%limits(1) = Cband-sig
            IF (.NOT. parinfo(ind%Cband(i))%limited(2)) &
              parinfo(ind%Cband(i))%limits(2) = Cband+sig
            parinfo(ind%Cband(i))%limited = .TRUE.
          END IF
          
          !! WSband
          IF (parinfo(ind%WSband(i))%fixed) THEN
            par(:,:,ind%WSband(i),:) = parinfo(ind%WSband(i))%value
          ELSE
            limi = MERGE( parinfo(ind%WSband(i))%limits(1), &
                          WSband*0.5_DP, &
                          parinfo(ind%WSband(i))%limited(1) ) ! lim inf
            lims = MERGE( parinfo(ind%WSband(i))%limits(2), &
                          WSband*2._DP, &
                          parinfo(ind%WSband(i))%limited(2) ) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%WSband(i),j) = WSband
              ELSE
                !! Uniform distribution
                par(:,:,ind%WSband(i),j) = (lims - limi) * theta3(:,:,j) + limi
              END IF
            END DO
            IF (.NOT. parinfo(ind%WSband(i))%limited(1)) &
              parinfo(ind%WSband(i))%limits(1) = WSband*0.5_DP
            IF (.NOT. parinfo(ind%WSband(i))%limited(2)) &
              parinfo(ind%WSband(i))%limits(2) = WSband*2._DP
            parinfo(ind%WSband(i))%limited = .TRUE.
          END IF
     
          !! WLband
          IF (parinfo(ind%WLband(i))%fixed) THEN
            par(:,:,ind%WLband(i),:) = parinfo(ind%WLband(i))%value
          ELSE
            limi = MERGE( parinfo(ind%WLband(i))%limits(1), &
                          WLband*0.5_DP, &
                          parinfo(ind%WLband(i))%limited(1) ) ! lim inf
            lims = MERGE( parinfo(ind%WLband(i))%limits(2), &
                          WLband*2._DP, &
                          parinfo(ind%WLband(i))%limited(2) ) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%WLband(i),j) = WLband
              ELSE
                !! Uniform distribution
                par(:,:,ind%WLband(i),j) = (lims - limi) * theta3(:,:,j) + limi
              END IF
            END DO
            IF (.NOT. parinfo(ind%WLband(i))%limited(1)) &
              parinfo(ind%WLband(i))%limits(1) = limi
            IF (.NOT. parinfo(ind%WLband(i))%limited(2)) &
              parinfo(ind%WLband(i))%limits(2) = lims
            parinfo(ind%WLband(i))%limited = .TRUE.
          END IF
          
          !! Check wavelength mask to decide if mask par
          DO x=1,Nx
            DO y=1,Ny
              IF (maskxy(x,y)) THEN
                !! Update par
                Cband = par(x,y,ind%Cband(i),1)
                WSband = par(x,y,ind%WSband(i),1)
                WLband = par(x,y,ind%WLband(i),1)
                limi = Cband - WSband
                lims = Cband + WLband
                IF ( ALL( wOBS(iwFIT(x,y)%arr1d(:)).LE.limi &
                          .OR. wOBS(iwFIT(x,y)%arr1d(:)).GE.lims ) ) THEN
                  mask(x,y,ind%lnRband(i)) = .FALSE.
                  mask(x,y,ind%Cband(i)) = .FALSE.
                  mask(x,y,ind%WSband(i)) = .FALSE.
                  mask(x,y,ind%WLband(i)) = .FALSE.
                  !! Update maskxy
                  maskxy(x,y) = ANY(mask(x,y,:))
                ELSE
                  CALL IWHERE( wOBS(:)>limi .AND. wOBS(:)<lims, indw )
                  DO iw=1,SIZE(indw)
                    CALL INCRARR( iwBL(x,y)%arr1d, indw(iw) )
                  END DO
                  DEALLOCATE(indw)
                END IF
              END IF
            END DO
          END DO
          
          !! lnRband
          IF (parinfo(ind%lnRband(i))%fixed) THEN
            par(:,:,ind%lnRband(i),:) = parinfo(ind%lnRband(i))%value
          ELSE
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                FORALL (x=1:Nx,y=1:Ny,mask(x,y,ind%lnRband(i)))
                  iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cband(i),j) )
                  !! Estimated via feature peak
                  par(x,y,ind%lnRband(i),j) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) / &
                    MAXVAL(LORENTZBAND(wOBS(:),par(x,y,ind%Cband(i),j), &
                                       par(x,y,ind%WSband(i),j), &
                                       par(x,y,ind%WLband(i),j),.TRUE.)) )
  
                END FORALL
              ELSE
                FORALL (x=1:Nx,y=1:Ny,mask(x,y,ind%lnRband(i)))
                  iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cband(i),j) )
                  !! Estimated via feature peak
                  val2(x,y) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) / &
                    MAXVAL(LORENTZBAND(wOBS(:),par(x,y,ind%Cband(i),j), &
                                       par(x,y,ind%WSband(i),j), &
                                       par(x,y,ind%WLband(i),j),.TRUE.)) )
                  limi2(x,y) = MERGE( parinfo(ind%lnRband(i))%limits(1), &
                                      val2(x,y)-.1_DP, & ! EXP(factor - 3)
                                      parinfo(ind%lnRband(i))%limited(1) ) ! lim inf
                  lims2(x,y) = MERGE( parinfo(ind%lnRband(i))%limits(2), &
                                      val2(x,y), &
                                      parinfo(ind%lnRband(i))%limited(2) ) ! lim sup
                  !! Uniform distribution
                  par(x,y,ind%lnRband(i),j) = (lims2(x,y)-limi2(x,y)) * &
                                              theta3(x,y,j) + limi2(x,y)
                  
                END FORALL
              END IF
            END DO
            
            IF (i == indref) lnIref(:,:,:) = par(:,:,ind%lnRband(i),:)
          END IF
          
        END DO

        DO i=1,Nband
          IF (i /= indref) &
            par(:,:,ind%lnRband(i),:) = par(:,:,ind%lnRband(i),:) - lnIref(:,:,:)
        END DO
        
        DO i=1,Nline
          CALL IWHERE( TRIMEQ(TABLine(:)%label, labL(i)), itab ) ! itab = ind of TABLine
          Cline = TABLine(itab)%wave
          Wline = degradeRes(Cline, .01_DP, 'NG-SL-LL')
          
          !! Cline
          IF (parinfo(ind%Cline(i))%fixed) THEN
            par(:,:,ind%Cline(i),:) = parinfo(ind%Cline(i))%value
          ELSE
            !! FWHM = 2*sqrt(2*log(2))*sigma ~ 2.355 (suppose PSF is gaussian)
            !! R = lambda / FWHM
            !! => sigma = lambda / R / 2.355
            sig = Wline / 2.355_DP
            limi = MERGE( parinfo(ind%Cline(i))%limits(1), &
                          Cline-sig, &
                          parinfo(ind%Cline(i))%limited(1)) ! lim inf
            lims = MERGE( parinfo(ind%Cline(i))%limits(2), &
                          Cline+sig, &
                          parinfo(ind%Cline(i))%limited(2)) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%Cline(i),j) = Cline
              ELSE
                !! Uniform distribution
                par(:,:,ind%Cline(i),j) = (lims - limi) * theta3(:,:,j) + limi
              END IF
            END DO
            IF (.NOT. parinfo(ind%Cline(i))%limited(1)) &
              parinfo(ind%Cline(i))%limits(1) = Cline-sig
            IF (.NOT. parinfo(ind%Cline(i))%limited(2)) &
              parinfo(ind%Cline(i))%limits(2) = Cline+sig
            parinfo(ind%Cline(i))%limited = .TRUE.
          END IF
          
          !! Wline
          IF (parinfo(ind%Wline(i))%fixed) THEN
            par(:,:,ind%Wline(i),:) = parinfo(ind%Wline(i))%value
          ELSE
            limi = MERGE( parinfo(ind%Wline(i))%limits(1), &
                          Wline*0.5_DP, &
                          parinfo(ind%Wline(i))%limited(1)) ! lim inf
            lims = MERGE( parinfo(ind%Wline(i))%limits(2), &
                          Wline*2._DP, &
                          parinfo(ind%Wline(i))%limited(2)) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%Wline(i),j) = Wline
              ELSE
                !! Uniform distribution
                par(:,:,ind%Wline(i),j) = (lims - limi) * theta3(:,:,j) + limi
              END IF
            END DO
            IF (.NOT. parinfo(ind%Wline(i))%limited(1)) &
              parinfo(ind%Wline(i))%limits(1) = Wline*0.5_DP
            IF (.NOT. parinfo(ind%Wline(i))%limited(2)) &
              parinfo(ind%Wline(i))%limits(2) = Wline*2._DP
            parinfo(ind%Wline(i))%limited = .TRUE.
          END IF
          
          !! Check wavelength mask to decide if mask par
          DO x=1,Nx
            DO y=1,Ny
              IF (maskxy(x,y)) THEN
                !! Update param
                Cline = par(x,y,ind%Cline(i),1)
                Wline = par(x,y,ind%Wline(i),1)
                limi = Cline - Wline
                lims = Cline + Wline
                IF ( ALL( wOBS(iwFIT(x,y)%arr1d(:)).LE.limi &
                          .OR. wOBS(iwFIT(x,y)%arr1d(:)).GE.lims ) ) THEN
                  mask(x,y,ind%lnRline(i)) = .FALSE.
                  mask(x,y,ind%Cline(i)) = .FALSE.
                  mask(x,y,ind%Wline(i)) = .FALSE.
                  !! Update maskxy
                  maskxy(x,y) = ANY(mask(x,y,:))
                ELSE
                  CALL IWHERE( wOBS(:)>limi .AND. wOBS(:)<lims, indw )
                  DO iw=1,SIZE(indw)
                    CALL INCRARR( iwBL(x,y)%arr1d, indw(iw) )
                  END DO
                  DEALLOCATE(indw)
                END IF
              END IF
            END DO
          END DO
          
          !! lnRline
          IF (parinfo(ind%lnRline(i))%fixed) THEN
            par(:,:,ind%lnRline(i),:) = parinfo(ind%lnRline(i))%value
          ELSE
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                FORALL (x=1:Nx,y=1:Ny,mask(x,y,ind%lnRline(i)))
                  iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cline(i),j) )
                  !! Estimated via feature peak
                  par(x,y,ind%lnRline(i),j) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) / &
                    MAXVAL(GAUSSLINE(wOBS(:),par(x,y,ind%Cline(i),j), &
                                     par(x,y,ind%Wline(i),j),.TRUE.)) ) - &
                    lnIref(x,y,j)

                END FORALL
              ELSE
                FORALL (x=1:Nx,y=1:Ny,mask(x,y,ind%lnRline(i)))
                  iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cline(i),j) )
                  !! Estimated via feature peak
                  val2(x,y) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) / &
                    MAXVAL(GAUSSLINE(wOBS(:),par(x,y,ind%Cline(i),j), &
                                     par(x,y,ind%Wline(i),j),.TRUE.)) ) - &
                    lnIref(x,y,j)
                  limi2(x,y) = MERGE( parinfo(ind%lnRline(i))%limits(1), &
                                      val2(x,y)-.1_DP, & ! EXP(factor - 3)
                                      parinfo(ind%lnRline(i))%limited(1)) ! lim inf
                  lims2(x,y) = MERGE( parinfo(ind%lnRline(i))%limits(2), &
                                      val2(x,y), &
                                      parinfo(ind%lnRline(i))%limited(2)) ! lim sup
                  !! Uniform distribution
                  par(x,y,ind%lnRline(i),j) = (lims2(x,y)-limi2(x,y)) * &
                                              theta3(x,y,j) + limi2(x,y)
                  
                END FORALL
              END IF
            END DO
          END IF
          
        END DO

        !! Create wvl grid with indices for band & line free spectra
        wcontx_mc: DO x=1,Nx
          wconty_mc: DO y=1,Ny
            IF (maskxy(x,y)) THEN
              NwCONT = Nw - SIZE(iwBL(x,y)%arr1d(:))
              DO iw=1,Nw
                IF ( ALL(iwBL(x,y)%arr1d(:).NE.iw) ) THEN
                  CALL INCRARR( iwCONT(x,y)%arr1d, iw )
                END IF
              END DO
              !! Larger T corresponds to smaller wvl
              iwCONT(x,y)%arr1d(:) = REVERSE(iwCONT(x,y)%arr1d(:))
        
              DO i=1,Ncont
                !! lnT LOG[K]
                IF (parinfo(ind%lnT(i))%fixed) THEN
                  par(x,y,ind%lnT(i),:) = parinfo(ind%lnT(i))%value
                ELSE
                  !! Generate a random value on lnT grid: LOG(50,500)
                  IF (ind%ordQ(i)>0) THEN
                    !! lnT(i) becomes lndT(i) = LOG( EXP(lnT(i)) - EXP(lnT(i-1)) )
                    limi = MERGE(parinfo(ind%lnT(i))%limits(1), &
                                 0._DP, &
                                 parinfo(ind%lnT(i))%limited(1))
                    lims = MERGE(parinfo(ind%lnT(i))%limits(2), &
                                 LOG(500._DP-50._DP), &
                                 parinfo(ind%lnT(i))%limited(2))
                  ELSE
                    limi = MERGE(parinfo(ind%lnT(i))%limits(1), &
                                 LOG(50._DP), &
                                 parinfo(ind%lnT(i))%limited(1))
                    lims = MERGE(parinfo(ind%lnT(i))%limits(2), &
                                 LOG(500._DP), &
                                 parinfo(ind%lnT(i))%limited(2))
                  END IF
                  parinfo(ind%lnT(i))%limits(1) = limi
                  parinfo(ind%lnT(i))%limits(2) = lims
                  parinfo(ind%lnT(i))%limited = .TRUE.
                  CALL RANDOM_NUMBER(theta(:))
                  IF (ind%ordQ(i)>0) THEN
                    FORALL (j=1:NiniMC) &
                      par(x,y,ind%lnT(i),j) = LOG( &
                        (EXP(lims)-EXP(limi))/Ncont * theta(j) + EXP(limi) )
                  ELSE
                    FORALL (j=1:NiniMC) &
                      par(x,y,ind%lnT(i),j) = LOG( &
                        (EXP(lims)-EXP(limi))/Ncont * theta(j) + EXP(limi) &
                        + (EXP(lims)-EXP(limi))* REAL(i-1)/Ncont ) !!! i/Ncont
                  END IF
                END IF
          
                !! lnFcont LOG[W/m2/sr]
                IF (parinfo(ind%lnFcont(i))%fixed) THEN
                  par(x,y,ind%lnFcont(i),:) = parinfo(ind%lnFcont(i))%value
                ELSE
                  CALL RANDOM_NUMBER(theta(:))
                  !! lnFcont is auto determined by modifBB(lnT) & the spectrum to fit
                  IF (Nw/Ncont.GE.1) THEN
                    iw = INT(NwCONT/Ncont) * i ! iw is NOT wOBS's index, but iwCONT(x,y)%arr1d
                  ELSE
                    !! Case where there are more CONT compo than its wvl grid (Ncont>NwCONT)
                    IF (i>NwCONT) THEN
                      iw = NwCONT
                    ELSE
                      iw = i
                    END IF
                  END IF
                  k = ind%ordQ(i)
                  lnT(x,y) = par(x,y,ind%lnT(i-k),1)
                  IF (k>0) &
                    lnT(x,y) = LOG( EXP(lnT(x,y)) &
                                    +SUM(EXP(par(x,y,ind%lnT(i-k+1:i),1))) )
                  DO j=1,NiniMC
                    IF (j == 1) THEN
                        val3(x,y,:) = MODIFBB(wOBS(:),EXP(lnT(x,y)),Qabs(i),.TRUE.,ind%refw)
                        val2(x,y) = LOG( ABS(FnuOBS(x,y,iwCONT(x,y)%arr1d(iw))) &
                                         / val3(x,y,iwCONT(x,y)%arr1d(iw)) )
                        par(x,y,ind%lnFcont(i),j) = val2(x,y)
                
                    ELSE
                      val3(x,y,:) = MODIFBB(wOBS(:),EXP(lnT(x,y)),Qabs(i),.TRUE.,ind%refw)
                      val2(x,y) = LOG( ABS(FnuOBS(x,y,iwCONT(x,y)%arr1d(iw))) &
                                       / val3(x,y,iwCONT(x,y)%arr1d(iw)) )
                      limi2(x,y) = MERGE( parinfo(ind%lnFcont(i))%limits(1), &
                                          val2(x,y), &
                                          parinfo(ind%lnFcont(i))%limited(1) ) ! lim inf
                      lims2(x,y) = MERGE( parinfo(ind%lnFcont(i))%limits(2), &
                                          val2(x,y), &
                                          parinfo(ind%lnFcont(i))%limited(2) ) ! lim sup
                      !! Uniform distribution
                      par(x,y,ind%lnFcont(i),j) = (lims2(x,y) - limi2(x,y)) * theta3(x,y,j) &
                                                   + limi2(x,y)
                    END IF
                  END DO
                END IF
              END DO
              
              DO i=1,Nstar
                !! lnFstar LOG[W/m2/sr]
                IF (parinfo(ind%lnFstar(i))%fixed) THEN
                  par(x,y,ind%lnFstar(i),:) = parinfo(ind%lnFstar(i))%value
                ELSE
                  CALL RANDOM_NUMBER(theta(:))
                  !! lnFstar is auto determined by BB & the spectrum to fit
                  iw = NwCONT - 10 ! avoid wvl edge
                  DO j=1,NiniMC
                    IF (j == 1) THEN
                      val3(x,y,:) = BLACKBODY(nuOBS(:),Tstar) / MKS%stefan/Tstar**4
                      val2(x,y) = LOG( ABS(FnuOBS(x,y,iwCONT(x,y)%arr1d(iw))) &
                                       /val3(x,y,iwCONT(x,y)%arr1d(iw)) /pi )
                      par(x,y,ind%lnFstar(i),j) = val2(x,y)
                    ELSE
                      val3(x,y,:) = BLACKBODY(nuOBS(:),Tstar) / MKS%stefan/Tstar**4
                      val2(x,y) = LOG( ABS(FnuOBS(x,y,iwCONT(x,y)%arr1d(iw))) &
                                       /val3(x,y,iwCONT(x,y)%arr1d(iw)) /pi )
                      limi2(x,y) = MERGE( parinfo(ind%lnFstar(i))%limits(1), &
                                          val2(x,y)-.1_DP, & ! EXP(factor - 3)
                                          parinfo(ind%lnFstar(i))%limited(1) ) ! lim inf
                      lims2(x,y) = MERGE( parinfo(ind%lnFstar(i))%limits(2), &
                                          val2(x,y)+.1_DP, &
                                          parinfo(ind%lnFstar(i))%limited(2) ) ! lim sup
                      !! Uniform distribution
                      par(x,y,ind%lnFstar(i),j) = (lims2(x,y)-limi2(x,y)) * &
                                                  theta3(x,y,j) + limi2(x,y)
                    END IF
                  END DO
                END IF
              END DO
                   
              END IF
            END DO wconty_mc
          END DO wcontx_mc
        
        DO i=1,Nextc
          !! lnAv LOG[mag]
          IF (parinfo(ind%lnAv(i))%fixed) THEN
            par(:,:,ind%lnAv(i),:) = parinfo(ind%lnAv(i))%value
          ELSE
            !! lnAv = 0 -> Av = 1 [mag] = no extinction
            val = 0._DP
            limi = MERGE( parinfo(ind%lnAv(i))%limits(1), &
                          val, &
                          parinfo(ind%lnAv(i))%limited(1) ) ! lim inf
            lims = MERGE( parinfo(ind%lnAv(i))%limits(2), &
                          1._DP, &
                          parinfo(ind%lnAv(i))%limited(2) ) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%lnAv(i),j) = val
              ELSE
                !! Uniform distribution
                par(:,:,ind%lnAv(i),j) = (lims - limi) * theta3(:,:,j) + limi
              END IF
            END DO
            IF (.NOT. parinfo(ind%lnAv(i))%limited(1)) &
              parinfo(ind%lnAv(i))%limits(1) = val
            IF (.NOT. parinfo(ind%lnAv(i))%limited(2)) &
              parinfo(ind%lnAv(i))%limits(2) = 5._DP
            parinfo(ind%lnAv(i))%limited = .TRUE.
          END IF

        END DO

      END IF
      DEALLOCATE(val3)
      
    ELSE
      
      !!----------------------------------------------------
      !! II. Read initial parameters from input file
      !!----------------------------------------------------

      !! Read the initial parameters from the file
      CALL READ_HDF5(STRARR1D=iniparname, FILE=filobs, &
                     NAME='Initial parameter label')
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filobs, NAME='Chi2init')
      chi2ini = TRIMEQ(strarr1d(1),'T')
      
      IF (chi2ini) THEN
        CALL READ_HDF5(DBLARR3D=iniparval, FILE=filobs, &
                       NAME='Chi2 fitted parameter value')
      ELSE
        CALL READ_HDF5(DBLARR3D=iniparval, FILE=filobs, &
                       NAME='Initial parameter value')
      END IF

      !! Set hard limits for (intensive) parameters
      !!------------------------------------------
      DO i=1,Nband
        CALL IWHERE( TRIMEQ(TABand(:)%label, labB(i)), itab ) ! itab = ind of TABand
        Cband = TABand(itab)%wave
        WSband = TABand(itab)%sigmaS + degradeRes(Cband, .01_DP, 'NG-SL-LL')
        WLband = TABand(itab)%sigmaL + degradeRes(Cband, .01_DP, 'NG-SL-LL')
        
        !! Cband
        sig = Cband_sig
        IF (.NOT. parinfo(ind%Cband(i))%limited(1)) &
          parinfo(ind%Cband(i))%limits(1) = Cband-sig
        IF (.NOT. parinfo(ind%Cband(i))%limited(2)) &
          parinfo(ind%Cband(i))%limits(2) = Cband+sig
        parinfo(ind%Cband(i))%limited = .TRUE.
        
        !! WSband
        IF (.NOT. parinfo(ind%WSband(i))%limited(1)) &
          parinfo(ind%WSband(i))%limits(1) = WSband*0.5_DP
        IF (.NOT. parinfo(ind%WSband(i))%limited(2)) &
          parinfo(ind%WSband(i))%limits(2) = WSband*2._DP
        parinfo(ind%WSband(i))%limited = .TRUE.
   
        !! WLband
        IF (.NOT. parinfo(ind%WLband(i))%limited(1)) &
          parinfo(ind%WLband(i))%limits(1) = WLband*0.5_DP
        IF (.NOT. parinfo(ind%WLband(i))%limited(2)) &
          parinfo(ind%WLband(i))%limits(2) = WLband*2._DP
        parinfo(ind%WLband(i))%limited = .TRUE.

        !! lnRband
        IF (i /= indref) THEN
          IF (.NOT. parinfo(ind%lnRband(i))%limited(1)) &
            parinfo(ind%lnRband(i))%limits(1) = -10._DP
          IF (.NOT. parinfo(ind%lnRband(i))%limited(2)) &
            parinfo(ind%lnRband(i))%limits(2) = 10._DP
          parinfo(ind%lnRband(i))%limited = .TRUE.
        END IF
        
        !! Check wavelength mask to decide if mask par
        DO x=1,Nx
          DO y=1,Ny
            IF (maskxy(x,y)) THEN
              !! Update param
              Cband = iniparval(x,y,ind%Cband(i))
              WSband = iniparval(x,y,ind%WSband(i))
              WLband = iniparval(x,y,ind%WLband(i))
              limi = Cband - WSband
              lims = Cband + WLband
              IF ( ALL( wOBS(iwFIT(x,y)%arr1d(:)).LE.limi &
                        .OR. wOBS(iwFIT(x,y)%arr1d(:)).GE.lims) ) THEN
                mask(x,y,ind%lnRband(i)) = .FALSE.
                mask(x,y,ind%Cband(i)) = .FALSE.
                mask(x,y,ind%WSband(i)) = .FALSE.
                mask(x,y,ind%WLband(i)) = .FALSE.
                !! Update maskxy
                maskxy(x,y) = ANY(mask(x,y,:))
              END IF
            END IF
          END DO
        END DO
        
      END DO
      
      DO i=1,Nline
        CALL IWHERE( TRIMEQ(TABLine(:)%label, labL(i)), itab ) ! itab = ind of TABLine
        !! FWHM = 2*sqrt(2*log(2))*sigma ~ 2.355 (suppose PSF is gaussian)
        !! R = lambda / FWHM
        !! => sigma = lambda / R / 2.355
        Cline = TABLine(itab)%wave
        Wline = degradeRes(Cline, .01_DP, 'NG-SL-LL')
        
        !! Cline
        sig = Wline / 2.355_DP
        IF (.NOT. parinfo(ind%Cline(i))%limited(1)) &
          parinfo(ind%Cline(i))%limits(1) = Cline-sig
        IF (.NOT. parinfo(ind%Cline(i))%limited(2)) &
          parinfo(ind%Cline(i))%limits(2) = Cline+sig
        parinfo(ind%Cline(i))%limited = .TRUE.
        
        !! Wline
        IF (.NOT. parinfo(ind%Wline(i))%limited(1)) &
          parinfo(ind%Wline(i))%limits(1) = Wline*0.5_DP
        IF (.NOT. parinfo(ind%Wline(i))%limited(2)) &
          parinfo(ind%Wline(i))%limits(2) = Wline*2._DP
        parinfo(ind%Wline(i))%limited = .TRUE.

        !! lnRline
        IF (.NOT. parinfo(ind%lnRline(i))%limited(1)) &
          parinfo(ind%lnRline(i))%limits(1) = -10._DP
        IF (.NOT. parinfo(ind%lnRline(i))%limited(2)) &
          parinfo(ind%lnRline(i))%limits(2) = 10._DP
        parinfo(ind%lnRline(i))%limited = .TRUE.

        !! Check wavelength mask to decide if mask par
        DO x=1,Nx
          DO y=1,Ny
            IF (maskxy(x,y)) THEN
              !! Update param
              Cline = iniparval(x,y,ind%Cline(i))
              Wline = iniparval(x,y,ind%Wline(i))
              limi = Cline - Wline
              lims = Cline + Wline
              IF ( ALL( wOBS(iwFIT(x,y)%arr1d(:)).LE.limi &
                        .OR. wOBS(iwFIT(x,y)%arr1d(:)).GE.lims) ) THEN
                mask(x,y,ind%lnRband(i)) = .FALSE.
                mask(x,y,ind%Cband(i)) = .FALSE.
                mask(x,y,ind%WSband(i)) = .FALSE.
                mask(x,y,ind%WLband(i)) = .FALSE.
                !! Update maskxy
                maskxy(x,y) = ANY(mask(x,y,:))
              END IF
            END IF
          END DO
        END DO
        
      END DO
      
      DO i=1,Ncont
        !! lnT LOG[K]
        limi = MERGE(parinfo(ind%lnT(i))%limits(1), &
                     LOG(50._DP), &
                     parinfo(ind%lnT(i))%limited(1))
        lims = MERGE(parinfo(ind%lnT(i))%limits(2), &
                     LOG(500._DP), &
                     parinfo(ind%lnT(i))%limited(2))
        IF (ind%ordQ(i)>0) THEN
          !! lnT(i) becomes lndT(i) = LOG( EXP(lnT(i)) - EXP(lnT(i-1)) )
          IF (.NOT. parinfo(ind%lnT(i))%limited(1)) &
            parinfo(ind%lnT(i))%limits(1) = 0._DP
          IF (.NOT. parinfo(ind%lnT(i))%limited(2)) &
            parinfo(ind%lnT(i))%limits(2) = LOG( EXP(lims)-EXP(limi) )
          parinfo(ind%lnT(i))%limited = .TRUE.
        ELSE
          IF (.NOT. parinfo(ind%lnT(i))%limited(1)) &
            parinfo(ind%lnT(i))%limits(1) = limi
          IF (.NOT. parinfo(ind%lnT(i))%limited(2)) &
            parinfo(ind%lnT(i))%limits(2) = lims
          parinfo(ind%lnT(i))%limited = .TRUE.
        END IF
      END DO

      DO i=1,Nextc
        !! lnAv LOG[mag]
        IF (.NOT. parinfo(ind%lnAv(i))%limited(1)) &
          parinfo(ind%lnAv(i))%limits(1) = 0._DP
        IF (.NOT. parinfo(ind%lnAv(i))%limited(2)) &
          parinfo(ind%lnAv(i))%limits(2) = 5._DP
        parinfo(ind%lnAv(i))%limited = .TRUE.
      END DO
      
      !! Replace the automatic initial values by those
      DO i=1,Npar
        IF ( ANY(TRIMEQ(iniparname(:),parinfo(i)%name)) ) THEN
          CALL IWHERE(TRIMEQ(iniparname(:),parinfo(i)%name), ipar)
          WHERE ( maskxy(:,:) ) par(:,:,i,1) = iniparval(:,:,ipar)
          
        END IF
      END DO

      !! Free memory space
      DEALLOCATE (iniparname,iniparval)
    
    END IF newinipar

    !! Check that the initial guesses satisfy the parameter constraints
    FORALL (i=1:Npar,j=1:SIZE(par,4))
      par(:,:,i,j) = MERGE(parinfo(i)%value,par(:,:,i,j),parinfo(i)%fixed)
      WHERE ( parinfo(i)%limited(1) .AND. (par(:,:,i,j) < parinfo(i)%limits(1)) &
             .AND. maskxy(:,:) )
        par(:,:,i,j) = parinfo(i)%limits(1)
      ELSEWHERE (parinfo(i)%limited(2) &
                 .AND. (par(:,:,i,j) > parinfo(i)%limits(2)) &
                 .AND. maskxy(:,:)) 
        par(:,:,i,j) = parinfo(i)%limits(2)
      END WHERE
    END FORALL

    !! Check there is no NaN
    IF (ANY(isNaN(par(:,:,:,:)))) THEN
      Nmc = MERGE(NiniMC,1,NiniMC > 0)
      ALLOCATE (mu(Npar,Nmc))
      FORALL (i=1:Npar,j=1:Nmc) &
        mu(i,j) = MEAN(par(:,:,i,j),MASK=(maskxy(:,:) &
                                          .AND. .NOT. isNaN(par(:,:,i,j))))
      FORALL (i=1:Npar,j=1:SIZE(par,4))
        WHERE ( isNaN(par(:,:,i,j)) .AND. maskxy(:,:) ) par(:,:,i,j) = mu(i,j)
      END FORALL
    END IF
    
    !! Enforce tying of parameters
    FORALL (i=1:Npar,itied(i) > 0) par(:,:,itied(i),:) = par(:,:,i,:)
    
  END SUBROUTINE initparam

  !!-------------------------------------------------------
  !!
  !!             Read the input analysis settings
  !!
  !!-------------------------------------------------------
  SUBROUTINE read_analysis( dirIN, dirOUT, ACF, verbose, t_burnin, t_end )

    USE utilities, ONLY: DP, warning, trimeq, trimlr, pring!, verbatim
    USE inout, ONLY: lenpar, lenpath, read_input_line, read_hdf5, h5ext
    USE arrays, ONLY: reallocate, iwhere
    USE interpolation, ONLY: interp_lin_sorted
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN), OPTIONAL :: dirIN ! default: ./
 
    CHARACTER(lenpath) :: filanal
    CHARACTER(lenpath), DIMENSION(:), ALLOCATABLE :: path1d
    INTEGER, DIMENSION(:), ALLOCATABLE :: intarr1d
    CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: strarr1d
    
    CHARACTER(lenpath), INTENT(OUT), OPTIONAL :: dirOUT
    INTEGER, DIMENSION(2), INTENT(OUT), OPTIONAL :: t_burnin, t_end
    LOGICAL, INTENT(OUT), OPTIONAL :: verbose
    LOGICAL, DIMENSION(2), INTENT(OUT), OPTIONAL :: ACF

    !! Read the input master file
    !!----------------------------
    IF (PRESENT(dirIN)) THEN
      filanal = TRIMLR(dirIN)//'input_analysis'//h5ext
    ELSE
      CALL WARNING('READ_ANALYSIS','Default input path')
      filanal = './input_analysis'//h5ext
    END IF

    !! Output path
    IF (PRESENT(dirOUT)) THEN
      CALL READ_HDF5(STRARR1D=path1d, FILE=filanal, NAME='output dir')
      dirOUT = path1d(1)
    END IF

    !! input_analysis
    IF (PRESENT(ACF)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filanal, NAME='ACF')
      ACF = [ TRIMEQ(strarr1d(1),'T'), TRIMEQ(strarr1d(2),'T') ]
    END IF
    IF (PRESENT(verbose)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filanal, NAME='verbose')
      verbose = TRIMEQ(strarr1d(1),'T')
    END IF
    IF (PRESENT(t_burnin)) THEN
      CALL READ_HDF5(INTARR1D=intarr1d, FILE=filanal, NAME='t_burnin')
      t_burnin = intarr1d
    END IF
    IF (PRESENT(t_end)) THEN
      CALL READ_HDF5(INTARR1D=intarr1d, FILE=filanal, NAME='t_end')
      t_end = intarr1d
    END IF

  END SUBROUTINE read_analysis
  
  !!-------------------------------------------------------
  !!
  !!     Sherman-Morrison approach to invert matrices
  !!
  !!-------------------------------------------------------
  SUBROUTINE check_SM( A, A_prev, pos, delta )
    
    USE utilities, ONLY: DP, strike
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(IN) :: A, A_prev
    INTEGER :: x, y, Nx, Ny
    LOGICAL :: flag
    REAL(DP), INTENT(OUT), OPTIONAL :: delta
    INTEGER, DIMENSION(2), INTENT(OUT), OPTIONAL :: pos
    
    Nx = SIZE(A,1)
    Ny = SIZE(A,2)
    flag = .FALSE.
    
    !! Find the mole
    DO x=1,Nx
      DO y=1,Ny
        IF (A(x,y)/=A_prev(x,y)) THEN
          IF (.NOT. flag) THEN
            flag = .TRUE.
            IF (PRESENT(pos)) pos(:) = [x,y]
            IF (PRESENT(delta)) delta = A(x,y) - A_prev(x,y)
          ELSE
            CALL STRIKE('invert_SM','Cannot use Sherman-Morrison approach.')
          END IF
        END IF
      END DO
    END DO
    
  END SUBROUTINE check_SM
  
  FUNCTION invert_SM( A, invA, pos, delta, Ndim, &
                      detA, determinant )

    USE utilities, ONLY: DP, strike
    IMPLICIT NONE

    INTEGER, DIMENSION(2), INTENT(IN) :: pos
    REAL(DP), INTENT(IN) :: delta
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: A, invA
    INTEGER, INTENT(IN), OPTIONAL :: Ndim
    REAL(DP), INTENT(IN), OPTIONAL :: detA
    INTEGER :: y, Ny
    REAL(DP) :: c0, c1
    REAL(DP), INTENT(OUT), OPTIONAL :: determinant
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: invert_SM

    IF (PRESENT(Ndim)) THEN
      Ny = Ndim
    ELSE
      Ny = SIZE(A,2)
    END IF

    !! Calculate invert matrix
    invert_SM(:,:) = 0._DP
    c0 = 1 + invA(pos(2),pos(1)) * delta

    IF (c0/=0._DP) THEN
      invert_SM(pos(2),:) = invA(pos(2),:) / c0
      invert_SM(:,pos(1)) = invA(:,pos(1)) / c0

      c1 = delta / c0
      FORALL (y=1:Ny,y/=pos(1)) &
        invert_SM(:,y) = invA(:,y) - invA(:,pos(1))*invA(pos(2),y)*c1

    END IF    

    IF (PRESENT(determinant)) THEN
      IF (PRESENT(detA)) THEN
        determinant = detA * (1+invA(pos(2),pos(1))*delta)
      ELSE
        CALL STRIKE('INVERT_SM','Input detA!')
      END IF
    END IF

  END FUNCTION invert_SM

  !!-------------------------------------------------------
  !!
  !! Modified Sherman-Morrison approach to invert matrices
  !!
  !!-------------------------------------------------------
  FUNCTION invert_mSM( A, invA, pos, delta, Ndim, &
                       detA, determinant, noposdef )

    USE utilities, ONLY: DP, isNaN, strike, warning
    ! USE matrices, ONLY: determinant_matrix
    IMPLICIT NONE
    
    INTEGER, DIMENSION(2), INTENT(IN) :: pos
    REAL(DP), INTENT(IN) :: delta
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: A, invA
    INTEGER, INTENT(IN), OPTIONAL :: Ndim
    REAL(DP), INTENT(IN), OPTIONAL :: detA
    INTEGER :: N
    INTEGER, DIMENSION(2) :: pos_trans
    REAL(DP) :: detA1, detA2
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: A1, invA1!, A2
    REAL(DP), INTENT(OUT), OPTIONAL :: determinant
    LOGICAL, INTENT(OUT), OPTIONAL :: noposdef
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: invert_mSM

    IF (PRESENT(Ndim)) THEN
      N = Ndim
    ELSE
      N = SIZE(A,2)
    END IF
    
    IF (PRESENT(noposdef)) noposdef = .FALSE.

    detA2 = 0._DP
    
    !! Create intermediate matrix A1
    A1(:,:) = A(:,:)
    A1(pos(1),pos(2)) = A1(pos(1),pos(2)) + delta

    pos_trans(:) = [pos(2), pos(1)]

    !! Calculate invert matrix
    IF (PRESENT(detA)) THEN
      invA1(:,:) = invert_SM(A, invA, pos, delta, N, &
                             detA, detA1)
      invert_mSM(:,:) = invert_SM(A1, invA1, pos_trans, delta, N, &
                                  detA1, detA2)
    ELSE
      invA1(:,:) = invert_SM(A, invA, pos, delta, N)
      invert_mSM(:,:) = invert_SM(A1, invA1, pos_trans, delta, N)
      IF (PRESENT(determinant)) &
        CALL STRIKE('INVERT_MSM','No detA input!')
        ! CALL WARNING('INVERT_MSM','No detA input!')
        ! A2(:,:) = A1(:,:)
        ! A2(pos(2),pos(1)) = A2(pos(2),pos(1)) + delta
        ! detA2 = determinant_matrix(A2)
    END IF
    
    IF (PRESENT(determinant)) determinant = detA2

    IF (PRESENT(noposdef)) &
      noposdef = ANY( isNaN(invert_mSM(:,:)) )

  END FUNCTION invert_mSM

  !!-------------------------------------------------------
  !!
  !!               Convert NaNs to zeros
  !!
  !!-------------------------------------------------------
  FUNCTION NaN2zero_int1D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(IN) :: array
    INTEGER :: i1
    INTEGER, DIMENSION(SIZE(array)) :: NaN2zero_int1D

    NaN2zero_int1D(:) = array(:)
    DO i1=1,SIZE(array)
      IF (isNaN(REAL(array(i1)))) NaN2zero_int1D(i1) = 0
    END DO

  END FUNCTION NaN2zero_int1D

  !!==============================

  FUNCTION NaN2zero_int2D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN) :: array
    INTEGER :: i1, i2
    INTEGER, DIMENSION(SIZE(array,1),SIZE(array,2)) :: NaN2zero_int2D

    NaN2zero_int2D(:,:) = array(:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        IF (isNaN(REAL(array(i1,i2)))) NaN2zero_int2D(i1,i2) = 0
      END DO
    END DO

  END FUNCTION NaN2zero_int2D

  !!==============================

  FUNCTION NaN2zero_int3D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: array
    INTEGER :: i1, i2, i3
    INTEGER, DIMENSION(SIZE(array,1),SIZE(array,2),SIZE(array,3)) :: NaN2zero_int3D

    NaN2zero_int3D(:,:,:) = array(:,:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        DO i3=1,SIZE(array,3)
          IF (isNaN(REAL(array(i1,i2,i3)))) NaN2zero_int3D(i1,i2,i3) = 0
        END DO
      END DO
    END DO

  END FUNCTION NaN2zero_int3D

  !!==============================

  FUNCTION NaN2zero_int4D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    INTEGER, DIMENSION(:,:,:,:), INTENT(IN) :: array
    INTEGER :: i1, i2, i3, i4
    INTEGER, DIMENSION(SIZE(array,1),SIZE(array,2),SIZE(array,3),SIZE(array,4)) :: &
      NaN2zero_int4D

    NaN2zero_int4D(:,:,:,:) = array(:,:,:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        DO i3=1,SIZE(array,3)
          DO i4=1,SIZE(array,4)
            IF (isNaN(REAL(array(i1,i2,i3,i4)))) NaN2zero_int4D(i1,i2,i3,i4) = 0
          END DO
        END DO
      END DO
    END DO

  END FUNCTION NaN2zero_int4D

  !!==============================

  FUNCTION NaN2zero_dbl1D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: array
    INTEGER :: i1
    REAL(DP), DIMENSION(SIZE(array)) :: NaN2zero_dbl1D

    NaN2zero_dbl1D(:) = array(:)
    DO i1=1,SIZE(array)
      IF (isNaN(array(i1))) NaN2zero_dbl1D(i1) = 0._DP
    END DO

  END FUNCTION NaN2zero_dbl1D

  !!==============================

  FUNCTION NaN2zero_dbl2D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(IN) :: array
    INTEGER :: i1, i2
    REAL(DP), DIMENSION(SIZE(array,1),SIZE(array,2)) :: NaN2zero_dbl2D

    NaN2zero_dbl2D(:,:) = array(:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        IF (isNaN(array(i1,i2))) NaN2zero_dbl2D(i1,i2) = 0._DP
      END DO
    END DO

  END FUNCTION NaN2zero_dbl2D

  !!==============================

  FUNCTION NaN2zero_dbl3D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: array
    INTEGER :: i1, i2, i3
    REAL(DP), DIMENSION(SIZE(array,1),SIZE(array,2),SIZE(array,3)) :: NaN2zero_dbl3D

    NaN2zero_dbl3D(:,:,:) = array(:,:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        DO i3=1,SIZE(array,3)
          IF (isNaN(array(i1,i2,i3))) NaN2zero_dbl3D(i1,i2,i3) = 0._DP
        END DO
      END DO
    END DO

  END FUNCTION NaN2zero_dbl3D

  !!==============================

  FUNCTION NaN2zero_dbl4D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: array
    INTEGER :: i1, i2, i3, i4
    REAL(DP), DIMENSION(SIZE(array,1),SIZE(array,2),SIZE(array,3),SIZE(array,4)) :: &
      NaN2zero_dbl4D

    NaN2zero_dbl4D(:,:,:,:) = array(:,:,:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        DO i3=1,SIZE(array,3)
          DO i4=1,SIZE(array,4)
            IF (isNaN(array(i1,i2,i3,i4))) NaN2zero_dbl4D(i1,i2,i3,i4) = 0._DP
          END DO
        END DO
      END DO
    END DO

  END FUNCTION NaN2zero_dbl4D

  !!-------------------------------------------------------
  !!
  !!               Convert zeros to NaNs
  !!
  !!-------------------------------------------------------
  FUNCTION zero2NaN_int1D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(IN) :: array
    INTEGER :: i1
    INTEGER, DIMENSION(SIZE(array)) :: zero2NaN_int1D

    zero2NaN_int1D(:) = array(:)
    DO i1=1,SIZE(array)
      IF (array(i1)==0) zero2NaN_int1D(i1) = INT(NaN())
    END DO

  END FUNCTION zero2NaN_int1D

  !!==============================

  FUNCTION zero2NaN_int2D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN) :: array
    INTEGER :: i1, i2
    INTEGER, DIMENSION(SIZE(array,1),SIZE(array,2)) :: zero2NaN_int2D

    zero2NaN_int2D(:,:) = array(:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        IF (array(i1,i2)==0) zero2NaN_int2D(i1,i2) = INT(NaN())
      END DO
    END DO

  END FUNCTION zero2NaN_int2D

  !!==============================

  FUNCTION zero2NaN_int3D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: array
    INTEGER :: i1, i2, i3
    INTEGER, DIMENSION(SIZE(array,1),SIZE(array,2),SIZE(array,3)) :: zero2NaN_int3D

    zero2NaN_int3D(:,:,:) = array(:,:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        DO i3=1,SIZE(array,3)
          IF (array(i1,i2,i3)==0) zero2NaN_int3D(i1,i2,i3) = INT(NaN())
        END DO
      END DO
    END DO

  END FUNCTION zero2NaN_int3D

  !!==============================

  FUNCTION zero2NaN_int4D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    INTEGER, DIMENSION(:,:,:,:), INTENT(IN) :: array
    INTEGER :: i1, i2, i3, i4
    INTEGER, DIMENSION(SIZE(array,1),SIZE(array,2),SIZE(array,3),SIZE(array,4)) :: &
      zero2NaN_int4D

    zero2NaN_int4D(:,:,:,:) = array(:,:,:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        DO i3=1,SIZE(array,3)
          DO i4=1,SIZE(array,4)
            IF (array(i1,i2,i3,i4)==0) zero2NaN_int4D(i1,i2,i3,i4) = INT(NaN())
          END DO
        END DO
      END DO
    END DO

  END FUNCTION zero2NaN_int4D

  !!==============================

  FUNCTION zero2NaN_dbl1D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: array
    INTEGER :: i1
    REAL(DP), DIMENSION(SIZE(array)) :: zero2NaN_dbl1D

    zero2NaN_dbl1D(:) = array(:)
    DO i1=1,SIZE(array)
      IF (array(i1)==0) zero2NaN_dbl1D(i1) = NaN()
    END DO

  END FUNCTION zero2NaN_dbl1D

  !!==============================

  FUNCTION zero2NaN_dbl2D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(IN) :: array
    INTEGER :: i1, i2
    REAL(DP), DIMENSION(SIZE(array,1),SIZE(array,2)) :: zero2NaN_dbl2D

    zero2NaN_dbl2D(:,:) = array(:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        IF (array(i1,i2)==0) zero2NaN_dbl2D(i1,i2) = NaN()
      END DO
    END DO

  END FUNCTION zero2NaN_dbl2D

  !!==============================

  FUNCTION zero2NaN_dbl3D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: array
    INTEGER :: i1, i2, i3
    REAL(DP), DIMENSION(SIZE(array,1),SIZE(array,2),SIZE(array,3)) :: zero2NaN_dbl3D

    zero2NaN_dbl3D(:,:,:) = array(:,:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        DO i3=1,SIZE(array,3)
          IF (array(i1,i2,i3)==0) zero2NaN_dbl3D(i1,i2,i3) = NaN()
        END DO
      END DO
    END DO

  END FUNCTION zero2NaN_dbl3D

  !!==============================

  FUNCTION zero2NaN_dbl4D( array )

    USE utilities, ONLY: DP, isNaN, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: array
    INTEGER :: i1, i2, i3, i4
    REAL(DP), DIMENSION(SIZE(array,1),SIZE(array,2),SIZE(array,3),SIZE(array,4)) :: &
      zero2NaN_dbl4D

    zero2NaN_dbl4D(:,:,:,:) = array(:,:,:,:)
    DO i1=1,SIZE(array,1)
      DO i2=1,SIZE(array,2)
        DO i3=1,SIZE(array,3)
          DO i4=1,SIZE(array,4)
            IF (array(i1,i2,i3,i4)==0) zero2NaN_dbl4D(i1,i2,i3,i4) = NaN()
          END DO
        END DO
      END DO
    END DO

  END FUNCTION zero2NaN_dbl4D
  
  !!-------------------------------------------------------
  !!
  !! Automatize the degradation of the spectral resolution
  !!
  !!-------------------------------------------------------
  FUNCTION degradeRes( wc, dw, instr )

    USE auxil, ONLY: res
    USE utilities, ONLY: DP
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: instr
    REAL(DP), INTENT(IN) :: wc, dw
    REAL(DP) :: degradeRes
    
    !! Line width fixed by the instrumental resolution
    SELECT CASE(instr)
      CASE('SH')
        degradeRes = res%dw_w_SH * wc
      CASE('LH')
        degradeRes = res%dw_w_LH * wc
      CASE('SL')
        degradeRes = res%dw_w_SL * wc
      CASE('LL')
        degradeRes = res%dw_w_LL * wc
      CASE('SL-SH')
        IF (wc.LT.10._DP) THEN
          degradeRes = res%dw_w_SL * wc
        ELSE
          degradeRes = res%dw_w_SH * wc
        END IF
      CASE('SL-LL')
        IF (wc.LT.14.29_DP) THEN
          degradeRes = res%dw_w_SL * wc
        ELSE
          degradeRes = res%dw_w_LL * wc
        END IF
      CASE('NG-SL-LL')
        IF (wc.LT.5.15_DP) THEN
          degradeRes = res%dw_w_AKARI_NG * wc
        ELSE IF (wc.LT.14.29_DP) THEN
          degradeRes = res%dw_w_SL * wc
        ELSE
          degradeRes = res%dw_w_LL * wc
        END IF
      CASE('CAM')
        degradeRes = res%dw_w_CAM * wc
      CASE('SWS')
        degradeRes = res%dw_w_SWS * wc
      CASE('SWSfine')
        degradeRes = res%dw_w_SWSfine * wc
      CASE DEFAULT
        degradeRes = dw

    END SELECT

  END FUNCTION degradeRes

  !!-------------------------------------------------------
  !!
  !!    Analytical functions of the individual features
  !!
  !!-------------------------------------------------------
  
  !!-----------------------
  !! Dust contimuum (MBB)
  !!-----------------------
  PURE FUNCTION modifBB_gen( x, temp, Qabs, wi, inorm )
    !! temp [K]; BLACKBODY [W/m2/sr/Hz]; 
    !! output array [W/sr/Hz/kg] or [Hz-1] if inorm presents
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: temp
    TYPE(Qabs_type), INTENT(IN), DIMENSION(:) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL :: wi ! input grid is in wavelength
    INTEGER, INTENT(IN), OPTIONAL :: inorm
    INTEGER :: Ntemp, itemp
    REAL(DP), DIMENSION(SIZE(x)) :: nu
    REAL(DP), DIMENSION(SIZE(temp),SIZE(x)) :: modifBB_gen

    Ntemp = SIZE(temp)

    !! nu should be the same with Qabs%nua
    IF (PRESENT(wi)) THEN
      IF (wi) nu(:) = MKS%clight/MKS%micron / x(:)
    ELSE
      nu(:) = x(:)
    END IF
    
    !! modifBB [W/sr/Hz/kg] = kappa [m2/kg] * BLACKBODY [W/m2/sr/Hz]
    FORALL (itemp=1:Ntemp) &
      modifBB_gen(itemp,:) = Qabs(itemp)%kappa(:)*BLACKBODY(nu(:),temp(itemp))

    IF (PRESENT(inorm)) THEN
      !! Normalised at nu(inorm)
      FORALL (itemp=1:Ntemp) &
        modifBB_gen(itemp,:) = modifBB_gen(itemp,:) / modifBB_gen(itemp,inorm)

    END IF
  
  END FUNCTION modifBB_gen

  PURE FUNCTION modifBB_0( x, temp, Qabs, wi, inorm )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: temp
    TYPE(Qabs_type), INTENT(IN) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL :: wi
    INTEGER, INTENT(IN), OPTIONAL :: inorm
    REAL(DP), DIMENSION(SIZE(x)) :: modifBB_0
    IF (PRESENT(wi)) THEN
      IF (PRESENT(inorm)) THEN
        modifBB_0(:) = RESHAPE( modifBB_gen(x(:),[temp],[Qabs],wi,inorm), &
                                [SIZE(x(:))] )
      ELSE
        modifBB_0(:) = RESHAPE(modifBB_gen(x(:),[temp],[Qabs],wi), [SIZE(x(:))])
      END IF
    ELSE
      IF (PRESENT(inorm)) THEN
        modifBB_0(:) = RESHAPE( modifBB_gen(x(:),[temp],[Qabs],INORM=inorm), &
                                [SIZE(x(:))] )
      ELSE
        modifBB_0(:) = RESHAPE(modifBB_gen(x(:),[temp],[Qabs]), [SIZE(x(:))])
      END IF
    END IF
  END FUNCTION modifBB_0

  PURE FUNCTION modifBB_1( x, temp, Qabs, wi, inorm )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: temp
    TYPE(Qabs_type), INTENT(IN), DIMENSION(:) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL :: wi
    INTEGER, INTENT(IN), OPTIONAL :: inorm
    REAL(DP), DIMENSION(SIZE(temp),SIZE(x)) :: modifBB_1
    IF (PRESENT(wi)) THEN
      IF (PRESENT(inorm)) THEN
        modifBB_1(:,:) = RESHAPE( modifBB_gen(x(:),temp(:),Qabs(:),wi,inorm), &
                                  [SIZE(temp(:)),SIZE(x(:))] )
      ELSE
        modifBB_1(:,:) = RESHAPE( modifBB_gen(x(:),temp(:),Qabs(:),wi), &
                                  [SIZE(temp(:)),SIZE(x(:))] )
      END IF
    ELSE
      IF (PRESENT(inorm)) THEN
        modifBB_1(:,:) = RESHAPE( modifBB_gen(x(:),temp(:),Qabs(:),INORM=inorm), &
                                  [SIZE(temp(:)),SIZE(x(:))] )
      ELSE
        modifBB_1(:,:) = RESHAPE( modifBB_gen(x(:),temp(:),Qabs(:)), &
                                  [SIZE(temp(:)),SIZE(x(:))] )
      END IF
    END IF
  END FUNCTION modifBB_1

  !!-----------------------------------------------------
  !! Atomic & molecular unresolved lines (Gauss profile)
  !!-----------------------------------------------------
  PURE FUNCTION gaussLine_gen( x, mu, sig, wi )
    !! Output array [Hz^-1]
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE distributions, ONLY: dist_gauss
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: mu, sig
    LOGICAL, INTENT(IN), OPTIONAL :: wi
    INTEGER :: N, i
    REAL(DP), DIMENSION(SIZE(x)) :: nu
    REAL(DP), DIMENSION(SIZE(mu)) :: numu, nusig
    REAL(DP), DIMENSION(SIZE(mu),SIZE(x)) :: gaussLine_gen

    N = SIZE(mu)

    nu(:) = x(:)
    numu(:) = mu(:)
    nusig(:) = sig(:)
    IF (PRESENT(wi)) THEN
      IF (wi) THEN
        nu(:) = MKS%clight/MKS%micron / x(:)
        numu(:) = MKS%clight/MKS%micron / mu(:)
        nusig(:) = MKS%clight/MKS%micron * &
                   (1./(mu(:)-sig(:)) - 1./(mu(:)+sig(:))) / 2._DP
      
      END IF
    END IF

    FORALL (i=1:N) &
      gaussLine_gen(i,:) = DIST_GAUSS(nu(:), numu(i), nusig(i))
    
  END FUNCTION gaussLine_gen

  PURE FUNCTION gaussLine_0( x, mu, sig, wi )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: mu, sig
    LOGICAL, INTENT(IN), OPTIONAL :: wi
    REAL(DP), DIMENSION(SIZE(x)) :: gaussLine_0
    IF (PRESENT(wi)) THEN
      gaussLine_0(:) = RESHAPE(gaussLine_gen(x(:),[mu],[sig],wi), [SIZE(x(:))])
    ELSE
      gaussLine_0(:) = RESHAPE(gaussLine_gen(x(:),[mu],[sig]), [SIZE(x(:))])
    END IF
  END FUNCTION gaussLine_0

  PURE FUNCTION gaussLine_1( x, mu, sig, wi )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: mu, sig
    LOGICAL, INTENT(IN), OPTIONAL :: wi
    REAL(DP), DIMENSION(SIZE(mu),SIZE(x)) :: gaussLine_1
    IF (PRESENT(wi)) THEN
      gaussLine_1(:,:) = RESHAPE(gaussLine_gen(x(:),mu(:),sig(:),wi), &
                                 [SIZE(mu(:)),SIZE(x(:))])
    ELSE
      gaussLine_1(:,:) = RESHAPE(gaussLine_gen(x(:),mu(:),sig(:)), &
                                 [SIZE(mu(:)),SIZE(x(:))])
    END IF
  END FUNCTION gaussLine_1
  
  !!------------------------------------------------------
  !! Resolved aromatic bands (Asymmetric Lorentz profile)
  !!------------------------------------------------------
  PURE FUNCTION lorentzBand_gen( x, mu, sigS, sigL, wi )
    !! Output array [Hz^-1]
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE distributions, ONLY: dist_splitlorentz
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: mu, sigS, sigL
    LOGICAL, INTENT(IN), OPTIONAL :: wi
    INTEGER :: N, i
    REAL(DP), DIMENSION(SIZE(x)) :: nu
    REAL(DP), DIMENSION(SIZE(mu)) :: numu, nusigS, nusigL
    REAL(DP), DIMENSION(SIZE(mu)) :: lambda, tau
    REAL(DP), DIMENSION(SIZE(mu),SIZE(x)) :: lorentzBand_gen

    N = SIZE(mu)

    nu(:) = x(:)
    numu(:) = mu(:)
    nusigS(:) = sigS(:)
    nusigL(:) = sigL(:)
    IF (PRESENT(wi)) THEN
      IF (wi) THEN
        nu(:) = MKS%clight/MKS%micron / x(:)
        numu(:) = MKS%clight/MKS%micron / mu(:)
        nusigS(:) = MKS%clight/MKS%micron * (1./(mu(:)-sigS(:)) - 1./mu(:))
        nusigL(:) = MKS%clight/MKS%micron * (1./mu(:) - 1./(mu(:)+sigL(:)))
        
      END IF
    END IF
    
    !! The SwING library defines the split Lorentzian distribution with
    !! the parameterization (lambda, tau) instead of (sigS, sigL).
    !! In order to make the larger/smaller wavelength/frequency width
    !! be at the longer wavelength side, we define (in frequency grid):
    !! lambda = 2sigL & tau = sigS/sigL,
    !! where S corresponds to short wavelength side and L long wavelength
    !! Note that mu is not the mean and that neither sigS nor sigL is
    !! the std_dev of a Lorentzian, 'cause it is undefined (infinite).
    !! Indeed, lorentzian has no moment.
    !! This is an analogue of (split) normal distribution.
    !! In addition, we denote lambda=2*delta_nu_long (analogue of FWHM)
    !! and tau=delta_wave_long/delta_wave_short.
    lambda(:) = 2*nusigL(:) ! lambda -> width param
    tau(:) = nusigS(:)/nusigL(:) ! tau -> shape param
    !! norm = lambda / (1+tau) / pi
    FORALL (i=1:N) &
      lorentzBand_gen(i,:) = DIST_SPLITLORENTZ(nu(:), numu(i), lambda(i), tau(i))

  END FUNCTION lorentzBand_gen
  
  PURE FUNCTION lorentzBand_0( x, mu, sigS, sigL, wi )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: mu, sigS, sigL
    LOGICAL, INTENT(IN), OPTIONAL :: wi
    REAL(DP), DIMENSION(SIZE(x)) :: lorentzBand_0
    IF (PRESENT(wi)) THEN
      lorentzBand_0(:) = RESHAPE(lorentzBand_gen(x(:),[mu],[sigS],[sigL],wi), &
                                 [SIZE(x(:))])
    ELSE
      lorentzBand_0(:) = RESHAPE(lorentzBand_gen(x(:),[mu],[sigS],[sigL]), &
                                 [SIZE(x(:))])
    END IF
  END FUNCTION lorentzBand_0

  PURE FUNCTION lorentzBand_1( x, mu, sigS, sigL, wi )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: mu, sigS, sigL
    LOGICAL, INTENT(IN), OPTIONAL :: wi
    REAL(DP), DIMENSION(SIZE(mu),SIZE(x)) :: lorentzBand_1
    IF (PRESENT(wi)) THEN
      lorentzBand_1(:,:) = RESHAPE(lorentzBand_gen(x(:),mu(:),sigS(:),sigL(:),wi), &
                                   [SIZE(mu(:)), SIZE(x(:))])
    ELSE
      lorentzBand_1(:,:) = RESHAPE(lorentzBand_gen(x(:),mu(:),sigS(:),sigL(:)), &
                                   [SIZE(mu(:)), SIZE(x(:))])
    END IF
  END FUNCTION lorentzBand_1

  !!------------------
  !! Extinction curve
  !!------------------
  FUNCTION extCurve( curve, x, wi )
    !! https://www.astro.princeton.edu/~draine/dust/dustmix.html
    !! Diffuse (MW) Rv = 3.1
    !! nu2w=.TRUE. if input is in frequency [Hz]
    USE utilities, ONLY: DP, trimlr, trimeq, warning
    USE constants, ONLY: MKS
    USE arrays, ONLY: closest
    USE inout, ONLY: lenpath, read_hdf5, h5ext
    USE interpolation, ONLY: interp_lin_sorted
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: curve
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: wi
    CHARACTER(lenpath) :: dir, fil
    REAL(DP) :: wave_V, Cext_V
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave0, Cext0
    REAL(DP), DIMENSION(SIZE(x)) :: wave
    REAL(DP), DIMENSION(SIZE(x)) :: extCurve

    dir = '../lib/'
    IF (TRIMEQ(curve,'D03')) THEN
      fil = TRIMLR(dir)//'Adust_Draine03'//h5ext
    ELSE IF (TRIMEQ(curve,'H2O')) THEN
      fil = TRIMLR(dir)//'Aice_H2O'//h5ext
    ELSE IF (TRIMEQ(curve,'CO')) THEN
      fil = TRIMLR(dir)//'Aice_CO'//h5ext
    ELSE IF (TRIMEQ(curve,'CO2')) THEN
      fil = TRIMLR(dir)//'Aice_CO2'//h5ext
    ELSE
      fil = 'noextc'
      CALL WARNING('EXTCURVE','Unknown component! Av set to 0')
    END IF
    
    IF (PRESENT(wi)) THEN
      IF (wi) wave(:) = x(:)
    ELSE
      wave(:) = MKS%clight/MKS%micron / x(:)
    END IF

    IF (TRIMEQ(fil,'noextc')) THEN
      extCurve(:) = 0._DP
    ELSE
      CALL READ_HDF5(DBLARR1D=wave0, NAME='lambda (micron)', FILE=fil)
      CALL READ_HDF5(DBLARR1D=Cext0, NAME='C_extovH (cm^2ovH)', FILE=fil)
      
      wave_V = 0.5470_DP
      Cext_V = Cext0( CLOSEST(wave0, wave_V) )
      !! interpolated with respect to wavelength
      extCurve(:) = INTERP_LIN_SORTED(Cext0/Cext_V, wave0, wave, &
                                      XLOG=.TRUE., YLOG=.TRUE., FORCE=.TRUE.)
    END IF
      
  END FUNCTION extCurve

  !!---------------------------------------
  !! Total model function for chi2 calling
  !!---------------------------------------

  !! 3D version
  !!============
  FUNCTION specModel_3D( wvl, parval, indpar, mask, Qabs, extinct, verbose, &
                         FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, maskpar, &
                         FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab )

    USE utilities, ONLY: DP, trimeq, trimlr, pring, strike, NaN, &!verbatim, &
                         initiate_clock, time_type, timinfo!, ustd
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: reallocate, iwhere, incrarr
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: wvl
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: extinct
    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: parval ! 3D (Nx,Ny,*Npar)
    TYPE(indpar_type), INTENT(IN) :: indpar
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL :: verbose
    LOGICAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: mask ! 3D (Nx,Ny,Nw)

    INTEGER :: Nw, Nx, Ny, Npar, Ncont, Nline, Nband, Nextc, Nstar
    INTEGER :: x, y, i, k, indref
    INTEGER, DIMENSION(:), ALLOCATABLE :: indw
    REAL(DP) :: Cband, WSband, WLband, Cline, Wline
    REAL(DP) :: Tstar
    REAL(DP), DIMENSION(SIZE(wvl)) :: nu
    REAL(DP), DIMENSION(SIZE(parval,1),SIZE(parval,2)) :: lnIref, lnT
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: dblarr4d ! (Nx,Ny,Nw,Ncomp)
    REAL(DP), DIMENSION(SIZE(parval,1),SIZE(parval,2),SIZE(wvl)) :: &
      FnuCONT0, FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: maskxy
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: mask3D

    LOGICAL, DIMENSION(SIZE(parval,1),SIZE(parval,2),SIZE(parval,3)), INTENT(OUT), OPTIONAL :: &
      maskpar
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
    REAL(DP), DIMENSION(SIZE(parval,1),SIZE(parval,2),SIZE(wvl)) :: specModel_3D

    !! Clock module
    LOGICAL :: printimer
    TYPE(time_type) :: timestr

    CALL INITIATE_CLOCK(timestr)
    
    !! Preliminaries
    !!---------------
    Nw = SIZE(wvl(:))
    Nx = SIZE(parval(:,:,:),1)
    Ny = SIZE(parval(:,:,:),2)
    Npar = SIZE(parval(:,:,:),3)
    Ncont = SIZE(Qabs(:))
    Nband = COUNT(indpar%lnRband(:) .NE. -1)
    Nline = COUNT(indpar%lnRline(:) .NE. -1)
    Nextc = COUNT(indpar%lnAv(:) .NE. -1)
    Nstar = COUNT(indpar%lnFstar(:) .NE. -1)
    Tstar = 5.E4_DP ! [K] high enough to stay in Rayleigh-Jeans limit (T>>5000K for 3um)
    nu(:) = MKS%clight / wvl(:) /MKS%micron ! [Hz] in BLACKBODY

    !! labB index of the reference band (default: 1)
    indref = indpar%refB
    !! (LOG) Intensity of the reference band
    lnIref(:,:) = parval(:,:,indpar%lnRband(indref))

    !! Mask
    ALLOCATE(maskxy(Nx,Ny),mask3D(Nx,Ny,Nw))
    IF (PRESENT(mask)) THEN
      mask3D(:,:,:) = mask(:,:,:)
    ELSE
      mask3D(:,:,:) = .TRUE.
    END IF
    FORALL (x=1:Nx,y=1:Ny) maskxy(x,y) = ANY(mask3D(x,y,:))

    !! mask -> maskpar (see also initparam)
    paramask: IF (PRESENT(maskpar)) THEN
      FORALL (x=1:Nx,y=1:Ny) maskpar(x,y,:) = maskxy(x,y)
      
      IF (PRESENT(mask)) THEN
        !! Find (indices of) all masked wavelengths
        DO x=1,Nx
          DO y=1,Ny
            IF (maskxy(x,y)) THEN
              !! indw is NOT-NaN index list
              CALL IWHERE( mask3D(x,y,:), indw )
              !! Check wavelength mask to decide if mask par
              
              !! Bands
              DO i=1,Nband
                Cband = parval(x,y,indpar%Cband(i))
                WSband = parval(x,y,indpar%WSband(i))
                WLband = parval(x,y,indpar%WLband(i))
                IF ( ALL( wvl(indw(:)).LE.Cband-WSband &
                          .OR. wvl(indw(:)).GE.Cband+WLband ) ) THEN
                  maskpar(x,y,indpar%lnRband(i)) = .FALSE.
                  maskpar(x,y,indpar%Cband(i)) = .FALSE.
                  maskpar(x,y,indpar%WSband(i)) = .FALSE.
                  maskpar(x,y,indpar%WLband(i)) = .FALSE.
                END IF
              END DO
              
              !! Lines
              DO i=1,Nline
                Cline = parval(x,y,indpar%Cline(i))
                Wline = parval(x,y,indpar%Wline(i))
                Wline = parval(x,y,indpar%Wline(i))
                IF ( ALL( wvl(indw(:)).LE.Cline-Wline &
                          .OR. wvl(indw(:)).GE.Cline+Wline ) ) THEN
                  maskpar(x,y,indpar%lnRline(i)) = .FALSE.
                  maskpar(x,y,indpar%Cline(i)) = .FALSE.
                  maskpar(x,y,indpar%Wline(i)) = .FALSE.
                END IF
              END DO
              
            END IF
          END DO
        END DO
      END IF
      
    END IF paramask
    
    !! Print timer
    printimer = .FALSE.
    IF (PRESENT(verbose)) printimer = verbose

    !! Initialization
    specModel_3D(:,:,:) = 0._DP
    
    !! 1. Continua
    !!-------------
    CALL REALLOCATE(dblarr4d,Nx,Ny,Nw,Ncont)
    dblarr4d(:,:,:,:) = 0._DP
    !! The actual unit depends on the input spectrum.
    !! Here as an example, suppose the input is in MKS.
    !! FnuCONT [W/m2/sr/Hz] = Fcont [W/m2/sr] * modifBB [Hz-1]
    DO i=1,Ncont
      k = indpar%ordQ(i)
      lnT(:,:) = parval(:,:,indpar%lnT(i-k))
      IF (k>0) &
        lnT(:,:) = LOG( EXP(lnT(:,:)) &
                        +SUM(EXP(parval(:,:,indpar%lnT(i-k+1:i))),DIM=3) )
      FORALL (x=1:Nx,y=1:Ny,maskxy(x,y)) &
        dblarr4d(x,y,:,i) = EXP(parval(x,y,indpar%lnFcont(i))) * &
                            MODIFBB(wvl(:),EXP(lnT(x,y)),Qabs(i),.TRUE.,indpar%refw)
    END DO
    FnuCONT0(:,:,:) = SUM(dblarr4d(:,:,:,:),DIM=4)
    IF (PRESENT(FnuCONT)) FnuCONT = FnuCONT0(:,:,:)
    IF (PRESENT(FnuCONT_tab)) FnuCONT_tab = dblarr4d(:,:,:,:)
    ! IF (PRESENT(FnuCONT)) FnuCONT = zero2NaN(FnuCONT0(:,:,:))
    ! IF (PRESENT(FnuCONT_tab)) FnuCONT_tab = zero2NaN(dblarr4d(:,:,:,:))
    
    !! 2. Bands
    !!----------
    CALL REALLOCATE(dblarr4d,Nx,Ny,Nw,Nband)
    dblarr4d(:,:,:,:) = 0._DP
    !! FnuBAND [W/m2/sr/Hz] = Iband [W/m2/sr] * lorentzBand [Hz-1]
    DO i=1,Nband
      IF (i == indref) THEN
        FORALL (x=1:Nx,y=1:Ny,maskxy(x,y)) &
          dblarr4d(x,y,:,i) = EXP( parval(x,y,indpar%lnRband(i)) ) * &
                              lorentzBand( wvl(:), parval(x,y,indpar%Cband(i)), &
                                parval(x,y,indpar%WSband(i)), &
                                parval(x,y,indpar%WLband(i)),.TRUE. )
      ELSE
        FORALL (x=1:Nx,y=1:Ny,maskxy(x,y)) &
          dblarr4d(x,y,:,i) = EXP( parval(x,y,indpar%lnRband(i))+lnIref(x,y) ) * &
                              lorentzBand( wvl(:), parval(x,y,indpar%Cband(i)), &
                                parval(x,y,indpar%WSband(i)), &
                                parval(x,y,indpar%WLband(i)),.TRUE. )
      END IF
    END DO
    FnuBAND0(:,:,:) = SUM(dblarr4d(:,:,:,:),DIM=4)
    IF (PRESENT(FnuBAND)) FnuBAND = FnuBAND0(:,:,:)
    IF (PRESENT(FnuBAND_tab)) FnuBAND_tab = dblarr4d(:,:,:,:)
    ! IF (PRESENT(FnuBAND)) FnuBAND = zero2NaN(FnuBAND0(:,:,:))
    ! IF (PRESENT(FnuBAND_tab)) FnuBAND_tab = zero2NaN(dblarr4d(:,:,:,:))
    
    !! 3. Stellar Continuum
    !!----------------------
    CALL REALLOCATE(dblarr4d,Nx,Ny,Nw,Nstar)
    dblarr4d(:,:,:,:) = 0._DP
    !! FnuSTAR [W/m2/sr/Hz] = Fstar [W/m2] * BLACKBODY [W/m2/sr/Hz] /
    !!                        stefan [W/m2/K4] / Tstar4 [K4]
    FORALL (x=1:Nx,y=1:Ny,i=1:Nstar,maskxy(x,y)) &
      dblarr4d(x,y,:,i) = EXP(parval(x,y,indpar%lnFstar(i))) * &
                          pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
    FnuSTAR0(:,:,:) = SUM(dblarr4d(:,:,:,:),DIM=4)
    IF (PRESENT(FnuSTAR)) FnuSTAR = FnuSTAR0(:,:,:)
    IF (PRESENT(FnuSTAR_tab)) FnuSTAR_tab = dblarr4d(:,:,:,:)
    ! IF (PRESENT(FnuSTAR)) FnuSTAR = zero2NaN(FnuSTAR0(:,:,:))
    ! IF (PRESENT(FnuSTAR_tab)) FnuSTAR_tab = zero2NaN(dblarr4d(:,:,:,:))
    
    !! 4. Screen extinction
    !!----------------------
    CALL REALLOCATE(dblarr4d,Nx,Ny,Nw,Nextc)
    dblarr4d(:,:,:,:) = 0._DP
    !! Pabs [/] = EXP( Av [mag] * extinction [/] )
    FORALL (x=1:Nx,y=1:Ny,i=1:Nextc,maskxy(x,y)) &
      dblarr4d(x,y,:,i) = -EXP(parval(x,y,indpar%lnAv(i)))/1.086_DP * extinct(i,:)
    Pabs0(:,:,:) = EXP( SUM(dblarr4d(:,:,:,:),DIM=4) )
    IF (PRESENT(Pabs)) Pabs = Pabs0(:,:,:)
    IF (PRESENT(Pabs_tab)) Pabs_tab = EXP( dblarr4d(:,:,:,:) )
    ! IF (PRESENT(Pabs)) Pabs = zero2NaN(Pabs0(:,:,:))
    ! IF (PRESENT(Pabs_tab)) Pabs_tab = EXP( zero2NaN(dblarr4d(:,:,:,:)) )
    
    !! 5. Lines
    !!----------
    CALL REALLOCATE(dblarr4d,Nx,Ny,Nw,Nline)
    dblarr4d(:,:,:,:) = 0._DP
    !! FnuLINE [W/m2/sr/Hz] = Iline [W/m2/sr] * gaussLine [Hz-1]
    DO i=1,Nline
      FORALL (x=1:Nx,y=1:Ny,maskxy(x,y)) &
        dblarr4d(x,y,:,i) = EXP( parval(x,y,indpar%lnRline(i))+lnIref(x,y) ) * &
                            gaussLine(wvl(:), parval(x,y,indpar%Cline(i)), &
                              parval(x,y,indpar%Wline(i)),.TRUE.)
    END DO
    FnuLINE0(:,:,:) = SUM(dblarr4d(:,:,:,:),DIM=4)
    IF (PRESENT(FnuLINE)) FnuLINE = FnuLINE0(:,:,:)
    IF (PRESENT(FnuLINE_tab)) FnuLINE_tab = dblarr4d(:,:,:,:)
    ! IF (PRESENT(FnuLINE)) FnuLINE = zero2NaN(FnuLINE0(:,:,:))
    ! IF (PRESENT(FnuLINE_tab)) FnuLINE_tab = zero2NaN(dblarr4d(:,:,:,:))

    !! Total model
    !!-------------
    specModel_3D(:,:,:) = (FnuCONT0(:,:,:) + FnuBAND0(:,:,:) &
                           + FnuSTAR0(:,:,:) + FnuLINE0(:,:,:)) * Pabs0(:,:,:)
    ! specModel_3D(:,:,:) = zero2NaN(specModel_3D(:,:,:))

    IF (printimer) &
      PRINT*, "[specModel] EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
    
  END FUNCTION specModel_3D

  !! 2D version
  !!============
  FUNCTION specModel_2D( wvl, parval, indpar, mask, Qabs, extinct, verbose, &
                         FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, &
                         FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: wvl
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: parval ! 2D (Nx,*Npar)
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: extinct
    TYPE(indpar_type), INTENT(IN) :: indpar
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL :: verbose
    LOGICAL, DIMENSION(:,:), INTENT(IN), OPTIONAL :: mask ! 2D (Nx,Nw)
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
      FnuCONT0, FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: &
      FnuCONT_tab0, FnuBAND_tab0, FnuSTAR_tab0, Pabs_tab0, FnuLINE_tab0
    LOGICAL :: verbose0
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: mask3D
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
    REAL(DP), DIMENSION(SIZE(parval,1),SIZE(wvl)) :: specModel_2D
    INTEGER :: x, Nw, Nx, Npar
    Nw = SIZE(wvl(:))
    Nx = SIZE(parval(:,:),1)
    Npar = SIZE(parval(:,:),2)
    verbose0 = .FALSE.
    IF (PRESENT(verbose)) verbose0 = verbose
    IF (PRESENT(mask)) THEN
      ALLOCATE(mask3D(Nx,1,Nw))
      FORALL (x=1:Nx) mask3D(x,1,:) = mask(x,:)
      specModel_2D(:,:) = RESHAPE( specModel_3D(wvl(:), INDPAR=indpar, &
                                     PARVAL=RESHAPE(parval(:,:),[Nx,1,Npar]), &
                                     MASK=mask3D(:,:,:), &
                                     QABS=Qabs, EXTINCT=extinct, VERBOSE=verbose0, &
                                     FNUCONT=FnuCONT0, FNUBAND=FnuBAND0, FNUSTAR=FnuSTAR0, &
                                     PABS=Pabs0, FNULINE=FnuLINE0, &
                                     FNUCONT_TAB=FnuCONT_tab0, FNUBAND_TAB=FnuBAND_tab0, &
                                     FNUSTAR_TAB=FnuSTAR_tab0, PABS_TAB=Pabs_tab0, &
                                     FNULINE_TAB=FnuLINE_tab0), &
                                   [Nx,SIZE(wvl(:))] )
    ELSE
      specModel_2D(:,:) = RESHAPE( specModel_3D(wvl(:), INDPAR=indpar, &
                                     PARVAL=RESHAPE(parval(:,:),[Nx,1,Npar]), &
                                     QABS=Qabs, EXTINCT=extinct, VERBOSE=verbose0, &
                                     FNUCONT=FnuCONT0, FNUBAND=FnuBAND0, FNUSTAR=FnuSTAR0, &
                                     PABS=Pabs0, FNULINE=FnuLINE0, &
                                     FNUCONT_TAB=FnuCONT_tab0, FNUBAND_TAB=FnuBAND_tab0, &
                                     FNUSTAR_TAB=FnuSTAR_tab0, PABS_TAB=Pabs_tab0, &
                                     FNULINE_TAB=FnuLINE_tab0), &
                                   [Nx,SIZE(wvl(:))] )
    END IF
    IF (PRESENT(FnuCONT)) FnuCONT = FnuCONT0(:,1,:)
    IF (PRESENT(FnuBAND)) FnuBAND = FnuBAND0(:,1,:)
    IF (PRESENT(FnuSTAR)) FnuSTAR = FnuSTAR0(:,1,:)
    IF (PRESENT(Pabs)) Pabs = Pabs0(:,1,:)
    IF (PRESENT(FnuLINE)) FnuLINE = FnuLINE0(:,1,:)
    IF (PRESENT(FnuCONT_tab)) FnuCONT_tab = FnuCONT_tab0(:,1,:,:)
    IF (PRESENT(FnuBAND_tab)) FnuBAND_tab = FnuBAND_tab0(:,1,:,:)
    IF (PRESENT(FnuSTAR_tab)) FnuSTAR_tab = FnuSTAR_tab0(:,1,:,:)
    IF (PRESENT(Pabs_tab)) Pabs_tab = Pabs_tab0(:,1,:,:)
    IF (PRESENT(FnuLINE_tab)) FnuLINE_tab = FnuLINE_tab0(:,1,:,:)
  END FUNCTION specModel_2D
  
  !! 1D version
  !!============
  FUNCTION specModel_1D( wvl, parval, indpar, mask, Qabs, extinct, verbose, &
                         FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, &
                         FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: wvl
    REAL(DP), DIMENSION(:), INTENT(IN) :: parval ! 1D (*Npar)
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: extinct
    TYPE(indpar_type), INTENT(IN) :: indpar
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL :: verbose
    LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: mask ! 1D (Nw)
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
      FnuCONT0, FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: &
      FnuCONT_tab0, FnuBAND_tab0, FnuSTAR_tab0, Pabs_tab0, FnuLINE_tab0
    LOGICAL :: verbose0
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: mask3D
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
    REAL(DP), DIMENSION(SIZE(wvl)) :: specModel_1D
    INTEGER :: Nw, Npar
    Nw = SIZE(wvl(:))
    Npar = SIZE(parval(:))
    verbose0 = .FALSE.
    IF (PRESENT(verbose)) verbose0 = verbose
    IF (PRESENT(mask)) THEN
      ALLOCATE(mask3D(1,1,Nw))
      mask3D(1,1,:) = mask(:)
      specModel_1D(:) = RESHAPE( specModel_3D(wvl(:), INDPAR=indpar, &
                                   PARVAL=RESHAPE(parval(:),[1,1,Npar]), &
                                   MASK=mask3D(:,:,:), &
                                   QABS=Qabs, EXTINCT=extinct, VERBOSE=verbose0, &
                                   FNUCONT=FnuCONT0, FNUBAND=FnuBAND0, FNUSTAR=FnuSTAR0, &
                                   PABS=Pabs0, FNULINE=FnuLINE0, &
                                   FNUCONT_TAB=FnuCONT_tab0, FNUBAND_TAB=FnuBAND_tab0, &
                                   FNUSTAR_TAB=FnuSTAR_tab0, PABS_TAB=Pabs_tab0, &
                                   FNULINE_TAB=FnuLINE_tab0), &
                                 [SIZE(wvl(:))] )
    ELSE
      specModel_1D(:) = RESHAPE( specModel_3D(wvl(:), INDPAR=indpar, &
                                   PARVAL=RESHAPE(parval(:),[1,1,Npar]), &
                                   QABS=Qabs, EXTINCT=extinct, VERBOSE=verbose0, &
                                   FNUCONT=FnuCONT0, FNUBAND=FnuBAND0, FNUSTAR=FnuSTAR0, &
                                   PABS=Pabs0, FNULINE=FnuLINE0, &
                                   FNUCONT_TAB=FnuCONT_tab0, FNUBAND_TAB=FnuBAND_tab0, &
                                   FNUSTAR_TAB=FnuSTAR_tab0, PABS_TAB=Pabs_tab0, &
                                   FNULINE_TAB=FnuLINE_tab0), &
                                 [SIZE(wvl(:))] )
    END IF
    IF (PRESENT(FnuCONT)) FnuCONT = FnuCONT0(1,1,:)
    IF (PRESENT(FnuBAND)) FnuBAND = FnuBAND0(1,1,:)
    IF (PRESENT(FnuSTAR)) FnuSTAR = FnuSTAR0(1,1,:)
    IF (PRESENT(Pabs)) Pabs = Pabs0(1,1,:)
    IF (PRESENT(FnuLINE)) FnuLINE = FnuLINE0(1,1,:)
    IF (PRESENT(FnuCONT_tab)) FnuCONT_tab = FnuCONT_tab0(1,1,:,:)
    IF (PRESENT(FnuBAND_tab)) FnuBAND_tab = FnuBAND_tab0(1,1,:,:)
    IF (PRESENT(FnuSTAR_tab)) FnuSTAR_tab = FnuSTAR_tab0(1,1,:,:)
    IF (PRESENT(Pabs_tab)) Pabs_tab = Pabs_tab0(1,1,:,:)
    IF (PRESENT(FnuLINE_tab)) FnuLINE_tab = FnuLINE_tab0(1,1,:,:)
  END FUNCTION specModel_1D
  
  !! Generic intertface
  !!====================
  FUNCTION specModel_gen( wvl, parvec, parname, parinfo, parval, &
                          indpar, Qabs, extinct )

    USE utilities, ONLY: DP, trimeq
    USE constants, ONLY: pi, MKS
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: wvl
    REAL(DP), DIMENSION(:), INTENT(IN) :: parvec, parval
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: extinct
    CHARACTER(*), INTENT(IN) :: parname
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo
    TYPE(indpar_type), INTENT(IN) :: indpar
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN) :: Qabs
    
    INTEGER :: Nw, Ngrid, Ncont, Nline, Nband, Nextc, Nstar
    INTEGER :: i, j, k, igrid, indref, iparef
    REAL(DP) :: Tstar, lnT
    REAL(DP), DIMENSION(SIZE(wvl)) :: nu
    LOGICAL :: gridlnFcont, gridlnT, gridlnRline, gridCline, gridWline
    LOGICAL :: gridlnRband, gridCband, gridWSband, gridWLband, gridlnAv, gridlnFstar
    LOGICAL :: gridCONT, gridBAND, gridLINE, gridSTAR, gridEXTC
    REAL(DP), DIMENSION(SIZE(wvl)) :: dblarr1D, spec1D, extc1D
    REAL(DP), DIMENSION(SIZE(parvec),SIZE(wvl)) :: specModel_gen

    
    !! Preliminaries
    !!---------------
    Nw = SIZE(wvl(:))
    Ngrid = SIZE(parvec(:))
    ! Ncont = SIZE(Qabs(:))
    Ncont = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'CONT')) ) / NparCONT
    Nline = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'LINE')) ) / NparLINE
    Nband = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'BAND')) ) / NparBAND
    Nextc = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'EXTC')) ) / NparEXTC
    ! Nstar = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'STAR')) ) / NparSTAR
    Nstar = 1 ! Stellar contamination
    Tstar = 5.E4_DP ! [K] high enough to stay in Rayleigh-Jeans limit (T>>5000K for 3um)
    nu(:) = MKS%clight / wvl(:) /MKS%micron
    
    !! labB index of the reference band (default: 1)
    indref = indpar%refB
    !! par index of the reference band
    iparef = indpar%lnRband(indref)
    
    !! Sampled parameter
    gridlnFcont = ( parname(1:7) == 'lnFcont' )
    gridlnT = ( parname(1:3) == 'lnT' )
    gridlnRband = ( parname(1:7) == 'lnRband' )
    gridCband = ( parname(1:5) == 'Cband' )
    gridWSband = ( parname(1:6) == 'WSband' )
    gridWLband = ( parname(1:6) == 'WLband' )
    gridlnRline = ( parname(1:7) == 'lnRline' )
    gridCline = ( parname(1:5) == 'Cline' )
    gridWline = ( parname(1:5) == 'Wline' )
    gridlnFstar = ( parname(1:7) == 'lnFstar' )
    gridlnAv = ( parname(1:4) == 'lnAv' )
    
    !! Constant components
    gridCONT = ( gridlnFcont .OR. gridlnT )
    gridBAND = ( gridlnRband .OR. gridCband .OR. gridWSband .OR. gridWLband )
    gridLINE = ( gridlnRline .OR. gridCline .OR. gridWline )
    gridSTAR = gridlnFstar
    gridEXTC = gridlnAv
    
    !! Initialization
    specModel_gen(:,:) = 0._DP
    spec1D(:) = 0._DP
    extc1D(:) = 1._DP

    !! 0. Extinction
    !!---------------
    sampEXTC: IF (gridEXTC) THEN
      READ (parname(5:6),*) j
      FORALL (igrid=1:Ngrid) &
        specModel_gen(igrid,:) = EXP( -EXP(parvec(igrid))/1.086_DP*extinct(j,:) )

      DO i=1,Nextc
        IF (i /= j) &
          extc1D(:) = extc1D(:) * &
                      EXP( -EXP(parval(indpar%lnAv(i)))/1.086_DP*extinct(i,:) )
      END DO
    ELSE
      DO i=1,Nextc
        extc1D(:) = extc1D(:) * &
                    EXP( -EXP(parval(indpar%lnAv(i)))/1.086_DP*extinct(i,:) )
      END DO
    END IF sampEXTC
    
    
    !! 1. Continuum
    !!--------------
    sampCONT: IF (gridCONT) THEN
      IF (gridlnFcont) THEN
        READ(parname(8:9),*) j
        k = indpar%ordQ(j)
        lnT = parval(indpar%lnT(j-k))
        IF (k>0) &
          lnT = LOG( EXP(lnT)+SUM(EXP(parval(indpar%lnT(j-k+1:j)))) )
        dblarr1d(:) = MODIFBB( wvl(:),EXP(lnT),Qabs(j),.TRUE.,indpar%refw )
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP(parvec(igrid)) * dblarr1d(:)
      ELSE IF (gridlnT) THEN
        READ(parname(4:5),*) j
        k = indpar%ordQ(j)
        IF (k==0) THEN
          FORALL (igrid=1:Ngrid) &
            specModel_gen(igrid,:) = EXP(parval(indpar%lnFcont(j))) &
              * MODIFBB( wvl(:),EXP(parvec(igrid)),Qabs(j),.TRUE.,indpar%refw )
        ELSE
          lnT = parval(indpar%lnT(j-k))
          IF (k>1) &
            lnT = LOG( EXP(lnT)+SUM(EXP(parval(indpar%lnT(j-k+1:j-1)))) )
          FORALL (igrid=1:Ngrid) &
            specModel_gen(igrid,:) = EXP(parval(indpar%lnFcont(j))) &
              * MODIFBB( wvl(:),EXP(lnT)+EXP(parvec(igrid)), &
                         Qabs(j),.TRUE.,indpar%refw )
        END IF
      END IF
      
      DO i=1,Ncont
        IF (i /= j) THEN
          k = indpar%ordQ(i)
          !! Fais gaffe !!
          IF (gridlnT .AND. j>=i-k .AND. j<i) THEN
            lnT = parval(indpar%lnT(i))
            IF (j==i-k .AND. k>1) THEN
              !! lnTj is lnT
              lnT = LOG( EXP(lnT)+SUM(EXP(parval(indpar%lnT(i-k+1:i-1)))) )
            ELSE IF (j>i-k .AND. j<i-1) THEN ! intrinsicly k>2
              !! lnTj is lndT
              lnT = LOG( EXP(lnT)+SUM(EXP(parval(indpar%lnT(i-k:j-1)))) &
                         +SUM(EXP(parval(indpar%lnT(j+1:i-1)))) )
            ELSE IF (j==i-1 .AND. k>1) THEN
              !! lnTj is lndT
              lnT = LOG( EXP(lnT)+SUM(EXP(parval(indpar%lnT(i-k:i-2)))) )
            END IF
            FORALL (igrid=1:Ngrid) &
              specModel_gen(igrid,:) = specModel_gen(igrid,:) &
                + EXP(parval(indpar%lnFcont(i))) &
                  * MODIFBB( wvl(:),EXP(lnT)+EXP(parvec(igrid)), &
                             Qabs(j),.TRUE.,indpar%refw )
          ELSE
            lnT = parval(indpar%lnT(i-k))
            IF (k>0) &
              lnT = LOG( EXP(lnT)+SUM(EXP(parval(indpar%lnT(i-k+1:i)))) )
            spec1D(:) = spec1D(:) &
                        + EXP(parval(indpar%lnFcont(i))) &
                          * MODIFBB( wvl(:),EXP(lnT),Qabs(i),.TRUE.,indpar%refw )
          END IF
        END IF
      END DO
    ELSE
      DO i=1,Ncont
        k = indpar%ordQ(i)
        lnT = parval(indpar%lnT(i-k))
        IF (k>0) &
          lnT = LOG( EXP(lnT)+SUM(EXP(parval(indpar%lnT(i-k+1:i)))) )
        spec1D(:) = spec1D(:) &
                    + EXP(parval(indpar%lnFcont(i))) &
                      * MODIFBB( wvl(:),EXP(lnT),Qabs(i),.TRUE.,indpar%refw )
      END DO
    END IF sampCONT
    
    
    !! 2. Bands
    !!----------
    sampBAND: IF (gridBAND) THEN
      IF (gridlnRband) THEN
        READ(parname(8:9),*) j
        dblarr1d(:) = LORENTZBAND( wvl(:),parval(indpar%Cband(j)), &
                                   parval(indpar%WSband(j)), &
                                   parval(indpar%WLband(j)),.TRUE. )
        IF (j == indref) THEN
          FORALL (igrid=1:Ngrid) &
            specModel_gen(igrid,:) = EXP( parvec(igrid) ) * dblarr1d(:)
        ELSE
          FORALL (igrid=1:Ngrid) &
            specModel_gen(igrid,:) = EXP( parvec(igrid)+parval(iparef) ) * dblarr1d(:)
        END IF
      ELSE IF (gridCband) THEN
        READ(parname(6:7),*) j
        IF (j == indref) THEN
          FORALL (igrid=1:Ngrid) &
            specModel_gen(igrid,:) = EXP( parval(indpar%lnRband(j)) ) &
                                     * LORENTZBAND( wvl(:),parvec(igrid), &
                                                    parval(indpar%WSband(j)), &
                                                    parval(indpar%WLband(j)),.TRUE. )
        ELSE
          FORALL (igrid=1:Ngrid) &
            specModel_gen(igrid,:) = EXP( parval(indpar%lnRband(j))+parval(iparef) ) &
                                     * LORENTZBAND( wvl(:),parvec(igrid), &
                                                    parval(indpar%WSband(j)), &
                                                    parval(indpar%WLband(j)),.TRUE. )
        END IF
      ELSE IF (gridWSband) THEN
        READ(parname(7:8),*) j
        IF (j == indref) THEN
          FORALL (igrid=1:Ngrid) &
            specModel_gen(igrid,:) = EXP( parval(indpar%lnRband(j)) ) &
                                     * LORENTZBAND( wvl(:),parval(indpar%Cband(j)), &
                                                    parvec(igrid), &
                                                    parval(indpar%WLband(j)),.TRUE. )
        ELSE
          FORALL (igrid=1:Ngrid) &
            specModel_gen(igrid,:) = EXP( parval(indpar%lnRband(j))+parval(iparef) ) &
                                     * LORENTZBAND( wvl(:),parval(indpar%Cband(j)), &
                                                    parvec(igrid), &
                                                    parval(indpar%WLband(j)),.TRUE. )
        END IF
      ELSE IF (gridWLband) THEN
        READ(parname(7:8),*) j
        IF (j == indref) THEN
          FORALL (igrid=1:Ngrid) &
            specModel_gen(igrid,:) = EXP( parval(indpar%lnRband(j)) ) &
                                     * LORENTZBAND( wvl(:),parval(indpar%Cband(j)), &
                                                    parval(indpar%WSband(j)), &
                                                    parvec(igrid),.TRUE. )
        ELSE
          FORALL (igrid=1:Ngrid) &
            specModel_gen(igrid,:) = EXP( parval(indpar%lnRband(j))+parval(iparef) ) &
                                     * LORENTZBAND( wvl(:),parval(indpar%Cband(j)), &
                                                    parval(indpar%WSband(j)), &
                                                    parvec(igrid),.TRUE. )
        END IF
      END IF
     
      DO i=1,Nband
        IF (i /= j) THEN
          IF (i == indref) THEN
            spec1D(:) = spec1D(:) + EXP( parval(indpar%lnRband(i)) ) &
                        * LORENTZBAND( wvl(:),parval(indpar%Cband(i)), &
                                       parval(indpar%WSband(i)), &
                                       parval(indpar%WLband(i)),.TRUE. )
          ELSE
            spec1D(:) = spec1D(:) + EXP( parval(indpar%lnRband(i))+parval(iparef) ) &
                        * LORENTZBAND( wvl(:),parval(indpar%Cband(i)), &
                                       parval(indpar%WSband(i)), &
                                       parval(indpar%WLband(i)),.TRUE. )
          END IF
        END IF
      END DO
    ELSE
      DO i=1,Nband
        IF (i == indref) THEN
          spec1D(:) = spec1D(:) + EXP( parval(indpar%lnRband(i)) ) &
                      * LORENTZBAND( wvl(:),parval(indpar%Cband(i)), &
                                     parval(indpar%WSband(i)), &
                                     parval(indpar%WLband(i)),.TRUE. )
        ELSE
          spec1D(:) = spec1D(:) + EXP( parval(indpar%lnRband(i))+parval(iparef) ) &
                      * LORENTZBAND( wvl(:),parval(indpar%Cband(i)), &
                                     parval(indpar%WSband(i)), &
                                     parval(indpar%WLband(i)),.TRUE. )
        END IF
      END DO
    END IF sampBAND
    
    
    !! 3. Lines
    !!----------
    sampLINE: IF (gridLINE) THEN
      IF (gridlnRline) THEN
        READ(parname(8:9),*) j
        dblarr1d(:) = GAUSSLINE( wvl(:),parval(indpar%Cline(j)), &
                                 parval(indpar%Wline(j)),.TRUE. )
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP( parvec(igrid)+parval(iparef) ) * dblarr1d(:)
      ELSE IF (gridCline) THEN
        READ(parname(6:7),*) j
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP( parval(indpar%lnRline(j))+parval(iparef) ) &
                                   * GAUSSLINE( wvl(:),parvec(igrid), &
                                                parval(indpar%Wline(j)),.TRUE. )
      ELSE IF (gridWline) THEN
        READ(parname(6:7),*) j
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP( parval(indpar%lnRline(j))+parval(iparef) ) &
                                   * GAUSSLINE( wvl(:),parval(indpar%Cline(j)), &
                                                parvec(igrid),.TRUE. )
      END IF
      
      DO i=1,Nline
        IF (i /= j) &
          spec1D(:) = spec1D(:) + EXP( parval(indpar%lnRline(i))+parval(iparef) ) &
                      * GAUSSLINE( wvl(:), parval(indpar%Cline(i)), &
                                   parval(indpar%Wline(i)),.TRUE. )
      END DO
    ELSE
      DO i=1,Nline
        spec1D(:) = spec1D(:) + EXP( parval(indpar%lnRline(i))+parval(iparef) ) &
                    * GAUSSLINE( wvl(:),parval(indpar%Cline(i)), &
                                 parval(indpar%Wline(i)),.TRUE. )
      END DO
    END IF sampLINE
    
    
    !! 4. Stars
    !!----------
    sampSTAR: IF (gridSTAR) THEN
      FORALL (igrid=1:Ngrid) &
        specModel_gen(igrid,:) = EXP(parvec(igrid)) &
                                 * pi * BLACKBODY( nu(:),Tstar ) / MKS%stefan / Tstar**4
    ELSE
      spec1D(:) = spec1D(:) + EXP(parval(indpar%lnFstar)) &
                  * pi * BLACKBODY( nu(:),Tstar ) / MKS%stefan / Tstar**4
    END IF sampSTAR
    
    
    !! 5. Sum all the components
    !!---------------------------
    IF (gridEXTC) THEN
      FORALL (igrid=1:Ngrid) &
        specModel_gen(igrid,:) = spec1D(:) * specModel_gen(igrid,:) * extc1D(:)
    ELSE
      FORALL (igrid=1:Ngrid) &
        specModel_gen(igrid,:) = (specModel_gen(igrid,:)+spec1D(:)) * extc1D(:)
    END IF
      
  END FUNCTION specModel_gen

  
  FUNCTION specModel_scl( wvl, parinfo, parval, indpar, Qabs, extinct)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: wvl
    REAL(DP), DIMENSION(:), INTENT(IN) :: parval ! 1D
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: extinct
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo
    TYPE(indpar_type), INTENT(IN) :: indpar
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN) :: Qabs
    REAL(DP), DIMENSION(SIZE(wvl)) :: specModel_scl
    specModel_scl(:) = RESHAPE( specModel_gen(wvl(:), PARVEC=[0._DP],PARNAME='', &
                                  PARINFO=parinfo, INDPAR=indpar, PARVAL=parval(:), &
                                  QABS=Qabs(:), EXTINCT=extinct(:,:)), &
                                [SIZE(wvl(:))] )
  END FUNCTION specModel_scl
  
END MODULE core
