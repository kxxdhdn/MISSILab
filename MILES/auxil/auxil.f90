!******************************************************************************
!*
!*                              MILES Auxiliaries
!*
!******************************************************************************


  !==========================================================================
  ! 1) AUTHOR: D. HU
  ! 
  ! 2) DESCRIPTION: NONE
  !
  ! 3) HISTORY: 
  !    - 20210305: Simplified specModel_gen (Fred)
  !    - 20210210: Moved extCurve from specModel to read_master
  !    - 20210120: Added description.
  !==========================================================================

MODULE auxil

  USE datable, ONLY: Ncont_max, Nline_max, Nband_max, Nextc_max, Nstar_max, Cband_sig
  USE utilities, ONLY: DP
  USE constants, ONLY: 
  USE inout, ONLY: lenpar
  IMPLICIT NONE
  PRIVATE

  !! only in set_indpar
  INTEGER, PARAMETER, PUBLIC :: NparCONT=2, NparLINE=3, NparBAND=4, NparEXTC=1, NparSTAR=1
  
  PUBLIC :: set_indpar, set_indref, make_Qabs, read_master, initparam
  PUBLIC :: check_SM, invert_SM, invert_mSM
  PUBLIC :: degradeRes, modifBB, gaussLine, lorentzBand, extCurve, specModel

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
    INTEGER, DIMENSION(Ncont_max) :: lnMovd2 = -1 ! dust MBB coeff LOG[Msun/pc2]
    !! If there are identical cont composants, the second lnT becomes delta lnT, etc.
    INTEGER, DIMENSION(Ncont_max) :: lnT = -1 ! dust temperature (MBB) LOG[K]
    INTEGER, DIMENSION(Nline_max) :: lnRline = -1 ! line intensity over ref band Ratio LOG[-]
    INTEGER, DIMENSION(Nline_max) :: Cline = -1 ! Center [um]
    INTEGER, DIMENSION(Nline_max) :: Wline = -1 ! Width [um]
    !! lnRband of ref band is LOG(intensity [W/m2/sr])
    INTEGER, DIMENSION(Nband_max) :: lnRband = -1 ! band intensity over ref band Ratio LOG[-]
    INTEGER, DIMENSION(Nband_max) :: Cband = -1 ! Center (peak) [um]
    INTEGER, DIMENSION(Nband_max) :: WSband = -1 ! Width Short nu side [um]
    INTEGER, DIMENSION(Nband_max) :: WLband = -1 ! Width Long nu side [um]
    INTEGER, DIMENSION(Nextc_max) :: lnAv = -1 ! Extinction LOG[mag]
    !! Total star surface brightness with dilution factor Omega = r/d LOG[W/m2/sr]
    INTEGER, DIMENSION(Nstar_max) :: lnFstar = -1
    INTEGER, DIMENSION(:), ALLOCATABLE :: extra
    INTEGER :: refB = 1 ! index of ref in labB
  END TYPE indpar_type

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
    INTEGER :: i, Nspec, Nr, Nw
    INTEGER, DIMENSION(SIZE(label)) :: indrad ! radius index
    REAL(DP), PARAMETER :: a0 = 1.E-2_DP ! grain radius
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave0 ! READ_OPTICS default grid
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: radius
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Qabsall ! (Nspec, Nr, Nw)
    TYPE(Qabs_type), DIMENSION(SIZE(label)), INTENT(OUT) :: Qabs

    CALL READ_OPTICS(label, WAVE=wave0, RADIUSALL=radius, QABSALL=Qabsall)    

    wave = wave0(:)
    IF (PRESENT(wavAll)) wave = wavAll(:)
    Nw = SIZE(wave(:))

    Nspec = SIZE(label)
    DO i=1,Nspec
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
  SUBROUTINE set_indpar( indpar, parinfo, refB, labB )
    
    USE utilities, ONLY: trimeq, trimlr, pring, strike
    USE arrays, ONLY: iwhere
    IMPLICIT NONE

    TYPE(indpar_type), INTENT(OUT) :: indpar
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo
    CHARACTER(*), INTENT(IN), OPTIONAL :: refB
    CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: labB
    
    INTEGER :: i, Ncont, Nline, Nband, Nextc, Nstar, Nextra
    LOGICAL :: CONTset, LINEset, BANDset, EXTCset, STARset, extraset

    CONTset = ( ANY(TRIMEQ(parinfo(:)%comp, 'CONT')) )
    LINEset = ( ANY(TRIMEQ(parinfo(:)%comp, 'LINE')) )
    BANDset = ( ANY(TRIMEQ(parinfo(:)%comp, 'BAND')) )
    EXTCset = ( ANY(TRIMEQ(parinfo(:)%comp, 'EXTC')) )
    STARset = ( ANY(TRIMEQ(parinfo(:)%comp, 'STAR')) )
    extraset = ( ANY(TRIMEQ(parinfo(:)%comp, 'extra')) )

    IF (CONTset) THEN
      Ncont = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'CONT')) ) / NparCONT
      DO i=1,Ncont
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'lnMovd2'//TRIMLR(PRING(i))), indpar%lnMovd2(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'lnT'//TRIMLR(PRING(i))), indpar%lnT(i) )
        
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
  !! Read the input master file for the Chi2/Bayesian run
  !!
  !!-------------------------------------------------------

  SUBROUTINE read_master( wavAll, dirIN, dirOUT, spec_unit, &
                          Nmcmc, verbose, NiniMC, &!robust_RMS, robust_cal, skew_RMS, &
                          resume, indresume, calib, newseed, newinit, dostop, nohi, &
                          labL, labB, refB, labQ, Qabs, labE, extinct, &
                          Ncont, Nband, Nline, Nextc, Nstar, Nextra, &
                          ! corrhypname, corrname, &
                          parinfo, parmodinfo, parhypinfo, parextinfo, &
                          indpar, Npar, Nparmod, Nparhyp, Ncorrhyp, Ncorr)

    USE utilities, ONLY: DP, strike, warning, trimeq, trimlr, pring!, verbatim
    USE inout, ONLY: lenpar, lenpath, read_input_line, read_hdf5, h5ext
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: incrarr, iwhere, reallocate, closest
    USE interpolation, ONLY: interp_lin_sorted
    USE statistics, ONLY: correl_parlist, N_corr
    USE grain_optics, ONLY: lendustQ, rho_grain, read_optics
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN), OPTIONAL :: dirIN ! default: ./
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: wavAll ! Interp Qabs if presents

    INTEGER, PARAMETER :: unit = 1

    CHARACTER(lenpath) :: filmas, filmod, filext
    ! CHARACTER(lenpath) :: dirOUT0
    ! CHARACTER(lenpar) :: spec_unit0
    CHARACTER(lendustQ) :: refB0
    CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ0, labL0, labB0, labE0
    ! INTEGER :: Nmcmc0, NiniMC0
    INTEGER :: Ncont0, Nband0, Nline0, Nextc0, Nstar0, Nextra0, Npar0
    INTEGER :: Nparmod0, Nparhyp0, Ncorrhyp0, Ncorr0, Nparall
    ! INTEGER :: i, iostat, ipar
    TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE :: parinfall, parinfo0
    ! TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE :: parinfoextra
    ! LOGICAL :: verbose0!, robust_RMS0, robust_cal0, skew_RMS0
    ! LOGICAL :: calib0, newseed0, newinit0, dostop0
    LOGICAL, DIMENSION(:), ALLOCATABLE :: boolparmod, boolpar

    CHARACTER(lenpath), DIMENSION(:), ALLOCATABLE :: Parr1d
    CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: strarr1d, Tarr1d
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
    ! CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
    !   corrhypname, corrname
    INTEGER, INTENT(OUT), OPTIONAL :: Nmcmc, NiniMC, indresume, &
      Ncont, Nband, Nline, Nextc, Nstar, Nextra
    INTEGER, INTENT(OUT), OPTIONAL :: Npar, Nparmod, Nparhyp, Ncorrhyp, Ncorr
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: extinct
    TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Qabs
    TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      parinfo, parmodinfo, parhypinfo, parextinfo
    TYPE(indpar_type), INTENT(OUT), OPTIONAL :: indpar
    LOGICAL, INTENT(OUT), OPTIONAL :: verbose!, robust_RMS, robust_cal, skew_RMS
    LOGICAL, INTENT(OUT), OPTIONAL :: calib, newseed, newinit, dostop, resume, nohi

    !! Read the input master file
    !!----------------------------
    IF (PRESENT(dirIN)) THEN
      filmas = TRIMLR(dirIN)//'input_fitMIR_master'//h5ext
      filmod = TRIMLR(dirIN)//'input_fitMIR_model'//h5ext
      filext = TRIMLR(dirIN)//'input_fitMIR_extra'//h5ext
    ELSE
      CALL WARNING('READ_MASTER','Default input path')
      filmas = './input_fitMIR_master'//h5ext
      filmod = './input_fitMIR_model'//h5ext
      filext = './input_fitMIR_extra'//h5ext
    END IF

    !! Output path
    IF (PRESENT(dirOUT)) THEN
      CALL READ_HDF5(STRARR1D=Parr1d, FILE=filmas, NAME='output dir')
      dirOUT = Parr1d(1)
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
    
    !! input_fitMIR_master
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
    ! IF (PRESENT(robust_RMS)) THEN
    !   CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='robust_RMS')
    !   robust_RMS = TRIMEQ(strarr1d(1),'T')
    ! END IF
    ! IF (PRESENT(robust_cal)) THEN
    !   CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='robust_cal')
    !   robust_cal = TRIMEQ(strarr1d(1),'T')
    ! END IF
    ! IF (PRESENT(skew_RMS)) THEN
    !   CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='skew_RMS')
    !   skew_RMS = TRIMEQ(strarr1d(1),'T')
    ! END IF
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
    
    !! input_fitMIR_model
    IF (PRESENT(spec_unit)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmas, NAME='Spectral unit')
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

    !! input_fitMIR_extra
    ! IF (PRESENT(Nextra)) THEN
    !   CALL READ_HDF5(INTARR1D=intarr1d, FILE=filext, NAME='Nextra')
    !   Nextra0 = intarr1d(1)
    !   Nextra = Nextra0
    ! END IF

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

    ALLOCATE(Tarr1d(Nparall))
    Tarr1d(:) = 'T'

    !! Read input model parameters
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo name')
    parinfall(1:Nparall-Nextra0)%name = strarr1d(:)
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo comp')
    parinfall(1:Nparall-Nextra0)%comp = strarr1d(:)
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo fixed')
    parinfall(1:Nparall-Nextra0)%fixed = TRIMEQ(strarr1d(:),Tarr1d(:))
    CALL READ_HDF5(STRARR2D=strarr2d, FILE=filmod, NAME='parinfo limited')
    CALL READ_HDF5(DBLARR2D=dblarr2d, FILE=filmod, NAME='parinfo limits')
    DO i=1,2
      parinfall(1:Nparall-Nextra0)%limited(i) = TRIMEQ(strarr2d(i,:),Tarr1d(:))
      parinfall(1:Nparall-Nextra0)%limits(i) = dblarr2d(i,:)
    END DO
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo model')
    parinfall(1:Nparall-Nextra0)%model = TRIMEQ(strarr1d(:),Tarr1d(:))
    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo hyper')
    parinfall(1:Nparall-Nextra0)%hyper = TRIMEQ(strarr1d(:),Tarr1d(:))
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
    WHERE (.NOT. TRIMEQ(parinfall(:)%tied, ""))
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
    ! IF (PRESENT(corrhypname)) THEN
    !   ALLOCATE (parhypname(Nparhyp0))
    !   IF (COUNT(boolpar(:)) > 0) THEN
    !     parhypname(:) = PACK(parinfall(:)%name,boolpar(:))
    !     CALL CORREL_PARLIST(parhypname(:),corrhypname)
    !   ELSE
    !     parhypname(:) = PACK(parinfall(:)%name, &
    !                          .NOT. parinfall(:)%fixed .AND. boolparmod(:))
    !     CALL CORREL_PARLIST(parhypname(:),corrhypname)
    !   END IF
    !   DEALLOCATE (parhypname)
    ! END IF

    !! Correlation names for all parameters
    ! IF (PRESENT(corrname)) CALL CORREL_PARLIST(parinfall(:)%name,corrname)
    
    !! Indices
    IF (PRESENT(indpar)) &
      CALL SET_INDPAR(INDPAR=indpar, PARINFO=parinfo0(:), REFB=refB0, LABB=labB0)

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
    
    USE datable, ONLY: TABLine, TABand
    USE utilities, ONLY: DP, trimeq, trimlr, pring, isNaN
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: iwhere, closest, reallocate, incrarr, reverse
    USE inout, ONLY: read_hdf5, lenpar
    USE statistics, ONLY: mean
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NiniMC
    TYPE(indpar_type), INTENT(IN) :: ind
    REAL(DP), DIMENSION(:,:,:,:), INTENT(INOUT) :: par ! (Nx,Ny,Npar,NiniMC)
    TYPE(parinfo_type), DIMENSION(:), INTENT(INOUT) :: parinfo
    INTEGER, DIMENSION(:), INTENT(IN) :: itied
    LOGICAL, DIMENSION(:,:,:), INTENT(IN) :: mask
    LOGICAL, INTENT(IN), OPTIONAL :: newinit
    CHARACTER(*), INTENT(IN), OPTIONAL :: filobs
    CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: labB, labL
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN), OPTIONAL :: Qabs
    
    INTEGER :: i, j, ipar, iw, x, y, Nx, Ny, Nw, Nmc, Nwc, itab, indref
    INTEGER :: Npar, Ncont, Nline, Nband, Nextc, Nstar, Nextra
    INTEGER, DIMENSION(:), ALLOCATABLE :: indwi, indw, indwc
    INTEGER, DIMENSION(SIZE(par,1),SIZE(par,2)) :: iw2
    REAL(DP) :: limi, lims, wcen, dw, val, sig, theta, Tstar
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wOBS, nuOBS
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mu
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: iniparval, FnuOBS, val3
    REAL(DP), DIMENSION(SIZE(par,1),SIZE(par,2)) :: limi2, lims2, val2
    REAL(DP), DIMENSION(SIZE(par,1),SIZE(par,2),MAX(NiniMC,1)) :: theta3, lnIref
    LOGICAL :: newini, chi2ini
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: maskxy
    CHARACTER(lenpar) :: spec_unit
    CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: iniparname, strarr1d

    !! Preliminaries
    Nx = SIZE(par(:,:,:,:),1)
    Ny = SIZE(par(:,:,:,:),2)
    Npar = SIZE(parinfo(:))
    Ncont = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'CONT')) ) / NparCONT
    Nline = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'LINE')) ) / NparLINE
    Nband = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'BAND')) ) / NparBAND
    Nextc = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'EXTC')) ) / NparEXTC
    Nstar = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'STAR')) ) / NparSTAR
    IF (PRESENT(newinit)) THEN
      newini = newinit
    ELSE
      newini = .FALSE.
    END IF
    ALLOCATE(maskxy(Nx,Ny))
    FORALL (x=1:Nx,y=1:Ny) maskxy(x,y) = ALL(mask(x,y,:))
    !! ANY is too strict, e.g. for spectra interpolated to the common grid

    Tstar = 5.E4_DP ! [K] high enough to stay in Rayleigh-Jeans limit (T>>5000K for 3um)

    CALL READ_HDF5(STRARR1D=strarr1d, FILE=filobs, NAME='Spectral unit')
    spec_unit = strarr1d(1)
    CALL READ_HDF5(DBLARR1D=wOBS, FILE=filobs, NAME='Wavelength (microns)')
    CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filobs, NAME='FnuOBS ('//TRIMLR(spec_unit)//')')

    Nw = SIZE(wOBS(:))
    ALLOCATE(nuOBS(Nw), val3(Nx,Ny,Nw))
    nuOBS(:) = MKS%clight / wOBS(:) /MKS%micron

    indref = ind%refB
    lnIref(:,:,:) = 1._DP
    
    newinipar: IF (.NOT. newini) THEN
      
      !!----------------------------------------------------
      !! I. Automatically generate initial parameter values
      !!----------------------------------------------------

      IF (NiniMC==0) THEN

        !! a. Simple initilization
        !!-------------------------

        DO i=1,Nband
          CALL IWHERE( TRIMEQ(TABand(:)%label, labB(i)), itab ) ! itab = ind of TABand
  
          !! Cband
          IF (parinfo(ind%Cband(i))%fixed) THEN
            par(:,:,ind%Cband(i),:) = parinfo(ind%Cband(i))%value
          ELSE
            val = TABand(itab)%wave
            par(:,:,ind%Cband(i),:) = val
            !! If Cband is not fixed, limits are indispensable
            sig = Cband_sig
            IF (.NOT. parinfo(ind%Cband(i))%limited(1)) &
              parinfo(ind%Cband(i))%limits(1) = val-sig
            IF (.NOT. parinfo(ind%Cband(i))%limited(2)) &
              parinfo(ind%Cband(i))%limits(2) = val+sig
            parinfo(ind%Cband(i))%limited = .TRUE.
          END IF
          
          !! WSband
          IF (parinfo(ind%WSband(i))%fixed) THEN
            par(:,:,ind%WSband(i),:) = parinfo(ind%WSband(i))%value
          ELSE
            val = TABand(itab)%sigmaS
            par(:,:,ind%WSband(i),:) = val
            IF (.NOT. parinfo(ind%WSband(i))%limited(1)) &
              parinfo(ind%WSband(i))%limits(1) = val*0.5_DP
            IF (.NOT. parinfo(ind%WSband(i))%limited(2)) &
              parinfo(ind%WSband(i))%limits(2) = val*2._DP
            parinfo(ind%WSband(i))%limited = .TRUE.
          END IF
     
          !! WLband
          IF (parinfo(ind%WLband(i))%fixed) THEN
            par(:,:,ind%WLband(i),:) = parinfo(ind%WLband(i))%value
          ELSE
            val = TABand(itab)%sigmaL
            par(:,:,ind%WLband(i),:) = val
            IF (.NOT. parinfo(ind%WLband(i))%limited(1)) &
              parinfo(ind%WLband(i))%limits(1) = val*0.5_DP
            IF (.NOT. parinfo(ind%WLband(i))%limited(2)) &
              parinfo(ind%WLband(i))%limits(2) = val*2._DP
            parinfo(ind%WLband(i))%limited = .TRUE.
          END IF
          
          !! lnRband
          IF (parinfo(ind%lnRband(i))%fixed) THEN
            par(:,:,ind%lnRband(i),:) = parinfo(ind%lnRband(i))%value
          ELSE
            FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
              iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cband(i),1) )
              !! Same estimation as lnRline
              par(x,y,ind%lnRband(i),:) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) / &
                MAXVAL(LORENTZBAND(wOBS(:),par(x,y,ind%Cband(i),1), &
                                   par(x,y,ind%WSband(i),1), &
                                   par(x,y,ind%WLband(i),1),.TRUE.)) )
            END FORALL
            
            IF (i == indref) lnIref(:,:,:) = par(:,:,ind%lnRband(i),:)
          END IF

          !! Find indw
          limi = par(1,1,ind%Cband(i),1) - par(1,1,ind%WLband(i),1)
          lims = par(1,1,ind%Cband(i),1) + par(1,1,ind%WSband(i),1)
          IF (ANY(wOBS(:)>limi .AND. wOBS(:)<lims)) THEN
            CALL IWHERE( wOBS(:)>limi .AND. wOBS(:)<lims, indwi )
            DO iw=1,SIZE(indwi)
              CALL INCRARR( indw, indwi(iw) )
            END DO
          END IF

        END DO

        DO i=1,Nband
          IF (i /= indref) &
            par(:,:,ind%lnRband(i),:) = par(:,:,ind%lnRband(i),:) - lnIref(:,:,:)
        END DO
        
        
        DO i=1,Nline
          CALL IWHERE( TRIMEQ(TABLine(:)%label, labL(i)), itab ) ! itab = ind of TABLine

          !! Cline
          IF (parinfo(ind%Cline(i))%fixed) THEN
            par(:,:,ind%Cline(i),:) = parinfo(ind%Cline(i))%value
          ELSE
            val = TABLine(itab)%wave
            par(:,:,ind%Cline(i),:) = val
            !! If Cline is not fixed, limits are indispensable
            !! FWHM = 2*sqrt(2*log(2))*sigma ~ 2.355 (suppose PSF is gaussian)
            !! R = lambda / FWHM
            !! => sigma = lambda / R / 2.355
            sig = degradeRes(val, 0.01_DP, 'SL-LL') / 2.355_DP
            IF (.NOT. parinfo(ind%Cline(i))%limited(1)) &
              parinfo(ind%Cline(i))%limits(1) = val-sig
            IF (.NOT. parinfo(ind%Cline(i))%limited(2)) &
              parinfo(ind%Cline(i))%limits(2) = val+sig
            parinfo(ind%Cline(i))%limited = .TRUE.
          END IF
          
          !! Wline
          IF (parinfo(ind%Wline(i))%fixed) THEN
            par(:,:,ind%Wline(i),:) = parinfo(ind%Wline(i))%value
          ELSE
            val = degradeRes(TABLine(itab)%wave, .01_DP, 'SL-LL')
            par(:,:,ind%Wline(i),:) = val
            IF (.NOT. parinfo(ind%Wline(i))%limited(1)) &
              parinfo(ind%Wline(i))%limits(1) = val*0.5_DP
            IF (.NOT. parinfo(ind%Wline(i))%limited(2)) &
              parinfo(ind%Wline(i))%limits(2) = val*2._DP
            parinfo(ind%Wline(i))%limited = .TRUE.
          END IF
          
          !! lnRline
          IF (parinfo(ind%lnRline(i))%fixed) THEN
            par(:,:,ind%lnRline(i),:) = parinfo(ind%lnRline(i))%value
          ELSE
            FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
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
          
          !! Find indw
          !! Take (1,1) Cline and Wline since they are supposed to vary little
          limi = par(1,1,ind%Cline(i),1) - par(1,1,ind%Wline(i),1)
          lims = par(1,1,ind%Cline(i),1) + par(1,1,ind%Wline(i),1)
          IF (ANY(wOBS(:)>limi .AND. wOBS(:)<lims)) THEN
            CALL IWHERE( wOBS(:)>limi .AND. wOBS(:)<lims, indwi )
            DO iw=1,SIZE(indwi)
              CALL INCRARR( indw, indwi(iw) )
            END DO
          END IF
          
        END DO

        !! Create wvl grid with indices (indwc) for band & line free spectra (cont only)
        Nwc = Nw - SIZE(indw(:)) ! Nb of wvl for the cont grid
        DO iw=1,Nw
          IF ( ALL(indw(:).NE.iw) ) THEN
            CALL INCRARR( indwc, iw )

          END IF
        END DO
        !! Larger T corresponds to smaller wvl
        indwc(:) = REVERSE(indwc(:))
        ! CALL REALLOCATE(indw, Nwc)
        ! indw(:) = indwc(:)
        ! indwc(:) = indw(Nwc:1:-1)
        
        DO i=1,Ncont
          !! lnT LOG[K]
          IF (parinfo(ind%lnT(i))%fixed) THEN
            par(:,:,ind%lnT(i),:) = parinfo(ind%lnT(i))%value
          ELSE
            !! Create lnT grid: (ln50,ln500) - (Wien Law: 58-5.8um)
            par(:,:,ind%lnT(i),:) = (LOG(500._DP)-LOG(50._DP)) * &
                                     REAL(i)/Ncont + LOG(50._DP) !!! i/Ncont
          END IF
          
          !! lnMovd2 LOG[Msun/pc2]
          IF ( parinfo(ind%lnMovd2(i))%fixed ) THEN
            par(:,:,ind%lnMovd2(i),:) = parinfo(ind%lnMovd2(i))%value
          ELSE
            !! lnMovd2 is auto determined by modifBB(lnT) & the spectrum to fit
            IF (Nwc/Ncont .GE. 1) THEN
              iw = INT(Nwc/Ncont) * i ! iw is NOT wOBS's index, but is indwc's index
            ELSE
              !! Case where there are more CONT compo than its wvl grid (Ncont>Nwc)
              IF (i>Nwc) THEN
                iw = Nwc
              ELSE
                iw = i
              END IF
            END IF
            FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
              val3(x,y,:) = modifBB(wOBS(:),EXP(par(x,y,ind%lnT(i),1)),Qabs(i),.TRUE.) * &
                            MKS%Msun/MKS%pc**2
              par(x,y,ind%lnMovd2(i),:) = LOG( ABS(FnuOBS(x,y,indwc(iw))) / &
                                                   val3(x,y,indwc(iw)) / Ncont )
            END FORALL
          END IF

        END DO

        DO i=1,Nextc
          !! lnAv LOG[mag]
          IF (parinfo(ind%lnAv(i))%fixed) THEN
            par(:,:,ind%lnAv(i),:) = parinfo(ind%lnAv(i))%value
          ELSE
            par(:,:,ind%lnAv(i),:) = 0._DP ! 1 [mag] = no extinction
            IF (.NOT. parinfo(ind%lnAv(i))%limited(1)) &
              parinfo(ind%lnAv(i))%limits(1) = 0._DP
            IF (.NOT. parinfo(ind%lnAv(i))%limited(2)) &
              parinfo(ind%lnAv(i))%limits(2) = 1._DP
            parinfo(ind%lnAv(i))%limited = .TRUE.
          END IF

        END DO

        DO i=1,Nstar
          !! lnFstar LOG[Lsun/pc2]
          IF (parinfo(ind%lnFstar(i))%fixed) THEN
            par(:,:,ind%lnFstar(i),:) = parinfo(ind%lnFstar(i))%value
          ELSE
            !! lnFstar is auto determined by BB & the spectrum to fit
            iw = Nwc - 10 ! avoid wvl edge
            FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
              val3(x,y,:) = BLACKBODY(nuOBS(:),Tstar) / MKS%stefan/Tstar**4
              val2(x,y) = LOG( ABS(FnuOBS(x,y,indwc(iw)))/val3(x,y,indwc(iw)) /pi )
              par(x,y,ind%lnFstar(i),:) = val2(x,y)
            END FORALL
          END IF

        END DO
        
      ELSE

        !! b. MC initialization
        !!----------------------

        DO i=1,Nband
          CALL IWHERE( TRIMEQ(TABand(:)%label, labB(i)), itab ) ! itab = ind of TABand
  
          !! Cband
          IF (parinfo(ind%Cband(i))%fixed) THEN
            par(:,:,ind%Cband(i),:) = parinfo(ind%Cband(i))%value
          ELSE
            val = TABand(itab)%wave
            ! sig = degradeRes(val, 0.01_DP, 'SL-LL') / 2.355_DP
            sig = Cband_sig
            limi = MERGE( parinfo(ind%Cband(i))%limits(1), &
                          val-sig, &
                          parinfo(ind%Cband(i))%limited(1) ) ! lim inf
            lims = MERGE( parinfo(ind%Cband(i))%limits(2), &
                          val+sig, &
                          parinfo(ind%Cband(i))%limited(2) ) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%Cband(i),j) = val
              ELSE
                !! Uniform distribution
                par(:,:,ind%Cband(i),j) = (lims - limi) * theta3(:,:,j) + limi
              END IF
            END DO
            IF (.NOT. parinfo(ind%Cband(i))%limited(1)) &
              parinfo(ind%Cband(i))%limits(1) = val-sig
            IF (.NOT. parinfo(ind%Cband(i))%limited(2)) &
              parinfo(ind%Cband(i))%limits(2) = val+sig
            parinfo(ind%Cband(i))%limited = .TRUE.
          END IF
          
          !! WSband
          IF (parinfo(ind%WSband(i))%fixed) THEN
            par(:,:,ind%WSband(i),:) = parinfo(ind%WSband(i))%value
          ELSE
            val = TABand(itab)%sigmaS
            limi = MERGE( parinfo(ind%WSband(i))%limits(1), &
                          val*0.5_DP, &
                          parinfo(ind%WSband(i))%limited(1) ) ! lim inf
            lims = MERGE( parinfo(ind%WSband(i))%limits(2), &
                          val*2._DP, &
                          parinfo(ind%WSband(i))%limited(2) ) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%WSband(i),j) = val
              ELSE
                !! Uniform distribution
                par(:,:,ind%WSband(i),j) = (lims - limi) * theta3(:,:,j) + limi
              END IF
            END DO
            IF (.NOT. parinfo(ind%WSband(i))%limited(1)) &
              parinfo(ind%WSband(i))%limits(1) = val*0.5_DP
            IF (.NOT. parinfo(ind%WSband(i))%limited(2)) &
              parinfo(ind%WSband(i))%limits(2) = val*2._DP
            parinfo(ind%WSband(i))%limited = .TRUE.
          END IF
     
          !! WLband
          IF (parinfo(ind%WLband(i))%fixed) THEN
            par(:,:,ind%WLband(i),:) = parinfo(ind%WLband(i))%value
          ELSE
            val = TABand(itab)%sigmaL
            limi = MERGE( parinfo(ind%WLband(i))%limits(1), &
                          val*0.5_DP, &
                          parinfo(ind%WLband(i))%limited(1) ) ! lim inf
            lims = MERGE( parinfo(ind%WLband(i))%limits(2), &
                          val*2._DP, &
                          parinfo(ind%WLband(i))%limited(2) ) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%WLband(i),j) = val
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
  
          !! lnRband
          IF (parinfo(ind%lnRband(i))%fixed) THEN
            par(:,:,ind%lnRband(i),:) = parinfo(ind%lnRband(i))%value
          ELSE
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
                  iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cband(i),j) )
                  !! Estimated via feature peak
                  par(x,y,ind%lnRband(i),j) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) / &
                    MAXVAL(LORENTZBAND(wOBS(:),par(x,y,ind%Cband(i),j), &
                                       par(x,y,ind%WSband(i),j), &
                                       par(x,y,ind%WLband(i),j),.TRUE.)) )
  
                END FORALL
              ELSE
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
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
                    
          !! Find indw
          limi = par(1,1,ind%Cband(i),1) - par(1,1,ind%WLband(i),1)
          lims = par(1,1,ind%Cband(i),1) + par(1,1,ind%WSband(i),1)
          IF (ANY(wOBS(:)>limi .AND. wOBS(:)<lims)) THEN
            CALL IWHERE( wOBS(:)>limi .AND. wOBS(:)<lims, indwi )
            DO iw=1,SIZE(indwi)
              CALL INCRARR( indw, indwi(iw) )
            END DO
          END IF
          
        END DO

        DO i=1,Nband
          IF (i /= indref) &
            par(:,:,ind%lnRband(i),:) = par(:,:,ind%lnRband(i),:) - lnIref(:,:,:)
        END DO
        
        DO i=1,Nline
          CALL IWHERE( TRIMEQ(TABLine(:)%label, labL(i)), itab ) ! itab = ind of TABLine
          
          !! Cline
          IF (parinfo(ind%Cline(i))%fixed) THEN
            par(:,:,ind%Cline(i),:) = parinfo(ind%Cline(i))%value
          ELSE
            !! FWHM = 2*sqrt(2*log(2))*sigma ~ 2.355 (suppose PSF is gaussian)
            !! R = lambda / FWHM
            !! => sigma = lambda / R / 2.355
            val = TABLine(itab)%wave
            sig = degradeRes(val, 0.01_DP, 'SL-LL') / 2.355_DP
            limi = MERGE( parinfo(ind%Cline(i))%limits(1), &
                          val-sig, &
                          parinfo(ind%Cline(i))%limited(1)) ! lim inf
            lims = MERGE( parinfo(ind%Cline(i))%limits(2), &
                          val+sig, &
                          parinfo(ind%Cline(i))%limited(2)) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%Cline(i),j) = val
              ELSE
                !! Uniform distribution
                par(:,:,ind%Cline(i),j) = (lims - limi) * theta3(:,:,j) + limi
              END IF
            END DO
            IF (.NOT. parinfo(ind%Cline(i))%limited(1)) &
              parinfo(ind%Cline(i))%limits(1) = val-sig
            IF (.NOT. parinfo(ind%Cline(i))%limited(2)) &
              parinfo(ind%Cline(i))%limits(2) = val+sig
            parinfo(ind%Cline(i))%limited = .TRUE.
          END IF
          
          !! Wline
          IF (parinfo(ind%Wline(i))%fixed) THEN
            par(:,:,ind%Wline(i),:) = parinfo(ind%Wline(i))%value
          ELSE
            val = degradeRes(TABLine(itab)%wave, .01_DP, 'SL-LL')
            limi = MERGE( parinfo(ind%Wline(i))%limits(1), &
                          val*0.5_DP, &
                          parinfo(ind%Wline(i))%limited(1)) ! lim inf
            lims = MERGE( parinfo(ind%Wline(i))%limits(2), &
                          val*2._DP, &
                          parinfo(ind%Wline(i))%limited(2)) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%Wline(i),j) = val
              ELSE
                !! Uniform distribution
                par(:,:,ind%Wline(i),j) = (lims - limi) * theta3(:,:,j) + limi
              END IF
            END DO
            IF (.NOT. parinfo(ind%Wline(i))%limited(1)) &
              parinfo(ind%Wline(i))%limits(1) = val*0.5_DP
            IF (.NOT. parinfo(ind%Wline(i))%limited(2)) &
              parinfo(ind%Wline(i))%limits(2) = val*2._DP
            parinfo(ind%Wline(i))%limited = .TRUE.
          END IF
          
          !! lnRline
          IF (parinfo(ind%lnRline(i))%fixed) THEN
            par(:,:,ind%lnRline(i),:) = parinfo(ind%lnRline(i))%value
          ELSE
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
                  iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cline(i),j) )
                  !! Estimated via feature peak
                  par(x,y,ind%lnRline(i),j) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) / &
                    MAXVAL(GAUSSLINE(wOBS(:),par(x,y,ind%Cline(i),j), &
                                     par(x,y,ind%Wline(i),j),.TRUE.)) ) - &
                    lnIref(x,y,j)

                END FORALL
              ELSE
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
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
          
          !! Find indw
          limi = par(1,1,ind%Cline(i),1) - par(1,1,ind%Wline(i),1)
          lims = par(1,1,ind%Cline(i),1) + par(1,1,ind%Wline(i),1)
          IF (ANY(wOBS(:)>limi .AND. wOBS(:)<lims)) THEN
            CALL IWHERE( wOBS(:)>limi .AND. wOBS(:)<lims, indwi )
            DO iw=1,SIZE(indwi)
              CALL INCRARR( indw, indwi(iw) )
            END DO
          END IF
          
        END DO

        !! Create wvl grid with indices (indwc) for band & line free spectra (cont only)
        Nwc = Nw - SIZE(indw(:)) ! Nb of wvl for the cont grid
        DO iw=1,Nw
          IF ( ALL(indw(:).NE.iw) ) THEN
            CALL INCRARR( indwc, iw )

          END IF
        END DO
        !! Larger T corresponds to smaller wvl
        indwc(:) = REVERSE(indwc(:))
        ! CALL REALLOCATE(indw, Nwc)
        ! indw(:) = indwc(:)
        ! indwc(:) = indw(Nwc:1:-1)
        
        DO i=1,Ncont
          !! lnT LOG[K]
          IF (parinfo(ind%lnT(i))%fixed) THEN
            par(:,:,ind%lnT(i),:) = parinfo(ind%lnT(i))%value
          ELSE
            !! Generate a random value on lnT grid: (ln50,ln500)
            CALL RANDOM_NUMBER(theta)
            val = (LOG(500._DP)-LOG(50._DP))*theta + LOG(50._DP) !!! i/Ncont
            limi = MERGE( parinfo(ind%lnT(i))%limits(1), &
                          val, &
                          parinfo(ind%lnT(i))%limited(1)) ! lim inf
            lims = MERGE( parinfo(ind%lnT(i))%limits(2), &
                          val, &
                          parinfo(ind%lnT(i))%limited(2)) ! lim sup
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                par(:,:,ind%lnT(i),j) = (LOG(500._DP)-LOG(50._DP)) * &
                                         REAL(i)/Ncont + LOG(50._DP) !!! i/Ncont
              ELSE
                !! Uniform distribution
                par(:,:,ind%lnT(i),j) = (lims - limi) * theta3(:,:,j) + limi
              END IF
            END DO
          END IF
          
          !! lnMovd2 LOG[Msun/pc2]
          IF (parinfo(ind%lnMovd2(i))%fixed) THEN
            par(:,:,ind%lnMovd2(i),:) = parinfo(ind%lnMovd2(i))%value
          ELSE
            CALL RANDOM_NUMBER(theta3(:,:,:))
            !! lnMovd2 is auto determined by modifBB(lnT) & the spectrum to fit
            IF (Nwc/Ncont .GE. 1) THEN
              iw = INT(Nwc/Ncont) * i ! iw is NOT wOBS's index, but is indwc's index
            ELSE
              !! Case where there are more CONT compo than its wvl grid (Ncont>Nwc)
              IF (i>Nwc) THEN
                iw = Nwc
              ELSE
                iw = i
              END IF
            END IF
            DO j=1,NiniMC
              IF (j == 1) THEN
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
                  val3(x,y,:) = modifBB(wOBS(:),EXP(par(x,y,ind%lnT(i),1)),Qabs(i),.TRUE.) &
                                * MKS%Msun/MKS%pc**2
                  val2(x,y) = LOG( ABS(FnuOBS(x,y,indwc(iw)))/val3(x,y,indwc(iw)) )
                  par(x,y,ind%lnMovd2(i),j) = val2(x,y)

                END FORALL
              ELSE
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
                  val3(x,y,:) = modifBB(wOBS(:),EXP(par(x,y,ind%lnT(i),1)),Qabs(i),.TRUE.) &
                                * MKS%Msun/MKS%pc**2
                  val2(x,y) = LOG( ABS(FnuOBS(x,y,indwc(iw)))/val3(x,y,indwc(iw)) )
                  limi2(x,y) = MERGE( parinfo(ind%lnMovd2(i))%limits(1), &
                                      val2(x,y), &
                                      parinfo(ind%lnMovd2(i))%limited(1) ) ! lim inf
                  lims2(x,y) = MERGE( parinfo(ind%lnMovd2(i))%limits(2), &
                                      val2(x,y), &
                                      parinfo(ind%lnMovd2(i))%limited(2) ) ! lim sup
                  !! Uniform distribution
                  par(x,y,ind%lnMovd2(i),j) = (lims2(x,y) - limi2(x,y)) * theta3(x,y,j) &
                                              + limi2(x,y)
                  
                END FORALL
              END IF
            END DO
          END IF
          
        END DO
        
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
              parinfo(ind%lnAv(i))%limits(2) = 1._DP
            parinfo(ind%lnAv(i))%limited = .TRUE.
          END IF

        END DO

        DO i=1,Nstar
          !! lnFstar LOG[Lsun/pc2]
          IF (parinfo(ind%lnFstar(i))%fixed) THEN
            par(:,:,ind%lnFstar(i),:) = parinfo(ind%lnFstar(i))%value
          ELSE
            CALL RANDOM_NUMBER(theta3(:,:,:))
            !! lnFstar is auto determined by BB & the spectrum to fit
            iw = Nwc - 10 ! avoid wvl edge
            DO j=1,NiniMC
              IF (j == 1) THEN
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
                  val3(x,y,:) = BLACKBODY(nuOBS(:),Tstar) / MKS%stefan/Tstar**4
                  val2(x,y) = LOG( ABS(FnuOBS(x,y,indwc(iw)))/val3(x,y,indwc(iw)) /pi )
                  par(x,y,ind%lnFstar(i),j) = val2(x,y)
                END FORALL
              ELSE
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
                  val3(x,y,:) = BLACKBODY(nuOBS(:),Tstar) / MKS%stefan/Tstar**4
                  val2(x,y) = LOG( ABS(FnuOBS(x,y,indwc(iw)))/val3(x,y,indwc(iw)) /pi )
                  limi2(x,y) = MERGE( parinfo(ind%lnFstar(i))%limits(1), &
                                      val2(x,y)-.1_DP, & ! EXP(factor - 3)
                                      parinfo(ind%lnFstar(i))%limited(1) ) ! lim inf
                  lims2(x,y) = MERGE( parinfo(ind%lnFstar(i))%limits(2), &
                                      val2(x,y)+.1_DP, &
                                      parinfo(ind%lnFstar(i))%limited(2) ) ! lim sup
                  !! Uniform distribution
                  par(x,y,ind%lnFstar(i),j) = (lims2(x,y)-limi2(x,y)) * &
                                              theta3(x,y,j) + limi2(x,y)
    
                END FORALL

              END IF
            END DO
          END IF
        END DO

      END IF
      DEALLOCATE(indwi, indw, indwc, val3)
      
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

        !! Set hard limits for intensive parameters
        !!------------------------------------------
        DO i=1,Nline
          CALL IWHERE( TRIMEQ(TABLine(:)%label, labL(i)), itab ) ! itab = ind of TABLine
          !! FWHM = 2*sqrt(2*log(2))*sigma ~ 2.355 (suppose PSF is gaussian)
          !! R = lambda / FWHM
          !! => sigma = lambda / R / 2.355
          wcen = TABLine(itab)%wave
          dw = degradeRes(TABLine(itab)%wave, .01_DP, 'SL-LL')
          sig = dw / 2.355_DP

          !! Cline
          par(:,:,ind%Cline(i),:) = wcen
          IF (.NOT. parinfo(ind%Cline(i))%limited(1)) &
            parinfo(ind%Cline(i))%limits(1) = wcen-sig
          IF (.NOT. parinfo(ind%Cline(i))%limited(2)) &
            parinfo(ind%Cline(i))%limits(2) = wcen+sig
          parinfo(ind%Cline(i))%limited = .TRUE.
          
          !! Wline
          par(:,:,ind%Wline(i),:) = dw
          IF (.NOT. parinfo(ind%Wline(i))%limited(1)) &
            parinfo(ind%Wline(i))%limits(1) = dw*0.5_DP
          IF (.NOT. parinfo(ind%Wline(i))%limited(2)) &
            parinfo(ind%Wline(i))%limits(2) = dw*2._DP
          parinfo(ind%Wline(i))%limited = .TRUE.
          
        END DO

        DO i=1,Nband
          CALL IWHERE( TRIMEQ(TABand(:)%label, labB(i)), itab ) ! itab = ind of TABand
          wcen = TABand(itab)%wave
          sig = Cband_sig
          
          !! Cband
          par(:,:,ind%Cband(i),:) = wcen
          IF (.NOT. parinfo(ind%Cband(i))%limited(1)) &
            parinfo(ind%Cband(i))%limits(1) = wcen-sig
          IF (.NOT. parinfo(ind%Cband(i))%limited(2)) &
            parinfo(ind%Cband(i))%limits(2) = wcen+sig
          parinfo(ind%Cband(i))%limited = .TRUE.
          
          !! WSband
          val = TABand(itab)%sigmaS
          par(:,:,ind%WSband(i),:) = val
          IF (.NOT. parinfo(ind%WSband(i))%limited(1)) &
            parinfo(ind%WSband(i))%limits(1) = val*0.5_DP
          IF (.NOT. parinfo(ind%WSband(i))%limited(2)) &
            parinfo(ind%WSband(i))%limits(2) = val*2._DP
          parinfo(ind%WSband(i))%limited = .TRUE.
   
          !! WLband
          val = TABand(itab)%sigmaL
          par(:,:,ind%WLband(i),:) = val
          IF (.NOT. parinfo(ind%WLband(i))%limited(1)) &
            parinfo(ind%WLband(i))%limits(1) = val*0.5_DP
          IF (.NOT. parinfo(ind%WLband(i))%limited(2)) &
            parinfo(ind%WLband(i))%limits(2) = val*2._DP
          parinfo(ind%WLband(i))%limited = .TRUE.

        END DO
        
      ELSE
        CALL READ_HDF5(DBLARR3D=iniparval, FILE=filobs, &
                       NAME='Initial parameter value')
      END IF
      
      !! Replace the automatic initial values by those
      DO i=1,Npar
        IF ( ANY(TRIMEQ(iniparname(:),parinfo(i)%name)) ) THEN
          CALL IWHERE(TRIMEQ(iniparname(:),parinfo(i)%name), ipar)
          WHERE (maskxy(:,:)) par(:,:,i,1) = iniparval(:,:,ipar)
          
        END IF
      END DO

      !! Free memory space
      DEALLOCATE (iniparname,iniparval)
    
    END IF newinipar

    !! Check that the initial guesses satisfy the parameter constraints
    FORALL (i=1:Npar,j=1:SIZE(par,4))
      par(:,:,i,j) = MERGE(parinfo(i)%value,par(:,:,i,j),parinfo(i)%fixed)
      WHERE (parinfo(i)%limited(1) .AND. (par(:,:,i,j) < parinfo(i)%limits(1)) &
             .AND. maskxy(:,:))
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
        WHERE (isNaN(par(:,:,i,j)) .AND. maskxy(:,:)) par(:,:,i,j) = mu(i,j)
      END FORALL
    END IF
    
    !! Enforce tying of parameters
    FORALL (i=1:Npar,itied(i) > 0) par(:,:,itied(i),:) = par(:,:,i,:)
    
  END SUBROUTINE initparam
  
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
  
  FUNCTION invert_SM( A, B_prev, pos, delta, Ndim, determinant )

    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, DIMENSION(2), INTENT(IN) :: pos
    REAL(DP), INTENT(IN) :: delta
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: A, B_prev
    INTEGER, INTENT(IN), OPTIONAL :: Ndim
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
    c0 = 1 + B_prev(pos(2),pos(1)) * delta

    IF (c0/=0._DP) THEN
      invert_SM(pos(2),:) = B_prev(pos(2),:) / c0
      invert_SM(:,pos(1)) = B_prev(:,pos(1)) / c0

      c1 = delta / c0
      FORALL (y=1:Ny,y/=pos(1)) &
        invert_SM(:,y) = B_prev(:,y) - B_prev(:,pos(1))*B_prev(pos(2),y)*c1

    END IF    
    
    IF (PRESENT(determinant)) PRINT*, 'COUCOU'

  END FUNCTION invert_SM

  !!-------------------------------------------------------
  !!
  !! Modified Sherman-Morrison approach to invert matrices
  !!
  !!-------------------------------------------------------
  FUNCTION invert_mSM( A, B_prev, pos, delta, Ndim, determinant )

    USE utilities, ONLY: DP, strike
    IMPLICIT NONE
    
    INTEGER, DIMENSION(2), INTENT(IN) :: pos
    REAL(DP), INTENT(IN) :: delta
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: A, B_prev
    INTEGER, INTENT(IN), OPTIONAL :: Ndim
    INTEGER :: Ny
    INTEGER, DIMENSION(2) :: pos_trans
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: A1, B1
    REAL(DP), INTENT(OUT), OPTIONAL :: determinant
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: invert_mSM

    IF (PRESENT(Ndim)) THEN
      Ny = Ndim
    ELSE
      Ny = SIZE(A,2)
    END IF
    
    !! Create intermediate matrix A1
    A1(:,:) = A(:,:)
    A1(pos(1),pos(2)) = A(pos(1),pos(2)) - delta

    pos_trans(:) = [pos(2), pos(1)]

    !! Calculate invert matrix
    B1(:,:) = invert_SM(A1, B_prev, pos_trans, delta, Ny)
    invert_mSM(:,:) = invert_SM(A, B1, pos, delta, Ny)
    
    IF (PRESENT(determinant)) PRINT*, 'COUCOU'

  END FUNCTION invert_mSM
  
  !!-------------------------------------------------------
  !!
  !! Automatize the degradation of the spectral resolution
  !!
  !!-------------------------------------------------------
  FUNCTION degradeRes( wc, dw, instr )

    USE datable, ONLY: res
    USE utilities, ONLY: DP
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)         :: instr
    REAL(DP), INTENT(IN)             :: wc, dw
    REAL(DP)                         :: degradeRes
    
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
        IF (wc.LT.12._DP) THEN
          degradeRes = res%dw_w_SL * wc
        ELSE
          degradeRes = res%dw_w_LL * wc
        END IF
      CASE('Ns-SL')
        IF (wc.LT.5._DP) THEN
          degradeRes = res%dw_w_AKARI_Ns * wc
        ELSE
          degradeRes = res%dw_w_SL * wc
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
  PURE FUNCTION modifBB_gen( x, temp, Qabs, wi )
    !! temp [K]; BLACKBODY [W/m2/sr/Hz]; 
    !! output array [W/sr/Hz/kg]
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: temp
    TYPE(Qabs_type), INTENT(IN), DIMENSION(:) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL :: wi
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
      modifBB_gen(itemp,:) = Qabs(itemp)%kappa * BLACKBODY(nu(:), temp(itemp))
    
  END FUNCTION modifBB_gen

  PURE FUNCTION modifBB_0( x, temp, Qabs, wi )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: temp
    TYPE(Qabs_type), INTENT(IN) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL :: wi
    REAL(DP), DIMENSION(SIZE(x)) :: modifBB_0
    IF (PRESENT(wi)) THEN
      modifBB_0(:) = RESHAPE(modifBB_gen(x(:),[temp],[Qabs],wi), [SIZE(x(:))])
    ELSE
      modifBB_0(:) = RESHAPE(modifBB_gen(x(:),[temp],[Qabs]), [SIZE(x(:))])
    END IF
  END FUNCTION modifBB_0

  PURE FUNCTION modifBB_1( x, temp, Qabs, wi )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: temp
    TYPE(Qabs_type), INTENT(IN), DIMENSION(:) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL :: wi
    REAL(DP), DIMENSION(SIZE(temp),SIZE(x)) :: modifBB_1
    IF (PRESENT(wi)) THEN
      modifBB_1(:,:) = RESHAPE(modifBB_gen(x(:),temp(:),Qabs(:),wi), &
                               [SIZE(temp(:)),SIZE(x(:))])
    ELSE
      modifBB_1(:,:) = RESHAPE(modifBB_gen(x(:),temp(:),Qabs(:)), &
                               [SIZE(temp(:)),SIZE(x(:))])
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
    !! In order to make the larger wavelength / the smaller frequency width
    !! be at the longer wavelength side, we define (frequency grid):
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

    dir = '../data/'
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
      !! interpolated with respect to Wavelength
      extCurve(:) = INTERP_LIN_SORTED(Cext0/Cext_V, wave0, wave, &
                                      XLOG=.TRUE., YLOG=.TRUE., FORCE=.TRUE.)
    END IF
      
  END FUNCTION extCurve

  !!---------------------------------------
  !! Total model function for chi2 calling
  !!---------------------------------------

  !! 3D version
  !!============
  FUNCTION specModel_3D( wvl, parval, indpar, Qabs, extinct, verbose, &
                         FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, &
                         FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab )

    USE utilities, ONLY: DP, trimeq, trimlr, pring, &!verbatim, &
                         initiate_clock, time_type, timinfo!, ustd
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: reallocate
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: wvl
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: extinct
    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: parval ! 3D (Nx,Ny,*Npar)
    TYPE(indpar_type), INTENT(IN) :: indpar
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL :: verbose

    INTEGER :: Nw, Nx, Ny, Ncont, Nline, Nband, Nextc, Nstar
    INTEGER :: x, y, i, indref
    REAL(DP) :: Tstar
    REAL(DP), DIMENSION(SIZE(wvl)) :: nu
    REAL(DP), DIMENSION(SIZE(parval,1),SIZE(parval,2)) :: lnIref
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: dblarr4d ! (Nx,Ny,Nw,Ncomp)
    REAL(DP), DIMENSION(SIZE(parval,1),SIZE(parval,2),SIZE(wvl)) :: &
      FnuCONT0, FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
    REAL(DP), DIMENSION(SIZE(parval,1),SIZE(parval,2),SIZE(wvl)) :: specModel_3D

    LOGICAL :: printimer
    TYPE(time_type) :: timestr

    CALL INITIATE_CLOCK(timestr)
    
    !! Preliminaries
    !!---------------
    Nw = SIZE(wvl(:))
    Nx = SIZE(parval(:,:,:),1)
    Ny = SIZE(parval(:,:,:),2)
    Ncont = SIZE(Qabs(:))
    Nband = COUNT(indpar%lnRband(:) .NE. -1)
    Tstar = 5.E4_DP ! [K] high enough to stay in Rayleigh-Jeans limit (T>>5000K for 3um)
    Nline = COUNT(indpar%lnRline(:) .NE. -1)
    Nextc = COUNT(indpar%lnAv(:) .NE. -1)
    Nstar = COUNT(indpar%lnFstar(:) .NE. -1)
    nu(:) = MKS%clight / wvl(:) /MKS%micron ! [Hz] in BLACKBODY

    !! labB index of the reference band (default: 1)
    indref = indpar%refB
    !! (LOG) Intensity of the reference band
    lnIref(:,:) = parval(:,:,indpar%lnRband(indref))
    
    !! Print timer
    printimer = .FALSE.
    IF (PRESENT(verbose)) printimer = verbose

    !! Initialization
    specModel_3D(:,:,:) = 0._DP
    
    !! 1. Continuum
    !!--------------
    CALL REALLOCATE(dblarr4d,Nx,Ny,Nw,Ncont)
    !! The actual unit depends on the input spectrum.
    !! Here as an example, suppose the input is in MKS.
    !! FnuCONT [W/m2/sr/Hz] = Movd2 [Msun/pc2] * [kg/Msun] / [m2/pc2] * modifBB [W/sr/Hz/kg]
    FORALL (x=1:Nx,y=1:Ny,i=1:Ncont) &
      dblarr4d(x,y,:,i) = EXP(parval(x,y,indpar%lnMovd2(i))) * MKS%Msun/MKS%pc**2 * &
                          modifBB(wvl(:),EXP(parval(x,y,indpar%lnT(i))),Qabs(i),.TRUE.)
    FnuCONT0(:,:,:) = SUM(dblarr4d(:,:,:,:),DIM=4)
    IF (PRESENT(FnuCONT)) FnuCONT = FnuCONT0(:,:,:)
    IF (PRESENT(FnuCONT_tab)) FnuCONT_tab = dblarr4d(:,:,:,:)

    ! IF (printimer) &
    !   PRINT*, "[specModel] CAL CONT IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! 2. Bands
    !!----------
    CALL REALLOCATE(dblarr4d,Nx,Ny,Nw,Nband)
    !! FnuBAND [W/m2/sr/Hz] = Iband [W/m2/sr] * lorentzBand [Hz-1]
    DO i=1,Nband
      IF (i == indref) THEN
        FORALL (x=1:Nx,y=1:Ny) &
          dblarr4d(x,y,:,i) = EXP( parval(x,y,indpar%lnRband(i)) ) * &
                              lorentzBand( wvl(:), parval(x,y,indpar%Cband(i)), &
                                parval(x,y,indpar%WSband(i)), &
                                parval(x,y,indpar%WLband(i)),.TRUE. )
      ELSE
        FORALL (x=1:Nx,y=1:Ny) &
          dblarr4d(x,y,:,i) = EXP( parval(x,y,indpar%lnRband(i))+lnIref(x,y) ) * &
                              lorentzBand( wvl(:), parval(x,y,indpar%Cband(i)), &
                                parval(x,y,indpar%WSband(i)), &
                                parval(x,y,indpar%WLband(i)),.TRUE. )
      END IF
    END DO
    FnuBAND0(:,:,:) = SUM(dblarr4d(:,:,:,:),DIM=4)
    IF (PRESENT(FnuBAND)) FnuBAND = FnuBAND0(:,:,:)
    IF (PRESENT(FnuBAND_tab)) FnuBAND_tab = dblarr4d(:,:,:,:)

    ! IF (printimer) &
    !   PRINT*, "[specModel] CAL BAND IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! 3. Stellar Continuum
    !!----------------------
    CALL REALLOCATE(dblarr4d,Nx,Ny,Nw,Nstar)
    !! FnuSTAR [W/m2/sr/Hz] = Fstar [W/m2] * BLACKBODY [W/m2/sr/Hz] /
    !!                        stefan [W/m2/K4] / Tstar4 [K4]
    FORALL (x=1:Nx,y=1:Ny,i=1:Nstar) &
      dblarr4d(x,y,:,i) = EXP(parval(x,y,indpar%lnFstar(i))) * &
                          pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
    FnuSTAR0(:,:,:) = SUM(dblarr4d(:,:,:,:),DIM=4)
    IF (PRESENT(FnuSTAR)) FnuSTAR = FnuSTAR0(:,:,:)
    IF (PRESENT(FnuSTAR_tab)) FnuSTAR_tab = dblarr4d(:,:,:,:)

    ! IF (printimer) &
    !   PRINT*, "[specModel] CAL STAR IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! 4. Screen extinction
    !!----------------------
    CALL REALLOCATE(dblarr4d,Nx,Ny,Nw,Nextc)
    !! Pabs [/] = EXP( Av [mag] * extinction [/] )
    FORALL (x=1:Nx,y=1:Ny,i=1:Nextc) &
      dblarr4d(x,y,:,i) = -EXP(parval(x,y,indpar%lnAv(i)))/1.086_DP * extinct(i,:)
    Pabs0(:,:,:) = EXP( SUM(dblarr4d(:,:,:,:),DIM=4) )
    IF (PRESENT(Pabs)) Pabs = Pabs0(:,:,:)
    IF (PRESENT(Pabs_tab)) Pabs_tab = EXP( dblarr4d(:,:,:,:) )
      
    ! IF (printimer) &
    !   PRINT*, "[specModel] CAL PABS IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! 5. Lines
    !!----------
    CALL REALLOCATE(dblarr4d,Nx,Ny,Nw,Nline)
    !! FnuLINE [W/m2/sr/Hz] = Iline [W/m2/sr] * gaussLine [Hz-1]
    FORALL (x=1:Nx,y=1:Ny,i=1:Nline) &
      dblarr4d(x,y,:,i) = EXP( parval(x,y,indpar%lnRline(i))+lnIref(x,y) ) * &
                          gaussLine(wvl(:), parval(x,y,indpar%Cline(i)), &
                            parval(x,y,indpar%Wline(i)),.TRUE.)
    FnuLINE0(:,:,:) = SUM(dblarr4d(:,:,:,:),DIM=4)
    IF (PRESENT(FnuLINE)) FnuLINE = FnuLINE0(:,:,:)
    IF (PRESENT(FnuLINE_tab)) FnuLINE_tab = dblarr4d(:,:,:,:)

    ! IF (printimer) &
    !   PRINT*, "[specModel] CAL LINE IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! Total model
    !!-------------
    specModel_3D(:,:,:) = (FnuCONT0(:,:,:) + FnuBAND0(:,:,:) + FnuSTAR0(:,:,:) + &
                           FnuLINE0(:,:,:)) * Pabs0(:,:,:)

    IF (printimer) &
      PRINT*, "[specModel] EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
    
  END FUNCTION specModel_3D

  !! 2D version
  !!============
  FUNCTION specModel_2D( wvl, parval, indpar, Qabs, extinct, verbose, &
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
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
      FnuCONT0, FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: &
      FnuCONT_tab0, FnuBAND_tab0, FnuSTAR_tab0, Pabs_tab0, FnuLINE_tab0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
    REAL(DP), DIMENSION(SIZE(parval,1),SIZE(wvl)) :: specModel_2D
    INTEGER :: Nx, Npar
    Nx = SIZE(parval(:,:),1)
    Npar = SIZE(parval(:,:),2)
    specModel_2D(:,:) = RESHAPE( specModel_3D(wvl(:), INDPAR=indpar, &
                                   PARVAL=RESHAPE(parval(:,:),[Nx,1,Npar]), &
                                   QABS=Qabs, EXTINCT=extinct, VERBOSE=verbose, &
                                   FNUCONT=FnuCONT0, FNUBAND=FnuBAND0, FNUSTAR=FnuSTAR0, &
                                   PABS=Pabs0, FNULINE=FnuLINE0, &
                                   FNUCONT_TAB=FnuCONT_tab0, FNUBAND_TAB=FnuBAND_tab0, &
                                   FNUSTAR_TAB=FnuSTAR_tab0, PABS_TAB=Pabs_tab0, &
                                   FNULINE_TAB=FnuLINE_tab0), &
                                 [Nx,SIZE(wvl(:))] )
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
  FUNCTION specModel_1D( wvl, parval, indpar, Qabs, extinct, verbose, &
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
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
      FnuCONT0, FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: &
      FnuCONT_tab0, FnuBAND_tab0, FnuSTAR_tab0, Pabs_tab0, FnuLINE_tab0
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
    REAL(DP), DIMENSION(SIZE(wvl)) :: specModel_1D
    INTEGER :: Npar
    Npar = SIZE(parval(:))
    specModel_1D(:) = RESHAPE( specModel_3D(wvl(:), INDPAR=indpar, &
                                 PARVAL=RESHAPE(parval(:),[1,1,Npar]), &
                                 QABS=Qabs, EXTINCT=extinct, VERBOSE=verbose, &
                                 FNUCONT=FnuCONT0, FNUBAND=FnuBAND0, FNUSTAR=FnuSTAR0, &
                                 PABS=Pabs0, FNULINE=FnuLINE0, &
                                 FNUCONT_TAB=FnuCONT_tab0, FNUBAND_TAB=FnuBAND_tab0, &
                                 FNUSTAR_TAB=FnuSTAR_tab0, PABS_TAB=Pabs_tab0, &
                                 FNULINE_TAB=FnuLINE_tab0), &
                               [SIZE(wvl(:))] )
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
    INTEGER :: i, igrid, indref, iparef, j
    REAL(DP) :: Tstar
    REAL(DP), DIMENSION(SIZE(wvl)) :: nu
    LOGICAL :: gridlnMovd2, gridlnT, gridlnRline, gridCline, gridWline
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
    gridlnMovd2 = ( parname(1:7) == 'lnMovd2' )
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
    gridCONT = ( gridlnMovd2 .OR. gridlnT )
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
        specModel_gen(igrid,:) = EXP( -EXP(parvec(igrid)) / 1.086_DP * extinct(j,:) )

      DO i=1,Nextc
        IF (i /= j) &
          extc1D(:) = extc1D(:) * &
                      EXP( -EXP(parval(indpar%lnAv(i))) / 1.086_DP * extinct(i,:) )
      END DO
    ELSE
      DO i=1,Nextc
        extc1D(:) = extc1D(:) * &
                    EXP( -EXP(parval(indpar%lnAv(i))) / 1.086_DP * extinct(i,:) )
      END DO
    END IF sampEXTC
    
    
    !! 1. Continuum
    !!--------------
    sampCONT: IF (gridCONT) THEN
      IF (gridlnMovd2) THEN
        READ(parname(8:9),*) j
        dblarr1d(:) = MKS%Msun/MKS%pc**2 &
                      * MODIFBB( wvl(:),EXP(parval(indpar%lnT(j))),Qabs(j),.TRUE. )
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP(parvec(igrid)) * dblarr1d(:)
      ELSE IF (gridlnT) THEN
        READ(parname(4:5),*) j
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP(parval(indpar%lnMovd2(j))) &
                                   * MKS%Msun/MKS%pc**2 &
                                   * MODIFBB( wvl(:),EXP(parvec(igrid)),Qabs(j),.TRUE. )
      END IF
      
      DO i=1,Ncont
        IF (i /= j) &
          spec1D(:) = spec1D(:) &
                      + EXP(parval(indpar%lnMovd2(i))) * MKS%Msun/MKS%pc**2 &
                      * MODIFBB( wvl(:),EXP(parval(indpar%lnT(i))),Qabs(i),.TRUE. )
      END DO
    ELSE
      DO i=1,Ncont
        spec1D(:) = spec1D(:) &
                    + EXP(parval(indpar%lnMovd2(i))) * MKS%Msun/MKS%pc**2 &
                    * MODIFBB( wvl(:),EXP(parval(indpar%lnT(i))),Qabs(i),.TRUE. )
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
  
END MODULE auxil
