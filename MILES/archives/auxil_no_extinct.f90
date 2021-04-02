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

  USE datable, ONLY: Ncont_max, Nline_max, Nband_max, Npabs_max, Nstar_max, Cband_sig
  USE utilities, ONLY: DP, STRIKE, WARNING
  USE constants, ONLY: 
  USE inout, ONLY: lenpar
  IMPLICIT NONE
  PRIVATE

  !! only in set_indpar
  INTEGER, PARAMETER, PUBLIC :: NparCONT=2, NparLINE=3, NparBAND=4, NparPABS=1, NparSTAR=1
  
  PUBLIC :: set_indpar, read_master, initparam
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
  
  !! [OBSOLETE] repalced by indpar_type
  TYPE, PUBLIC :: par_type
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lnMovd2
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lnT
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lnIline
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Cline
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Wline
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lnIband
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Cband
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WSband
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WLband
    REAL(DP) :: lnAv
    REAL(DP) :: lnFstar
  END TYPE par_type
 
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
    INTEGER, DIMENSION(Ncont_max) :: lnT = -1 ! dust temperature (MBB) LOG[K]
    INTEGER, DIMENSION(Nline_max) :: lnIline = -1 ! line Intensity LOG[W/m2/sr]
    INTEGER, DIMENSION(Nline_max) :: Cline = -1 ! Center [um]
    INTEGER, DIMENSION(Nline_max) :: Wline = -1 ! Width [um]
    INTEGER, DIMENSION(Nband_max) :: lnIband = -1 ! band Intensity LOG[W/m2/sr]
    INTEGER, DIMENSION(Nband_max) :: Cband = -1 ! Center (peak) [um]
    INTEGER, DIMENSION(Nband_max) :: WSband = -1 ! Width Short nu side [um]
    INTEGER, DIMENSION(Nband_max) :: WLband = -1 ! Width Long nu side [um]
    INTEGER, DIMENSION(Npabs_max) :: lnAv = -1 ! Attenuation in the V band LOG[mag]
    !! Total star surface brightness with dilution factor Omega = r/d LOG[W/m2/sr]
    INTEGER, DIMENSION(Nstar_max) :: lnFstar = -1
    INTEGER, DIMENSION(:), ALLOCATABLE :: extra
  END TYPE indpar_type

  TYPE, PUBLIC :: Qabs_type
    !! rho(kg/m3); wave(micron); nu(Hz); kappa(m2/kg)
    REAL(DP)                            :: rho
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave, nu
    REAL(DP), DIMENSION(:), ALLOCATABLE :: kappa
  END TYPE Qabs_type
 
CONTAINS

  !!-------------------------------------------------------
  !!
  !! Fill the INDPAR_TYPE structure, from a PARINFO_TYPE structure
  !!
  !!-------------------------------------------------------
  SUBROUTINE set_indpar( indpar, parinfo )
    
    USE utilities, ONLY: trimeq, trimlr, pring
    USE arrays, ONLY: iwhere
    IMPLICIT NONE

    TYPE(indpar_type), INTENT(OUT) :: indpar
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo

    INTEGER :: i, Ncont, Nline, Nband, Npabs, Nstar, Nextra
    LOGICAL :: CONTset, LINEset, BANDset, PABSset, STARset, extraset

    CONTset = ( ANY(TRIMEQ(parinfo(:)%comp, 'CONT')) )
    LINEset = ( ANY(TRIMEQ(parinfo(:)%comp, 'LINE')) )
    BANDset = ( ANY(TRIMEQ(parinfo(:)%comp, 'BAND')) )
    PABSset = ( ANY(TRIMEQ(parinfo(:)%comp, 'PABS')) )
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
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'lnIline'//TRIMLR(PRING(i))), indpar%lnIline(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'Cline'//TRIMLR(PRING(i))), indpar%Cline(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'Wline'//TRIMLR(PRING(i))), indpar%Wline(i) )

      END DO
    END IF

    IF (BANDset) THEN
      Nband = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'BAND')) ) / NparBAND
      DO i=1,Nband
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'lnIband'//TRIMLR(PRING(i))), indpar%lnIband(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'Cband'//TRIMLR(PRING(i))), indpar%Cband(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'WSband'//TRIMLR(PRING(i))), indpar%WSband(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), 'WLband'//TRIMLR(PRING(i))), indpar%WLband(i) )
        
      END DO
    END IF

    IF (PABSset) THEN
      Npabs = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'PABS')) ) / NparPABS
      DO i=1,Npabs
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
    
  END SUBROUTINE set_indpar

  !!-------------------------------------------------------
  !!
  !! Read the input master file for the Chi2/Bayesian run
  !!
  !!-------------------------------------------------------

  SUBROUTINE read_master( wavAll, dirIN, dirOUT, spec_unit, &
                          Nmcmc, verbose, NiniMC, &!robust_RMS, robust_cal, skew_RMS, &
                          calib, newseed, newinit, dostop, nohi, &
                          labQ, labL, labB, Qabs, extinct, &
                          Ncont, Nband, Nline, Npabs, Nstar, Nextra, &
                          ! corrhypname, corrname, &
                          parinfo, parmodinfo, parhypinfo, parextinfo, &
                          indpar, Npar, Nparmod, Nparhyp, Ncorrhyp, Ncorr)

    USE utilities, ONLY: DP, strike, trimeq, trimlr, pring!, verbatim
    USE inout, ONLY: lenpar, lenpath, read_input_line, read_hdf5, h5ext
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: incrarr, iwhere, reallocate, closest
    USE interpolation, ONLY: interp_lin_sorted
    USE statistics, ONLY: correl_parlist, N_corr
    USE grain_optics, ONLY: lendustQ, rho_grain, read_optics
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN), OPTIONAL :: dirIN ! default: ./
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: wavAll ! Interpolate Qabs if wavAll presents

    INTEGER, PARAMETER :: unit = 1

    CHARACTER(lenpath) :: filmas, filmod, filext
    ! CHARACTER(lenpath) :: dirOUT0
    ! CHARACTER(lenpar) :: spec_unit0
    CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ0, labL0, labB0
    ! INTEGER :: Nmcmc0, NiniMC0
    INTEGER :: Ncont0, Nband0, Nline0, Npabs0, Nstar0, Nextra0, Npar0
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
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: radius ! (Nr, Nw)
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Qabsall ! (Ncont, Nr, Nw)
    
    CHARACTER(lenpath), INTENT(OUT), OPTIONAL :: dirOUT
    CHARACTER(lenpar), INTENT(OUT), OPTIONAL :: spec_unit
    CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: labQ, labB, labL
    ! CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
    !   corrhypname, corrname
    INTEGER, INTENT(OUT), OPTIONAL :: Nmcmc, NiniMC, Ncont, Nband, Nline, Npabs, Nstar, Nextra
    INTEGER, INTENT(OUT), OPTIONAL :: Npar, Nparmod, Nparhyp, Ncorrhyp, Ncorr
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: extinct
    TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Qabs
    TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      parinfo, parmodinfo, parhypinfo, parextinfo
    TYPE(indpar_type), INTENT(OUT), OPTIONAL :: indpar
    LOGICAL, INTENT(OUT), OPTIONAL :: verbose!, robust_RMS, robust_cal, skew_RMS
    LOGICAL, INTENT(OUT), OPTIONAL :: calib, newseed, newinit, dostop, nohi
    
    !! Read the input master file
    !!----------------------------
    IF (PRESENT(dirIN)) THEN
      filmas = TRIMLR(dirIN)//'input_fitMIR_master'//h5ext
      filmod = TRIMLR(dirIN)//'input_fitMIR_model'//h5ext
      filext = TRIMLR(dirIN)//'input_fitMIR_extra'//h5ext
    ELSE
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
    ALLOCATE(intarr1d(1))
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

    !! Get Npabs
    Npabs0 = 1

    !! Get Nstar
    Nstar0 = 1
    
    IF (PRESENT(Ncont)) Ncont = Ncont0
    IF (PRESENT(Nband)) Nband = Nband0
    IF (PRESENT(Nline)) Nline = Nline0
    IF (PRESENT(Npabs)) Npabs = Npabs0
    IF (PRESENT(Nstar)) Nstar = Nstar0

    !! Read optical properties (Qabs)
    !!--------------------------------
    IF (PRESENT(Qabs)) THEN
      !! Qabs % rho(kg/m3); wave(micron); nu(Hz); kappa(m2/kg)
      CALL READ_OPTICS(labQ0, WAVE=wave0, RADIUSALL=radius, QABSALL=Qabsall)    
      
      Nw = SIZE(wave0(:))

      ALLOCATE(indrad(Ncont0), Qabs(Ncont0))
      
      DO i=1,Ncont0
        
        ALLOCATE(Qabs(i)%kappa(Nw))
        
        Nr = SIZE(radius(i,:))
        Qabs(i)%rho = rho_grain(labQ0(i))
        Qabs(i)%wave = wave0(:)
        Qabs(i)%nu = MKS%clight / wave0(:) /MKS%micron
        indrad(i) = CLOSEST(radius(i,:), a0)
        ! PRINT*, 'Radius of ', labQ0(i), ': ', radius(i,indrad(i)), '->', indrad(i)
        Qabs(i)%kappa = 3._DP/4._DP/rho_grain(labQ0(i)) * Qabsall(i,indrad(i),:)/radius(i,indrad(i))
        
        !! [Optional] Interpolate Qabs to input wave grid
        IF (PRESENT(wavAll)) THEN
          Qabs(i)%wave = wavAll(:)
          Qabs(i)%nu = MKS%clight / wavAll(:) /MKS%micron
          Qabs(i)%kappa = interp_lin_sorted(Qabs(i)%kappa, wave0(:), wavAll(:), &
                                            XLOG=.TRUE., YLOG=.TRUE., FORCE=.TRUE.)
          
        END IF
      END DO
      !! Free memory space
      DEALLOCATE(wave0, radius, Qabsall)
      
    END IF

    !! Read extinction curve
    !!-----------------------
    IF (PRESENT(extinct)) THEN
      IF (PRESENT(wavAll)) THEN
        extinct = extCurve(wavAll(:))
      ELSE
        CALL STRIKE('READ_MASTER', 'Please input wavelength grid. ')

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
    Nparall = 2*Ncont0 + 4*Nband0 + 3*Nline0 + Npabs0 + Nstar0 + Nextra0

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
    IF (PRESENT(indpar)) CALL SET_INDPAR(indpar, parinfo0(:))

    !! Free memory space
    !!-------------------
    DEALLOCATE (labQ0, labL0, labB0, boolpar, boolparmod, parinfo0, parinfall)
    
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
    
    INTEGER :: i, j, ipar, iw, x, y, Nx, Ny, Nw, Nmc, Nwc, itab
    INTEGER :: Npar, Ncont, Nline, Nband, Npabs, Nstar, Nextra
    INTEGER, DIMENSION(:), ALLOCATABLE :: indwi, indw, indwc
    INTEGER, DIMENSION(SIZE(par,1),SIZE(par,2)) :: iw2
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: iniparval, FnuOBS, val3
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mu
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wOBS, nuOBS
    REAL(DP), DIMENSION(SIZE(par,1),SIZE(par,2),MAX(NiniMC,1)) :: theta3
    REAL(DP), DIMENSION(SIZE(par,1),SIZE(par,2)) :: limi2, lims2, val2
    REAL(DP) :: limi, lims, wcen, dw, val, sig, theta, Tstar
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
    Npabs = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'PABS')) ) / NparPABS
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
    ! FnuOBS(:,:,:) = FnuOBS(:,:,:) * MKS%Jy*1.E6 ! convert MJy/sr to W/m2/Hz/sr (MKS)

    newinipar: IF (.NOT. newini) THEN
      
      !!----------------------------------------------------
      !! I. Automatically generate initial parameter values
      !!----------------------------------------------------

      IF (NiniMC==0) THEN

        !! a. Simple initilization
        !!-------------------------

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
            parinfo(ind%Cline(i))%limited(1) = .TRUE.
            IF (.NOT. parinfo(ind%Cline(i))%limited(2)) &
              parinfo(ind%Cline(i))%limits(2) = val+sig
            parinfo(ind%Cline(i))%limited(2) = .TRUE.
          END IF
          
          !! Wline
          IF (parinfo(ind%Wline(i))%fixed) THEN
            par(:,:,ind%Wline(i),:) = parinfo(ind%Wline(i))%value
          ELSE
            val = degradeRes(TABLine(itab)%wave, .01_DP, 'SL-LL')
            par(:,:,ind%Wline(i),:) = val
            IF (.NOT. parinfo(ind%Wline(i))%limited(1)) &
              parinfo(ind%Wline(i))%limits(1) = val*0.5_DP
            parinfo(ind%Wline(i))%limited(1) = .TRUE.
            IF (.NOT. parinfo(ind%Wline(i))%limited(2)) &
              parinfo(ind%Wline(i))%limits(2) = val*2._DP
            parinfo(ind%Wline(i))%limited(2) = .TRUE.
          END IF
          
          !! lnIline
          IF (parinfo(ind%lnIline(i))%fixed) THEN
            par(:,:,ind%lnIline(i),:) = parinfo(ind%lnIline(i))%value
          ELSE
            FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
              iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cline(i),1) )
              !! lnIline is auto determined by the spectrum to fit (idem. for lnIband)
              !! Can be overestimated (should subtract other compo from FnuOBS)
              !! or underestimated (due to extinction Pabs).
              !! Here we just need a reasonable order of magnitude.
              par(x,y,ind%lnIline(i),:) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) )
            END FORALL
          END IF
          
          !! Find indw
          !! Take (1,1) Cline and Wline since they are supposed to vary little (idem. below)
          limi = par(1,1,ind%Cline(i),1) - par(1,1,ind%Wline(i),1)
          lims = par(1,1,ind%Cline(i),1) + par(1,1,ind%Wline(i),1)
          IF (ANY(wOBS(:)>limi .AND. wOBS(:)<lims)) THEN
            CALL IWHERE( wOBS(:)>limi .AND. wOBS(:)<lims, indwi )
            DO iw=1,SIZE(indwi)
              CALL INCRARR( indw, indwi(iw) )
            END DO
          END IF
          
        END DO

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
            parinfo(ind%Cband(i))%limited(1) = .TRUE.
            IF (.NOT. parinfo(ind%Cband(i))%limited(2)) &
              parinfo(ind%Cband(i))%limits(2) = val+sig
            parinfo(ind%Cband(i))%limited(2) = .TRUE.
          END IF
          
          !! WSband
          IF (parinfo(ind%WSband(i))%fixed) THEN
            par(:,:,ind%WSband(i),:) = parinfo(ind%WSband(i))%value
          ELSE
            val = TABand(itab)%sigmaS
            par(:,:,ind%WSband(i),:) = val
            IF (.NOT. parinfo(ind%WSband(i))%limited(1)) &
              parinfo(ind%WSband(i))%limits(1) = val*0.5_DP
            parinfo(ind%WSband(i))%limited(1) = .TRUE.
            IF (.NOT. parinfo(ind%WSband(i))%limited(2)) &
              parinfo(ind%WSband(i))%limits(2) = val*2._DP
            parinfo(ind%WSband(i))%limited(2) = .TRUE.
          END IF
   
          !! WLband
          IF (parinfo(ind%WLband(i))%fixed) THEN
            par(:,:,ind%WLband(i),:) = parinfo(ind%WLband(i))%value
          ELSE
            val = TABand(itab)%sigmaL
            par(:,:,ind%WLband(i),:) = val
            IF (.NOT. parinfo(ind%WLband(i))%limited(1)) &
              parinfo(ind%WLband(i))%limits(1) = val*0.5_DP
            parinfo(ind%WLband(i))%limited(1) = .TRUE.
            IF (.NOT. parinfo(ind%WLband(i))%limited(2)) &
              parinfo(ind%WLband(i))%limits(2) = val*2._DP
            parinfo(ind%WLband(i))%limited(2) = .TRUE.
          END IF
          
          !! lnIband
          IF (parinfo(ind%lnIband(i))%fixed) THEN
            par(:,:,ind%lnIband(i),:) = parinfo(ind%lnIband(i))%value
          ELSE
            FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
              iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cband(i),1) )
              !! Same estimation as lnIline
              par(x,y,ind%lnIband(i),:) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) )
            END FORALL
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
              val3(x,y,:) = modifBB(wOBS(:), EXP(par(x,y,ind%lnT(i),1)), Qabs(i)) * MKS%Msun/MKS%pc**2
              par(x,y,ind%lnMovd2(i),:) = LOG( ABS(FnuOBS(x,y,indwc(iw)))/val3(x,y,indwc(iw)) )
            END FORALL
          END IF

        END DO

        DO i=1,Npabs
          !! lnAv LOG[mag]
          IF (parinfo(ind%lnAv(i))%fixed) THEN
            par(:,:,ind%lnAv(i),:) = parinfo(ind%lnAv(i))%value
          ELSE
            par(:,:,ind%lnAv(i),:) = 0. ! 1 [mag] = no extinction
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
          END IF
          
          !! lnIline
          IF (parinfo(ind%lnIline(i))%fixed) THEN
            par(:,:,ind%lnIline(i),:) = parinfo(ind%lnIline(i))%value
          ELSE
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
                  iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cline(i),j) )
                  !! Estimated via feature peak
                  par(x,y,ind%lnIline(i),j) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) )

                END FORALL
              ELSE
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
                  iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cline(i),j) )
                  !! Estimated via feature peak
                  val2(x,y) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) )
                  limi2(x,y) = MERGE( parinfo(ind%lnIline(i))%limits(1), &
                                      val2(x,y)-.1_DP, & ! EXP(factor - 3)
                                      parinfo(ind%lnIline(i))%limited(1)) ! lim inf
                  lims2(x,y) = MERGE( parinfo(ind%lnIline(i))%limits(2), &
                                      val2(x,y), &
                                      parinfo(ind%lnIline(i))%limited(2)) ! lim sup
                  !! Uniform distribution
                  par(x,y,ind%lnIline(i),j) = (lims2(x,y) - limi2(x,y)) * theta3(x,y,j) + limi2(x,y)
                  
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
          END IF

          !! lnIband
          IF (parinfo(ind%lnIband(i))%fixed) THEN
            par(:,:,ind%lnIband(i),:) = parinfo(ind%lnIband(i))%value
          ELSE
            CALL RANDOM_NUMBER(theta3(:,:,:))
            DO j=1,NiniMC
              IF (j == 1) THEN
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
                  iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cband(i),j) )
                  !! Estimated via feature peak
                  par(x,y,ind%lnIband(i),j) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) )

                END FORALL
              ELSE
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
                  iw2(x,y) = CLOSEST( wOBS(:), par(x,y,ind%Cband(i),j) )
                  !! Estimated via feature peak
                  val2(x,y) = LOG( ABS(FnuOBS(x,y,iw2(x,y))) )
                  limi2(x,y) = MERGE( parinfo(ind%lnIband(i))%limits(1), &
                                      val2(x,y)-.1_DP, & ! EXP(factor - 3)
                                      parinfo(ind%lnIband(i))%limited(1) ) ! lim inf
                  lims2(x,y) = MERGE( parinfo(ind%lnIband(i))%limits(2), &
                                      val2(x,y), &
                                      parinfo(ind%lnIband(i))%limited(2) ) ! lim sup
                  !! Uniform distribution
                  par(x,y,ind%lnIband(i),j) = (lims2(x,y) - limi2(x,y)) * theta3(x,y,j) + limi2(x,y)
                  
                END FORALL
              END IF
            END DO
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
                  val3(x,y,:) = modifBB(wOBS(:), EXP(par(x,y,ind%lnT(i),1)), Qabs(i)) &
                                * MKS%Msun/MKS%pc**2
                  val2(x,y) = LOG( ABS(FnuOBS(x,y,indwc(iw)))/val3(x,y,indwc(iw)) )
                  par(x,y,ind%lnMovd2(i),j) = val2(x,y)

                END FORALL
              ELSE
                FORALL (x=1:Nx,y=1:Ny,maskxy(x,y))
                  val3(x,y,:) = modifBB(wOBS(:), EXP(par(x,y,ind%lnT(i),1)), Qabs(i)) &
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
        
        DO i=1,Npabs
          !! lnAv LOG[mag]
          IF (parinfo(ind%lnAv(i))%fixed) THEN
            par(:,:,ind%lnAv(i),:) = parinfo(ind%lnAv(i))%value
          ELSE
            !! lnAv = 0 -> Av = 1 [mag] = no extinction
            val = 0.
            limi = MERGE( parinfo(ind%lnAv(i))%limits(1), &
                          val-.1_DP, & ! EXP(factor - 3)
                          parinfo(ind%lnAv(i))%limited(1) ) ! lim inf
            lims = MERGE( parinfo(ind%lnAv(i))%limits(2), &
                          val+.1_DP, &
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
                  par(x,y,ind%lnFstar(i),j) = (lims2(x,y) - limi2(x,y)) * theta3(x,y,j) + limi2(x,y)
    
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
          parinfo(ind%Cline(i))%limited(1) = .TRUE.
          IF (.NOT. parinfo(ind%Cline(i))%limited(2)) &
            parinfo(ind%Cline(i))%limits(2) = wcen+sig
          parinfo(ind%Cline(i))%limited(2) = .TRUE.
          
          !! Wline
          par(:,:,ind%Wline(i),:) = dw
          IF (.NOT. parinfo(ind%Wline(i))%limited(1)) &
            parinfo(ind%Wline(i))%limits(1) = dw*0.5_DP
          parinfo(ind%Wline(i))%limited(1) = .TRUE.
          IF (.NOT. parinfo(ind%Wline(i))%limited(2)) &
            parinfo(ind%Wline(i))%limits(2) = dw*2._DP
          parinfo(ind%Wline(i))%limited(2) = .TRUE.
          
        END DO

        DO i=1,Nband
          CALL IWHERE( TRIMEQ(TABand(:)%label, labB(i)), itab ) ! itab = ind of TABand
          wcen = TABand(itab)%wave
          sig = Cband_sig
          
          !! Cband
          par(:,:,ind%Cband(i),:) = wcen
          IF (.NOT. parinfo(ind%Cband(i))%limited(1)) &
            parinfo(ind%Cband(i))%limits(1) = wcen-sig
          parinfo(ind%Cband(i))%limited(1) = .TRUE.
          IF (.NOT. parinfo(ind%Cband(i))%limited(2)) &
            parinfo(ind%Cband(i))%limits(2) = wcen+sig
          parinfo(ind%Cband(i))%limited(2) = .TRUE.
          
          !! WSband
          val = TABand(itab)%sigmaS
          par(:,:,ind%WSband(i),:) = val
          IF (.NOT. parinfo(ind%WSband(i))%limited(1)) &
            parinfo(ind%WSband(i))%limits(1) = val*0.5_DP
          parinfo(ind%WSband(i))%limited(1) = .TRUE.
          IF (.NOT. parinfo(ind%WSband(i))%limited(2)) &
            parinfo(ind%WSband(i))%limits(2) = val*2._DP
          parinfo(ind%WSband(i))%limited(2) = .TRUE.
   
          !! WLband
          val = TABand(itab)%sigmaL
          par(:,:,ind%WLband(i),:) = val
          IF (.NOT. parinfo(ind%WLband(i))%limited(1)) &
            parinfo(ind%WLband(i))%limits(1) = val*0.5_DP
          parinfo(ind%WLband(i))%limited(1) = .TRUE.
          IF (.NOT. parinfo(ind%WLband(i))%limited(2)) &
            parinfo(ind%WLband(i))%limits(2) = val*2._DP
          parinfo(ind%WLband(i))%limited(2) = .TRUE.

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
  !! Dust contimuum (N BB)
  !!-----------------------
  PURE FUNCTION modifBB_gen( dblarr, temp, Qabs )
    !! temp [K]; BLACKBODY [W/m2/sr/Hz]; 
    !! output array [W/sr/Hz/kg]
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:), INTENT(IN)           :: dblarr
    REAL(DP), DIMENSION(:), INTENT(IN)           :: temp
    TYPE(Qabs_type), INTENT(IN), DIMENSION(:)    :: Qabs
    INTEGER :: Nt, it
    REAL(DP), DIMENSION(SIZE(dblarr))            :: nu
    REAL(DP), DIMENSION(SIZE(temp),SIZE(dblarr)) :: modifBB_gen

    Nt = SIZE(temp)
    nu(:) = MKS%clight / dblarr(:) /MKS%micron ! should be the same with Qabs%nu

    !! modifBB [W/sr/Hz/kg] = kappa [m2/kg] * BLACKBODY [W/m2/sr/Hz]
    FORALL (it=1:Nt) &
      modifBB_gen(it,:) = Qabs(it)%kappa * BLACKBODY(nu(:), temp(it))
    
  END FUNCTION modifBB_gen

  PURE FUNCTION modifBB_0(dblarr, temp, Qabs)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN)           :: dblarr
    REAL(DP), INTENT(IN)                         :: temp
    TYPE(Qabs_type), INTENT(IN)                  :: Qabs
    REAL(DP), DIMENSION(SIZE(dblarr))            :: modifBB_0
    modifBB_0(:) = RESHAPE(modifBB_gen(dblarr(:), [temp], [Qabs]), [SIZE(dblarr(:))])
  END FUNCTION modifBB_0

  PURE FUNCTION modifBB_1(dblarr, temp, Qabs)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN)           :: dblarr
    REAL(DP), DIMENSION(:), INTENT(IN)           :: temp
    TYPE(Qabs_type), INTENT(IN), DIMENSION(:)    :: Qabs
    REAL(DP), DIMENSION(SIZE(temp),SIZE(dblarr)) :: modifBB_1
    modifBB_1(:,:) = RESHAPE(modifBB_gen(dblarr(:), temp(:), Qabs(:)), &
                             [SIZE(temp(:)), SIZE(dblarr(:))])
  END FUNCTION modifBB_1

  !!-----------------------------------------------------
  !! Atomic & molecular unresolved lines (Gauss profile)
  !!-----------------------------------------------------
  PURE FUNCTION gaussLine_gen( dblarr, ref, sig, w2nu )
    !! w2nu=.TRUE. if input is in wavelength
    !! Output array [Hz^-1]
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE distributions, ONLY: dist_gauss
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN)          :: dblarr
    REAL(DP), DIMENSION(:), INTENT(IN)          :: ref, sig
    LOGICAL, INTENT(IN), OPTIONAL               :: w2nu
    INTEGER :: N, i
    REAL(DP), DIMENSION(SIZE(dblarr))           :: nuIN
    REAL(DP), DIMENSION(SIZE(ref))              :: nuref, nusig
    REAL(DP), DIMENSION(SIZE(ref),SIZE(dblarr)) :: gaussLine_gen

    N = SIZE(ref)
    nuIN(:) = dblarr(:)
    nuref(:) = ref(:)
    nusig(:) = sig(:)
    IF (PRESENT(w2nu)) THEN
      IF (w2nu) THEN
        nuIN(:) = MKS%clight / dblarr(:) /MKS%micron
        nuref(:) = MKS%clight / ref(:) /MKS%micron
        nusig(:) = MKS%clight * (1./(ref(:)-sig(:)) - 1./(ref(:)+sig(:))) / 2._DP /MKS%micron

      END IF
    END IF
    FORALL (i=1:N) &
      gaussLine_gen(i,:) = DIST_GAUSS(nuIN(:), nuref(i), nusig(i))

  END FUNCTION gaussLine_gen

  PURE FUNCTION gaussLine_0( dblarr, ref, sig, w2nu )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: dblarr
    REAL(DP), INTENT(IN)               :: ref, sig
    LOGICAL, INTENT(IN), OPTIONAL      :: w2nu
    REAL(DP), DIMENSION(SIZE(dblarr))  :: gaussLine_0
    IF (PRESENT(w2nu)) THEN
      gaussLine_0(:) = RESHAPE(gaussLine_gen(dblarr(:), [ref], [sig], w2nu), &
                               [SIZE(dblarr(:))])
    ELSE
      gaussLine_0(:) = RESHAPE(gaussLine_gen(dblarr(:), [ref], [sig]), [SIZE(dblarr(:))])
    END IF
  END FUNCTION gaussLine_0


  PURE FUNCTION gaussLine_1( dblarr, ref, sig, w2nu )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN)          :: dblarr
    REAL(DP), DIMENSION(:), INTENT(IN)          :: ref, sig
    LOGICAL, INTENT(IN), OPTIONAL               :: w2nu
    REAL(DP), DIMENSION(SIZE(ref),SIZE(dblarr)) :: gaussLine_1
    IF (PRESENT(w2nu)) THEN
      gaussLine_1(:,:) = RESHAPE(gaussLine_gen(dblarr(:), ref(:), sig(:), w2nu), &
                                 [SIZE(ref(:)), SIZE(dblarr(:))])
    ELSE
      gaussLine_1(:,:) = RESHAPE(gaussLine_gen(dblarr(:), ref(:), sig(:)), &
                                 [SIZE(ref(:)), SIZE(dblarr(:))])
    END IF
  END FUNCTION gaussLine_1
  
  !!------------------------------------------------------
  !! Resolved aromatic bands (Asymmetric Lorentz profile)
  !!------------------------------------------------------
  PURE FUNCTION lorentzBand_gen( dblarr, ref, sigL, sigS, w2nu )
    !! w2nu=.TRUE. if input is in wavelength
    !! Output array [Hz^-1]
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE distributions, ONLY: dist_splitlorentz
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN)          :: dblarr
    REAL(DP), DIMENSION(:), INTENT(IN)          :: ref, sigS, sigL
    LOGICAL, INTENT(IN), OPTIONAL               :: w2nu
    INTEGER :: N, i
    REAL(DP), DIMENSION(SIZE(dblarr))           :: nuIN
    REAL(DP), DIMENSION(SIZE(ref))              :: nuref, nusigS, nusigL
    REAL(DP), DIMENSION(SIZE(ref))              :: lambda, tau
    REAL(DP), DIMENSION(SIZE(ref),SIZE(dblarr)) :: lorentzBand_gen

    N = SIZE(ref)
    
    nusigS = 1._DP
    nusigL = 1._DP
    IF (PRESENT(w2nu)) THEN
      IF (w2nu) THEN
        nuIN(:) = MKS%clight / dblarr(:) /MKS%micron
        nuref(:) = MKS%clight / ref(:) /MKS%micron
        nusigS(:) = MKS%clight * (1./(ref(:)-sigS(:)) - 1./ref(:)) /MKS%micron
        nusigL(:) = MKS%clight * (1./ref(:) - 1./(ref(:)+sigL(:))) /MKS%micron
      ELSE
        nuIN(:) = dblarr(:)
        nuref(:) = ref(:)
        nusigS(:) = sigS(:)
        nusigL(:) = sigL(:)
        
      END IF
    ELSE
      nuIN(:) = dblarr(:)
      nuref(:) = ref(:)
      nusigS(:) = sigS(:)
      nusigL(:) = sigL(:)
      
    END IF
    !! The SwING library adopted the arbitary parameters such that 
    !! dnu_short corresponds to longer wavelength side of the curve,
    !! Note that neither dnu_S nor dnu_L is std_dev of a lorentzian,
    !! 'cause it is undefined (infinite).
    !! This is analogue of (split) normal distribution.
    !! In addition, we denote lambda=2*dnu_long (analogue of FWHM)
    !! and tau=dnu_short/dnu_long (shape param).
    lambda(:) = 2*nusigL(:) ! nusigL = lambda/2 (lambda -> width param)
    tau(:) = nusigS(:)/nusigL(:) ! nusigS = tau * nusigL (tau -> shape param)
    !! norm = lambda / (1+tau) / pi
    FORALL (i=1:N) &
      lorentzBand_gen(i,:) = DIST_SPLITLORENTZ(nuIN(:), nuref(i), lambda(i), tau(i))

  END FUNCTION lorentzBand_gen

  PURE FUNCTION lorentzBand_0( dblarr, ref, sigS, sigL, w2nu )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: dblarr
    REAL(DP), INTENT(IN)               :: ref, sigS, sigL
    LOGICAL, INTENT(IN), OPTIONAL      :: w2nu
    REAL(DP), DIMENSION(SIZE(dblarr))  :: lorentzBand_0
    IF (PRESENT(w2nu)) THEN
      lorentzBand_0(:) = RESHAPE(lorentzBand_gen(dblarr(:), [ref], [sigS], [sigL], w2nu), &
                                 [SIZE(dblarr(:))])
    ELSE
      lorentzBand_0(:) = RESHAPE(lorentzBand_gen(dblarr(:), [ref], [sigS], [sigL]), &
                                 [SIZE(dblarr(:))])
    END IF
  END FUNCTION lorentzBand_0

  PURE FUNCTION lorentzBand_1( dblarr, ref, sigS, sigL, w2nu )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN)          :: dblarr
    REAL(DP), DIMENSION(:), INTENT(IN)          :: ref, sigS, sigL
    LOGICAL, INTENT(IN), OPTIONAL               :: w2nu
    REAL(DP), DIMENSION(SIZE(ref),SIZE(dblarr)) :: lorentzBand_1
    IF (PRESENT(w2nu)) THEN
      lorentzBand_1(:,:) = RESHAPE(lorentzBand_gen(dblarr(:), ref(:), sigS(:), sigL(:), w2nu), &
                                   [SIZE(ref(:)), SIZE(dblarr(:))])
    ELSE
      lorentzBand_1(:,:) = RESHAPE(lorentzBand_gen(dblarr(:), ref(:), sigS(:), sigL(:)), &
                                   [SIZE(ref(:)), SIZE(dblarr(:))])
    END IF
  END FUNCTION lorentzBand_1

  !!------------------
  !! Extinction curve
  !!------------------
  FUNCTION extCurve( dblarr, nu2w )
    !! https://www.astro.princeton.edu/~draine/dust/dustmix.html
    !! Diffuse (MW) Rv = 3.1
    !! nu2w=.TRUE. if input is in frequency [Hz]
    USE utilities, ONLY: DP
    USE constants, ONLY: MKS
    USE arrays, ONLY: closest
    USE inout, ONLY: read_hdf5, h5ext
    USE interpolation, ONLY: interp_lin_sorted
    IMPLICIT NONE

    LOGICAL, INTENT(IN), OPTIONAL       :: nu2w
    REAL(DP), DIMENSION(:), INTENT(IN)  :: dblarr
    INTEGER                             :: Nw0
    REAL(DP)                            :: wave_V, Cext_V
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave, Cext
    REAL(DP), DIMENSION(SIZE(dblarr))   :: waveIN
    REAL(DP), DIMENSION(SIZE(dblarr))   :: extCurve

    IF (PRESENT(nu2w)) THEN
      IF (nu2w) THEN
        waveIN = MKS%clight / dblarr /MKS%micron ! input [Hz]
      END IF
    ELSE
      waveIN = dblarr
      
    END IF
    
    CALL READ_HDF5(DBLARR1D=wave, NAME='lambda (micron)', &
                   N1=Nw0, FILE='../data/extcurve'//h5ext)
    CALL READ_HDF5(DBLARR1D=Cext, NAME='C_extovH (cm^2ovH)', &
                   FILE='../data/extcurve'//h5ext)
    wave_V = 0.5470_DP
    Cext_V = Cext(closest(wave, wave_V))
    !! interpolated with respect to Wavelength
    extCurve(:) = interp_lin_sorted(Cext/Cext_V, wave, waveIN, &
                                 XLOG=.TRUE., YLOG=.TRUE., FORCE=.TRUE.)

    DEALLOCATE(wave, Cext)

  END FUNCTION extCurve

  !!---------------------------------------
  !! Total model function for chi2 calling
  !!---------------------------------------

  !! 3D version
  !!============
  FUNCTION specModel_3D( wvl, indpar, parval, Qabs, extinct, verbose, &
                         FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, &
                         FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab )

    USE utilities, ONLY: DP, trimeq, trimlr, pring, &!verbatim, &
                         initiate_clock, time_type, timinfo!, ustd
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: reallocate
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: wvl, extinct
    TYPE(indpar_type), INTENT(IN) :: indpar
    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: parval ! 3D (Nx,Ny,*Npar)
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL :: verbose

    INTEGER :: Nw, Nx, Ny, Ncont, Nline, Nband, Npabs, Nstar
    INTEGER :: x, y, i
    REAL(DP) :: Tstar
    REAL(DP), DIMENSION(SIZE(wvl)) :: nu
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: Const4D ! (Nx,Ny,Nw,Ncomp)
    REAL(DP), DIMENSION(SIZE(parval,1),SIZE(parval,2),SIZE(wvl)) :: FnuCONT0, &
      FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuCONT, &
      FnuBAND, FnuSTAR, Pabs, FnuLINE
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuCONT_tab, &
      FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
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
    Nband = COUNT(indpar%lnIband(:) .NE. -1)
    Tstar = 5.E4_DP ! [K] high enough to stay in Rayleigh-Jeans limit (T>>5000K for 3um)
    Nline = COUNT(indpar%lnIline(:) .NE. -1)
    Npabs = COUNT(indpar%lnAv(:) .NE. -1)
    Nstar = COUNT(indpar%lnFstar(:) .NE. -1)
    nu(:) = MKS%clight / wvl(:) /MKS%micron ! [Hz] in BLACKBODY
    
    !! Print timer
    printimer = .FALSE.
    IF (PRESENT(verbose)) printimer = verbose

    !! Initialization
    specModel_3D(:,:,:) = 0._DP
    
    !! 1. Continuum
    !!--------------
    CALL REALLOCATE(Const4D,Nx,Ny,Nw,Ncont)
    !! The actual unit depends on the input spectrum.
    !! Here as an example, suppose the input is in MKS.
    !! FnuCONT [W/m2/sr/Hz] = Movd2 [Msun/pc2] * [kg/Msun] / [m2/pc2] * modifBB [W/sr/Hz/kg]
    FORALL (x=1:Nx,y=1:Ny,i=1:Ncont) &
      Const4D(x,y,:,i) = EXP(parval(x,y,indpar%lnMovd2(i))) * MKS%Msun/MKS%pc**2 * &
                         modifBB(wvl(:), EXP(parval(x,y,indpar%lnT(i))), Qabs(i))
    FnuCONT0(:,:,:) = SUM(Const4D(:,:,:,:),DIM=4)
    IF (PRESENT(FnuCONT)) FnuCONT = FnuCONT0(:,:,:)
    IF (PRESENT(FnuCONT_tab)) FnuCONT_tab = Const4D(:,:,:,:)

    ! IF (printimer) &
    !   PRINT*, "[specModel] CAL CONT IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! 2. Bands
    !!----------
    CALL REALLOCATE(Const4D,Nx,Ny,Nw,Nband)
    !! FnuBAND [W/m2/sr/Hz] = Iband [W/m2/sr] * lorentzBand [Hz-1]
    FORALL (x=1:Nx,y=1:Ny,i=1:Nband) &
      Const4D(x,y,:,i) = EXP(parval(x,y,indpar%lnIband(i))) * &
                         lorentzBand(wvl(:), parval(x,y,indpar%Cband(i)), &
                           parval(x,y,indpar%WSband(i)), parval(x,y,indpar%WLband(i)))
    FnuBAND0(:,:,:) = SUM(Const4D(:,:,:,:),DIM=4)
    IF (PRESENT(FnuBAND)) FnuBAND = FnuBAND0(:,:,:)
    IF (PRESENT(FnuBAND_tab)) FnuBAND_tab = Const4D(:,:,:,:)

    ! IF (printimer) &
    !   PRINT*, "[specModel] CAL BAND IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! 3. Stellar Continuum
    !!----------------------
    CALL REALLOCATE(Const4D,Nx,Ny,Nw,Nstar)
    !! FnuSTAR [W/m2/sr/Hz] = Fstar [W/m2] * BLACKBODY [W/m2/sr/Hz] / stefan [W/m2/K4] / Tstar4 [K4]
    FORALL (x=1:Nx,y=1:Ny,i=1:Nstar) &
      Const4D(x,y,:,i) = EXP(parval(x,y,indpar%lnFstar(i))) * &
                         pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
    FnuSTAR0(:,:,:) = SUM(Const4D(:,:,:,:),DIM=4)
    IF (PRESENT(FnuSTAR)) FnuSTAR = FnuSTAR0(:,:,:)
    IF (PRESENT(FnuSTAR_tab)) FnuSTAR_tab = Const4D(:,:,:,:)

    ! IF (printimer) &
    !   PRINT*, "[specModel] CAL STAR IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! 4. Screen extinction
    !!----------------------
    CALL REALLOCATE(Const4D,Nx,Ny,Nw,Npabs)
    !! Pabs [/] = EXP( Av [mag] * extinction [/] )
    FORALL (x=1:Nx,y=1:Ny,i=1:Npabs) &
      Const4D(x,y,:,i) = -EXP(parval(x,y,indpar%lnAv(i))) / 1.086_DP * extinct(:)
    Pabs0(:,:,:) = EXP( SUM(Const4D(:,:,:,:),DIM=4) )
    IF (PRESENT(Pabs)) Pabs = Pabs0(:,:,:)
    IF (PRESENT(Pabs_tab)) THEN
      Const4D(:,:,:,:) = EXP( Const4D(:,:,:,:) )
      Pabs_tab = Const4D(:,:,:,:)
    END IF
    
    ! IF (printimer) &
    !   PRINT*, "[specModel] CAL PABS IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! 5. Lines
    !!----------
    CALL REALLOCATE(Const4D,Nx,Ny,Nw,Nline)
    !! FnuLINE [W/m2/sr/Hz] = Iline [W/m2/sr] * gaussLine [Hz-1]
    FORALL (x=1:Nx,y=1:Ny,i=1:Nline) &
      Const4D(x,y,:,i) = EXP(parval(x,y,indpar%lnIline(i))) * &
                         gaussLine(wvl(:), parval(x,y,indpar%Cline(i)), &
                           parval(x,y,indpar%Wline(i)))
    FnuLINE0(:,:,:) = SUM(Const4D(:,:,:,:),DIM=4)
    IF (PRESENT(FnuLINE)) FnuLINE = FnuLINE0(:,:,:)
    IF (PRESENT(FnuLINE_tab)) FnuLINE_tab = Const4D(:,:,:,:)
! print*, "COUCOU"
! print*, parval(6,6,indpar%Cline(1)),parval(6,6,indpar%Wline(1))
! print*, MKS%clight/MKS%micron/parval(6,6,indpar%Cline(1)),MKS%clight/MKS%micron/parval(6,6,indpar%Wline(1))
! print*, gaussLine(wvl(:),parval(6,6,indpar%Cline(1)),parval(6,6,indpar%Wline(1)))

    ! IF (printimer) &
    !   PRINT*, "[specModel] CAL LINE IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! Total model
    !!-------------
    specModel_3D(:,:,:) = (FnuCONT0(:,:,:) + FnuBAND0(:,:,:) + FnuSTAR0(:,:,:)) * &
                          Pabs0(:,:,:) + FnuLINE0(:,:,:)

    IF (printimer) &
      PRINT*, "[specModel] EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
    
  END FUNCTION specModel_3D

  !! 2D version
  !!============
  FUNCTION specModel_2D( wvl, indpar, parval, Qabs, extinct, verbose, &
                         FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, &
                         FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN)        :: wvl, extinct
    TYPE(indpar_type), INTENT(IN)             :: indpar
    REAL(DP), DIMENSION(:,:), INTENT(IN)      :: parval ! 2D (Nx,*Npar)
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL             :: verbose
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT0, &
      FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FnuCONT_tab0, &
      FnuBAND_tab0, FnuSTAR_tab0, Pabs_tab0, FnuLINE_tab0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuCONT, &
      FnuBAND, FnuSTAR, Pabs, FnuLINE
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuCONT_tab, &
      FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
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
  FUNCTION specModel_1D( wvl, indpar, parval, Qabs, extinct, verbose, &
                         FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, &
                         FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN)        :: wvl, extinct
    TYPE(indpar_type), INTENT(IN)             :: indpar
    REAL(DP), DIMENSION(:), INTENT(IN)        :: parval ! 1D (*Npar)
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN) :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL             :: verbose
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT0, &
      FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FnuCONT_tab0, &
      FnuBAND_tab0, FnuSTAR_tab0, Pabs_tab0, FnuLINE_tab0
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuCONT, &
      FnuBAND, FnuSTAR, Pabs, FnuLINE
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuCONT_tab, &
      FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
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
  FUNCTION specModel_gen( wvl, parvec, parname, parinfo, indpar, parval, Qabs, extinct)

    USE utilities, ONLY: DP, trimeq
    USE constants, ONLY: pi, MKS
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN)           :: wvl, extinct
    REAL(DP), DIMENSION(:), INTENT(IN)           :: parvec, parval
    CHARACTER(*), INTENT(IN)                     :: parname
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo
    TYPE(indpar_type), INTENT(IN)                :: indpar
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN)    :: Qabs
    
    INTEGER :: Nw, Ngrid, Ncont, Nline, Nband, Npabs, Nstar
    INTEGER :: i, igrid, j
    REAL(DP)                                     :: Tstar
    REAL(DP), DIMENSION(SIZE(wvl))               :: nu
    LOGICAL :: gridlnMovd2, gridlnT, gridlnIline, gridCline, gridWline
    LOGICAL :: gridlnIband, gridCband, gridWSband, gridWLband, gridlnAv, gridlnFstar
    LOGICAL :: gridCONT, gridBAND, gridLINE, gridSTAR, gridPABS
    REAL(DP), DIMENSION(SIZE(wvl)) :: dblarr1D, spec1D
    REAL(DP), DIMENSION(SIZE(parvec),SIZE(wvl)) :: specModel_gen

    
    !! Preliminaries
    !!---------------
    Nw = SIZE(wvl(:))
    Ngrid = SIZE(parvec(:))
    ! Ncont = SIZE(Qabs(:))
    Ncont = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'CONT')) ) / NparCONT
    Nline = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'LINE')) ) / NparLINE
    Nband = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'BAND')) ) / NparBAND
    Npabs = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'PABS')) ) / NparPABS
    Nstar = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'STAR')) ) / NparSTAR
    Tstar = 5.E4_DP ! [K] high enough to stay in Rayleigh-Jeans limit (T>>5000K for 3um)
    nu(:) = MKS%clight / wvl(:) /MKS%micron
    
    ! Sampled parameter
    gridlnMovd2 = ( parname(1:7) == 'lnMovd2' )
    gridlnT = ( parname(1:3) == 'lnT' )
    gridlnIband = ( parname(1:7) == 'lnIband' )
    gridCband = ( parname(1:5) == 'Cband' )
    gridWSband = ( parname(1:6) == 'WSband' )
    gridWLband = ( parname(1:6) == 'WLband' )
    gridlnIline = ( parname(1:7) == 'lnIline' )
    gridCline = ( parname(1:5) == 'Cline' )
    gridWline = ( parname(1:5) == 'Wline' )
    gridlnFstar = ( parname(1:7) == 'lnFstar' )
    gridlnAv = ( parname(1:4) == 'lnAv' )
    
    ! Constant components
    gridCONT = ( gridlnMovd2 .OR. gridlnT )
    gridBAND = ( gridlnIband .OR. gridCband .OR. gridWSband .OR. gridWLband )
    gridLINE = ( gridlnIline .OR. gridCline .OR. gridWline )
    gridSTAR = gridlnFstar
    gridPABS = gridlnAv
    
    !! Initialization
    specModel_gen(:,:) = 0._DP
    spec1D(:) = 0._DP
    
    
    ! 1. Continuum
    !-------------
    sampCONT: IF (gridCONT) THEN
      IF (gridlnMovd2) THEN
        READ(parname(8:9),*) j
        dblarr1d(:) = MKS%Msun/MKS%pc**2 &
                      * MODIFBB( wvl(:),EXP(parval(indpar%lnT(j))),Qabs(j) )
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP(parvec(igrid)) * dblarr1d(:)
      ELSE IF (gridlnT) THEN
        READ(parname(4:5),*) j
        FORALL (igrid=1:Ngrid) &
        specModel_gen(igrid,:) = EXP(parval(indpar%lnMovd2(j))) &
                                 * MKS%Msun/MKS%pc**2 &
                                 * MODIFBB( wvl(:),EXP(parvec(igrid)),Qabs(j) )
      END IF
      
      DO i=1,Ncont
        IF (i /= j) &
          spec1D(:) = spec1D(:) &
                      + EXP(parval(indpar%lnMovd2(i))) * MKS%Msun/MKS%pc**2 &
                      * MODIFBB( wvl(:),EXP(parval(indpar%lnT(i))),Qabs(i) )
      END DO
    ELSE
      DO i=1,Ncont
        spec1D(:) = spec1D(:) &
                    + EXP(parval(indpar%lnMovd2(i))) * MKS%Msun/MKS%pc**2 &
                    * MODIFBB( wvl(:),EXP(parval(indpar%lnT(i))),Qabs(i) )
      END DO
    END IF sampCONT
    
    
    ! 2. Bands
    !---------
    sampband: IF (gridBAND) THEN
      IF (gridlnIband) THEN
        READ(parname(8:9),*) j
        dblarr1d(:) = LORENTZBAND( wvl(:),parval(indpar%Cband(j)), &
                                   parval(indpar%WSband(j)), &
                                   parval(indpar%WLband(j)) )
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP(parvec(igrid)) * dblarr1d(:)
      ELSE IF (gridCband) THEN
        READ(parname(6:7),*) j
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP(parval(indpar%lnIband(j))) &
                                   * LORENTZBAND( wvl(:),parvec(igrid), &
                                                  parval(indpar%WSband(j)), &
                                                  parval(indpar%WLband(j)) )
      ELSE IF (gridWSband) THEN
        READ(parname(7:8),*) j
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP(parval(indpar%lnIband(j))) &
                                   * LORENTZBAND( wvl(:),parval(indpar%Cband(j)), &
                                                  parvec(igrid), &
                                                  parval(indpar%WLband(j)) )
      ELSE IF (gridWLband) THEN
        READ(parname(7:8),*) j
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP(parval(indpar%lnIband(j))) &
                                   * LORENTZBAND( wvl(:),parval(indpar%Cband(j)), &
                                                  parval(indpar%WSband(j)), &
                                                  parvec(igrid) )
      END IF
     
      DO i=1,Nband
        IF (i /= j) &
          spec1D(:) = spec1D(:) + EXP(parval(indpar%lnIband(i))) &
                      * LORENTZBAND( wvl(:),parval(indpar%Cband(i)), &
                                     parval(indpar%WSband(i)), &
                                     parval(indpar%WLband(i)) )
      END DO
    ELSE
      DO i=1,Nband
        spec1D(:) = spec1D(:) + EXP(parval(indpar%lnIband(i))) &
                    * LORENTZBAND( wvl(:),parval(indpar%Cband(i)), &
                                   parval(indpar%WSband(i)), &
                                   parval(indpar%WLband(i)) )
      END DO
    END IF sampband
    
    
    ! 3. Lines
    !---------
    sampline: IF (gridLINE) THEN
      IF (gridlnIline) THEN
        READ(parname(8:9),*) j
        dblarr1d(:) = GAUSSLINE( wvl(:),parval(indpar%Cline(j)), &
                                 parval(indpar%Wline(j)) )
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP(parvec(igrid)) * dblarr1d(:)
      ELSE IF (gridCline) THEN
        READ(parname(6:7),*) j
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP(parval(indpar%lnIline(j))) &
                                   * GAUSSLINE( wvl(:),parvec(igrid), &
                                                parval(indpar%Wline(j)) )
      ELSE IF (gridWline) THEN
        READ(parname(6:7),*) j
        FORALL (igrid=1:Ngrid) &
          specModel_gen(igrid,:) = EXP(parval(indpar%lnIline(j))) &
                                   * GAUSSLINE( wvl(:),parval(indpar%Cline(j)), &
                                                parvec(igrid) )
      END IF
      
      DO i=1,Nline
        IF (i /= j) &
          spec1D(:) = spec1D(:) + EXP(parval(indpar%lnIline(i))) &
                      * GAUSSLINE( wvl(:), parval(indpar%Cline(i)), &
                                   parval(indpar%Wline(i)) )
      END DO
    ELSE
      DO i=1,Nline
        spec1D(:) = spec1D(:) + EXP(parval(indpar%lnIline(i))) &
                    * GAUSSLINE( wvl(:),parval(indpar%Cline(i)), &
                                 parval(indpar%Wline(i)) )
      END DO
    END IF sampline
    
    
    ! 4. Stars
    !---------
    sampstar: IF (gridSTAR) THEN
      FORALL (igrid=1:Ngrid) &
        specModel_gen(igrid,:) = EXP(parvec(igrid)) &
                                 * pi * BLACKBODY(nu(:), Tstar) / MKS%stefan / Tstar**4
    ELSE
      spec1D(:) = spec1D(:) + EXP(parval(indpar%lnFstar)) &
                  * pi * BLACKBODY(nu(:), Tstar) / MKS%stefan / Tstar**4
    END IF sampstar
    
    
    ! 5. Av
    !------
    
    
    ! 6. Sum all the components
    !--------------------------
    FORALL (igrid=1:Ngrid) &
      specModel_gen(igrid,:) = specModel_gen(igrid,:) + spec1D(:)

  END FUNCTION specModel_gen

  
  FUNCTION specModel_scl( wvl, parinfo, indpar, parval, Qabs, extinct)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN)           :: wvl, extinct
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo
    TYPE(indpar_type), INTENT(IN)                :: indpar
    REAL(DP), DIMENSION(:), INTENT(IN)           :: parval ! 1D
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN)    :: Qabs
    REAL(DP), DIMENSION(SIZE(wvl))               :: specModel_scl
    specModel_scl(:) = RESHAPE( specModel_gen(wvl(:), PARVEC=[0._DP],PARNAME='', &
                                  PARINFO=parinfo, INDPAR=indpar, PARVAL=parval(:), &
                                  QABS=Qabs(:), EXTINCT=extinct(:)), &
                                [SIZE(wvl(:))] )
  END FUNCTION specModel_scl
  
END MODULE auxil
