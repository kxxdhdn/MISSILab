MODULE auxil

  USE utilities, ONLY: DP, STRIKE, WARNING
  USE constants, ONLY: 
  USE inout, ONLY: lenpar
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: NparCONT = 2, NparLINE = 3, NparBAND = 4 ! only in set_indpar
  INTEGER, PARAMETER, PUBLIC :: Ncont_max = 31, Nline_max = 46, Nband_max = 33
  ! INTEGER, PARAMETER :: Npabs_max = 
  ! INTEGER, PARAMETER :: Nstar_max =
  
  PUBLIC :: set_indpar, read_master!, initparam
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
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lnMovd2 ! dust MBB coeff LOG[Msun/pc2]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lnT ! dust temperature (MBB) LOG[K]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lnIline ! line Intensity LOG[W/m2/sr]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Cline ! Center [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Wline ! Width [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lnIband ! band Intensity LOG[W/m2/sr]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Cband ! Center (peak) [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WSband ! Width Short nu side [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WLband ! Width Long nu side [um]
    REAL(DP) :: lnAv ! Attenuation in the V band LOG[mag]
    REAL(DP) :: lnFstar ! Total star surface brightness with dilution factor Omega = r/d LOG[W/m2/sr]
  END TYPE par_type
 
  TYPE, PUBLIC :: parinfo_type
    CHARACTER(lenpar) :: name = ''
    CHARACTER(lenpar) :: comp = ''
    LOGICAL :: fixed = .TRUE.
    LOGICAL, DIMENSION(2) :: limited = [.FALSE., .FALSE.]
    REAL(DP), DIMENSION(2) :: limits = [0._DP, 0._DP]
    LOGICAL :: hyper = .FALSE.
    CHARACTER(lenpar) :: tied = ''
    REAL(DP) :: value = 0._DP
    INTEGER :: itied = 0
    LOGICAL :: model = .TRUE.
    REAL(DP) :: mean = 0._DP
    REAL(DP) :: sigma = 0._DP
    INTEGER :: ind = -1
  END TYPE parinfo_type

  TYPE, PUBLIC :: indpar_type
    INTEGER, DIMENSION(Ncont_max) :: lnMovd2 = -1 ! dust MBB coeff LOG[Msun/pc2]
    INTEGER, DIMENSION(Ncont_max) :: lnT = -1 ! dust temperature (MBB) LOG[K]
    INTEGER, DIMENSION(Nline_max) :: lnIline = -1 ! line Intensity LOG[W/m2/sr]
    INTEGER, DIMENSION(Nline_max) :: Cline = -1 ! Center [um]
    INTEGER, DIMENSION(Nline_max) :: Wline = -1 ! Width [um]
    INTEGER, DIMENSION(Nband_max) :: lnIband = -1 ! band Intensity LOG[W/m2/sr]
    INTEGER, DIMENSION(Nband_max) :: Cband = -1 ! Center (peak) [um]
    INTEGER, DIMENSION(Nband_max) :: WSband = -1 ! Width Short nu side [um]
    INTEGER, DIMENSION(Nband_max) :: WLband = -1 ! Width Long nu side [um]
    INTEGER :: lnAv = -1 ! Attenuation in the V band LOG[mag]
    INTEGER :: lnFstar = -1 ! Total star surface brightness with dilution factor Omega = r/d LOG[W/m2/sr]
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

    INTEGER :: i, Ncont, Nline, Nband, Nextra
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
      CALL IWHERE(TRIMEQ(parinfo%name(:),'lnAv'), indpar%lnAv)
      
    END IF
    
    IF (STARset) THEN
      CALL IWHERE(TRIMEQ(parinfo%name(:),'lnFstar'),indpar%lnFstar)
      
    END IF

    IF (extraset) THEN
      Nextra = COUNT(TRIMEQ(parinfo(:)%comp, 'extra'))
      IF (Nextra > 0) CALL IWHERE(TRIMEQ(parinfo(:)%comp,'extra'),indpar%extra)
      
    END IF
    
  END SUBROUTINE set_indpar

  !!-------------------------------------------------------
  !!
  !! Read the input master file for the Chi2/Bayesian run
  !!
  !!-------------------------------------------------------

  SUBROUTINE read_master( wavALL, &
                          Nmcmc, verbose, NiniMC, &!robust_RMS, robust_cal, skew_RMS, &
                          calib, newseed, newinit, &
                          labQ, labL, labB, Qabs, &
                          Ncont, Nband, Nline, Nextra, dostop, &
                          parinfo, parhypinfo, parextinfo, &
                          indpar, Npar, &!Nparmod, Nparhyp, Ncorrhyp, Ncorr, &
                          spec_unit)

    ! USE datable, ONLY: TABLine, TABand
    USE utilities, ONLY: DP, strike, trimeq, trimlr, pring!, verbatim
    USE inout, ONLY: lenpar, lenline, read_input_line, read_hdf5, h5ext
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: incrarr, iwhere, reallocate, closest
    USE interpolation, ONLY: interp_lin_sorted
    USE grain_optics, ONLY: lendustQ, rho_grain, read_optics
    IMPLICIT NONE

    CHARACTER(lenpar), INTENT(INOUT), OPTIONAL :: spec_unit ! ??
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: wavALL ! Interpolate Qabs if wavALL presents
    
    INTEGER, PARAMETER :: unit = 1
    CHARACTER(*), PARAMETER :: filmas = './input_fitMIR_master'//h5ext
    CHARACTER(*), PARAMETER :: filmod = './input_fitMIR_model'//h5ext
    CHARACTER(*), PARAMETER :: filext = './input_fitMIR_extra'//h5ext

    ! CHARACTER(lenpar) :: spec_unit0
    CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ0, labL0, labB0
    ! INTEGER :: Nmcmc0, NiniMC0
    INTEGER :: Ncont0, Nband0, Nline0, Nextra0, Npar0
    ! INTEGER :: Nparall, Nparmod0, Nparhyp0,  Ncorrhyp0, Ncorr0
    ! INTEGER :: i, iostat, ipar
    ! TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE :: parinfo0
    ! TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE :: parinfoextra, parinfall
    ! LOGICAL :: verbose0!, robust_RMS0, robust_cal0, skew_RMS0
    ! LOGICAL :: calib0, newseed0, newinit0, dostop0

    CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE   :: strarr1d, Tarr1d
    CHARACTER(lenpar), DIMENSION(:,:), ALLOCATABLE :: strarr2d
    INTEGER, DIMENSION(:), ALLOCATABLE             :: intarr1d
    REAL(DP), DIMENSION(:), ALLOCATABLE            :: dblarr1d
    REAL(DP), DIMENSION(:,:), ALLOCATABLE          :: dblarr2d
    INTEGER                                 :: i, Nr, Nw
    INTEGER, DIMENSION(:), ALLOCATABLE      :: indrad ! radius index
    REAL(DP), PARAMETER                     :: a0 = 1.E-2_DP ! grain radius
    REAL(DP), DIMENSION(:), ALLOCATABLE     :: wave0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: radius ! (Nr, Nw)
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Qabsall ! (Ncont, Nr, Nw)
    
    CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL ::labQ, labB, labL
    INTEGER, INTENT(OUT), OPTIONAL :: Nmcmc, NiniMC, Ncont, Nband, Nline, Nextra
    INTEGER, INTENT(OUT), OPTIONAL :: Npar!, Nparmod, Nparhyp, Ncorrhyp, Ncorr
    TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Qabs
    TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: &
      parinfo, parhypinfo, parextinfo!, parmodinfo
    TYPE(indpar_type), INTENT(OUT), OPTIONAL :: indpar
    LOGICAL, INTENT(OUT), OPTIONAL :: verbose!, robust_RMS, robust_cal, skew_RMS
    LOGICAL, INTENT(OUT), OPTIONAL :: calib, newseed, newinit, dostop
    
    !! Read the input master file
    !!----------------------------

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

    !! input_fitMIR_model
    IF (PRESENT(spec_unit)) THEN
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='spec_unit')
      spec_unit = strarr1d(1)
    END IF
    
    IF (PRESENT(labQ)) THEN
      CALL READ_HDF5(STRARR1D=labQ0, FILE=filmod, NAME='labQ')
      labQ = labQ0(:)
    ELSE
      labQ0 = ['ACH2_Z96             ']
    END IF
    Ncont0 = SIZE(labQ0(:))
    
    IF (PRESENT(labL)) THEN
      CALL READ_HDF5(STRARR1D=labL0, FILE=filmod, NAME='labL')
    ELSE
      ALLOCATE(labL0(0))
    END IF
    Nline0 = SIZE(labL0(:))
    IF (PRESENT(labL)) labL = labL0(:)

    IF (PRESENT(labB)) THEN
      CALL READ_HDF5(STRARR1D=labB0, FILE=filmod, NAME='labB')
    ELSE
      ALLOCATE(labB0(0))
    END IF
    Nband0 = SIZE(labB0(:))
    IF (PRESENT(labB)) labB = labB0(:)

    IF (PRESENT(Ncont)) Ncont = Ncont0
    IF (PRESENT(Nband)) Nband = Nband0
    IF (PRESENT(Nline)) Nline = Nline0

    !! Read optical properties (Qabs)
    !!--------------------------------
    IF (PRESENT(Qabs)) THEN
      !! Qabs % rho(kg/m3); wave(micron); nu(Hz); kappa(m2/kg)
      CALL READ_OPTICS(labQ0, WAVE=wave0, RADIUSALL=radius, QABSALL=Qabsall)    
      
      Nw = SIZE(wave0)

      ALLOCATE(indrad(Ncont0), Qabs(Ncont0))
      
      DO i=1,Ncont0
        
        ALLOCATE(Qabs(i)%kappa(Nw))
        
        Nr = SIZE(radius(i,:))
        Qabs(i)%rho = rho_grain(labQ0(i))
        Qabs(i)%wave = wave0
        Qabs(i)%nu = MKS%clight/MKS%micron / wave0
        indrad(i) = closest(radius(i,:), a0)
        ! PRINT*, 'Radius of ', labQ0(i), ': ', radius(i,indrad(i)), '->', indrad(i)
        Qabs(i)%kappa = 3._DP/4._DP/rho_grain(labQ0(i)) * Qabsall(i,indrad(i),:)/radius(i,indrad(i))
        
        !! [Optional] Interpolate Qabs to input wave grid
        IF (PRESENT(wavALL)) THEN
          Qabs(i)%wave = wavALL
          Qabs(i)%nu = MKS%clight/MKS%micron / wavALL
          Qabs(i)%kappa = interp_lin_sorted(Qabs(i)%kappa, wave0, wavALL, &
                                                XLOG=.TRUE., YLOG=.TRUE., FORCE=.TRUE.)
          
        END IF
      END DO
      !! Free memory space
      DEALLOCATE(wave0, radius, Qabsall)
      
    END IF
    
    !! Extra parameters
    !!------------------
    Nextra0 = 0

    !! input_fitMIR_extra
    IF (PRESENT(Nextra)) THEN
      CALL READ_HDF5(INTARR1D=intarr1d, FILE=filext, NAME='Nextra')
      Nextra0 = intarr1d(1)
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
    Npar0 = 2*Ncont0 + 4*Nband0 + 3*Nline0 + 2 + Nextra0
    IF (PRESENT(Npar)) Npar = Npar0
    
    IF (PRESENT(parinfo)) THEN
      !! General parameter structure
      ALLOCATE(parinfo(Npar0), Tarr1d(Npar0))
      Tarr1d(:) = 'T'
      parinfo(:)%ind = [(i,i=1,Npar0)]
      
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo name')
      print*, Npar0, Ncont0
      parinfo(:)%name = strarr1d(:)
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo comp')
      parinfo(:)%comp = strarr1d(:)
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo fixed')
      parinfo(:)%fixed = TRIMEQ(strarr1d(:),Tarr1d(:))
      CALL READ_HDF5(STRARR2D=strarr2d, FILE=filmod, NAME='parinfo limited')
      CALL READ_HDF5(DBLARR2D=dblarr2d, FILE=filmod, NAME='parinfo limits')
      DO i=1,2
        parinfo(:)%limited(i) = TRIMEQ(strarr2d(i,:),Tarr1d(:))
        parinfo(:)%limits(i) = dblarr2d(i,:)
      END DO
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo hyper')
      parinfo(:)%hyper = TRIMEQ(strarr1d(:),Tarr1d(:))
      CALL READ_HDF5(STRARR1D=strarr1d, FILE=filmod, NAME='parinfo tied')
      parinfo(:)%tied = strarr1d(:)
      CALL READ_HDF5(DBLARR1D=dblarr1d, FILE=filmod, NAME='parinfo value')
      parinfo(:)%value = dblarr1d(:)
      
      !! Indices
      IF (PRESENT(indpar)) CALL SET_INDPAR(indpar, parinfo(:))
      
    ELSE
      IF (PRESENT(indpar)) CALL STRIKE("READ_MASTER", "need parinfo to set indpar")
      
    END IF
    
    !! Free memory space
    !!-------------------
    DEALLOCATE (labQ0, labL0, labB0)
    
  END SUBROUTINE read_master

  !!-------------------------------------------------------
  !!
  !!     Automatic initialization of model parameters
  !!
  !!-------------------------------------------------------
  SUBROUTINE initparam( NiniMC, ind, par, parinfo, itied, mask, &
                        newinit, filobs )
    
    USE utilities, ONLY: DP, trimeq, trimlr, pring, isNaN
    USE arrays, ONLY: iwhere, reallocate
    USE inout, ONLY: read_hdf5, lenpar
    USE constants, ONLY: MKS
    USE statistics, ONLY: mean
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:,:,:), INTENT(INOUT) :: par
    INTEGER, INTENT(IN) :: NiniMC
    TYPE(indpar_type), INTENT(IN) :: ind
    INTEGER, DIMENSION(:), INTENT(IN) :: itied
    LOGICAL, DIMENSION(:,:,:), INTENT(IN) :: mask
    LOGICAL, INTENT(IN), OPTIONAL :: newinit
    CHARACTER(*), INTENT(IN), OPTIONAL :: filobs
    
    INTEGER :: i, j, ipar, x, y, Nx, Ny, Npar, Nextra, Nmc
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: iniparval
    REAL(DP), DIMENSION(:), ALLOCATABLE :: parval
    REAL(DP), DIMENSION(SIZE(par,1),SIZE(par,2),MAX(NiniMC,1)) :: theta
    REAL(DP), DIMENSION(SIZE(par,1),SIZE(par,2)) :: limi2, lims2
    LOGICAL :: newini
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: maskxy
    CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: iniparname
    
    !! Preliminaries
    Nx = SIZE(par(:,:,:,:),1)
    Ny = SIZE(par(:,:,:,:),2)
    Npar = SIZE(parinfo(:))
    IF (PRESENT(newinit)) THEN
      newini = newinit
    ELSE
      newini = .FALSE.
    END IF
    ALLOCATE(parval(Npar))
    ALLOCATE(maskxy(Nx,Ny))
    FORALL (x=1:Nx,y=1:Ny) maskxy(x,y) = ANY(mask(x,y,:))

    newinipar: IF (.NOT. newini) THEN
      
      !!----------------------------------------------------
      !! I. Automatically generate initial parameter values
      !!----------------------------------------------------

      IF (NiniMC==0) THEN

        !! a. Simple initilization
        !!-------------------------
        
        !! lnMovd2
        IF (parinfo(ind%lnMovd2)%fixed) THEN
          par(:,:,ind%lnMovd2,:) = parinfo(ind%lnMovd2)%value
        ELSE
          par(:,:,ind%lnMovd2,1) =
        END IF

        !! lnT
        IF (parinfo(ind%lnT)%fixed) THEN
          par(:,:,ind%lnT,:) = parinfo(ind%lnT)%value
        ELSE
          par(:,:,ind%lnT,1) =
        END IF

      ELSE

        !! b. MC initilization
        !!---------------------

        !! lnMovd2
        IF (parinfo(ind%lnMovd2)%fixed) THEN
          par(:,:,ind%lnMovd2,:) = parinfo(ind%lnMovd2)%value
        ELSE
          limi = MERGE( parinfo(ind%lnMovd2)%limits(1), &
                        MINVAL(), &
                        parinfo(ind%lnMovd2)%limits(1)) ! lim inf
          lims = MERGE( parinfo(ind%lnMovd2)%limits(1), &
                        MINVAL(), &
                        parinfo(ind%lnMovd2)%limits(1)) ! lim sup
          CALL RANDOM_NUMBER(theta(:,:,:))
          par(:,:,ind%lnMovd2,:) = (lims - limi) * theta(:,:,:) + limi
        END IF

    ELSE
      
      !!----------------------------------------------------
      !! II. Read initial parameters from input file
      !!----------------------------------------------------

      !! Read the initial parameters from the file
      CALL READ_HDF5(STRARR1D=iniparname, FILE=filobs, &
                     NAME='Initial parameter label')
      CALL READ_HDF5(DBLARR3D=iniparval, FILE=filobs, &
                     NAME='Initial parameter value')

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
    !! temp [K]; blackbody [W/m2/sr/Hz]; 
    !! output array [MJy m2/sr/kg]
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
    nu(:) = MKS%clight/MKS%micron / dblarr(:) ! should be the same with Qabs%nu

    FORALL (it=1:Nt) &
      modifBB_gen(it,:) = Qabs(it)%kappa * BLACKBODY(nu(:), temp(it)) ! W/sr/Hz/kg
    !! Unit onversion (from [W/m2/Hz] to [MJy])
    modifBB_gen(:,:) = modifBB_gen(:,:) /MKS%Jy/1.E6
    
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
        nuIN(:) = MKS%clight/MKS%micron / dblarr(:)
        nuref(:) = MKS%clight/MKS%micron / ref(:)
        nusig(:) = MKS%clight/MKS%micron * (1./(ref(:)-sig(:)) - 1./(ref(:)+sig(:))) / 2._DP

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
  PURE FUNCTION lorentzBand_gen( dblarr, ref, sigS, sigL, w2nu )
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
        nuIN(:) = MKS%clight/MKS%micron / dblarr(:)
        nuref(:) = MKS%clight/MKS%micron / ref(:)
        nusigS(:) = MKS%clight/MKS%micron * (1./(ref(:)-sigS(:)) - 1./ref(:))
        nusigL(:) = MKS%clight/MKS%micron * (1./ref(:) - 1./(ref(:)+sigL(:)))
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
    !! In addition we denote lambda=2*dnu_long (analogue of FWHM)
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
    !! nu2w=.TRUE. if input is in frequency
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
        waveIN = MKS%clight/MKS%micron / dblarr
      END IF
    ELSE
      waveIN = dblarr
      
    END IF
    
    CALL read_hdf5(DBLARR1D=wave, NAME='lambda (micron)', &
                   N1=Nw0, FILE='../data/extcurve'//h5ext)
    CALL read_hdf5(DBLARR1D=Cext, NAME='C_extovH (cm^2ovH)', &
                   FILE='../data/extcurve'//h5ext)
    wave_V = 0.5470_DP
    Cext_V = Cext(closest(wave, wave_V))
    !! interpolated with respect to Wavelength
    extCurve = interp_lin_sorted(Cext/Cext_V, wave, waveIN, &
                                 XLOG=.TRUE., YLOG=.TRUE., FORCE=.TRUE.)

    DEALLOCATE(wave, Cext)

  END FUNCTION extCurve

  !!---------------------------------------
  !! Total model function for chi2 calling
  !!---------------------------------------

  !! 3D version
  !!============
  FUNCTION specModel_3D( wvl, indpar, parval, Qabs, &
                         FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE )

    USE utilities, ONLY: DP, trimeq, trimlr, pring
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: reallocate
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN)                           :: wvl
    TYPE(indpar_type), INTENT(IN)                                :: indpar
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)                       :: parval ! 3D (Nx,Ny,*Npar)
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN)                    :: Qabs
    INTEGER :: Nw, Nx, Ny, Ncont, Nline, Nband
    INTEGER :: x, y, i
    REAL(DP)                                                     :: Tstar
    REAL(DP), DIMENSION(SIZE(wvl))                               :: nu, extinction
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE                    :: Const4D ! (Nx,Ny,Nw,Ncomp)
    REAL(DP), DIMENSION(SIZE(parval,1),SIZE(parval,2),SIZE(wvl)) :: FnuCONT0, &
      FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuCONT, &
      FnuBAND, FnuSTAR, Pabs, FnuLINE
    REAL(DP), DIMENSION(SIZE(parval,1),SIZE(parval,2),SIZE(wvl)) :: specModel_3D

    !! Preliminaries
    !!---------------
    Nw = SIZE(wvl(:))
    Nx = SIZE(parval(:,:,:),1)
    Ny = SIZE(parval(:,:,:),2)
    Ncont = SIZE(Qabs(:))
    Nband = COUNT(indpar%lnIband(:) .NE. -1)
    Tstar = 5.E4_DP ! high enough to stay in Rayleigh-Jeans limit [K]
    Nline = COUNT(indpar%lnIline(:) .NE. -1)
    nu(:) = MKS%clight/MKS%micron / wvl(:)
    extinction(:) = extCurve(wvl(:))

    !! Initialization
    specModel_3D(:,:,:) = 0._DP
    
    !! 1. Continuum
    !!--------------
    CALL REALLOCATE(Const4D,Nx,Ny,Nw,Ncont)
    FORALL (x=1:Nx,y=1:Ny,i=1:Ncont) &
      Const4D(x,y,:,i) = EXP(parval(x,y,indpar%lnMovd2(i))) * MKS%Msun/MKS%pc**2 * &
                         modifBB(wvl(:), EXP(parval(x,y,indpar%lnT(i))), Qabs(i))
    FnuCONT0(:,:,:) = SUM(Const4D(:,:,:,:),DIM=4)

    IF (PRESENT(FnuCONT)) THEN
      ALLOCATE(FnuCONT(Nx,Ny,Nw))
      FnuCONT(:,:,:) = FnuCONT0(:,:,:)
    END IF
    
    !! 2. Bands
    !!----------
    CALL REALLOCATE(Const4D,Nx,Ny,Nw,Nband)
    FORALL (x=1:Nx,y=1:Ny,i=1:Nband) &
      Const4D(x,y,:,i) = EXP(parval(x,y,indpar%lnIband(i))) /MKS%Jy/1.E6 * &
                         lorentzBand(wvl(:), parval(x,y,indpar%Cband(i)), &
                           parval(x,y,indpar%WSband(i)), parval(x,y,indpar%WLband(i)), .TRUE.)
    FnuBAND0(:,:,:) = SUM(Const4D(:,:,:,:),DIM=4)

    IF (PRESENT(FnuBAND)) THEN
      ALLOCATE(FnuBAND(Nx,Ny,Nw))
      FnuBAND(:,:,:) = FnuBAND0(:,:,:)
    END IF
    
    !! 3. Stellar Continuum
    !!----------------------
    FORALL (x=1:Nx,y=1:Ny) &
      FnuSTAR0(x,y,:) = EXP(parval(x,y,indpar%lnFstar)) /MKS%Jy/1.E6 * &
                        pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

    IF (PRESENT(FnuSTAR)) THEN
      ALLOCATE(FnuSTAR(Nx,Ny,Nw))
      FnuSTAR(:,:,:) = FnuSTAR0(:,:,:)
    END IF
    
    !! 4. Screen extinction
    !!----------------------
    FORALL (x=1:Nx,y=1:Ny) &
      Pabs0(x,y,:) = EXP(-EXP(parval(x,y,indpar%lnAv))/1.086_DP * extinction(:))

    IF (PRESENT(Pabs)) THEN
      ALLOCATE(Pabs(Nx,Ny,Nw))
      Pabs(:,:,:) = Pabs0(:,:,:)
    END IF
    
    !! 5. Lines
    !!----------
    CALL REALLOCATE(Const4D,Nx,Ny,Nw,Nline)
    FORALL (x=1:Nx,y=1:Ny,i=1:Nline) &
      Const4D(x,y,:,i) = EXP(parval(x,y,indpar%lnIline(i))) /MKS%Jy/1.E6 * &
                         gaussLine(wvl(:), parval(x,y,indpar%Cline(i)), &
                           parval(x,y,indpar%Wline(i)), .TRUE.)
    FnuLINE0(:,:,:) = SUM(Const4D(:,:,:,:),DIM=4)

    IF (PRESENT(FnuLINE)) THEN
      ALLOCATE(FnuLINE(Nx,Ny,Nw))
      FnuLINE(:,:,:) = FnuLINE0(:,:,:)
    END IF
    
    !! Total model
    !!-------------
    specModel_3D(:,:,:) = (FnuCONT0(:,:,:) + FnuBAND0(:,:,:) + FnuSTAR0(:,:,:)) * &
                          Pabs0(:,:,:) + FnuLINE0(:,:,:)

  END FUNCTION specModel_3D

  !! 2D version
  !!============
  FUNCTION specModel_2D( wvl, indpar, parval, Qabs )

    USE utilities, ONLY: DP, trimeq, trimlr, pring
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: reallocate
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN)            :: wvl
    TYPE(indpar_type), INTENT(IN)                 :: indpar
    REAL(DP), DIMENSION(:,:), INTENT(IN)          :: parval ! 2D (*Npar,Npar)
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN)     :: Qabs
    INTEGER :: Nw, Npar, Ncont, Nline, Nband
    INTEGER :: ipar, i
    REAL(DP)                                      :: Tstar
    REAL(DP), DIMENSION(SIZE(wvl))                :: nu, extinction
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE       :: Const3D ! (Npar,Nw,Ncomp)
    REAL(DP), DIMENSION(SIZE(parval,2),SIZE(wvl)) :: FnuCONT0, &
      FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    REAL(DP), DIMENSION(SIZE(parval,2),SIZE(wvl)) :: specModel_2D

    !! Preliminaries
    !!---------------
    Nw = SIZE(wvl(:))
    Npar = SIZE(parval(:,:),2)
    Ncont = SIZE(Qabs(:))
    Nband = COUNT(indpar%lnIband(:) .NE. -1)
    Tstar = 5.E4_DP ! high enough to stay in Rayleigh-Jeans limit [K]
    Nline = COUNT(indpar%lnIline(:) .NE. -1)
    nu(:) = MKS%clight/MKS%micron / wvl(:)
    extinction(:) = extCurve(wvl(:))

    !! Initialization
    specModel_2D(:,:) = 0._DP
    
    !! 1. Continuum
    !!--------------
    CALL REALLOCATE(Const3D,Npar,Nw,Ncont)
    FORALL (ipar=1:Npar, i=1:Ncont) &
      Const3D(ipar,:,i) = EXP(parval(indpar%lnMovd2(i),ipar)) * MKS%Msun/MKS%pc**2 * &
                          modifBB(wvl(:), EXP(parval(indpar%lnT(i),ipar)), Qabs(i))
    FnuCONT0(:,:) = SUM(Const3D(:,:,:),DIM=3)
    
    !! 2. Bands
    !!----------
    CALL REALLOCATE(Const3D,Npar,Nw,Nband)
    FORALL (ipar=1:Npar, i=1:Nband) &
      Const3D(ipar,:,i) = EXP(parval(indpar%lnIband(i),ipar)) /MKS%Jy/1.E6 * &
                          lorentzBand(wvl(:), parval(indpar%Cband(i),ipar), &
                            parval(indpar%WSband(i),ipar), parval(indpar%WLband(i),ipar), .TRUE.)
    FnuBAND0(:,:) = SUM(Const3D(:,:,:),DIM=3)
    
    !! 3. Stellar Continuum
    !!----------------------
    FORALL (ipar=1:Npar) &
      FnuSTAR0(ipar,:) = EXP(parval(indpar%lnFstar,ipar)) /MKS%Jy/1.E6 * &
                         pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

    !! 4. Screen extinction
    !!----------------------
    FORALL (ipar=1:Npar) &
      Pabs0(ipar,:) = EXP(-EXP(parval(indpar%lnAv,ipar))/1.086_DP * extinction(:))

    !! 5. Lines
    !!----------
    CALL REALLOCATE(Const3D,Npar,Nw,Nline)
    FORALL (ipar=1:Npar, i=1:Nline) &
      Const3D(ipar,:,i) = EXP(parval(indpar%lnIline(i),ipar)) /MKS%Jy/1.E6 * &
                          gaussLine(wvl(:), parval(indpar%Cline(i),ipar), &
                            parval(indpar%Wline(i),ipar), .TRUE.)
    FnuLINE0(:,:) = SUM(Const3D(:,:,:),DIM=3)

    !! Total model
    !!-------------
    specModel_2D(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                        Pabs0(:,:) + FnuLINE0(:,:)

  END FUNCTION specModel_2D

  !! 1D version
  !!============
  ! FUNCTION specModel_1D( wvl, indpar, parval, Qabs, &
  !                        FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE )
  !   USE utilities, ONLY: DP
  !   IMPLICIT NONE
  !   REAL(DP), DIMENSION(:), INTENT(IN)        :: wvl
  !   TYPE(indpar_type), INTENT(IN)             :: indpar
  !   REAL(DP), DIMENSION(:), INTENT(IN)        :: parval ! 1D (*Npar)
  !   TYPE(Qabs_type), DIMENSION(:), INTENT(IN) :: Qabs
  !   REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuCONT, &
  !     FnuBAND, FnuSTAR, Pabs, FnuLINE
  !   REAL(DP), DIMENSION(SIZE(wvl))            :: specModel_1D
  !   specModel_1D(:) = RESHAPE( specModel_3D(wvl(:), INDPAR=indpar, &
  !                                PARVAL=parval, QABS=Qabs, &
  !                                FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
  !                                PABS=Pabs, FNULINE=FnuLINE), [SIZE(wvl(:))] )
  ! END FUNCTION specModel_1D
  
  FUNCTION specModel_1D( wvl, indpar, parval, Qabs, &
                         FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE )

    USE utilities, ONLY: DP, trimeq, trimlr, pring
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: reallocate
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE
  
    REAL(DP), DIMENSION(:), INTENT(IN)        :: wvl
    TYPE(indpar_type), INTENT(IN)             :: indpar
    REAL(DP), DIMENSION(:), INTENT(IN)        :: parval ! 1D (*Npar)
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN) :: Qabs
    INTEGER :: Nw, Ncont, Nline, Nband
    INTEGER :: i
    REAL(DP)                                  :: Tstar
    REAL(DP), DIMENSION(SIZE(wvl))            :: nu, extinction
    REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: Const2D ! (Nw,Ncomp)
    REAL(DP), DIMENSION(SIZE(wvl))            :: FnuCONT0, &
      FnuBAND0, FnuSTAR0, Pabs0, FnuLINE0
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuCONT, &
      FnuBAND, FnuSTAR, Pabs, FnuLINE
    REAL(DP), DIMENSION(SIZE(wvl))            :: specModel_1D

    !! Preliminaries
    !!---------------
    Nw = SIZE(wvl(:))
    Ncont = SIZE(Qabs(:))
    Nband = COUNT(indpar%lnIband(:) .NE. -1)
    Tstar = 5.E4_DP ! high enough to stay in Rayleigh-Jeans limit [K]
    Nline = COUNT(indpar%lnIline(:) .NE. -1)
    nu(:) = MKS%clight/MKS%micron / wvl(:)
    extinction(:) = extCurve(wvl(:))

    !! Initialization
    specModel_1D(:) = 0._DP
  
    !! 1. Continuum
    !!--------------
    CALL REALLOCATE(Const2D,Nw,Ncont)
    FORALL (i=1:Ncont) &
      Const2D(:,i) = EXP(parval(indpar%lnMovd2(i))) * MKS%Msun/MKS%pc**2 * &
                     modifBB(wvl(:), EXP(parval(indpar%lnT(i))), Qabs(i))
    FnuCONT0(:) = SUM(Const2D,DIM=2)
    
    IF (PRESENT(FnuCONT)) THEN
      ALLOCATE(FnuCONT(Nw))
      FnuCONT(:) = FnuCONT0(:)
    END IF
  
    !! 2. Bands
    !!----------
    CALL REALLOCATE(Const2D,Nw,Nband)
    FORALL (i=1:Nband) &
      Const2D(:,i) = EXP(parval(indpar%lnIband(i))) /MKS%Jy/1.E6 * &
                     lorentzBand(wvl(:), parval(indpar%Cband(i)), &
                       parval(indpar%WSband(i)), parval(indpar%WLband(i)), .TRUE.)
    FnuBAND0(:) = SUM(Const2D,DIM=2)

    IF (PRESENT(FnuBAND)) THEN
      ALLOCATE(FnuBAND(Nw))
      FnuBAND(:) = FnuBAND0(:)
    END IF
  
    !! 3. Stellar Continuum
    !!----------------------
    FnuSTAR0(:) = EXP(parval(indpar%lnFstar)) /MKS%Jy/1.E6 * &
                  pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

    IF (PRESENT(FnuSTAR)) THEN
      ALLOCATE(FnuSTAR(Nw))
      FnuSTAR(:) = FnuSTAR0(:)
    END IF

    !! 4. Screen extinction
    !!----------------------
    Pabs0(:) = EXP(-EXP(parval(indpar%lnAv))/1.086_DP * extinction(:))

    IF (PRESENT(Pabs)) THEN
      ALLOCATE(Pabs(Nw))
      Pabs(:) = Pabs0(:)
    END IF

    !! 5. Lines
    !!----------
    CALL REALLOCATE(Const2D,Nw,Nline)
    FORALL (i=1:Nline) &
      Const2D(:,i) = EXP(parval(indpar%lnIline(i))) /MKS%Jy/1.E6 * &
                     gaussLine(wvl(:), parval(indpar%Cline(i)), &
                       parval(indpar%Wline(i)), .TRUE.)
    FnuLINE0(:) = SUM(Const2D,DIM=2)

    IF (PRESENT(FnuLINE)) THEN
      ALLOCATE(FnuLINE(Nw))
      FnuLINE(:) = FnuLINE0(:)
    END IF

    !! Total model
    !!-------------
    specModel_1D(:) = (FnuCONT0(:) + FnuBAND0(:) + FnuSTAR0(:)) * &
                      Pabs0(:) + FnuLINE0(:)

  END FUNCTION specModel_1D
  
  !! Generic intertface
  !!====================
  FUNCTION specModel_gen( wvl, parvec, parname, parinfo, indpar, parval, Qabs )

    USE utilities, ONLY: DP, trimeq, trimlr, pring
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: reallocate
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN)           :: wvl
    REAL(DP), DIMENSION(:), INTENT(IN)           :: parvec
    CHARACTER(*), INTENT(IN)                     :: parname
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo
    TYPE(indpar_type), INTENT(IN)                :: indpar
    REAL(DP), DIMENSION(:), INTENT(IN)           :: parval ! 1D
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN)    :: Qabs
    INTEGER :: Nw, Ngrid, Ncont, Nline, Nband
    INTEGER :: i, igrid, iw, i1, i2, i3!, itied
    REAL(DP)                                     :: Tstar
    REAL(DP), DIMENSION(SIZE(wvl))               :: nu, extinction
    LOGICAL :: gridlnMovd2, gridlnT, gridlnIline, gridCline, gridWline
    LOGICAL :: gridlnIband, gridCband, gridWSband, gridWLband, gridlnAv, gridlnFstar
    REAL(DP), DIMENSION(SIZE(parvec),SIZE(wvl)) :: FnuCONT0, FnuBAND0, FnuSTAR0
    REAL(DP), DIMENSION(SIZE(parvec),SIZE(wvl)) :: Pabs0, FnuLINE0
    REAL(DP) :: Const
    REAL(DP), DIMENSION(SIZE(wvl)) :: Const1D
    ! REAL(DP), DIMENSION(SIZE(parvec),SIZE(wvl)) :: Const2D
    ! REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL   :: Pabs
    ! REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuCONT
    ! REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuBAND
    ! REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL   :: FnuSTAR
    ! REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: FnuLINE
    REAL(DP), DIMENSION(SIZE(parvec),SIZE(wvl)) :: specModel_gen

    !! Preliminaries
    !!---------------
    Nw = SIZE(wvl(:))
    Ngrid = SIZE(parvec(:))
    ! Ncont = SIZE(Qabs(:))
    Ncont = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'CONT')) ) / NparCONT
    Nline = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'LINE')) ) / NparLINE
    Nband = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, 'BAND')) ) / NparBAND
    Tstar = 5.E4_DP ! high enough to stay in Rayleigh-Jeans limit [K]
    nu(:) = MKS%clight/MKS%micron / wvl(:)
    extinction(:) = extCurve(wvl(:))

    !! Is the parameter tied to another one?
    ! tied = ANY(TRIMEQ(parinfo(:)%tied, parname))
    ! IF (tied) CALL IWHERE( TRIMEQ(parinfo(:)%tied, parname), itied )

    !! Initialization
    specModel_gen(:,:) = 0._DP
    
    !! 1. Continuum
    !!--------------
    loopCONT: DO i=1,Ncont
      gridlnMovd2 = TRIMEQ(parname, 'lnMovd2'//TRIMLR(PRING(i)))
      gridlnT = TRIMEQ(parname, 'lnT'//TRIMLR(PRING(i)))
      
      IF (gridlnMovd2) THEN
        !! CONT *lnMovd2
        Const1D(:) = MKS%Msun/MKS%pc**2 * modifBB(wvl(:), EXP(parval(indpar%lnT(i))), Qabs(i))
        FORALL (igrid=1:Ngrid) &
          FnuCONT0(igrid,:) = EXP(parvec(igrid)) * Const1D(:)
        !! CONT
        Const1D(:) = 0._DP
        DO i1=1,Ncont
          IF (i1.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(i1))) * MKS%Msun/MKS%pc**2 * &
                         modifBB(wvl(:), EXP(parval(indpar%lnT(i1))), Qabs(i1))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = FnuCONT0(:,iw) + Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO i2=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(i2))) /MKS%Jy/1.E6 * &
                       lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                         parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR
        Const1D(:) = EXP(parval(indpar%lnFstar)) /MKS%Jy/1.E6 * &
                     pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO i3=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(i3))) /MKS%Jy/1.E6 * &
                       gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                         parval(indpar%Wline(i3)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-EXP(parval(indpar%lnAv))/1.086_DP * extinction(iw))
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      ELSE IF (gridlnT) THEN
        !! CONT *lnT
        Const = EXP(parval(indpar%lnMovd2(i))) * MKS%Msun/MKS%pc**2
        FORALL (igrid=1:Ngrid) &
          FnuCONT0(igrid,:) = Const * modifBB(wvl(:), EXP(parvec(igrid)), Qabs(i))
        !! CONT
        Const1D(:) = 0._DP
        DO i1=1,Ncont
          IF (i1.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(i1))) * MKS%Msun/MKS%pc**2 * &
                         modifBB(wvl(:), EXP(parval(indpar%lnT(i1))), Qabs(i1))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = FnuCONT0(:,iw) + Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO i2=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(i2))) /MKS%Jy/1.E6 * &
                       lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                         parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR
        Const1D(:) = EXP(parval(indpar%lnFstar)) /MKS%Jy/1.E6 * &
                     pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO i3=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(i3))) /MKS%Jy/1.E6 * &
                       gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                         parval(indpar%Wline(i3)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-EXP(parval(indpar%lnAv))/1.086_DP * extinction(iw))
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      END IF
    END DO loopCONT
    
    !! 2. Bands
    !!----------
    loopBAND: DO i=1,Nband
      gridlnIband = TRIMEQ(parname, 'lnIband'//TRIMLR(PRING(i)))
      gridCband = TRIMEQ(parname, 'Cband'//TRIMLR(PRING(i)))
      gridWSband = TRIMEQ(parname, 'WSband'//TRIMLR(PRING(i)))
      gridWLband = TRIMEQ(parname, 'WLband'//TRIMLR(PRING(i)))

      IF (gridlnIband) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO i1=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(i1))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(i1))), Qabs(i1))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND *lnIband
        Const1D(:) = 1._DP /MKS%Jy/1.E6 * &
                     lorentzBand(wvl(:), parval(indpar%Cband(i)), &
                       parval(indpar%WSband(i)), parval(indpar%WLband(i)), .TRUE.)
        FORALL (igrid=1:Ngrid) &
          FnuBAND0(igrid,:) = EXP(parvec(igrid)) * Const1D(:)
        !! BAND
        Const1D(:) = 0._DP
        DO i2=1,Nband
          IF (i2.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(i2))) /MKS%Jy/1.E6 * &
                         lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                           parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = FnuBAND0(:,iw) + Const1D(iw)
        !! STAR
        Const1D(:) = EXP(parval(indpar%lnFstar)) /MKS%Jy/1.E6 * &
                     pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO i3=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(i3))) /MKS%Jy/1.E6 * &
                       gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                         parval(indpar%Wline(i3)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-EXP(parval(indpar%lnAv))/1.086_DP * extinction(iw))
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      ELSE IF (gridCband) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO i1=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(i1))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(i1))), Qabs(i1))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND *Cband
        Const = EXP(parval(indpar%lnIband(i))) /MKS%Jy/1.E6
        FORALL (igrid=i:Ngrid) &
          FnuBAND0(igrid,:) = Const * lorentzBand(wvl(:), parvec(igrid), &
                                        parval(indpar%WSband(i)), parval(indpar%WLband(i)), .TRUE.)
        !! BAND
        Const1D(:) = 0._DP
        DO i2=1,Nband
          IF (i2.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(i2))) /MKS%Jy/1.E6 * &
                         lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                           parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = FnuBAND0(:,iw) + Const1D(iw)
        !! STAR
        Const1D(:) = EXP(parval(indpar%lnFstar)) /MKS%Jy/1.E6 * &
                     pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO i3=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(i3))) /MKS%Jy/1.E6 * &
                       gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                         parval(indpar%Wline(i3)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-EXP(parval(indpar%lnAv))/1.086_DP * extinction(iw))
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      ELSE IF (gridWSband) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO i1=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(i1))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(i1))), Qabs(i1))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND *WSband
        Const = EXP(parval(indpar%lnIband(i))) /MKS%Jy/1.E6
        FORALL (igrid=i:Ngrid) &
          FnuBAND0(igrid,:) = Const * lorentzBand(wvl(:), parval(indpar%Cband(i)), &
                                        parvec(igrid), parval(indpar%WLband(i)), .TRUE.)
        !! BAND
        Const1D(:) = 0._DP
        DO i2=1,Nband
          IF (i2.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(i2))) /MKS%Jy/1.E6 * &
                         lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                           parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = FnuBAND0(:,iw) + Const1D(iw)
        !! STAR
        Const1D(:) = EXP(parval(indpar%lnFstar)) /MKS%Jy/1.E6 * &
                     pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO i3=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(i3))) /MKS%Jy/1.E6 * &
                       gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                         parval(indpar%Wline(i3)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-EXP(parval(indpar%lnAv))/1.086_DP * extinction(iw))
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)        

      ELSE IF (gridWLband) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO i1=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(i1))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(i1))), Qabs(i1))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND *WLband
        Const = EXP(parval(indpar%lnIband(i))) /MKS%Jy/1.E6
        FORALL (igrid=i:Ngrid) &
          FnuBAND0(igrid,:) = Const * lorentzBand(wvl(:), parval(indpar%Cband(i)), &
                                       parval(indpar%WSband(i)), parvec(igrid), .TRUE.)
        !! BAND
        Const1D(:) = 0._DP
        DO i2=1,Nband
          IF (i2.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(i2))) /MKS%Jy/1.E6 * &
                         lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                           parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = FnuBAND0(:,iw) + Const1D(iw)
        !! STAR
        Const1D(:) = EXP(parval(indpar%lnFstar)) /MKS%Jy/1.E6 * &
                     pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO i3=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(i3))) /MKS%Jy/1.E6 * &
                       gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                         parval(indpar%Wline(i3)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-EXP(parval(indpar%lnAv))/1.086_DP * extinction(iw))
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      END IF
    END DO loopBAND
    
    !! 3. Stellar Continuum
    !!----------------------
    gridlnFstar = TRIMEQ(parname, 'lnFstar')
    
    IF (gridlnFstar) THEN
      !! CONT
      Const1D(:) = 0._DP
      DO i1=1,Ncont
        Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(i1))) * MKS%Msun/MKS%pc**2 * &
                     modifBB(wvl(:), EXP(parval(indpar%lnT(i1))), Qabs(i1))

      END DO
      FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
      !! BAND
      Const1D(:) = 0._DP
      DO i2=1,Nband
        Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(i2))) /MKS%Jy/1.E6 * &
                     lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                       parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
        
      END DO
      FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
      !! STAR *lnFstar
      FORALL (igrid=1:Ngrid) &
        FnuSTAR0(igrid,:) = EXP(parvec(igrid)) /MKS%Jy/1.E6 * &
                            pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
      !! LINE
      Const1D(:) = 0._DP
      DO i3=1,Nline
        Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(i3))) /MKS%Jy/1.E6 * &
                     gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                       parval(indpar%Wline(i3)), .TRUE.)

      END DO
      FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
      !! PABS
      FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-EXP(parval(indpar%lnAv))/1.086_DP * extinction(iw))
      !! Total model
      specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                           Pabs0(:,:) + FnuLINE0(:,:)

    END IF
    
    !! 4. Screen extinction
    !!----------------------
    gridlnAv = TRIMEQ(parname, 'lnAv')
    
    IF (gridlnAv) THEN
      !! CONT
      Const1D(:) = 0._DP
      DO i1=1,Ncont
        Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(i1))) * MKS%Msun/MKS%pc**2 * &
                     modifBB(wvl(:), EXP(parval(indpar%lnT(i1))), Qabs(i1))

      END DO
      FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
      !! BAND
      Const1D(:) = 0._DP
      DO i2=1,Nband
        Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(i2))) /MKS%Jy/1.E6 * &
                     lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                       parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
        
      END DO
      FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
      !! STAR
      Const1D(:) = EXP(parval(indpar%lnFstar)) /MKS%Jy/1.E6 * &
                   pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
      FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
      !! LINE
      Const1D(:) = 0._DP
      DO i3=1,Nline
        Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(i3))) /MKS%Jy/1.E6 * &
                     gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                       parval(indpar%Wline(i3)), .TRUE.)

      END DO
      FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
      !! PABS *lnAv
      FORALL (igrid=1:Ngrid) &
        Pabs0(igrid,:) = EXP(-EXP(parvec(igrid))/1.086_DP * extinction(:))
      !! Total model
      specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                           Pabs0(:,:) + FnuLINE0(:,:)

    END IF
    
    !! 5. Lines
    !!----------
    loopLINE: DO i=1,Nline
      gridlnIline = TRIMEQ(parname, 'lnIline'//TRIMLR(PRING(i)))
      gridCline = TRIMEQ(parname, 'Cline'//TRIMLR(PRING(i)))
      gridWline = TRIMEQ(parname, 'Wline'//TRIMLR(PRING(i)))
      
      IF (gridlnIline) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO i1=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(i1))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(i1))), Qabs(i1))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO i2=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(i2))) /MKS%Jy/1.E6 * &
                       lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                         parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR
        Const1D(:) = EXP(parval(indpar%lnFstar)) /MKS%Jy/1.E6 * &
                     pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE *lnIline
        Const1D(:) =  1._DP /MKS%Jy/1.E6 * &
                      gaussLine(wvl(:), parval(indpar%Cline(i)), &
                        parval(indpar%Wline(i)), .TRUE.)
        FORALL (igrid=1:Ngrid) &
          FnuLINE0(igrid,:) = EXP(parvec(igrid)) * Const1D(:)
        !! LINE
        Const1D(:) = 0._DP
        DO i3=1,Nline
          IF (i3.NE.i) &
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(i3))) /MKS%Jy/1.E6 * &
                       gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                         parval(indpar%Wline(i3)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = FnuLINE0(:,iw) + Const1D(iw)
        !! PABS
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-EXP(parval(indpar%lnAv))/1.086_DP * extinction(iw))
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      ELSE IF (gridCline) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO i1=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(i1))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(i1))), Qabs(i1))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO i2=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(i2))) /MKS%Jy/1.E6 * &
                       lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                         parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR
        Const1D(:) = EXP(parval(indpar%lnFstar)) /MKS%Jy/1.E6 * &
                     pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE *Cline
        Const = EXP(parval(indpar%lnIline(i))) /MKS%Jy/1.E6
        FORALL (igrid=1:Ngrid) &
          FnuBAND0(igrid,:) = Const * gaussLine(wvl(:), parvec(igrid), &
                                       parval(indpar%Wline(i)), .TRUE.)
        !! LINE
        Const1D(:) = 0._DP
        DO i3=1,Nline
          IF (i3.NE.i) &
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(i3))) /MKS%Jy/1.E6 * &
                       gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                         parval(indpar%Wline(i3)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = FnuLINE0(:,iw) + Const1D(iw)
        !! PABS
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-EXP(parval(indpar%lnAv))/1.086_DP * extinction(iw))
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      ELSE IF (gridWline) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO i1=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(i1))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(i1))), Qabs(i1))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO i2=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(i2))) /MKS%Jy/1.E6 * &
                       lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                         parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR
        Const1D(:) = EXP(parval(indpar%lnFstar)) /MKS%Jy/1.E6 * &
                     pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE *Wline
        Const = EXP(parval(indpar%lnIline(i))) /MKS%Jy/1.E6
        FORALL (igrid=1:Ngrid) &
          FnuBAND0(igrid,:) = Const * gaussLine(wvl(:), parvec(igrid), &
                                       parval(indpar%Wline(i)), .TRUE.)
        !! LINE
        Const1D(:) = 0._DP
        DO i3=1,Nline
          IF (i3.NE.i) &
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(i3))) /MKS%Jy/1.E6 * &
                       gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                         parval(indpar%Wline(i3)), .TRUE.)
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = FnuLINE0(:,iw) + Const1D(iw)
        !! PABS
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-EXP(parval(indpar%lnAv))/1.086_DP * extinction(iw))
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      END IF
    END DO loopLINE

  END FUNCTION specModel_gen

  FUNCTION specModel_scl( wvl, parinfo, indpar, parval, Qabs )
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN)           :: wvl
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo
    TYPE(indpar_type), INTENT(IN)                :: indpar
    REAL(DP), DIMENSION(:), INTENT(IN)           :: parval ! 1D
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN)    :: Qabs
    REAL(DP), DIMENSION(SIZE(wvl))               :: specModel_scl
    specModel_scl(:) = RESHAPE( specModel_gen(wvl(:), PARVEC=[0._DP],PARNAME='', &
                                  PARINFO=parinfo, INDPAR=indpar, PARVAL=parval(:), &
                                  QABS=Qabs(:)), [SIZE(wvl(:))] )
  END FUNCTION specModel_scl
  
END MODULE auxil
