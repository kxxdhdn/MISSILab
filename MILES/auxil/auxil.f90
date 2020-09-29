MODULE auxil

  USE utilities, ONLY: DP, STRIKE, WARNING
  USE constants, ONLY: 
  USE inout, ONLY: lenpar
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: NparCONT = 2, NparLINE = 3, NparBAND = 4 ! set_indpar
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
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Mcont ! dust mass (Modified BlackBody) [Msun]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Tcont ! dust temperature (MBB) [K]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Iline ! line Intensity [W/m2/sr]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Cline ! Center [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Wline ! Width [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Iband ! band Intensity [W/m2/sr]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Cband ! Center (pic) [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WSband ! Width Short nu side [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WLband ! Width Long nu side [um]
    REAL(DP) :: Av ! Attenuation in the V band [mag]
    REAL(DP) :: Fstar ! Total star luminosity [Lsun/pc2]
  END TYPE par_type
 
  TYPE, PUBLIC :: parinfo_type
    CHARACTER(lenpar) :: name = ""
    CHARACTER(lenpar) :: comp = ""
    REAL(DP) :: value = 0._DP
    REAL(DP), DIMENSION(2) :: limits = [0._DP, 0._DP]
    LOGICAL, DIMENSION(2) :: limited = [.FALSE., .FALSE.]
    LOGICAL :: fixed = .TRUE.
    LOGICAL :: hyper = .FALSE.
    LOGICAL :: model = .TRUE.
    CHARACTER(lenpar) :: tied = ""
    INTEGER :: itied = 0
    REAL(DP) :: mean = 0._DP
    REAL(DP) :: sigma = 0._DP
    INTEGER :: ind = -1
  END TYPE parinfo_type

  TYPE, PUBLIC :: indpar_type
    INTEGER, DIMENSION(Ncont_max) :: Mcont = -1
    INTEGER, DIMENSION(Ncont_max) :: Tcont = -1
    INTEGER, DIMENSION(Nline_max) :: Iline = -1
    INTEGER, DIMENSION(Nline_max) :: Cline = -1
    INTEGER, DIMENSION(Nline_max) :: Wline = -1
    INTEGER, DIMENSION(Nband_max) :: Iband = -1
    INTEGER, DIMENSION(Nband_max) :: Cband = -1
    INTEGER, DIMENSION(Nband_max) :: WSband = -1
    INTEGER, DIMENSION(Nband_max) :: WLband = -1
    INTEGER :: Av = -1
    INTEGER :: Fstar = -1
    INTEGER, DIMENSION(:), ALLOCATABLE :: extra
  END TYPE indpar_type

  TYPE, PUBLIC :: Qabs_type
    !! rho(kg/m3); wave(micron); nu(Hz); Qova(m-2)
    REAL(DP)                            :: rho
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave, nu
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Qova ! a0=1.E-2_DP micron
    REAL(DP), DIMENSION(:), ALLOCATABLE :: coeffMBB
  END TYPE Qabs_type
 
CONTAINS

  !!-------------------------------------------------------
  !!
  !! Fill the INDPAR_TYPE structure, from a PARINFO_TYPE structure
  !!
  !!-------------------------------------------------------
  SUBROUTINE set_indpar(indpar, parinfo)
    
    USE utilities, ONLY: trimeq, trimlr, pring
    USE arrays, ONLY: iwhere
    IMPLICIT NONE

    TYPE(indpar_type), INTENT(OUT) :: indpar
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo

    INTEGER :: i, Ncont, Nline, Nband, Nextra
    LOGICAL :: CONTset, LINEset, BANDset, PABSset, STARset, extraset

    CONTset = ( ANY(TRIMEQ(parinfo(:)%comp, "CONT")) )
    LINEset = ( ANY(TRIMEQ(parinfo(:)%comp, "LINE")) )
    BANDset = ( ANY(TRIMEQ(parinfo(:)%comp, "BAND")) )
    PABSset = ( ANY(TRIMEQ(parinfo(:)%comp, "PABS")) )
    STARset = ( ANY(TRIMEQ(parinfo(:)%comp, "STAR")) )
    extraset = ( ANY(TRIMEQ(parinfo(:)%comp, "extra")) )

    IF (CONTset) THEN
      Ncont = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, "CONT")) ) / NparCONT
      DO i=1,Ncont
        CALL IWHERE( TRIMEQ(parinfo%name(:), "Mcont"//TRIMLR(PRING(i))), indpar%Mcont(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), "Tcont"//TRIMLR(PRING(i))), indpar%Tcont(i) )
        
      END DO
    END IF
    
    IF (LINEset) THEN
      Nline = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, "LINE")) ) / NparLINE
      DO i=1,Nline
        CALL IWHERE( TRIMEQ(parinfo%name(:), "Iline"//TRIMLR(PRING(i))), indpar%Iline(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), "Cline"//TRIMLR(PRING(i))), indpar%Cline(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), "Wline"//TRIMLR(PRING(i))), indpar%Wline(i) )

      END DO
    END IF

    IF (BANDset) THEN
      Nband = SIZE( PACK(parinfo(:), TRIMEQ(parinfo(:)%comp, "BAND")) ) / NparBAND
      DO i=1,Nband
        CALL IWHERE( TRIMEQ(parinfo%name(:), "Iband"//TRIMLR(PRING(i))), indpar%Iband(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), "Cband"//TRIMLR(PRING(i))), indpar%Cband(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), "WSband"//TRIMLR(PRING(i))), indpar%WSband(i) )
        CALL IWHERE( TRIMEQ(parinfo%name(:), "WLband"//TRIMLR(PRING(i))), indpar%WLband(i) )
        
      END DO
    END IF

    IF (PABSset) THEN
      CALL IWHERE(TRIMEQ(parinfo%name(:),"Av"), indpar%Av)
      
    END IF
    
    IF (STARset) THEN
      CALL IWHERE(TRIMEQ(parinfo%name(:),"Fstar"),indpar%Fstar)
      
    END IF

    IF (extraset) THEN
      Nextra = COUNT(TRIMEQ(parinfo(:)%comp, "extra"))
      IF (Nextra > 0) CALL IWHERE(TRIMEQ(parinfo(:)%comp,"extra"),indpar%extra)
      
    END IF
    
  END SUBROUTINE set_indpar

  !!-------------------------------------------------------
  !!
  !! Read the input master file for the Chi2/Bayesian run
  !!
  !!-------------------------------------------------------
  SUBROUTINE read_master(labQ, labBIN, labLIN, waveall, &
                         Qabs, indBIN, indLIN)
                        ! Nmcmc, verbose)

    USE datable, ONLY: BIN, LIN
    USE utilities, ONLY: DP, strike, trimeq, trimlr, pring
    USE inout, ONLY: lenpar, lenline
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: closest, iwhere, reallocate
    USE interpolation, ONLY: interp_lin_sorted
    USE grain_optics, ONLY: rho_grain, read_optics
    IMPLICIT NONE

    !! Temporary inputs
    CHARACTER(*), DIMENSION(:), INTENT(IN)               :: labQ
    CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL     :: labBIN, labLIN
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL         :: waveall
    
    ! INTEGER, PARAMETER :: unit = 1
    INTEGER :: i!, iostat
    ! CHARACTER(lenline) :: inpn, inpv
    ! CHARACTER(lenline), DIMENSION(:), ALLOCATABLE :: inpname, inpval

    INTEGER                                              :: Nspec, Nr, Nw, Nband, Nline
    INTEGER, DIMENSION(SIZE(labQ))                       :: indr ! radius index
    REAL(DP), PARAMETER                                  :: a0 = 1.E-2_DP ! grain radius
    REAL(DP), DIMENSION(:), ALLOCATABLE                  :: wave0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE                :: radius
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE              :: Q ! (Nspec, Nr, Nw)
    
    ! INTEGER, INTENT(OUT), OPTIONAL :: Nmcmc
    ! LOGICAL, INTENT(OUT), OPTIONAL :: verbose
    TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Qabs
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL         :: indBIN
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL         :: indLIN

    !!----------------------------
    !! Read the input master file
    !!----------------------------
    ! OPEN (unit, FILE="./input_fitMIR_master.txt", &
    !       STATUS="OLD", ACTION="READ", IOSTAT=iostat)
    !   DO 
    !     IF (iostat /= 0) EXIT
    !     CALL READ_INPUT_LINE(unit,inpn,inpv,IOSTAT=iostat)
    !     CALL INCRARR(inpname,inpn)
    !     CALL INCRARR(inpval,inpv)
        
    !   END DO
    ! CLOSE (unit)

    !!--------------------------------
    !! Read optical properties (Qabs)
    !!--------------------------------
    IF (PRESENT(Qabs)) THEN
      !! Qabs % wave [um], nu [Hz], Qova [m^-1], rho [kg/m^3]
      CALL read_optics(labQ, WAVE=wave0, RADIUSALL=radius, QABSALL=Q)    
      
      Nspec = SIZE(labQ)
      Nw = SIZE(wave0)

      ALLOCATE(Qabs(SIZE(labQ)))
      
      DO i=1,Nspec
        
        ALLOCATE(Qabs(i)%Qova(Nw), Qabs(i)%coeffMBB(Nw))
        
        Nr = SIZE(radius(i,:))
        Qabs(i)%rho = rho_grain(labQ(i))
        Qabs(i)%wave = wave0
        Qabs(i)%nu = MKS%clight/MKS%micron / wave0
        indr(i) = closest(radius(i,:), a0)
        Qabs(i)%Qova = Q(i,indr(i),:)/radius(i,indr(i))
        ! PRINT*, 'Radius of ', labQ(i), ': ', radius(i,indr(i)), '->', indr(i)
        Qabs(i)%coeffMBB = 3._DP*pi/4._DP/rho_grain(labQ(i)) * Q(i,indr(i),:)/radius(i,indr(i))
        
        !! [Optional] Interpolate Qabs to input wave grid
        IF (PRESENT(waveall)) THEN
          Qabs(i)%wave = waveall
          Qabs(i)%nu = MKS%clight/MKS%micron / waveall
          Qabs(i)%Qova = interp_lin_sorted(Qabs(i)%Qova, wave0, waveall, &
                                            XLOG=.TRUE., YLOG=.TRUE., FORCE=.TRUE.)
          Qabs(i)%coeffMBB = interp_lin_sorted(Qabs(i)%coeffMBB, wave0, waveall, &
                                                XLOG=.TRUE., YLOG=.TRUE., FORCE=.TRUE.)
          
        END IF
      END DO
      !! Free memory space
      DEALLOCATE(wave0, radius, Q)
      
    END IF

    !!----------------------
    !! Read bands and lines
    !!----------------------
    IF (PRESENT(indBIN)) THEN
      IF (PRESENT(labBIN)) THEN
        Nband = SIZE(labBIN)
        ALLOCATE(indBIN(Nband))
        DO i=1,Nband
          CALL IWHERE( TRIMEQ(BIN(:)%label, labBIN(i)), indBIN(i) )
          
        END DO
      ELSE
        Nband = Nband_max
        ALLOCATE(indBIN(Nband))
        indBIN(:) = (/ (i, i=1,Nband) /)
      END IF
    END IF

    IF (PRESENT(indLIN)) THEN
      IF (PRESENT(labLIN)) THEN
        Nline = SIZE(labLIN)
        ALLOCATE(indLIN(Nline))
        DO i=1,Nline
          CALL IWHERE( TRIMEQ(LIN(:)%label, labLIN(i)), indLIN(i) )
          
        END DO
      ELSE
        Nline = Nline_max
        ALLOCATE(indLIN(Nline))
        indLIN(:) = (/ (i, i=1,Nline) /)
      END IF
    END IF

  END SUBROUTINE read_master

  !!-------------------------------------------------------
  !!
  !!     Automatic initialization of model parameters
  !!
  !!-------------------------------------------------------
  ! SUBROUTINE initparam(NiniMC, par, mask)
    
  !   USE utilities, ONLY: DP, trimeq, trimlr, pring, isNaN
  !   IMPLICIT NONE

  !   INTEGER, INTENT(IN) :: NiniMC
  !   LOGICAL, DIMENSION(:,:,:), INTENT(IN) :: mask
  !   REAL(DP), DIMENSION(:,:,:,:), INTENT(INOUT) :: par
  !   INTEGER :: i, j, ipar, x, y, Nx, Ny, Npar, Nextra, Nmc
  !   LOGICAL, DIMENSION(:,:), ALLOCATABLE :: maskxy
    
  !   !! Preliminary
  !   Nx = SIZE(par(:,:,:,:),1)
  !   Ny = SIZE(par(:,:,:,:),2)
  !   ALLOCATE(maskxy(Nx,Ny))
  !   FORALL (x=1:Nx,y=1:Ny) maskxy(x,y) = ANY(mask(x,y,:))

  ! END SUBROUTINE initparam
  
  !!-------------------------------------------------------
  !!
  !! Automatize the degradation of the spectral resolution
  !!
  !!-------------------------------------------------------
  FUNCTION degradeRes(wc, dw, instr)

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
  PURE FUNCTION modifBB_gen(dblarr, temp, Qabs)
    !! temp [K]; blackbody [W/m2/sr/Hz]; 
    !! output array [MJy/sr]
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
      modifBB_gen(it,:) = Qabs(it)%coeffMBB * BLACKBODY(nu(:), temp(it)) ! in W/m2/sr/Hz
    !! Unit onversion ([W/m2/sr/Hz] to [MJy/sr])
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
  PURE FUNCTION gaussLine_gen(dblarr, ref, sig, w2nu)
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

  PURE FUNCTION gaussLine_0(dblarr, ref, sig, w2nu)
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


  PURE FUNCTION gaussLine_1(dblarr, ref, sig, w2nu)
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
  PURE FUNCTION lorentzBand_gen(dblarr, ref, sigS, sigL, w2nu)
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

  PURE FUNCTION lorentzBand_0(dblarr, ref, sigS, sigL, w2nu)
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

  PURE FUNCTION lorentzBand_1(dblarr, ref, sigS, sigL, w2nu)
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
  FUNCTION extCurve(dblarr, nu2w)
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
    
    CALL read_hdf5(DBLARR1D=wave, NAME="lambda (micron)", &
                   N1=Nw0, FILE="../data/extcurve"//h5ext)
    CALL read_hdf5(DBLARR1D=Cext, NAME="C_extovH (cm^2ovH)", &
                   FILE="../data/extcurve"//h5ext)
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
  FUNCTION specModel_3D(wvl, indpar, parval, Qabs, &
                        FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE)

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
    Nband = COUNT(indpar%Iband(:) .NE. -1)
    Tstar = 5.E4_DP ! high enough to stay in Rayleigh-Jeans limit [K]
    Nline = COUNT(indpar%Iline(:) .NE. -1)
    nu(:) = MKS%clight/MKS%micron / wvl(:)
    extinction(:) = extCurve(wvl(:))

    !! Initialization
    specModel_3D(:,:,:) = 0._DP
    
    !! 1. Continuum
    !!--------------
    CALL REALLOCATE(Const4D,Nx,Ny,Nw,Ncont)
    FORALL (x=1:Nx,y=1:Ny,i=1:Ncont) &
      Const4D(x,y,:,i) = parval(x,y,indpar%Mcont(i)) * MKS%Msun/MKS%pc**2 * &
                         modifBB(wvl(:), parval(x,y,indpar%Tcont(i)), Qabs(i))
    FnuCONT0(:,:,:) = SUM(Const4D(:,:,:,:),DIM=4)

    IF (PRESENT(FnuCONT)) THEN
      ALLOCATE(FnuCONT(Nx,Ny,Nw))
      FnuCONT(:,:,:) = FnuCONT0(:,:,:)
    END IF
    
    !! 2. Bands
    !!----------
    CALL REALLOCATE(Const4D,Nx,Ny,Nw,Nband)
    FORALL (x=1:Nx,y=1:Ny,i=1:Nband) &
      Const4D(x,y,:,i) = parval(x,y,indpar%Iband(i)) /MKS%Jy/1.E6 * &
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
      FnuSTAR0(x,y,:) = parval(x,y,indpar%Fstar) * MKS%Lsun/4._DP/pi/MKS%pc**2 * &
                         pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4 /MKS%Jy/1.E6

    IF (PRESENT(FnuSTAR)) THEN
      ALLOCATE(FnuSTAR(Nx,Ny,Nw))
      FnuSTAR(:,:,:) = FnuSTAR0(:,:,:)
    END IF
    
    !! 4. Screen extinction
    !!----------------------
    FORALL (x=1:Nx,y=1:Ny) &
      Pabs0(x,y,:) = EXP(-parval(x,y,indpar%Av)/1.086_DP * extinction(:))

    IF (PRESENT(Pabs)) THEN
      ALLOCATE(Pabs(Nx,Ny,Nw))
      Pabs(:,:,:) = Pabs0(:,:,:)
    END IF
    
    !! 5. Lines
    !!----------
    CALL REALLOCATE(Const4D,Nx,Ny,Nw,Nline)
    FORALL (x=1:Nx,y=1:Ny,i=1:Nline) &
      Const4D(x,y,:,i) = parval(x,y,indpar%Iline(i)) /MKS%Jy/1.E6 * &
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
  FUNCTION specModel_2D(wvl, indpar, parval, Qabs)

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
    Nband = COUNT(indpar%Iband(:) .NE. -1)
    Tstar = 5.E4_DP ! high enough to stay in Rayleigh-Jeans limit [K]
    Nline = COUNT(indpar%Iline(:) .NE. -1)
    nu(:) = MKS%clight/MKS%micron / wvl(:)
    extinction(:) = extCurve(wvl(:))

    !! Initialization
    specModel_2D(:,:) = 0._DP
    
    !! 1. Continuum
    !!--------------
    CALL REALLOCATE(Const3D,Npar,Nw,Ncont)
    FORALL (ipar=1:Npar, i=1:Ncont) &
      Const3D(ipar,:,i) = parval(indpar%Mcont(i),ipar) * MKS%Msun/MKS%pc**2 * &
                          modifBB(wvl(:), parval(indpar%Tcont(i),ipar), Qabs(i))
    FnuCONT0(:,:) = SUM(Const3D(:,:,:),DIM=3)
    
    !! 2. Bands
    !!----------
    CALL REALLOCATE(Const3D,Npar,Nw,Nband)
    FORALL (ipar=1:Npar, i=1:Nband) &
      Const3D(ipar,:,i) = parval(indpar%Iband(i),ipar) /MKS%Jy/1.E6 * &
                          lorentzBand(wvl(:), parval(indpar%Cband(i),ipar), &
                            parval(indpar%WSband(i),ipar), parval(indpar%WLband(i),ipar), .TRUE.)
    FnuBAND0(:,:) = SUM(Const3D(:,:,:),DIM=3)
    
    !! 3. Stellar Continuum
    !!----------------------
    FORALL (ipar=1:Npar) &
      FnuSTAR0(ipar,:) = parval(indpar%Fstar,ipar) * MKS%Lsun/4._DP/pi/MKS%pc**2 * &
                         pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4 /MKS%Jy/1.E6

    !! 4. Screen extinction
    !!----------------------
    FORALL (ipar=1:Npar) &
      Pabs0(ipar,:) = EXP(-parval(indpar%Av,ipar)/1.086_DP * extinction(:))

    !! 5. Lines
    !!----------
    CALL REALLOCATE(Const3D,Npar,Nw,Nline)
    FORALL (ipar=1:Npar, i=1:Nline) &
      Const3D(ipar,:,i) = parval(indpar%Iline(i),ipar) /MKS%Jy/1.E6 * &
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
  ! FUNCTION specModel_1D(wvl, indpar, parval, Qabs, &
  !                       FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE)
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
  
  FUNCTION specModel_1D(wvl, indpar, parval, Qabs, &
                        FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE)

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
    Nband = COUNT(indpar%Iband(:) .NE. -1)
    Tstar = 5.E4_DP ! high enough to stay in Rayleigh-Jeans limit [K]
    Nline = COUNT(indpar%Iline(:) .NE. -1)
    nu(:) = MKS%clight/MKS%micron / wvl(:)
    extinction(:) = extCurve(wvl(:))

    !! Initialization
    specModel_1D(:) = 0._DP
  
    !! 1. Continuum
    !!--------------
    CALL REALLOCATE(Const2D,Nw,Ncont)
    FORALL (i=1:Ncont) &
      Const2D(:,i) = parval(indpar%Mcont(i)) * MKS%Msun/MKS%pc**2 * &
                     modifBB(wvl(:), parval(indpar%Tcont(i)), Qabs(i))
    FnuCONT0(:) = SUM(Const2D,DIM=2)

    IF (PRESENT(FnuCONT)) THEN
      ALLOCATE(FnuCONT(Nw))
      FnuCONT(:) = FnuCONT0(:)
    END IF
  
    !! 2. Bands
    !!----------
    CALL REALLOCATE(Const2D,Nw,Nband)
    FORALL (i=1:Nband) &
      Const2D(:,i) = parval(indpar%Iband(i)) /MKS%Jy/1.E6 * &
                     lorentzBand(wvl(:), parval(indpar%Cband(i)), &
                       parval(indpar%WSband(i)), parval(indpar%WLband(i)), .TRUE.)
    FnuBAND0(:) = SUM(Const2D,DIM=2)

    IF (PRESENT(FnuBAND)) THEN
      ALLOCATE(FnuBAND(Nw))
      FnuBAND(:) = FnuBAND0(:)
    END IF
  
    !! 3. Stellar Continuum
    !!----------------------
    FnuSTAR0(:) = parval(indpar%Fstar) * MKS%Lsun/4._DP/pi/MKS%pc**2 * &
                  pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4 /MKS%Jy/1.E6

    IF (PRESENT(FnuSTAR)) THEN
      ALLOCATE(FnuSTAR(Nw))
      FnuSTAR(:) = FnuSTAR0(:)
    END IF

    !! 4. Screen extinction
    !!----------------------
    Pabs0(:) = EXP(-parval(indpar%Av)/1.086_DP * extinction(:))

    IF (PRESENT(Pabs)) THEN
      ALLOCATE(Pabs(Nw))
      Pabs(:) = Pabs0(:)
    END IF

    !! 5. Lines
    !!----------
    CALL REALLOCATE(Const2D,Nw,Nline)
    FORALL (i=1:Nline) &
      Const2D(:,i) = parval(indpar%Iline(i)) /MKS%Jy/1.E6 * &
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
  FUNCTION specModel_gen(wvl, parvec, parname, parinfo, indpar, &
                         parval, Qabs, indBIN, indLIN)
    !! Nstar = 0 or 1
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
    INTEGER, DIMENSION(:), INTENT(IN)            :: indBIN, indLIN
    INTEGER :: Nw, Npar, Ncont, Nline, Nband
    INTEGER :: i, ipar, iw, i1, i2, i3!, itied
    REAL(DP)                                     :: Tstar
    REAL(DP), DIMENSION(SIZE(wvl))               :: nu, extinction
    LOGICAL :: CONTset, LINEset, BANDset, PABSset, STARset
    LOGICAL :: gridMcont, gridTcont, gridIline, gridCline, gridWline
    LOGICAL :: gridIband, gridCband, gridWSband, gridWLband, gridAv, gridFstar
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
    Npar = SIZE(parvec(:))
    Ncont = SIZE(Qabs(:))
    Nband = SIZE(indBIN(:))
    Tstar = 5.E4_DP ! high enough to stay in Rayleigh-Jeans limit [K]
    Nline = SIZE(indLIN(:))
    nu(:) = MKS%clight/MKS%micron / wvl(:)
    extinction(:) = extCurve(wvl(:))

    !! Component settings
    CONTset = ANY(TRIMEQ(parinfo(:)%comp, "CONT"))
    LINEset = ANY(TRIMEQ(parinfo(:)%comp, "LINE"))
    BANDset = ANY(TRIMEQ(parinfo(:)%comp, "BAND"))
    PABSset = ANY(TRIMEQ(parinfo(:)%comp, "PABS"))
    STARset = ANY(TRIMEQ(parinfo(:)%comp, "STAR"))

    !! Is the parameter tied to another one?
    ! tied = ANY(TRIMEQ(parinfo(:)%tied, parname))
    ! IF (tied) CALL IWHERE( TRIMEQ(parinfo(:)%tied, parname), itied )

    !! Initialization
    specModel_gen(:,:) = 0._DP
    
    !! 1. Continuum
    !!--------------
    loopCONT: DO i=1,Ncont
      gridMcont = TRIMEQ(parname, "Mcont"//TRIMLR(PRING(i)))
      gridTcont = TRIMEQ(parname, "Tcont"//TRIMLR(PRING(i)))
      IF (gridMcont) THEN
      !! CONT *Mcont
        Const1D(:) = MKS%Msun/MKS%pc**2 * modifBB(wvl(:), parval(indpar%Tcont(i)), Qabs(i))
        FORALL (ipar=1:Npar) &
          FnuCONT0(ipar,:) = parvec(ipar) * Const1D(:)
      ELSE IF (gridTcont) THEN
      !! CONT *Tcont
        Const = parval(indpar%Mcont(i)) * MKS%Msun/MKS%pc**2
        FORALL (ipar=1:Npar) &
             FnuCONT0(ipar,:) = Const * modifBB(wvl(:), parvec(ipar), Qabs(i))

      END IF
      !! CONT
      Const1D(:) = 0._DP
      DO i1=1,Ncont
        IF (i1.NE.i) &
          Const1D(:) = Const1D(:) + parval(indpar%Mcont(i1)) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), parval(indpar%Tcont(i1)), Qabs(i1))

      END DO
      FORALL (iw=1:Nw) FnuCONT0(:,iw) = FnuCONT0(:,iw) + Const1D(iw)
      !! BAND
      Const1D(:) = 0._DP
      DO i2=1,Nband
        Const1D(:) = Const1D(:) + parval(indpar%Iband(i2)) /MKS%Jy/1.E6 * &
                     lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                       parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
        
      END DO
      FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
      !! STAR
      Const1D(:) = parval(indpar%Fstar) * MKS%Lsun/4._DP/pi/MKS%pc**2 * &
                   pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4 /MKS%Jy/1.E6
      FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
      !! LINE
      Const1D(:) = 0._DP
      DO i3=1,Nline
        Const1D(:) = Const1D(:) + parval(indpar%Iline(i3)) /MKS%Jy/1.E6 * &
                     gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                       parval(indpar%Wline(i3)), .TRUE.)
      
      END DO
      FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
      !! PABS
      FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-parval(indpar%Av)/1.086_DP * extinction(iw))
      !! Total model
      specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                           Pabs0(:,:) + FnuLINE0(:,:)

    END DO loopCONT
    
    !! 2. Bands
    !!----------
    loopBAND: DO i=1,Nband
      gridIband = TRIMEQ(parname, "Iband"//TRIMLR(PRING(i)))
      gridCband = TRIMEQ(parname, "Cband"//TRIMLR(PRING(i)))
      gridWSband = TRIMEQ(parname, "WSband"//TRIMLR(PRING(i)))
      gridWLband = TRIMEQ(parname, "WLband"//TRIMLR(PRING(i)))
      !! CONT
      Const1D(:) = 0._DP
      DO i1=1,Ncont
        Const1D(:) = Const1D(:) + parval(indpar%Mcont(i1)) * MKS%Msun/MKS%pc**2 * &
                     modifBB(wvl(:), parval(indpar%Tcont(i1)), Qabs(i1))
      
      END DO
      FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
      IF (gridIband) THEN
      !! BAND *Iband
        Const1D(:) = 1._DP /MKS%Jy/1.E6 * &
                     lorentzBand(wvl(:), parval(indpar%Cband(i)), &
                       parval(indpar%WSband(i)), parval(indpar%WLband(i)), .TRUE.)
        FORALL (ipar=1:Npar) &
          FnuBAND0(ipar,:) = parvec(ipar) * Const1D(:)
      ELSE IF (gridCband) THEN
      !! BAND *Cband
        Const = parval(indpar%Iband(i)) /MKS%Jy/1.E6
        FORALL (ipar=i:Npar) &
          FnuBAND0(ipar,:) = Const * lorentzBand(wvl(:), parvec(ipar), &
                                       parval(indpar%WSband(i)), parval(indpar%WLband(i)), .TRUE.)
      ELSE IF (gridWSband) THEN
      !! BAND *WSband
        Const = parval(indpar%Iband(i)) /MKS%Jy/1.E6
        FORALL (ipar=i:Npar) &
          FnuBAND0(ipar,:) = Const * lorentzBand(wvl(:), parval(indpar%Cband(i)), &
                                       parvec(ipar), parval(indpar%WLband(i)), .TRUE.)
      ELSE IF (gridWLband) THEN
      !! BAND *WLband
        Const = parval(indpar%Iband(i)) /MKS%Jy/1.E6
        FORALL (ipar=i:Npar) &
          FnuBAND0(ipar,:) = Const * lorentzBand(wvl(:), parval(indpar%Cband(i)), &
                                       parval(indpar%WSband(i)), parvec(ipar), .TRUE.)

      END IF
      !! BAND
      Const1D(:) = 0._DP
      DO i2=1,Nband
        IF (i2.NE.i) &
          Const1D(:) = Const1D(:) + parval(indpar%Iband(i2)) /MKS%Jy/1.E6 * &
                       lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                         parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)

      END DO
      FORALL (iw=1:Nw) FnuBAND0(:,iw) = FnuBAND0(:,iw) + Const1D(iw)
      !! STAR
      Const1D(:) = parval(indpar%Fstar) * MKS%Lsun/4._DP/pi/MKS%pc**2 * &
                   pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4 /MKS%Jy/1.E6
      FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
      !! LINE
      Const1D(:) = 0._DP
      DO i3=1,Nline
        Const1D(:) = Const1D(:) + parval(indpar%Iline(i3)) /MKS%Jy/1.E6 * &
                     gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                       parval(indpar%Wline(i3)), .TRUE.)
      
      END DO
      FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
      !! PABS
      FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-parval(indpar%Av)/1.086_DP * extinction(iw))
      !! Total model
      specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                           Pabs0(:,:) + FnuLINE0(:,:)

    END DO loopBAND
    
    !! 3. Stellar Continuum
    !!----------------------
    gridFstar = TRIMEQ(parname, "Fstar")
    IF (gridFstar) THEN
      !! CONT
      Const1D(:) = 0._DP
      DO i1=1,Ncont
        Const1D(:) = Const1D(:) + parval(indpar%Mcont(i1)) * MKS%Msun/MKS%pc**2 * &
                     modifBB(wvl(:), parval(indpar%Tcont(i1)), Qabs(i1))

      END DO
      FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
      !! BAND
      Const1D(:) = 0._DP
      DO i2=1,Nband
        Const1D(:) = Const1D(:) + parval(indpar%Iband(i2)) /MKS%Jy/1.E6 * &
                     lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                       parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
        
      END DO
      FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
      !! STAR *Fstar
      FORALL (ipar=1:Npar) &
        FnuSTAR0(ipar,:) = parvec(ipar) * MKS%Lsun/4._DP/pi/MKS%pc**2 * &
                           pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4 /MKS%Jy/1.E6
      !! LINE
      Const1D(:) = 0._DP
      DO i3=1,Nline
        Const1D(:) = Const1D(:) + parval(indpar%Iline(i3)) /MKS%Jy/1.E6 * &
                     gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                       parval(indpar%Wline(i3)), .TRUE.)

      END DO
      FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
      !! PABS
      FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-parval(indpar%Av)/1.086_DP * extinction(iw))
      !! Total model
      specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                           Pabs0(:,:) + FnuLINE0(:,:)

    END IF
    
    !! 4. Screen extinction
    !!----------------------
    gridAv = TRIMEQ(parname, "Av")
    IF (gridAv) THEN
      !! CONT
      Const1D(:) = 0._DP
      DO i1=1,Ncont
        Const1D(:) = Const1D(:) + parval(indpar%Mcont(i1)) * MKS%Msun/MKS%pc**2 * &
                     modifBB(wvl(:), parval(indpar%Tcont(i1)), Qabs(i1))

      END DO
      FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
      !! BAND
      Const1D(:) = 0._DP
      DO i2=1,Nband
        Const1D(:) = Const1D(:) + parval(indpar%Iband(i2)) /MKS%Jy/1.E6 * &
                     lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                       parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
        
      END DO
      FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
      !! STAR
      Const1D(:) = parval(indpar%Fstar) * MKS%Lsun/4._DP/pi/MKS%pc**2 * &
                   pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4 /MKS%Jy/1.E6
      FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
      !! LINE
      Const1D(:) = 0._DP
      DO i3=1,Nline
        Const1D(:) = Const1D(:) + parval(indpar%Iline(i3)) /MKS%Jy/1.E6 * &
                     gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                       parval(indpar%Wline(i3)), .TRUE.)

      END DO
      FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
      !! PABS *Av
      FORALL (ipar=1:Npar) &
        Pabs0(ipar,:) = EXP(-parvec(ipar)/1.086_DP * extinction(:))
      !! Total model
      specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                           Pabs0(:,:) + FnuLINE0(:,:)

    END IF
    
    !! 5. Lines
    !!----------
    loopLINE: DO i=1,Nline
      gridIline = TRIMEQ(parname, "Iline"//TRIMLR(PRING(i)))
      gridCline = TRIMEQ(parname, "Cline"//TRIMLR(PRING(i)))
      gridWline = TRIMEQ(parname, "Wline"//TRIMLR(PRING(i)))
      !! CONT
      Const1D(:) = 0._DP
      DO i1=1,Ncont
        Const1D(:) = Const1D(:) + parval(indpar%Mcont(i1)) * MKS%Msun/MKS%pc**2 * &
                     modifBB(wvl(:), parval(indpar%Tcont(i1)), Qabs(i1))
      
      END DO
      FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
      !! BAND
      Const1D(:) = 0._DP
      DO i2=1,Nband
        Const1D(:) = Const1D(:) + parval(indpar%Iband(i2)) /MKS%Jy/1.E6 * &
                     lorentzBand(wvl(:), parval(indpar%Cband(i2)), &
                       parval(indpar%WSband(i2)), parval(indpar%WLband(i2)), .TRUE.)
        
      END DO
      FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
      !! STAR
      Const1D(:) = parval(indpar%Fstar) * MKS%Lsun/4._DP/pi/MKS%pc**2 * &
                   pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4 /MKS%Jy/1.E6
      FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
      IF (gridIline) THEN
      !! LINE *Iline
        Const1D(:) =  1._DP /MKS%Jy/1.E6 * &
                      gaussLine(wvl(:), parval(indpar%Cline(i)), &
                        parval(indpar%Wline(i)), .TRUE.)
        FORALL (ipar=1:Npar) &
          FnuLINE0(ipar,:) = parvec(ipar) * Const1D(:)
      ELSE IF (gridCline) THEN
      !! LINE *Cline
        Const = parval(indpar%Iline(i)) /MKS%Jy/1.E6
        FORALL (ipar=1:Npar) &
          FnuBAND0(ipar,:) = Const * gaussLine(wvl(:), parvec(ipar), &
                                       parval(indpar%Wline(i)), .TRUE.)
      ELSE IF (gridWline) THEN
      !! LINE *Wline
        Const = parval(indpar%Iline(i)) /MKS%Jy/1.E6
        FORALL (ipar=1:Npar) &
          FnuBAND0(ipar,:) = Const * gaussLine(wvl(:), parvec(ipar), &
                                       parval(indpar%Wline(i)), .TRUE.)

      END IF
      !! LINE
      Const1D(:) = 0._DP
      DO i3=1,Nline
        IF (i3.NE.i) &
        Const1D(:) = Const1D(:) + parval(indpar%Iline(i3)) /MKS%Jy/1.E6 * &
                     gaussLine(wvl(:), parval(indpar%Cline(i3)), &
                       parval(indpar%Wline(i3)), .TRUE.)
      
      END DO
      FORALL (iw=1:Nw) FnuLINE0(:,iw) = FnuLINE0(:,iw) + Const1D(iw)
      !! PABS
      FORALL (iw=1:Nw) Pabs0(:,iw) = EXP(-parval(indpar%Av)/1.086_DP * extinction(iw))
      !! Total model
      specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                           Pabs0(:,:) + FnuLINE0(:,:)

    END DO loopLINE

  END FUNCTION specModel_gen

  FUNCTION specModel_scl(wvl, parinfo, indpar, parval, Qabs, indBIN, indLIN)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN)           :: wvl
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo
    TYPE(indpar_type), INTENT(IN)                :: indpar
    REAL(DP), DIMENSION(:), INTENT(IN)           :: parval ! 1D
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN)    :: Qabs
    INTEGER, DIMENSION(:), INTENT(IN)            :: indBIN, indLIN
    REAL(DP), DIMENSION(SIZE(wvl))               :: specModel_scl
    specModel_scl(:) = RESHAPE( specModel_gen(wvl(:), PARVEC=[0._DP],PARNAME="", &
                                  PARINFO=parinfo, INDPAR=indpar, PARVAL=parval(:), &
                                  QABS=Qabs(:), INDBIN=indBIN(:), INDLIN=indLIN(:)), [SIZE(wvl(:))] )
  END FUNCTION specModel_scl
  
END MODULE auxil
