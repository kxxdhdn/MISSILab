MODULE auxil

  USE utilities, ONLY: DP, STRIKE, WARNING
  USE constants, ONLY: pi, MKS
  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: par_str
    REAL(DP), DIMENSION(:), ALLOCATABLE :: massBB ! dust mass (Modified BlackBody) [Msun]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: tempBB ! dust temperature (MBB) [K]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Iline ! line Intensity [W/m2/Hz]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Cline ! Center [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Wline ! Width [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Iband ! band Intensity [W/m2/Hz]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Cband ! Center (pic) [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WSband ! Width Short nu side [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WLband ! Width Long nu side [um]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Av ! Attenuation in the V band [mag]
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Fstar ! Total star luminosity [Lsun]
  END TYPE par_str
 
  TYPE(par_str), SAVE, PUBLIC :: par

  TYPE, PUBLIC :: parinfo
    CHARACTER(30) :: name = "" !!!
    CHARACTER(30) :: comp = "" !!!
    REAL(DP) :: value = 0._DP
    REAL(DP), DIMENSION(2) :: limits = [0._DP,0._DP]
    LOGICAL, DIMENSION(2) :: limited = [.False.,.False.]
    LOGICAL :: fixed = .False.
    LOGICAL :: hyper = .True.
    LOGICAL :: asis = .True.
    LOGICAL :: model = .True.
    CHARACTER(30) :: tied = ""
    REAL(DP) :: mean = 0._DP
    REAL(DP) :: sigma = 0._DP
    INTEGER :: ind = -1
  END TYPE parinfo

  TYPE, PUBLIC :: Qabs_str
    !! rho(kg/m3); wave(micron); nu(Hz); Qova(m-2)
    REAL(DP)                            :: rho
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave, nu
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Qova ! a0=1.E-2_DP micron
  END TYPE Qabs_str

  PUBLIC :: make_par, make_Qabs
  PUBLIC :: degradeRes, modifBB, gaussLine, lorentzBand, extCurve, specModel
 
CONTAINS

  !!-------------------------------------------------------
  !!
  !!             Create the parameter structure
  !!
  !!-------------------------------------------------------
  SUBROUTINE make_par(parr, Nbb, Nline, Nband, NAv, Nstar, par)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nbb, Nline, Nband, NAv, Nstar
    REAL(DP), DIMENSION(:), INTENT(IN) :: parr ! par array
    INTEGER i, i0
    TYPE(par_str), INTENT(OUT) :: par

    ALLOCATE(par%massBB(Nbb), par%tempBB(Nbb)) ! Nbb*2
    ALLOCATE(par%Iline(Nline), par%Cline(Nline), par%Wline(Nline)) ! Nline*3
    ALLOCATE(par%Iband(Nband), par%Cband(Nband))
    ALLOCATE(par%WSband(Nband), par%WLband(Nband)) ! Nband*4
    ALLOCATE(par%Av(NAv)) ! NAv*1
    ALLOCATE(par%Fstar(Nstar))
    !! Nstar*1

    !! 1) Continuum
    i0 = 0
    DO i=1,Nbb
      par%massBB(i) = parr(i0+2*i-1) ! i0 + 2*(i-1) + 1
      par%tempBB(i) = parr(i0+2*i)
      
    END DO

    !! 2) Lines
    i0 = i0 + 2*Nbb
    DO i=1,Nline
      par%Iline(i) = parr(i0+3*i-2)
      par%Cline(i) = parr(i0+3*i-1)
      par%Wline(i) = parr(i0+3*i)
      
    END DO

    !! 3) Bands
    i0 = i0 + 3*Nline
    DO i=1,Nband
      par%Iband(i) = parr(i0+4*i-3)
      par%Cband(i) = parr(i0+4*i-2)
      par%WSband(i) = parr(i0+4*i-1)
      par%WLband(i) = parr(i0+4*i)
      
    END DO

    !! 4) Av
    i0 = i0 + 4*Nband
    DO i=1,NAv
      par%Av(i) = parr(i0+i)

    END DO

    !! 5) Star
    i0 = i0 + NAv
    IF (Nstar .GT. 0) THEN
      DO i=1,Nstar
        par%Fstar(i) = parr(i0+1)
        
      END DO
    END IF
  
  END SUBROUTINE make_par

  !!-------------------------------------------------------
  !!
  !!                 Read optical properties
  !!
  !!-------------------------------------------------------
  SUBROUTINE make_Qabs(label, nQabs)
    !! wave [um], nu [Hz], Qova [m^-1], rho [kg/m^3]
    USE utilities, ONLY: DP
    USE constants, ONLY: MKS
    USE arrays, ONLY: closest
    USE grain_optics, ONLY: rho_grain, read_optics
    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:), INTENT(IN)  :: label
    INTEGER                                 :: i, Nspec, Nr, Nw
    INTEGER, DIMENSION(SIZE(label))         :: rind ! radius index
    REAL(DP), PARAMETER                     :: a0 = 1.E-2_DP ! grain radius
    REAL(DP), DIMENSION(:), ALLOCATABLE     :: wave0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: radius
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Q ! (Nspec, Nr, Nw)
    TYPE(Qabs_str), DIMENSION(SIZE(label)), INTENT(OUT) :: nQabs
    
    CALL read_optics(label, WAVE=wave0, RADIUSALL=radius, QABSALL=Q)

    Nspec = SIZE(label)
    Nw = SIZE(wave0)

    DO i=1,Nspec
      
      ALLOCATE(nQabs(i)%Qova(Nw))
      
      Nr = SIZE(radius(i,:))
      nQabs(i)%rho = rho_grain(label(i))
      nQabs(i)%wave = wave0
      nQabs(i)%nu = MKS%clight/MKS%micron / wave0
      rind(i) = closest(radius(i,:), a0)
      nQabs(i)%Qova = Q(i,rind(i),:)/radius(i,rind(i))
      PRINT*, 'Radius of ', label(i), ': ', radius(i,rind(i)), '->', rind(i)

    END DO
    
    !! Free memory space
    DEALLOCATE(wave0, radius, Q)

  END SUBROUTINE make_Qabs

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
          degradeRes = RES%dw_w_SL * wc
        ELSE
          degradeRes = RES%dw_w_LL * wc
        END IF
      CASE('Ns-SL')
        IF (wc.LT.5._DP) THEN
          degradeRes = RES%dw_w_AKARI_Ns * wc
        ELSE
          degradeRes = RES%dw_w_SL * wc
        END IF
      CASE('CAM')
        degradeRes = RES%dw_w_CAM * wc
      CASE('SWS')
        degradeRes = RES%dw_w_SWS * wc
      CASE('SWSfine')
        degradeRes = RES%dw_w_SWSfine * wc
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
  PURE FUNCTION modifBB(temp, Qabs)
    !! temp [K]; blackbody [W/m2/sr/Hz]; 
    !! output array [MJy/sr]
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE
    
    TYPE(Qabs_str), INTENT(IN)          :: Qabs
    REAL(DP), INTENT(IN)                :: temp
    INTEGER                             :: Nw
    REAL(DP), DIMENSION(:), ALLOCATABLE :: modifBB

    Nw = SIZE(Qabs%nu)
    
    ALLOCATE(modifBB(Nw))

    modifBB = 3._DP*pi/4._DP/Qabs%rho * Qabs%Qova &
              * blackbody(Qabs%nu, temp) ! in W/m2/sr/Hz
    !! Unit onversion ([W/m2/sr/Hz] to [MJy/sr])
    modifBB = modifBB /MKS%Jy/1.E6
    
  END FUNCTION modifBB

  !!-----------------------------------------------------
  !! Atomic & molecular unresolved lines (Gauss profile)
  !!-----------------------------------------------------
  PURE FUNCTION gaussLine(dblarr, ref, sig, w2nu)
    !! w2nu=.TRUE. if input is in wavelength
    !! Output array [Hz^-1]
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE statistics, ONLY: mean, median, sigma
    USE distributions, ONLY: dist_gauss
    IMPLICIT NONE

    LOGICAL, INTENT(IN), OPTIONAL      :: w2nu
    REAL(DP), INTENT(IN)               :: ref, sig
    REAL(DP), DIMENSION(:), INTENT(IN) :: dblarr
    REAL(DP)                           :: nuref, nusig
    REAL(DP), DIMENSION(SIZE(dblarr))  :: nuIN
    REAL(DP), DIMENSION(SIZE(dblarr))  :: gaussLine

    IF (PRESENT(w2nu)) THEN
      IF (w2nu) THEN
        nuIN = MKS%clight/MKS%micron / dblarr
        nuref = MKS%clight/MKS%micron / ref
        nusig = MKS%clight/MKS%micron * (1./(ref-sig) - 1./(ref+sig)) / 2._DP

      END IF
    ELSE
      nuIN = dblarr
      nuref = ref
      nusig = sig
      ! print*, dblarr(582), SIZE(dblarr) ! remove PURE for printing tests
      
    END IF
    gaussLine = dist_gauss(nuIN, nuref, nusig)

  END FUNCTION gaussLine

  !!------------------------------------------------------
  !! Resolved aromatic bands (Asymmetric Lorentz profile)
  !!------------------------------------------------------
  PURE FUNCTION lorentzBand(dblarr, ref, sigS, sigL, w2nu)
    !! w2nu=.TRUE. if input is in wavelength
    !! Output array [Hz^-1]
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE distributions, ONLY: dist_splitlorentz
    IMPLICIT NONE

    LOGICAL, INTENT(IN), OPTIONAL      :: w2nu
    REAL(DP), INTENT(IN)               :: ref, sigS, sigL
    REAL(DP), DIMENSION(:), INTENT(IN) :: dblarr
    REAL(DP)                           :: nuref, nusigS, nusigL
    REAL(DP), DIMENSION(SIZE(dblarr))  :: nuIN
    REAL(DP)                           :: lambda, tau
    REAL(DP), DIMENSION(SIZE(dblarr))  :: lorentzBand

    IF (PRESENT(w2nu)) THEN
      IF (w2nu) THEN
        nuIN = MKS%clight/MKS%micron / dblarr
        nuref = MKS%clight/MKS%micron / ref
        nusigS = MKS%clight/MKS%micron * (1./(ref-sigS) - 1./ref)
        nusigL = MKS%clight/MKS%micron * (1./ref - 1./(ref+sigL))
      END IF
    ELSE
      nuIN = dblarr
      nuref = ref
      nusigS = sigS
      nusigL = sigL
      
    END IF
    !! The SwING library adopted the arbitary parameters such that 
    !! dnu_short corresponds to longer wavelength side of the curve,
    !! Note that neither dnu_S nor dnu_L is std_dev of a lorentzian,
    !! 'cause it is undefined (infinite).
    !! This is analogue of (split) normal distribution.
    !! In addition we denote lambda=2*dnu_long (analogue of FWHM)
    !! and tau=dnu_short/dnu_long (shape param).
    lambda = 2*nusigL ! nusigL = lambda/2 (lambda -> width param)
    tau = nusigS/nusigL ! nusigS = tau * nusigL (tau -> shape param)
    !! norm = lambda / (1+tau) / pi
    lorentzBand = dist_splitlorentz(nuIN, nuref, lambda, tau)

  END FUNCTION lorentzBand

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
  FUNCTION specModel(waveIN, parr, nQabs, Nbb, Nline, Nband, NAv, Nstar, &
                     tau_tab, Fnu_cont_tab, Fnu_line_tab, Fnu_band_tab, &
                     Pabs, Fnu_cont, Fnu_line, Fnu_band, Fnu_star)
    !! Nstar = 0 or 1
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: REALLOCATE
    USE interpolation, ONLY: interp_lin_sorted
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN)       :: waveIN
    REAL(DP), DIMENSION(:), INTENT(IN)       :: parr
    TYPE(Qabs_str), DIMENSION(:), INTENT(IN) :: nQabs
    INTEGER, INTENT(IN)                      :: Nbb, Nline, Nband, NAv, Nstar
    TYPE(par_str)                                                :: par
    INTEGER                                                      :: i, Nw
    REAL(DP)                                                     :: Tstar
    REAL(DP), DIMENSION(SIZE(waveIN))                            :: nu
    REAL(DP), DIMENSION(:,:), ALLOCATABLE                        :: tau_tab0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE                        :: Fnu_cont_tab0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE                        :: Fnu_line_tab0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE                        :: Fnu_band_tab0
    REAL(DP), DIMENSION(:), ALLOCATABLE                          :: Pabs0
    REAL(DP), DIMENSION(:), ALLOCATABLE                          :: Fnu_cont0
    REAL(DP), DIMENSION(:), ALLOCATABLE                          :: Fnu_line0
    REAL(DP), DIMENSION(:), ALLOCATABLE                          :: Fnu_band0
    REAL(DP), DIMENSION(:), ALLOCATABLE                          :: Fnu_star0
    REAL(DP), DIMENSION(NAv,SIZE(waveIN)), INTENT(OUT), OPTIONAL :: tau_tab
    REAL(DP), DIMENSION(Nbb,SIZE(waveIN)), INTENT(OUT), OPTIONAL :: Fnu_cont_tab
    REAL(DP), DIMENSION(Nline,SIZE(waveIN)), INTENT(OUT), OPTIONAL :: Fnu_line_tab
    REAL(DP), DIMENSION(Nband,SIZE(waveIN)), INTENT(OUT), OPTIONAL :: Fnu_band_tab
    REAL(DP), DIMENSION(SIZE(waveIN)), INTENT(OUT), OPTIONAL   :: Pabs
    REAL(DP), DIMENSION(SIZE(waveIN)), INTENT(OUT), OPTIONAL   :: Fnu_cont
    REAL(DP), DIMENSION(SIZE(waveIN)), INTENT(OUT), OPTIONAL   :: Fnu_line
    REAL(DP), DIMENSION(SIZE(waveIN)), INTENT(OUT), OPTIONAL   :: Fnu_band
    REAL(DP), DIMENSION(SIZE(waveIN)), INTENT(OUT), OPTIONAL   :: Fnu_star
    REAL(DP), DIMENSION(:), ALLOCATABLE                          :: specModel

    !! Preliminaries
    !!---------------
    Nw = SIZE(waveIN)
    nu = MKS%clight/MKS%micron / waveIN

    CALL make_par(parr, Nbb, Nline, Nband, NAv, Nstar, par)

    ALLOCATE(Pabs0(Nw), Fnu_cont0(Nw), Fnu_line0(Nw), Fnu_band0(Nw), Fnu_star0(Nw))
    ALLOCATE(specModel(Nw))

    !! Screen extinction
    !!-------------------
    IF (NAv.GT.0) THEN
      ALLOCATE(tau_tab0(NAv, Nw))
      DO i=1,NAv
        tau_tab0(i,:) = par%Av(i)/1.086_DP * extCurve(waveIN)

      END DO
    ELSE
      CALL REALLOCATE(tau_tab0, 1, Nw)

    END IF
    Pabs0 = EXP(-SUM(tau_tab0(:,:), DIM=1))
    !! output option
    IF (PRESENT(tau_tab)) tau_tab = tau_tab0
    IF (PRESENT(Pabs)) Pabs = Pabs0

    !! Continuum
    !!-----------
    IF (Nbb.GT.0) THEN
      ALLOCATE(Fnu_cont_tab0(Nbb, Nw))
      DO i=1,Nbb
        !! massBB in unit of [Msun/pc2]
        Fnu_cont_tab0(i,:) = par%massBB(i) * MKS%Msun/MKS%pc**2 * &
                             interp_lin_sorted(modifBB(par%tempBB(i), nQabs(i)), &
                             nQabs(i)%nu, nu, xlog=.True., ylog=.True.)

      END DO
    ELSE
      CALL REALLOCATE(Fnu_cont_tab0, 1, Nw)

    END IF
    Fnu_cont0 = SUM(Fnu_cont_tab0(:,:), DIM=1)
    !! output option
    IF (PRESENT(Fnu_cont_tab)) Fnu_cont_tab = Fnu_cont_tab0
    IF (PRESENT(Fnu_cont)) Fnu_cont = Fnu_cont0

    !! Lines
    !!-------
    IF (Nline.GT.0) THEN
      ALLOCATE(Fnu_line_tab0(Nline, Nw))
      DO i=1,Nline
        Fnu_line_tab0(i,:) = par%Iline(i) /MKS%Jy/1.E6 * &
                             gaussLine(waveIN, par%Cline(i), par%Wline(i), .TRUE.)

      END DO
    ELSE
      CALL REALLOCATE(Fnu_line_tab0, 1, Nw)

    END IF
    Fnu_line0 = SUM(Fnu_line_tab0(:,:), DIM=1)
    !! output option
    IF (PRESENT(Fnu_line_tab)) Fnu_line_tab = Fnu_line_tab0
    IF (PRESENT(Fnu_line)) Fnu_line = Fnu_line0

    !! Bands
    !!-------
    IF (Nband.GT.0) THEN
      ALLOCATE(Fnu_band_tab0(Nband, Nw))
      DO i=1,Nband
        Fnu_band_tab0(i,:) = par%Iband(i) /MKS%Jy/1.E6 * &
                             lorentzBand(waveIN, par%Cband(i), &
                             par%WSband(i), par%WLband(i), .TRUE.)

      END DO
    ELSE
      CALL REALLOCATE(Fnu_band_tab0, 1, Nw)

    END IF
    Fnu_band0 = SUM(Fnu_band_tab0(:,:), DIM=1)
    !! output option
    IF (PRESENT(Fnu_band_tab)) Fnu_band_tab = Fnu_band_tab0
    IF (PRESENT(Fnu_band)) Fnu_band = Fnu_band0

    !! Stellar Continuum
    !!-------------------
    IF (Nstar.GT.0) THEN
      DO i=1,Nstar
        Tstar = 5.E4_DP ! high enough to stay in Rayleigh-Jeans limit [K]
        !! Fstar in unit of [Lsun/pc2] & Fnu_star [MJy/sr] with notations as follows: 
        !! Fnu [W/m2/Hz] = pi * Bnu [W/m2/sr/Hz] (isotropically emitting surface)
        !! F [W/m2] = integ(Fnudnu) = pi * integ(Bnudnu) = pi * B [W/m2] = sig * T^4
        !! cf. energy flux (luminosity) [W]; surface brightness [W/m2] (or [mag]); 
        !! flux density (irradiance) [W/m2]; radiant flux density (intensity) [W/m2/sr]
        Fnu_star0 = par%Fstar(i) * MKS%Lsun/4._DP/pi/MKS%pc**2 * &
                    pi * blackbody(nu, Tstar) /MKS%stefan/Tstar**4 /MKS%Jy/1.E6

      END DO
    END IF
    !! output option
    IF (PRESENT(Fnu_star)) Fnu_star = Fnu_star0

    !! Total
    !!-------
    specModel(:) = (Fnu_cont0+Fnu_band0+Fnu_star0)*Pabs0 + Fnu_line0

    DEALLOCATE(tau_tab0, Fnu_cont_tab0, Fnu_line_tab0, Fnu_band_tab0)
    DEALLOCATE(Pabs0, Fnu_cont0, Fnu_line0, Fnu_band0, Fnu_star0)

  END FUNCTION specModel

END MODULE auxil
