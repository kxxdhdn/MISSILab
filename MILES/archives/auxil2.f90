MODULE auxil

  USE utilities, ONLY: DP
  USE constants, ONLY: pi, MKS
  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: par_STR

    REAL(DP), DIMENSION(:), ALLOCATABLE :: massBB
    REAL(DP), DIMENSION(:), ALLOCATABLE :: tempBB
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Iline ! Intensity
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Cline ! Center
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Wline ! Width
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Iband
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Cband
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WSband
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WLband
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Av
    REAL(DP)                            :: massStar = 0._DP

  END TYPE par_STR
  TYPE(par_STR), SAVE, PUBLIC :: par

  TYPE, PUBLIC :: parinfo_TYPE

    CHARACTER(30) :: name = ""
    CHARACTER(30) :: comp = ""
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

  END TYPE parinfo_TYPE

  TYPE, PUBLIC :: Qabs_STR
    !! rho(kg/m3); wave(micron); nu(Hz); Qova(m-2)
    REAL(DP)                            :: rho
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave, nu
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Qova ! a0=1.E-2_DP micron

  END TYPE Qabs_STR

  PUBLIC :: par_STRUCT
  PUBLIC :: optics_INIT, spec_MODEL, chi2_INIT
  PUBLIC :: degradeRes, modifBB, gaussLine, lorentzBand
 
CONTAINS

  !!-------------------------------------------------------
  !!
  !!             Create the parameter structure
  !!
  !!-------------------------------------------------------
  SUBROUTINE par_STRUCT(parm, Nbb, Nline, Nband, NAv, Nstar, par)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nbb, Nline, Nband, NAv, Nstar
    REAL(DP), DIMENSION(:), INTENT(IN) :: parm
    INTEGER i, i0
    TYPE(par_STR), INTENT(OUT) :: par

    ALLOCATE(par%massBB(Nbb), par%tempBB(Nbb))
    ALLOCATE(par%Iline(Nline), par%Cline(Nline), par%Wline(Nline))
    ALLOCATE(par%Iband(Nband), par%Cband(Nband))
    ALLOCATE(par%WSband(Nband), par%WLband(Nband))
    ALLOCATE(par%Av(NAv))

    !! 1) Continuum
    i0 = 0
    DO i=1,Nbb
      par%massBB(i) = parm(i0+2*i-1) ! i0 + 2*(i-1) + 1
      par%tempBB(i) = parm(i0+2*i)

    END DO

    !! 2) Lines
    i0 = i0 + 2*Nbb
    DO i=1,Nline
      par%Iline(i) = parm(i0+3*i-2)
      par%Cline(i) = parm(i0+3*i-1)
      par%Wline(i) = parm(i0+3*i)

    END DO

    !! 3) Bands
    i0 = i0 + 3*Nline
    DO i=1,Nband
      par%Iband(i) = parm(i0+4*i-3)
      par%Cband(i) = parm(i0+4*i-2)
      par%WSband(i) = parm(i0+4*i-1)
      par%WLband(i) = parm(i0+4*i)

    END DO

    !! 4) Av
    i0 = i0 + 4*Nband
    DO i=1,NAv
      par%Av(i) = parm(i0+i)

    END DO

    !! 5) Star
    i0 = i0 + NAv
    IF (Nstar .GT. 0) par%massStar = parm(i0+1)
  
  END SUBROUTINE par_STRUCT

  !!-------------------------------------------------------
  !!
  !!                 Read optical properties
  !!
  !!-------------------------------------------------------
  SUBROUTINE optics_INIT(label, nQ_str)

    USE utilities, ONLY: DP
    USE constants, ONLY: MKS
    USE arrays, ONLY: closest
    USE grain_optics, ONLY: rho_grain, read_optics
    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:), INTENT(IN)  :: label
    INTEGER                                 :: i, Nspec, Nw, Nr
    INTEGER, DIMENSION(SIZE(label))        :: ind ! radius index
    REAL(DP), PARAMETER                     :: a0=1.E-2_DP ! grain radius
    REAL(DP), DIMENSION(:), ALLOCATABLE     :: waveIN
    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: radius
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Qabs
    TYPE(Qabs_STR), DIMENSION(SIZE(label)), INTENT(OUT) :: nQ_str
    
    CALL read_optics(label, WAVE=waveIN, &
                     RADIUSALL=radius, QABSALL=Qabs)

    Nspec = SIZE(label)
    Nw = SIZE(waveIN)

    DO i=1,Nspec
      
      ALLOCATE(nQ_str(i)%Qova(Nw))
      
      Nr = SIZE(radius(i,:))
      nQ_str(i)%rho = rho_grain(label(i))
      nQ_str(i)%wave = waveIN
      nQ_str(i)%nu = MKS%clight/MKS%micron / waveIN
      ind(i) = closest(radius(i,:), a0)
      nQ_str(i)%Qova = Qabs(i,ind(i),:)/radius(i,ind(i))
      PRINT*, 'Radius of ', label(i), ': ', radius(i,ind(i)), '->', ind(i)

    END DO
    
    !! Free memory space
    DEALLOCATE(waveIN, radius, Qabs)

  END SUBROUTINE optics_INIT

  !!-------------------------------------------------------
  !!
  !! Automatize the degradation of the spectral resolution
  !!
  !!-------------------------------------------------------
  FUNCTION degradeRes(wc, dw, instr)

    USE factable, ONLY: RES
    USE utilities, ONLY: DP
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)         :: instr
    REAL(DP), INTENT(IN)             :: wc, dw
    REAL(DP)                         :: degradeRes
    
    !! Line width fixed by the instrumental resolution
    SELECT CASE(instr)
      CASE('SH')
        degradeRes = RES%SH * wc
      CASE('LH')
        degradeRes = RES%LH * wc
      CASE('SL')
        degradeRes = RES%SL * wc
      CASE('LL')
        degradeRes = RES%LL * wc
      CASE('SL-SH')
        IF (wc.LT.10._DP) THEN
          degradeRes = RES%SL * wc
        ELSE
          degradeRes = RES%SH * wc
        END IF
      CASE('SL-LL')
        IF (wc.LT.12._DP) THEN
          degradeRes = RES%SL * wc
        ELSE
          degradeRes = RES%LL * wc
        END IF
      CASE('Ns-SL')
        IF (wc.LT.5._DP) THEN
          degradeRes = RES%AKARI_Ns * wc
        ELSE
          degradeRes = RES%SL * wc
        END IF
      CASE('CAM')
        degradeRes = RES%CAM * wc
      CASE('SWS')
        degradeRes = RES%SWS * wc
      CASE('SWSfine')
        degradeRes = RES%SWSfine * wc
      CASE DEFAULT
        degradeRes = dw

    END SELECT

  END FUNCTION degradeRes

  !!-------------------------------------------------------
  !!
  !!    Analytical functions of the individual features
  !!
  !!-------------------------------------------------------

  !! Dust contimuum (N BB)
  !!-----------------------
  PURE FUNCTION modifBB(temp, Q_str)
    !! temp(K); blackbody(W/m2/sr/Hz=Jy/sr); 
    !! Fnu(MJy/sr) [in Msun/kpc^2]
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE
    
    TYPE(Qabs_STR), INTENT(IN)          :: Q_str
    REAL(DP), INTENT(IN)                :: temp
    INTEGER                             :: Nw
    REAL(DP), DIMENSION(:), ALLOCATABLE :: modifBB

    Nw = SIZE(Q_str%nu)
    
    ALLOCATE(modifBB(Nw))

    modifBB = 3._DP*pi/4._DP/Q_str%rho * Q_str%Qova &
              * blackbody(Q_str%nu, temp)
    !! unit conversion (MJy/sr) [in Msun/kpc^2]
    modifBB = modifBB * MKS%Msun/MKS%kpc**2 / 1.E6_DP

  END FUNCTION modifBB

  !! Atomic & molecular unresolved lines (Gauss profile)
  !!-----------------------------------------------------
  PURE FUNCTION gaussLine(nuIN, nucen, nusig)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    USE statistics, ONLY: mean, median, sigma
    USE distributions, ONLY: dist_gauss
    IMPLICIT NONE
    
    REAL(DP), INTENT(IN)               :: nucen, nusig
    REAL(DP), DIMENSION(:), INTENT(IN) :: nuIN
    REAL(DP), DIMENSION(SIZE(nuIN))    :: gaussLine

    gaussLine = dist_gauss(nuIN, nucen, nusig)

  END FUNCTION gaussLine

  PURE FUNCTION gaussLine_w(waveIN, wref, sigma)
    !! Input wave (obsolete)
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE statistics, ONLY: mean, median, sigma
    USE distributions, ONLY: dist_gauss
    IMPLICIT NONE
    
    REAL(DP), INTENT(IN)                :: wref, sigma
    REAL(DP), DIMENSION(:), INTENT(IN)  :: waveIN
    REAL(DP)                            :: nucen, nusig
    REAL(DP), DIMENSION(SIZE(waveIN))   :: nuIN
    REAL(DP), DIMENSION(:), ALLOCATABLE :: gaussLine_w

    nuIN = MKS%clight/MKS%micron / waveIN

    nucen = MKS%clight/MKS%micron / wref
    nusig = MKS%clight/MKS%micron * (1./(wref-sigma) - 1./(wref+sigma)) / 2.

    ALLOCATE(gaussLine_w(SIZE(waveIN)))

    gaussLine_w = dist_gauss(nuIN, nucen, nusig)

  END FUNCTION gaussLine_w

  !! Resolved aromatic bands (Asymmetric Lorentz profile)
  !!------------------------------------------------------
  PURE FUNCTION lorentzBand(nuIN, nuref, nusigS, nusigL)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    USE distributions, ONLY: dist_splitlorentz
    IMPLICIT NONE
    
    REAL(DP), INTENT(IN)               :: nuref, nusigS, nusigL
    REAL(DP), DIMENSION(:), INTENT(IN) :: nuIN
    REAL(DP)                           :: lambda, tau
    REAL(DP), DIMENSION(SIZE(nuIN))    :: lorentzBand

    lambda = 2*nusigL ! nusigL = lambda/2 (lambda -> width param)
    tau = nusigS/nusigL ! nusigS = tau * nusigL (tau -> shape param)
    !! norm = lambda / (1+tau) / pi
    lorentzBand = dist_splitlorentz(nuIN, nuref, lambda, tau)

  END FUNCTION lorentzBand

  PURE FUNCTION lorentzBand_w(waveIN, wref, sigmaS, sigmaL)
    !! Input wave (obsolete)
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE distributions, ONLY: dist_splitlorentz
    IMPLICIT NONE
    
    REAL(DP), INTENT(IN)                :: wref, sigmaS, sigmaL
    REAL(DP), DIMENSION(:), INTENT(IN)  :: waveIN
    REAL(DP)                            :: lambda, tau
    REAL(DP)                            :: nucen, nusigS, nusigL
    REAL(DP), DIMENSION(SIZE(waveIN))   :: nuIN
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lorentzBand_w

    nuIN = MKS%clight/MKS%micron / waveIN

    nucen = MKS%clight/MKS%micron / wref
    nusigS = MKS%clight/MKS%micron * (1./(wref-sigmaS) - 1./wref)
    nusigL = MKS%clight/MKS%micron * (1./wref - 1./(wref+sigmaL))

    ALLOCATE(lorentzBand_w(SIZE(waveIN)))

    lambda = 2*nusigL
    tau = nusigS/nusigL
    !! norm = lambda / (1+tau) / pi
    lorentzBand_w = dist_splitlorentz(nuIN, nucen, lambda, tau)

  END FUNCTION lorentzBand_w

  !! Total model
  !!-------------
  SUBROUTINE spec_MODEL(label, parm, NAv, Nbb, Nline, Nband, Nstar, &
                        Pabs, Fnu_cont, Fnu_line, Fnu_band, Fnu_star, &
                        FnuOUT, waveIN)

    USE utilities, ONLY: DP
    USE constants, ONLY: MKS
    USE arrays, ONLY: REALLOCATE
    USE statistical_physics, ONLY: blackbody
    USE inout, ONLY: write_hdf5, h5ext
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(IN), OPTIONAL :: waveIN
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(IN)           :: parm
    CHARACTER(*), DIMENSION(:), ALLOCATABLE, INTENT(IN)       :: label
    INTEGER, INTENT(IN) :: NAv, Nbb, Nline, Nband, Nstar
    INTEGER                                   :: i, Nw
    TYPE(par_STR)                             :: par
    TYPE(Qabs_STR), DIMENSION(:), ALLOCATABLE :: nQ_str
    REAL(DP), DIMENSION(:), ALLOCATABLE       :: nu, wave
    REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: tau_tab
    REAL(DP), DIMENSION(:,:), ALLOCATABLE              :: Fnu_cont_tab
    REAL(DP), DIMENSION(:,:), ALLOCATABLE              :: Fnu_line_tab
    REAL(DP), DIMENSION(:,:), ALLOCATABLE              :: Fnu_band_tab
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)   :: Pabs
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)   :: Fnu_cont
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)   :: Fnu_line
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)   :: Fnu_band
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)   :: Fnu_star
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)   :: FnuOUT
    
    !! Preliminaries
    !!---------------
    !! Total model parameters
    CALL par_STRUCT(parm, Nbb, Nline, Nband, NAv, Nstar, par)
    !! Continuum (modifBB) parameters
    ALLOCATE(nQ_str(Nbb))
    CALL optics_INIT(label, nQ_str)

    !! I/O parameters
    IF (PRESENT(waveIN)) THEN
      Nw = SIZE(waveIN)
      ALLOCATE(wave(Nw))
      wave =waveIN ! interpolate grid
    ELSE
      Nw = SIZE(nQ_str(1)%nu)
      ALLOCATE(wave(Nw))
      wave = nQ_str(1)%wave ! init grid (often more fin)
      PRINT*, '>> Default (1) wavelength grid. <<'

    END IF
    ALLOCATE(nu(Nw))
    nu = MKS%clight/MKS%micron / wave
    ALLOCATE(Pabs(Nw), Fnu_cont(Nw), Fnu_line(Nw), Fnu_band(Nw))
    ALLOCATE(Fnu_star(Nw), FnuOUT(Nw))

    ! Screen extinction
    !-------------------
    IF (NAv.GT.0) THEN
      ALLOCATE(tau_tab(NAv, Nw))
      DO i=1,NAv
        tau_tab(i,:) = par%Av(i)/1.086_DP*1._DP ! tau0 ?

      END DO
    ELSE
      CALL REALLOCATE(tau_tab, 1, Nw)

    END IF
    Pabs = EXP(-SUM(tau_tab(:,:), DIM=1))

    !! Continuum
    !!-----------
    IF (Nbb.GT.0) THEN
      ALLOCATE(Fnu_cont_tab(Nbb, Nw))
      DO i=1,Nbb
        Fnu_cont_tab(i,:) = par%massBB(i) &
                            * modifBB(par%tempBB(i), nQ_str(i))

      END DO
    ELSE
      CALL REALLOCATE(Fnu_cont_tab, 1, Nw)

    END IF
    Fnu_cont = SUM(Fnu_cont_tab(:,:), DIM=1)

    !! Lines
    !!-------
    IF (Nline.GT.0) THEN
      ALLOCATE(Fnu_line_tab(Nline, Nw))
      DO i=1,Nline
        Fnu_line_tab(i,:) = gaussLine_w(wave, par%Cline(i), par%Wline(i))

      END DO
    ELSE
      CALL REALLOCATE(Fnu_line_tab, 1, Nw)

    END IF
    Fnu_line = SUM(Fnu_line_tab(:,:), DIM=1)
    ! print*, par%Wline

    !! Bands
    !!-------
    IF (Nband.GT.0) THEN
      ALLOCATE(Fnu_band_tab(Nband, Nw))
      DO i=1,Nband
        Fnu_band_tab(i,:) = lorentzBand_w(wave, par%Cband(i), &
                            par%WSband(i), par%WLband(i))

      END DO
    ELSE
      CALL REALLOCATE(Fnu_band_tab, 1, Nw)

    END IF
    Fnu_band = SUM(Fnu_band_tab(:,:), DIM=1)

    !! Stellar Continuum
    !!-------------------
    IF (Nstar.GT.0) THEN
      Fnu_star = par%massStar * blackbody(nu, 5.E4_DP)

    END IF

    !! Final fit & write
    !!-------------------
    FnuOUT = (Fnu_cont+Fnu_band+Fnu_star)*Pabs + Fnu_line

    CALL write_hdf5(wave, NAME="Wavelength (micron)", &
                    COMPRESS=.False., APPEND=.False., &
                    FILE="output/totalout"//h5ext)
    CALL write_hdf5(Fnu_cont*Pabs, NAME="FnuBB (MJy/sr)", &
                    COMPRESS=.False., APPEND=.True., &
                    FILE="output/totalout"//h5ext)
    CALL write_hdf5(Fnu_line+(Fnu_cont+Fnu_star)*Pabs, &
                    NAME="FnuLINE (MJy/sr)", &
                    COMPRESS=.False., APPEND=.True., &
                    FILE="output/totalout"//h5ext)
    CALL write_hdf5((Fnu_band+Fnu_cont+Fnu_star)*Pabs, &
                    NAME="FnuBAND (MJy/sr)", &
                    COMPRESS=.False., APPEND=.True., &
                    FILE="output/totalout"//h5ext)
    CALL write_hdf5(Fnu_star*Pabs, NAME="FnuSTAR (MJy/sr)", &
                    COMPRESS=.False., APPEND=.True., &
                    FILE="output/totalout"//h5ext)
    CALL write_hdf5(FnuOUT, NAME="Fnu (MJy/sr)", &
                    COMPRESS=.False., APPEND=.True., &
                    FILE="output/totalout"//h5ext)

    DEALLOCATE(nQ_str, tau_tab, Fnu_cont_tab, Fnu_line_tab, Fnu_band_tab)

  END SUBROUTINE spec_MODEL

  !! Initialization of parameters for chi2 method
  !!----------------------------------------------
  SUBROUTINE chi2_INIT(label, parm, NAv, Nbb, Nline, Nband, Nstar)

    USE factable, ONLY: LIN, BIN
    USE utilities, ONLY: DP
    USE constants, ONLY: MKS
    IMPLICIT NONE

    INTEGER :: Nparm
    INTEGER, INTENT(OUT) :: Nbb, Nline, Nband, NAv, Nstar
    CHARACTER(*), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: label
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)      :: parm

    Nbb = 3
    label = ['ACH2_Z96             ', &
             'BE_Z96               ', &
             'Sil_D03              ']
    ! Nbb = SIZE(label)
 
    NAv = 1
    Nline = 3
    Nband = 2
    Nstar = 1
    Nparm = 2*Nbb + 3*Nline + 4*Nband + NAv + Nstar

    ALLOCATE(parm(Nparm))
    !! parm = [massBB, tempBB, &
    !!         Iline, Cline, Wline, &
    !!         Iband, Cband, WSband, WLband, &
    !!         Av, &
    !!         massStar]
    parm = [1.E23_DP, 100._DP, .1_DP, 100._DP, .1_DP, 100._DP, &
            2._DP, MKS%clight/MKS%micron / LIN(5)%wave, &
            degradeRes(LIN(5)%wave, .01_DP, 'SL-LL'), &
            3._DP, MKS%clight/MKS%micron / LIN(10)%wave, &
            degradeRes(LIN(10)%wave, .01_DP, 'SL-LL'), &
            1._DP, MKS%clight/MKS%micron / LIN(25)%wave, &
            degradeRes(LIN(25)%wave, .01_DP, 'SL-LL'), &
            10._DP, MKS%clight/MKS%micron / BIN(10)%wave, &
            BIN(10)%sigmaS, BIN(10)%sigmaL, &
            5._DP, MKS%clight/MKS%micron / BIN(15)%wave, &
            BIN(15)%sigmaS, BIN(15)%sigmaL, &
            .1_DP, &
            5.E-16_DP]

  END SUBROUTINE chi2_INIT

END MODULE auxil
