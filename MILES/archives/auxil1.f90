MODULE auxil

  USE utilities, ONLY: DP
  USE constants, ONLY: pi, MKS
  IMPLICIT NONE
  !! A sr * sr2arcsec2 = B arcsec^2
  REAL(DP), PARAMETER, PUBLIC :: sr2arcsec2 = (180._DP/pi)**2
  !! A deg * deg2arcsec = B arcsec
  REAL(DP), PARAMETER, PUBLIC :: deg2arcsec = 3600._DP
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
    REAL(DP)                            :: rho
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Qova ! a0=1.E-2_DP micron

  END TYPE Qabs_STR

  PUBLIC :: par_STRUCT, optics_INIT, degradeRes
  PUBLIC :: modifBB, gaussLine, lorentzBand
 
CONTAINS

  !!--------------------------------
  !! Create the parameter structure
  !!--------------------------------
  SUBROUTINE par_STRUCT(parm, Nbb, Nline, Nband, NAv, Nstar, par)

    USE utilities, ONLY: DP

    INTEGER i, i0
    INTEGER, INTENT(IN) :: Nbb, Nline, Nband, NAv, Nstar
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: parm
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

  !!-------------------------
  !! Read optical properties
  !!-------------------------
  SUBROUTINE optics_INIT(labels, wave, nu, nQ_str)
    !! rho(kg/m3); wave(micron); nu(Hz); Qova(m-2)
    USE utilities, ONLY: DP
    USE constants, ONLY: MKS
    USE arrays, ONLY: closest
    USE grain_optics, ONLY: rho_grain, read_optics
    IMPLICIT NONE

    INTEGER i, Nspec, Nw, Nr
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind ! radius index
    CHARACTER(30), DIMENSION(:), ALLOCATABLE, INTENT(IN)       :: labels
    REAL(DP), PARAMETER :: a0=1.E-2_DP ! grain radius
    REAL(DP), DIMENSION(:,:), ALLOCATABLE                      :: radius
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: Qabs
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: wave, nu
    TYPE(Qabs_STR), DIMENSION(:), ALLOCATABLE, INTENT(OUT)     :: nQ_str
    
    CALL read_optics(labels, WAVE=wave, RADIUSALL=radius, QABSALL=Qabs)

    Nspec = SIZE(labels)
    Nw = SIZE(wave)

    ALLOCATE(nu(Nw))
    ALLOCATE(nQ_str(Nspec))
    ALLOCATE(ind(Nspec))
    
    nu = MKS%clight / MKS%micron / wave
    DO i=1,Nspec
      
      ALLOCATE(nQ_str(i)%Qova(Nw))
      
      Nr = SIZE(radius(i,:))
      ind(i) = closest(radius(i,:), a0)
      nQ_str(i)%Qova = Qabs(i,ind(i),:)/radius(i,ind(i))
      PRINT*, 'Radius of ', labels(i), ': ', radius(i,ind(i)), '->', ind(i)
      nQ_str(i)%rho = rho_grain(labels(i))

    END DO
    
    !! Free memory space
    DEALLOCATE(ind, radius, Qabs)

  END SUBROUTINE optics_INIT

  !!-------------------------------------------------------
  !! Automatize the degradation of the spectral resolution
  !!-------------------------------------------------------
  FUNCTION degradeRes(wc, dw, instr)

    USE factable, ONLY: RES
    USE utilities, ONLY: DP
    IMPLICIT NONE

    CHARACTER(30), INTENT(IN)        :: instr
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

  !!-------------------------------------------------
  !! Analytical functions of the individual features
  !!-------------------------------------------------

  !! Dust contimuum (N BB)
  !!-----------------------
  FUNCTION modifBB(nu, temp, Q_str)
    !! temp(K); rho(kg/m3); nu(Hz); Qova(m-2);
    !! blackbody(W/m2/sr/Hz=Jy/sr); Fnu(MJy/sr) [in Msun/kpc^2]
    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE
    
    INTEGER Nw
    TYPE(Qabs_STR), INTENT(IN)                      :: Q_str
    REAL(DP), INTENT(IN)                            :: temp
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: nu
    REAL(DP), DIMENSION(:), ALLOCATABLE             :: modifBB

    Nw = SIZE(nu)
    
    ALLOCATE(modifBB(Nw))

    modifBB = 3._DP*pi/4._DP/Q_str%rho * Q_str%Qova &
              * blackbody(nu, temp)
    !! unit conversion (MJy/sr) [in Msun/kpc^2]
    modifBB = modifBB * MKS%Msun/MKS%kpc**2 / 1.E6_DP

  END FUNCTION modifBB

  !! Atomic & molecular unresolved lines (Gauss profile)
  !!-----------------------------------------------------
  FUNCTION gaussLine(nu, wref, sigma)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    USE statistics, ONLY: mean, median, sigma
    USE distributions, ONLY: dist_gauss
    IMPLICIT NONE
    
    INTEGER Nw
    ! INTEGER, INTENT(IN)  :: Nspec
    REAL(DP), INTENT(IN)                            :: wref, sigma
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: nu
    REAL(DP), DIMENSION(:), ALLOCATABLE             :: gaussLine
    
    Nw = SIZE(nu)

    ALLOCATE(gaussLine(Nw))

    gaussLine = dist_gauss(nu, wref, sigma)

  END FUNCTION gaussLine

  !! Resolved aromatic bands (Asymmetric Lorentz profile)
  !!------------------------------------------------------
  FUNCTION lorentzBand(nu, wref, sigmaS, sigmaL)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi
    USE distributions, ONLY: dist_splitlorentz
    IMPLICIT NONE
    
    INTEGER Nw
    ! INTEGER, INTENT(IN)  :: Nspec
    REAL(DP) :: lambda, tau
    REAL(DP), INTENT(IN)                            :: wref, sigmaS, sigmaL
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: nu
    REAL(DP), DIMENSION(:), ALLOCATABLE             :: lorentzBand
    
    Nw = SIZE(nu)

    ALLOCATE(lorentzBand(Nw))

    lambda = 2*sigmaL
    tau = sigmaS/sigmaL
    !! norm = lambda / (1+tau) / pi
    lorentzBand = dist_splitlorentz(nu, wref, lambda, tau)

  END FUNCTION lorentzBand

END MODULE auxil
