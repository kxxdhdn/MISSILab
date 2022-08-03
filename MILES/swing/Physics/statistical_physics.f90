!******************************************************************************
!*
!*                             STATISTICAL PHYSICS
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 08/2007 
  !    - 01/2015: add Maxwell-Boltzmann distribution.
  ! 
  ! 3) DESCRIPTION: Various functions of statistical physics.
  !==========================================================================

MODULE statistical_physics

  USE utilities, ONLY: DP
  USE constants, ONLY:
  USE special_functions, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: blackbody, maxwell_boltzmann

  INTERFACE blackbody
    MODULE PROCEDURE blackbody_s, blackbody_v
  END INTERFACE blackbody

  INTERFACE maxwell_boltzmann
    MODULE PROCEDURE maxwell_boltzmann_s, maxwell_boltzmann_v
  END INTERFACE maxwell_boltzmann

  ! Peak of the black body function
  REAL(DP), PARAMETER, PUBLIC :: wT_nu = 5.0996E-3_DP ! [m/K] wave*T for Bnu_max
  REAL(DP), PARAMETER, PUBLIC :: wT_w = 3.6698E-3_DP ! [m/K] wave*T for Bw_max


CONTAINS


  !==========================================================================
  ! Bnu[N,M] = BLACKBODY(nu[N],Temperature[M])
  !
  !   Evaluates the planck function at wavelengths nu [Hz], in W/m2/sr/Hz.
  !==========================================================================

  PURE FUNCTION blackbody_v (nu,temp)

    USE utilities, ONLY: DP
    USE constants, ONLY: MKS
    USE special_functions, ONLY: expm1
    IMPLICIT NONE 

    REAL(DP), INTENT(IN), DIMENSION(:) :: nu ! Frequency grid in Hz
    REAL(DP), INTENT(IN), DIMENSION(:) :: temp ! Temperature grid in K
    REAL(DP), DIMENSION(SIZE(nu),SIZE(temp)) :: blackbody_v ! Bnu in W/m2/sr/Hz

    REAL(DP) :: twicehovc2, hovk ! Intermediate factors
    REAL(DP), DIMENSION(SIZE(temp)) :: hbeta
    INTEGER :: i, Ntemp 

    !-----------------------------------------------------------------------

    Ntemp = SIZE(temp(:))
    hovk = MKS%hplanck / MKS%kboltz
    twicehovc2 = 2*MKS%hplanck / MKS%clight**2
    hbeta(:) = hovk / temp(:)
 
    FORALL (i=1:Ntemp) &
      blackbody_v(:,i) = twicehovc2*nu**3 / EXPM1(hbeta(i)*nu)

    !-----------------------------------------------------------------------

  END FUNCTION blackbody_v

  !==========================================================================

  PURE FUNCTION blackbody_s (nu,temp)

    USE utilities, ONLY: DP
    IMPLICIT NONE 

    REAL(DP), INTENT(IN), DIMENSION(:) :: nu ! Frequency grid in Hz
    REAL(DP), INTENT(IN) :: temp ! Temperature in K
    REAL(DP), DIMENSION(SIZE(nu)) :: blackbody_s ! Bnu in W/m2/sr/Hz

    !-----------------------------------------------------------------------

    blackbody_s(:) = RESHAPE(BLACKBODY_V(nu(:),[temp]),[SIZE(nu(:))])

    !-----------------------------------------------------------------------

  END FUNCTION blackbody_s


  !==========================================================================
  ! f(E)[N] = MAXWELL_BOLTZMANN(E[N],T)
  !
  !   Returns the Maxwell-Boltzmann distribution in energy [J-1].
  !==========================================================================

  PURE FUNCTION maxwell_boltzmann_v (E,T)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: E
    REAL(DP), INTENT(IN) :: T
    REAL(DP), DIMENSION(SIZE(E)) :: maxwell_boltzmann_v

    REAL(DP) :: kT

    !-----------------------------------------------------------------------

    kT = MKS%kboltz * T ! [J]
    maxwell_boltzmann_v(:) = 2._DP / SQRT(pi) * SQRT(E(:)/kT**3) * EXP(-E(:)/kT)

    !-----------------------------------------------------------------------

  END FUNCTION maxwell_boltzmann_v

  !==========================================================================

  PURE FUNCTION maxwell_boltzmann_s (E,T)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: E
    REAL(DP), INTENT(IN) :: T
    REAL(DP) :: maxwell_boltzmann_s

    !-----------------------------------------------------------------------

    maxwell_boltzmann_s = MAXVAL(MAXWELL_BOLTZMANN_V([E],T))

    !-----------------------------------------------------------------------

  END FUNCTION maxwell_boltzmann_s


END MODULE statistical_physics
