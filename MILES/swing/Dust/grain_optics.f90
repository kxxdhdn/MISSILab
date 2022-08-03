!******************************************************************************
!*
!*                     VARIOUS GRAIN OPTICAL PROPERTIES
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 08/2007.
  !    - Add Mie scattering and other related routines (02/2014), and clean-up
  !      the obsolete stuff.
  !    - Update with the HDF5 format.
  !    - 01/2016: add THEMIS data.
  !    - 07/2017: add the revised THEMIS data.
  !    - 07/2019: make read_optics more flexible (interface, Qabs label 
  !      aliases, etc.).
  ! 
  ! 3) DESCRIPTION: Functions returning the X-ray to radio grain
  !                 optical properties.
  !==========================================================================

MODULE grain_optics

  USE constants, ONLY:
  USE arrays, ONLY:
  USE inout, ONLY:
  USE interpolation, ONLY:
  USE integration, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Mie_scattering, RayleighGans_scattering, Geometric_scattering
  PUBLIC :: grain_cross_section, rho_grain, coll_cross_section
  PUBLIC :: read_optics

  INTERFACE rho_grain
    MODULE PROCEDURE rho_grain_pure, rho_grain_gen
  END INTERFACE rho_grain
  INTERFACE label_alias
    MODULE PROCEDURE label_alias_1D, label_alias_scl
  END INTERFACE label_alias
  INTERFACE read_optics
    MODULE PROCEDURE read_optics_gen, read_optics_scl
  END INTERFACE read_optics

  ! Enthalpy Label
  INTEGER, PARAMETER, PUBLIC :: lendustQ = 30 ! label length


CONTAINS


  !==========================================================================
  ! CALL MIE_SCATTERING (x,refrel,Nang,Qext,Qsca,Qabs,Qback,gsca, &
  !                      S1[2*Nang-1],S2[2*Nang-1])
  !
  !   Returns scattering and absorption by a homogeneous isotropic sphere, 
  ! as a function of X = 2*pi*a/lambda, given REFREL, its relative refractive 
  ! index (complex refractive index of the sphere divided by the real part 
  ! index of the surrounding medium). Nang>1 is the number of scattering angles 
  ! between 0 and pi/2 (will calculate 2*Nang-1 directions from 0 to pi/2).
  ! m  = SQRT(eps1+i.eps2);
  ! Qext = C_ext/pi*a**2 = efficiency factor for extinction;
  ! Qsca = C_sca/pi*a**2 = efficiency factor for scattering;
  ! Qabs = Qext - Qsca = efficiency for absorption;
  ! Qback = 4.*pi*(dC_sca/domega)/pi*a**2 =  backscattering efficiency;
  ! gsca = <cos(theta)> = asymetry parameter for scattering.
  ! S1 and S2 are diagonal elements of the amplitude scattering matrix:
  ! S1(1:2*Nang-1) = -i*f_22 (incid. E perp. to scatt. plane,
  !                           scatt. E perp. to scatt. plane);
  ! S2(1:2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
  !                           scatt. E parr. to scatt. plane);
  !
  !   Cf. the book by Bohren & Huffman (Sects. 4.4, 4.8 and appendix A).  
  !==========================================================================

  SUBROUTINE Mie_Scattering (x,refrel,Nang,Qext,Qsca,Qabs,Qback,gsca,S1,S2)

    USE utilities, ONLY: DP, CDP, strike
    USE constants, ONLY: dpi => pi
    USE arrays, ONLY: ramp
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: x
    COMPLEX(CDP), INTENT(IN) :: refrel
    INTEGER, INTENT(INOUT) :: Nang
    REAL(DP), INTENT(OUT) :: Qext, Qsca, Qabs, Qback, gsca
    COMPLEX(CDP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: S1, S2

    INTEGER, PARAMETER :: maxNang = 1000, Nmaxmax = 1500000

    INTEGER :: j, jj, n, Nstop, Nmax
    REAL(DP) :: chi, chi0, chi1, En, Fn, p, psi, psi0, psi1
    REAL(DP) :: xstop, ymod
    REAL(DP), DIMENSION(:), ALLOCATABLE :: mu, pi, pi0, pi1, tau
    COMPLEX(CDP) :: an, an1, bn, bn1, xi, xi1, y
    COMPLEX(CDP), DIMENSION(:), ALLOCATABLE :: D, SS1, SS2
      
    !-----------------------------------------------------------------------

    ! Preliminary
    !------------
    ! Scattering angles (mu = COS(theta))
    IF (Nang > maxNang) CALL STRIKE("MIE_SCATTERING","Nang > maxNang")
    IF (Nang < 2) Nang = 2
    ALLOCATE (mu(Nang))
    mu(:) = COS(RAMP(Nang,0._DP,dpi/2._DP))

    ! Initialize the variables
    y = x * refrel
    ymod = ABS(y)

    ! Termination of the series expansion (appendix A)
    xstop = x + 4._DP*x**(1._DP/3._DP) + 2._DP
    Nmax = NINT(MAX(xstop,ymod)) + 15
    Nstop = NINT(xstop)
    IF (Nmax > Nmaxmax) CALL STRIKE("MIE_SCATTERING","Nmax too big")


    ! Series expansion for Mie development
    !-------------------------------------
    ALLOCATE (pi(Nang),pi0(Nang),pi1(Nang),tau(Nang))
    pi0(:) = 0._DP
    pi1(:) = 1._DP
    ALLOCATE (SS1(2*Nang-1),SS2(2*Nang-1))
    SS1(:) = CMPLX(0._DP,0._DP,CDP)
    SS2(:) = CMPLX(0._DP,0._DP,CDP)

    ! Logarithmic derivative D calculated by downward recurrence (Eq. 4.89)
    ALLOCATE (D(Nmax))
    D(Nmax) = CMPLX(0._DP,0._DP,CDP)
    DO n=1,Nmax-1
      En = Nmax - n + 1._DP
      D(Nmax-n) = En/y - 1._DP / ( D(Nmax-n+1) + En/y )
    END DO

    ! Useless initializations to shut up the compilation warnings
    An1 = (0._DP,0._DP)
    Bn1 = (0._DP,0._DP)    

    ! Riccati-Bessel functions with real argument X calculated by upward 
    ! recurrence
    psi0 = COS(x)
    psi1 = SIN(x)
    chi0 = -SIN(x)
    chi1 = COS(x)
    xi1 = CMPLX(psi1,-chi1,CDP)
    Qsca = 0._DP
    gsca = 0._DP
    p = -1._DP
    DO n=1,Nstop
      En = n
      Fn = (2._DP*En+1._DP) / (En*(En+1._DP))

      ! For given N, psi  = psi_n        chi  = chi_n
      !              psi1 = psi_{n-1}    chi1 = chi_{n-1}
      !              psi0 = psi_{n-2}    chi0 = chi_{n-2}
      psi = (2._DP*En-1._DP) * psi1/x - psi0
      chi = (2._DP*En-1._DP) * chi1/x - chi0
      xi = CMPLX(psi,-chi,CDP)

      ! Compute the scattering coefficients (Eq. 4.88)
      IF (n > 1) THEN
        An1 = An
        Bn1 = Bn
      END IF
      An = ( (D(n)/refrel + En/x) * psi - psi1 ) &
         / ( (D(n)/refrel + En/x) * xi - xi1 )
      Bn = ( (refrel*D(n) + En/x) * psi - psi1 ) &
         / ( (refrel*D(n) + En/x) * xi - xi1 )

      ! Build the scattering quantities
      Qsca = Qsca + REAL(( 2._DP*En + 1._DP ) * ( ABS(An)**2 + ABS(Bn)**2 ), DP)
      gsca = gsca + REAL( ( (2._DP*En + 1._DP) / (En * (En + 1._DP)) ) &
                        * ( REAL(An,DP)*REAL(Bn,DP) + AIMAG(An)*AIMAG(Bn) ), DP)
      IF (n > 1) THEN
        gsca = gsca &
             + REAL( ( (En-1._DP) * (En+1._DP) / En ) &
                   * ( REAL(An1,DP)*REAL(An,DP) + AIMAG(An1)*AIMAG(An) &
             + REAL(Bn1,DP)*REAL(Bn,DP) + AIMAG(Bn1)*AIMAG(Bn) ) )
      END IF

      ! Now calculate scattering intensity pattern
      ! First do angles from 0 to 90
      DO j=1,Nang
        jj = 2*Nang - j
        pi(j) = pi1(j)
        tau(j) = En*mu(j)*pi(j) - (En+1._DP)*pi0(j)
        SS1(j) = SS1(j) + Fn * ( An*pi(j) + Bn*tau(j) )
        SS2(j) = SS2(j) + Fn * ( An*tau(j) + Bn*pi(j) )
      END DO

      ! Now do angles greater than 90 using PI and TAU from angles less than 90.
      ! P=1 for N=1,3,...; P=-1 for N=2,4,...
      p = -p
      DO j=1,Nang-1
        jj = 2*Nang - j
        SS1(jj) = SS1(jj) + Fn*p*( An*pi(j) - Bn*tau(j) )
        SS2(jj) = SS2(jj) + Fn*p*( Bn*pi(j) - An*tau(j) )
      END DO
      psi0 = psi1
      psi1 = psi
      chi0 = chi1
      chi1 = chi
      xi1 = CMPLX(psi1,-chi1,CDP)

      ! Compute pi_n for next value of n
      ! For each angle J, compute pi_n+1 from PI = pi_n , PI0 = pi_n-1
      DO j=1,Nang
        pi1(j) = ( (2._DP*En+1._DP) * mu(j)*pi(j) - (En+1._DP)*pi0(j) ) / En
        pi0(j) = pi(j)
      END DO
    
    END DO


    ! Efficiencies
    !-------------
    gsca = 2._DP * gsca / Qsca
    Qsca = 2._DP / (x*x) * Qsca
    Qext = 4._DP / (x*x) * REAL(SS1(1),DP)
    Qback = 4._DP * (ABS(SS1(2*Nang-1))/x)**2
    Qabs = Qext - Qsca

    ! Optional
    IF (PRESENT(S1) .AND. PRESENT(S2)) THEN
      ALLOCATE (S1(2*Nang-1),S2(Nang-1))
      S1(:) = SS1(:)
      S2(:) = SS2(:)
    END IF
    DEALLOCATE (mu,pi,pi0,pi1,tau,D,SS1,SS2)

    RETURN
 
    !-----------------------------------------------------------------------

  END SUBROUTINE Mie_Scattering


  !==========================================================================
  ! CALL RayleighGans_Scattering(x,refrel,Qext,Qsca,Qabs,gsca)
  !
  ! Rayleigh-Gans Regime approximation. From Laor & Draine (1993, Sect. 2.2.1).
  !==========================================================================

  SUBROUTINE RayleighGans_scattering (x,refrel,Qext,Qsca,Qabs,gsca)

    USE utilities, ONLY: DP, CDP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: x
    COMPLEX(CDP), INTENT(IN) :: refrel
    REAL(DP), INTENT(OUT) :: Qext, Qsca, Qabs, gsca

    !-----------------------------------------------------------------------

    ! Eqs. (5,7,8) of Laor & Draine (1993)
    Qabs = 8._DP/3._DP * AIMAG(refrel) * x
    Qsca = 32._DP * ABS(refrel-1)**2 * x**4 / (27._DP+16._DP*x**2)
    gsca = 0.3_DP * x**2 / (1._DP+0.3_DP*x**2) 
    Qext = Qabs + Qsca

    !-----------------------------------------------------------------------

  END SUBROUTINE RayleighGans_scattering


  !==========================================================================
  ! CALL Geometric_Scattering (x,refrel,Qext,Qsca,Qabs,gsca)
  !
  !   Geometrical Optics Regime approximation. From Laor & Draine (1993, 
  ! Sect. 2.2.2).  
  !==========================================================================

  SUBROUTINE Geometric_Scattering (x,refrel,Qext,Qsca,Qabs,gsca)

    USE utilities, ONLY: DP, CDP
    USE constants, ONLY: pi
    USE arrays, ONLY: ramp
    USE integration, ONLY: integ_tabulated
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: x
    COMPLEX(CDP), INTENT(IN) :: refrel
    REAL(DP), INTENT(OUT) :: Qext, Qsca, Qabs, gsca

    INTEGER, PARAMETER :: Nang = 100, Nn = 2
    INTEGER :: i, n
    COMPLEX(CDP) :: epsilon
    REAL(DP) :: Qabs_RG, Qsca_RG, Qext_RG
    REAL(DP), DIMENSION(Nang) :: theta, ptheta, qtheta, psi, ppsi, qpsi
    REAL(DP), DIMENSION(Nang) :: R0per, R0par, R1per, R1par, A
    REAL(DP), DIMENSION(Nang,Nn) :: hn

    !-----------------------------------------------------------------------

    ! Angle grid
    theta(:) = RAMP(Nang,0._DP,pi/2._DP)


    ! Absorption efficiency
    !----------------------
    epsilon = refrel**2

    ! 1) Interface reflections
    !  - Eq. (9) of Laor & Draine (1993)
    ptheta(:) = SQRT(0.5_DP*( (SIN(theta(:)))**2 - REAL(epsilon,DP) &
                            + SQRT( (2*AIMAG(epsilon))**2 &
                                  + (REAL(epsilon,DP)-(SIN(theta(:)))**2)**2) ))
    !  - Eq. (10) of Laor & Draine (1993)
    qtheta(:) = SQRT(0.5_DP*( REAL(epsilon,DP) - (SIN(theta(:)))**2 &
                            + SQRT((2*AIMAG(epsilon))**2 &
                                  + (REAL(epsilon,DP)-(SIN(theta(:)))**2)**2) ))
    !  - Eq. (11) of Laor & Draine (1993)
    R0per(:) = ( (qtheta(:)-COS(theta(:)))**2 + ptheta(:)**2 ) &
            / ( (qtheta(:)+COS(theta(:)))**2 + ptheta(:)**2 )
    !  - Eq. (12) of Laor & Draine (1993)
    R0par(:) = R0per(:) * ( (qtheta(:)-SIN(theta(:))*TAN(theta(:)))**2 &
                            + ptheta(:)**2 ) &
                          / ( (qtheta(:)+SIN(theta(:))*TAN(theta(:)))**2 &
                            + ptheta(:)**2 )

    ! 2) Internal reflections (theta->psi,epsilon->1/epsilon)
    !  - Eq. (14) of Laor & Draine (1993)
    psi(:) = ASIN( SIN(theta(:)) / SQRT(qtheta(:)**2+(SIN(theta(:)))**2) )

    !  - Eq. (9) of Laor & Draine (1993)
    ppsi(:) = SQRT(0.5_DP*( (SIN(psi(:)))**2 - REAL(1._DP/epsilon,DP) & 
                          + SQRT( (2*AIMAG(1._DP/epsilon))**2 &
                              + (REAL(1._DP/epsilon,DP)-(SIN(psi(:)))**2)**2) ))
    !  - Eq. (10) of Laor & Draine (1993)
    qpsi(:) = SQRT(0.5_DP*( REAL(1._DP/epsilon,DP) - (SIN(psi(:)))**2 & 
                          + SQRT( (2*AIMAG(1._DP/epsilon))**2 &
                              + (REAL(1._DP/epsilon,DP)-(SIN(psi(:)))**2)**2) ))
    !  - Eq. (11) of Laor & Draine (1993)
    R1per(:) = ( (qpsi(:)-COS(psi(:)))**2 + ppsi(:)**2 ) &
             / ( (qpsi(:)+COS(psi(:)))**2 + ppsi(:)**2 )
    !  - Eq. (12) of Laor & Draine (1993)
    R1par(:) = R1per(:) * ((qpsi(:)-SIN(psi(:))*TAN(psi(:)))**2 + ppsi(:)**2) &
             / ( (qpsi(:)+SIN(psi(:))*TAN(psi(:)))**2 + ppsi(:)**2 )

    ! Eq. (16) of Laor & Draine (1993)
    A(:) = EXP( -4._DP*x * AIMAG(refrel) * COS(psi(:)) )

    ! Eq. (15) of Laor & Draine (1993)
    Qabs = INTEG_TABULATED( theta(:), &
                            SIN(theta(:))*COS(theta(:))*(1._DP-A(:)) &
                            *( (1._DP-R0per(:))/(1._DP-A(:)*R1per(:)) &
                             + (1._DP-R0par(:))/(1._DP-A(:)*R1par(:)) ), &
                            XRANGE=[0._DP,pi/2._DP] )


    ! Extinction efficiency
    !----------------------
    ! Rayleigh-Gans approximation
    Qabs_RG = 8._DP/3._DP * AIMAG(refrel)*x
    Qsca_RG = 32._DP*(ABS(refrel-1._DP))**2 * x**4 / (27._DP+16._DP*x**2)
    Qext_RG = Qabs_RG + Qsca_RG

    ! Eq. (17) of Laor & Draine (1993)
    Qext = Qext_RG / SQRT(1._DP+0.25_DP*Qext_RG**2)


    ! Scattering terms
    !-----------------
    ! Scattering efficiency
    Qsca = Qext - Qabs

    ! Asymmetry
    DO i=1,Nn
      n = i - 1
      hn(:,i) = A(:)**(n+1._DP) * COS(2*(theta(:)-psi(:)) + n*(pi-2*psi(:))) &
              * ( (1._DP-R0per(:))*R1per(:)**n*(1._DP-R1per(:)) &
                  + (1._DP-R0par(:))*R1par(:)**n*(1._DP-R1par(:)) )
    END DO
    gsca = INTEG_TABULATED( theta(:), &
                            SIN(theta(:))*COS(theta(:)) &
                            * ( (R0per(:)+R0par(:))*COS(2*theta(:)) &
                              + SUM(hn(:,:),DIM=2) ) ) / Qsca

    !-----------------------------------------------------------------------

  END SUBROUTINE Geometric_scattering


  !==========================================================================
  ! CALL GRAIN_CROSS_SECTION (radius,wave,refrel,Qext,Qsca,Qabs,gsca)
  !
  !   Compute the wavelength dependent cross-sections, for a set of sizes,
  ! using the method of Laor & Draine (1993; three different regimes).
  !==========================================================================

  SUBROUTINE grain_cross_section (radius,wave,refrel,Qext,Qsca,Qabs,gsca, &
                                  timestr)

    USE utilities, ONLY: DP, CDP, trimlr, pring, verbatim, timinfo, time_type
    USE constants, ONLY: pi
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: radius, wave
    COMPLEX(CDP), DIMENSION(:), INTENT(IN) :: refrel
    REAL(DP), DIMENSION(SIZE(radius),SIZE(wave)), INTENT(OUT) :: Qext, Qsca
    REAL(DP), DIMENSION(SIZE(radius),SIZE(wave)), INTENT(OUT) :: Qabs, gsca
    TYPE(time_type), INTENT(INOUT), OPTIONAL :: timestr

    INTEGER :: i, j, Nw, Nr, Nang
    REAL(DP) :: Qext_wr, Qsca_wr, Qabs_wr, Qback_wr, gsca_wr
    REAL(DP), DIMENSION(SIZE(wave)) :: x, mx, m_1x
    LOGICAL, PARAMETER :: onlyMie = .True.

    !-----------------------------------------------------------------------

    Nw = SIZE(wave(:))
    Nr = SIZE(radius(:))
    Nang = 2

    ! Loop on the wavelengths
    Cgrain_size: DO i=1,Nr
      IF (verbatim .AND. (MODULO(i,10) == 1 .OR. i == Nr)) THEN
        IF (PRESENT(timestr)) THEN
          PRINT*, "  Calculating grain cross-section for a=" &
                           //TRIMLR(PRING(radius(i),NDEC=6))//" microns" &
                //" ("//TRIMLR(TIMINFO(timestr))//")"
        ELSE
          PRINT*, "  Calculating grain cross-section for a=" &
                           //TRIMLR(PRING(radius(i),NDEC=6))//" microns"
        END IF
      END IF      

      ! Size parameter
      x(:) = 2._DP*pi*radius(i) / wave(:)
      mx(:) = ABS(refrel(:)) * x(:)
      m_1x(:) = ABS(refrel(:)-1._DP) * x(:)

      ! Loop on the sizes
      Cgrain_wave: DO j=1,Nw

        ! Actual calculation
        IF (onlyMie) THEN
          ! We implement Mie everywhere, even in the R-G and the geometric 
          ! optics regimes
          CALL MIE_SCATTERING(x(j),refrel(j),Nang,Qext_wr,Qsca_wr,Qabs_wr, &
                              Qback_wr,gsca_wr)
        ELSE
          ! Laor & Draine (1993) three regime method
          IF (mx(j) <= 1.E3_DP) THEN
            CALL MIE_SCATTERING(x(j),refrel(j),Nang,Qext_wr,Qsca_wr,Qabs_wr, &
                                Qback_wr,gsca_wr)
          ELSE IF (mx(j) > 1.E3_DP .AND. m_1x(j) <= 1.E-3_DP) THEN
            CALL RAYLEIGHGANS_SCATTERING(x(j),refrel(j),Qext_wr,Qsca_wr, &
                                         Qabs_wr,gsca_wr)
          ELSE IF (mx(j) > 1.E3_DP .AND. m_1x(j) > 1.E-3_DP) THEN
            CALL GEOMETRIC_SCATTERING(x(j),refrel(j),Qext_wr,Qsca_wr,Qabs_wr, &
                                      gsca_wr)
          END IF   
        END IF
        
        ! Final storage
        Qext(i,j) = Qext_wr
        Qsca(i,j) = Qsca_wr
        Qabs(i,j) = Qabs_wr
        gsca(i,j) = gsca_wr

      END DO Cgrain_wave

    END DO Cgrain_size

    !-----------------------------------------------------------------------

  END SUBROUTINE grain_cross_section


  !==========================================================================
  ! sigma = COLL_CROSS_SECTION(E,a)
  !
  !    Collision cross-section for electrons, following Eq. (4) of Dwek (1986).
  ! E must be in J, a in microns, sigma is returned in m2.
  !==========================================================================

  PURE FUNCTION coll_cross_section (E,a)

    USE utilities, ONLY: DP
    USE constants, ONLY: pi, MKS
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: E
    REAL(DP), INTENT(IN) :: a
    REAL(DP), DIMENSION(SIZE(E)) :: coll_cross_section

    REAL(DP) :: Estar

    !-----------------------------------------------------------------------

    Estar = 3.7E-15_DP * a**(2._DP/3._DP) ! barrier energy [J]
    WHERE (E(:) > Estar)
      coll_cross_section(:) = pi * (a*MKS%micron)**2 &
        * ( 1._DP - ( 1._DP - (Estar/E(:))**1.5_DP )**(2._DP/3._DP) )
    ELSEWHERE
      coll_cross_section(:) = pi * (a*MKS%micron)**2
    END WHERE

    !-----------------------------------------------------------------------

  END FUNCTION coll_cross_section


  !==========================================================================
  ! label = LABEL_ALIAS(lab)
  !
  !  Aliases for common optical properties.
  !==========================================================================

  PURE FUNCTION label_alias_1D(lab)

    USE utilities, ONLY: trimeq
    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:), INTENT(IN) :: lab
    CHARACTER(lendustQ), DIMENSION(SIZE(lab)) :: label_alias_1D
    INTEGER :: i, Nspec

    !-----------------------------------------------------------------------

    Nspec = SIZE(lab(:))
    label_alias_1D(:) = lab(:)
    DO i=1,Nspec
      IF (TRIMEQ(lab(i),"Sil") .OR. TRIMEQ(lab(i),"silicate")) THEN 
        label_alias_1D(i) = "Sil_D03"
      ELSE IF (TRIMEQ(lab(i),"SiC")) THEN
        label_alias_1D(i) = "SiC_LD93"
      ELSE IF (TRIMEQ(lab(i),"Gra") .OR. TRIMEQ(lab(i),"graphite")) THEN
        label_alias_1D(i) = "Gra_D03"
      ELSE IF (TRIMEQ(lab(i),"PAHn") .OR. TRIMEQ(lab(i),"PAH0")) THEN
        label_alias_1D(i) = "PAHn_DL07"
      ELSE IF (TRIMEQ(lab(i),"PAHi") .OR. TRIMEQ(lab(i),"PAH+")) THEN
        label_alias_1D(i) = "PAHi_DL07"
      ELSE IF (TRIMEQ(lab(i), "a-C(:H)") .OR. TRIMEQ(lab(i),"HAC")) THEN
        label_alias_1D(i) = "a-C_man20nm_small"
      END IF
    END DO

    !-----------------------------------------------------------------------

  END FUNCTION label_alias_1D

  !-------------------------------------------------------------------------

  PURE FUNCTION label_alias_scl(lab)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: lab
    CHARACTER(lendustQ) :: label_alias_scl
    CHARACTER(lendustQ), DIMENSION(1) :: tmp
    tmp = LABEL_ALIAS_1D([lab])
    label_alias_scl = tmp(1)
  END FUNCTION label_alias_scl

  !==========================================================================
  ! rho = rho_grain(species,name,reference)
  !
  !   Returns the mass per unit volume of material for grain species in kg/m3.
  ! The label of the species should be the same as for the cross-sections.
  !==========================================================================

  PURE FUNCTION rho_grain_pure (species)

    USE utilities, ONLY: DP, trimeq
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: species
    REAL(DP) :: rho_grain_pure
    CHARACTER(lendustQ) :: species0

    !-----------------------------------------------------------------------

    species0 = (LABEL_ALIAS(species))
    IF (TRIMEQ(species0,"Sil_LD93")) THEN
      rho_grain_pure = 3.3_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"Sil_WD01")) THEN
      rho_grain_pure = 3.5_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"Sil_LD01")) THEN
      rho_grain_pure = 3.5_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"Sil_D03")) THEN
      rho_grain_pure = 3.5_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"SiC_LD93")) THEN
      rho_grain_pure = 3.22_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"Gra_D03")) THEN
      rho_grain_pure = 2.24_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"Gra_LD93")) THEN
      rho_grain_pure = 2.24_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"ACAR_Z96")) THEN
      rho_grain_pure = 1.85_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"ACH2_Z96")) THEN
      rho_grain_pure = 1.85_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"BE_Z96")) THEN
      rho_grain_pure = 1.85_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"PAHn_LD01")) THEN
      rho_grain_pure = 2.24_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"PAHi_LD01")) THEN
      rho_grain_pure = 2.24_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"PAHn_DL07")) THEN
      rho_grain_pure = 2.24_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"PAHi_DL07")) THEN
      rho_grain_pure = 2.24_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"PAHn_DL07_G11")) THEN
      rho_grain_pure = 2.24_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"PAHi_DL07_G11")) THEN
      rho_grain_pure = 2.24_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"PAHn_DL07_C11")) THEN
      rho_grain_pure = 2.24_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"PAHi_DL07_C11")) THEN
      rho_grain_pure = 2.24_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"a-C_man20nm_J13")) THEN
      rho_grain_pure = 1.30_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"a-Forst_Fe_man5nm_J13")) THEN
      rho_grain_pure = 1.60_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"a-Enst_Fe_man5nm_J13")) THEN
      rho_grain_pure = 1.60_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"a-C_man20nm_small")) THEN
      rho_grain_pure = 1.60_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"a-C_man20nm_big")) THEN
      rho_grain_pure = 1.51_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"a-Forst_Fe_man5nm")) THEN
      rho_grain_pure = 2.19_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"a-Enst_Fe_man5nm")) THEN
      rho_grain_pure = 2.19_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"Sil_Mg07_J03")) THEN
      rho_grain_pure = 3.50_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"Sil_Mg10_J03")) THEN
      rho_grain_pure = 3.50_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"Sil_Mg15_J03")) THEN
      rho_grain_pure = 3.50_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"Sil_Mg20_J03")) THEN
      rho_grain_pure = 3.50_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"Sil_Mg24_J03")) THEN
      rho_grain_pure = 3.50_DP * 1.E3_DP ! [kg/m3]
    ELSE IF (TRIMEQ(species0,"FeO_H95")) THEN
      rho_grain_pure = 5.70_DP * 1.E3_DP ! [kg/m3]
    ELSE
      rho_grain_pure = 0._DP
    END IF

    !-----------------------------------------------------------------------

  END FUNCTION rho_grain_pure

  !==========================================================================

  FUNCTION rho_grain_gen (species,name,reference)

    USE utilities, ONLY: DP, trimeq
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: species
    CHARACTER(*), INTENT(OUT) :: name, reference
    REAL(DP) :: rho_grain_gen
    CHARACTER(lendustQ) :: species0

    !-----------------------------------------------------------------------

    species0 = LABEL_ALIAS(species)
    IF (TRIMEQ(species0,"Sil_LD93")) THEN
      name = "Original astronomical silicate"
      reference = "(Laor & Draine 1993, ApJ, 402, 441)"
    ELSE IF (TRIMEQ(species0,"Sil_WD01")) THEN
      name = "Smoothed UV astronomical silicate"
      reference = "(Weingartner & Draine 2001, ApJ, 548, 296)"
    ELSE IF (TRIMEQ(species0,"Sil_LD01")) THEN
      name = "Smoothed UV astronomical silicate & submm excess"
      reference = "(Li & Draine 2001, ApJ, 554, 778)"
    ELSE IF (TRIMEQ(species0,"Sil_D03")) THEN
      name = "Astronomical silicate"
      reference = "(Draine 2003, ApJ, 598, 1017)"
    ELSE IF (TRIMEQ(species0,"SiC_LD93")) THEN
      name = "Silicon carbide"
      reference = "(Laor & Draine 1993, ApJ, 402, 441)"
    ELSE IF (TRIMEQ(species0,"Gra_D03")) THEN
      name = "Graphitic carbon (1/3-2/3)"
      reference = "(Draine 2003, ApJ, 598, 1017)"
    ELSE IF (TRIMEQ(species0,"Gra_LD93")) THEN
      name = "Original graphitic carbon (1/3-2/3)"
      reference = "(Laor & Draine 1993, ApJ, 402, 441)"
    ELSE IF (TRIMEQ(species0,"ACAR_Z96")) THEN
      name = "ACAR amorphous carbon"
      reference = "(Zubko et al. 1996, MNRAS, 282, 1321)"
    ELSE IF (TRIMEQ(species0,"ACH2_Z96")) THEN
      name = "ACH2 amorphous carbon"
      reference = "(Zubko et al. 1996, MNRAS, 282, 1321)"
    ELSE IF (TRIMEQ(species0,"BE_Z96")) THEN
      name = "BE amorphous carbon"
      reference = "(Zubko et al. 1996, MNRAS, 282, 1321)"
    ELSE IF (TRIMEQ(species0,"PAHn_LD01")) THEN
      name = "Original neutral PAHs"
      reference = "(Li & Draine 2001, ApJ, 554, 778)"
    ELSE IF (TRIMEQ(species0,"PAHi_LD01")) THEN
      name = "Original ionized PAHs"
      reference = "(Li & Draine 2001, ApJ, 554, 778)"
    ELSE IF (TRIMEQ(species0,"PAHn_DL07")) THEN
      name = "Neutral PAHs"
      reference = "(Draine & Li 2007, ApJ, 657, 810)"
    ELSE IF (TRIMEQ(species0,"PAHi_DL07")) THEN
      name = "Ionized PAHs"
      reference = "(Draine & Li 2007, ApJ, 657, 810)"
    ELSE IF (TRIMEQ(species0,"PAHn_DL07_G11")) THEN
      name = "Neutral PAHs modified for the Z04 model"
      reference = "(Galliano et al. 2011, A&A, 536, A88)"
    ELSE IF (TRIMEQ(species0,"PAHi_DL07_G11")) THEN
      name = "Ionized PAHs modified for the Z04 model"
      reference = "(Galliano et al. 2011, A&A, 536, A88)"
    ELSE IF (TRIMEQ(species0,"PAHn_DL07_C11")) THEN
      name = "Neutral PAHs modified for the C11 model"
      reference = "(Compiegne et al. 2011, A&A, 525, A103)"
    ELSE IF (TRIMEQ(species0,"PAHi_DL07_C11")) THEN
      name = "Ionized PAHs modified for the C11 model"
      reference = "(Compiegne et al. 2011, A&A, 525, A103)"
    ELSE IF (TRIMEQ(species0,"a-C_man20nm_J13")) THEN
      name = "a-C(:H) + 20 nm a-C"
      reference = "(Jones et al. 2013, A&A, 558, A62)"
    ELSE IF (TRIMEQ(species0,"a-Forst_Fe_man5nm_J13")) THEN
      name = "a-Forsterite + 10 % Fe/FeS + 5 nm a-C"
      reference = "(Jones et al. 2013, A&A, 558, A62)"
    ELSE IF (TRIMEQ(species0,"a-Enst_Fe_man5nm_J13")) THEN
      name = "a-Enstatite + 10 % Fe/FeS + 5 nm a-C"
      reference = "(Jones et al. 2013, A&A, 558, A62)"
    ELSE IF (TRIMEQ(species0,"a-C_man20nm_small")) THEN
      name = "Small a-C(:H) + 20 nm a-C"
      reference = "(Jones et al. 2017, A&A, 602, A46)"
    ELSE IF (TRIMEQ(species0,"a-C_man20nm_big")) THEN
      name = "Big a-C(:H) + 20 nm a-C"
      reference = "(Jones et al. 2017, A&A, 602, A46)"
    ELSE IF (TRIMEQ(species0,"a-Forst_Fe_man5nm")) THEN
      name = "a-Forsterite + 10 % Fe/FeS + 5 nm a-C"
      reference = "(Jones et al. 2017, A&A, 602, A46)"
    ELSE IF (TRIMEQ(species0,"a-Enst_Fe_man5nm")) THEN
      name = "a-Enstatite + 10 % Fe/FeS + 5 nm a-C"
      reference = "(Jones et al. 2017, A&A, 602, A46)"
    ELSE IF (TRIMEQ(species0,"Sil_Mg07_J03")) THEN
      name = "Silicate Mg_0.7 Si O_2.7"
      reference = "(Jaeger et al. 2003, A&A, 408, 193)"
    ELSE IF (TRIMEQ(species0,"FeO_H95")) THEN
      name = "Iron Oxide FeO"
      reference = "(Henning et al. 1995, A&AS, 112, 143)"
    END IF
    rho_grain_gen = RHO_GRAIN_PURE(species0)

    !-----------------------------------------------------------------------

  END FUNCTION rho_grain_gen


  !==========================================================================
  ! CALL READ_OPTICS(lab,Nwall,Nrall,NTall,Nmax,Nrmax,NTmax, & 
  !                             waveall,radiusall,tempall, &
  !                             Qabsall,Qscaall,gscall,Qavall, & 
  !                             wrange,wave,nu,Nw)
  !
  !   Reads the wavelength, radius and temperature dependent precomputed
  ! optical properties of several species into one big array per quantity,
  ! and optionally reinterpolate them on a common wavelength grid, without
  ! loosing spectral resolution.
  !==========================================================================

  SUBROUTINE read_optics_gen (lab,Nwall,Nrall,NTall,Nwmax,Nrmax,NTmax, & 
                               waveall,radiusall,tempall, &
                               Qabsall,Qscaall,gscaall,Qavall, & 
                               wrange,wave,nu,Nw)
  
    USE utilities, ONLY: DP, trimeq, trimlr, lenmax, temproot, strreplace
    USE constants, ONLY: MKS
    USE arrays, ONLY: reallocate, closest
    USE inout, ONLY: read_hdf5, h5ext
    USE interpolation, ONLY: interp_lin_sorted
    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:), INTENT(IN) :: lab
    INTEGER, DIMENSION(SIZE(lab)), INTENT(OUT), OPTIONAL :: Nwall, Nrall, NTall
    INTEGER, INTENT(OUT), OPTIONAL :: Nwmax, Nrmax, NTmax
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: waveall
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: radiusall
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: tempall
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Qabsall
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Qscaall
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Qavall
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: gscaall
    REAL(DP), DIMENSION(2), INTENT(IN), OPTIONAL :: wrange
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: wave, nu
    INTEGER, INTENT(OUT), OPTIONAL :: Nw

    INTEGER :: i, j, Nspec, Nw1
    INTEGER :: Nwmax0, Nrmax0, NTmax0, Nwmax1, Nrmax1, NTmax1, Nw0, Nr0, NT0
    INTEGER, DIMENSION(SIZE(lab)) :: Nwall0, Nrall0, NTall0
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind
    REAL(DP) :: lnwmin, lnwmax, dlnwave, dlnwmax
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wave_r, radius_r, temp_r, wave0
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lnwave, lnwave0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dlnwaveall
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: waveall0, radiusall0, tempall0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: wave1, radius1, temp1
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Qabs1, Qsca1, gsca1, Qav1
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Qabs_r, Qsca_r, Qav_r, gsca_r
    CHARACTER(lenmax) :: dirQabs
    CHARACTER(lendustQ), DIMENSION(SIZE(lab)) :: lab0
    LOGICAL :: wavinterp
    LOGICAL, DIMENSION(:), ALLOCATABLE :: boolind

    !-----------------------------------------------------------------------

    ! 1) Preliminaries
    !-----------------
    Nspec = SIZE(lab(:))
    wavinterp = .True.
    IF (.NOT. PRESENT(wave)) wavinterp = .False.

    ! Solve label degeneracies
    lab0(:) = LABEL_ALIAS(lab(:))


    ! 2) Concatenate the optical properties
    !--------------------------------------
    dirQabs = TRIMLR(temproot())//"Cross_sections/Data/"
    readoptics: DO i=1,Nspec
      nocomp: IF (.NOT. TRIMEQ(lab0(i),"")) THEN

        ! Load the variables
        CALL READ_HDF5(DBLARR1D=wave_r,NAME="Wavelength (microns)", &
                       FILE=TRIMLR(dirQabs)//"optics_"//TRIMLR(lab0(i))//h5ext)
        CALL READ_HDF5(DBLARR1D=radius_r,NAME="Grain radius (microns)", &
                       FILE=TRIMLR(dirQabs)//"optics_"//TRIMLR(lab0(i))//h5ext)
        CALL READ_HDF5(DBLARR1D=temp_r,NAME="Grain temperature (K)", &
                       FILE=TRIMLR(dirQabs)//"optics_"//TRIMLR(lab0(i))//h5ext)
        IF (PRESENT(Qabsall)) &
          CALL READ_HDF5(DBLARR2D=Qabs_r, &
                        NAME="Grain absorption efficiency (Qabs)", &
                        FILE=TRIMLR(dirQabs)//"optics_"//TRIMLR(lab0(i))//h5ext)
        IF (PRESENT(Qscaall)) &
          CALL READ_HDF5(DBLARR2D=Qsca_r, &
                        NAME="Grain scattering efficiency (Qsca)", &
                        FILE=TRIMLR(dirQabs)//"optics_"//TRIMLR(lab0(i))//h5ext)
        IF (PRESENT(gscaall)) &
          CALL READ_HDF5(DBLARR2D=gsca_r, &
                        NAME="Grain asymmetry parameter (<cos(theta)>)", &
                        FILE=TRIMLR(dirQabs)//"optics_"//TRIMLR(lab0(i))//h5ext)
        IF (PRESENT(Qavall)) &
          CALL READ_HDF5(DBLARR2D=Qav_r,NAME="Planck average of Qabs", &
                        FILE=TRIMLR(dirQabs)//"optics_"//TRIMLR(lab0(i))//h5ext)
        Nw0 = SIZE(wave_r(:))
        Nr0 = SIZE(radius_r(:))
        NT0 = SIZE(temp_r(:))
        Nwall0(i) = Nw0
        Nrall0(i) = Nr0
        NTall0(i) = NT0

        ! Put everything in the same array
        IF (i == 1) THEN 
          ALLOCATE( waveall0(Nspec,Nw0), radiusall0(Nspec,Nr0), &
                    tempall0(Nspec,NT0) )
          waveall0(i,:) = wave_r(:)
          radiusall0(i,:) = radius_r(:)
          tempall0(i,:) = temp_r(:)
          IF (PRESENT(Qabsall)) THEN
            ALLOCATE(Qabsall(Nspec,Nr0,Nw0))
            Qabsall(i,:,:) = Qabs_r(:,:)
          END IF
          IF (PRESENT(Qscaall)) THEN
            ALLOCATE(Qscaall(Nspec,Nr0,Nw0))
            Qscaall(i,:,:) = Qsca_r(:,:)
          END IF
          IF (PRESENT(gscaall)) THEN
            ALLOCATE(gscaall(Nspec,Nr0,Nw0))
            gscaall(i,:,:) = gsca_r(:,:)
          END IF
          IF (PRESENT(Qavall)) THEN
            ALLOCATE(Qavall(Nspec,Nr0,Nw0))
            Qavall(i,:,:) = Qav_r(:,:)
          END IF
          Nwmax0 = Nw0
          Nrmax0 = Nr0
          NTmax0 = NT0
        ELSE
          ! Concatenate wave
          Nwmax1 = Nwmax0
          Nwmax0 = MAX(Nwmax0,Nw0)
          CALL REALLOCATE(wave1,Nspec,Nwmax0)
          wave1(:,:) = 0._DP
          wave1(1:i-1,1:Nwmax1) = waveall0(1:i-1,:)
          IF (Nwmax0 > Nwmax1) wave1(1:i-1,Nwmax1+1:Nwmax0) = 0._DP
          wave1(i,1:Nw0) = wave_r(:)
          IF (Nw0 > Nwmax0) wave1(i,Nw0+1:Nwmax0) = 0._DP
          CALL REALLOCATE(waveall0,Nspec,Nwmax0)
          waveall0(:,:) = wave1(:,:)
          DEALLOCATE (wave1)
          ! Concatenate radius
          Nrmax1 = Nrmax0
          Nrmax0 = MAX(Nrmax0,Nr0)
          CALL REALLOCATE(radius1,Nspec,Nrmax0)
          radius1(:,:) = 0._DP
          radius1(1:i-1,1:Nrmax1) = radiusall0(1:i-1,:)
          radius1(i,1:Nr0) = radius_r(:)
          CALL REALLOCATE(radiusall0,Nspec,Nrmax0)
          radiusall0(:,:) = radius1(:,:)
          DEALLOCATE(radius1)
          ! Concatenate temp
          NTmax1 = NTmax0
          NTmax0 = MAX(NTmax0,NT0)
          CALL REALLOCATE(temp1,Nspec,NTmax0)
          temp1(:,:) = 0._DP
          temp1(1:i-1,1:NTmax1) = tempall0(1:i-1,:)
          temp1(i,1:NT0) = temp_r(:)
          CALL REALLOCATE(tempall0,Nspec,NTmax0)
          tempall0(:,:) = temp1(:,:)
          DEALLOCATE(temp1)
          ! Concatenate Qabs
          IF (PRESENT(Qabsall)) THEN
            CALL REALLOCATE(Qabs1,Nspec,Nrmax0,Nwmax0)
            Qabs1(:,:,:) = 0._DP
            Qabs1(1:i-1,1:Nrmax1,1:Nwmax1) = Qabsall(1:i-1,:,:)
            Qabs1(i,1:Nr0,1:Nw0) = Qabs_r(:,:)
            CALL REALLOCATE(Qabsall,Nspec,Nrmax0,Nwmax0)
            Qabsall(:,:,:) = Qabs1(:,:,:)
          END IF
          ! Concatenate Qsca
          IF (PRESENT(Qscaall)) THEN
            CALL REALLOCATE(Qsca1,Nspec,Nrmax0,Nwmax0)
            Qsca1(:,:,:) = 0._DP
            Qsca1(1:i-1,1:Nrmax1,1:Nwmax1) = Qscaall(1:i-1,:,:)
            Qsca1(i,1:Nr0,1:Nw0) = Qsca_r(:,:)
            CALL REALLOCATE(Qscaall,Nspec,Nrmax0,Nwmax0)
            Qscaall(:,:,:) = Qsca1(:,:,:)
          END IF
          ! Concatenate gsca
          IF (PRESENT(gscaall)) THEN
            CALL REALLOCATE(gsca1,Nspec,Nrmax0,Nwmax0)
            gsca1(:,:,:) = 0._DP
            gsca1(1:i-1,1:Nrmax1,1:Nwmax1) = gscaall(1:i-1,:,:)
            gsca1(i,1:Nr0,1:Nw0) = gsca_r(:,:)
            CALL REALLOCATE(gscaall,Nspec,Nrmax0,Nwmax0)
            gscaall(:,:,:) = gsca1(:,:,:)
          END IF
          ! Concatenate Qav
          IF (PRESENT(Qavall)) THEN
            CALL REALLOCATE(Qav1,Nspec,Nrmax0,NTmax0)
            Qav1(:,:,:) = 0._DP
            Qav1(1:i-1,1:Nrmax1,1:NTmax1) = Qavall(1:i-1,:,:)
            Qav1(i,1:Nr0,1:NT0) = Qav_r(:,:)
            CALL REALLOCATE(Qavall,Nspec,Nrmax0,NTmax0)
            Qavall(:,:,:) = Qav1(:,:,:)
          END IF
        END IF

        ! Free memory space
        IF (ALLOCATED(wave_r)) DEALLOCATE (wave_r)
        IF (ALLOCATED(radius_r)) DEALLOCATE(radius_r)
        IF (ALLOCATED(temp_r)) DEALLOCATE(temp_r)
        IF (ALLOCATED(Qabs_r)) DEALLOCATE(Qabs_r)
        IF (ALLOCATED(Qsca_r)) DEALLOCATE(Qsca_r)
        IF (ALLOCATED(gsca_r)) DEALLOCATE(gsca_r)

      ELSE

        Nwall0(i) = 1
        Nrall0(i) = 1
        NTall0(i) = 1
        waveall0(i,:) = 0._DP
        radiusall0(i,:) = 0._DP
        tempall0(i,:) = 0._DP
        IF (PRESENT(Qabsall)) Qabsall(i,:,:) = 0._DP
        IF (PRESENT(Qscaall)) Qscaall(i,:,:) = 0._DP
        IF (PRESENT(gscaall)) gscaall(i,:,:) = 0._DP
        IF (PRESENT(Qavall)) Qavall(i,:,:) = 0._DP

      END IF nocomp
  
    END DO readoptics

    
    ! 3) Interpolate on a common wavelength grid
    !-------------------------------------------
    commonwave: IF (Nspec > 1 .AND. wavinterp) THEN
       
      ! Wavelength range
      IF (PRESENT(wrange)) THEN
        lnwmin = LOG(wrange(1))
        lnwmax = LOG(wrange(2))
      ELSE
        lnwmin = LOG(1.E-2_DP)
        lnwmax = LOG(3.E4_DP)
      END IF

      ! Design a common grid, without loosing resolution
      ALLOCATE (dlnwaveall(Nspec,Nwmax0))
      dlnwaveall(:,1:Nwmax0-1) = LOG(waveall0(:,2:Nwmax0)) &
                               - LOG(waveall0(:,1:Nwmax0-1))
      dlnwaveall(:,Nwmax0) = dlnwaveall(:,Nwmax0-1)
      Nw0 = 1
      ALLOCATE (lnwave(Nw0))
      ALLOCATE (ind(Nspec))
      ALLOCATE (boolind(Nspec))
      lnwave(1) = lnwmin
      dlnwmax = MAXVAL(dlnwaveall(:,:),MASK=(waveall0(:,:) > 0._DP))
      DO 
        Nw1 = Nw0
        FORALL (i=1:Nspec,.NOT. TRIMEQ(lab0(i),"")) &
          ind(i) = CLOSEST(waveall0(i,1:Nwall0(i)),EXP(lnwave(Nw1)))
        boolind(:) = ( ind(:) < Nwall0(:) )
        dlnwave = dlnwmax
        DO i=1,Nspec
          IF (boolind(i) .AND. (.NOT. TRIMEQ(lab0(i),""))) THEN
            IF (dlnwave > dlnwaveall(i,ind(i))) dlnwave = dlnwaveall(i,ind(i))
          END IF
        END DO
        Nw0 = Nw1 + 1
        CALL REALLOCATE(lnwave0,Nw0)
        lnwave0(1:Nw1) = lnwave(:)
        lnwave0(Nw0) = lnwave(Nw1) + dlnwave
        CALL REALLOCATE(lnwave,Nw0)
        lnwave(:) = lnwave0(:)
        IF (lnwave(Nw0) >= lnwmax-dlnwave/2._DP) THEN
          lnwave(Nw0) = lnwmax
          EXIT
        END IF
      END DO
      CALL REALLOCATE (wave0,Nw0)
      wave0(:) = EXP(lnwave(:))
      DEALLOCATE (ind,boolind)

      ! Interpolate the Qabs
      IF (PRESENT(Qabsall)) THEN
        CALL REALLOCATE(Qabs1,Nspec,Nrmax0,Nw0)
        Qabs1(:,:,:) = 0._DP
        FORALL ( i=1:Nspec, j=1:Nrmax0, &
                (j <= Nrall0(i)) .AND. (.NOT. TRIMEQ(lab0(i),"")) ) &
          Qabs1(i,j,:) = INTERP_LIN_SORTED( Qabsall(i,j,1:Nwall0(i)), &
                                            waveall0(i,1:Nwall0(i)), wave0(:), &
                                            XLOG=.True., YLOG=.True. )
        CALL REALLOCATE(Qabsall,Nspec,Nrmax0,Nw0)
        Qabsall(:,:,:) = Qabs1
      END IF

      ! Interpolate the Qsca
      IF (PRESENT(Qscaall)) THEN
        CALL REALLOCATE(Qsca1,Nspec,Nrmax0,Nw0)
        Qsca1(:,:,:) = 0._DP
        FORALL ( i=1:Nspec, j=1:Nrmax0, &
                (j <= Nrall0(i)) .AND. (.NOT. TRIMEQ(lab0(i),"")) ) &
          Qsca1(i,j,:) = INTERP_LIN_SORTED( Qscaall(i,j,1:Nwall0(i)), &
                                            waveall0(i,1:Nwall0(i)), wave0(:), &
                                            XLOG=.True., YLOG=.True. )
        CALL REALLOCATE(Qscaall,Nspec,Nrmax0,Nw0)
        Qscaall(:,:,:) = Qsca1
      END IF

      ! Interpolate the gsca
      IF (PRESENT(gscaall)) THEN
        CALL REALLOCATE(gsca1,Nspec,Nrmax0,Nw0)
        gsca1(:,:,:) = 0._DP
        FORALL ( i=1:Nspec, j=1:Nrmax0, &
                (j <= Nrall0(i)) .AND. (.NOT. TRIMEQ(lab0(i),"")) ) &
          gsca1(i,j,:) = INTERP_LIN_SORTED( gscaall(i,j,1:Nwall0(i)), &
                                            waveall0(i,1:Nwall0(i)), wave0(:), &
                                            XLOG=.True., YLOG=.False. )
        CALL REALLOCATE(gscaall,Nspec,Nrmax0,Nw0)
        gscaall(:,:,:) = gsca1
      END IF

      ! Free memory space 
      IF (ALLOCATED(Qabs1)) DEALLOCATE(Qabs1)
      IF (ALLOCATED(Qsca1)) DEALLOCATE(Qsca1)
      IF (ALLOCATED(gsca1)) DEALLOCATE(gsca1)
      IF (ALLOCATED(Qav1)) DEALLOCATE(Qav1)
      IF (ALLOCATED(lnwave)) DEALLOCATE(lnwave)
      IF (ALLOCATED(dlnwaveall)) DEALLOCATE(dlnwaveall)
      IF (ALLOCATED(lnwave0)) DEALLOCATE(lnwave0)

    ELSE IF (Nspec == 1 .AND. wavinterp) THEN

      CALL REALLOCATE (wave0,Nwmax0)
      wave0(:) = waveall0(1,:)

    END IF commonwave


    ! 4) Optional output
    !-------------------
    IF (PRESENT(Nwmax)) Nwmax = Nwmax0
    IF (PRESENT(Nrmax)) Nrmax = Nrmax0
    IF (PRESENT(NTmax)) NTmax = NTmax0
    IF (PRESENT(Nwall)) Nwall(:) = Nwall0(:)
    IF (PRESENT(Nrall)) Nrall(:) = Nrall0(:)
    IF (PRESENT(NTall)) NTall(:) = NTall0(:)
    IF (PRESENT(waveall)) THEN
      ALLOCATE(waveall(Nspec,Nwmax0))
      waveall(:,:) = waveall0(:,:)
    END IF
    IF (PRESENT(radiusall)) THEN
      ALLOCATE(radiusall(Nspec,Nrmax0))
      radiusall(:,:) = radiusall0(:,:)
    END IF
    IF (PRESENT(tempall)) THEN
      ALLOCATE(tempall(Nspec,NTmax0))
      tempall(:,:) = tempall0(:,:)
    END IF
    IF (PRESENT(Nw)) Nw = Nw0
    IF (PRESENT(wave)) THEN 
      ALLOCATE(wave(Nw0))
      wave(:) = wave0(:)
    END IF
    IF (PRESENT(nu)) THEN 
      ALLOCATE(nu(Nw0))
      nu(:) = MKS%clight / MKS%micron / wave0(:)
    END IF
    IF (ALLOCATED(wave0)) DEALLOCATE (wave0)
    IF (ALLOCATED(waveall0)) DEALLOCATE(waveall0)
    IF (ALLOCATED(radiusall0)) DEALLOCATE(radiusall0)
    IF (ALLOCATED(tempall0)) DEALLOCATE(tempall0)    

    !-----------------------------------------------------------------------

  END SUBROUTINE read_optics_gen

  !-------------------------------------------------------------------------

  SUBROUTINE read_optics_scl (lab,Qabs,Qsca,gsca,Qav, & 
                               wrange,wave,nu,radius,temp,Nw,Nr,NT)
    USE utilities, ONLY: DP
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: lab
    INTEGER, INTENT(OUT), OPTIONAL :: Nw, Nr, NT
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Qabs, Qsca
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Qav, gsca
    REAL(DP), DIMENSION(2), INTENT(IN), OPTIONAL :: wrange
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: wave, nu
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: radius, temp
    INTEGER :: Nr0, Nw0, NT0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: radius0, temp0
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Qabs0, Qsca0, Qav0, gsca0
    CALL READ_OPTICS_GEN([lab],NWMAX=Nw0,NRMAX=Nr0,NTMAX=NT0, & 
                           RADIUSALL=radius0,TEMPALL=temp0,QABSALL=Qabs0, &
                           QSCAALL=Qsca0,GSCAALL=gsca0,QAVALL=Qav0, & 
                           WRANGE=wrange,WAVE=wave,NU=nu)
    IF (PRESENT(Nr)) Nr = Nr0
    IF (PRESENT(Nw)) Nw = Nw0
    IF (PRESENT(NT)) NT = NT0
    IF (PRESENT(radius)) THEN
      ALLOCATE(radius(Nr0))
      radius(:) = radius0(1,:)
    END IF
    IF (PRESENT(temp)) THEN 
      ALLOCATE(temp(NT0))
      temp(:) = temp0(1,:)
    ENDIF
    IF (PRESENT(Qabs)) THEN 
      ALLOCATE(Qabs(Nr0,Nw0))
      Qabs(:,:) = Qabs0(1,:,:)
    END IF
    IF (PRESENT(Qsca)) THEN 
      ALLOCATE(Qsca(Nr0,Nw0))
      Qsca(:,:) = Qsca0(1,:,:)
    END IF
    IF (PRESENT(Qav)) THEN 
      ALLOCATE(Qav(Nr0,NT0))
      Qav(:,:) = Qav0(1,:,:)
    END IF
    IF (PRESENT(gsca)) THEN 
      ALLOCATE(gsca(Nr0,Nw0))
      gsca(:,:) = gsca0(1,:,:)
    END IF
  END SUBROUTINE read_optics_scl


END MODULE grain_optics
