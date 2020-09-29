    Grid settings
    dogrid = ( SIZE(parvec(:)) > 1 ) ! par grid for Gibbs sampling
    ! IF (tied .AND. dogrid) THEN
    IF (dogrid) THEN
      IF (CONTset) THEN
        ALLOCATE(CONTgrid(Ncont))
        DO i=1,Ncont
          CONTgrid(i) = ( TRIMEQ(parname, "Mcont"//TRIMLR(PRING(i))) &
                          .OR. TRIMEQ(parname, "Tcont"//TRIMLR(PRING(i))) )

        END DO
      ELSE
        ALLOCATE(CONTgrid(1))
        CONTgrid(:) = .FALSE.

      END IF
      IF (BANDset) THEN
        ALLOCATE(BANDgrid(Nband))
        DO i=1,Nband
          BANDgrid(i) = ( TRIMEQ(parname, "Iband"//TRIMLR(PRING(i))) &
                          .OR. TRIMEQ(parname, "Cband"//TRIMLR(PRING(i))) &
                          .OR. TRIMEQ(parname, "WSband"//TRIMLR(PRING(i))) &
                          .OR. TRIMEQ(parname, "WLband"//TRIMLR(PRING(i))) )

        END DO
      ELSE
        ALLOCATE(BANDgrid(1))
        BANDgrid(:) = .FALSE.

      END IF
      STARgrid = ( TRIMEQ(compname, "Fstar") )
      PABSgrid = ( TRIMEQ(compname, "Av") )
      IF (LINEset) THEN
        ALLOCATE(LINEgrid(Nline))
        DO i=1,Nline
          LINEgrid(i) = ( TRIMEQ(parname, "Iline"//TRIMLR(PRING(i))) &
                          .OR. TRIMEQ(parname, "Cline"//TRIMLR(PRING(i))) &
                          .OR. TRIMEQ(parname, "Wline"//TRIMLR(PRING(i))) )

        END DO
      ELSE
        ALLOCATE(LINEgrid(1))
        LINEgrid(:) = .FALSE.

      END IF
    END IF

!!-------------------------------------------------------
  !!
  !!            [obsolet] Read optical properties
  !!
  !!-------------------------------------------------------
  SUBROUTINE make_Qabs(label, nQabs, waveall)
    !! wave [um], nu [Hz], Qova [m^-1], rho [kg/m^3]
    USE utilities, ONLY: DP
    USE constants, ONLY: MKS
    USE arrays, ONLY: closest
    USE interpolation, ONLY: interp_lin_sorted
    USE grain_optics, ONLY: rho_grain, read_optics
    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:), INTENT(IN)               :: label
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL         :: waveall
    INTEGER                                              :: i, Nspec, Nr, Nw
    INTEGER, DIMENSION(SIZE(label))                      :: indr ! radius index
    REAL(DP), PARAMETER                                  :: a0 = 1.E-2_DP ! grain radius
    REAL(DP), DIMENSION(:), ALLOCATABLE                  :: wave0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE                :: radius
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE              :: Q ! (Nspec, Nr, Nw)
    TYPE(Qabs_type), DIMENSION(SIZE(label)), INTENT(OUT) :: nQabs
    
    CALL read_optics(label, WAVE=wave0, RADIUSALL=radius, QABSALL=Q)

    Nspec = SIZE(label)
    Nw = SIZE(wave0)

    DO i=1,Nspec
      
      ALLOCATE(nQabs(i)%Qova(Nw))
      
      Nr = SIZE(radius(i,:))
      nQabs(i)%rho = rho_grain(label(i))
      nQabs(i)%wave = wave0
      nQabs(i)%nu = MKS%clight/MKS%micron / wave0
      indr(i) = closest(radius(i,:), a0)
      nQabs(i)%Qova = Q(i,indr(i),:)/radius(i,indr(i))
      ! PRINT*, 'Radius of ', label(i), ': ', radius(i,indr(i)), '->', indr(i)

      !! [Optional] Interpolate nQabs to input wave grid
      IF (PRESENT(waveall)) THEN
        nQabs(i)%wave = waveall
        nQabs(i)%nu = MKS%clight/MKS%micron / waveall
        nQabs(i)%Qova = interp_lin_sorted(nQabs(i)%Qova, wave0, waveall, &
                                          XLOG=.TRUE., YLOG=.TRUE., FORCE=.TRUE.)

      END IF
      
    END DO

    !! Free memory space
    DEALLOCATE(wave0, radius, Q)

  END SUBROUTINE make_Qabs
  
  !!-------------------------------------------------------
  !!
  !!       [obsolet] Create the parameter structure
  !!
  !!-------------------------------------------------------
  SUBROUTINE make_par(parvec, Nbb, Nline, Nband, NAv, Nstar, &
                      par, Npar, parinfo)

    USE utilities, ONLY: DP, PRING, TRIMLR
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nbb, Nline, Nband, NAv, Nstar
    REAL(DP), DIMENSION(:), INTENT(IN) :: parvec ! param vector
    INTEGER i, i0, Npar0
    TYPE(par_type) :: par0
    TYPE(par_type), INTENT(OUT), OPTIONAL :: par
    INTEGER, INTENT(OUT), OPTIONAL :: Npar
    TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: parinfo

    ! Npar0 = 2*Nbb + 3*Nline + 4*Nband + NAv + Nstar
    Npar0 = SIZE(parvec)
    IF (PRESENT(Npar)) Npar = Npar0
    IF (PRESENT(parinfo)) ALLOCATE(parinfo(Npar0))

    ALLOCATE(par0%Mcont(Nbb), par0%Tcont(Nbb)) ! Nbb*2
    ALLOCATE(par0%Iline(Nline), par0%Cline(Nline), par0%Wline(Nline)) ! Nline*3
    ALLOCATE(par0%Iband(Nband), par0%Cband(Nband), &
             par0%WSband(Nband), par0%WLband(Nband)) ! Nband*4
    ALLOCATE(par0%Av(NAv)) ! NAv*1
    ALLOCATE(par0%Fstar(Nstar)) ! Nstar*1

    !! 1) Continuum
    i0 = 0
    DO i=1,Nbb
      par0%Mcont(i) = parvec(i0+2*i-1) ! i0 + 2*(i-1) + 1
      par0%Tcont(i) = parvec(i0+2*i)

      IF (PRESENT(parinfo)) THEN
        parinfo(i0+2*i-1)%name = 'Mcont'//TRIMLR(PRING(i))
        parinfo(i0+2*i-1)%limits = [0._DP,1._DP]
        parinfo(i0+2*i-1)%limited = [.True.,.False.]
        parinfo(i0+2*i)%name = 'Tcont'//TRIMLR(PRING(i))
        parinfo(i0+2*i)%limits = [50._DP,500._DP]
        parinfo(i0+2*i)%limited = [.True.,.True.]
        
      END IF
      
    END DO

    !! 2) Lines
    i0 = i0 + 2*Nbb
    DO i=1,Nline
      par0%Iline(i) = parvec(i0+3*i-2)
      par0%Cline(i) = parvec(i0+3*i-1)
      par0%Wline(i) = parvec(i0+3*i)

      IF (PRESENT(parinfo)) THEN
        parinfo(i0+3*i-2)%name = 'Iline'//TRIMLR(PRING(i))
        parinfo(i0+3*i-2)%limits = [0._DP,1._DP]
        parinfo(i0+3*i-2)%limited = [.True.,.False.]
        parinfo(i0+3*i-1)%name = 'Cline'//TRIMLR(PRING(i))
        parinfo(i0+3*i-1)%limits = &
          [ par0%Cline(i) - par0%Wline(i)/par0%Cline(i)/2._DP, &
            par0%Cline(i) + par0%Wline(i)/par0%Cline(i)/2._DP ]
        parinfo(i0+3*i-1)%limited = [.True.,.True.]
        parinfo(i0+3*i)%name = 'Wline'//TRIMLR(PRING(i))
        parinfo(i0+3*i)%fixed = .True.
        
      END IF
      
    END DO

    !! 3) Bands
    i0 = i0 + 3*Nline
    DO i=1,Nband
      par0%Iband(i) = parvec(i0+4*i-3)
      par0%Cband(i) = parvec(i0+4*i-2)
      par0%WSband(i) = parvec(i0+4*i-1)
      par0%WLband(i) = parvec(i0+4*i)

      IF (PRESENT(parinfo)) THEN
        parinfo(i0+4*i-3)%name = 'Iband'//TRIMLR(PRING(i))
        parinfo(i0+4*i-3)%limits = [0._DP,1._DP]
        parinfo(i0+4*i-3)%limited = [.True.,.False.]
        parinfo(i0+4*i-2)%name = 'Cband'//TRIMLR(PRING(i))
        parinfo(i0+4*i-2)%fixed = .True.
        parinfo(i0+4*i-1)%name = 'WSband'//TRIMLR(PRING(i))
        parinfo(i0+4*i-1)%fixed = .True.
        parinfo(i0+4*i)%name = 'WLband'//TRIMLR(PRING(i))
        parinfo(i0+4*i)%fixed = .True.

      END IF
      
    END DO

    !! 4) Av
    i0 = i0 + 4*Nband
    DO i=1,NAv
      par0%Av(i) = parvec(i0+i)

      IF (PRESENT(parinfo)) THEN
        parinfo(i0+i)%name = 'Av'
        parinfo(i0+i)%fixed = .True.

      END IF

    END DO

    !! 5) Star
    i0 = i0 + NAv
    IF (Nstar .GT. 0) THEN
      DO i=1,Nstar
        par0%Fstar(i) = parvec(i0+1)

        IF (PRESENT(parinfo)) THEN
          parinfo(i0+i)%name = 'Fstar'
          parinfo(i0+i)%limits = [0._DP,1._DP]
          parinfo(i0+i)%limited = [.True.,.False.]

        END IF
        
      END DO
    END IF

    IF (PRESENT(par)) par = par0
  
  END SUBROUTINE make_par
