  !!-------------------------------------------------------
  !!
  !!               [obsolete] Generic intertface
  !!
  !!-------------------------------------------------------
  FUNCTION specModel_gen( wvl, parvec, parname, parinfo, indpar, parval, Qabs, extinct, &
                          verbose, debug )

    USE utilities, ONLY: DP, trimeq, trimlr, pring, &!verbatim, &
                         initiate_clock, time_type, timinfo!, ustd
    USE constants, ONLY: pi, MKS
    USE arrays, ONLY: reallocate
    USE statistical_physics, ONLY: blackbody
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN)           :: wvl, extinct
    REAL(DP), DIMENSION(:), INTENT(IN)           :: parvec, parval
    CHARACTER(*), INTENT(IN)                     :: parname
    TYPE(parinfo_type), DIMENSION(:), INTENT(IN) :: parinfo
    TYPE(indpar_type), INTENT(IN)                :: indpar
    TYPE(Qabs_type), DIMENSION(:), INTENT(IN)    :: Qabs
    LOGICAL, INTENT(IN), OPTIONAL                :: verbose, debug
    
    INTEGER :: Nw, Ngrid, Ncont, Nline, Nband, Npabs, Nstar
    INTEGER :: i, igrid, iw, j!, itied
    REAL(DP)                                     :: Tstar
    REAL(DP), DIMENSION(SIZE(wvl))               :: nu
    LOGICAL :: gridlnMovd2, gridlnT, gridlnIline, gridCline, gridWline
    LOGICAL :: gridlnIband, gridCband, gridWSband, gridWLband, gridlnAv, gridlnFstar
    REAL(DP), DIMENSION(SIZE(parvec),SIZE(wvl)) :: FnuCONT0, FnuBAND0, FnuSTAR0
    REAL(DP), DIMENSION(SIZE(parvec),SIZE(wvl)) :: Pabs0, FnuLINE0
    REAL(DP) :: Const
    REAL(DP), DIMENSION(SIZE(wvl)) :: Const1D
    REAL(DP), DIMENSION(SIZE(parvec),SIZE(wvl)) :: specModel_gen
    
    LOGICAL :: printimer, bls
    TYPE(time_type) :: timestr

    CALL INITIATE_CLOCK(timestr)
    
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
    
    !! Verbose
    printimer = .FALSE. ! Print timer
    IF (PRESENT(verbose)) printimer = verbose

    !! Debug
    bls = .FALSE.
    IF (PRESENT(debug)) bls = debug
    
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
        DO j=1,Ncont
          IF (j.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(j))) * MKS%Msun/MKS%pc**2 * &
                         modifBB(wvl(:), EXP(parval(indpar%lnT(j))), Qabs(j))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = FnuCONT0(:,iw) + Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO j=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(j))) * &
                       lorentzBand(wvl(:), parval(indpar%Cband(j)), &
                         parval(indpar%WSband(j)), parval(indpar%WLband(j)))
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR
        Const1D(:) = 0._DP
        DO j=1,Nstar
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnFstar(j))) * &
                       pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

        END DO
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO j=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(j))) * &
                       gaussLine(wvl(:), parval(indpar%Cline(j)), &
                         parval(indpar%Wline(j)))
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS (May need to modify if consider different geometries)
        Const1D(:) = 0._DP
        DO j=1,Npabs
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnAv(j)))

        END DO
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP( -Const1D(iw)/1.086_DP * extinct(iw) )
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
        DO j=1,Ncont
          IF (j.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(j))) * MKS%Msun/MKS%pc**2 * &
                         modifBB(wvl(:), EXP(parval(indpar%lnT(j))), Qabs(j))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = FnuCONT0(:,iw) + Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO j=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(j))) * &
                       lorentzBand(wvl(:), parval(indpar%Cband(j)), &
                         parval(indpar%WSband(j)), parval(indpar%WLband(j)))
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR
        Const1D(:) = 0._DP
        DO j=1,Nstar
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnFstar(j))) * &
                       pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

        END DO
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO j=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(j))) * &
                       gaussLine(wvl(:), parval(indpar%Cline(j)), &
                         parval(indpar%Wline(j)))
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        Const1D(:) = 0._DP
        DO j=1,Npabs
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnAv(j)))

        END DO
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP( -Const1D(iw)/1.086_DP * extinct(iw) )
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      END IF

      !! BLS (Bug Localization System)
      IF (bls) THEN
        IF (ANY(isNaN(specModel_gen))) THEN
          IF (gridlnMovd2) THEN
            PRINT*, "lnMovd2"//TRIMLR(PRING(i))
          ELSE IF (gridlnT) THEN
            PRINT*, "lnT"//TRIMLR(PRING(i))
          END IF
          
          IF (ANY(isNaN(FnuCONT0))) THEN
            PRINT*, 'FnuCONT0'
          END IF
          IF (ANY(isNaN(FnuBAND0))) THEN
            PRINT*, 'FnuBAND0'
          END IF
          IF (ANY(isNaN(FnuSTAR0))) THEN
            PRINT*, 'FnuSTAR0'
          END IF
          IF (ANY(isNaN(PABS0))) THEN
            PRINT*, 'PABS0'
          END IF
          IF (ANY(isNaN(FnuLINE0))) THEN
            PRINT*, 'FnuLINE0'
          END IF
          STOP
        END IF
      END IF
      
    END DO loopCONT

    IF (printimer) &
      PRINT*, "[specModel] CAL CONT IN "//TRIMLR(TIMINFO(timestr))//"."
    
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
        DO j=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(j))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(j))), Qabs(j))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND *lnIband
        Const1D(:) = lorentzBand(wvl(:), parval(indpar%Cband(i)), &
                       parval(indpar%WSband(i)), parval(indpar%WLband(i)))
        FORALL (igrid=1:Ngrid) &
          FnuBAND0(igrid,:) = EXP(parvec(igrid)) * Const1D(:)
! if (any(isnan(Const1D))) print*, "*lnIband", &
  ! new_line('a'), parval(indpar%cband(i)), parval(indpar%wsband(i)), parval(indpar%wlband(i))

        !! BAND
        Const1D(:) = 0._DP
        DO j=1,Nband
          IF (j.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(j))) * &
                         lorentzBand(wvl(:), parval(indpar%Cband(j)), &
                           parval(indpar%WSband(j)), parval(indpar%WLband(j)))
! if (j.NE.i .AND. any(isnan(lorentzBand(wvl(:), parval(indpar%Cband(j)), &
!   parval(indpar%WSband(j)), parval(indpar%WLband(j)))))) &
!   print*, "lnIband", j, new_line('a'), &
!   new_line('a'), parval(indpar%cband(j)), parval(indpar%wsband(j)), parval(indpar%wlband(j))

        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = FnuBAND0(:,iw) + Const1D(iw)
        !! STAR
        Const1D(:) = 0._DP
        DO j=1,Nstar
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnFstar(j))) * &
                       pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

        END DO
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO j=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(j))) * &
                       gaussLine(wvl(:), parval(indpar%Cline(j)), &
                         parval(indpar%Wline(j)))
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        Const1D(:) = 0._DP
        DO j=1,Npabs
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnAv(j)))

        END DO
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP( -Const1D(iw)/1.086_DP * extinct(iw) )
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      ELSE IF (gridCband) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO j=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(j))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(j))), Qabs(j))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND *Cband
        Const = EXP(parval(indpar%lnIband(i)))
        FORALL (igrid=1:Ngrid) &
          FnuBAND0(igrid,:) = Const * lorentzBand(wvl(:), parvec(igrid), &
                                        parval(indpar%WSband(i)), parval(indpar%WLband(i)))
! if (any(isnan(fnuband0))) print*, "*Cband"//trimlr(pring(i)), &
  ! new_line('a'), parval(indpar%cband(i)), parval(indpar%wsband(i)), parval(indpar%wlband(i))

        !! BAND
        Const1D(:) = 0._DP
        DO j=1,Nband
          IF (j.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(j))) * &
                         lorentzBand(wvl(:), parval(indpar%Cband(j)), &
                           parval(indpar%WSband(j)), parval(indpar%WLband(j)))
! if (j.NE.i .AND. any(isnan(lorentzBand(wvl(:), parval(indpar%Cband(j)), &
!   parval(indpar%WSband(j)), parval(indpar%WLband(j)))))) &
!   print*, "Cband", j, new_line('a'), &
!   new_line('a'), parval(indpar%cband(j)), parval(indpar%wsband(j)), parval(indpar%wlband(j))

        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = FnuBAND0(:,iw) + Const1D(iw)
        !! STAR
        Const1D(:) = 0._DP
        DO j=1,Nstar
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnFstar(j))) * &
                       pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

        END DO
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO j=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(j))) * &
                       gaussLine(wvl(:), parval(indpar%Cline(j)), &
                         parval(indpar%Wline(j)))
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        Const1D(:) = 0._DP
        DO j=1,Npabs
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnAv(j)))

        END DO
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP( -Const1D(iw)/1.086_DP * extinct(iw) )
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      ELSE IF (gridWSband) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO j=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(j))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(j))), Qabs(j))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND *WSband
        Const = EXP(parval(indpar%lnIband(i)))
        FORALL (igrid=1:Ngrid) &
          FnuBAND0(igrid,:) = Const * lorentzBand(wvl(:), parval(indpar%Cband(i)), &
                                        parvec(igrid), parval(indpar%WLband(i)))
! if (any(isnan(fnuband0))) print*, "*WSband", &
!   new_line('a'), parval(indpar%cband(i)), parval(indpar%wsband(i)), parval(indpar%wlband(i))

        !! BAND
        Const1D(:) = 0._DP
        DO j=1,Nband
          IF (j.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(j))) * &
                         lorentzBand(wvl(:), parval(indpar%Cband(j)), &
                           parval(indpar%WSband(j)), parval(indpar%WLband(j)))
! if (j.NE.i .AND. any(isnan(lorentzBand(wvl(:), parval(indpar%Cband(j)), &
!   parval(indpar%WSband(j)), parval(indpar%WLband(j)))))) &
!   print*, "WSband", j, new_line('a'), &
!   new_line('a'), parval(indpar%cband(j)), parval(indpar%wsband(j)), parval(indpar%wlband(j))
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = FnuBAND0(:,iw) + Const1D(iw)
        !! STAR
        Const1D(:) = 0._DP
        DO j=1,Nstar
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnFstar(j))) * &
                       pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

        END DO
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO j=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(j))) * &
                       gaussLine(wvl(:), parval(indpar%Cline(j)), &
                         parval(indpar%Wline(j)))
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        Const1D(:) = 0._DP
        DO j=1,Npabs
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnAv(j)))

        END DO
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP( -Const1D(iw)/1.086_DP * extinct(iw) )
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)        

      ELSE IF (gridWLband) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO j=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(j))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(j))), Qabs(j))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND *WLband
        Const = EXP(parval(indpar%lnIband(i)))
        FORALL (igrid=1:Ngrid) &
          FnuBAND0(igrid,:) = Const * lorentzBand(wvl(:), parval(indpar%Cband(i)), &
                                        parval(indpar%WSband(i)), parvec(igrid))
! if (any(isnan(fnuband0))) print*, "*WLband", &
!   new_line('a'), parval(indpar%cband(i)), parval(indpar%wsband(i)), parval(indpar%wlband(i))

        !! BAND
        Const1D(:) = 0._DP
        DO j=1,Nband
          IF (j.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(j))) * &
                         lorentzBand(wvl(:), parval(indpar%Cband(j)), &
                           parval(indpar%WSband(j)), parval(indpar%WLband(j)))
! if (j.NE.i .AND. any(isnan(lorentzBand(wvl(:), parval(indpar%Cband(j)), &
!   parval(indpar%WSband(j)), parval(indpar%WLband(j)))))) &
!   print*, "WLband", j, new_line('a'), &
!   new_line('a'), parval(indpar%cband(j)), parval(indpar%wsband(j)), parval(indpar%wlband(j))

        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = FnuBAND0(:,iw) + Const1D(iw)
        !! STAR
        Const1D(:) = 0._DP
        DO j=1,Nstar
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnFstar(j))) * &
                       pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

        END DO
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO j=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(j))) * &
                       gaussLine(wvl(:), parval(indpar%Cline(j)), &
                         parval(indpar%Wline(j)))
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        Const1D(:) = 0._DP
        DO j=1,Npabs
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnAv(j)))

        END DO
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP( -Const1D(iw)/1.086_DP * extinct(iw) )
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      END IF

      !! BLS (Bug Localization System)
      IF (bls) THEN
        IF (ANY(isNaN(specModel_gen))) THEN
          IF (gridlnIband) THEN
            PRINT*, "lnIband"//TRIMLR(PRING(i))
          ELSE IF (gridCband) THEN
            PRINT*, "Cband"//TRIMLR(PRING(i))
          ELSE IF (gridWSband) THEN
            PRINT*, "WSband"//TRIMLR(PRING(i))
          ELSE IF (gridWLband) THEN
            PRINT*, "WLband"//TRIMLR(PRING(i))
          END IF
          
          IF (ANY(isNaN(FnuCONT0))) THEN
            PRINT*, 'FnuCONT0'
          END IF
          IF (ANY(isNaN(FnuBAND0))) THEN
            PRINT*, 'FnuBAND0'
          END IF
          IF (ANY(isNaN(FnuSTAR0))) THEN
            PRINT*, 'FnuSTAR0'
          END IF
          IF (ANY(isNaN(PABS0))) THEN
            PRINT*, 'PABS0'
          END IF
          IF (ANY(isNaN(FnuLINE0))) THEN
            PRINT*, 'FnuLINE0'
          END IF
          STOP
        END IF
      END IF
      
    END DO loopBAND

    IF (printimer) &
      PRINT*, "[specModel] CAL BAND IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! 3. Stellar Continuum
    !!----------------------
    loopSTAR: DO i=1,Nstar
      gridlnFstar = TRIMEQ(parname, 'lnFstar'//TRIMLR(PRING(i)))
      
      IF (gridlnFstar) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO j=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(j))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(j))), Qabs(j))
  
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO j=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(j))) * &
                       lorentzBand(wvl(:), parval(indpar%Cband(j)), &
                         parval(indpar%WSband(j)), parval(indpar%WLband(j)))
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR *lnFstar
        FORALL (igrid=1:Ngrid) &
          FnuSTAR0(igrid,:) = EXP(parvec(igrid)) * &
                              pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4
        !! STAR
        Const1D(:) = 0._DP
        DO j=1,Nstar
          IF (j.NE.i) &
            Const1D(:) = Const1D(:) + EXP(parval(indpar%lnFstar(j))) * &
                         pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

        END DO
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = FnuSTAR0(:,iw) + Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO j=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(j))) * &
                       gaussLine(wvl(:), parval(indpar%Cline(j)), &
                         parval(indpar%Wline(j)))
  
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS
        Const1D(:) = 0._DP
        DO j=1,Npabs
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnAv(j)))

        END DO
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP( -Const1D(iw)/1.086_DP * extinct(iw) )
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      END IF

      !! BLS (Bug Localization System)
      IF (bls) THEN
        IF (ANY(isNaN(specModel_gen))) THEN
          IF (gridlnFstar) PRINT*, "lnFstar"//TRIMLR(PRING(j))

          IF (ANY(isNaN(FnuCONT0))) THEN
            PRINT*, 'FnuCONT0'
          END IF
          IF (ANY(isNaN(FnuBAND0))) THEN
            PRINT*, 'FnuBAND0'
          END IF
          IF (ANY(isNaN(FnuSTAR0))) THEN
            PRINT*, 'FnuSTAR0'
          END IF
          IF (ANY(isNaN(PABS0))) THEN
            PRINT*, 'PABS0'
          END IF
          IF (ANY(isNaN(FnuLINE0))) THEN
            PRINT*, 'FnuLINE0'
          END IF
          STOP
        END IF
      END IF
      
    END DO loopSTAR
    
    IF (printimer) &
      PRINT*, "[specModel] CAL STAR IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! 4. Screen extinction
    !!----------------------
    loopPABS: DO i=1,Npabs
      gridlnAv = TRIMEQ(parname, 'lnAv'//TRIMLR(PRING(i)))
      
      IF (gridlnAv) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO j=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(j))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(j))), Qabs(j))
      
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO j=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(j))) * &
                       lorentzBand(wvl(:), parval(indpar%Cband(j)), &
                         parval(indpar%WSband(j)), parval(indpar%WLband(j)))
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR
        Const1D(:) = 0._DP
        DO j=1,Nstar
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnFstar(j))) * &
                       pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

        END DO
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE
        Const1D(:) = 0._DP
        DO j=1,Nline
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(j))) * &
                       gaussLine(wvl(:), parval(indpar%Cline(j)), &
                         parval(indpar%Wline(j)))
      
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = Const1D(iw)
        !! PABS *lnAv
        FORALL (igrid=1:Ngrid) &
          Pabs0(igrid,:) = EXP( -EXP(parvec(igrid))/1.086_DP * extinct(:) )
        !! PABS
        Const1D(:) = 0._DP
        DO j=1,Npabs
          IF (j.NE.i) Const1D(:) = Const1D(:) + EXP(parval(indpar%lnAv(j)))

        END DO
        FORALL (iw=1:Nw) Pabs0(:,iw) = Pabs0(:,iw) * EXP( -Const1D(:)/1.086_DP * extinct(iw) )
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      END IF

      !! BLS (Bug Localization System)
      IF (bls) THEN
        IF (ANY(isNaN(specModel_gen))) THEN
          IF (gridlnAv) PRINT*, "lnAv"//TRIMLR(PRING(i))
          
          IF (ANY(isNaN(FnuCONT0))) THEN
            PRINT*, 'FnuCONT0'
          END IF
          IF (ANY(isNaN(FnuBAND0))) THEN
            PRINT*, 'FnuBAND0'
          END IF
          IF (ANY(isNaN(FnuSTAR0))) THEN
            PRINT*, 'FnuSTAR0'
          END IF
          IF (ANY(isNaN(PABS0))) THEN
            PRINT*, 'PABS0'
          END IF
          IF (ANY(isNaN(FnuLINE0))) THEN
            PRINT*, 'FnuLINE0'
          END IF
          STOP
        END IF
      END IF
      
    END DO loopPABS

    IF (printimer) &
      PRINT*, "[specModel] CAL PABS IN "//TRIMLR(TIMINFO(timestr))//"."
    
    !! 5. Lines
    !!----------
    loopLINE: DO i=1,Nline
      gridlnIline = TRIMEQ(parname, 'lnIline'//TRIMLR(PRING(i)))
      gridCline = TRIMEQ(parname, 'Cline'//TRIMLR(PRING(i)))
      gridWline = TRIMEQ(parname, 'Wline'//TRIMLR(PRING(i)))
      
      IF (gridlnIline) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO j=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(j))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(j))), Qabs(j))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO j=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(j))) * &
                       lorentzBand(wvl(:), parval(indpar%Cband(j)), &
                         parval(indpar%WSband(j)), parval(indpar%WLband(j)))
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR
        Const1D(:) = 0._DP
        DO j=1,Nstar
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnFstar(j))) * &
                       pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

        END DO
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE *lnIline
        Const1D(:) =  gaussLine(wvl(:), parval(indpar%Cline(i)), &
                        parval(indpar%Wline(i)))
! if (isnan(EXP(parval(indpar%lnIline(j))))) print*, "exp*"
! if (any(isnan(gaussLine(wvl(:), parval(indpar%Cline(j)), parval(indpar%Wline(j)))))) print*, "gauss*"
        FORALL (igrid=1:Ngrid) &
          FnuLINE0(igrid,:) = EXP(parvec(igrid)) * Const1D(:)
        !! LINE
        Const1D(:) = 0._DP
        DO j=1,Nline
          IF (j.NE.i) &
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(j))) * &
                       gaussLine(wvl(:), parval(indpar%Cline(j)), &
                         parval(indpar%Wline(j)))
! if (isnan(EXP(parval(indpar%lnIline(j))))) print*, "exp"
! if (any(isnan(gaussLine(wvl(:), parval(indpar%Cline(j)), parval(indpar%Wline(j)))))) print*, "gauss"        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = FnuLINE0(:,iw) + Const1D(iw)
        !! PABS
        Const1D(:) = 0._DP
        DO j=1,Npabs
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnAv(j)))

        END DO
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP( -Const1D(iw)/1.086_DP * extinct(iw) )
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      ELSE IF (gridCline) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO j=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(j))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(j))), Qabs(j))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO j=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(j))) * &
                       lorentzBand(wvl(:), parval(indpar%Cband(j)), &
                         parval(indpar%WSband(j)), parval(indpar%WLband(j)))
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR
        Const1D(:) = 0._DP
        DO j=1,Nstar
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnFstar(j))) * &
                       pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

        END DO
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE *Cline
        Const = EXP(parval(indpar%lnIline(i)))
        FORALL (igrid=1:Ngrid) &
          FnuLINE0(igrid,:) = Const * gaussLine(wvl(:), parvec(igrid), &
                                       parval(indpar%Wline(i)))
        !! LINE
        Const1D(:) = 0._DP
        DO j=1,Nline
          IF (j.NE.i) &
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(j))) * &
                       gaussLine(wvl(:), parval(indpar%Cline(j)), &
                         parval(indpar%Wline(j)))
! if (isnan(EXP(parval(indpar%lnIline(j))))) print*, "exp"
! if (any(isnan(gaussLine(wvl(:), parval(indpar%Cline(j)), parval(indpar%Wline(j)))))) print*, "gauss"
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = FnuLINE0(:,iw) + Const1D(iw)
        !! PABS
        Const1D(:) = 0._DP
        DO j=1,Npabs
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnAv(j)))

        END DO
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP( -Const1D(iw)/1.086_DP * extinct(iw) )
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      ELSE IF (gridWline) THEN
        !! CONT
        Const1D(:) = 0._DP
        DO j=1,Ncont
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnMovd2(j))) * MKS%Msun/MKS%pc**2 * &
                       modifBB(wvl(:), EXP(parval(indpar%lnT(j))), Qabs(j))
        
        END DO
        FORALL (iw=1:Nw) FnuCONT0(:,iw) = Const1D(iw)
        !! BAND
        Const1D(:) = 0._DP
        DO j=1,Nband
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIband(j))) * &
                       lorentzBand(wvl(:), parval(indpar%Cband(j)), &
                         parval(indpar%WSband(j)), parval(indpar%WLband(j)))
          
        END DO
        FORALL (iw=1:Nw) FnuBAND0(:,iw) = Const1D(iw)
        !! STAR
        Const1D(:) = 0._DP
        DO j=1,Nstar
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnFstar(j))) * &
                       pi * BLACKBODY(nu(:), Tstar) /MKS%stefan/Tstar**4

        END DO
        FORALL (iw=1:Nw) FnuSTAR0(:,iw) = Const1D(iw)
        !! LINE *Wline
        Const = EXP(parval(indpar%lnIline(i)))
        FORALL (igrid=1:Ngrid) &
          FnuLINE0(igrid,:) = Const * gaussLine(wvl(:), parvec(igrid), &
                                       parval(indpar%Wline(i)))
        !! LINE
        Const1D(:) = 0._DP
        DO j=1,Nline
          IF (j.NE.i) &
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnIline(j))) * &
                       gaussLine(wvl(:), parval(indpar%Cline(j)), &
                         parval(indpar%Wline(j)))
        
        END DO
        FORALL (iw=1:Nw) FnuLINE0(:,iw) = FnuLINE0(:,iw) + Const1D(iw)
        !! PABS
        Const1D(:) = 0._DP
        DO j=1,Npabs
          Const1D(:) = Const1D(:) + EXP(parval(indpar%lnAv(j)))

        END DO
        FORALL (iw=1:Nw) Pabs0(:,iw) = EXP( -Const1D(iw)/1.086_DP * extinct(iw) )
        !! Total model
        specModel_gen(:,:) = (FnuCONT0(:,:) + FnuBAND0(:,:) + FnuSTAR0(:,:)) * &
                             Pabs0(:,:) + FnuLINE0(:,:)

      END IF

      !! BLS (Bug Localization System)
      IF (bls) THEN
        IF (ANY(isNaN(specModel_gen))) THEN
          IF (gridlnIline) THEN
            PRINT*, "lnIline"//TRIMLR(PRING(i))
          ELSE IF (gridCline) THEN
            PRINT*, "Cline"//TRIMLR(PRING(i))
          ELSE IF (gridWline) THEN
            PRINT*, "Wline"//TRIMLR(PRING(i))
          END IF
          
          IF (ANY(isNaN(FnuCONT0))) THEN
            PRINT*, 'FnuCONT0'
          END IF
          IF (ANY(isNaN(FnuBAND0))) THEN
            PRINT*, 'FnuBAND0'
          END IF
          IF (ANY(isNaN(FnuSTAR0))) THEN
            PRINT*, 'FnuSTAR0'
          END IF
          IF (ANY(isNaN(PABS0))) THEN
            PRINT*, 'PABS0'
          END IF
          IF (ANY(isNaN(FnuLINE0))) THEN
            PRINT*, 'FnuLINE0'
          END IF
          STOP
        END IF
      END IF
      
    END DO loopLINE

    IF (printimer) &
      PRINT*, "[specModel] CAL LINE IN "//TRIMLR(TIMINFO(timestr))//"."
    
  END FUNCTION specModel_gen

  !!-------------------------------------------------------
  !!
  !!                      Grid settings
  !!
  !!-------------------------------------------------------
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
