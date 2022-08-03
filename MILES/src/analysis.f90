!******************************************************************************
!*
!*                           MILES POST-PROCESSING
!*
!******************************************************************************


PROGRAM analysis

  USE core, ONLY: read_master, read_analysis, specModel
  USE utilities, ONLY: DP, pring, trimLR, &
                       banner_program, ustd, isNaN, warning, &
                       trimeq, time_type, timinfo, initiate_clock
  USE arrays, ONLY: iwhere, closest, reallocate, sort
  USE statistics, ONLY: mean, sigma, median, autocorrel, intautocorrtime
  USE grain_optics, ONLY: lendustQ
  USE inout, ONLY: read_hdf5, write_hdf5, check_hdf5, h5ext, ascext, &
                   lenpar, lenpath
  USE hb, ONLY: Nx, Ny, NwOBS, FnuOBS, Nparhyp, Ncorrhyp, &
                mask, ind, parinfo, parhypinfo, Qabs, extinct
  IMPLICIT NONE

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: textwid = 60
  CHARACTER(*), PARAMETER :: nampar = "Parameter values"
  CHARACTER(*), PARAMETER :: nammu = "Mean of hyperdistribution"
  CHARACTER(*), PARAMETER :: namsig = "Sigma of hyperdistribution"
  CHARACTER(*), PARAMETER :: namcorr = "Correlation of hyperdistribution"
  LOGICAL, PARAMETER :: compress = .True., forcetint = .True.
  ! LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Input variables
  INTEGER :: i, j, x, y, iw, ipar, ihyp, icorr
  INTEGER :: Npar, Nmcmc, lastind!, NiniMC
  ! INTEGER :: Ncont, Nband, Nline, Nextc, Nstar, Nextra
  INTEGER, DIMENSION(:), ALLOCATABLE :: indresume
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: maskint ! convert mask=0 to mask=T
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wOBS
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: parerr, par
  CHARACTER(lenpath) :: dirIN, dirOUT, filOBS, filMCMC, fiLOG, fileCHI2, fileHB
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lenpath), DIMENSION(:), ALLOCATABLE :: path1d
  !! Analysis variables
  INTEGER :: t_burnin, t_end, Nmcmc_eff
  INTEGER :: Neff_ref1, Neff_ref2, Nsou
  INTEGER, DIMENSION(:), ALLOCATABLE :: Neff, Neff2
  REAL(DP) :: t_int_ref1, t_int_ref2
  REAL(DP), DIMENSION(:), ALLOCATABLE :: xs
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meanmu, stdevmu, meansig, stdevsig
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meancorr, stdevcorr
  REAL(DP), DIMENSION(:), ALLOCATABLE :: autocorfun, autocorfun2
  REAL(DP), DIMENSION(:), ALLOCATABLE :: t_int, t_int2
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mumcmc, sigmcmc, corrmcmc
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: meanpar, stdevpar
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT, &
    FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FnuCONT_tab, &
    FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab, FnuMOD_mcmc, n_FnuMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parmcmc, qpar
  LOGICAL :: chi2OK, hbOK, ACF, verbose!, calib
  LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: mask4D

  !! Output variables
  ! INTEGER :: Ncorr
  INTEGER, DIMENSION(2) :: unitlog
  TYPE(time_type) :: timestr

  !! Preliminaires
  !!---------------
  CALL READ_HDF5(STRARR1D=path1d, FILE='../out/set_input'//h5ext, NAME='input dir')
  dirIN = path1d(1)
  filOBS = TRIMLR(dirIN)//'observation_MIR'//h5ext

  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='wavelength (microns)', N1=NwOBS)
  CALL READ_MASTER(WAVALL=wOBS(:), DIRIN=dirIN, DIROUT=dirOUT, VERBOSE=verbose, &
                   NMCMC=Nmcmc, QABS=Qabs, EXTINCT=extinct, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit, &
                   PARHYPINFO=parhypinfo, NPARHYP=Nparhyp, NCORRHYP=Ncorrhyp)
  
  fiLOG = TRIMLR(dirOUT)//'log_analysis'//ascext

  CALL INITIATE_CLOCK(timestr)
  OPEN (ulog,FILE=fiLOG,STATUS="REPLACE",ACTION="WRITE")
  unitlog(:) = [ ulog, ustd ]
  
  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, &
                 NAME='FnuOBS ('//TRIMLR(spec_unit)//')', N1=Nx, N2=Ny)
  CALL READ_HDF5(INTARR3D=maskint, FILE=filOBS, NAME='NaN mask')
  ALLOCATE(mask(Nx,Ny,NwOBS))
  mask(:,:,:) = ( maskint(:,:,:) == 0 )
  Nsou = COUNT(ANY(mask(:,:,:),DIM=3))

  !!------------------------------------------------------------------------
  !!                      LE MIROIR (Chi2 Fit) Analysis
  !!------------------------------------------------------------------------
  
  fileCHI2 = TRIMLR(dirOUT)//'fit_chi2'//h5ext
  chi2OK = CHECK_HDF5(FILE=fileCHI2)
  IF (chi2OK) THEN
    PRINT*, 'LE MIROIR: Chi2 fit analysis'

    CALL READ_HDF5(DBLARR3D=par, FILE=fileCHI2, NAME='Best fitted parameter value')
    CALL READ_HDF5(DBLARR3D=parerr, FILE=fileCHI2, NAME='Best fitted parameter error')

    CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME='Parameter label', &
                    FILE=fileCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
    CALL WRITE_HDF5(DBLARR3D=par(:,:,:), NAME='Best fitted parameter value', &
                    FILE=fileCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR3D=parerr(:,:,:), NAME='Best fitted parameter error', &
                    FILE=fileCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

    !! Calculate model
    !!-----------------
    ALLOCATE(FnuMOD(Nx,Ny,NwOBS))
    
    FnuMOD(:,:,:) = specModel( wOBS(:), INDPAR=ind, PARVAL=par(:,:,:), &
                               QABS=Qabs(:), EXTINCT=extinct(:,:), &
                               FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
                               PABS=Pabs, FNULINE=FnuLINE, &
                               FNUCONT_TAB=FnuCONT_tab, FNUBAND_TAB=FnuBAND_tab, &
                               FNUSTAR_TAB=FnuSTAR_tab, PABS_TAB=Pabs_tab, &
                               FNULINE_TAB=FnuLINE_tab )
    
    CALL WRITE_HDF5(DBLARR4D=FnuCONT_tab, NAME='FnuCONT ('//TRIMLR(spec_unit)//')', &
                    FILE=fileCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR4D=FnuLINE_tab, NAME='FnuLINE ('//TRIMLR(spec_unit)//')', &
                    FILE=fileCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR4D=FnuBAND_tab, NAME='FnuBAND ('//TRIMLR(spec_unit)//')', &
                    FILE=fileCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR4D=FnuSTAR_tab, NAME='FnuSTAR ('//TRIMLR(spec_unit)//')', &
                    FILE=fileCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR4D=Pabs_tab, NAME='PABS', &
                    FILE=fileCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR3D=FnuMOD, NAME='FnuMOD ('//TRIMLR(spec_unit)//')', &
                    FILE=fileCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    
    DO i=1,MERGE(2,1,verbose)
      WRITE(unitlog(i),*)
      WRITE(unitlog(i),*) 'LE MIROIR: fit_chi2 analysis [done]'
      WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
    END DO
    
    !! Free memory space
    DEALLOCATE(FnuMOD)
    
  END IF

  !!------------------------------------------------------------------------
  !!                        HISTOIRE (HB Fit) Analysis
  !!------------------------------------------------------------------------

  DO j=1,2
    IF (j==1) THEN
      fileHB = TRIMLR(dirOUT)//'fit_bb'//h5ext
      filMCMC = TRIMLR(dirOUT)//'parlog_fit_bb'//h5ext
    ELSE
      fileHB = TRIMLR(dirOUT)//'fit_hb'//h5ext
      filMCMC = TRIMLR(dirOUT)//'parlog_fit_hb'//h5ext
    END IF
    hbOK = CHECK_HDF5(FILE=filMCMC)
    IF (hbOK) THEN
      IF (j==1) THEN
        PRINT*, 'HISTOIRE: BB fit analysis'
      ELSE
        PRINT*, 'HISTOIRE: HB fit analysis'
      END IF

      CALL READ_ANALYSIS(DIRIN=dirIN, DIROUT=dirOUT, &
                         ACF=ACF, VERBOSE=verbose, &
                         T_BURNIN=t_burnin, T_END=t_end)
  
      !! Compute the statistics
      !!------------------------
      lastind = 0
      CALL READ_HDF5(INTARR1D=indresume, FILE=filMCMC, NAME='Last index')
      IF (lastind < indresume(1)) lastind = indresume(1)
      IF (lastind < Nmcmc) &
        CALL WARNING("HISTOIRE", "Unfinished fits! Nmcmc shrank to "//TRIMLR(PRING(lastind)))
      IF (t_end <= 0 .OR. t_end > Nmcmc) t_end = Nmcmc
      IF (t_end > lastind) t_end = lastind
      IF (t_burnin <= 0 .OR. t_burnin >= t_end) t_burnin = INT(t_end * 0.1_DP)

      Nmcmc_eff = t_end - t_burnin + 1
      
      ALLOCATE(meanpar(Nx,Ny,Npar), stdevpar(Nx,Ny,Npar))
      ALLOCATE(qpar(Nx,Ny,Npar,3), xs(Nmcmc_eff))
      
      CALL READ_HDF5(DBLARR4D=parmcmc, FILE=filMCMC, &
                     NAME=nampar, IND4=[1,t_end])
      meanpar(:,:,:) = MEAN(parmcmc(:,:,:,t_burnin:t_end), DIM=4)
      stdevpar(:,:,:) = SIGMA(parmcmc(:,:,:,t_burnin:t_end), DIM=4)
      qpar(:,:,:,2) = MEDIAN(parmcmc(:,:,:,t_burnin:t_end), DIM=4)
      !! Calculate 1st and 5th 6-quantiles (sextiles)
      DO x=1,Nx
        DO y=1,Ny
          DO ipar=1,Npar
            xs(:) = SORT(parmcmc(x,y,ipar,t_burnin:t_end))
            qpar(x,y,ipar,1) = xs(NINT(1._DP/6*Nmcmc_eff)) ! round
            qpar(x,y,ipar,3) = xs(NINT(5._DP/6*Nmcmc_eff))
          END DO
        END DO
      END DO

      CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME="Parameter label", &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
      CALL WRITE_HDF5(DBLARR3D=meanpar(:,:,:), NAME="Mean of parameter value", &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR3D=stdevpar(:,:,:), NAME="Sigma of parameter value", &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR4D=qpar(:,:,:,:), NAME="Quantiles of parameter value", &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

      !! Compute the average, standard deviations and correlations of each quantity
      !!----------------------------------------------------------------------------
      CALL READ_HDF5(DBLARR2D=mumcmc, FILE=filMCMC, &
                     NAME=nammu, IND2=[1,t_end])
      CALL READ_HDF5(DBLARR2D=sigmcmc, FILE=filMCMC, &
                     NAME=namsig, IND2=[1,t_end])
      CALL READ_HDF5(DBLARR2D=corrmcmc, FILE=filMCMC, &
                     NAME=namcorr, IND2=[1,t_end])
  
      ALLOCATE (meanmu(Nparhyp),stdevmu(Nparhyp),meansig(Nparhyp),stdevsig(Nparhyp))
      
      meanmu(:) = MEAN(mumcmc(:,t_burnin:t_end), DIM=2)
      stdevmu(:) = SIGMA(mumcmc(:,t_burnin:t_end), DIM=2)
      meansig(:) = MEAN(sigmcmc(:,t_burnin:t_end), DIM=2)
      stdevsig(:) = SIGMA(sigmcmc(:,t_burnin:t_end), DIM=2)
      
      ALLOCATE (meancorr(Ncorrhyp),stdevcorr(Ncorrhyp))
      
      meancorr(:) = MEAN(corrmcmc(:,t_burnin:t_end), DIM=2)
      stdevcorr(:) = SIGMA(corrmcmc(:,t_burnin:t_end), DIM=2)

      CALL WRITE_HDF5(STRARR1D=parhypinfo(:)%name, NAME='Hyperparameter label', &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR1D=meanmu(:),NAME="Mean of parameter average", &
                      FILE=fileHB,COMPRESS=compress,VERBOSE=debug,APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR1D=stdevmu(:),NAME="Sigma of parameter average", &
                      FILE=fileHB,COMPRESS=compress,VERBOSE=debug,APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR1D=meansig(:), &
                      NAME="Mean of parameter standard deviation", &
                      FILE=fileHB,COMPRESS=compress,VERBOSE=debug,APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR1D=stdevsig(:), &
                      NAME="Sigma of parameter standard deviation", &
                      FILE=fileHB,COMPRESS=compress,VERBOSE=debug,APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR1D=meancorr(:),NAME="Mean of parameter correlation", &
                      FILE=fileHB,COMPRESS=compress,VERBOSE=debug,APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR1D=stdevcorr(:),NAME="Sigma of parameter correlation", &
                      FILE=fileHB,COMPRESS=compress,VERBOSE=debug,APPEND=.TRUE.)

      !! ACFs
      !!------
      doacf: IF (ACF) THEN

        DO i=1,MERGE(2,1,verbose)
          WRITE(unitlog(i),*) " - Write the ACFs; also quote the integrated " &
                              //"autocorrelation times (t_int) "
          WRITE(unitlog(i),*) "   & the effective sample sizes (Neff; " &
                              //TRIMLR(TIMINFO(timestr))//")"
        END DO
        
        !! Means
        CALL REALLOCATE(t_int,Nparhyp)
        CALL REALLOCATE(t_int2,Nparhyp)
        CALL REALLOCATE(Neff,Nparhyp)
        CALL REALLOCATE(Neff2,Nparhyp)
        DO ihyp=1,Nparhyp
          CALL AUTOCORREL(mumcmc(ihyp,1:t_end),autocorfun)
          CALL AUTOCORREL(mumcmc(ihyp,t_burnin:t_end),autocorfun2)
          t_int(ihyp) = INTAUTOCORRTIME(autocorfun(:),t_end,NEFF=Neff(ihyp), &
                                        WARN=.False.,FORCE=forcetint)
          t_int2(ihyp) = INTAUTOCORRTIME(autocorfun2(:),Nmcmc_eff, &
                                         NEFF=Neff2(ihyp),WARN=.False., &
                                         FORCE=forcetint)
          ! DO i=1,MERGE(2,1,verbose)
          !   WRITE(unitlog(i),"(A50,F10.1,A3,A10,I10,A3)") &
          !     "mu("//TRIMLR(parhypinfo(ihyp)%name)//"): t_int = ", &
          !     ABS(t_int2(ihyp)), MERGE("   ","(?)",t_int2(ihyp)>0._DP), &
          !     "; Neff =", ABS(Neff2(ihyp)), MERGE("   ","(?)",t_int2(ihyp)>0._DP)
          ! END DO
          CALL WRITE_HDF5(DBLARR1D=autocorfun(:),FILE=fileHB,COMPRESS=compress, &
                          NAME="Autocorrelation function for mu(" &
                               //TRIMLR(parhypinfo(ihyp)%name)//")", &
                          VERBOSE=debug,APPEND=.True.)
          CALL WRITE_HDF5(DBLARR1D=autocorfun2(:),FILE=fileHB,COMPRESS=compress, &
                          NAME="Autocorrelation function for mu(" &
                               //TRIMLR(parhypinfo(ihyp)%name)//") after burn-in", &
                          VERBOSE=debug,APPEND=.True.)
        
        END DO
        t_int_ref1 = MAXVAL(PACK(t_int2(:),MASK=t_int2(:)>0._DP))
        t_int_ref2 = MAXVAL(ABS(t_int2(:)))
        Neff_ref1 = MINVAL(PACK(Neff2(:),MASK=Neff2(:)>0))
        Neff_ref2 = MINVAL(ABS(Neff2(:)))
        CALL WRITE_HDF5(DBLARR1D=t_int(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Integrated autocorrelation time for the means", &
                        VERBOSE=debug,APPEND=.True.)
        CALL WRITE_HDF5(DBLARR1D=t_int2(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Integrated autocorrelation time for the means " &
                           //"after burn-in", &
                        VERBOSE=debug,APPEND=.True.)
        CALL WRITE_HDF5(INTARR1D=Neff(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Effective sample size for the means", &
                        VERBOSE=debug,APPEND=.True.)
        CALL WRITE_HDF5(INTARR1D=Neff2(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Effective sample size for the means after burn-in", &
                        VERBOSE=debug,APPEND=.True.)
        
        !! Standard-deviations
        DO ihyp=1,Nparhyp
          CALL AUTOCORREL(sigmcmc(ihyp,1:t_end),autocorfun)
          CALL AUTOCORREL(sigmcmc(ihyp,t_burnin:t_end),autocorfun2)
          t_int(ihyp) = INTAUTOCORRTIME(autocorfun(:),t_end,NEFF=Neff(ihyp), &
                                        WARN=.False.,FORCE=forcetint)
          t_int2(ihyp) = INTAUTOCORRTIME(autocorfun2(:),Nmcmc_eff, &
                                         NEFF=Neff2(ihyp),WARN=.False., &
                                         FORCE=forcetint)
          ! IF (Nsou > 1) THEN
          !   DO i=1,MERGE(2,1,verbose)
          !     WRITE(unitlog(i),"(A50,F10.1,A3,A10,I10,A3)") &
          !       "sig("//TRIMLR(parhypinfo(ihyp)%name)//"): t_int = ", &
          !       ABS(t_int2(ihyp)), MERGE("   ","(?)",t_int2(ihyp)>0._DP), &
          !       "; Neff =", ABS(Neff2(ihyp)), MERGE("   ","(?)",t_int2(ihyp)>0._DP)
          !   END DO
          ! END IF
          CALL WRITE_HDF5(DBLARR1D=autocorfun(:),FILE=fileHB,COMPRESS=compress, &
                          NAME="Autocorrelation function for sig(" &
                               //TRIMLR(parhypinfo(ihyp)%name)//")", &
                          VERBOSE=debug,APPEND=.True.)
          CALL WRITE_HDF5(DBLARR1D=autocorfun2(:),FILE=fileHB,COMPRESS=compress, &
                          NAME="Autocorrelation function for sig(" &
                               //TRIMLR(parhypinfo(ihyp)%name)//") after burn-in", &
                          VERBOSE=debug,APPEND=.True.)
        END DO
        t_int_ref1 = MAX(MAXVAL(PACK(t_int2(:),MASK=t_int2(:)>0._DP)),t_int_ref1)
        t_int_ref2 = MAX(MAXVAL(ABS(t_int2(:))),t_int_ref2)
        Neff_ref1 = MIN(MINVAL(PACK(Neff2(:),MASK=Neff2(:)>0)),Neff_ref1)
        Neff_ref2 = MIN(MINVAL(ABS(Neff2(:))),Neff_ref2)
        CALL WRITE_HDF5(DBLARR1D=t_int(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Integrated autocorrelation time for the standard " &
                           //"deviations", &
                        VERBOSE=debug,APPEND=.True.)
        CALL WRITE_HDF5(DBLARR1D=t_int2(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Integrated autocorrelation time for the standard " &
                           //"deviations after burn-in", &
                        VERBOSE=debug,APPEND=.True.)
        CALL WRITE_HDF5(INTARR1D=Neff(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Effective sample size for the standard deviations", &
                        VERBOSE=debug,APPEND=.True.)
        CALL WRITE_HDF5(INTARR1D=Neff2(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Effective sample size for the standard deviations " &
                           //"after burn-in", &
                        VERBOSE=debug,APPEND=.True.)
        
        !! Correlations
        CALL REALLOCATE(t_int,Ncorrhyp)
        CALL REALLOCATE(t_int2,Ncorrhyp)
        CALL REALLOCATE(Neff,Ncorrhyp)
        CALL REALLOCATE(Neff2,Ncorrhyp)
        DO icorr=1,Ncorrhyp
          CALL AUTOCORREL(corrmcmc(icorr,1:t_end),autocorfun)
          CALL AUTOCORREL(corrmcmc(icorr,t_burnin:t_end),autocorfun2)
          t_int(icorr) = INTAUTOCORRTIME(autocorfun(:),t_end,NEFF=Neff(icorr), &
                                         WARN=.False.,FORCE=forcetint)
          t_int2(icorr) = INTAUTOCORRTIME(autocorfun2(:),Nmcmc_eff, &
                                          NEFF=Neff2(icorr),WARN=.False., &
                                          FORCE=forcetint)
          ! IF (Nsou > 1) THEN
          !   DO i=1,MERGE(2,1,verbose)
          !     WRITE(unitlog(i),"(A50,F10.1,A3,A10,I10,A3)") &
          !       "corr("//TRIMLR(PRING(icorr))//"): t_int = ", &
          !       ABS(t_int2(icorr)), MERGE("   ","(?)",t_int2(icorr)>0._DP), &
          !       "; Neff =", ABS(Neff2(icorr)), MERGE("   ","(?)",t_int2(icorr)>0._DP)
          !   END DO
          ! END IF
          CALL WRITE_HDF5(DBLARR1D=autocorfun(:),FILE=fileHB,COMPRESS=compress, &
                          NAME="Autocorrelation function for corr(" &
                               //TRIMLR(PRING(icorr))//")", &
                          VERBOSE=debug,APPEND=.True.)
          CALL WRITE_HDF5(DBLARR1D=autocorfun2(:),FILE=fileHB,COMPRESS=compress, &
                          NAME="Autocorrelation function for corr(" &
                               //TRIMLR(PRING(icorr))//") after burn-in", &
                          VERBOSE=debug,APPEND=.True.)
        END DO
        t_int_ref1 = MAX(MAXVAL(PACK(t_int2(:),MASK=t_int2(:)>0._DP)),t_int_ref1)
        t_int_ref2 = MAX(MAXVAL(ABS(t_int2(:))),t_int_ref2)
        Neff_ref1 = MIN(MINVAL(PACK(Neff2(:),MASK=Neff2(:)>0)),Neff_ref1)
        Neff_ref2 = MIN(MINVAL(ABS(Neff2(:))),Neff_ref2)
        CALL WRITE_HDF5(DBLARR1D=t_int(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Integrated autocorrelation time for the " &
                           //"correlations", &
                        VERBOSE=debug,APPEND=.True.)
        CALL WRITE_HDF5(DBLARR1D=t_int2(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Integrated autocorrelation time for the " &
                           //"correlations after burn-in", &
                        VERBOSE=debug,APPEND=.True.)
        CALL WRITE_HDF5(INTARR1D=Neff(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Effective sample size for the correlations", &
                        VERBOSE=debug,APPEND=.True.)
        CALL WRITE_HDF5(INTARR1D=Neff2(:),FILE=fileHB,COMPRESS=compress, &
                        NAME="Effective sample size for the correlations " &
                           //"after burn-in", &
                        VERBOSE=debug,APPEND=.True.)

        !! Worst variables
        DO i=1,MERGE(2,1,verbose)
          WRITE(unitlog(i),*)
          WRITE(unitlog(i),"(A50,F10.1,A10,I10)") &
            "=> Worst case (reliable): t_int = ", t_int_ref1, "; Neff =", Neff_ref1
          WRITE(unitlog(i),"(A50,F10.1,A10,I10)") &
            "=> Worst case (unreliable): t_int = ", t_int_ref2, "; Neff =", Neff_ref2
        END DO
        
      END IF doacf
        
      !! Calculate model
      !!-----------------
      ALLOCATE(FnuMOD(Nx,Ny,NwOBS))
      
      FnuMOD(:,:,:) = specModel( wOBS(:), INDPAR=ind, PARVAL=qpar(:,:,:,2), &
                                 QABS=Qabs(:), EXTINCT=extinct(:,:), &
                                 FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
                                 PABS=Pabs, FNULINE=FnuLINE, &
                                 FNUCONT_TAB=FnuCONT_tab, FNUBAND_TAB=FnuBAND_tab, &
                                 FNUSTAR_TAB=FnuSTAR_tab, PABS_TAB=Pabs_tab, &
                                 FNULINE_TAB=FnuLINE_tab )
      
      CALL WRITE_HDF5(DBLARR4D=FnuCONT_tab, NAME='FnuCONT ('//TRIMLR(spec_unit)//')', &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR4D=FnuLINE_tab, NAME='FnuLINE ('//TRIMLR(spec_unit)//')', &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR4D=FnuBAND_tab, NAME='FnuBAND ('//TRIMLR(spec_unit)//')', &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR4D=FnuSTAR_tab, NAME='FnuSTAR ('//TRIMLR(spec_unit)//')', &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR4D=Pabs_tab, NAME='PABS', &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR3D=FnuMOD, NAME='FnuMOD ('//TRIMLR(spec_unit)//')', &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

      !! Spectrum density
      !!------------------
      ALLOCATE(FnuMOD_mcmc(Nx,Ny,NwOBS,Nmcmc_eff),n_FnuMOD(Nx,Ny,NwOBS,3), &
               mask4D(Nx,Ny,NwOBS,Nmcmc_eff))!, &
               ! FnuCONT_mcmc(Nx,Ny,NwOBS,Nmcmc_eff),n_FnuCONT(Nx,Ny,NwOBS,3), &
               ! FnuBAND_mcmc(Nx,Ny,NwOBS,Nmcmc_eff),n_FnuBAND(Nx,Ny,NwOBS,3), &
               ! FnuLINE_mcmc(Nx,Ny,NwOBS,Nmcmc_eff),n_FnuLINE(Nx,Ny,NwOBS,3), &
               ! FnuSTAR_mcmc(Nx,Ny,NwOBS,Nmcmc_eff),n_FnuSTAR(Nx,Ny,NwOBS,3), &
               ! Pabs_mcmc(Nx,Ny,NwOBS,Nmcmc_eff),n_Pabs(Nx,Ny,NwOBS,3), )

      DO i=1,Nmcmc_eff
        FnuMOD_mcmc(:,:,:,i) = specModel( wOBS(:), INDPAR=ind, &
                                          PARVAL=parmcmc(:,:,:,t_burnin+i-1), &
                                          QABS=Qabs(:), EXTINCT=extinct(:,:) )!, &
                                          ! FNUCONT=FnuCONT_mcmc(:,:,:,i), &
                                          ! FNUBAND=FnuBAND_mcmc(:,:,:,i), &
                                          ! FNUSTAR=FnuSTAR_mcmc(:,:,:,i), &
                                          ! PABS=Pabs_mcmc(:,:,:,i), &
                                          ! FNULINE=FnuLINE_mcmc(:,:,:,i) )
      END DO

      mask4D(:,:,:,:) = .NOT. ISNAN(FnuMOD_mcmc(:,:,:,:))
      n_FnuMOD(:,:,:,2) = MEDIAN(FnuMOD_mcmc(:,:,:,:), DIM=4, MASK=mask4D(:,:,:,:))
      ! n_FnuCONT(:,:,:,2) = MEDIAN(FnuCONT_mcmc(:,:,:,:), DIM=4)
      ! n_FnuBAND(:,:,:,2) = MEDIAN(FnuBAND_mcmc(:,:,:,:), DIM=4)
      ! n_FnuLINE(:,:,:,2) = MEDIAN(FnuLINE_mcmc(:,:,:,:), DIM=4)
      ! n_FnuSTAR(:,:,:,2) = MEDIAN(FnuSTAR_mcmc(:,:,:,:), DIM=4)
      ! n_Pabs(:,:,:,2) = MEDIAN(Pabs_mcmc(:,:,:,:), DIM=4)
      
      !! Calculate 1st and 5th 6-quantiles (sextiles)
      DO x=1,Nx
        DO y=1,Ny
          DO iw=1,NwOBS
            xs(:) = SORT(PACK(FnuMOD_mcmc(x,y,iw,:),MASK=mask4D(x,y,iw,:)))
            n_FnuMOD(x,y,iw,1) = xs(NINT(1._DP/6*Nmcmc_eff)) ! round
            n_FnuMOD(x,y,iw,3) = xs(NINT(5._DP/6*Nmcmc_eff))
          END DO
        END DO
      END DO

      CALL WRITE_HDF5(DBLARR4D=n_FnuMOD, NAME='n_FnuMOD ('//TRIMLR(spec_unit)//')', &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      
      !! Free memory space
      DEALLOCATE(meanpar, stdevpar, qpar, meanmu, stdevmu, meansig, stdevsig, &
                 meancorr, stdevcorr, xs, FnuMOD, FnuMOD_mcmc, n_FnuMOD, mask4D)

    END IF
  END DO
    
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'HISTOIRE: fit_hb analysis [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO

END PROGRAM analysis
