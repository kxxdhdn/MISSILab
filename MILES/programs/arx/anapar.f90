PROGRAM anapar

  USE auxil, ONLY: read_master, specModel
  USE utilities, ONLY: DP, pring, trimLR, trimeq, timinfo, &
                       banner_program, ustd, isNaN
  USE arrays, ONLY: iwhere, closest
  USE statistics, ONLY: mean, sigma
  USE grain_optics, ONLY: lendustQ
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, lenpar
  USE HB_kit, ONLY: Nx, Ny, NwOBS, FnuOBS, Nparhyp, Ncorrhyp, &
                    ind, parinfo, parhypinfo, Qabs, extinct
  IMPLICIT NONE

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: textwid = 60
  CHARACTER(*), PARAMETER :: dirIN = '../out1/'
  CHARACTER(*), PARAMETER :: filOBS = dirIN//'galspec'//h5ext
  CHARACTER(*), PARAMETER :: filMCMC = dirOUT//'parlog_fitpar_HB'//h5ext
  CHARACTER(*), PARAMETER :: filHB = dirOUT//'fitpar_HB'//h5ext
  CHARACTER(*), PARAMETER :: filCHI2 = dirOUT//'fitpar_chi2'//h5ext
  CHARACTER(*), PARAMETER :: nampar = "Parameter values"
  CHARACTER(*), PARAMETER :: nammu = "Mean of hyperdistribution"
  CHARACTER(*), PARAMETER :: namsig = "Sigma of hyperdistribution"
  CHARACTER(*), PARAMETER :: namcorr = "Correlation of hyperdistribution"
  CHARACTER(*), PARAMETER :: namrat = "Band ratio values"
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .TRUE.
  
  !! Input variables
  INTEGER :: i, Npar, Nmcmc, NiniMC, Nrat, Ncorr
  INTEGER :: Ncont, Nband, Nline, Npabs, Nstar, Nextra
  INTEGER, DIMENSION(2) :: unitlog
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wOBS
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: parerr, par
  CHARACTER(lenpath) :: dirOUT, filMCMC, filOUT
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB, labL, ratname, namBcorr
  CHARACTER(lenpar) :: spec_unit
  LOGICAL :: verbose, calib, newseed, newinit, dostop
  !! Analysis variables
  INTEGER :: t_burnin, t_end, Nmcmc_eff
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: meanpar, stdevpar
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: meanrat, stdevrat
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meanmu, stdevmu, meansig, stdevsig
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meancorr, stdevcorr
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mumcmc, sigmcmc, corrmcmc
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT, &
    FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FnuCONT_tab, &
    FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parmcmc, ratmcmc, parpix, ratpix
  LOGICAL :: chi2ana, HBana

  !! Log output
  !!------------
  ! CALL INITIATE_CLOCK(timestr)
  ! OPEN (ulog,FILE=fiLOG,STATUS="REPLACE",ACTION="WRITE")
  unitlog(:) = [ ulog, ustd ]

  !! Select fits to analyse
  chi2ana = .FALSE.
  HBana = .FALSE.
  
  ! chi2ana = .TRUE.
  HBana = .TRUE.

  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='Wavelength (microns)', N1=NwOBS)
  CALL READ_MASTER(WAVALL=wOBS(:), DIRIN=dirIN, DIROUT=dirOUT, &
                   VERBOSE=verbose, NMCMC=Nmcmc, NINIMC=NiniMC, &
                   CALIB=calib, NEWSEED=newseed, NEWINIT=newinit, &
                   LABL=labL, LABB=labB, QABS=Qabs, EXTINCT=extinct, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, &
                   NPABS=Npabs, NSTAR=Nstar, NEXTRA=Nextra, &
                   DOSTOP=dostop, NCORR=Ncorr, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit, &
                   PARHYPINFO=parhypinfo, NPARHYP=Nparhyp, NCORRHYP=Ncorrhyp)
  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, &
                 NAME='FnuOBS ('//TRIMLR(spec_unit)//')', N1=Nx, N2=Ny)
  CALL READ_HDF5(STRARR1D=ratname, FILE=filOBS, NAME='Band ratio name', N1=Nrat)
  CALL READ_HDF5(STRARR1D=namBcorr, FILE=filOBS, NAME='Correlated band name')

  !!------------------------------------------------------------------------
  !!                             Chi2 Fit
  !!------------------------------------------------------------------------
  IF (chi2ana) THEN
    PRINT*, 'MIROIR (chi2 fit) analysis'

    filOUT = TRIMLR(dirOUT)

    CALL READ_HDF5(DBLARR3D=par, FILE=filCHI2, NAME='Best fitted parameter value')
    CALL READ_HDF5(DBLARR3D=parerr, FILE=filCHI2, NAME='Best fitted parameter error')

    !! Calculate model
    !!-----------------
    ALLOCATE(FnuMOD(Nx,Ny,NwOBS))
    
    FnuMOD(:,:,:) = specModel( wOBS(:), INDPAR=ind, PARVAL=par(:,:,:), &
                               QABS=Qabs(:), EXTINCT=extinct(:), &
                               FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
                               PABS=Pabs, FNULINE=FnuLINE, &
                               FNUCONT_TAB=FnuCONT_tab, FNUBAND_TAB=FnuBAND_tab, &
                               FNUSTAR_TAB=FnuSTAR_tab, PABS_TAB=Pabs_tab, &
                               FNULINE_TAB=FnuLINE_tab )
    
    CALL WRITE_HDF5(DBLARR4D=FnuCONT_tab, NAME='FnuCONT ('//TRIMLR(spec_unit)//')', &
                    FILE=filCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR4D=FnuLINE_tab, NAME='FnuLINE ('//TRIMLR(spec_unit)//')', &
                    FILE=filCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR4D=FnuBAND_tab, NAME='FnuBAND ('//TRIMLR(spec_unit)//')', &
                    FILE=filCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR4D=FnuSTAR_tab, NAME='FnuSTAR ('//TRIMLR(spec_unit)//')', &
                    FILE=filCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR4D=Pabs_tab, NAME='PABS', &
                    FILE=filCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR3D=FnuMOD, NAME='FnuMOD ('//TRIMLR(spec_unit)//')', &
                    FILE=filCHI2, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    
    DO i=1,MERGE(2,1,verbose)
      WRITE(unitlog(i),*)
      WRITE(unitlog(i),*) 'fitpar_chi2 analysis [done]'
      WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
    END DO
    
    !! Final
    !!-------
    ! DO i=1,MERGE(2,1,verbose)
    !   WRITE(unitlog(i),*)
    !   WRITE(unitlog(i),*) "PROGRAM EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
    !   WRITE(unitlog(i),*)
    ! END DO
    
    !! Free memory space
    DEALLOCATE(wOBS, FnuOBS)
    
  END IF

  !!------------------------------------------------------------------------
  !!                              HB Fit
  !!------------------------------------------------------------------------
  IF (HBana) THEN
    PRINT*, 'HIBARI (HB fit) analysis'

    !! Compute the statistics
    !!------------------------
    t_burnin = 2000000
    t_end = Nmcmc ! end of Markov chain Monte Carlo
    IF (t_burnin <= 0 .OR. t_burnin >= t_end) t_burnin = INT(t_end * 0.8_DP)
    Nmcmc_eff = t_end - t_burnin + 1
    
    ALLOCATE(parpix(Nx,Ny,Npar,t_end), meanpar(Nx,Ny,Npar), stdevpar(Nx,Ny,Npar))
    ALLOCATE(ratpix(Nx,Ny,Nrat,t_end), meanrat(Nx,Ny,Nrat), stdevrat(Nx,Ny,Nrat))
    
    CALL READ_HDF5(DBLARR4D=parmcmc, FILE=filMCMC, &
                   NAME=nampar, IND4=[1,t_end])
    parpix(:,:,:,:) = parmcmc(:,:,:,:)
    meanpar(:,:,:) = MEAN(parpix(:,:,:,t_burnin:t_end), DIM=4)
    stdevpar(:,:,:) = SIGMA(parpix(:,:,:,t_burnin:t_end), DIM=4)

    CALL READ_HDF5(DBLARR4D=ratmcmc, FILE=filMCMC, &
                   NAME=namrat, IND4=[1,t_end])
    ratpix(:,:,:,:) = ratmcmc(:,:,:,:)
    meanrat(:,:,:) = MEAN(ratpix(:,:,:,t_burnin:t_end), DIM=4)
    stdevrat(:,:,:) = SIGMA(ratpix(:,:,:,t_burnin:t_end), DIM=4)
    
    CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME="Parameter label", &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
    CALL WRITE_HDF5(DBLARR3D=meanpar(:,:,:), NAME="Mean of parameter value", &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR3D=stdevpar(:,:,:), NAME="Sigma of parameter value", &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(STRARR1D=namBcorr(:), NAME="Correlated band label", &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(STRARR1D=ratname(:), NAME="Band ratio label", &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR3D=meanrat(:,:,:), NAME="Mean of band ratio value", &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR3D=stdevrat(:,:,:), NAME="Sigma of band ratio value", &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

    !! Compute the average, standard deviations and correlations of each quantity
    !!----------------------------------------------------------------------------
    CALL READ_HDF5(DBLARR2D=mumcmc, FILE=filMCMC, &
                   NAME=nammu, IND2=[1,t_end])
    CALL READ_HDF5(DBLARR2D=sigmcmc, FILE=filMCMC, &
                   NAME=namsig, IND2=[1,t_end])
    CALL READ_HDF5(DBLARR2D=corrmcmc, FILE=filMCMC, &
                   NAME=namcorr, IND2=[1,t_end])
    
    ! DO i=1,MERGE(2,1,verbose)
    !   WRITE(unitlog(i),*) " - Compute the statistics ("//TRIMLR(TIMINFO(timestr)) &
    !                       //")"
    !   WRITE(unitlog(i),*)
    ! END DO
    
    ALLOCATE (meanmu(Nparhyp),stdevmu(Nparhyp),meansig(Nparhyp),stdevsig(Nparhyp))
    
    meanmu(:) = MEAN(mumcmc(:,t_burnin:t_end), DIM=2)
    stdevmu(:) = SIGMA(mumcmc(:,t_burnin:t_end), DIM=2)
    meansig(:) = MEAN(sigmcmc(:,t_burnin:t_end), DIM=2)
    stdevsig(:) = SIGMA(sigmcmc(:,t_burnin:t_end), DIM=2)
    
    ALLOCATE (meancorr(Ncorrhyp),stdevcorr(Ncorrhyp))
    
    meancorr(:) = MEAN(corrmcmc(:,t_burnin:t_end), DIM=2)
    stdevcorr(:) = SIGMA(corrmcmc(:,t_burnin:t_end), DIM=2)
    
    CALL WRITE_HDF5(DBLARR1D=meanmu(:),NAME="Mean of parameter average", &
                    FILE=filHB,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)
    CALL WRITE_HDF5(DBLARR1D=stdevmu(:),NAME="Sigma of parameter average", &
                    FILE=filHB,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)
    CALL WRITE_HDF5(DBLARR1D=meansig(:), &
                    NAME="Mean of parameter standard deviation", &
                    FILE=filHB,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)
    CALL WRITE_HDF5(DBLARR1D=stdevsig(:), &
                    NAME="Sigma of parameter standard deviation", &
                    FILE=filHB,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)
    CALL WRITE_HDF5(DBLARR1D=meancorr(:),NAME="Mean of parameter correlation", &
                    FILE=filHB,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)
    CALL WRITE_HDF5(DBLARR1D=stdevcorr(:),NAME="Sigma of parameter correlation", &
                    FILE=filHB,COMPRESS=compress,VERBOSE=debug,APPEND=.True.)
    
    !! Calculate model
    !!-----------------
    ALLOCATE(FnuMOD(Nx,Ny,NwOBS))
    
    FnuMOD(:,:,:) = specModel( wOBS(:), INDPAR=ind, PARVAL=meanpar(:,:,:), &
                               QABS=Qabs(:), EXTINCT=extinct(:), &
                               FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
                               PABS=Pabs, FNULINE=FnuLINE, &
                               FNUCONT_TAB=FnuCONT_tab, FNUBAND_TAB=FnuBAND_tab, &
                               FNUSTAR_TAB=FnuSTAR_tab, PABS_TAB=Pabs_tab, &
                               FNULINE_TAB=FnuLINE_tab )
    
    DO i=1,Npabs
      FnuCONT_tab(:,:,:,i) = FnuCONT_tab(:,:,:,i) * Pabs(:,:,:)
    END DO 
    CALL WRITE_HDF5(DBLARR4D=FnuCONT_tab, NAME='FnuCONT ('//TRIMLR(spec_unit)//')', &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    DO i=1,Nline
      FnuLINE_tab(:,:,:,i) = FnuLINE_tab(:,:,:,i) + (FnuCONT(:,:,:)+FnuSTAR(:,:,:))*Pabs(:,:,:)
    END DO
    CALL WRITE_HDF5(DBLARR4D=FnuLINE_tab, NAME='FnuLINE ('//TRIMLR(spec_unit)//')', &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    DO i=1,Nband
      FnuBAND_tab(:,:,:,i) = FnuBAND_tab(:,:,:,i) + (FnuCONT(:,:,:)+FnuSTAR(:,:,:))*Pabs(:,:,:)
    END DO
    CALL WRITE_HDF5(DBLARR4D=FnuBAND_tab, NAME='FnuBAND ('//TRIMLR(spec_unit)//')', &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    DO i=1,Nstar
      FnuSTAR_tab(:,:,:,i) = FnuSTAR_tab(:,:,:,i) * Pabs(:,:,:)
    END DO
    CALL WRITE_HDF5(DBLARR4D=FnuSTAR_tab, NAME='FnuSTAR ('//TRIMLR(spec_unit)//')', &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR3D=FnuMOD, NAME='FnuMOD ('//TRIMLR(spec_unit)//')', &
                    FILE=filHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    
    DO i=1,MERGE(2,1,verbose)
      WRITE(unitlog(i),*)
      WRITE(unitlog(i),*) 'fitpar_HB analysis [done]'
      WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
    END DO
    
    !! Final
    !!-------
    ! DO i=1,MERGE(2,1,verbose)
    !   WRITE(unitlog(i),*)
    !   WRITE(unitlog(i),*) "PROGRAM EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
    !   WRITE(unitlog(i),*)
    ! END DO
    
    !! Free memory space
    DEALLOCATE(wOBS, FnuMOD)

  END IF
  
END PROGRAM anapar
