!******************************************************************************
!*
!*                      Fit Simulation Parameters (HB)
!*
!******************************************************************************


  !==========================================================================
  ! 1) AUTHOR: D. HU
  ! 
  ! 2) DESCRIPTION: NONE
  !
  ! 3) HISTORY: 
  !    - 20210120: Created.
  !==========================================================================
PROGRAM anapar

  USE auxil, ONLY: read_master, specModel
  USE utilities, ONLY: DP, pring, trimLR, trimeq, &
                       banner_program, ustd, isNaN
  USE arrays, ONLY: iwhere, closest, reallocate
  USE statistics, ONLY: mean, sigma
  USE grain_optics, ONLY: lendustQ
  USE inout, ONLY: read_hdf5, write_hdf5, check_hdf5, h5ext, lenpar, lenpath
  USE HB_kit, ONLY: Nx, Ny, NwOBS, FnuOBS, Nparhyp, Ncorrhyp, &
                    ind, parinfo, Qabs, extinct
  IMPLICIT NONE

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: textwid = 60
  CHARACTER(*), PARAMETER :: nampar = "Parameter values"
  CHARACTER(*), PARAMETER :: nammu = "Mean of hyperdistribution"
  CHARACTER(*), PARAMETER :: namsig = "Sigma of hyperdistribution"
  CHARACTER(*), PARAMETER :: namcorr = "Correlation of hyperdistribution"
  CHARACTER(*), PARAMETER :: dirIN = '../out1/'
  CHARACTER(*), PARAMETER :: filOBS = dirIN//'galspec'//h5ext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Input variables
  INTEGER :: i, j, Npar, Nmcmc!, NiniMC, Ncorr
  ! INTEGER :: Ncont, Nband, Nline, Nextc, Nstar, Nextra
  INTEGER, DIMENSION(2) :: unitlog
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wOBS
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: parerr, par
  CHARACTER(lenpath) :: dirOUT, filMCMC, fileCHI2, fileHB
  CHARACTER(lenpar) :: spec_unit
  LOGICAL :: verbose!, calib, newseed, newinit, dostop
  !! Analysis variables
  INTEGER :: t_burnin, t_end, Nmcmc_eff
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: meanpar, stdevpar
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meanmu, stdevmu, meansig, stdevsig
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meancorr, stdevcorr
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mumcmc, sigmcmc, corrmcmc
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT, &
    FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FnuCONT_tab, &
    FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parmcmc
  LOGICAL :: chi2OK, HBOK

  !! Preliminaires
  !!---------------
  unitlog(:) = [ ulog, ustd ]

  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='Wavelength (microns)', N1=NwOBS)
  CALL READ_MASTER(WAVALL=wOBS(:), DIRIN=dirIN, DIROUT=dirOUT, &
                   VERBOSE=verbose, NMCMC=Nmcmc, &!NINIMC=NiniMC, &
                   ! CALIB=calib, NEWSEED=newseed, NEWINIT=newinit, &
                   QABS=Qabs, EXTINCT=extinct, &
                   ! NCONT=Ncont, NBAND=Nband, NLINE=Nline, &
                   ! NEXTC=Nextc, NSTAR=Nstar, NEXTRA=Nextra, &
                   ! DOSTOP=dostop, NCORR=Ncorr, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit, &
                   NPARHYP=Nparhyp, NCORRHYP=Ncorrhyp)
  
  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, &
                 NAME='FnuOBS ('//TRIMLR(spec_unit)//')', N1=Nx, N2=Ny)

  !!------------------------------------------------------------------------
  !!                             Chi2 Fit
  !!------------------------------------------------------------------------
  
  fileCHI2 = TRIMLR(dirOUT)//'fitpar_chi2'//h5ext
  chi2OK = CHECK_HDF5(FILE=fileCHI2)
  IF (chi2OK) THEN
    PRINT*, 'LE MIROIR (chi2 fit) analysis'

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
      WRITE(unitlog(i),*) 'fitpar_chi2 analysis [done]'
      WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
    END DO
    
    !! Free memory space
    DEALLOCATE(FnuMOD)
    
  END IF

  !!------------------------------------------------------------------------
  !!                              HB Fit
  !!------------------------------------------------------------------------

  DO j=1,2
    IF (j==1) THEN
      fileHB = TRIMLR(dirOUT)//'fitpar_BB'//h5ext
      filMCMC = TRIMLR(dirOUT)//'parlog_fitpar_BB'//h5ext
    ELSE
      fileHB = TRIMLR(dirOUT)//'fitpar_HB'//h5ext
      filMCMC = TRIMLR(dirOUT)//'parlog_fitpar_HB'//h5ext
    END IF
    HBOK = CHECK_HDF5(FILE=filMCMC)
    IF (HBOK) THEN
      PRINT*, 'HIBARI (HB fit) analysis'
  
      !! Compute the statistics
      !!------------------------
      t_burnin = 2000000
      t_end = Nmcmc ! end of Markov chain Monte Carlo
      IF (t_burnin <= 0 .OR. t_burnin >= t_end) t_burnin = INT(t_end * 0.5_DP)
      Nmcmc_eff = t_end - t_burnin + 1
      
      ALLOCATE(meanpar(Nx,Ny,Npar), stdevpar(Nx,Ny,Npar))
      
      CALL READ_HDF5(DBLARR4D=parmcmc, FILE=filMCMC, &
                     NAME=nampar, IND4=[1,t_end])
      meanpar(:,:,:) = MEAN(parmcmc(:,:,:,t_burnin:t_end), DIM=4)
      stdevpar(:,:,:) = SIGMA(parmcmc(:,:,:,t_burnin:t_end), DIM=4)
      
      CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME="Parameter label", &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
      CALL WRITE_HDF5(DBLARR3D=meanpar(:,:,:), NAME="Mean of parameter value", &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR3D=stdevpar(:,:,:), NAME="Sigma of parameter value", &
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
  
      !! Calculate model
      !!-----------------
      ALLOCATE(FnuMOD(Nx,Ny,NwOBS))
      
      FnuMOD(:,:,:) = specModel( wOBS(:), INDPAR=ind, PARVAL=meanpar(:,:,:), &
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

      !! Free memory space
      DEALLOCATE(meanpar, stdevpar, meanmu, stdevmu, meansig, stdevsig, &
                 meancorr, stdevcorr, FnuMOD)

    END IF
  END DO
    
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'fitpar_HB analysis [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO

END PROGRAM anapar
