!******************************************************************************
!*
!*                           MILES POST-PROCESSING
!*
!******************************************************************************


PROGRAM presimu_an

  USE core, ONLY: read_master, read_analysis, specModel
  USE utilities, ONLY: DP, pring, trimLR, &
                       banner_program, ustd, isNaN, NaN, warning, &
                       trimeq, time_type, timinfo, initiate_clock
  USE arrays, ONLY: closest, reallocate, sort
  USE statistics, ONLY: mean, sigma, median
  USE grain_optics, ONLY: lendustQ
  USE inout, ONLY: read_hdf5, write_hdf5, check_hdf5, h5ext, ascext, &
                   lenpar, lenpath
  USE ext_hb, ONLY: Nx, Ny, NwOBS, FnuOBS, Nparhyp, Ncorrhyp, &
                    mask, ind, parinfo, parhypinfo, Qabs, extinct, &
                    specOBS, calib
  IMPLICIT NONE

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: textwid = 60
  CHARACTER(*), PARAMETER :: nampar = "Parameter values"
  CHARACTER(*), PARAMETER :: namln1pd = "ln(1+delta)"
  CHARACTER(*), PARAMETER :: nammu = "Mean of hyperdistribution"
  CHARACTER(*), PARAMETER :: namsig = "Sigma of hyperdistribution"
  CHARACTER(*), PARAMETER :: namcorr = "Correlation of hyperdistribution"
  LOGICAL, PARAMETER :: compress = .TRUE., forcetint = .TRUE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Input variables
  INTEGER :: i, x, y, iw, ical, ipar, it
  INTEGER :: Npar, Nmcmc, lastind, Ncalib
  ! INTEGER :: Ncont, Nband, Nline, Nextc, Nstar, Nextra
  INTEGER, DIMENSION(:), ALLOCATABLE :: indresume
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: maskint ! convert mask=0 to mask=T
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wOBS
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: parerr, par
  CHARACTER(lenpath) :: dirIN, dirOUT, filOBS, filMCMC, fiLOG, fileCHI2, fileHB
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lenpath), DIMENSION(:), ALLOCATABLE :: path1d
  !! Analysis variables
  INTEGER :: Nmcmc_eff, Nsou
  INTEGER, DIMENSION(2) :: t_burnin, t_end
  REAL(DP), DIMENSION(:), ALLOCATABLE :: xs
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meanmu, stdevmu, meansig, stdevsig
  REAL(DP), DIMENSION(:), ALLOCATABLE :: meancorr, stdevcorr
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mumcmc, sigmcmc, corrmcmc
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: FnuMOD_mcmc, n_FnuMOD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: meanln1pd, stdevln1pd
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: meanpar, stdevpar
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT, &
    FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FnuCONT_tab, &
    FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: ln1pdmcmc, parmcmc, qpar, qcal
  LOGICAL :: chi2OK, hbOK, verbose
  LOGICAL, DIMENSION(2) :: ACF
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: mask2D
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
                   NMCMC=Nmcmc, QABS=Qabs, EXTINCT=extinct, CALIB=calib, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit, &
                   PARHYPINFO=parhypinfo, NPARHYP=Nparhyp, NCORRHYP=Ncorrhyp)
  
  fiLOG = TRIMLR(dirOUT)//'log_presimu_an'//ascext

  CALL INITIATE_CLOCK(timestr)
  OPEN (ulog,FILE=fiLOG,STATUS="REPLACE",ACTION="WRITE")
  unitlog(:) = [ ulog, ustd ]
  
  CALL READ_HDF5(STRARR1D=specOBS, FILE=filOBS, &
                 NAME='spectroscopic module labels')
  Ncalib = SIZE(specOBS)
  
  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, &
                 NAME='FnuOBS ('//TRIMLR(spec_unit)//')', N1=Nx, N2=Ny)
  CALL READ_HDF5(INTARR3D=maskint, FILE=filOBS, NAME='NaN mask')
  ALLOCATE(mask(Nx,Ny,NwOBS))
  mask(:,:,:) = ( maskint(:,:,:) == 0 )
  Nsou = COUNT(ANY(mask(:,:,:),DIM=3))

  !!------------------------------------------------------------------------
  !!                      LE MIROIR (Chi2 Fit) Analysis
  !!------------------------------------------------------------------------
  
  fileCHI2 = TRIMLR(dirOUT)//'presimu_chi2'//h5ext
  chi2OK = CHECK_HDF5(FILE=fileCHI2)
  IF (chi2OK) THEN
    PRINT*, 'LE MIROIR: presimu_chi2 analysis'

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
                               MASK=mask(:,:,:), QABS=Qabs(:), EXTINCT=extinct(:,:), &
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
      WRITE(unitlog(i),*) 'LE MIROIR: presimu_chi2 analysis [done]'
      WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
    END DO
    
    !! Free memory space
    DEALLOCATE(FnuMOD)
    
  END IF

  !!------------------------------------------------------------------------
  !!                        HISTOIRE (HB Fit) Analysis
  !!------------------------------------------------------------------------

  fileHB = TRIMLR(dirOUT)//'presimu_hb'//h5ext
  filMCMC = TRIMLR(dirOUT)//'parlog_presimu_hb'//h5ext
  hbOK = CHECK_HDF5(FILE=filMCMC)
  IF (hbOK) THEN

    PRINT*, 'HISTOIRE: presimu_hb analysis'

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
    IF (t_end(1) <= 0 .OR. t_end(1) > Nmcmc) t_end(1) = Nmcmc
    IF (t_end(1) > lastind) t_end(1) = lastind
    IF (t_burnin(1) <= 0 .OR. t_burnin(1) >= t_end(1)) &
      t_burnin(1) = INT(t_end(1) * 0.2_DP) ! By default, burn 20 percents

    Nmcmc_eff = t_end(1) - t_burnin(1) + 1
    
    ALLOCATE(meanpar(Nx,Ny,Npar), stdevpar(Nx,Ny,Npar), qpar(Nx,Ny,Npar,3))
    
    CALL READ_HDF5(DBLARR4D=parmcmc, FILE=filMCMC, &
                   NAME=nampar, IND4=[1,t_end(1)])
    meanpar(:,:,:) = MEAN(parmcmc(:,:,:,t_burnin(1):t_end(1)), DIM=4)
    stdevpar(:,:,:) = SIGMA(parmcmc(:,:,:,t_burnin(1):t_end(1)), DIM=4)

    !! Calculate 2nd and 4th 6-quantiles (sextiles)
    ALLOCATE(xs(Nmcmc_eff), mask4D(Nx,Ny,Npar,Nmcmc_eff))
    mask4D(:,:,:,:) = .NOT. ISNAN(parmcmc(:,:,:,t_burnin(1):t_end(1)))
    ! qpar(:,:,:,2) = MEDIAN(parmcmc(:,:,:,t_burnin(1):t_end(1)), DIM=4, &
    !                        MASK=mask4D(:,:,:,:))
    DO x=1,Nx
      DO y=1,Ny
        DO ipar=1,Npar
          xs(:) = SORT(PACK(parmcmc(x,y,ipar,t_burnin(1):t_end(1)), &
                            MASK=mask4D(x,y,ipar,:)))
          qpar(x,y,ipar,1) = xs(NINT(2._DP/6*Nmcmc_eff)) ! round
          qpar(x,y,ipar,2) = xs(NINT(3._DP/6*Nmcmc_eff))
          qpar(x,y,ipar,3) = xs(NINT(4._DP/6*Nmcmc_eff))
        END DO
      END DO
    END DO
    DEALLOCATE(mask4D, xs)

    CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME="Parameter label", &
                    FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
    CALL WRITE_HDF5(DBLARR3D=meanpar(:,:,:), NAME="Mean of parameter value", &
                    FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR3D=stdevpar(:,:,:), NAME="Sigma of parameter value", &
                    FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    CALL WRITE_HDF5(DBLARR4D=qpar(:,:,:,:), NAME="Quantiles of parameter value", &
                    FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

    !! Calibration errors
    !!--------------------
    IF (calib) THEN
      CALL READ_HDF5(DBLARR4D=ln1pdmcmc, FILE=filMCMC, &
                     NAME=namln1pd, IND4=[1,t_end(1)])

      ALLOCATE (meanln1pd(Nx,Ny,Ncalib),stdevln1pd(Nx,Ny,Ncalib),qcal(Nx,Ny,Ncalib,3))
      
      meanln1pd(:,:,:) = MEAN(ln1pdmcmc(:,:,:,t_burnin(1):t_end(1)), DIM=4)
      stdevln1pd(:,:,:) = SIGMA(ln1pdmcmc(:,:,:,t_burnin(1):t_end(1)), DIM=4)

      !! Calculate 2nd and 4th 6-quantiles (sextiles)
      ALLOCATE(xs(Nmcmc_eff), mask4D(Nx,Ny,Npar,Nmcmc_eff))
      mask4D(:,:,:,:) = .NOT. ISNAN(ln1pdmcmc(:,:,:,t_burnin(1):t_end(1)))
      ! qcal(:,:,:,2) = MEDIAN(ln1pdmcmc(:,:,:,t_burnin(1):t_end(1)), DIM=4, &
      !                        MASK=mask4D(:,:,:,:))
      DO x=1,Nx
        DO y=1,Ny
          DO ical=1,Ncalib
            xs(:) = SORT(PACK(ln1pdmcmc(x,y,ical,t_burnin(1):t_end(1)), &
                              MASK=mask4D(x,y,ical,:)))
            qcal(x,y,ical,1) = xs(NINT(2._DP/6*Nmcmc_eff)) ! round
            qcal(x,y,ical,2) = xs(NINT(3._DP/6*Nmcmc_eff))
            qcal(x,y,ical,3) = xs(NINT(4._DP/6*Nmcmc_eff))
          END DO
        END DO
      END DO
      DEALLOCATE(mask4D, xs)

      CALL WRITE_HDF5(DBLARR3D=meanln1pd(:,:,:), NAME="Mean of ln(1+delta)", &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR3D=stdevln1pd(:,:,:),NAME="Sigma of ln(1+delta)", &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
      CALL WRITE_HDF5(DBLARR4D=qcal(:,:,:,:), NAME="Quantiles of ln(1+delta)", &
                      FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
    ELSE
      ALLOCATE(meanln1pd(Nx,Ny,NwOBS),stdevln1pd(Nx,Ny,NwOBS),qcal(Nx,Ny,NwOBS,3))
      meanln1pd(:,:,:) = 0._DP
      stdevln1pd(:,:,:) = 0._DP
      qcal(:,:,:,:) = 0._DP
    END IF
    
    !! Compute the average, standard deviations and correlations of each quantity
    !!----------------------------------------------------------------------------
    CALL READ_HDF5(DBLARR2D=mumcmc, FILE=filMCMC, &
                   NAME=nammu, IND2=[1,t_end(1)])
    CALL READ_HDF5(DBLARR2D=sigmcmc, FILE=filMCMC, &
                   NAME=namsig, IND2=[1,t_end(1)])
    CALL READ_HDF5(DBLARR2D=corrmcmc, FILE=filMCMC, &
                   NAME=namcorr, IND2=[1,t_end(1)])
  
    ALLOCATE (meanmu(Nparhyp),stdevmu(Nparhyp),meansig(Nparhyp),stdevsig(Nparhyp))
    
    meanmu(:) = MEAN(mumcmc(:,t_burnin(1):t_end(1)), DIM=2)
    stdevmu(:) = SIGMA(mumcmc(:,t_burnin(1):t_end(1)), DIM=2)
    meansig(:) = MEAN(sigmcmc(:,t_burnin(1):t_end(1)), DIM=2)
    stdevsig(:) = SIGMA(sigmcmc(:,t_burnin(1):t_end(1)), DIM=2)
    
    ALLOCATE (meancorr(Ncorrhyp),stdevcorr(Ncorrhyp))
    
    meancorr(:) = MEAN(corrmcmc(:,t_burnin(1):t_end(1)), DIM=2)
    stdevcorr(:) = SIGMA(corrmcmc(:,t_burnin(1):t_end(1)), DIM=2)

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

    DEALLOCATE(meanln1pd, stdevln1pd)
    DEALLOCATE(meanmu, stdevmu, meansig, stdevsig, meancorr, stdevcorr)
      
    !! Calculate model
    !!-----------------
    ALLOCATE(FnuMOD(Nx,Ny,NwOBS))
    
    FnuMOD(:,:,:) = specModel( wOBS(:), INDPAR=ind, PARVAL=qpar(:,:,:,2), &
                               MASK=mask(:,:,:), QABS=Qabs(:), EXTINCT=extinct(:,:), &
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
    ALLOCATE(FnuMOD_mcmc(NwOBS,Nmcmc_eff),n_FnuMOD(NwOBS,3))

    CALL WRITE_HDF5(INITDBLARR=[Nx,Ny,NwOBS,3], NAME='n_FnuMOD ('//TRIMLR(spec_unit)//')', &
                    FILE=fileHB, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

    DO x=1,Nx
      DO y=1,Ny
        FnuMOD_mcmc(:,:) = NaN()
        n_FnuMOD(:,:) = NaN()
        IF (ANY(mask(x,y,:))) THEN
          DO i=1,MERGE(2,1,verbose)
            WRITE(unitlog(i),*) "CALCULATING sextiles (" &
                 //TRIMLR(PRING(x))//","//TRIMLR(PRING(y))//") - " &
                 //TRIMLR(TIMINFO(timestr))//"."!//NEW_LINE('')
          END DO
          
          !! Calculate MCMC models
          DO it=1,Nmcmc_eff
            FnuMOD_mcmc(:,it) = specModel( wOBS(:), INDPAR=ind, &
                                           PARVAL=parmcmc(x,y,:,t_burnin(1)+it-1), &
                                           ! MASK=mask(:,:,:), &
                                           QABS=Qabs(:), EXTINCT=extinct(:,:) )
          END DO
          
          !! Calculate 2nd, 3rd (median) and 4th 6-quantiles (sextiles)
          ALLOCATE(mask2D(NwOBS,Nmcmc_eff), xs(Nmcmc_eff))      
          mask2D(:,:) = .NOT. ISNAN(FnuMOD_mcmc(:,:))
          DO iw=1,NwOBS
            xs(:) = SORT(PACK(FnuMOD_mcmc(iw,:),MASK=mask2D(iw,:)))
            n_FnuMOD(iw,1) = xs(NINT(2._DP/6*Nmcmc_eff)) ! round
            n_FnuMOD(iw,2) = xs(NINT(3._DP/6*Nmcmc_eff))
            n_FnuMOD(iw,3) = xs(NINT(4._DP/6*Nmcmc_eff))
          END DO
          DEALLOCATE(xs, mask2D)
        END IF
        
        CALL WRITE_HDF5(DBLARR4D=RESHAPE(n_FnuMOD(:,:),[1,1,NwOBS,3]), &
                        NAME='n_FnuMOD ('//TRIMLR(spec_unit)//')', &
                        FILE=fileHB, COMPRESS=compress, VERBOSE=debug, &
                        IND1=[x,x],IND2=[y,y])          
      END DO
    END DO

    !! Free memory space
    DEALLOCATE(meanpar, stdevpar, qpar, qcal, FnuMOD, FnuMOD_mcmc, n_FnuMOD)

  END IF
    
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'HISTOIRE: presimu_hb analysis [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO

END PROGRAM presimu_an
