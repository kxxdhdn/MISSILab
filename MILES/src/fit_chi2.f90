!******************************************************************************
!*
!*                          FIT MIR SPECTRA (CHI2)
!*
!******************************************************************************


PROGRAM fit_chi2

  USE utilities, ONLY: DP, pring, trimLR, trimeq, timinfo, &
                       banner_program, ustd, isNaN, initiate_clock, time_type
  USE arrays, ONLY: iwhere, closest
  USE constants, ONLY: MKS
  USE statistics, ONLY: MEDIAN
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, ascext, lenpar, lenpath
  USE grain_optics, ONLY: lendustQ
  USE random, ONLY: generate_newseed
  USE chi2_minimization, ONLY: chi2min_LM
  USE core, ONLY: read_master, specModel, initparam
  USE ext_chi2, ONLY: Nx, Ny, NwOBS, wOBS, nuOBS, FnuOBS, dFnuOBS, &
                      residuals, Fnu_mod, resid, invLcovarOBS, &
                      ind, parinfo, Qabs, extinct, xOBS, yOBS, mask, iwfree
  IMPLICIT NONE

  !! Parameters
  INTEGER, PARAMETER :: ulog = 2
  INTEGER, PARAMETER :: textwid = 60
  ! REAL(DP), PARAMETER :: tol = 1.E-10_DP
  LOGICAL, PARAMETER :: compress = .TRUE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Input variables
  INTEGER :: i, j, x, y, Npar, Nparfree, Nwfree, NiniMC
  INTEGER :: Ncont, Nband, Nline, Nextc, Nstar, Nextra
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: maskint ! convert mask=0 to mask=T
  LOGICAL :: verbose, calib, newseed, newinit, dostop
  CHARACTER(lenpath) :: dirIN, dirOUT, filOBS, filOUT, fiLOG
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lenpath), DIMENSION(:), ALLOCATABLE :: path1d
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB, labL

  !! Monte-Carlo parameters
  INTEGER :: ibest
  INTEGER, DIMENSION(:), ALLOCATABLE :: statusMC, NiterMC
  REAL(DP), DIMENSION(:), ALLOCATABLE :: chi2redMC
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: parerrMC
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: covarMC
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: parini
  
  !! chi2min_LM parameters
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: status
  INTEGER, DIMENSION(:), ALLOCATABLE :: itied
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: limits!, medSovN
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: chi2red
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: parerr, par
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: covarOBS, covpar
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: limited, maskS
  
  !! Output variables
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Niter
  INTEGER, DIMENSION(2) :: unitlog
  TYPE(time_type) :: timestr

  !! Analysis variables
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT, &
    FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FnuCONT_tab, &
    FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab

  
  
  !!------------------------------------------------------------------------
  !!                            I. Read the inputs
  !!------------------------------------------------------------------------

  CALL READ_HDF5(STRARR1D=path1d, FILE='../out/set_input'//h5ext, NAME='input dir')
  dirIN = path1d(1)
  filOBS = TRIMLR(dirIN)//'observation_MIR'//h5ext
  
  !! a. Main settings
  !!------------------
  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='wavelength (microns)', N1=NwOBS)
  CALL READ_MASTER(WAVALL=wOBS(:), DIRIN=dirIN, DIROUT=dirOUT, &
                   VERBOSE=verbose, NINIMC=NiniMC, &
                   CALIB=calib, NEWSEED=newseed, NEWINIT=newinit, &
                   LABL=labL, LABB=labB, QABS=Qabs, EXTINCT=extinct, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, &
                   NEXTC=Nextc, NSTAR=Nstar, NEXTRA=Nextra, DOSTOP=dostop, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, SPEC_UNIT=spec_unit)

  !! Output settings
  filOUT = TRIMLR(dirOUT)//'fit_chi2'//h5ext
  fiLOG = TRIMLR(dirOUT)//'log_fit_chi2'//ascext
  
  CALL INITIATE_CLOCK(timestr)
  OPEN (ulog,FILE=fiLOG,STATUS="REPLACE",ACTION="WRITE")
  unitlog(:) = [ ulog, ustd ]
  
  IF (newseed) CALL GENERATE_NEWSEED()
  IF (verbose) PRINT*
  DO i=1,MERGE(2,1,verbose)
    CALL BANNER_PROGRAM("LE MIROIR: LEast-squares fitting of Mid-IR emission " &
                        //"OptImized Routine", UNIT=unitlog(i), SWING=.TRUE.)
  END DO

  !! b. Observation sample
  !!-----------------------
  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, &
                 NAME='FnuOBS ('//TRIMLR(spec_unit)//')', N1=Nx, N2=Ny)
  CALL READ_HDF5(DBLARR3D=dFnuOBS, FILE=filOBS, &
                 NAME='dFnuOBS ('//TRIMLR(spec_unit)//')')

  !! For the input, the convention is mask=1 => blocks; mask=0 => passes.
  !! However, within Fortran, the convention is the Fortran mask=T => passes.
  CALL READ_HDF5(INTARR3D=maskint, FILE=filOBS, NAME='NaN mask')
  ALLOCATE(mask(Nx,Ny,NwOBS))
  mask(:,:,:) = ( maskint(:,:,:) == 0 )

  ALLOCATE(nuOBS(NwOBS))
  nuOBS(:) = MKS%clight/MKS%micron / wOBS(:)
  
  Nwfree = COUNT( ALL(ALL(.NOT. mask(:,:,:) &
                            .AND. dFnuOBS(:,:,:) > 0._DP,DIM=1),DIM=1) )
  IF (Nwfree > 0) &
    CALL IWHERE(ALL(ALL(.NOT. mask(:,:,:) &
                        .AND. dFnuOBS(:,:,:) > 0._DP,DIM=1),DIM=1), iwfree)

  !! Constraints
  ALLOCATE(itied(Npar))
  itied(:) = 0
  DO i=1,Npar
    IF (.NOT. TRIMEQ(parinfo(i)%tied,"")) &
      CALL IWHERE(TRIMEQ(parinfo(:)%name,parinfo(i)%tied),itied(i))
  END DO
  
  !! c. Compute the covariance matrix of the uncertainties
  !!-------------------------------------------------------
  !! Mask for the RMS covariance (diagonal)
  ALLOCATE(maskS(NwOBS,NwOBS))
  maskS(:,:) = .FALSE.
  FORALL (i=1:NwOBS) maskS(i,i) = .TRUE.

  !! 4D mask for the covariance matrix
  ! ALLOCATE(mask4D(Nx,Ny,NwOBS,NwOBS))
  ! mask4D(:,:,:,:) = .TRUE.
  ! FORALL (x=1:Nx,y=1:Ny,i=1:NwOBS, .NOT. mask(x,y,i))
  !   mask4D(x,y,i,:) = .FALSE.
  !   mask4D(x,y,:,i) = .FALSE.
  ! END FORALL

  !! Covariance matrix for each pixel and its inverse
  ALLOCATE(covarOBS(Nx,Ny,NwOBS,NwOBS), invLcovarOBS(Nx,Ny,NwOBS,NwOBS))

  covarOBS(:,:,:,:) = 0._DP
  !! invLcovarOBS is a lower triangle matrix L (Likelihood) results from 
  !! the Cholesky decomposition of the covariance matrix (Appendix C2, Galliano18a)
  invLcovarOBS(:,:,:,:) = 0._DP
  FORALL (x=1:Nx,y=1:Ny,ALL(mask(x,y,:)))
    covarOBS(x,y,:,:) = UNPACK(dFnuOBS(x,y,:)**2, maskS, FIELD=0._DP)
    invLcovarOBS(x,y,:,:) = UNPACK(1._DP/dFnuOBS(x,y,:), maskS, FIELD=0._DP)

  END FORALL

  !! No wave calib (between diff instr/module)
  ! DO x=1,Nx
  !   DO y=1,Ny
  !     DO i=1,NwOBS
  !       IF (.NOT. mask(x,y,i)) THEN
  !         invLcovarOBS(x,y,i,:) = 0._DP
  !         invLcovarOBS(x,y,:,i) = 0._DP

  !       END IF
  !     END DO
  !   END DO
  ! END DO

  ! DEALLOCATE(covarOBS)

  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'LE MIROIR: Read the inputs [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO
  
  !!------------------------------------------------------------------------
  !!                          II. Initial parameters
  !!------------------------------------------------------------------------
  
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*) "ESTIMATING THE INITIAL VALUES OF THE PARAMETERS (" &
                        //TRIMLR(TIMINFO(timestr))//")"
  END DO

  !! a. Rough spectrum estimators
  !!------------------------------

  !! b. Automatic init param estimates
  !!-----------------------------------
  ALLOCATE(parini(Nx, Ny, Npar, MAX(NiniMC,1)))
  parini(:,:,:,:) = 0._DP
  CALL INITPARAM(NiniMC, IND=ind, PAR=parini(:,:,:,:), PARINFO=parinfo(:), &
                 ITIED=itied(:), MASK=mask(:,:,:), &
                 NEWINIT=newinit, FILOBS=filOBS, &
                 LABB=labB(:), LABL=labL(:), QABS=Qabs(:))
  
  Nparfree = COUNT((.NOT. parinfo(:)%fixed) .AND. (itied(:) <= 0))

  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*) 'Number of free param: '//TRIMLR(PRING(Nparfree))
    
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'LE MIROIR: Initial parameters [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO

  !!------------------------------------------------------------------------
  !!                            Run the fitter
  !!------------------------------------------------------------------------
  
  ALLOCATE(par(Nx,Ny,Npar), status(Nx,Ny), resid(NwOBS), chi2red(Nx,Ny), &
           parerr(Nx,Ny,Npar), covpar(Nx,Ny,Npar,Npar), Niter(Nx,Ny), &
           Fnu_mod(NwOBS), limits(Npar,2), limited(Npar,2))

  FORALL (i=1:Npar)
    limits(i,:) = parinfo(i)%limits(:)
    limited(i,:) = parinfo(i)%limited(:)
  END FORALL
  
  ALLOCATE (statusMC(MAX(NiniMC,1)),NiterMC(MAX(NiniMC,1)), &
            parerrMC(Npar,MAX(NiniMC,1)),chi2redMC(MAX(NiniMC,1)), &
            covarMC(Npar,Npar,MAX(NiniMC,1)))

  par(:,:,:) = parini(:,:,:,1)
  CALL WRITE_HDF5(STRARR1D=parinfo(:)%name, NAME='Parameter label', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)
  CALL WRITE_HDF5(INITDBLARR=[Nx,Ny,Npar], NAME='Best fitted parameter value', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(INITDBLARR=[Nx,Ny,Npar], NAME='Best fitted parameter error', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

  chi2fitx: DO xOBS=1,Nx
    chi2fity: DO yOBS=1,Ny
      notmasked: IF (ANY( mask(xOBS,yOBS,:) )) THEN
        DO i=1,MERGE(2,1,verbose)
          WRITE(unitlog(i),*) 'pos: ('//TRIMLR(PRING(xOBS))//', '//TRIMLR(PRING(yOBS))//'): '
          WRITE(unitlog(i),*)
        END DO

        !! a. Levenberg-Marquardt method
        !!-------------------------------
        pariniMC: IF (NiniMC == 0) THEN

          CALL CHI2MIN_LM (residuals, NwOBS, PAR=par(xOBS,yOBS,:), VERBOSE=debug, &
                           STATUS=status(xOBS,yOBS), PARNAME=parinfo(:)%name, &
                           LIMITED=limited(:,:), LIMITS=limits(:,:), &
                           FIXED=parinfo(:)%fixed, ITIED=itied(:), &
                           CHI2RED=chi2red(xOBS,yOBS), NITER=Niter(xOBS,yOBS), &
                           PARERR=parerr(xOBS,yOBS,:), COVAR=covpar(xOBS,yOBS,:,:))

          DO i=1,MERGE(2,1,verbose)
            WRITE(unitlog(i),*) "PROGRAM EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
            WRITE(unitlog(i),*)
          END DO
          
        ELSE

          !! Vary the init param and keep the best fit
          DO j=1,NiniMC
            DO i=1,MERGE(2,1,verbose)
              WRITE(unitlog(i),*) 'iniMC iteration: '//TRIMLR(PRING(j))//'/'//TRIMLR(PRING(NiniMC))
            END DO
            
            CALL CHI2MIN_LM (residuals, NwOBS, PAR=parini(xOBS,yOBS,:,j), VERBOSE=debug, &
                             STATUS=statusMC(j), PARNAME=parinfo(:)%name, &
                             LIMITED=limited(:,:), LIMITS=limits(:,:), &
                             FIXED=parinfo(:)%fixed, ITIED=itied(:), &
                             CHI2RED=chi2redMC(j), NITER=NiterMC(j), &
                             PARERR=parerrMC(:,j), COVAR=covarMC(:,:,j))

            DO i=1,MERGE(2,1,verbose)
              WRITE(unitlog(i),*) "PROGRAM EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
              WRITE(unitlog(i),*)
            END DO
          
          END DO
          ibest = MAXVAL(MINLOC(chi2redMC(:)))
          par(xOBS,yOBS,:) = parini(xOBS,yOBS,:,ibest)
          parerr(xOBS,yOBS,:) = parerrMC(:,ibest)
          covpar(xOBS,yOBS,:,:) = covarMC(:,:,ibest)
          status(xOBS,yOBS) = statusMC(ibest)
          Niter(xOBS,yOBS) = NiterMC(ibest)
          chi2red(xOBS,yOBS) = chi2redMC(ibest)
          
          DO i=1,MERGE(2,1,verbose)
            WRITE(unitlog(i),*) "  MC chi2red is between " &
              //TRIMLR(PRING(MINVAL(chi2redMC(:)),NDEC=6))//" and " &
              //TRIMLR(PRING(MAXVAL(chi2redMC(:)),NDEC=6))
          
          END DO

        END IF pariniMC

        !! b. Derived quantities
        !!-----------------------

        !! c. Correlation between parameters
        !!-----------------------------------

        !! d. Summary
        !!------------
        DO i=1,MERGE(2,1,verbose)
          WRITE(unitlog(i),*) '>>>'
          WRITE(unitlog(i),*) "Checking output ("//TRIMLR(PRING(xOBS))//', '//TRIMLR(PRING(yOBS))//'): '
          WRITE(unitlog(i),*) "  status = "//TRIMLR(PRING(status(xOBS,yOBS)))
          WRITE(unitlog(i),*) "  chi2red = "//TRIMLR(PRING(chi2red(xOBS,yOBS),NDEC=10))
          WRITE(unitlog(i),*) "  Niter = "//TRIMLR(PRING(Niter(xOBS,yOBS)))
          WRITE(unitlog(i),*) '<<<'
        END DO

        CALL WRITE_HDF5(DBLARR3D=par(xOBS:xOBS,yOBS:yOBS,:), &
                        NAME='Best fitted parameter value', &
                        FILE=filOUT, COMPRESS=compress, VERBOSE=debug, &
                        IND1=[xOBS,xOBS], IND2=[yOBS,yOBS])
        CALL WRITE_HDF5(DBLARR3D=parerr(xOBS:xOBS,yOBS:yOBS,:), &
                        NAME='Best fitted parameter error', &
                        FILE=filOUT, COMPRESS=compress, VERBOSE=debug, &
                        IND1=[xOBS,xOBS], IND2=[yOBS,yOBS])

      END IF notmasked
    END DO chi2fity
  END DO chi2fitx

  !! Extending init param for HB inputs (alternatively by Python)
  ! CALL WRITE_HDF5(DBLARR3D=par(:,:,:), NAME='Chi2 fitted parameter value', &
  !                 FILE=filOBS, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)

  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'LE MIROIR: Chi2 fit [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO
  
  !!------------------------------------------------------------------------
  !!          Quick analysis (full post-processing see analysis.f90)
  !!------------------------------------------------------------------------

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
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR4D=FnuLINE_tab, NAME='FnuLINE ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR4D=FnuBAND_tab, NAME='FnuBAND ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR4D=FnuSTAR_tab, NAME='FnuSTAR ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR4D=Pabs_tab, NAME='PABS', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=FnuMOD, NAME='FnuMOD ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) 'LE MIROIR: fit_chi2 analysis [done]'
    WRITE(unitlog(i),*) REPEAT('=',textwid)//NEW_LINE('')
  END DO
  
  !! Final
  !!-------
  DO i=1,MERGE(2,1,verbose)
    WRITE(unitlog(i),*)
    WRITE(unitlog(i),*) "PROGRAM EXECUTED IN "//TRIMLR(TIMINFO(timestr))//"."
    WRITE(unitlog(i),*)
  END DO
  
  !! Free memory space
  DEALLOCATE(wOBS, FnuOBS, dFnuOBS, maskint, mask)

END PROGRAM fit_chi2
