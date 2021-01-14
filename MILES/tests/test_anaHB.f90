PROGRAM test_anaHB

  USE auxil, ONLY: parinfo_type, indpar_type, Qabs_type, read_master, specModel
  USE utilities, ONLY: DP, pring, trimLR, trimeq, timinfo, &
                       banner_program, ustd, isNaN
  USE arrays, ONLY: iwhere, closest
  USE statistics, ONLY: MEDIAN
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, lenpar
  IMPLICIT NONE

  !! Parameters
  CHARACTER(*), PARAMETER :: filOBS = './dat/observations_fitMIR'//h5ext
  CHARACTER(*), PARAMETER :: dirOUT = './out/'
  CHARACTER(*), PARAMETER :: filOUT = dirOUT//'test_fitHB'//h5ext
  LOGICAL, PARAMETER :: compress = .FALSE.
  LOGICAL, PARAMETER :: debug = .FALSE.
  
  !! Input variables
  INTEGER :: i, Nx, Ny, NwOBS
  INTEGER :: Ncont, Nband, Nline, Npabs, Nstar
  CHARACTER(lenpar) :: spec_unit
  LOGICAL :: verbose
  ! LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: mask
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wOBS
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuOBS, dFnuOBS
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: par
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: FnuCONT, &
    FnuBAND, FnuSTAR, Pabs, FnuLINE, FnuMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: FnuCONT_tab, &
    FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab
  TYPE(indpar_type) :: ind
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE :: Qabs
  
  !!------------------------------------------------------------------------
  !!                            I. Read the inputs
  !!------------------------------------------------------------------------
  CALL READ_HDF5(wOBS, FILE=filOBS, NAME='Wavelength (microns)', N1=NwOBS)
  CALL READ_MASTER(WAVALL=wOBS(:), VERBOSE=verbose, QABS=Qabs, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, NPABS=Npabs, NSTAR=Nstar, &
                   INDPAR=ind, SPEC_UNIT=spec_unit)
  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE=filOBS, &
                 NAME='FnuOBS ('//TRIMLR(spec_unit)//')', N1=Nx, N2=Ny)
  CALL READ_HDF5(DBLARR3D=dFnuOBS, FILE=filOBS, &
                 NAME='dFnuOBS ('//TRIMLR(spec_unit)//')')
  CALL READ_HDF5(DBLARR3D=par, FILE=filOUT, &
                 NAME='Best fitted parameter value')

  CALL WRITE_HDF5(DBLARR3D=par, NAME='Best fitted parameter value', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.FALSE.)

  print*, 'Read fitted parameters [done]'//NEW_LINE('')

  print*, '=================================================='
  
  !!------------------------------------------------------------------------
  !!                                 Analysis
  !!------------------------------------------------------------------------

  !! Calculate model
  !!-----------------
  FnuMOD = specModel( wOBS(:), INDPAR=ind, PARVAL=par(:,:,:), QABS=Qabs(:), &
                      FNUCONT=FnuCONT, FNUBAND=FnuBAND, FNUSTAR=FnuSTAR, &
                      PABS=Pabs, FNULINE=FnuLINE, &
                      FNUCONT_TAB=FnuCONT_tab, FNUBAND_TAB=FnuBAND_tab, &
                      FNUSTAR_TAB=FnuSTAR_tab, PABS_TAB=Pabs_tab, &
                      FNULINE_TAB=FnuLINE_tab )

  DO i=1,Npabs
    FnuCONT_tab(:,:,:,i) = FnuCONT_tab(:,:,:,i) * Pabs(:,:,:)
  END DO 
  CALL WRITE_HDF5(DBLARR4D=FnuCONT_tab, NAME='FnuCONT ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  DO i=1,Nline
    FnuLINE_tab(:,:,:,i) = FnuLINE_tab(:,:,:,i) + (FnuCONT(:,:,:)+FnuSTAR(:,:,:))*Pabs(:,:,:)
  END DO
  CALL WRITE_HDF5(DBLARR4D=FnuLINE_tab, NAME='FnuLINE ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  DO i=1,Nband
    FnuBAND_tab(:,:,:,i) = FnuBAND_tab(:,:,:,i) + (FnuCONT(:,:,:)+FnuSTAR(:,:,:))*Pabs(:,:,:)
  END DO
  CALL WRITE_HDF5(DBLARR4D=FnuBAND_tab, NAME='FnuBAND ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  DO i=1,Nstar
    FnuSTAR_tab(:,:,:,i) = FnuSTAR_tab(:,:,:,i) * Pabs(:,:,:)
  END DO
  CALL WRITE_HDF5(DBLARR4D=FnuSTAR_tab, NAME='FnuSTAR ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  CALL WRITE_HDF5(DBLARR3D=FnuMOD, NAME='FnuMOD ('//TRIMLR(spec_unit)//')', &
                  FILE=filOUT, COMPRESS=compress, VERBOSE=debug, APPEND=.TRUE.)
  
  print*, 'fitHB analysis [done]'//NEW_LINE('')

  print*, '=================================================='


  !! Free memory space
  DEALLOCATE(par, wOBS, FnuOBS, dFnuOBS, FnuMOD, &
             FnuCONT, FnuBAND, FnuSTAR, Pabs, FnuLINE, &
             FnuCONT_tab, FnuBAND_tab, FnuSTAR_tab, Pabs_tab, FnuLINE_tab)

END PROGRAM test_anaHB
