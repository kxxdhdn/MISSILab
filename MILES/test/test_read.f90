PROGRAM test_read

  USE core, ONLY: read_master, Qabs_type, parinfo_type, indpar_type, lencorr
  USE utilities, ONLY: DP, trimLR
  USE arrays, ONLY: ramp
  USE inout, ONLY: lenpar, lenpath, read_hdf5, h5ext
  USE grain_optics, ONLY: lendustQ
  IMPLICIT NONE

  INTEGER, PARAMETER :: Ngrid = 1000
  INTEGER :: Nmcmc, NiniMC, Npar, Nparhyp, Ncorr, Ncorrhyp, indres0
  INTEGER :: Ncont, Nband, Nline, Nextc, Nstar, Nextra
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wave
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: extinct
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB, labL
  CHARACTER(lenpath) :: dirIN, dirOUT, filOBS
  CHARACTER(lenpath), DIMENSION(:), ALLOCATABLE :: path1d
  CHARACTER(lencorr), DIMENSION(:), ALLOCATABLE :: corrhypname
  LOGICAL :: verbose, calib, newseed, newinit, dostop, resume, nohi
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE :: Qabs
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE :: parinfo, parhypinfo
  TYPE(indpar_type) :: ind

  wave = RAMP(Ngrid, 1._DP, 40._DP, XLOG=.TRUE.)

  CALL READ_HDF5(STRARR1D=path1d, FILE='./dat/set_input'//h5ext, NAME='input dir')
  dirIN = path1d(1)
  filOBS = TRIMLR(dirIN)//'observation_MIR'//h5ext
  CALL READ_MASTER(WAVALL=wave(:), DIRIN=dirIN, DIROUT=dirOUT, &
                   VERBOSE=verbose, NMCMC=Nmcmc, NINIMC=NiniMC, &
                   CALIB=calib, NEWSEED=newseed, NEWINIT=newinit, &
                   LABL=labL, LABB=labB, QABS=Qabs, EXTINCT=extinct, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, &
                   NEXTC=Nextc, NSTAR=Nstar, NEXTRA=Nextra, NCORR=Ncorr, &
                   DOSTOP=dostop, RESUME=resume, INDRESUME=indres0, NOHI=nohi, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, &
                   CORRHYPNAME=corrhypname, SPEC_UNIT=spec_unit, &
                   PARHYPINFO=parhypinfo, NPARHYP=Nparhyp, NCORRHYP=Ncorrhyp)

  !! Print outputs
  !!---------------
  PRINT*, 'Nmcmc = ', Nmcmc
  PRINT*, 'verbose = ', verbose
  ! PRINT*, 'Qabs%wave = ', Qabs(1)%wave
  ! PRINT*, 'Qabs%rho = ', Qabs(1)%rho
  ! PRINT*, 'Qabs%kappa = ', Qabs(2)%kappa
  PRINT*, 'parinfo%name = ', parinfo%name
  PRINT*, 'parinfo%fixed = ', parinfo%fixed
  PRINT*, 'parinfo(1)%limited = ', parinfo(1)%limited

END PROGRAM test_read
