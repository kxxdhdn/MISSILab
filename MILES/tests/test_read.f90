PROGRAM test_read

  USE auxil, ONLY: read_master, Qabs_type, parinfo_type, indpar_type
  USE utilities, ONLY: DP
  USE arrays, ONLY: ramp
  USE inout, ONLY: lenpar
  USE grain_optics, ONLY: lendustQ
  IMPLICIT NONE

  INTEGER, PARAMETER :: Ngrid = 1000
  INTEGER :: Nmcmc, NiniMC, Npar
  INTEGER :: Ncont, Nband, Nline
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wave
  CHARACTER(lenpar) :: spec_unit
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labB, labL
  CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ
  LOGICAL :: verbose, calib, newseed, newinit, dostop
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE :: Qabs
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE :: parinfo
  TYPE(indpar_type) :: ind

  wave = RAMP(Ngrid, 1._DP, 40._DP, XLOG=.TRUE.)

  CALL READ_MASTER(WAVALL=wave(:), &
                   NMCMC=Nmcmc, VERBOSE=verbose, NiniMC=NiniMC, &
                   CALIB=calib, NEWSEED=newseed, NEWINIT=newinit, &
                   LABQ=labQ, LABL=labL, LABB=labB, QABS=Qabs, &
                   NCONT=Ncont, NBAND=Nband, NLINE=Nline, DOSTOP=dostop, &
                   PARINFO=parinfo, INDPAR=ind, NPAR=Npar, &
                   SPEC_UNIT=spec_unit)

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
