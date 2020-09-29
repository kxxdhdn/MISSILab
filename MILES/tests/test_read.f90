PROGRAM test_read

  USE auxil, ONLY: read_master, Qabs_type
  USE utilities, ONLY: DP
  USE arrays, ONLY: ramp
  USE inout, ONLY: lenpar
  USE grain_optics, ONLY: lendustQ
  IMPLICIT NONE

  INTEGER, PARAMETER :: Ngrid = 1000
  INTEGER :: Nband, Nline
  INTEGER, DIMENSION(:), ALLOCATABLE :: indBIN, indLIN
  REAL(DP), DIMENSION(:), ALLOCATABLE :: wave
  CHARACTER(lenpar), DIMENSION(:), ALLOCATABLE :: labBIN, labLIN
  CHARACTER(lendustQ), DIMENSION(:), ALLOCATABLE :: labQ
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE :: Qabs

  wave = RAMP(Ngrid, 1._DP, 40._DP, XLOG=.TRUE.)

  ! ALLOCATE(labQ(2), Qabs(2))
  labQ = (/'Sil_D03 ', 'ACAR_Z96'/)
  labBIN = (/'Main 6.2 (1) ', & ! 7
             'Main 6.2 (2) ', & ! 8
             'Plateau 7.7  ', & ! 12
             'Main 7.7 (1) ', & ! 13
             'Main 7.7 (2) ', & ! 14
             'Main 8.6     ', & ! 16
             'Small 11.0   ', & ! 19
             'Main 11.2    ', & ! 20
             'Plateau 11.3 ', & !21
             'Main 12.7 (1)', & ! 24
             'Main 12.7 (2)', & ! 25
             'Plateau 17.0 '/) ! 30
  labLIN = (/'H2S7  ', & ! 3
             'H2S5  ', & ! 9
             'ArII  ', & ! 11
             'ArIII1', & ! 20
             'H2S3  ', & ! 23
             'SIV   ', & ! 26
             'H2S2  ', & ! 27
             'NeII  ', & ! 29
             'NeIII1', & ! 33
             'H2S1  ', & ! 34
             'SIII1 '/) ! 36
  Nband = SIZE(labBIN)
  Nline = SIZE(labLIN)

  CALL READ_MASTER(LABQ=labQ(:), LABBIN=labBIN(:), LABLIN=labLIN(:), WAVEALL=wave(:), &
                   QABS=Qabs, INDBIN=indBIN, INDLIN=indLIN)

  !! Print outputs
  !!---------------
  ! PRINT*, 'Qabs1%wave = ', Qabs(1)%wave
  ! PRINT*, 'Qabs1%rho = ', Qabs(1)%rho
  ! PRINT*, 'Qabs1%Qova = ', Qabs(1)%Qova
  ! PRINT*, 'Qabs2%coeffMBB = ', Qabs(2)%coeffMBB
  PRINT*, 'Band indices: ', indBIN
  PRINT*, 'Line indices: ', indLIN

END PROGRAM test_read
