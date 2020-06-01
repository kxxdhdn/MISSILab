MODULE modifBB_external

  USE utilities, ONLY: DP
  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PUBLIC :: Npar, NwOBS, Nx, Ny, xOBS, yOBS
  CHARACTER(30), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: label
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC      :: wOBS, nuOBS
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC  :: LnuOBS, dLnuOBS

  PUBLIC :: funcResid

CONTAINS

  !! Main function for L-M
  !!-----------------------
  FUNCTION funcesid(par, NwOBS0)
    
    USE auxil, ONLY: Qabs_str, modifBB
    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER i
    INTEGER, INTENT(IN)                      :: NwOBS0
    REAL(DP), DIMENSION(:), INTENT(IN)       :: par
    TYPE(Qabs_str), DIMENSION(:), INTENT(IN) :: nQabs
    REAL(DP), DIMENSION(Nx,Ny,NwOBS)         :: LnuMOD
    REAL(DP)                                 :: M, T
    REAL(DP), DIMENSION(NwOBS0)              :: funcResid

    DO i=2,Npar,2
      M = par(i-1)
      T = par(i)
      LnuMOD(i,:) = M * modifBB(T, nQabs(i))

    END DO
    
    WHERE (mask(xOBS,yOBS,:))
      funcResid = LnuOBS(xOBS,yOBS,:) - SUM(LnuMOD(i,:), DIM=1)
    ELSE WHERE
      funcResid(:) = 0._DP

    END WHERE

  END FUNCTION funcResid
  
END MODULE modifBB_external

!!==========================================================================

PROGRAM test_modifBB

  USE auxil, ONLY: Qabs_str, make_Qabs 
  USE utilities, ONLY: DP, TRIMLR
  USE inout, ONLY: write_hdf5, h5ext, write_ascii, ascext
  USE chi2_minimization, ONLY: chi2min_LM
  USE modifBB_external, ONLY: funcResid, Npar, NwOBS, wOBS, nuOBS, LnuOBS, dLnuOBS
  IMPLICIT NONE

  INTERGER :: i, status, Nparfree
  INTEGER, DIMENSION(Npar)                  :: itied
  REAL(DP), PARAMETER                       :: tol = 1.E-10_DP
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: chi2red
  REAL(DP), DIMENSION(Npar)                 :: parerr
  REAL(DP), DIMENSION(Npar,Npar)            :: covar
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: LnuMOD, residuals
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: par
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: limits
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: fixedINT, itiedINT
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: limitedINT
  LOGICAL, DIMENSION(Npar)                  :: fixed
  LOGICAL, DIMENSION(Npar,2)                :: limited
  CHARACTER(20), DIMENSION(Npar)            :: parname

  CHARACTER(30)                             :: path
  CHARACTER(30), DIMENSION(:), ALLOCATABLE  :: label
  TYPE(Qabs_str), DIMENSION(:), ALLOCATABLE :: nQabs

  !! INPUTS
  CALL READ_HDF5(LnuOBS, FILE='dat/M83'//h5ext, NAME='LnuOBS (MJyovsr)')
  CALL READ_HDF5(wOBS, FILE='dat/M83'//h5ext, NAME='Wavelength (microns)')
  CALL READ_HDF5(dLnuOBS, FILE='dat/M83'//h5ext, NAME='dLnuOBS (MJyovsr)')
  print*, 'import obs spectra [done]'

  ALLOCATE(label(n), nQabs(n), M(n), T(n))
  label = (/'Sil_D03 ', 'ACAR_Z96'/)
  M = 3._DP
  T = 100._DP

  PRINT*, "Number of components: ", Nbb

  !! Get par
  CALL make_Qabs(label, nQabs)
  parname(:) = (/'M1', 'T1', 'M2', 'T2'/)

  !! Run the fitter
  !!----------------
  ALLOCATE (LnuMOD(NwOBS), residuals(NwOBS))
  CALL chi2min_LM(funcResid, NwOBS, par, tol, residuals, status, .TRUE., &
                  limited, limits, fixed, itied, parname, &
                  PARERR=parerr, COVAR=covar)
  !! Calculate Fnu
  LnuOBS(:) = 0._DP
  dLnuOBS(:) = 1._DP
  LnuMOD(:) = funcResid(par, NwOBS)
  path = 'out/test_modifBB'
  PRINT*
  PRINT*, "Checking output:"
  PRINT*, "  status = "//TRIMLR(PRING(Status))
  chi2red = SUM(deviate(:)**2) / (Nobs-Nparfree)
  PRINT*, "  chi2 = "//TRIMLR(PRING(chi2red,NDEC=10))
  CALL WRITE_ASCII(VEC1=par,VEC2=parerr,FILE="out/chi2min.txt")
  PRINT*
  PRINT*, "Covariance matrix:"
  DO i=1,Npar
    PRINT*, REAL(covar(:,i),KIND(0.))
  END DO
  PRINT*
  CALL write_hdf5(nQabs(1)%wave, FILE=TRIMLR(path)//h5ext, &
                  COMPRESS=.False., NAME="Wavelength (micron)")
  CALL write_hdf5(nQabs(1)%nu, FILE=TRIMLR(path)//h5ext, &
                  COMPRESS=.False., NAME="Wavelength (Hz)", &
                  APPEND=.True.)
  CALL write_hdf5(Fnu, FILE=TRIMLR(path)//h5ext, &
                  COMPRESS=.False., NAME="Fnu (MJy/sr)", &
                  APPEND=.True.)
  CALL write_ascii(VEC1=nQabs(1)%wave, VEC2=Fnu, &
                   FILE=TRIMLR(path)//ascext)

  !! Free memory space
  DEALLOCATE(label, M, T, Fnu)

END PROGRAM test_modifBB
