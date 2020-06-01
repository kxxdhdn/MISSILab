MODULE miles_external

  ! USE auxil
  ! USE factable
  USE utilities, ONLY: DP
  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PUBLIC :: NwOBS, Nx, Ny, xOBS, yOBS
  INTEGER, SAVE, PUBLIC :: NAv, Nbb, Nline, Nband, Nstar
  CHARACTER(30), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: label
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC      :: parm
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC      :: wavOBS, nuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC      :: Lnutot
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC  :: LnuOBS, dLnuOBS
  REAl(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC      :: Pabs, Fnu_cont
  REAl(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC      :: Fnu_line, Fnu_band
  REAl(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC      :: Fnu_star, Fnu
  !! A sr * sr2arcsec2 = B arcsec^2
  ! REAL(DP), PARAMETER, PUBLIC :: sr2arcsec2 = (180._DP/pi)**2
  !! A deg * deg2arcsec = B arcsec
  ! REAL(DP), PARAMETER, PUBLIC :: deg2arcsec = 3600._DP

  PUBLIC :: residuals

CONTAINS

  !! Main function for L-M
  !!-----------------------
  FUNCTION residuals(par, NwOBS0)
    
    USE auxil, ONLY: specMODEL
    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN)                :: NwOBS0
    REAL(DP), DIMENSION(:), INTENT(IN) :: par
    REAL(DP), DIMENSION(NwOBS0)        :: residuals

    Lnutot = specModel(label, par, &
             NAv, Nbb, Nline, Nband, Nstar, &
             wavOBS)

    residuals = LnuOBS(xOBS,yOBS,:) - Lnutot

  END FUNCTION residuals
  
END MODULE miles_external

  !!-------------------------------------------------------

PROGRAM miles
  
  USE auxil
  USE utilities, ONLY: DP
  USE constants, ONLY: MKS
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext
  USE interpolation, ONLY: interp_lin_sorted
  USE chi2_minimization, ONLY: chi2min_LM
  USE miles_external, ONLY: NwOBS, Nx, Ny, xOBS, yOBS, &
                            wavOBS, nuOBS, LnuOBS, dLnuOBS, &
                            label, parm, NAv, Nbb, Nline, Nband, Nstar, &
                            Pabs, Fnu_cont, Fnu_line, Fnu_band, Fnu_star, &
                            Fnu, Lnutot, residuals
  IMPLICIT NONE

  INTEGER :: i, Nparm
  ! REAL(DP), DIMENSION(:), ALLOCATABLE :: wave, nu
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: chi2red
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: parerr, par
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: covpar
  
  !!-------------------------------------------------------
  !!                         INUPTS
  !!-------------------------------------------------------

  !! Banners
  !!---------
  !! Info
  !! Manual

  !! Input properties
  !!------------------
  !! Observation
  CALL read_hdf5(DBLARR1D=wavOBS, NAME="Wavelength (micron)", &
                 N1=NwOBS, FILE="input/test_sample"//h5ext)
  CALL read_hdf5(DBLARR3D=LnuOBS, NAME="LnuOBS (MJy/sr)", &
                 N1=Nx, N2=Ny, FILE="input/test_sample"//h5ext)
  CALL read_hdf5(DBLARR3D=dLnuOBS, NAME="dLnuOBS (MJy/sr)", &
                 FILE="input/test_sample"//h5ext)
  nuOBS = MKS%clight/MKS%micron / wavOBS
  print*, '>> read_DATA done <<'
  print*, 'Nx, Ny = ', Nx, Ny

  ! wavinterp = .True.
  ! IF (.NOT. PRESENT(nu)) THEN
  !   Nw = SIZE(nu0)
  !   wavinterp = .False.
  ! ELSE
  !   Nw = SIZE(nu)

  ! END IF
  ! IF (wavinterp.EQ..True.) THEN
  !   Fnu = interp_lin_sorted(nu0, Fnu, nu) ! ???

  ! END IF
  print*, '>> wave_INTERP undo <<'

  !!-------------------------------------------------------
  !!       Spectral Feature properties and Templates
  !!-------------------------------------------------------
  
  !!-------------------------------------------------------
  !!        Initial Guesses and Parameter Variations
  !!-------------------------------------------------------
  CALL chi2_INIT(label, parm, NAv, Nbb, Nline, Nband, Nstar)
  print*, '>> chi2_INIT done <<'
  !! Continuum (modifBB) parameters
  ALLOCATE(nQ_str(Nbb))
  CALL optics_INIT(label, nQ_str)
  print*, '>> optics_INIT done <<'
  !! Parameters
  !!------------
  Nparm = 2*Nbb + 3*Nline + 4*Nband + NAv + Nstar
  ALLOCATE(chi2red(Nx,Ny), par(Nx,Ny,Nparm), parerr(Nx,Ny,Nparm), &
           covpar(Nx,Ny,Nparm,Nparm))
  ! Output arrays

  ! Big loop
  !----------
  chi2fitx: DO xOBS=1,Nx
    chi2fity: DO yOBS=1,Ny
      CALL chi2min_LM(residuals, NwOBS,PAR=par(xOBS,yOBS,:), &
                      VERBOSE=.False., &
                      ! LIMITED=limited(:,:), LIMITS=limits(:,:), &
                      ! ITIED=itied(:), FIXED=parinfo(:)%fixed, &
                      ! PARNAME=parinfo(:)%name, STATUS=status(xOBS,yOBS), &
                      PARERR=parerr(xOBS,yOBS,:), &
                      COVAR=covpar(xOBS,yOBS,:,:), &
                      ! NITER=Niter(xOBS,yOBS), &
                      CHI2=chi2red(xOBS,yOBS))
                      ! STEP=step(:), RELSTEP=relstep(:), TWOSIDE=twoside(:))

  
  !!-------------------------------------------------------
  !!                       Complete Fit
  !!-------------------------------------------------------
    END DO chi2fity

  END DO chi2fitx
  print*, '>> chi2min_LM done <<'

  
  CALL spec_MODEL(label, parm, &
                  NAv, Nbb, Nline, Nband, Nstar, &
                  Pabs, Fnu_cont, Fnu_line, Fnu_band, &
                  Fnu_star, Fnu, wavOBS)
  print*, '--------------- FIN ---------------'

  !!-------------------------------------------------------
  !!                         Outputs
  !!-------------------------------------------------------
  CALL write_hdf5(wavOBS, NAME="Wavelength (micron)", &
                  COMPRESS=.False., APPEND=.False., &
                  FILE="output/chi2fit"//h5ext)
  CALL write_hdf5(Fnu_cont*Pabs, NAME="FnuBB (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/chi2fit"//h5ext)
  CALL write_hdf5(Fnu_line+(Fnu_cont+Fnu_star)*Pabs, NAME="FnuLINE (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/chi2fit"//h5ext)
  CALL write_hdf5((Fnu_band+Fnu_cont+Fnu_star)*Pabs, NAME="FnuBAND (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/chi2fit"//h5ext)
  CALL write_hdf5(Fnu_star*Pabs, NAME="FnuSTAR (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/chi2fit"//h5ext)
  CALL write_hdf5(Fnu, NAME="Fnu (MJy/sr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE="output/chi2fit"//h5ext)
  print*, 'write_RESULTS done'

END PROGRAM miles
