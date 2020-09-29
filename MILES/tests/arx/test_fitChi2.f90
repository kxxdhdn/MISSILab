MODULE fitChi2_external

  USE auxil, ONLY: parinfo_type, Qabs_type
  USE utilities, ONLY: DP
  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PUBLIC :: Nx, Ny, NwOBS, xOBS, yOBS, jw
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC        :: ifreefilt
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC       :: wOBS, nuOBS
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC       :: FnuTOT, resid
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC   :: Pabs
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC   :: Fnu_cont
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC   :: Fnu_line
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC   :: Fnu_band
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC   :: Fnu_star
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC   :: FnuOBS, dFnuOBS
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC   :: FnuMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE, PUBLIC :: invLcovarOBS
  LOGICAL, SAVE, PUBLIC                                   :: iid
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE, PUBLIC    :: mask

  !! Init param
  !!------------
  INTEGER, SAVE, PUBLIC :: Nbb, Nline, Nband, NAv, Nstar, Npar
  TYPE(parinfo_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: parinfo
  CHARACTER(30), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC     :: compdust
  TYPE(Qabs_type), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC    :: nQabs
  
  PUBLIC :: funcResid

CONTAINS

  !! Main function for L-M
  !!-----------------------
  FUNCTION funcResid(parr0, NwOBS0)
    
    USE auxil, ONLY: specModel
    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN)                      :: NwOBS0
    REAL(DP), DIMENSION(:), INTENT(IN)       :: parr0
    INTEGER                                  :: i
    
    REAL(DP), DIMENSION(NwOBS0)              :: funcResid

    FnuTOT = specModel(wOBS, parr0, nQabs, Nbb, Nline, Nband, NAv, Nstar)

    !! Unweighted residuals    
    WHERE (mask(xOBS,yOBS,:))
      resid(:) = FnuOBS(xOBS,yOBS,:) - FnuTOT
    ELSEWHERE
      resid(:) = 0._DP

    END WHERE
    
    !! Weighted residuals
    ! funcResid(:) = MATMUL(invLcovarOBS(xOBS,yOBS,:,:), resid(:))
    
    !! With keyword iid
    IF (iid) THEN
      FORALL (i=1:NwOBS0) &
        funcResid(i) = invLcovarOBS(xOBS,yOBS,i,i) * resid(i)
    ELSE
      funcResid(:) = MATMUL(invLcovarOBS(xOBS,yOBS,:,:), resid(:))

    END IF

  END FUNCTION funcResid

END MODULE fitChi2_external


!!==========================================================================
!!                           Main Program
!!==========================================================================


PROGRAM test_fitChi2

  USE auxil, ONLY: make_par, make_Qabs, specModel, degradeRes
  USE datable, ONLY: LIN, BIN
  USE utilities, ONLY: DP, PRING, TRIMLR, isNaN
  USE arrays, ONLY: IWHERE, CLOSEST
  USE constants, ONLY: MKS
  USE statistics, ONLY: MEDIAN
  USE inout, ONLY: read_hdf5, write_hdf5, h5ext, write_ascii, ascext
  USE chi2_minimization, ONLY: chi2min_LM
  USE fitChi2_external, ONLY: funcResid, Nx, Ny, NwOBS, wOBS, nuOBS, &
                             FnuOBS, dFnuOBS, FnuMOD, resid, invLcovarOBS, &
                             Pabs, Fnu_cont, Fnu_line, Fnu_band, Fnu_star, &
                             Nbb, Nline, Nband, NAv, Nstar, Npar, parinfo, &
                             compdust, nQabs, xOBS, yOBS, mask, ifreefilt, iid
  IMPLICIT NONE

  INTEGER, DIMENSION(:,:,:), ALLOCATABLE    :: maskint
  REAL(DP), DIMENSION(:), ALLOCATABLE       :: parr0
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: parr

  INTEGER :: i, x, y, Niter, Nparfree, Nfreefilt
  INTEGER, DIMENSION(:,:), ALLOCATABLE      :: status
  REAL(DP), PARAMETER                       :: tol = 1.E-10_DP
  ! REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: medSovN
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: chi2red
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: parerr
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: limits
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: covarOBS, covpar
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: Fnu_cont_tab
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: Fnu_line_tab
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: Fnu_band_tab
  LOGICAL, DIMENSION(:,:), ALLOCATABLE      :: limited, maskS
  ! LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: mask4D

  CHARACTER(30)                             :: filOUT

  !! wavelengths indpdt odentically dist (diagonal covar matrix)
  iid = .True.
  
  !!-----------------
  !! Read the inputs
  !!-----------------
  CALL READ_HDF5(wOBS, FILE='dat/M83'//h5ext, NAME='Wavelength (microns)')
  CALL READ_HDF5(DBLARR3D=FnuOBS, FILE='dat/M83'//h5ext, NAME='FnuOBS (MJyovsr)')
  CALL READ_HDF5(DBLARR3D=dFnuOBS, FILE='dat/M83'//h5ext, NAME='dFnuOBS (MJyovsr)')
  CALL READ_HDF5(INTARR3D=maskint, FILE='dat/M83'//h5ext, NAME='Mask')

  Nx = SIZE(FnuOBS, 1)
  Ny = SIZE(FnuOBS, 2)    
  NwOBS = SIZE(wOBS)
  
  !! Fortran mask convention is mask=T => passes; mask=F => blocks
  ALLOCATE(mask(Nx,Ny,NwOBS))
  mask(:,:,:) = ( maskint(:,:,:) == 0 )
  
  WHERE (isNaN(FnuOBS(:,:,:))) FnuOBS(:,:,:) = 0._DP
  WHERE (isNaN(dFnuOBS(:,:,:))) dFnuOBS(:,:,:) = 0._DP
  
  Nfreefilt = COUNT(ALL(ALL(.NOT. mask(:,:,:) &
                            .AND. dFnuOBS(:,:,:) > 0._DP,DIM=1),DIM=1))
  IF (Nfreefilt > 0) &
    CALL IWHERE(ALL(ALL(.NOT. mask(:,:,:) &
                        .AND. dFnuOBS(:,:,:) > 0._DP,DIM=1),DIM=1),ifreefilt)

  !! Median signal-to-noise ratio
  ! ALLOCATE(medSovN(Nx,Ny))
  ! medSovN(:,:) = 0._DP
  ! DO x=1,Nx
  !   DO y=1,Ny
  !     IF (ANY(mask(x,y,:))) THEN
  !       medSovN(x,y) = MEDIAN(FnuOBS(x,y,:)/dFnuOBS(x,y,:),MASK=mask(x,y,:))
  !       IF (medSovN(x,y) < 0._DP) medSovN(x,y) = 0._DP

  !     END IF
  !   END DO
  ! END DO

  !! Read the instrumental covariance matrix (TBD)

  !! Compute the covariance matrix of the uncertainties
  !!----------------------------------------------------
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

  print*, 'Import obs spectra [done]'//NEW_LINE('')
  
  !!------------
  !! Init param
  !!------------
  Nbb=4
  Nline=11
  Nband=31
  NAv=1
  Nstar=1
  
  ALLOCATE(nuOBS(NwOBS), compdust(Nbb), nQabs(Nbb))
  
  nuOBS = MKS%clight/MKS%micron / wOBS
  
  !! Modified blackbody components
  compdust = [&
              'ACH2_Z96             ', &
              'Sil_D03              ', &
              'ACH2_Z96             ', &
              'ACH2_Z96             ' &
              ]
  CALL make_Qabs(compdust, nQabs, WAVEALL=wOBS)

  !! parr0 = [massBB [Msun/pc2], tempBB [K], &
  parr0 = [10._DP, 80._DP, 10._DP, 150._DP, 10._DP, 180._DP, 10._DP, 50._DP, &
  !!          Iline, Cline, Wline, &
              1._DP, LIN(3)%wave, degradeRes(LIN(3)%wave, .01_DP, 'SL-LL'), & 
              1._DP, LIN(9)%wave, degradeRes(LIN(9)%wave, .01_DP, 'SL-LL'), & 
              1._DP, LIN(11)%wave, degradeRes(LIN(11)%wave, .01_DP, 'SL-LL'), &
              1._DP, LIN(20)%wave, degradeRes(LIN(20)%wave, .01_DP, 'SL-LL'), &
              1._DP, LIN(23)%wave, degradeRes(LIN(23)%wave, .01_DP, 'SL-LL'), &
              ! 1._DP, LIN(24)%wave, degradeRes(LIN(24)%wave, .01_DP, 'SL-LL'), &
              1._DP, LIN(26)%wave, degradeRes(LIN(26)%wave, .01_DP, 'SL-LL'), &
              1._DP, LIN(27)%wave, degradeRes(LIN(27)%wave, .01_DP, 'SL-LL'), &
              1._DP, LIN(29)%wave, degradeRes(LIN(29)%wave, .01_DP, 'SL-LL'), &
              1._DP, LIN(33)%wave, degradeRes(LIN(33)%wave, .01_DP, 'SL-LL'), &
              1._DP, LIN(34)%wave, degradeRes(LIN(34)%wave, .01_DP, 'SL-LL'), &
              1._DP, LIN(36)%wave, degradeRes(LIN(36)%wave, .01_DP, 'SL-LL'), &
  !!          Iband, Cband, WSband, WLband, &
              ! 1._DP, BIN(1)%wave, BIN(1)%sigmaS, BIN(1)%sigmaL, & 
              ! 1._DP, BIN(2)%wave, BIN(2)%sigmaS, BIN(2)%sigmaL, & 
              1._DP, BIN(7)%wave, BIN(7)%sigmaS, BIN(7)%sigmaL, & 
              1._DP, BIN(8)%wave, BIN(8)%sigmaS, BIN(8)%sigmaL, & 
              1._DP, BIN(12)%wave, BIN(12)%sigmaS, BIN(12)%sigmaL, & 
              1._DP, BIN(13)%wave, BIN(13)%sigmaS, BIN(13)%sigmaL, & 
              1._DP, BIN(14)%wave, BIN(14)%sigmaS, BIN(14)%sigmaL, & 
              1._DP, BIN(16)%wave, BIN(16)%sigmaS, BIN(16)%sigmaL, & 
              1._DP, BIN(19)%wave, BIN(19)%sigmaS, BIN(19)%sigmaL, & 
              1._DP, BIN(20)%wave, BIN(20)%sigmaS, BIN(20)%sigmaL, &
              1._DP, BIN(21)%wave, BIN(21)%sigmaS, BIN(21)%sigmaL, & 
              1._DP, BIN(24)%wave, BIN(24)%sigmaS, BIN(24)%sigmaL, & 
              1._DP, BIN(25)%wave, BIN(25)%sigmaS, BIN(25)%sigmaL, & 
              1._DP, BIN(30)%wave, BIN(30)%sigmaS, BIN(30)%sigmaL, &

              1._DP, BIN(3)%wave, BIN(3)%sigmaS, BIN(3)%sigmaL, & 
              1._DP, BIN(4)%wave, BIN(4)%sigmaS, BIN(4)%sigmaL, & 
              1._DP, BIN(5)%wave, BIN(5)%sigmaS, BIN(5)%sigmaL, & 
              1._DP, BIN(6)%wave, BIN(6)%sigmaS, BIN(6)%sigmaL, & 
              1._DP, BIN(9)%wave, BIN(9)%sigmaS, BIN(9)%sigmaL, & 
              1._DP, BIN(10)%wave, BIN(10)%sigmaS, BIN(10)%sigmaL, & 
              1._DP, BIN(11)%wave, BIN(11)%sigmaS, BIN(11)%sigmaL, & 
              1._DP, BIN(15)%wave, BIN(15)%sigmaS, BIN(15)%sigmaL, & 
              1._DP, BIN(17)%wave, BIN(17)%sigmaS, BIN(17)%sigmaL, & 
              1._DP, BIN(18)%wave, BIN(18)%sigmaS, BIN(18)%sigmaL, & 
              1._DP, BIN(22)%wave, BIN(22)%sigmaS, BIN(22)%sigmaL, & 
              1._DP, BIN(23)%wave, BIN(23)%sigmaS, BIN(23)%sigmaL, & 
              1._DP, BIN(26)%wave, BIN(26)%sigmaS, BIN(26)%sigmaL, & 
              1._DP, BIN(27)%wave, BIN(27)%sigmaS, BIN(27)%sigmaL, & 
              1._DP, BIN(28)%wave, BIN(28)%sigmaS, BIN(28)%sigmaL, & 
              1._DP, BIN(29)%wave, BIN(29)%sigmaS, BIN(29)%sigmaL, & 
              1._DP, BIN(31)%wave, BIN(31)%sigmaS, BIN(31)%sigmaL, & 
              1._DP, BIN(32)%wave, BIN(32)%sigmaS, BIN(32)%sigmaL, & 
              1._DP, BIN(33)%wave, BIN(33)%sigmaS, BIN(33)%sigmaL, & 

  !!          Av [mag], &
              0._DP, & 
  !!          Fstar [Lsun/pc2]]
              1.E4_DP]
  
  !! Get Npar & parinfo
  CALL make_par(parr0, Nbb, Nline, Nband, NAv, Nstar, NPAR=Npar, PARINFO=parinfo)
  Nparfree = COUNT((.NOT. parinfo(:)%fixed) .AND. (parinfo(:)%itied <= 0))
  PRINT*, 'Number of free param: '//TRIMLR(PRING(Nparfree))

  !! The same parr for all pixels
  ALLOCATE(parr(Nx,Ny,Npar))

  FORALL (xOBS=1:Nx,yOBS=1:Ny) parr(xOBS,yOBS,:) = parr0(:)

  !! Multi dim type elements (limits & limited) not converted to multi dim real arrays
  ALLOCATE(limits(Npar,2), limited(Npar,2))
  
  FORALL (i=1:Npar)
    limits(i,:) = parinfo(i)%limits(:)
    limited(i,:) = parinfo(i)%limited(:)
    
  END FORALL
  
  !!----------------
  !! Run the fitter
  !!----------------
  ALLOCATE(resid(NwOBS),status(Nx,Ny), chi2red(Nx,Ny))
  ALLOCATE(parerr(Nx,Ny,Npar), covpar(Nx,Ny,Npar,Npar))

  ! chi2fitx: DO xOBS=6,7
    ! chi2fity: DO yOBS=6,7
  chi2fitx: DO xOBS=1,Nx
    chi2fity: DO yOBS=1,Ny
      notmasked: IF (ALL( mask(xOBS,yOBS,:) )) THEN
        PRINT*
        PRINT*, 'pos: ('//TRIMLR(PRING(xOBS))//', '//TRIMLR(PRING(yOBS))//'): '
        
        CALL chi2min_LM (funcResid, NwOBS, parr(xOBS,yOBS,:), tol, &
                         STATUS=status(xOBS,yOBS), VERBOSE=.False., &
                         LIMITED=limited, LIMITS=limits, &
                         FIXED=parinfo(:)%fixed, ITIED=parinfo(:)%itied, &
                         PARNAME=parinfo(:)%name, CHI2RED=chi2red(xOBS,yOBS), &
                         NITER=Niter, PARERR=parerr, COVAR=covpar)

        !! Compute model result
        PRINT*, '>>>'
        PRINT*, "Checking output ("//TRIMLR(PRING(xOBS))//', '//TRIMLR(PRING(yOBS))//'): '
        PRINT*, "  status = "//TRIMLR(PRING(status(xOBS,yOBS)))
        PRINT*, "  chi2 = "//TRIMLR(PRING(chi2red(xOBS,yOBS),NDEC=10))
        PRINT*, "  Niter = "//TRIMLR(PRING(Niter))
        ! CALL write_ascii(VEC1=parr(xOBS,yOBS,:), &
        !                  VEC2=parerr(xOBS,yOBS,:), &
        !                  FILE="out/chi2min_fitChi2.txt")
        PRINT*
        PRINT*, "Covariance matrix:"
        ! DO i=1,Npar
        !   PRINT*, REAL(covpar(xOBS,yOBS,:,i), KIND(0.))
        ! END DO
        PRINT*, '<<<'
        PRINT*
       
      END IF notmasked
    END DO chi2fity
  END DO chi2fitx

  ALLOCATE(Pabs(Nx,Ny,NwOBS), Fnu_cont(Nx,Ny,NwOBS), &
           Fnu_line(Nx,Ny,NwOBS), Fnu_band(Nx,Ny,NwOBS), Fnu_star(Nx,Ny,NwOBS))
  ALLOCATE(Fnu_cont_tab(Nx,Ny,Nbb,NwOBS), &
           Fnu_line_tab(Nx,Ny,Nline,NwOBS), Fnu_band_tab(Nx,Ny,Nband,NwOBS))

  FnuMOD = specModel(wOBS, parr, nQabs, Nbb, Nline, Nband, NAv, Nstar, &
                     PABS=Pabs, FNU_CONT=Fnu_cont, FNU_LINE=Fnu_line, &
                     FNU_BAND=Fnu_band, FNU_STAR=Fnu_star, &
                     FNU_CONT_TAB=Fnu_cont_tab, &
                     FNU_LINE_TAB=Fnu_line_tab, FNU_BAND_TAB=Fnu_band_tab)

  filOUT = 'out/test_fitChi2'
  CALL write_hdf5(wOBS, NAME="Wavelength (microns)", &
                  COMPRESS=.False., APPEND=.False., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(DBLARR3D=dFnuOBS, NAME="FnuUNC (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(DBLARR3D=FnuOBS, NAME="FnuOBS (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(DBLARR3D=Fnu_cont*Pabs, NAME="FnuBB (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(DBLARR3D=Fnu_line+(Fnu_cont+Fnu_star)*Pabs, NAME="FnuLINE (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(DBLARR3D=(Fnu_band+Fnu_cont+Fnu_star)*Pabs, NAME="FnuBAND (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(DBLARR3D=Fnu_star*Pabs, NAME="FnuSTAR (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  CALL write_hdf5(DBLARR3D=FnuMOD, NAME="FnuMOD (MJyovsr)", &
                  COMPRESS=.False., APPEND=.True., &
                  FILE=TRIMLR(filOUT)//h5ext)
  DO i=1,Nbb
    CALL write_hdf5(DBLARR3D=Fnu_cont_tab(:,:,i,:)*Pabs, &
                    NAME='FnuBB'//TRIMLR(PRING(i))//' (MJyovsr)', &
                    COMPRESS=.False., APPEND=.True., &
                    FILE=TRIMLR(filOUT)//h5ext)
  END DO
  DO i=1,Nline
    CALL write_hdf5(DBLARR3D=Fnu_line_tab(:,:,i,:)+(Fnu_cont+Fnu_star)*Pabs, &
                    NAME='FnuLINE'//TRIMLR(PRING(i))//' (MJyovsr)', &
                    COMPRESS=.False., APPEND=.True., &
                    FILE=TRIMLR(filOUT)//h5ext)
  END DO
  DO i=1,Nband
    CALL write_hdf5(DBLARR3D=Fnu_band_tab(:,:,i,:)+(Fnu_cont+Fnu_star)*Pabs, &
                    NAME='FnuBAND'//TRIMLR(PRING(i))//' (MJyovsr)', &
                    COMPRESS=.False., APPEND=.True., &
                    FILE=TRIMLR(filOUT)//h5ext)

  END DO
  
  print*, 'Chi2 fit M83 IRS (15 pix * 15 pix) [done]'//NEW_LINE('')

  print*, '=================================================='

  !! Free memory space
  DEALLOCATE(wOBS, FnuOBS, dFnuOBS, FnuMOD, maskint, &
             mask, compdust, nQabs, nuOBS, &
             parr, limits, limited, status, chi2red)!, parerr, covpar)
  DEALLOCATE(Pabs,Fnu_cont,Fnu_line,Fnu_band,Fnu_star)

END PROGRAM test_fitChi2
