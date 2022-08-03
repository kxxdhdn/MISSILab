!******************************************************************************
!*
!*                      VARIOUS APPLICATIONS OF FFT
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 07/2013 
  ! 
  ! 3) DESCRIPTION: implements FFTs, and various relaed functions.
  !==========================================================================

MODULE FFT_specials

  USE utilities, ONLY:
  USE constants, ONLY:
  USE arrays, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: realft, correl


CONTAINS


  !==========================================================================
  ! CALL FOURROW(data,isign)
  !
  !   Replaces each row (constant first index) of data(1:M,1:N) by its discrete
  ! Fourier transform (transform on second index), if isign is input as 1; or 
  ! replaces each row of data by N times its inverse discrete Fourier 
  ! transform, if isign is input as −1. N must be an integer power of 2. 
  ! Parallelism is M-fold on the first index of data.
  !==========================================================================

  SUBROUTINE fourrow (data,isign)

    USE utilities, ONLY: DP, CDP, swap, strike
    USE constants, ONLY: pi
    IMPLICIT NONE
    COMPLEX(CDP), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: isign
    INTEGER :: N, i, istep, j, m, Mmax, N2
    REAL(DP) :: theta
    COMPLEX(CDP), DIMENSION(SIZE(data,1)) :: temp
    COMPLEX(CDP) :: w,wp
    COMPLEX(CDP) :: ws

    !-------------------------------------------------------------------------

    N = SIZE(data,2)
    IF (IAND(N,N-1) /= 0) CALL STRIKE("FOOURROW","n must be a power of 2")
    N2 = N / 2
    j = N2
    DO i=1,N-2
      IF (j > i) CALL SWAP(data(:,j+1),data(:,i+1))
      m = N2
      DO
        IF (m < 2 .OR. j < m) EXIT
        j = j - m
        m = m / 2
      END DO
      j = j + m
    END DO
    Mmax = 1
    DO
      IF (N <= Mmax) EXIT
      istep = 2 * Mmax
      theta = pi / (isign*Mmax)
      wp = CMPLX(-2._DP*SIN(0.5_DP*theta)**2,SIN(theta),KIND=CDP)
      w = CMPLX(1._DP,0._DP,KIND=CDP)
      DO m=1,Mmax
        ws = w
        DO i=m,n,istep
          j = i + Mmax
          temp = ws * data(:,j)
          data(:,j) = data(:,i) - temp
          data(:,i) = data(:,i) + temp
        END DO
        w = w * wp + w
      END DO
      Mmax = istep
    END DO

    !-------------------------------------------------------------------------

  END SUBROUTINE fourrow


  !==========================================================================
  ! CALL FOUR1(data,isign)
  !
  !   Replaces a complex array data by its discrete Fourier transform, if 
  ! isign is input as 1; or replaces data by its inverse discrete Fourier 
  ! transform times the size of data, if isign is input as −1. The size of 
  ! data must be an integer power of 2. Parallelism is achieved by internally 
  ! reshaping the input array to two dimensions. (Use this version if fourrow 
  ! is faster than fourcol on your machine.)
  !==========================================================================

  SUBROUTINE four1 (data,isign)

    USE utilities, ONLY: DP, CDP, strike, arth
    USE constants, ONLY: twopi
    IMPLICIT NONE
    COMPLEX(CDP), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: isign
    COMPLEX(CDP), DIMENSION(:,:), ALLOCATABLE :: dat, temp
    COMPLEX(CDP), DIMENSION(:), ALLOCATABLE :: w, wp
    REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
    INTEGER :: N, m1, m2, j

    !-------------------------------------------------------------------------

    N = SIZE(data)
    IF (IAND(N,N-1) /=0) CALL STRIKE("FOUR1","N must be a power of 2")
    m1 = 2**CEILING( 0.5_DP * LOG(REAL(N,DP)) / 0.693147_DP )
    m2 = n / m1
    ALLOCATE (dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat = RESHAPE(data,SHAPE(dat))
    CALL FOURROW(dat,isign)
    theta = ARTH(0,isign,m1) * TWOPI / N
    wp = CMPLX(-2._DP*SIN(0.5_DP*theta)**2,SIN(theta),KIND=CDP)
    w = CMPLX(1._DP,0._DP,KIND=CDP)
    DO j=2,m2
      w = w * wp + w
      dat(:,j) = dat(:,j) * w
    END DO
    temp = TRANSPOSE(dat)
    CALL FOURROW(temp,isign)
    data = RESHAPE(temp,SHAPE(data))
    DEALLOCATE (dat,w,wp,theta,temp)

    !-------------------------------------------------------------------------

  END SUBROUTINE four1


  !==========================================================================
  ! CALL REALFT(data,isign,zdata)
  !
  !   When isign = 1, calculates the Fourier transform of a set of N 
  ! real-valued data points, input in the array data. If the optional argument 
  ! zdata is not present, the data are replaced by the positive frequency half 
  ! of its complex Fourier transform. The real-valued first and last components 
  ! of the complex transform are returned as elements data(1) and data(2),
  ! respectively. If the complex array zdata of length N/2 is present, data is 
  ! unchanged and the transform is returned in zdata. N must be a power of 2. 
  ! If isign = 1, this routine calculates the inverse transform of a complex 
  ! data array if it is the transform of real data. (Result in this case must 
  ! be multiplied by 2/N.) The data can be supplied either in data, with zdata 
  ! absent, or in zdata.
  !==========================================================================

  SUBROUTINE REALFT(data,isign,zdata)

    USE utilities, ONLY: DP, CDP, strike, zroots_unity
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: isign
    COMPLEX(CDP), DIMENSION(:), OPTIONAL, TARGET :: zdata
    INTEGER :: N, Ndum, Nh, Nq
    COMPLEX(CDP), DIMENSION(SIZE(data)/4) :: w
    COMPLEX(CDP), DIMENSION(SIZE(data)/4-1) :: h1, h2
    COMPLEX(CDP), DIMENSION(:), POINTER :: cdata
    COMPLEX(CDP) :: z
    REAL(DP) :: c1 = 0.5_DP, c2

    !-------------------------------------------------------------------------

    N = SIZE(data)
    IF (IAND(N,N-1) /= 0) CALL STRIKE("REALFT","N must be a power of 2")
    Nh = N / 2
    Nq = N / 4
    IF (PRESENT(zdata)) THEN
      Ndum = N / 2 
      IF (SIZE(zdata) /= N/2) CALL STRIKE("REALFT","Wonrg size of input")
      cdata => zdata
      IF (isign == 1) cdata = CMPLX(data(1:N-1:2),data(2:N:2),KIND=CDP)
    ELSE
      ALLOCATE (cdata(N/2))
      cdata = CMPLX(data(1:N-1:2),data(2:N:2),KIND=CDP)
    END IF

    IF (isign == 1) THEN
      c2 = - 0.5_DP
      CALL four1(cdata,+1)
    ELSE
      c2 = 0.5_DP
    END IF
    w = ZROOTS_UNITY(SIGN(N,isign),N/4)
    w = CMPLX(-AIMAG(w),REAL(w),KIND=CDP)
    h1 = c1 * (cdata(2:Nq)+CONJG(cdata(Nh:Nq+2:-1)))
    h2 = c2 * (cdata(2:Nq)-CONJG(cdata(Nh:Nq+2:-1)))
    cdata(2:Nq) = h1 + w(2:Nq) * h2
    cdata(Nh:Nq+2:-1) = CONJG(h1-w(2:Nq)*h2)
    z = cdata(1)
    IF (isign == 1) THEN
      cdata(1) = CMPLX(REAL(z)+AIMAG(z),REAL(z)-AIMAG(z),KIND=CDP)
    ELSE
      cdata(1) = CMPLX(c1*(REAL(z)+AIMAG(z)),c1*(REAL(z)-AIMAG(z)),KIND=CDP)
      CALL FOUR1(cdata,-1)
    END IF
    IF (PRESENT(zdata)) THEN
      IF (isign /= 1) THEN
        data(1:N-1:2) = REAL(cdata)
         data(2:N:2) = AIMAG(cdata)
      END IF
    ELSE
      data(1:N-1:2) = REAL(cdata)
      data(2:N:2) = AIMAG(cdata)
      DEALLOCATE(cdata)
    END IF

    !-------------------------------------------------------------------------

  END SUBROUTINE realft


  !==========================================================================
  ! f[N] = CORREL(data1[N],data2[N])
  !
  !   Computes the correlation of two arrays for a lag: f(N) to f(N/2+1) are
  ! negative lags, f(1) is lag 0, and f(1) to f(N/2) are positive lags.
  !==========================================================================

  FUNCTION correl (data1,data2)
   
    USE utilities, ONLY: DP, CDP, strike
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: data1, data2
    REAL(DP), DIMENSION(SIZE(data1)) :: correl

    COMPLEX(CDP), DIMENSION(SIZE(data1)/2) :: cdat1, cdat2
    INTEGER :: No2, N

    !-------------------------------------------------------------------------

    ! Preliminaries
    N = SIZE(data1)
    IF (SIZE(data2) /= N) CALL STRIKE("CORREL","Wrong size of input")
    IF (IAND(N,N-1) /= 0) CALL STRIKE("CORREL","N must be a power of 2")
    No2 = N / 2

    ! Transform both data vectors
    CALL REALFT(data1,1,cdat1) 
    CALL REALFT(data2,1,cdat2)
    
    ! Multiply to find FFT of their correlation
    cdat1(1) = CMPLX( REAL(cdat1(1)) * REAL(cdat2(1)) / No2, & 
                      AIMAG(cdat1(1)) * AIMAG(cdat2(1)) / No2, KIND=CDP ) 
    cdat1(2:) = cdat1(2:) * CONJG(cdat2(2:)) / No2

    ! Inverse transform gives correlation
    CALL REALFT(correl,-1,cdat1) 

    !-------------------------------------------------------------------------

  END FUNCTION correl


END MODULE FFT_specials
