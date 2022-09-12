!===============================================================================
! Copyright 2006-2020 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

! Content:
!       Intel(R) Math Kernel Library (Intel(R) MKL) DFTI implementation
!       through FFTW interface (via wrappers) example support functions
!
!*****************************************************************************

      SUBROUTINE INIT_COMPLEX_VECTOR_Z(X,N)

!
!     Input parameters
!
      COMPLEX*16 X(*)
      INTEGER N

!
!     Local parameters
!
      INTEGER I
      REAL*8 F_STEP

!
!     Body
!
      DO I=1,N
          F_STEP=DFLOAT(I)
          X(I)=DCMPLX((DSIN(F_STEP)*DSQRT(3.0D0))/2.0D0,
     *                DSIN(F_STEP)/DSQRT(3.0D0))
      END DO

      END SUBROUTINE

      SUBROUTINE INIT_COMPLEX_VECTOR_C(X,N)

!
!     Input parameters
!
      COMPLEX*8 X(*)
      INTEGER N

!
!     Local parameters
!
      INTEGER I
      REAL*4 F_STEP

!
!     Body
!
      DO I=1,N
          F_STEP=FLOAT(I)
          X(I)=CMPLX((SIN(F_STEP)*SQRT(3.0))/2.0,
     *               SIN(F_STEP)/SQRT(3.0))
      END DO

      END SUBROUTINE

      SUBROUTINE CHECK_RESULT_Z(IN,EXP_X,N,ERR)

!
!     Input parameters
!
      COMPLEX*16 IN(*),EXP_X(*)
      INTEGER N

!
!     Output parameters
!
      REAL*8 ERR

!
!     Local parameters
!
      COMPLEX*16 D
      REAL*8 E
      INTEGER I

!
!     Body
!
      E = 0.0D0
      DO I=1,N
        D=EXP_X(I)-IN(I)
        IF (CDABS(D) .GT. E) THEN
            E=CDABS(D)
        END IF
      END DO
      ERR=E

      END SUBROUTINE

      SUBROUTINE CHECK_RESULT_C(IN,EXP_X,N,ERR)

!
!     Input parameters
!
      COMPLEX*8 IN(*),EXP_X(*)
      INTEGER N

!
!     Output parameters
!
      REAL*4 ERR

!
!     Local parameters
!
      COMPLEX*8 D
      REAL*4 E
      INTEGER I

!
!     Body
!
      E = 0.0
      DO I=1,N
        D=EXP_X(I)-IN(I)
        IF (CABS(D) .GT. E) THEN
            E=CABS(D)
        END IF
      END DO
      ERR=E

      END SUBROUTINE

      SUBROUTINE CHECK_RESULT_D(IN,EXP_X,N,ERR)
!
!     Input parameters
!
      REAL*8 IN(*),EXP_X(*)
      INTEGER N

!
!     Output parameters
!
      REAL*8 ERR

!
!     Local parameters
!
      REAL*8 E,D
      INTEGER I

!
!     Body
!
      E = 0.0D0
      DO I=1,N
        D=EXP_X(I)-IN(I)
        IF (DABS(D) .GT. E) THEN
            E=DABS(D)
        END IF
      END DO
      ERR=E

      END SUBROUTINE


      SUBROUTINE CHECK_RESULT_MULTIPLE_D(IN,EXP_X,N,HOWMANY,
     *                                   ISTRIDE,IDIST,ERR)
!
!     Input parameters
!
      REAL*8 IN(*),EXP_X(*)
      INTEGER N,HOWMANY,ISTRIDE,IDIST

!
!     Output parameters
!
      REAL*8 ERR

!
!     Local parameters
!
      REAL*8 E,D
      INTEGER I, IS, ID, J

!
!     Body
!
      E = 0.0D0
      ID=1
      DO J=1,HOWMANY
        DO I=1,N
          IS = ID + (I-1)*ISTRIDE
          D=EXP_X(IS)-IN(IS)
          IF (DABS(D) .GT. E) THEN
            E=DABS(D)
          END IF
        END DO
        ID=ID+IDIST
      END DO
      ERR=E
      END SUBROUTINE

      SUBROUTINE CHECK_RESULT_MULTIPLE_S(IN,EXP_X,N,HOWMANY,
     *                                   ISTRIDE,IDIST,ERR)
!
!     Input parameters
!
      REAL*4 IN(*),EXP_X(*)
      INTEGER N,HOWMANY,ISTRIDE,IDIST

!
!     Output parameters
!
      REAL*4 ERR

!
!     Local parameters
!
      REAL*4 E,D
      INTEGER I, IS, ID, J

!
!     Body
!
      E = 0.0
      ID=1
      DO J=1,HOWMANY
        DO I=1,N
          IS = ID + (I-1)*ISTRIDE
          D=EXP_X(IS)-IN(IS)
          IF (ABS(D) .GT. E) THEN
            E=ABS(D)
          END IF
        END DO
        ID=ID+IDIST
      END DO
      ERR=E
      END SUBROUTINE

      SUBROUTINE CHECK_RESULT_S(IN,EXP_X,N,ERR)

!
!     Input parameters
!
      REAL*4 IN(*),EXP_X(*)
      INTEGER N

!
!     Output parameters
!
      REAL*4 ERR

!
!     Local parameters
!
      REAL*4 E,D
      INTEGER I

!
!     Body
!
      E = 0.0
      DO I=1,N
        D=EXP_X(I)-IN(I)
        IF (ABS(D) .GT. E) THEN
            E=ABS(D)
        END IF
      END DO
      ERR=E

      END SUBROUTINE

      SUBROUTINE SCALE_WITH_STRIDES_Z(IN,SCALE,HOWMANY,N,ISTRIDE,IDIST)

!
!     Input parameters
!
      COMPLEX*16 IN(*)
      INTEGER HOWMANY,N,ISTRIDE,IDIST
      REAL*8 SCALE

!
!     Local parameters
!
      INTEGER I,J,K

!
!     Body
!
      DO I=1,HOWMANY
        DO J=1,N
          K=(I-1)*IDIST+(J-1)*ISTRIDE+1
          IN(K)=SCALE*IN(K)
        END DO
      END DO

      END SUBROUTINE

      SUBROUTINE SCALE_WITH_STRIDES_C(IN,SCALE,HOWMANY,N,ISTRIDE,IDIST)

!
!     Input parameters
!
      COMPLEX*8 IN(*)
      INTEGER HOWMANY,N,ISTRIDE,IDIST
      REAL*4 SCALE

!
!     Local parameters
!
      INTEGER I,J,K

!
!     Body
!
      DO I=1,HOWMANY
        DO J=1,N
          K=(I-1)*IDIST+(J-1)*ISTRIDE+1
          IN(K)=SCALE*IN(K)
        END DO
      END DO

      END SUBROUTINE

      SUBROUTINE SCALE_WITH_STRIDES_D(IN,SCALE,HOWMANY,N,ISTRIDE,IDIST)

!
!     Input parameters
!
      REAL*8 IN(*)
      INTEGER HOWMANY,N,ISTRIDE,IDIST
      REAL*8 SCALE

!
!     Local parameters
!
      INTEGER I,J,K

!
!     Body
!
      DO I=1,HOWMANY
        DO J=1,N
          K=(I-1)*IDIST+(J-1)*ISTRIDE+1
          IN(K)=SCALE*IN(K)
        END DO
      END DO

      END SUBROUTINE

      SUBROUTINE SCALE_WITH_STRIDES_S(IN,SCALE,HOWMANY,N,ISTRIDE,IDIST)

!
!     Input parameters
!
      REAL*4 IN(*)
      INTEGER HOWMANY,N,ISTRIDE,IDIST
      REAL*4 SCALE

!
!     Local parameters
!
      INTEGER I,J,K

!
!     Body
!
      DO I=1,HOWMANY
        DO J=1,N
          K=(I-1)*IDIST+(J-1)*ISTRIDE+1
          IN(K)=SCALE*IN(K)
        END DO
      END DO

      END SUBROUTINE
