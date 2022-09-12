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
!       through FFTW interface (via wrappers) example program (Fortran-interface)
!
! Real-to-complex and complex-to-real 2D transform for REAL*8 and COMPLEX*16
! data not inplace.
!
! Configuration parameters for Intel MKL DFTI:
!           DFTI_FORWARD_DOMAIN = DFTI_REAL        (obligatory)
!           DFTI_PRECISION      = DFTI_DOUBLE      (obligatory)
!           DFTI_DIMENSION      = 2                (obligatory)
!           DFTI_LENGTHS        = {n,m}            (obligatory)
!           DFTI_PLACEMENT      = DFTI_NOT_INPLACE (default=DFTI_INPLACE)
!           DFTI_FORWARD_SCALE  = 1.0              (default)
!           DFTI_BACKWARD_SCALE = 1.0/(n*m)        (default=1.0)
!
! Other default configuration parameters are in the mkl_dfti.f90 interface file
!*****************************************************************************

      PROGRAM REAL_2D_DOUBLE_EX2

      INCLUDE 'fftw_f77.i'
      INCLUDE 'mkl_fftw_examples.fi'

      INTEGER N,M,OUT_N,I,J
      PARAMETER (M=10)
      PARAMETER (N=15)
      PARAMETER (OUT_N=INT(N/2.0)+1)
      INTEGER*8 MY_PLAN

      REAL*8 IN(N,M),EXP_X(N,M)
      COMPLEX*16 OUT(OUT_N,M)
      REAL*8 ERR,SCALE

!
!     Initialize IN and copy to expected EXP_X
!
      PRINT *, 'Initialize data array'
      DO J=1,M
        DO I=1,N
         IN(I,J)=DSIN(DFLOAT(I*J))
         EXP_X(I,J)=IN(I,J)
        END DO
      END DO

!
!     Create FFTW plan for 2D real to complex transform
!
      PRINT *, 'Create FFTW plan for 2D real to complex transform'
      CALL RFFTW2D_F77_CREATE_PLAN(MY_PLAN,N,M,FFTW_REAL_TO_COMPLEX,
     *                 FFTW_ESTIMATE)

!
!     Compute 2D real to complex transform
!
      PRINT *, 'Compute 2D real to complex transform'
      CALL RFFTWND_F77_ONE_REAL_TO_COMPLEX(MY_PLAN,IN,OUT)

!
!     Destroy FFTW plan
!
      PRINT *, 'Destroy FFTW plan'
      CALL RFFTWND_F77_DESTROY_PLAN(MY_PLAN)

!
!     Create FFTW plan for 2D complex to real transform
!
      PRINT *, 'Create FFTW plan for 2D complex to real transform'
      CALL RFFTW2D_F77_CREATE_PLAN(MY_PLAN,N,M,FFTW_COMPLEX_TO_REAL,
     *                 FFTW_ESTIMATE)

!
!     Compute 2D complex to real transform
!
      PRINT *, 'Compute 2D complex to real transform'
      CALL RFFTWND_F77_ONE_COMPLEX_TO_REAL(MY_PLAN,OUT,IN)

!
!     Destroy FFTW plan
!
      PRINT *, 'Destroy FFTW plan'
      CALL RFFTWND_F77_DESTROY_PLAN(MY_PLAN)

!
!     Scale result. FFTW can't do this itself.
!
      PRINT *, 'Scale result by 1/(N*M)'
      SCALE=1.0D0/(N*M)
      DO J=1,M
        DO I=1,N
         IN(I,J)=SCALE*IN(I,J)
        END DO
      END DO

!
!     Check results
!
      PRINT *, 'Check results'
      CALL CHECK_RESULT_D(IN,EXP_X,N*M,ERR)
      PRINT *, 'Accuracy=',ERR
      IF (ERR .GT. MAX_DOUBLE_ERR) THEN
       PRINT *, 'TEST FAILED'
       STOP 1
      END IF
      PRINT *, 'TEST PASSED'

      PRINT *, 'END OF TEST'

      END PROGRAM
