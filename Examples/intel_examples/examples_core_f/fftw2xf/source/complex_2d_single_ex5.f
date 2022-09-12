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
! Complex-to-complex 2D transform for COMPLEX*8 data inplace with threads.
!
! Configuration parameters for Intel MKL DFTI:
!           DFTI_FORWARD_DOMAIN = DFTI_COMPLEX (obligatory)
!           DFTI_PRECISION      = DFTI_SINGLE  (obligatory)
!           DFTI_DIMENSION      = 2            (obligatory)
!           DFTI_LENGTHS        = {n,m}        (obligatory)
!           DFTI_PLACEMENT      = DFTI_INPLACE (default)
!           DFTI_FORWARD_SCALE  = 1.0          (default)
!           DFTI_BACKWARD_SCALE = 1.0/(n*m)    (default=1.0)
!
! Other default configuration parameters are in the mkl_dfti.f90 interface file
!*****************************************************************************

      PROGRAM COMPLEX_2D_SINGLE_EX5

      INCLUDE 'fftw_f77.i'
      INCLUDE 'mkl_fftw_examples.fi'

      INTEGER N,M,NTHREADS,STATUS,I,J
      PARAMETER (M=10)
      PARAMETER (N=15)
      PARAMETER (NTHREADS=2)
      INTEGER*8 MY_PLAN

      COMPLEX*8 IN(N,M),EXP_X(N,M)
      REAL*4 ERR,SCALE

!
!     Initialize IN and copy to expected EXP_X
!
      PRINT *, 'Initialize data array'
      CALL INIT_COMPLEX_VECTOR_C(IN,N*M)
      DO J=1,M
        DO I=1,N
         EXP_X(I,J)=IN(I,J)
        END DO
      END DO

!
!     FFTW requires user to call this routine before calling any other FFTW routines
!
      CALL FFTW_F77_THREADS_INIT(STATUS)

!
!     Create FFTW plan for 2D forward transform
!
      PRINT *, 'Create FFTW plan for 2D forward transform'
      CALL FFTW2D_F77_CREATE_PLAN(MY_PLAN,N,M,FFTW_FORWARD,
     *                FFTW_ESTIMATE+FFTW_IN_PLACE)

!
!     Compute Forward
!
      PRINT *, 'Compute Forward'
      CALL FFTWND_F77_THREADS_ONE(NTHREADS,MY_PLAN,IN,0)

!
!     Destroy FFTW plan
!
      PRINT *, 'Destroy FFTW plan'
      CALL FFTWND_F77_DESTROY_PLAN(MY_PLAN)

!
!     Create FFTW plan for 2D backward transform
!
      PRINT *, 'Create FFTW plan for 2D backward transform'
      CALL FFTW2D_F77_CREATE_PLAN(MY_PLAN,N,M,FFTW_BACKWARD,
     *                FFTW_ESTIMATE+FFTW_IN_PLACE)

!
!     Compute Backward
!
      PRINT *, 'Compute Backward'
      CALL FFTWND_F77_THREADS_ONE(NTHREADS,MY_PLAN,IN,0)

!
!     Destroy FFTW plan
!
      PRINT *, 'Destroy FFTW plan'
      CALL FFTWND_F77_DESTROY_PLAN(MY_PLAN)

!
!     Scale result. FFTW can't do this itself.
!
      PRINT *, 'Scale result by 1/(N*M)'
      SCALE=1.0/REAL((N*M),4)
      DO J=1,M
        DO I=1,N
         IN(I,J)=SCALE*IN(I,J)
        END DO
      END DO

!
!     Check results
!
      PRINT *, 'Check results'
      CALL CHECK_RESULT_C(IN,EXP_X,N*M,ERR)
      PRINT *, 'Accuracy=',ERR
      IF (ERR .GT. MAX_SINGLE_ERR) THEN
       PRINT *, 'TEST FAILED'
       STOP 1
      END IF
      PRINT *, 'TEST PASSED'

      PRINT *, 'END OF TEST'

      END PROGRAM
