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
! Complex-to-complex 1D transform for COMPLEX*16 data inplace with strides
! and threads
!
! Configuration parameters for Intel MKL DFTI:
!           DFTI_FORWARD_DOMAIN = DFTI_COMPLEX (obligatory)
!           DFTI_PRECISION      = DFTI_DOUBLEE (obligatory)
!           DFTI_DIMENSION      = 1            (obligatory)
!           DFTI_LENGTHS        = n            (obligatory)
!           DFTI_PLACEMENT      = DFTI_INPLACE (default)
!           DFTI_FORWARD_SCALE  = 1.0          (default)
!           DFTI_BACKWARD_SCALE = 1.0/n        (default=1.0)
!
! Other default configuration parameters are in the mkl_dfti.f90 interface file
!*****************************************************************************

      PROGRAM COMPLEX_1D_DOUBLE_EX7

      INCLUDE 'fftw_f77.i'
      INCLUDE 'mkl_fftw_examples.fi'

      INTEGER SIZE,N,IDIST,ISTRIDE,HOWMANY,STATUS,NTHREADS,I
      PARAMETER (N=10)
      PARAMETER (ISTRIDE=2)
      PARAMETER (IDIST=N*ISTRIDE+3)
      PARAMETER (HOWMANY=2)
      PARAMETER (SIZE=HOWMANY*IDIST)
      PARAMETER (NTHREADS=2)
      INTEGER*8 MY_PLAN

      COMPLEX*16 IN(SIZE),OUT(N),EXP_X(SIZE)
      REAL*8 ERR,SCALE

!
!     Initialize IN and copy to expected EXP_X
!
      PRINT *, 'Initialize data array'
      CALL INIT_COMPLEX_VECTOR_Z(IN,SIZE)
      DO I=1,SIZE
       EXP_X(I)=IN(I)
      END DO

!
!     FFTW requires user to call this routine before calling any other FFTW routines
!
      CALL FFTW_F77_THREADS_INIT(STATUS)

!
!     Create FFTW plan for 1D forward transform
!
      PRINT *, 'Create FFTW plan for 1D forward transform'
      CALL FFTW_F77_CREATE_PLAN(MY_PLAN,N,FFTW_FORWARD,
     *              FFTW_ESTIMATE+FFTW_IN_PLACE)

!
!     Compute Forward
!
      PRINT *, 'Compute Forward'
      CALL FFTW_F77_THREADS(NTHREADS,MY_PLAN,HOWMANY,IN,ISTRIDE,IDIST,
     *              OUT,0,0)

!
!     Destroy FFTW plan
!
      PRINT *, 'Destroy FFTW plan'
      CALL FFTW_F77_DESTROY_PLAN(MY_PLAN)

!
!     Create FFTW plan for 1D backward transform
!
      PRINT *, 'Create FFTW plan for 1D backward transform'
      CALL FFTW_F77_CREATE_PLAN(MY_PLAN,N,FFTW_BACKWARD,
     *              FFTW_ESTIMATE+FFTW_IN_PLACE)

!
!     Compute Backward
!
      PRINT *, 'Compute Backward'
      CALL FFTW_F77_THREADS(NTHREADS,MY_PLAN,HOWMANY,IN,ISTRIDE,IDIST,
     *              OUT,0,0)

!
!     Destroy FFTW plan
!
      PRINT *, 'Destroy FFTW plan'
      CALL FFTW_F77_DESTROY_PLAN(MY_PLAN)

!
!     Scale result. FFTW can't do this itself.
!
      PRINT *, 'Scale result by 1/N'
      SCALE=1.0D0/N
      CALL SCALE_WITH_STRIDES_Z(IN,SCALE,HOWMANY,N,ISTRIDE,IDIST)

!
!     Check results
!
      PRINT *, 'Check results'
      CALL CHECK_RESULT_Z(IN,EXP_X,SIZE,ERR)
      PRINT *, 'Accuracy=',ERR
      IF (ERR .GT. MAX_DOUBLE_ERR) THEN
       PRINT *, 'TEST FAILED'
       STOP 1
      END IF
      PRINT *, 'TEST PASSED'

      PRINT *, 'END OF TEST'

      END PROGRAM
