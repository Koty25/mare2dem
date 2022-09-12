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
! Real-to-complex and complex-to-real 2D transform for REAL*4 and COMPLEX*8
! data not inplace with strides and threads.
!
! Configuration parameters for Intel MKL DFTI:
!           DFTI_FORWARD_DOMAIN = DFTI_REAL        (obligatory)
!           DFTI_PRECISION      = DFTI_SINGLE      (obligatory)
!           DFTI_DIMENSION      = 2                (obligatory)
!           DFTI_LENGTHS        = {n,m}            (obligatory)
!           DFTI_PLACEMENT      = DFTI_NOT_INPLACE (default=DFTI_INPLACE)
!           DFTI_FORWARD_SCALE  = 1.0              (default)
!           DFTI_BACKWARD_SCALE = 1.0/(n*m)        (default=1.0)
!
! Other default configuration parameters are in the mkl_dfti.f90 interface file
!*****************************************************************************

      PROGRAM REAL_2D_SINGLE_EX8

      INCLUDE 'fftw_f77.i'
      INCLUDE 'mkl_fftw_examples.fi'

      INTEGER SIZE,N,M,OUT_N,NTHREADS,STATUS,IDIST,ODIST,ISTRIDE,OSTRIDE
      INTEGER HOWMANY,I
      PARAMETER (M=5)
      PARAMETER (N=8)
      PARAMETER (OUT_N=INT(N/2.0)+1)
      PARAMETER (ISTRIDE=2)
      PARAMETER (OSTRIDE=3)
      PARAMETER (IDIST=N*M*ISTRIDE+3)
      PARAMETER (ODIST=OUT_N*M*OSTRIDE+5)
      PARAMETER (HOWMANY=2)
      PARAMETER (SIZE=HOWMANY*(IDIST+ODIST))
      PARAMETER (NTHREADS=2)
      INTEGER*8 MY_PLAN

      REAL*4 IN(SIZE),EXP_X(SIZE)
      COMPLEX*8 OUT(SIZE)
      REAL*4 ERR,SCALE

!
!     Initialize IN and copy to expected EXP_X
!
      PRINT *, 'Initialize data array'
      DO I=1,SIZE
        IN(I)=SIN(FLOAT(I))
        EXP_X(I)=IN(I)
      END DO

!
!     FFTW requires user to call this routine before calling any other FFTW routines
!
      CALL FFTW_F77_THREADS_INIT(STATUS)

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
      CALL RFFTWND_F77_THREADS_REAL_TO_COMPLEX(NTHREADS,MY_PLAN,
     *                 HOWMANY,IN,ISTRIDE,IDIST,OUT,OSTRIDE,ODIST)

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
      CALL RFFTWND_F77_THREADS_COMPLEX_TO_REAL(NTHREADS,MY_PLAN,
     *                 HOWMANY,OUT,OSTRIDE,ODIST,IN,ISTRIDE,IDIST)

!
!     Destroy FFTW plan
!
      PRINT *, 'Destroy FFTW plan'
      CALL RFFTWND_F77_DESTROY_PLAN(MY_PLAN)

!
!     Scale result. FFTW can't do this itself.
!
      PRINT *, 'Scale result by 1/(N*M)'
      SCALE=1.0/REAL((N*M),4)
      CALL SCALE_WITH_STRIDES_S(IN,SCALE,HOWMANY,N*M,ISTRIDE,IDIST)

!
!     Check results
!
      PRINT *, 'Check results'
      CALL CHECK_RESULT_S(IN,EXP_X,SIZE,ERR)
      PRINT *, 'Accuracy=',ERR
      IF (ERR .GT. MAX_SINGLE_ERR) THEN
       PRINT *, 'TEST FAILED'
       STOP 1
      END IF
      PRINT *, 'TEST PASSED'

      PRINT *, 'END OF TEST'

      END PROGRAM
