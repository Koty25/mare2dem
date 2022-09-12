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
! Real-to-complex and complex-to-real 3D transform for REAL*8 and COMPLEX*16
! data not inplace with strides and threads.
!
! Configuration parameters for Intel MKL DFTI:
!           DFTI_FORWARD_DOMAIN = DFTI_REAL        (obligatory)
!           DFTI_PRECISION      = DFTI_DOUBLE      (obligatory)
!           DFTI_DIMENSION      = 3                (obligatory)
!           DFTI_LENGTHS        = {n,m,l}          (obligatory)
!           DFTI_PLACEMENT      = DFTI_NOT_INPLACE (default=DFTI_INPLACE)
!           DFTI_FORWARD_SCALE  = 1.0              (default)
!           DFTI_BACKWARD_SCALE = 1.0/(n*m*l)      (default=1.0)
!
! Other default configuration parameters are in the mkl_dfti.f90 interface file
!*****************************************************************************

      PROGRAM REAL_3D_DOUBLE_EX8

      INCLUDE 'fftw_f77.i'
      INCLUDE 'mkl_fftw_examples.fi'

      INTEGER SIZE,N,M,L,OUT_N,NTHREADS,STATUS,IDIST,ODIST,ISTRIDE
      INTEGER OSTRIDE,HOWMANY,I
      PARAMETER (L=3)
      PARAMETER (M=5)
      PARAMETER (N=8)
      PARAMETER (OUT_N=INT(N/2.0)+1)
      PARAMETER (ISTRIDE=2)
      PARAMETER (OSTRIDE=3)
      PARAMETER (IDIST=N*M*L*ISTRIDE+3)
      PARAMETER (ODIST=OUT_N*M*L*OSTRIDE+5)
      PARAMETER (HOWMANY=2)
      PARAMETER (SIZE=HOWMANY*(IDIST+ODIST))
      PARAMETER (NTHREADS=2)
      INTEGER*8 MY_PLAN

      REAL*8 IN(SIZE),EXP_X(SIZE)
      COMPLEX*16 OUT(SIZE)
      REAL*8 ERR,SCALE

!
!     Initialize IN and copy to expected EXP_X
!
      PRINT *, 'Initialize data array'
      DO I=1,SIZE
        IN(I)=DSIN(DFLOAT(I))
        EXP_X(I)=IN(I)
      END DO

!
!     FFTW requires user to call this routine before calling any other FFTW routines
!
      CALL FFTW_F77_THREADS_INIT(STATUS)

!
!     Create FFTW plan for 3D real to complex transform
!
      PRINT *, 'Create FFTW plan for 3D real to complex transform'
      CALL RFFTW3D_F77_CREATE_PLAN(MY_PLAN,N,M,L,FFTW_REAL_TO_COMPLEX,
     *                 FFTW_ESTIMATE)

!
!     Compute 3D real to complex transform
!
      PRINT *, 'Compute 3D real to complex transform'
      CALL RFFTWND_F77_THREADS_REAL_TO_COMPLEX(NTHREADS,MY_PLAN,
     *                 HOWMANY,IN,ISTRIDE,IDIST,OUT,OSTRIDE,ODIST)

!
!     Destroy FFTW plan
!
      PRINT *, 'Destroy FFTW plan'
      CALL RFFTWND_F77_DESTROY_PLAN(MY_PLAN)

!
!     Create FFTW plan for 3D complex to real transform
!
      PRINT *, 'Create FFTW plan for 3D complex to real transform'
      CALL RFFTW3D_F77_CREATE_PLAN(MY_PLAN,N,M,L,FFTW_COMPLEX_TO_REAL,
     *                 FFTW_ESTIMATE)

!
!     Compute 3D complex to real transform
!
      PRINT *, 'Compute 3D complex to real transform'
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
      PRINT *, 'Scale result by 1/(N*M*L)'
      SCALE=1.0D0/(N*M*L)
      CALL SCALE_WITH_STRIDES_D(IN,SCALE,HOWMANY,N*M*L,ISTRIDE,IDIST)

!
!     Check results
!
      PRINT *, 'Check results'
      CALL CHECK_RESULT_D(IN,EXP_X,SIZE,ERR)
      PRINT *, 'Accuracy=',ERR
      IF (ERR .GT. MAX_DOUBLE_ERR) THEN
       PRINT *, 'TEST FAILED'
       STOP 1
      END IF
      PRINT *, 'TEST PASSED'

      PRINT *, 'END OF TEST'

      END PROGRAM
