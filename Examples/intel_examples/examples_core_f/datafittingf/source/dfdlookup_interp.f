!===============================================================================
! Copyright 2010-2020 Intel Corporation.
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

!  Content:
!    Interpolation with Look up method
!    Example Program Text
!*******************************************************************************

      include 'mkl_df.f90'
      include "errcheck.inc"
      include "rescheck.inc"
      include "generatedata.inc"

      PROGRAM MKL_DF_TEST

      USE MKL_DF_TYPE
      USE MKL_DF
      USE DF_GENERATE_DATA
      USE DF_EXAMPLE_RESCHECK

      ! number of break points
      INTEGER,PARAMETER :: N          = 7
      ! number of functions
      INTEGER,PARAMETER :: NNY        = 1
      ! total number of interpolation sites
      INTEGER,PARAMETER :: NNSITE     = N
      ! size of array describing derivative orders to compute
      INTEGER,PARAMETER :: NNDERORDER = 1
      ! number of derivatives to compute
      INTEGER,PARAMETER :: NNDER      = 1

      ! left limit of interpolation interval
      INTEGER,PARAMETER :: LEFT_LIMIT  = 1.0d0
      ! right limit of interpolation interval
      INTEGER,PARAMETER :: RIGHT_LIMIT = 2.0d0

      INTEGER,PARAMETER :: FFREQ = 1.7d0

      ! Data Fitting task descriptor
      TYPE(DF_TASK) task

      ! spline order
      INTEGER :: sorder
      ! spline type
      INTEGER :: stype
      ! number of break points
      INTEGER :: nx
      ! additional info about break points
      INTEGER :: xhint
      ! number of functions
      INTEGER :: ny
      ! total number of interpolation sites
      INTEGER :: nsite
      ! additional info about interpolation sites
      INTEGER :: sitehint
      ! size of array describing derivative orders
      INTEGER :: ndorder
      ! array describing derivative orders
      INTEGER :: derorder(1)
      ! interpolation results storage format
      INTEGER :: rhint

      ! array of break points
      REAL(8) :: x(N)
      ! function values
      REAL(8) :: y(NNY*N)
      ! array of interpolation sites
      REAL(8) :: site(NNSITE)
      ! spline evaluation results
      REAL(8) :: r(NNSITE)
      ! limits of the interpolation interval
      REAL(8) :: left,right
      ! type of calculations
      INTEGER :: type
      ! method that is used to perform calculations
      INTEGER :: method

      REAL(8) :: freq
      INTEGER :: i,j
      INTEGER :: errnums
      INTEGER(4) :: errcode

      errcode = 0
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****

      sorder = DF_PP_STD
      stype = DF_LOOKUP_INTERPOLANT

      !***** Parameters describing interpolation interval *****
      left = LEFT_LIMIT
      right = RIGHT_LIMIT
      nx = N
      xhint = DF_NON_UNIFORM_PARTITION

      !***** Parameters describing function *****
      ny = NNY

      !***** Parameters describing interpolation sites *****
      nsite      = NNSITE;
      sitehint   = DF_NON_UNIFORM_PARTITION

      !***** Parameter describing interpolation results storage *****
      rhint = DF_MATRIX_STORAGE_ROWS

      !***** Parameter describing array for derivative orders *****
      ndorder = NNDERORDER
      derorder = (/1/)

      !***** Number of derivatives to compute *****
      nder = NNDER

      !***** Generate array of uniformly distributed break points *****
      errcode = dUniformRandSortedData( x, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Generate function y = sin(2 * Pi * freq * x) *****
      freq = FFREQ
      errcode = dSinDataNotUniformGrid( y, x, freq , nx )
      CALL CheckDfError(errcode)

      !***** Generate interpolation sites ******
      DO j = 1, nx
        site(j) = x(j)
      END DO

      !***** Create Data Fitting task *****
      errcode = dfdNewTask1D( task, nx, x, xhint, ny, y )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for look up interpolant *****
      errcode = dfdEditPPSpline1D( task, sorder, stype )
      CALL CheckDfError(errcode)

      !***** Interpolate using look up method *****
      type = DF_INTERP
      method = DF_METHOD_PP

      errcode = dfdInterpolate1D( task, type, method, nsite,             &
     &    site, sitehint, ndorder, derorder, r=r, rhint=rhint )
      CALL CheckDfError(errcode)

      !***** Delete Data Fitting task *****
      errcode = dfDeleteTask( task )
      CALL CheckDfError(errcode)

      !***** Check results of interpolation *****
      DO j = 1, nsite
        IF ( ABS( r(j) - y(j) ) > EPSILON_DOUBLE ) errnums = errnums+1
      END DO

      !***** Print results *****
      WRITE (*,901) "Number of break points : ",nx
      WRITE (*,901) "Number of sites : ",nsite

      !***** Print given function and computed results *****
      WRITE (*,902) "    X             Y(X)          R(X) "
      DO j = 1,nx
        WRITE (*,903) " ",x(j),"   ",y(j),"   ",r(j)
      END DO

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,907) "","Error: Computed interpolation results "
        WRITE (*,904) " are incorrect"
        STOP 1
      ELSE
        WRITE (*,907) "","Computed interpolation results "
        WRITE (*,904) " are correct"
      END IF
      STOP 0

  901 FORMAT (A,I0)
  902 FORMAT (/A)
  903 FORMAT (A,SP,F11.6,A,SP,F11.6,A,SP,F11.6)
  904 FORMAT (99A)
  905 FORMAT (A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,    &
     & A,SP,F11.6)
  906 FORMAT (A,SP,F11.6,A,SP,F11.6)
  907 FORMAT (//A)
      END PROGRAM
