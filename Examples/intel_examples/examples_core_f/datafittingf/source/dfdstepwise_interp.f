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
!    Continuous left step-wise constant interpolation
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
      ! total number of interpolation
      INTEGER,PARAMETER :: NNSITE     = (N-1)
      ! size of array describing derivative orders to compute
      INTEGER,PARAMETER :: NNDERORDER = 1
      ! number of derivatives to compute
      INTEGER,PARAMETER :: NNDER      = NNDERORDER

      ! left limit of interpolation interval
      REAL(8),PARAMETER :: LEFT_LIMIT  = 1.0d0
      ! right limit of interpolation interval
      REAL(8),PARAMETER :: RIGHT_LIMIT = 3.0d0

      REAL(8),PARAMETER :: FFREQ = 0.5d0

      ! Data Fitting task descriptor
      TYPE (DF_TASK) task

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
      ! array of break points
      REAL(8) :: x(N)
      ! function values
      REAL(8) :: y(NNY*N)

      ! total number of interpolation sites
      INTEGER :: nsite
      ! additional info about interpolation sites
      INTEGER :: sitehint
      ! array of interpolation sites
      REAL(8) :: site(NNSITE)
      ! size of array describing derivative orders
      INTEGER :: ndorder
      ! array describing derivative orders
      INTEGER :: derorder(NNDERORDER)
      ! interpolation results storage format
      INTEGER :: rhint
      ! type of the computations
      INTEGER :: type
      ! spline evaluation results
      REAL(8) :: r(NNSITE)
      ! indices of cells containing interpolation sites
      INTEGER :: cell(NNSITE)
      ! limits of the interpolation interval
      REAL(8) :: left,right

      REAL(8) :: freq

      INTEGER :: test_cell(NNSITE)
      INTEGER :: i,j
      INTEGER :: errnums
      INTEGER(4) :: errcode

      errcode = 0
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****

      sorder = DF_PP_STD
      stype = DF_CR_STEPWISE_CONST_INTERPOLANT

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

      !***** Parameter describing array for derivative orders *****
      ndorder = NNDERORDER
      derorder = (/1/)

      !***** Number of derivatives to compute *****
      nder = NNDER

      !***** Parameter describing interpolation results storage *****
      rhint = DF_MATRIX_STORAGE_ROWS

      !***** Generate array of uniformly distributed break points *****
      errcode = dUniformRandSortedData( x, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Generate function y = sin(2 * Pi * freq * x) *****
      freq = FFREQ
      errcode = dSinDataNotUniformGrid( y, x, freq , nx )
      CALL CheckDfError(errcode)

      !***** Generate interpolation sites ******
      DO j = 1,nsite
        site(j) = (x(j)+x(j+1))/2
      END DO

      !***** Create Data Fitting task *****
      errcode = dfdnewtask1d( task, nx, x, xhint, ny, y )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for look up interpolant *****
      errcode = dfdeditppspline1d( task, sorder, stype )
      CALL CheckDfError(errcode)

      !***** Interpolate using look up method *****
      type = IOR(DF_CELL,DF_INTERP)
      errcode = dfdinterpolate1d( task, type, DF_METHOD_PP, nsite, site, &
     &      sitehint, ndorder, derorder, r=r, rhint=rhint, cell=cell )
        CALL CheckDfError(errcode)

      !***** Check search results ******
      errcode = dFindCells( nx, x, nsite, site, test_cell )
      test_cell = test_cell - 1
      DO i = 1, nsite
        IF ( test_cell(i) /= cell(i) ) errnums = errnums+1
      END DO

      !***** Check results of interpolation *****
      DO j = 1,nsite
        IF ( ABS( r(j) - y(j) ) > EPSILON_DOUBLE ) errnums = errnums+1
      END DO

      !***** Print results *****
      WRITE (*,901) "Number of break points : ",nx
      WRITE (*,901) "Number of sites : ",nsite

      !***** Print given function and computed results *****
      WRITE (*,902) "    X             Y(X)"
      DO j = 1,nx
        WRITE (*,903) " ",x(j),"   ",y(j)
      END DO
      WRITE (*,902) "Results of cell search and interpolate:"
      WRITE (*,904) "    Site          Computed idx    Expected idx      &
     &  Computed res      Expected res"
      DO i = 1, nsite
       WRITE (*,905) " ",site(i),"     ",cell(i),"     ",                &
     &  test_cell(i),"        ",r(i),"        ",y(i)
      END DO

      !***** Delete Data Fitting task *****
      errcode = dfdeletetask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,907) "","Error: Computed interpolation results "
        WRITE (*,904) " are incorrect"
        STOP 1
      ELSE
        WRITE (*,907) "","Computed interpolation results are correct"
      END IF
      STOP 0

  901 FORMAT (A,I0)
  902 FORMAT (/A)
  903 FORMAT (A,SP,F11.6,A,SP,F11.6)
  904 FORMAT (99A)
  905 FORMAT (A,SP,F11.6,A,I12,A,I12,A,SP,F11.6,A,SP,F11.6,              &
     & A,SP,F11.6)
  907 FORMAT (//A)
      END PROGRAM
