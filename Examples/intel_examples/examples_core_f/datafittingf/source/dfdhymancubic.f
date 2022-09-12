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
!    Construction of co-monotone cubic spline, accordingly Hyman algorithm
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
      INTEGER,PARAMETER :: N           = 11
      ! number of functions
      INTEGER,PARAMETER :: NNY         = 1
      ! size of array describing derivative orders to compute
      INTEGER,PARAMETER :: NNDERORDER  = 2
      ! number of derivatives to compute
      INTEGER,PARAMETER :: NNDER       = 2
      ! total number of spline coefficients
      INTEGER,PARAMETER :: NNSCOEFF    = NNY*(N-1)*DF_PP_CUBIC
      ! total number of interpolation sites
      INTEGER,PARAMETER :: NNSITE      = 10

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
      ! additional info about functions
      INTEGER :: yhint
      ! number of spline coefficients
      INTEGER :: nscoeff
      ! spline coefficients storage format
      INTEGER :: scoeffhint
      ! boundary conditions type
      INTEGER :: bc_type
      ! number of interpolation sites
      INTEGER :: nsite
      ! additional info about interpolation sites
      INTEGER :: sitehint
      ! size of array describing derivative orders
      INTEGER :: ndorder
      ! array of derivative orders
      INTEGER :: derorder(2)
      ! interpolation results storage format
      INTEGER :: rhint
      ! breakpoints
      REAL(8) :: x(N)
      ! function values
      REAL(8) :: y(NNY*N)
      ! array of spline coefficients
      REAL(8) :: scoeff(NNSCOEFF)
      ! array of interpolation sites
      REAL(8),TARGET :: site(NNSITE)
      ! spline evaluation results
      REAL(8),TARGET :: r(NNDER,NNSITE)
      ! boundary conditions
      REAL(8) :: bc(2)
      ! spline format
      INTEGER :: s_format
      ! method that is used to perform calculations
      INTEGER :: s_method
      ! type of calculations
      INTEGER :: intepr_type
      ! method that is used to perform calculations
      INTEGER :: intepr_method

      INTEGER :: xi,xxi
      INTEGER(4) :: errcode
      INTEGER :: comonot_status,comonot_status_1interval
      REAL(8) :: slope,di,dip1,alpha,beta,v1,v2

      errcode = 0

      ! **** Initializing parameters for Data Fitting task ****

      sorder = DF_PP_CUBIC
      stype  = DF_PP_HYMAN

      !***** Parameters describing construction type and method *****
      s_format = DF_PP_SPLINE
      s_method = DF_METHOD_STD

      !***** Parameters describing interpolation type and method *****
      intepr_type = DF_INTERP
      intepr_method = DF_METHOD_PP

      !***** Parameters describing interpolation interval *****
      nx    = N
      xhint = DF_NON_UNIFORM_PARTITION

      !***** Breakpoints *****
      x( 1) =  0.0D0
      x( 2) =  2.0D0
      x( 3) =  3.0D0
      x( 4) =  5.0D0
      x( 5) =  6.0D0
      x( 6) =  8.0D0
      x( 7) =  9.0D0
      x( 8) = 11.0D0
      x( 9) = 12.0D0
      x(10) = 14.0D0
      x(11) = 15.0D0

      !***** Function values *****
      y( 1) = 10.0D0
      y( 2) = 10.0D0
      y( 3) = 10.0D0
      y( 4) = 10.0D0
      y( 5) = 10.0D0
      y( 6) = 10.0D0
      y( 7) = 10.5D0
      y( 8) = 15.0D0
      y( 9) = 50.0D0
      y(10) = 60.0D0
      y(11) = 85.0D0

      !***** Parameters describing function *****
      ny = NNY
      yhint = DF_MATRIX_STORAGE_ROWS

      !***** Parameters describing spline coefficients storage *****
      nscoeff = NNSCOEFF
      scoeffhint = DF_NO_HINT

      !***** Parameters describing boundary conditions type *****
      bc_type = IOR(DF_BC_1ST_LEFT_DER , DF_BC_1ST_RIGHT_DER)
      bc = (/0.0D0,0.0D0/)

      !***** Parameters describing interpolation sites *****
      nsite      = NNSITE
      sitehint   = DF_NON_UNIFORM_PARTITION

      !***** Parameter describing interpolation results storage *****
      rhint = DF_MATRIX_STORAGE_ROWS

      !***** Parameter describing array for derivative orders *****
      ndorder = NNDERORDER
      ! spline values and 1st derivatives will be computed
      derorder = (/1,1/)

      !***** Number of derivatives to compute *****
      nder = NNDER


      !***** Create Data Fitting task *****
      errcode = dfdNewTask1D( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for natural cubic spline construction *****
      errcode = dfdEditPPSpline1D( task, sorder, stype, bc_type, bc=bc,  &
     & scoeff=scoeff, scoeffhint=scoeffhint )
      CALL CheckDfError(errcode)

      !***** Construct co-monotone cubic spline *****
      errcode = dfdConstruct1D( task, s_format, s_method )
      CALL CheckDfError(errcode)

      !***** Co-monotonicity test *****
      ! Test accordingly Fritch-Carlson criteria described in
      ! Fritsch, F. N.; Carlson, R. E. (1980). "Monotone Piecewise Cubic
      ! Interpolation". SIAM Journal on Numerical Analysis
      ! (SIAM) 17 (2): 238-246.
      comonot_status = 1
      DO xi = 1, nx-1
        ! For sub-interval number xi
        WRITE (*,901) "Interval ",xi
        WRITE (*,902) ", x=[",x(xi),";",x(xi+1),"]:"

        ! Interpolate on current interval
        ! and get 1st derivatives on its ends
        DO xxi = 1, nsite
            site(xxi) = x(xi) + (x(xi+1)-x(xi))/(nsite-1)*(xxi-1)
        END DO
        errcode = dfdInterpolate1D( task, intepr_type, intepr_method,    &
     &      nsite, site, sitehint, ndorder, derorder, r=r, rhint=rhint )
        CALL CheckDfError(errcode)

        comonot_status_1interval = 0

        !***** Co-monotonicity test on current sub-interval started *****
        ! Slope of the input data on current sub-interval
        slope = (y(xi+1)-y(xi)) / (x(xi+1)-x(xi))
        ! 1st derivative on left  side of current sub-interval
        di    = r(2,1)
        ! 1st derivative on right side of current sub-interval
        dip1  = r(2,nsite)

        ! Using Fritsch-Carlson criteria
        IF ( slope .EQ. 0.0 ) THEN
          ! Here, when input function is constant
          IF ( di .EQ. 0.0 .AND. dip1 .EQ. 0.0 ) THEN
            comonot_status_1interval = 1
          END IF
        ELSE
          ! Here, when input function is increasing or decreasing
          alpha = di   / slope
          beta  = dip1 / slope
          v1 = 2.0*alpha+beta-3.0
          v1 = v1 * v1
          v2 = 3.0*alpha*(alpha+beta-2.0)

          IF ( di    .GE. 0.0 .AND.                                      &
     &         dip1  .GE. 0.0 .AND.                                      &
     &         slope .GE. 0.0 .OR.                                       &
     &         di    .LE. 0.0 .AND.                                      &
     &         dip1  .LE. 0.0 .AND.                                      &
     &         slope .LE. 0.0 ) THEN
            IF ( alpha .GE. 0.0 .AND. alpha .LE. 3.0 .AND.               &
     &           beta  .GE. 0.0 .AND. beta  .LE. 3.0 .OR.                &
     &           v1 .LE. v2 ) THEN
              comonot_status_1interval = 1
            END IF
          END IF
        END IF
        !***** Co-monotonicity test on current sub-interval finished *****

        IF ( comonot_status_1interval == 1 ) THEN
          WRITE (*,*) " spline is co-monotone here."
        ELSE
          WRITE (*,*) " spline is not co-monotone here."
          WRITE (*,903) "  slope=",slope,", di=",di,", dip1=",dip1
          WRITE (*,904) "  alpha=",alpha,", beta=",beta

          DO xxi = 1, nsite
            WRITE (*,905) "  site(",xxi,")=",site(xxi)
            WRITE (*,906) " func=" ,r(1,xxi)
            WRITE (*,907) " deriv=",r(2,xxi)
          END DO
        END IF
        comonot_status = IAND(comonot_status,comonot_status_1interval)
      END DO


      !***** Delete Data Fitting task *****
      errcode = dfDeleteTask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      IF (comonot_status == 1) THEN
        WRITE (*,908) "Constructed spline"
        WRITE (*,*) " is co-monotone with input data."
      ELSE
        WRITE (*,908) "Error: constructed spline"
        WRITE (*,*) " is not co-monotone with input data."
        STOP 1
      END IF
      STOP 0


  901 FORMAT (A,I0,$)
  902 FORMAT (A,F9.6,A,F9.6,A,$)
  903 FORMAT (A,F9.6,A,F9.6,A,F9.6)
  904 FORMAT (A,F9.6,A,F9.6)
  905 FORMAT (A,I0,A,F9.6,$)
  906 FORMAT (A,F9.6,$)
  907 FORMAT (A,F9.6)
  908 FORMAT (A,$)

      END PROGRAM
