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
!    Construction of cubic spline with given second derivative
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
      INTEGER,PARAMETER :: N    = 7
      ! number of functions
      INTEGER,PARAMETER :: NNY  = 1
      ! number of internal conditions
      INTEGER,PARAMETER :: NNIC = N-2
      ! number of boundary conditions
      INTEGER,PARAMETER :: NNBC = 2

      ! total number of spline coefficients
      INTEGER,PARAMETER :: NNSCOEFF=(NNY*(N-1)*DF_PP_CUBIC)

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
      ! number of spline coefficients
      INTEGER :: nscoeff
      ! spline coefficients storage format
      INTEGER :: scoeffhint
      ! boundary conditions type
      INTEGER :: bc_type
      ! internal conditions type
      INTEGER :: ic_type
      ! number of internal conditions
      INTEGER :: nic
      ! array of break points
      REAL(8),TARGET  :: x(N)
      REAL(8),POINTER :: x_ptr(:)
      ! function values
      REAL(8) :: y(NNY*N)
      ! array of internal conditions
      REAL(8) :: ic(NNIC)
      ! array of boundary conditions
      REAL(8),TARGET  :: bc(NNBC)
      REAL(8),POINTER :: bc_ptr(:)
      ! array of spline coefficients
      REAL(8) :: scoeff(NNSCOEFF)
      ! limits of the interpolation interval
      REAL(8) :: left,right
      ! type of calculations
      INTEGER :: type
      ! method that is used to perform calculations
      INTEGER :: method

      REAL(8) :: freq
      REAL(8) :: left_val(N-1),right_val(N-1)
      REAL(8) :: left_der2(N-1),right_der2(N-1)
      INTEGER :: i,j
      INTEGER :: errnums
      INTEGER(4) :: errcode

      errcode = 0
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****

      sorder = DF_PP_CUBIC
      stype  = DF_PP_DEFAULT

      !***** Parameters describing interpolation interval *****
      left = LEFT_LIMIT
      right = RIGHT_LIMIT
      nx = N
      xhint = DF_NON_UNIFORM_PARTITION

      !***** Parameters describing function *****
      ny = NNY

      !***** Parameters describing spline coefficients storage *****
      nscoeff = NNSCOEFF
      scoeffhint = DF_MATRIX_STORAGE_ROWS

      !***** Parameters describing boundary conditions type *****
      bc_type = IOR(DF_BC_2ND_LEFT_DER, DF_BC_2ND_RIGHT_DER)
      nbc = NNBC

      !***** Parameters describing internal conditions type *****
      ic_type = DF_IC_2ND_DER
      nic = NNIC

      !***** Generate array of uniformly distributed break points *****
      errcode = dUniformRandSortedData( x, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Generate function y = sin(2 * Pi * freq * x) *****
      freq = FFREQ
      errcode = dSinDataNotUniformGrid( y, x, freq , nx )
      CALL CheckDfError(errcode)

      !***** Generate internal conditions needed for the spline construction *****
      x_ptr => x(2:)
      errcode = dSinDer2DataNotUniformGrid( ic, x_ptr, freq , nic )
      CALL CheckDfError(errcode)

      !***** Generate boundary conditions *****
      x_ptr => x(1:)
      bc_ptr => bc(1:)
      errcode = dSinDer2DataNotUniformGrid( bc_ptr, x_ptr, freq , 1 )
      CALL CheckDfError(errcode)
      x_ptr => x(nx:)
      bc_ptr => bc(2:)
      errcode = dSinDer2DataNotUniformGrid( bc_ptr, x_ptr, freq , 1 )
      CALL CheckDfError(errcode)

      !***** Create Data Fitting task *****
      errcode = dfdnewtask1d( task, nx, x, xhint, ny, y )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for cubic spline with provided
      !      2nd derivatives construction *****
      errcode = dfdeditppspline1d( task, sorder, stype, bc_type, bc,     &
     &    ic_type, ic, scoeff, scoeffhint )
      CALL CheckDfError(errcode)

      !***** Construct cubic spline with provided 2nd derivatives
      !      using STD method *****
      type = DF_PP_SPLINE
      method = DF_METHOD_STD

      errcode = dfdconstruct1d( task, type, method )
      CALL CheckDfError(errcode)

      !***** Check computed coefficients *****

      !***** Check spline values in break points *****
      errcode = dCheckCubBreakPoints( nx, x, ny, y, scoeff,              &
     &    left_val, right_val )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Check that spline 2nd derivatives are equal for left and right
      !      piece of the spline for each break point *****
      errcode = dCheckCub2ndDerConsistency( nx, x, ny, scoeff,           &
     &    left_der2, right_der2 )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Check internal conditions *****
      DO j = 1, nic
        IF ( ABS( ic(j) - left_der2(j) ) > EPSILON_DOUBLE )              &
     &       errnums = errnums+1
      END DO

      !***** Check boundary conditions *****
      errcode = dCheckCubBC( nx, x, ny, scoeff, bc_type, bc )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Print results *****
      WRITE (*,901) "Number of break points : ",nx

      !***** Print given function *****
      WRITE (*,902) "    X             Y(X)          Y''(X)"
      WRITE (*,903) " ",x(1),"   ",y(1),"   ",bc(1)
      DO j = 2, nx-1
        WRITE (*,903) " ",x(j),"   ",y(j),"   ",ic(j-1)
      END DO
      WRITE (*,903) " ",x(nx),"   ",y(nx),"   ",bc(2)

      !***** Print computed spline coefficients *****
      WRITE (*,904) "Coefficients are calculated for a polynomial of the &
     & form:",""
      WRITE (*,905) "Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 + Di*( &
     & x - x(i))^3"
      WRITE (*,905) "    where x(i) <= x < x(i+1)"
      WRITE (*,902) "Spline coefficients for Y:"
      WRITE (*,915) " i    Ai            Bi            Ci            Di"
      WRITE (*,915) "            P(X[i])       P(X[i+1]) "
      WRITE (*,905) "    P''(X[i])     P''(X[i+1])"
      DO j = 1, nx-1
        WRITE (*,906,ADVANCE='NO') " ",j," ",scoeff(sorder*(j-1) + 1),   &
     &   "   ",scoeff(sorder*(j-1) + 2),"   ",scoeff(sorder*(j-1) + 3),  &
     &   "   ",scoeff(sorder*(j-1) + 4),"   ",right_val(j),"   ",        &
     &   left_val(j)
        WRITE (*,907) "   ",right_der2(j),"   ",left_der2(j)
      END DO

      !***** Delete Data Fitting task *****
      errcode = dfdeletetask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      WRITE (*,905) ""
      IF (errnums /= 0) THEN
        WRITE (*,915) "Error: "
        WRITE (*,915) "Computed default cubic spline coefficients"
        WRITE (*,905) " are incorrect"
        STOP 1
      ELSE
        WRITE (*,915) "Computed default cubic spline coefficients"
        WRITE (*,905) " are correct"
      END IF
      STOP 0

901   FORMAT (A,I0)
902   FORMAT (/A)
903   FORMAT (A,SP,F11.6,A,SP,F11.6,A,SP,F11.6)
904   FORMAT (/A/A)
905   FORMAT (99A)
915   FORMAT (A,$)
906   FORMAT (A,I1,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,     &
     &        F11.6,A,SP,F11.6)
907   FORMAT (A,SP,F11.6,A,SP,F11.6)
908   FORMAT (//A)
918   FORMAT (A,$)

      END PROGRAM
