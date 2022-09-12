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
!    Construction of Subbotin spline Example Program Text
!*******************************************************************************

      include 'mkl_df.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "rescheck.inc"

      PROGRAM MKL_DF_TEST

      USE MKL_DF_TYPE
      USE MKL_DF
      USE DF_GENERATE_DATA
      USE DF_EXAMPLE_RESCHECK

      ! number of break points
      INTEGER,PARAMETER :: N       = 6
      ! number of functions
      INTEGER,PARAMETER :: NNY     = 1
      ! number of boundary conditions
      INTEGER,PARAMETER :: NNBC    = 2
      ! number of Subbotin spline knots
      INTEGER,PARAMETER :: NNIC    = (N+1)
      ! total number of spline coefficients
      INTEGER,PARAMETER :: NNSCOEFF=NNY*N*DF_PP_QUADRATIC

      ! left limit of interpolation interval
      REAL(8),PARAMETER :: LEFT_LIMIT  = -2.0d0
      ! right limit of interpolation interval
      REAL(8),PARAMETER :: RIGHT_LIMIT =  2.0d0

      REAL(8),PARAMETER :: FFREQ = 1.5d0

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
      ! additional info about spline coefficients
      INTEGER :: scoeffhint
      ! boundary conditions type
      INTEGER :: bc_type
      ! internal conditions type
      INTEGER :: ic_type
      ! number of internal conditions
      INTEGER :: nic
      ! array of break points
      REAL(8) :: x(N)
      ! array of Subbotin spline knots
      REAL(8) :: ic(NNIC)
      ! function values
      REAL(8) :: y(N*NNY)
      ! array of spline coefficients
      REAL(8) :: scoeff(NNSCOEFF)
      ! array of boundary conditions
      REAL(8) :: bc(NNBC)

      REAL(8) :: freq,left,right
      REAL(8) :: spline_val(N),left_val(N),right_val(N)
      REAL(8) :: left_der(N),right_der(N)
      INTEGER :: i,j,errnums
      INTEGER(4) :: errcode

      errcode = 0
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****

      sorder = DF_PP_QUADRATIC
      stype = DF_PP_SUBBOTIN

      !***** Parameters describing interpolation interval *****
      nx = N
      xhint = DF_NON_UNIFORM_PARTITION

      !***** Parameters describing function *****
      ny = NNY
      yhint = DF_NO_HINT

      !***** Parameters describing spline coefficients storage *****
      scoeffhint = DF_NO_HINT

      !***** Parameters describing boundary conditions type *****
      bc_type = IOR(DF_BC_1ST_LEFT_DER, DF_BC_1ST_RIGHT_DER)
      bc = (/1.0D0,-1.0D0/)

      !***** Parameters describing internal conditions type *****
      ic_type = DF_IC_Q_KNOT
      nic = NNIC

      !***** Generate uniformly distributed break points *****
      left = LEFT_LIMIT
      right = RIGHT_LIMIT
      errcode = dUniformRandSortedData( x, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Generate array of Subbotin spline knots *****
      ic(1)    = x(1)
      ic(nx+1) = x(nx)
      DO i = 1, nx-1
        ic(i+1) = ( x(i+1) + x(i) ) / 2.0
      END DO

      !***** Generate function y = sin(2 * Pi * freq * x) *****
      freq = FFREQ
      errcode = dSinDataNotUniformGrid( y, x, freq, nx )
      CALL CheckDfError(errcode)

      !***** Create Data Fitting task *****
      errcode = dfdnewtask1d( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for Subbotin spline coefficients
      !      computation *****
      errcode = dfdeditppspline1d( task, sorder, stype, bc_type, bc,     &
     &    ic_type, ic, scoeff, scoeffhint )
      CALL CheckDfError(errcode)

      !***** Construct Subbotin spline using STD method *****
      errcode = dfdconstruct1d( task, DF_PP_SPLINE, DF_METHOD_STD )
      CALL CheckDfError(errcode)

      !***** Delete Data Fitting task *****
      errcode = dfdeletetask( task )
      CALL CheckDfError(errcode)

      !***** Check computed coefficients *****

      !***** Check Subbotin spline values in break points *****
      errcode = dCheckSubbBreakPoints( nx, x, ny, y,                     &
     &    scoeff, spline_val )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Check spline values and 1st derivatives consistency
      !      in Subbotin spline knots *****
      errcode = dCheckSubbQNodes( nx, x, nic, ic, ny,                    &
     &    scoeff, left_val, right_val )
      IF ( errcode < 0 ) errnums = errnums+1

      errcode = dCheckQuadSubb1stDerConsistency( nic, x, ic, ny,         &
     &     scoeff, left_der, right_der )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Check boundary conditions *****
      errcode = dCheckQuadBC( nx, ic, ny,  scoeff, bc_type, bc )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Print results *****
      WRITE (*,901) "Number of break points : ",nx

      !***** Print given function *****
      WRITE (*,902) "    X             Y"
      DO j = 1, nx
        WRITE (*,903) " ",x(j),"   ",y(j)
      END DO

      !***** Print array of Subbotin spline knots *****
      WRITE (*,902) "  Subbotin spline knots:"
      DO j = 1, nic
        WRITE (*,904) " ",ic(j)
      END DO

      WRITE (*,902) "Boundary conditions(1st derivatives on both ends):"
      WRITE (*,903) " ",bc(1),"   ",bc(2)

      !***** Print computed spline coefficients *****
      WRITE (*,902) "Spline coefficients for Y :"
      WRITE (*,915) "    X^0           X^1           X^2       "
      WRITE (*,915) "    P(X[i])   "
      WRITE (*,915) "    P(Q[i])       P(Q[i+1]) "
      WRITE (*,905) "    P'(Q[i])      P'(Q[i+1])"
      DO j = 0, nx-1
        WRITE (*,906,ADVANCE='NO') " ",scoeff(sorder*j + 1),"   ",           &
     &  scoeff(sorder*j + 2),"   ",scoeff(sorder*j + 3),"   ",               &
     &  spline_val(j+1), "   ",right_val(j+1),"   ",left_val(j+1)
        WRITE (*,907) "   ",right_der(j+1),"   ",left_der(j+1)
      END DO

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,908) "","Error: Computed Subbotin spline coefficients"
        WRITE (*,905) " are incorrect"
        STOP 1
      ELSE
        WRITE (*,908) "","Computed Subbotin spline coefficients"
        WRITE (*,905) " are correct"
      END IF
      STOP 0

901   FORMAT (A,I0)
902   FORMAT (/A)
903   FORMAT (A,SP,F11.6,A,SP,F11.6)
904   FORMAT (A,SP,F11.6)
905   FORMAT (99A)
915   FORMAT (A,$)
906   FORMAT (A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,  &
     &        SP,F11.6)
907   FORMAT (A,SP,F11.6,A,SP,F11.6)
908   FORMAT (//A)

      END PROGRAM
