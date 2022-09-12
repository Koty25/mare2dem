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
!    Construction and evaluation of quadratic spline Example Program Text
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
      INTEGER,PARAMETER :: N          = 9
      ! number of functions
      INTEGER,PARAMETER :: NNY        = 1
      ! number of pairs of integration limits
      INTEGER,PARAMETER :: NNSITE     = 10
      ! size of array describing derivative orders to compute
      INTEGER,PARAMETER :: NNDORDER   = 1
      ! total number of spline coefficients
      INTEGER,PARAMETER :: NNSCOEFF   = (NNY*(N-1)*DF_PP_QUADRATIC)

      ! left limit of interpolation interval
      REAL(8),PARAMETER :: LLIM_X    = -2.0d0
      ! right limit of interpolation interval
      REAL(8),PARAMETER :: RLIM_X    =  2.0d0
      ! left limit of interpolation sites
      REAL(8),PARAMETER :: LLIM_SITE = -3.0d0
      ! right limit of interpolation site
      REAL(8),PARAMETER :: RLIM_SITE =  3.0d0
      REAL(8),PARAMETER :: FFREQ     =  0.5d0

      ! type of calculations
      INTEGER :: type
      ! method that is used to perform calculations
      INTEGER :: method

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
      ! additional info about function
      INTEGER :: yhint
      ! spline coefficients storage format
      INTEGER :: scoeffhint
      ! boundary conditions type
      INTEGER :: bc_type
      ! internal conditions type
      INTEGER :: ic_type
      ! number of interpolation sites
      INTEGER :: nsite
      ! additional info about interpolation sites
      INTEGER :: sitehint
      ! size of array describing derivative orders
      INTEGER :: ndorder
      ! array describing derivative orders
      INTEGER :: dorder(NNDORDER)
      ! interpolation results storage format
      INTEGER :: rhint

      ! array of break points
      REAL(8) :: x(N)
      ! function values
      REAL(8) :: y(NNY*N)
      ! boundary condition
      REAL(8) :: bc(1)
      ! array of spline coefficients
      REAL(8) :: scoeff(NNSCOEFF)
      ! limits of interpolation sites
      REAL(8) :: site(2)
      ! array of interpolation sites in full format
      REAL(8) :: site_full(NNSITE)
      ! spline evaluation results
      REAL(8) :: r(NNDORDER*NNSITE)

      REAL(8) :: left,right
      REAL(8) :: freq
      REAL(8) :: left_val(N-1), right_val(N-1)
      REAL(8) :: left_der(N-1), right_der(N-1)
      REAL(8) :: ref_r(NNDORDER*NNSITE)

      INTEGER :: i,j,errnums
      INTEGER(4) :: errcode
      errcode = DF_STATUS_OK
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****
      sorder     = DF_PP_QUADRATIC
      stype      = DF_PP_DEFAULT

      !***** Parameters describing interpolation interval *****
      left  = LLIM_X
      right = RLIM_X
      nx         = N
      xhint      = DF_NON_UNIFORM_PARTITION

      !***** Parameters describing function *****
      ny         = NNY
      yhint      = DF_NO_HINT

      !***** Parameters describing spline coefficients storage *****
      scoeffhint = DF_NO_HINT

      !***** Parameters describing boundary conditions type *****
      bc_type    = DF_BC_Q_VAL
      bc(1)      = 1.0

      !***** Parameters describing internal conditions type *****
      ! No internal conditions are provided for quadratic spline
      ic_type    = DF_NO_IC

      !***** Parameters describing interpolation sites *****
      nsite      = NNSITE
      sitehint   = DF_UNIFORM_PARTITION
      ! Limits of interpolation interval are provided if the sites
      ! are uniform
      site(1)    = LLIM_SITE
      site(2)    = RLIM_SITE

      !**** Parameter describing interpolation results storage *****
      rhint      = DF_NO_HINT

      !**** Parameters describing derivative orders *****
      ndorder    = NNDORDER
      dorder(1)  = 1

      !***** Generate uniformly distributed break points *****
      errcode = dUniformRandSortedData( x, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Generate function y = sin(2 * Pi * freq * x) *****
      freq = FFREQ
      errcode = dSinDataNotUniformGrid( y, x, freq, nx )
      CALL CheckDfError(errcode)

      !***** Create Data Fitting task *****
      errcode = dfdnewtask1d( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for quadratic spline construction *****
      errcode = dfdeditppspline1d( task, sorder, stype, bc_type, bc,     &
     &     ic_type, scoeff=scoeff, scoeffhint=scoeffhint )
      CALL CheckDfError(errcode)

      !***** Construct quadratic spline using STD method *****
      errcode = dfdconstruct1d( task, DF_PP_SPLINE, DF_METHOD_STD )
      CALL CheckDfError(errcode)

      !***** Interpolate using PP method *****
      type = DF_INTERP
      method = DF_METHOD_PP
      errcode = dfdinterpolate1d( task, type, method, nsite,             &
     &      site, sitehint, ndorder, dorder, r=r, rhint=rhint )
      CALL CheckDfError(errcode)

      !***** Check computed coefficients *****/

      !***** Check spline values in break points *****/
      errcode = dCheckQuadBreakPoints( nx, x, ny, y, scoeff,             &
     &      left_val, right_val )
      IF ( errcode < 0 ) errnums = errnums + 1

      !***** Check 1st derivatives in break points *****/
      errcode = dCheckQuad1stDerConsistency( nx, x, ny, scoeff,          &
     &      left_der, right_der )
      IF ( errcode < 0 ) errnums = errnums + 1

      !***** Check boundary conditions *****
      errcode = dCheckQuadBC( nx, x, ny, scoeff, bc_type, bc )
      IF ( errcode < 0 ) errnums = errnums + 1

      !***** Check results of interpolation *****/
      errcode = dUniformData( site_full, LLIM_SITE, RLIM_SITE, nsite )
      CALL CheckDfError(errcode)

      errcode = dCheckQuadInterpRes( nx, x, ny, scoeff,                  &
     &                               nsite, site_full, ndorder, dorder,  &
     &                               r, ref_r )
      CALL CheckDfError(errcode)

      !***** Print results *****
      WRITE (*,901) "Number of break points : ",nx

      !***** Print given function *****
      WRITE (*,902) " i    x(i)          y(i)"
      DO j = 1, nx
        WRITE (*,907) " ",j," ",x(j),"   ",y(j)
      END DO

      !***** Print computed spline coefficients *****
      WRITE (*,904) "Coefficients are calculated for a polynomial of the &
     & form:",""
      WRITE (*,902) "Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2"
      WRITE (*,904) "    where x(i) <= x < x(i+1)",""
      WRITE (*,902) "Spline coefficients for Y:"
      WRITE (*,902) " i    Ai          Bi          Ci          P(x(i))   &
     &   P(x(i+1))   P'(x(i))    P'(x(i+1))"

      DO j = 1, nx-1
        WRITE (*,910) " ", j, " ", scoeff(sorder*(j-1) + 1), " ",        &
     &      scoeff(sorder*(j-1) + 2), " ",  scoeff(sorder*(j-1) + 3)
        WRITE (*,911) " ",right_val(j)," ",left_val(j)," ",              &
     &      right_der(j)," ",left_der(j)
      END DO

      !***** Print interpolation results ******
      WRITE (*,902) "Results of interpolation:"
      PRINT  *,     "    Sites         Function value"
      PRINT  *,     "               Obtained    Expected"
      DO j = 1, nsite
        WRITE (*,909) " ", site_full(j), " ", r(j), " ", ref_r(j)
      END DO

      !***** Delete Data Fitting task *****
      errcode = dfdeletetask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,904) "","Error: Computed quadratic spline coefficients  &
     & or values are incorrect"
        STOP 1
      ELSE
        WRITE (*,904) "","Computed quadratic spline coefficients and     &
     & values are correct"
      END IF
      STOP 0

901   FORMAT (A,I0)
902   FORMAT (/A)
904   FORMAT (/A/A)
907   FORMAT (A,I1,A,SP,F11.6,A,SP,F11.6)
909   FORMAT (A,SP,F11.6,A,SP,F11.6,A,SP,F11.6)
910   FORMAT (A,I1,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,$)
911   FORMAT (A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6)

      END PROGRAM
