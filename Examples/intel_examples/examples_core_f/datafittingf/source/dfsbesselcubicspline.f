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
!   Calculation of Bessel cubic spline coefficients and spline evaluation
!   with user-defined extrapolation Example Program Text
!*******************************************************************************

      include 'mkl_df.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "rescheck.inc"

      MODULE MKL_DF_TEST_TYPE

       ! Structure containing linear extrapolation parameters
       TYPE DF_TEST_EXTRAP_PARAMS
          REAL(4) :: x
          REAL(4) :: y
          REAL(4) :: secant
       END TYPE DF_TEST_EXTRAP_PARAMS

       INTEGER,PARAMETER :: DF_TEST_EXTRAP_PARAMS_4BYTE_SIZE = 3

      END MODULE MKL_DF_TEST_TYPE


      PROGRAM MKL_DF_TEST

      USE MKL_DF_TYPE
      USE MKL_DF
      USE MKL_DF_TEST_TYPE
      USE DF_EXAMPLE_RESCHECK
      USE DF_GENERATE_DATA

      ! Data Fitting task descriptor
      TYPE (DF_TASK) task

      ! number of breakpoints
      INTEGER,PARAMETER :: N         =  7
      ! number of datasets to interpolate
      INTEGER,PARAMETER :: NNY       =  1
      ! number of interpolation sites
      INTEGER,PARAMETER :: NNSITE    =  10
      ! size of array describing derivative orders to compute
      INTEGER,PARAMETER :: NNDORDER  =  1
      ! number of spline coefficients
      INTEGER,PARAMETER :: NSCOEFF   =  NNY*(N-1)*DF_PP_CUBIC
      ! left  limit of interpolation interval
      REAL(4),PARAMETER :: LLIM_X    =  0.0
      ! right limit of interpolation interval
      REAL(4),PARAMETER :: RLIM_X    =  3.0
      ! left  limit of interpolation sites
      REAL(4),PARAMETER :: LLIM_SITE = -1.0
      ! right limit of interpolation sites
      REAL(4),PARAMETER :: RLIM_SITE =  4.0
      REAL(4),PARAMETER :: FFREQ     =  0.7

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
      ! functions storage format
      INTEGER :: yhint
      ! additional info about spline coefficients
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
      INTEGER :: dorder(1)
      ! interpolation results storage format
      INTEGER :: rhint
      ! indices of cells containing interpolation sites
      INTEGER,POINTER :: cell_idx
      ! construction method
      INTEGER :: method
      ! spline format
      INTEGER :: s_format
      ! type of calculations
      INTEGER :: intepr_type
      ! method that is used to perform calculations
      INTEGER :: intepr_method

      ! limits of interpolation interval
      REAL(4) :: left,right
      ! limits of interpolation sites
      REAL(4) :: lsite,rsite
      ! break points
      REAL(4) :: x(N)
      ! function values
      REAL(4) :: y(N*NNY)
      ! Bessel spline coefficients
      REAL(4) :: scoeff(NSCOEFF)
      ! array of interpolation sites
      REAL(4),TARGET :: site(NNSITE)
      ! spline evaluation results
      REAL(4) :: r(NNDORDER*NNSITE)
      ! interpolation call backs parameters
      INTEGER :: le_params(DF_TEST_EXTRAP_PARAMS_4BYTE_SIZE)
      INTEGER :: re_params(DF_TEST_EXTRAP_PARAMS_4BYTE_SIZE)
      TYPE (DF_TEST_EXTRAP_PARAMS) :: le_test_params,re_test_params

      REAL(4) :: left_val(NNY*(N-1)),right_val(NNY*(N-1))
      REAL(4) :: left_der1(NNY*(N-1)),right_der1(NNY*(N-1))
      REAL(4),TARGET :: ref_r(NNDORDER*NNSITE)
      REAL(4),POINTER :: y_cur,scoeff_cur
      REAL(4) :: freq
      INTEGER(8) :: tmp_cell(1)
      REAL(4),POINTER :: site_ptr(:)
      REAL(4),POINTER :: ref_r_ptr(:)

      INTEGER :: i,j,errnums
      INTEGER(4) :: errcode

      INTEGER(8) :: sites_num

      INTERFACE
        INTEGER FUNCTION linear_extrap(n,cell,site,r,params)
           INTEGER(8),INTENT(IN)  :: n
           INTEGER(8),INTENT(IN)  :: cell(*)
           INTEGER,INTENT(IN)  :: params(3)
           REAL(4),INTENT(IN)  :: site(*)
           REAL(4),INTENT(OUT) :: r(*)
        END FUNCTION
      END INTERFACE

      errcode = DF_STATUS_OK
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****

      !***** Parameters describing order and type of the spline *****
      sorder = DF_PP_CUBIC
      stype = DF_PP_BESSEL

      !***** Parameters describing interpolation interval *****
      left = LLIM_X
      right = RLIM_X
      nx = N
      xhint = DF_NON_UNIFORM_PARTITION

      !***** Parameters describing functions *****
      ny = NNY
      yhint = DF_MATRIX_STORAGE_ROWS

      !***** Parameters describing spline coefficients storage *****
      scoeffhint = DF_NO_HINT

      !***** Parameters describing boundary conditions type *****
      bc_type = DF_BC_NOT_A_KNOT

      !***** Parameters describing internal conditions type *****
      ic_type = DF_NO_IC

      !***** Parameters describing interpolation sites *****
      lsite = LLIM_SITE
      rsite = RLIM_SITE
      nsite = NNSITE
      sitehint = DF_SORTED_DATA

      !**** Parameter describing interpolation results storage *****
      rhint = DF_NO_HINT


      !**** Parameter describing array for derivative orders *****
      ndorder = NNDORDER
      dorder(1) = 1

      !***** Generate independent variables array with quasi-uniform
      !      break points *****
      errcode = sQuasiUniformData( x, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Generate interpolation sites *****
      errcode = sUniformRandSortedData( site, lsite, rsite, nsite )
      CALL CheckDfError(errcode)
      freq = FFREQ

      !***** Generate function y = sin(2 * Pi * freq * x) *****
      errcode = sSinDataNotUniformGrid( y, x, freq, nx )
      CALL CheckDfError(errcode)

      !***** Set left call back parameters ******
      le_test_params%x = x(1)
      le_test_params%y = y(1)
      le_test_params%secant = (y(2)  - y(1))    / (x(2)  - x(1))
      le_params = transfer(le_test_params, le_params)

      !***** Set right call back parameters ******
      re_test_params%x = x(nx-1)
      re_test_params%y = y(nx-1)
      re_test_params%secant = (y(nx) - y(nx-1)) / (x(nx) - x(nx-1))
      re_params = transfer(re_test_params, re_params)

      !***** Create Data Fitting task *****
      errcode = dfsnewtask1d( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for Bessel cubic spline coefficients
      !      computation *****
      errcode = dfseditppspline1d( task, sorder, stype, bc_type,         &
     & ic_type=ic_type, scoeff=scoeff, scoeffhint=scoeffhint )
      CALL CheckDfError(errcode)

      !***** Compute Bessel cubic spline coefficients using STD method *****
      s_format = DF_PP_SPLINE
      method = DF_METHOD_STD
      errcode = dfsconstruct1d( task, s_format, method )
      CALL CheckDfError(errcode)

      !***** Interpolate and use call backs for the left and right
      !      extrapolation *****
      intepr_type = DF_INTERP
      intepr_method = DF_METHOD_PP
      errcode = dfsinterpolateex1d( task, intepr_type, intepr_method,    &
     & nsite, site, sitehint, ndorder, dorder, r=r, rhint=rhint,         &
     & le_cb=linear_extrap, le_params=le_params, re_cb=linear_extrap,    &
     & re_params=re_params)
      CALL CheckDfError(errcode)

      !***** Check computed coefficients *****

      !***** Check spline values in break points *****
      errcode = sCheckCubBreakPoints( nx, x, ny, y, scoeff,              &
     &     left_val, right_val )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Check that spline 1st derivatives are equal for left
      !      and right piece of the spline for each break point *****
      errcode = sCheckCub1stDerConsistency( nx, x, ny, scoeff,           &
     &      left_der1, right_der1 )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Check computed interpolation results *****
      errcode = sCheckCubInterpRes( nx, x, ny, scoeff, nsite,            &
     & site, ndorder, dorder, r, ref_r )
      IF ( errcode < 0 ) errnums = errnums+1

      sites_num = 1

      !***** Check computed extrapolation results *****
      DO i = 1, nsite
        IF ( site(i) < x(1) ) THEN
          ! Check left  extrapolation results
          site_ptr  => site(i:)
          ref_r_ptr => ref_r(i:)

          errcode = linear_extrap( sites_num, tmp_cell, site_ptr,        &
     &                             ref_r_ptr, le_params )
          dif = ref_r(i) - r(i)
          IF ( dif < 0 ) dif = -dif
          IF ( dif > EPSILON_DOUBLE ) errnums = errnums+1
        END IF
        IF ( site(i) > x(nx) ) THEN
          ! Check right extrapolation results
          site_ptr  => site(i:)
          ref_r_ptr => ref_r(i:)

          errcode = linear_extrap( sites_num, tmp_cell, site_ptr,        &
     &                             ref_r_ptr, re_params )
          dif = ref_r(i) - r(i)
          IF ( dif < 0 ) dif = -dif
          IF ( dif > EPSILON_DOUBLE ) errnums = errnums+1
        END IF
      END DO

      !***** Print results *****
      WRITE (*,901) "Number of break points : ",nx

      !***** Print given tabular function *****
      WRITE (*,902) " i    x(i)          y(i)"
      DO j = 1, nx
        WRITE (*,903) " ",j," ",x(j),"   ",y(j)
      END DO

      !***** Print computed spline coefficients *****
      WRITE (*,904) "Coefficients are calculated for a polynomial of     &
     & the form:",""
      WRITE (*,905) "Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 +      &
     & Di*(x - x(i))^3"
      WRITE (*,905) "    where x(i) <= x < x(i+1)"
      WRITE (*,902) "Spline coefficients for Y:"
      WRITE (*,905) " i    Ai            Bi            Ci            Di  &
     &           Pi(x(i))      Pi(x(i+1))    Pi'(x(i))    Pi'(x(i+1))"
      DO j = 1, nx-1
        WRITE (*,906,ADVANCE='NO') " ",j," ",scoeff(sorder*(j-1) + 1),   &
     &  "   ",scoeff(sorder*(j-1) + 2),"   ",scoeff(sorder*(j-1) + 3),   &
     &  "   ",scoeff(sorder*(j-1) + 4),"   ",right_val(j),"   ",         &
     &  left_val(j)
        WRITE (*,907) "   ",right_der1(j),"   ",left_der1(j)
      END DO

      !***** Print interpolation results ******
      WRITE (*,902) "Results of interpolation:"
      WRITE (*,905) "    Sites         Obtained     Expected"
      DO i = 1, nsite
        WRITE (*,908) " ",site(i),"  ",r(i),"  ",ref_r(i)
      END DO

      !***** Delete Data Fitting task *****
      errcode = dfDeleteTask( task )
      CALL CheckDfError(errcode)

      !***** Printing summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,909) "","Error: Computed Bessel cubic spline            &
     &   coefficients or spline values are incorrect"
        STOP 1
      ELSE
        WRITE (*,909) "","Computed Bessel cubic spline coefficients"
        WRITE (*,905) " and spline values are correct"
      END IF
      STOP 0

901   FORMAT (A,I0)
902   FORMAT (/A)
903   FORMAT (A,I0,A,SP,F11.6,A,SP,F11.6)
904   FORMAT (/A/A)
905   FORMAT (99A)
906   FORMAT (A,I1,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,     &
     &        F11.6,A,SP,F11.6)
907   FORMAT (A,SP,F11.6,A,SP,F11.6)
908   FORMAT (A,SP,F11.6,A,SP,F11.6,A,SP,F11.6)
909   FORMAT (//A)

      END PROGRAM

!*******************************************************************************
!   Call back for linear extrapolation.
!
! API
!  FUNCTION linear_extrap(n,cell,site,r,params)  RESULT (errcode)
!
! INPUT PARAMETERS:
!  n      - number of interpolation sites
!  cell   - array of size n containing indices of the cells to which the
!           interpolation sites in array 'site' belong
!  site   - array of size n that holds the interpolation sites
!  params - structure contatining extrapolation parameters
!
! OUTPUT PARAMETERS:
!  r      - array of integration results
!
! RETURN VALUE:
!  The status returned by the callback function:
!   0  - indicates successful completion of the callback operation
!   <0 - error
!   >0 - warning
!*******************************************************************************
      FUNCTION linear_extrap(n,cell,site,r,params)  RESULT (errcode)

        USE MKL_DF_TYPE
        USE MKL_DF_TEST_TYPE

        INTEGER :: errcode
        INTEGER(8),INTENT(IN) :: n
        INTEGER(8),INTENT(IN) :: cell(*)
        INTEGER,INTENT(IN) :: params(DF_TEST_EXTRAP_PARAMS_4BYTE_SIZE)
        REAL(4),INTENT(IN) :: site(*)
        REAL(4),INTENT(OUT) :: r(*)

        INTEGER(8) :: i
        TYPE(DF_TEST_EXTRAP_PARAMS) :: p

        p = transfer(params, p)

        DO i = 1, n
          r(i) = p%y + p%secant * ( site(i) - p%x )
        END DO
        errcode = 0
        RETURN

      END FUNCTION
