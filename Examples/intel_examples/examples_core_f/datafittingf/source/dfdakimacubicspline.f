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
!    Calculation of Akima cubic spline coefficients and integral
!    computation with callback function  Example Program Text
!******************************************************************************/


      include 'mkl_df.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "rescheck.inc"

      MODULE MKL_DF_TEST_VARS

        ! Number of pairs of integration limits
        INTEGER,PARAMETER   :: NNLIM  = 2
        ! size of array describing derivative orders to compute
        INTEGER,PARAMETER   :: NNDORDER  =  1

        ! Integration intervals limits
        REAL(8)     :: llim(NNLIM),rlim(NNLIM)
        ! Corresponding cell search results
        INTEGER     :: llim_cell(NNLIM),rlim_cell(NNLIM)
        INTEGER(8)  :: llim_cell8(NNLIM),rlim_cell8(NNLIM)

      END MODULE MKL_DF_TEST_VARS




      MODULE MKL_DF_TEST_TYPES

        USE MKL_DF_TYPE

        TYPE MY_PARAMS
          REAL(8) :: v1
          REAL(8) :: v2
          REAL(8) :: v3
        END TYPE

        ! Size of MY_PARAMS in 32-bit blocks
        INTEGER,PARAMETER :: MKL_DF_TEST_MY_PARAMS_SIZE = 6

        INTEGER,PARAMETER :: DF_TEST_PARAMS_4BYTE_SIZE = 2

      END MODULE MKL_DF_TEST_TYPES




      MODULE MKL_DF_TEST_FUNCS

        INTERFACE
          INTEGER FUNCTION integr_search_cb(n,site,cell,flag,            &
     &          user_params,libparams)
            USE MKL_DF_TYPE
            USE MKL_DF_TEST_TYPES
            INTEGER(8),INTENT(IN)  :: n
            REAL(8),INTENT(IN)     :: site(*)
            INTEGER(8),INTENT(OUT) :: cell(*)
            INTEGER(4),INTENT(OUT) :: flag(*)
            TYPE(MY_PARAMS),INTENT(IN) :: user_params
            TYPE(DF_SEARCH_CALLBACK_LIBRARY_PARAMS),INTENT(IN)           &
     &                                 :: libparams
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION left_akima_integr(n,lcell,llim,rcell,rlim,    &
     &          r,x,libparams)
            USE MKL_DF_TYPE
           INTEGER(8),INTENT(IN)  :: n
           INTEGER(8),INTENT(IN)  :: lcell(*)
           REAL(8),INTENT(IN)     :: llim(*)
           INTEGER(8),INTENT(IN)  :: rcell(*)
           REAL(8),INTENT(IN)     :: rlim(*)
           REAL(8),INTENT(OUT)    :: r(*)
            INTEGER,INTENT(IN),OPTIONAL :: x(*)
            TYPE(DF_INTEGR_CALLBACK_LIBRARY_PARAMS),INTENT(IN),          &
     &      OPTIONAL                    :: libparams
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION right_akima_integr(n,lcell,llim,rcell,rlim,   &
     &          r,x,libparams)
            USE MKL_DF_TYPE
           INTEGER(8),INTENT(IN)  :: n
           INTEGER(8),INTENT(IN)  :: lcell(*)
           REAL(8),INTENT(IN)     :: llim(*)
           INTEGER(8),INTENT(IN)  :: rcell(*)
           REAL(8),INTENT(IN)     :: rlim(*)
           REAL(8),INTENT(OUT)    :: r(*)
            INTEGER,INTENT(IN),OPTIONAL :: x(*)
            TYPE(DF_INTEGR_CALLBACK_LIBRARY_PARAMS),INTENT(IN),          &
     &      OPTIONAL                    :: libparams
          END FUNCTION
        END INTERFACE

      END MODULE MKL_DF_TEST_FUNCS




      PROGRAM MKL_DF_TEST

      USE MKL_DF_TYPE
      USE MKL_DF
      USE DF_GENERATE_DATA
      USE DF_EXAMPLE_RESCHECK
      USE MKL_DF_TEST_TYPES
      USE MKL_DF_TEST_FUNCS
      USE MKL_DF_TEST_VARS

      ! Data Fitting task descriptor
      TYPE (DF_TASK) task

      ! number of breakpoints
      INTEGER,PARAMETER :: N           =  10
      ! number of datasets to integrate
      INTEGER,PARAMETER :: NNY         =  2
      INTEGER,PARAMETER :: NNSCOEFF    =  (N-1)*DF_PP_CUBIC*NNY
      ! left  limit of interpolation interval
      REAL(8),PARAMETER :: LLIM_X      =  0.0d0
      ! right limit of interpolation interval
      REAL(8),PARAMETER :: RLIM_X      =  3.0d0
      ! left  limit of integration interval
      REAL(8),PARAMETER :: LLIM_INTEGR = -0.5d0
      ! right limit of integration interval
      REAL(8),PARAMETER :: RLIM_INTEGR =  3.5d0
      REAL(8),PARAMETER :: FFREQ       =  0.75d0

      ! spline order
      INTEGER :: sorder
      ! spline type
      INTEGER :: stype
      ! number of break points
      INTEGER :: nx
      ! additional info about break points
      INTEGER :: xhint
      ! additional info about sites
      INTEGER :: sitehint
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
      ! number of pairs of integration limits
      INTEGER :: nlim
      ! integration limits storage formats
      INTEGER :: llimhint,rlimhint
      ! integration results storage format
      INTEGER :: rhint
      ! break points
      REAL(8) :: x(N)
      ! function values
      REAL(8) :: y(N*NNY)
      ! boundary conditions
      REAL(8) :: bc(2)

      ! Akima spline coefficients
      REAL(8) :: scoeff(NNSCOEFF)
      ! limits of interpolation interval
      REAL(8) :: left,right
      ! integration results
      REAL(8) :: r(NNLIM*NNY)

      INTEGER :: le_params(DF_TEST_PARAMS_4BYTE_SIZE)
      INTEGER :: re_params(DF_TEST_PARAMS_4BYTE_SIZE)

      REAL(8) :: left_val(NNY*(N-1)),right_val(NNY*(N-1))
      REAL(8) :: left_der1(NNY*(N-1)),right_der1(NNY*(N-1))
      REAL(8) :: freq,r_ref

      INTEGER :: i,j,errnums
      INTEGER(4) :: errcode

      INTEGER :: s_params(MKL_DF_TEST_MY_PARAMS_SIZE)
      TYPE(MY_PARAMS) :: s_p

      errcode = DF_STATUS_OK
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****

      !***** Parameters describing order and type of the spline *****
      sorder = DF_PP_CUBIC
      stype  = DF_PP_AKIMA

      !***** Parameters describing interpolation interval *****
      left  = LLIM_X
      right = RLIM_X
      nx    = N
      xhint = DF_NON_UNIFORM_PARTITION
      sitehint = DF_NON_UNIFORM_PARTITION

      !***** Parameters describing function *****
      ny = NNY
      yhint = DF_NO_HINT

      !***** Parameter describing additional info about spline
      !      coefficients *****
      scoeffhint = DF_NO_HINT

      !***** Parameters describing boundary conditions type *****
      bc_type = IOR(DF_BC_2ND_LEFT_DER, DF_BC_1ST_RIGHT_DER)
      bc = (/0.0D0,-0.5D0/)

      !***** Parameters describing internal conditions type *****
      ! No internal conditions are provided for Akima cubic spline
      ic_type = DF_NO_IC

      !***** Parameters decsribing integration limits *****
      nlim = NNLIM
      llim(1) = LLIM_INTEGR
      llim(2) = LLIM_INTEGR+0.5
      rlim(1) = RLIM_INTEGR
      rlim(2) = RLIM_INTEGR+0.5
      llimhint = DF_NO_HINT
      rlimhint = DF_NO_HINT

      !***** Parameter dascribing integration results storage format *****
      rhint = DF_NO_HINT

      !***** Generate partition with uniformly distributed break points *****
      errcode = dUniformRandSortedData( x, left, right, nx )
      CALL CheckDfError(errcode)

      le_params = transfer(left, le_params)
      re_params = transfer(right, re_params)

      !***** Set user parameters for search callback *****
      s_p%v1 = 1.0
      s_p%v2 = 2.0
      s_p%v3 = 3.0
      !***** ... and copy them to array at address s_params *****
      s_params=transfer(s_p,s_params)

      !***** Generate function y1 = sin(2 * Pi * freq * x) *****
      freq = FFREQ
      errcode = dSinDataNotUniformGrid( y(1:N), x, freq, nx )
      CALL CheckDfError(errcode)
      !***** Generate function y2 = sin(2 * Pi * freq * x) *****
      freq = FFREQ*2.0
      errcode = dSinDataNotUniformGrid( y(N+1:2*N), x, freq, nx )
      CALL CheckDfError(errcode)

      !***** Create Data Fitting task *****
      errcode = dfdnewtask1d( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Search for left integration limits, save result to llim_cell *****
      errcode = dfdsearchcells1d( task, DF_METHOD_STD, NLIM,llim,       &
     &      sitehint, cell=llim_cell )
      CALL CheckDfError(errcode)

      !***** Search for right integration limits, save result to rlim_cell *****
      errcode = dfdsearchcells1d( task, DF_METHOD_STD, NLIM,rlim,       &
     &      sitehint, cell=rlim_cell )
      CALL CheckDfError(errcode)

        llim_cell8 = llim_cell
        rlim_cell8 = rlim_cell


      !***** Edit task parameters for Akima cubic spline construction *****
      errcode = dfdeditppspline1d( task, sorder, stype, bc_type, bc,     &
     &     ic_type, scoeff=scoeff, scoeffhint=scoeffhint )
      CALL CheckDfError(errcode)

      !***** Construct Akima cubic spline using STD method *****
      errcode = dfdconstruct1d( task, DF_PP_SPLINE, DF_METHOD_STD )
      CALL CheckDfError(errcode)

      !***** Compute integral for the spline on the interval (llim, rlim) *****
      errcode = dfdintegrateex1d( task, DF_METHOD_PP, nlim, llim,        &
     &   llimhint, rlim, rlimhint, r=r, rhint=rhint,                     &
     &   le_cb=left_akima_integr, le_params=le_params,                   &
     &   re_cb=right_akima_integr, re_params=re_params,                  &
     &   search_cb=integr_search_cb, search_params=s_params              &
     &   )
      CALL CheckDfError(errcode)

      !***** Check computed coefficients *****

      !***** Check spline values in break points *****
      errcode = dCheckCubBreakPoints( nx, x, ny, y, scoeff,              &
     &    left_val, right_val )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Check that spline 1st derivatives are equal for left
      !      and right piece of the spline for each break point *****
      errcode = dCheckCub1stDerConsistency( nx, x, ny, scoeff,           &
     &    left_der1, right_der1 )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Print results *****

      !***** Print given function *****
      WRITE (*,901) "Number of break points : ",nx
      WRITE (*,913) " i    x(i)"
      DO i = 1, ny
        WRITE (*,911) "          y",i
        WRITE (*,913) "(i)"
      END DO
      WRITE (*,910) ""
      DO j = 1, nx
        WRITE (*,911) " ",j
        WRITE (*,912) " ",x(j)
        DO i = 1, ny
          WRITE (*,912) "   ",y(j+(i-1)*N)
        END DO
        WRITE (*,910) ""
      END DO

      !***** Print computed spline coefficients *****
      WRITE (*,904) "Coefficients are calculated for a polynomial of the &
     & form:",""
      WRITE (*,905) "Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 +      &
     & Di*(x - x(i))^3"
      WRITE (*,905) "    where x(i) <= x < x(i+1)"
      WRITE (*,902) "Spline coefficients for Y:"
      WRITE (*,915) " i      Ai              Bi              Ci          &
     &     Di      "
      WRITE (*,905) "      Pi(x(i))      Pi(x(i+1))    Pi'(x(i))         &
     & Pi'(x(i+1))"

      DO i = 1, ny
        WRITE (*,911) "   for function y",i
        WRITE (*,910) ""
        DO j = 1, nx-1
          WRITE (*,906,ADVANCE='NO') " ",j,                              &
     &      " "  ,scoeff((i-1)*(nx-1)*sorder + sorder*(j-1) + 1),        &
     &      "   ",scoeff((i-1)*(nx-1)*sorder + sorder*(j-1) + 2),        &
     &      "   ",scoeff((i-1)*(nx-1)*sorder + sorder*(j-1) + 3),        &
     &      "   ",scoeff((i-1)*(nx-1)*sorder + sorder*(j-1) + 4),        &
     &      "   ",right_val((i-1)*(nx-1) + j),                           &
     &      "   ",left_val( (i-1)*(nx-1) + j)
          WRITE (*,907) "   ",right_der1((i-1)*(nx-1) + j),              &
     &      "   ",left_der1((i-1)*(nx-1) + j)
        END DO
      END DO

      !***** Print computed integration results *****
      WRITE (*,911) "Spline-based integrals for ",ny
      WRITE (*,911)     " functions on ",nlim
      WRITE (*,910)     " intervals:"
      DO j = 1, nlim
        WRITE (*,911) "interval_",j
        WRITE (*,912)   "=[ ",llim(j)
        WRITE (*,912)   ", ",rlim(j)
        WRITE (*,913)   " ):"
        DO i = 1, ny
          WRITE (*,911) "  integral_of_y",i
          WRITE (*,912) "=",r(j+(i-1)*nlim)
        END DO
        WRITE (*,910)   ""
      END DO


      !***** Delete Data Fitting task *****
      errcode = dfdeletetask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,909) "","Error: Computed Akima cubic spline             &
     &   coefficients"
        WRITE (*,905) " are incorrect"
        STOP 1
      ELSE
        WRITE (*,909) "","Computed Akima cubic spline coefficients"
        WRITE (*,905) " are correct"
      END IF
      STOP 0

901    FORMAT (A,I0)
902    FORMAT (/A)
903    FORMAT (A,I1,A,SP,F11.6,A,SP,F11.6)
904    FORMAT (/A/A)
905    FORMAT (99A)
915    FORMAT (A,$)
906    FORMAT (A,I1,A,SP,F13.6,A,SP,F13.6,A,SP,F13.6,A,SP,F13.6,A,SP,    &
     &         F11.6,A,SP,F11.6)
907    FORMAT (A,SP,F11.6,A,SP,F11.6)
908    FORMAT (/A,F4.1,A,F4.1,A,F4.1)
909    FORMAT (//A)

910    FORMAT (A)         ! Text, new line.
911    FORMAT (A,I0,$)    ! Text, integer.
912    FORMAT (A,F11.6,$) ! Text, float.
913    FORMAT (A,$)       ! Text.

      END PROGRAM


      FUNCTION right_akima_integr(n,lcell,llim,rcell,rlim,r,xN,libpars)  &
     & RESULT (output_4)
      USE MKL_DF_TEST_TYPES

        INTEGER    :: output_4
        INTEGER(8) :: n
        INTEGER(8) :: lcell(*),rcell(*)
        INTEGER    :: xN(DF_TEST_PARAMS_4BYTE_SIZE)
        REAL(8)    :: llim(*),rlim(*),r(*)

        INTEGER(8) :: i
        REAL(8)    :: x

        x = transfer(xN, x)
        DO i = 0,n-1
          r(i+1) = x * ( rlim(i+1) - x )
        END DO
        output_4 = 0
        RETURN

      END FUNCTION


      FUNCTION left_akima_integr(n,lcell,llim,rcell,rlim,r,x0,libpars)   &
     & RESULT (output_4)
      USE MKL_DF_TEST_TYPES

        INTEGER    :: output_4
        INTEGER(8) :: n
        INTEGER(8) :: lcell(*),rcell(*)
        INTEGER    :: x0(DF_TEST_PARAMS_4BYTE_SIZE)
        REAL(8)    :: llim(*),rlim(*),r(*)

        INTEGER(8) :: i
        REAL(8)    :: x

        x = transfer(x0, x)
        DO i = 0,n-1
          r(i+1) = x * ( x - llim(i+1) )
        END DO
        output_4 = 0
        RETURN

      END FUNCTION



      FUNCTION integr_search_cb(n,site,cell,flag,                        &
     &  user_params,libparams)       RESULT (res)

        USE MKL_DF_TYPE
        USE MKL_DF_TEST_TYPES
        USE MKL_DF_TEST_VARS

        ! FUNCTION ARGUMENTS DECLARATION
        INTEGER(8),INTENT(IN)                   :: n
        REAL(8),INTENT(IN)                      :: site(*)
        INTEGER(8),INTENT(OUT)                  :: cell(*)
        INTEGER(4),INTENT(OUT)                  :: flag(*)
        TYPE(MY_PARAMS),INTENT(IN)              :: user_params
        TYPE(DF_SEARCH_CALLBACK_LIBRARY_PARAMS),INTENT(IN)::libparams

        ! LOCAL VARIABLES DECLARATION
        INTEGER(8)                              :: k
        INTEGER(8)                              :: nlim1
        INTEGER                                 :: res
        INTEGER                                 :: i
        INTEGER                                 :: j
        INTEGER(4)                              :: lim_type
        INTEGER,PARAMETER                       :: NNLIM1 = 6

        lim_type = libparams%limit_type_flag

        nlim1 = n

        IF (user_params%v1 .NE. 1.0  .OR.                                &
     &      user_params%v2 .NE. 2.0  .OR.                                &
     &      user_params%v3 .NE. 3.0) THEN
          WRITE (*,700) "Error: one of user_params is incorrect"
          res = DF_STATUS_EXACT_RESULT
          RETURN
        END IF

        IF (lim_type == DF_INTEGR_SEARCH_CB_LLIM_FLAG) THEN
          ! Here if search cb called for left integration limits
          DO i = 1, nlim1, 1
            cell(i) = llim_cell8(NNLIM)
            DO j = 1, NNLIM1, 1
              IF (site(i) .EQ. llim(j)) THEN
                cell(i) = llim_cell8(j)
                EXIT
              END IF
            END DO
            flag(i) = 1
          END DO
        ELSE
          ! Here if search cb called for right integration limits
          DO i = 1, nlim1, 1
            cell(i) = rlim_cell8(NNLIM)
            DO j = 1, NNLIM1, 1
              IF (site(i) .EQ. rlim(j)) THEN
                cell(i) = rlim_cell8(j)
                EXIT
              END IF
            END DO
            flag(i) = 1
          END DO
        END IF

        res = DF_STATUS_EXACT_RESULT

        RETURN

700    FORMAT (A)

      END FUNCTION
