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
!    Construction of natural cubic spline and spline-based integration
!  Example Program Text
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
      INTEGER,PARAMETER :: N          = 6
      ! number of functions
      INTEGER,PARAMETER :: NNY        = 1
      ! number of pairs of integration limits
      INTEGER,PARAMETER :: NNLIM      = 4
      ! total number of spline coefficients
      INTEGER,PARAMETER :: NNSCOEFF   = (NNY*(N-1)*DF_PP_CUBIC)
      ! left  limit of interpolation interval
      REAL(8),PARAMETER :: LLIM_X     =  0.0d0
      ! right limit of interpolation interval
      REAL(8),PARAMETER :: RLIM_X     =  7.0d0
      REAL(8),PARAMETER :: LEFT_LLIM  =  1.0d0
      REAL(8),PARAMETER :: RIGHT_LLIM =  3.0d0
      REAL(8),PARAMETER :: LEFT_RLIM  =  4.0d0
      REAL(8),PARAMETER :: RIGHT_RLIM =  5.5d0
      REAL(8),PARAMETER :: FFREQ      =  0.15d0

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
      ! number of pairs of integration limits
      INTEGER :: nlim
      ! additional info about the structure of left integration limits
      INTEGER :: llimhint
      ! additional info about the structure of right integration limits
      INTEGER :: rlimhint
      ! integration results storage format
      INTEGER :: rhint
      ! limits of the interpolation interval
      REAL(8) :: x(2)
      ! break points in full format
      REAL(8) :: xx(N)
      ! function values
      REAL(8) :: y(NNY*N)
      ! array of spline coefficients
      REAL(8) :: scoeff(NNSCOEFF)
      ! left integration limits
      REAL(8) :: llim(NNLIM)
      ! right integration limits
      REAL(8) :: rlim(NNLIM)
      ! integration results
      REAL(8) :: r(NNY*NNLIM)
      ! reference integration results
      REAL(8) :: ref_r(NNY*NNLIM)

      REAL(8) :: left,right
      REAL(8) :: l_llim, r_llim, l_rlim, r_rlim
      REAL(8) :: freq
      REAL(8) :: left_val(N-1), right_val(N-1)
      REAL(8) :: left_der1(N-1), right_der1(N-1)
      REAL(8) :: left_der2(N-1), right_der2(N-1)

      INTEGER :: i,j,errnums
      INTEGER(4) :: errcode

      errcode = DF_STATUS_OK
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****
      sorder = DF_PP_CUBIC
      stype  = DF_PP_NATURAL

      !***** Parameters describing interpolation interval *****
      nx = N
      xhint = DF_UNIFORM_PARTITION
      ! Limits of interpolation interval are provided in case of uniform grid
      left = LLIM_X
      right = RLIM_X
      x(1) = left
      x(2) = right

      !***** Parameters describing function *****
      ny = NNY
      yhint = DF_NO_HINT

      !***** Parameters describing spline coefficients storage *****
      scoeffhint = DF_NO_HINT

      !***** Parameters describing boundary conditions type *****
      bc_type = DF_BC_FREE_END

      !***** Parameters describing internal conditions type *****
      !***** No internal conditions are provided for natural cubic spline *****
      ic_type = DF_NO_IC

      !***** Parameters decsribing integration limits *****
      nlim = NNLIM
      llimhint = DF_NON_UNIFORM_PARTITION
      rlimhint = DF_NON_UNIFORM_PARTITION

      !***** Parameter dascribing integration results storage format *****
      rhint = DF_NO_HINT

      freq = FFREQ
      errcode = dSinDataUniformGrid( y, freq, left, freq, nx )
      CALL CheckDfError(errcode)

      !***** Generate limits of integration intervals *****
      l_llim = LEFT_LLIM
      r_llim = RIGHT_LLIM
      errcode = dUniformRandSortedData( llim, l_llim, r_llim, nlim )
      CALL CheckDfError(errcode)

      l_rlim = LEFT_RLIM
      r_rlim = RIGHT_RLIM
      errcode = dUniformRandSortedData( rlim, l_rlim, r_rlim, nlim )
      CALL CheckDfError(errcode)

      !***** Create Data Fitting task *****
      errcode = dfdnewtask1d( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for natural cubic spline construction *****
      errcode = dfdeditppspline1d( task, sorder, stype, bc_type,         &
     &     ic_type=ic_type, scoeff=scoeff, scoeffhint=scoeffhint )
      CALL CheckDfError(errcode)

      !***** Construct natural cubic spline using STD method *****
      errcode = dfdconstruct1d( task, DF_PP_SPLINE, DF_METHOD_STD )
      CALL CheckDfError(errcode)

      !***** Compute integrals *****
      errcode = dfdintegrate1d( task, DF_METHOD_PP, nlim, llim,          &
     &    llimhint, rlim, rlimhint, r=r, rhint=rhint )
      CALL CheckDfError(errcode)

      !***** Check computed coefficients *****
      errcode = dUniformData( xx, left, right, nx );
      CALL CheckDfError(errcode)

      !***** Check spline values in break points *****
      errcode = dCheckCubBreakPoints( nx, xx, ny, y, scoeff,             &
     &      left_val, right_val )
      CALL CheckDfError(errcode)

      !***** Check that spline 1st derivatives are equal for left
      !      and right piece of the spline for each break point *****
      errcode = dCheckCub1stDerConsistency( nx, xx, ny, scoeff,          &
     &      left_der1, right_der1 )
      CALL CheckDfError(errcode)

      !***** Check that spline 2nd derivatives are equal for left
      !     and right piece of the spline for each break point *****/
      errcode = dCheckCub2ndDerConsistency( nx, xx, ny, scoeff,          &
     &      left_der2, right_der2 )
      CALL CheckDfError(errcode)

      !***** Check boundary conditions *****
      errcode = dCheckCubBC( nx, xx, ny, scoeff, bc_type )
      CALL CheckDfError(errcode)

      !***** Check results of integration *****
      errcode = dCheckCubIntegrRes( nx, xx, ny, scoeff,                  &
     & nlim, llim, rlim, r, ref_r )
      IF (errcode < 0) errnums = errnums+1

      !***** Print results *****
      WRITE (*,901) "Number of break points : ",nx

      !***** Print given function *****
      WRITE (*,902) " i    x(i)          y(i)"
      DO j = 1, nx
        WRITE (*,907) " ",j," ",xx(j),"   ",y(j)
      END DO

      !***** Print computed spline coefficients *****
      WRITE (*,904) "Coefficients are calculated for a polynomial of the &
     & form:",""
      WRITE (*,902) "Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 + Di*(x&
     & - x(i))^3"
      WRITE (*,904) "    where x(i) <= x < x(i+1)",""
      WRITE (*,902) "Spline coefficients for Y:"
      WRITE (*,902) " i    Ai          Bi          Ci          Di        &
     &   P(x(i))     P(x(i+1))   P'(x(i))    P'(x(i+1))  P''(x(i))    P''&
     &(x(i+1))"

      DO j = 1, nx-1
        WRITE (*,910) " ", j, " ", scoeff(sorder*(j-1) + 1), " ",        &
     &      scoeff(sorder*(j-1) + 2), " ", scoeff(sorder*(j-1) + 3)," ", &
     &      scoeff(sorder*(j-1) + 4)
        WRITE (*,911) " ",right_val(j)," ",left_val(j)," ",              &
     &      right_der1(j)," ",left_der1(j)
        WRITE (*,912) " ",right_der2(j)," ",left_der2(j)
      END DO

      WRITE (*,904) "","Integration results for Y:"
      WRITE (*,902) " Integration interval               Result"

      DO j = 1, nlim
        WRITE (*,913) " ( ",llim(j)," , ",rlim(j)," )    ",r(j)
      END DO

      !***** Delete Data Fitting task *****
      errcode = dfdeletetask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,904) "","Error: Computed natural cubic spline           &
     &    coefficients or integrals are incorrect"
        STOP 1
      ELSE
        WRITE (*,904) "","Computed natural cubic spline coefficients and &
     &  integrals are correct"
      END IF
      STOP 0

901   FORMAT (A,I0)
902   FORMAT (/A)
904   FORMAT (/A/A)
907   FORMAT (A,I1,A,SP,F11.6,A,SP,F11.6)
910   FORMAT (A,I1,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,$)
911   FORMAT (A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,$)
912   FORMAT (A,SP,F11.6,A,SP,F11.6)
913   FORMAT (A,SP,F11.6,A,SP,F11.6,A,SP,F11.6)

      END PROGRAM
