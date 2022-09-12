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
!    Construction and integration of linear spline Example Program Text
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
      INTEGER,PARAMETER :: N          = 10
      ! number of functions
      INTEGER,PARAMETER :: NNY        = 2
      ! number of pairs of integration limits
      INTEGER,PARAMETER :: NNLIM      = 3
      ! total number of spline coefficients
      INTEGER,PARAMETER :: NNSCOEFF   = (NNY*(N-1)*DF_PP_LINEAR)
      ! left  limit of interpolation interval
      REAL(4),PARAMETER :: LLIM_X     = -2.0
      ! right limit of interpolation interval
      REAL(4),PARAMETER :: RLIM_X     =  2.0
      REAL(4),PARAMETER :: LEFT_LLIM  = -1.5
      REAL(4),PARAMETER :: RIGHT_LLIM = -0.5
      REAL(4),PARAMETER :: LEFT_RLIM  =  0.5
      REAL(4),PARAMETER :: RIGHT_RLIM =  1.5
      REAL(4),PARAMETER :: FFREQ      =  1.3

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
      ! array of break points
      REAL(4) :: x(N)
      ! function values
      REAL(4),TARGET  :: y(NNY*N)
      REAL(4),POINTER :: y_ptr(:)
      ! array of spline coefficients
      REAL(4),TARGET :: scoeff(NNSCOEFF)
      ! left and right limits of the intervals are provided
      ! if integration limits are uniform
      REAL(4) :: llim(2),rlim(2)
      ! integration results
      REAL(4) :: r(NNY*NNLIM)
      ! left integration limits in full format
      REAL(4) :: llim_full(NNLIM)
      ! right integration limits in full format
      REAL(4) :: rlim_full(NNLIM)
      REAL(4),POINTER :: scoeff_ptr(:)

      REAL(4) :: left,right
      REAL(4) :: freq
      INTEGER :: i,j,errnums
      INTEGER(4) :: errcode

      errcode = DF_STATUS_OK
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****
      sorder = DF_PP_LINEAR
      stype = DF_PP_DEFAULT

      !***** Parameters describing interpolation interval *****
      left = LLIM_X
      right = RLIM_X
      nx = N
      xhint = DF_NON_UNIFORM_PARTITION

      !***** Parameters describing function *****
      ny = NNY
      yhint = DF_NO_HINT

      !***** Parameters describing spline coefficients storage *****
      scoeffhint = DF_NO_HINT

      !***** Parameters describing boundary conditions type *****
      !***** No boundary conditions are provided for linear spline *****
      bc_type = DF_NO_BC

      !***** Parameters describing internal conditions type *****
      !***** No internal conditions are provided for linear spline *****
      ic_type = DF_NO_IC

      !***** Parameters decsribing integration limits *****
      nlim = NNLIM
      llimhint = DF_UNIFORM_PARTITION
      rlimhint = DF_UNIFORM_PARTITION

      !***** Left and right points of interval are provided in case of uniform partition *****
      llim(1) = LEFT_LLIM
      llim(2) = RIGHT_LLIM
      rlim(1) = LEFT_RLIM
      rlim(2) = RIGHT_RLIM

      !***** Parameter dascribing integration results storage format *****
      rhint = DF_NO_HINT

      !***** Generate array of uniformly distributed break points *****
      errcode = sUniformRandSortedData( x, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Generate functions y = sin(2 * Pi * freq * x) *****

      DO i = 0, ny-1
        freq = (i + 1) * FFREQ
        y_ptr => y(i*nx+1:)
        errcode = sSinDataNotUniformGrid( y_ptr, x, freq, nx )
        CALL CheckDfError(errcode)
      END DO

      !***** Create Data Fitting task *****
      errcode = dfsnewtask1d( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for linear spline construction *****
      errcode = dfseditppspline1d( task, sorder, stype, bc_type,         &
     &     ic_type=ic_type, scoeff=scoeff, scoeffhint=scoeffhint )
      CALL CheckDfError(errcode)

      !***** Construct linear spline using STD method *****
      errcode = dfsconstruct1d( task, DF_PP_SPLINE, DF_METHOD_STD )
      CALL CheckDfError(errcode)

      !***** Compute integrals *****
      errcode = dfsintegrate1d( task, DF_METHOD_PP, nlim, llim,          &
     &    llimhint, rlim, rlimhint, r=r, rhint=rhint )
      CALL CheckDfError(errcode)

      !***** Delete Data Fitting task *****
      errcode = dfdeletetask( task )
      CALL CheckDfError(errcode)

      !***** Check computed coefficients *****
      errcode = sCheckLinSplineBreakPoints( nx, x, ny, y, scoeff )
      IF (errcode < 0) errnums = errnums+1

      !***** Print results *****
      WRITE (*,901) "Number of break points : ",nx

      !***** Print given function *****
      WRITE (*,902) "    X             Y1            Y2"
      DO j = 1, nx
        WRITE (*,903) " ",x(j),"   ",y(j),"   ",y(j+nx)
      END DO

      !***** Print computed spline coefficients *****
      WRITE (*,904) "Coefficients are calculated for a polynomial of the &
     & form:",""
      WRITE (*,905) "Pi(x) = Ai + Bi*(x - x(i))"
      WRITE (*,905) "    where x(i) <= x < x(i+1)"
      scoeff_ptr => scoeff
      DO i = 0, ny-1
        WRITE (*,906) "Spline coefficients for Y",i+1 ," :"
        WRITE (*,905) " i    Ai            Bi"
        DO j = 0, nx-2
          WRITE (*,907) " ",j," ",scoeff_ptr(sorder*j + 1),"   ",        &
     &                  scoeff_ptr(sorder*j + 2)
        END DO
        scoeff_ptr => scoeff_ptr(( sorder*(nx-1))+1:)
      END DO

      !***** Print integration results *****
      errcode = sUniformData(llim_full, llim(1), llim(2), nlim)
      CALL CheckDfError(errcode)
      errcode = sUniformData(rlim_full, rlim(1), rlim(2), nlim)
      CALL CheckDfError(errcode)
      DO i = 0, ny-1
        WRITE (*,906) "Integration results for Y",i+1," :"
        WRITE (*,905) "Integration interval      Result"
        DO j = 1, nlim
          WRITE (*,908) " ( ",llim_full(j),", ",rlim_full(j),            &
     &                  " )        ",r(i*nlim + j)
        END DO
      END DO

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,909) "",                                                &
     &    "Error: Computed linear spline coefficients or integrals"
        WRITE (*,905) " are incorrect"
        STOP 1
      ELSE
        WRITE (*,909) "",                                                &
     &  "Computed linear spline coefficients and integrals"
        WRITE (*,905) " are correct"
      END IF
      STOP 0

901   FORMAT (A,I0)
902   FORMAT (/A)
903   FORMAT (A,SP,F11.6,A,SP,F11.6,A,SP,F11.6)
904   FORMAT (/A/A)
905   FORMAT (99A)
906   FORMAT (/A,I0,A)
907   FORMAT (A,I1,A,SP,F11.6,A,SP,F11.6)
908   FORMAT (A,SP,F4.1,A,SP,F4.1,A,SP,F11.6)
909   FORMAT (//A)

      END PROGRAM
