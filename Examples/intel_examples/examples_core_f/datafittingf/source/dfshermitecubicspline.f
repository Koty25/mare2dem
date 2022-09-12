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
!    Construction of Hermite cubic spline and interpolation
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
      INTEGER,PARAMETER :: N        = 7
      ! number of functions
      INTEGER,PARAMETER :: NNY      = 1
      ! number of boundary conditions
      INTEGER,PARAMETER :: NNBC     = 2
      ! number of internal conditions
      INTEGER,PARAMETER :: NNIC     = N-2
      ! number of interpolation sites
      INTEGER,PARAMETER :: NNSITE   = 15
      ! size of array describing derivative orders to compute
      INTEGER,PARAMETER :: NNDORDER = 1
      ! total number of spline coefficients
      INTEGER,PARAMETER :: NNSCOEFF = (NNY*(N-1)*DF_PP_CUBIC)
      ! left  limit of interpolation interval
      REAL(4),PARAMETER :: LLIM_X   = 0.0
      ! right limit of interpolation interval
      REAL(4),PARAMETER :: RLIM_X   = 3.0
      REAL(4),PARAMETER :: FFREQ    = 0.3

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
      ! number of spline coefficients
      INTEGER :: nscoeff
      ! spline coefficients storage format
      INTEGER :: scoeffhint
      ! boundary conditions type
      INTEGER :: bc_type
      ! number of boundary conditions
      INTEGER :: nbc
      ! internal conditions type
      INTEGER :: ic_type
      ! number of internal conditions
      INTEGER :: nic
      ! number of interpolation sites
      INTEGER :: nsite
      ! additional info about interpolation sites
      INTEGER :: sitehint
      !  size of array describing derivative orders
      INTEGER :: ndorder
      ! interpolation results storage format
      INTEGER :: rhint
      ! left limit of the interpolation interval
      REAL(4) :: left
      ! right limit of the interpolation interval
      REAL(4) :: right
      ! limits of the interpolation interval
      REAL(4) :: x(2)
      ! partition in full format
      REAL(4) :: xx(N)
      ! function values
      REAL(4) :: y(NNY*N)
      ! array of internal conditions
      REAL(4) :: ic(NNIC)
      ! array of boundary conditions
      REAL(4) :: bc(NNBC)
      ! array of spline coefficients
      REAL(4) :: scoeff(NNSCOEFF)
      ! array of interpolation sites
      REAL(4) :: site(NNSITE)
      ! spline evaluation results
      REAL(4),TARGET :: r(NNDORDER*NNSITE)
      ! array of derivatives orders
      INTEGER :: dorder(NNDORDER)
      ! type of the computations
      INTEGER :: type

      REAL(4) :: llim(NNSITE)
      REAL(4) :: rlim(NNSITE)
      REAL(4) :: left_val(N-1),right_val(N-1)
      REAL(4) :: left_der(N-1),right_der(N-1)
      REAL(4),TARGET :: ref_r(NNDORDER*NNSITE)
      REAL(4),POINTER :: r_ptr(:), ref_r_ptr(:)

      REAL(4) :: freq,delta,dif
      INTEGER :: i,j,errnums
      INTEGER(4) :: errcode

      errcode = 0
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****
      sorder = DF_PP_CUBIC
      stype  = DF_PP_HERMITE

      !***** Parameters describing interpolation interval *****
      left = LLIM_X
      right = RLIM_X
      nx = N
      xhint = DF_UNIFORM_PARTITION

      !***** Limits of interpolation interval are provided in case of uniform partition *****
      x(1) = left
      x(2) = right

      !***** Parameters describing function *****
      ny = NNY
      yhint = DF_NO_HINT

      !***** Parameters describing spline coefficients storage *****
      nscoeff = NNSCOEFF

      !Row-major storage will be user in case no additional info about storage is provided
      scoeffhint = DF_NO_HINT

      !***** Parameters describing boundary conditions type *****
      bc_type = IOR(DF_BC_1ST_LEFT_DER,DF_BC_1ST_RIGHT_DER)

      !***** Parameters describing internal conditions *****
      ic_type = DF_IC_1ST_DER
      nic = NNIC
      !***** Parameters describing interpolation sites *****
      nsite = NNSITE
      sitehint = DF_NON_UNIFORM_PARTITION

      !**** Parameter describing interpolation results storage *****
      rhint = DF_MATRIX_STORAGE_COLS

      !**** Parameter describing array for derivative orders *****
      ndorder = NNDORDER
      dorder = (/1/)

      !***** Generate partition in full format *****
      errcode = sUniformData( xx, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Generate function y = sin(2 * Pi * freq * x) *****
      freq = FFREQ
      errcode = sSinDataNotUniformGrid( y, xx, freq, nx )
      CALL CheckDfError(errcode)

      !***** Generate array of 1st derivatives needed for Hermite spline construction *****
      errcode = sSinDerDataNotUniformGrid( ic, xx(2), freq, nic )
      CALL CheckDfError(errcode)

      !***** Generate boundary conditions *****
      errcode = sSinDerDataNotUniformGrid( bc(1), xx(1),  freq, 1 )
      CALL CheckDfError(errcode)
      errcode = sSinDerDataNotUniformGrid( bc(2), xx(nx), freq, 1 )
      CALL CheckDfError(errcode)

      !***** Generate interpolation sites *****
      errcode = sUniformRandData( site, left, right, nsite )
      CALL CheckDfError(errcode)

      !***** Create Data Fitting task *****
      errcode = dfsnewtask1d( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for Hermite cubic spline construction *****
      errcode = dfseditppspline1d( task, sorder, stype, bc_type, bc,     &
     &      ic_type, ic, scoeff, scoeffhint )
      CALL CheckDfError(errcode)

      !***** Construct Hermite cubic spline using STD method *****
      errcode = dfsconstruct1d( task, DF_PP_SPLINE, DF_METHOD_STD )
      CALL CheckDfError(errcode)

      !***** Interpolate using PP method *****
      type = DF_INTERP
      errcode = dfsinterpolate1d( task, type, DF_METHOD_PP, nsite, site, &
     &      sitehint, ndorder, dorder, r=r, rhint=rhint)
      CALL CheckDfError(errcode)

      !***** Check computed coefficients *****

      !***** Check spline values in break points *****
      errcode = sCheckCubBreakPoints( nx, xx, ny, y, scoeff,             &
     &      left_val, right_val )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Check that spline 1st derivatives are equal for left *****
      !***** and right piece of the spline for each break point *****
      errcode = sCheckCub1stDerConsistency( nx, xx, ny, scoeff,          &
     &      left_der, right_der )
      IF ( errcode < 0 ) errnums = errnums+1

      !***** Check internal conditions *****
      DO i = 1, nic
        IF ( ABS(ic(i) - left_der(i)) > EPSILON_SINGLE )                 &
     &      errnums = errnums + 1
      END DO

      !***** Check boundary conditions *****
      errcode = sCheckCubBC( nx, xx, ny, scoeff, bc_type, bc )
      IF ( errcode < 0 ) errnums = errnums + 1

      !***** Check results of interpolation *****
      errcode = sCheckCubInterpRes( nx, xx, ny, scoeff, nsite, site,     &
     &      ndorder, dorder, r, ref_r )
      IF ( errcode < 0 ) errnums = errnums + 1

      !***** Print results *****
      WRITE (*,901) "Number of break points : ",nx

      !***** Print given function *****
      WRITE (*,902) " i  x(i)           y(i)          y'(i)"
      WRITE (*,903) " ",1," ",xx(1),"   ",y(1),"   ",bc(1)
      DO j = 2, nx-1
        WRITE (*,903) " ",j," ",xx(j),"   ",y(j),"   ",ic(j-1)
      END DO
      WRITE (*,903) " ",nx," ",xx(nx),"   ",y(nx),"   ", bc(2)

      !***** Print computed spline coefficients *****
      WRITE (*,904) "Coefficients are calculated for a polynomial of     &
     &the form:",""
      WRITE (*,905) "Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 +      &
     &Di*(x - x(i))^3"
      WRITE (*,905) "    where x(i) <= x < x(i+1)"
      WRITE (*,902) "Spline coefficients for Y:"
      WRITE (*,910) " i    Ai            Bi            Ci"
      WRITE (*,910) "            Di"
      WRITE (*,910) "            P(x(i))       P(x(i+1)) "
      WRITE (*,905) "    P'(x(i))      P'(x(i+1))"
      DO j = 1, nx-1
        WRITE (*,906,ADVANCE='NO') " ",j," ",scoeff(sorder*(j-1) + 1),   &
     & "   ",scoeff(sorder*(j-1) + 2),"   ",scoeff(sorder*(j-1) + 3),    &
     & "   ",scoeff(sorder*(j-1) + 4),"   ",right_val(j),"   ",          &
     & left_val(j)
        WRITE (*,907) "   ",right_der(j),"   ",left_der(j)
      END DO
      WRITE (*,902) "  i       Sites          Spline value"
      WRITE (*,905) "                         Computed       Expected"
      DO j = 1,nsite
        WRITE (*,908) " ",j,"    ",site(j),"    ",r(j),"    ",           &
     &  ref_r(j)
      END DO

      !***** Delete Data Fitting task *****
      errcode = dfdeletetask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,909) "","Error: Computed Hermite cubic spline           &
     &  coefficients, spline values are incorrect"
        STOP 1
      ELSE
        WRITE (*,909) "","Computed Hermite cubic spline coefficients     &
     &  spline values are correct"
      END IF
      STOP 0

901    FORMAT (A,I0)
902    FORMAT (/A)
912    FORMAT (/A,$)
903    FORMAT (A,I1,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6)
904    FORMAT (/A/A)
905    FORMAT (99A)
906    FORMAT (A,I1,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,    &
     &         F11.6,A,SP,F11.6)
907    FORMAT (A,SP,F11.6,A,SP,F11.6)
908    FORMAT (A,I2,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,       &
     &         SP,F11.6)
909    FORMAT (/A)
910    FORMAT(A,$)

      END PROGRAM
