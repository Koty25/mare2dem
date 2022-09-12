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
!    Construction of natural cubic spline, interpolation, and computation
!    of the 2nd derivative Example Program Text
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
      INTEGER,PARAMETER :: N           = 9
      ! number of functions
      INTEGER,PARAMETER :: NNY         = 1
      ! number of sites for spline-based interpolation in one block
      INTEGER,PARAMETER :: NSITE_BLOCK = 10
      ! number of blocks
      INTEGER,PARAMETER :: NNBLOCKS    = 4
      ! size of array describing derivative orders to compute
      INTEGER,PARAMETER :: NNDERORDER  = 3
      ! number of derivatives to compute
      INTEGER,PARAMETER :: NNDER       = 2
      ! total number of spline coefficients
      INTEGER,PARAMETER :: NNSCOEFF    = NNY*(N-1)*DF_PP_CUBIC
      ! total number of interpolation sites
      INTEGER,PARAMETER :: NNSITE      = NSITE_BLOCK*NNBLOCKS

      ! left  limit of interpolation interval
      REAL(4),PARAMETER :: LEFT_LIMIT  = -1.0
      ! right limit of interpolation interval
      REAL(4),PARAMETER :: RIGHT_LIMIT =  2.0
      REAL(4),PARAMETER :: FFREQ = 1.7

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
      ! number of boundary conditions
      INTEGER :: nbc
      ! number of interpolaton sites in one block
      INTEGER :: nsite_bl
      ! number of blocks
      INTEGER :: nblocks
      ! total number of interpolation sites
      INTEGER :: nsite
      ! additional info about interpolation sites
      INTEGER :: sitehint
      ! size of array describing derivative orders
      INTEGER :: ndorder
      ! array of derivative orders
      INTEGER :: derorder(3)
      ! interpolation results storage format
      INTEGER :: rhint
      ! left limit of the interpolation interval
      REAL(4) :: left
      ! right limit of the interpolation interval
      REAL(4) :: right
      ! limits of the interpolation interval
      REAL(4) :: x(2)
      ! function values
      REAL(4) :: y(NNY*N)
      ! array of spline coefficients
      REAL(4) :: scoeff(NNSCOEFF)
      ! array of interpolation sites
      REAL(4),TARGET :: site(NNSITE)
      ! spline evaluation results
      REAL(4),TARGET :: r(NNDER*NNSITE)
      ! type of calculations
      INTEGER :: type
      ! method that is used to perform calculations
      INTEGER :: method

      REAL(4),POINTER :: site_ptr(:)
      REAL(4),POINTER :: r_ptr(:)

      INTEGER :: i,j,nder, yhint
      INTEGER :: errnum
      INTEGER(4) :: errcode

      REAL(4) :: freq

      REAL(4) :: xx(N)

      REAL(4) :: ref_r(NNDER*NNSITE)

      errcode = 0
      errnums = 0

      ! **** Initializing parameters for Data Fitting task ****

      sorder = DF_PP_CUBIC
      stype  = DF_PP_NATURAL

      !***** Parameters describing interpolation interval *****
      nx    = N
      xhint = DF_UNIFORM_PARTITION
      ! Limits of interpolation interval are provided in case
      ! of uniform partition
      left  = LEFT_LIMIT
      right = RIGHT_LIMIT
      x(1) = left
      x(2) = right

      !***** Parameters describing function *****
      ny = NNY
      yhint = DF_NO_HINT

      !***** Parameters describing spline coefficients storage *****
      nscoeff = NNSCOEFF
      scoeffhint = DF_NO_HINT

      !***** Parameters describing boundary conditions type *****
      !  No additional parameters are provided for Free-End boundary conditions
      bc_type = DF_BC_FREE_END

      !***** Parameters describing interpolation sites *****
      nsite_bl   = NSITE_BLOCK
      nsite      = NNSITE
      sitehint   = DF_NON_UNIFORM_PARTITION

      !***** Parameter describing interpolation results storage *****
      rhint = DF_MATRIX_STORAGE_ROWS

      !***** Parameter describing array for derivative orders *****
      ndorder = NNDERORDER
      ! spline values and 2nd derivatives will be computed
      derorder = (/1,0,1/)

      !***** Number of derivatives to compute *****
      nder = NNDER

      !***** Generate function y = sin(2 * Pi * freq * x) *****
      freq = FFREQ
      errcode = sSinDataUniformGrid( y, freq, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Generate interpolation sites *****
      errcode = sUniformRandData( site, left, right, nsite )
      CALL CheckDfError(errcode)

      !***** Create Data Fitting task *****
      errcode = dfsNewTask1D( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for natural cubic spline construction *****
      errcode = dfsEditPPSpline1D( task, sorder, stype, bc_type,         &
     & scoeff=scoeff, scoeffhint=scoeffhint )
      CALL CheckDfError(errcode)

      !***** Construct natural cubic spline using STD method *****
      type = DF_PP_SPLINE
      method = DF_METHOD_STD

      errcode = dfsConstruct1D( task, type, method )
      CALL CheckDfError(errcode)

      !***** Interpolate, and compute 2nd order derivatives
      !      using STD method *****
      nblocks = NNBLOCKS

      type = DF_INTERP
      method = DF_METHOD_PP

      DO i = 1, nblocks
        site_ptr => site((i-1)*nsite_bl+1:)
        r_ptr    => r((i-1)*nder*nsite_bl+1:)
        errcode = dfsInterpolate1D( task, type, method,                  &
     &      nsite_bl, site_ptr, sitehint, ndorder, derorder, r=r_ptr,    &
     &      rhint=rhint )
        CALL CheckDfError(errcode)
      END DO

      !***** Check results of interpolation *****
      errcode = sUniformData( xx, left, right, nx )
      CALL CheckDfError(errcode)

      errcode = sCheckCubInterpRes( nx, xx, ny, scoeff,                  &
     & nsite, site, ndorder, derorder, r, ref_r )
      IF (errcode < 0) errnums = errnums+1

      !***** Print results *****

      WRITE (*,901) "Number of break points : ", nx

      !***** Print function ****
      WRITE (*,902) " i     x(i)          y(i)"
      DO j = 1, nx
        WRITE (*,903) " ",j,"  ",xx(j),"   ",y(j)
      END DO

      !***** Print spline coefficients *****
      WRITE (*,902) "Spline coefficients for Y:"
      WRITE (*,904) "  i   Ai           Bi           Ci           Di"
      DO j = 1, nx-1
        WRITE (*,905) " ",scoeff(sorder*(j-1)+1),"   ",                  &
     &    scoeff(sorder*(j-1) + 2),"   ",scoeff(sorder*(j-1) + 3),"   ", &
     &    scoeff(sorder*(j-1) + 4)
      END DO

      !***** Print interpolation results *****
      WRITE (*,902) "Results of interpolation :"
      WRITE (*,914) "                    Function value"
      WRITE (*,904) "             Second derivative"
      WRITE (*,904) "    Sites         Obtained     Expected    Obtained &
     &     Expected"
      DO i = 0,nsite-1
        WRITE (*,906,ADVANCE='NO') " ",site(i+1)
        DO j = 0,nder-1
           WRITE (*,907,ADVANCE='NO') "   ",r(i*nder + j+1),"  ",        &
     &      ref_r(i*nder + j+1)
        END DO
        WRITE (*,904) ""
      END DO

      !***** Delete Data Fitting task *****
      errcode = dfDeleteTask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,908) "","Error: Computed interpolation results"
        WRITE (*,904) " are incorrect"
        STOP 1
      ELSE
        WRITE (*,908) "","Computed interpolation results"
        WRITE (*,904) " are correct"
      END IF
      STOP 0


  901 FORMAT (A,I0)
  902 FORMAT (/A)
  903 FORMAT (A,I0,A,SP,F11.6,A,SP,F11.6)
  904 FORMAT (99A)
  914 FORMAT (A,$)
  905 FORMAT (A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6)
  906 FORMAT (A,SP,F11.6)
  907 FORMAT (A,SP,F11.6,A,SP,F11.6)
  908 FORMAT (//A)

      END PROGRAM
