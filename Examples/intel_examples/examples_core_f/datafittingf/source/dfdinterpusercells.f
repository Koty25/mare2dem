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
!    Construction and interpolation of linear spline using user provided
!  cell indices Example Program Text
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
      INTEGER,PARAMETER :: N           = 7
      ! number of functions
      INTEGER,PARAMETER :: NNY         = 1
      ! number of interpolation sites
      INTEGER,PARAMETER :: NNSITE      = 10
      ! number of spline coefficients
      INTEGER,PARAMETER :: NSCOEFF     = DF_PP_LINEAR*(N-1)
      ! left  limit of interpolation interval
      REAL(8),PARAMETER :: LLIM_X      = -10.0d0
      ! right limit of interpolation interval
      REAL(8),PARAMETER :: RLIM_X      = 10.0d0

      REAL(8),PARAMETER :: FFREQ       = 1.75d0

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
      ! number of derivatives to calculate
      INTEGER :: ndorder
      ! derivatives to calculate
      INTEGER :: dorder(1)
      ! interpolation results storage format
      INTEGER :: rhint
      ! type of computations
      INTEGER :: comp_type
      ! computation method
      INTEGER :: method

      ! limits of interpolation interval
      REAL(8) :: left, right
      ! array of break points
      REAL(8) :: x(N)
      ! function values
      REAL(8) :: y(N)
      ! array of spline coefficients
      REAL(8) :: scoeff(NSCOEFF)
      ! interpolation sites
      REAL(8) :: site(NNSITE)
      ! interpolation results
      REAL(8) :: r(NNSITE)
      ! reference interpolation results
      REAL(8) :: ref_r(NNSITE)

      REAL(8) :: freq

      ! cell search results
      INTEGER :: cell(NNSITE)
      ! reference cell search results
      INTEGER :: ref_cell(NNSITE)

      INTEGER :: i, errnums
      INTEGER(4) :: errcode

      errcode = DF_STATUS_OK
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****
      sorder      = DF_PP_LINEAR
      stype       = DF_PP_DEFAULT

      !***** Parameters describing interpolation interval *****
      left        = LLIM_X
      right       = RLIM_X
      nx          = N
      xhint       = DF_NON_UNIFORM_PARTITION

      !***** Parameters describing function *****
      ny = NNY
      yhint = DF_NO_HINT

      !***** Parameters describing spline coefficients storage *****
      scoeffhint  = DF_NO_HINT

      !***** Parameters describing boundary conditions type *****
      !***** No boundary conditions are provided for linear spline *****
      bc_type     = DF_NO_BC

      !***** Parameters describing internal conditions type *****
      !***** No internal conditions are provided for linear spline *****
      ic_type     = DF_NO_IC

      !***** Parameters describing interpolation sites *****
      nsite       = NNSITE
      sitehint    = DF_NON_UNIFORM_PARTITION

      !***** Parameter describing interpolation results storage *****
      rhint       = DF_NO_HINT

      !***** Generate array of uniformly distributed break points *****
      errcode = dUniformRandSortedData( x, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Generate interpolation sites *****
      errcode = dUniformRandData( site, left, right, nsite )
      CALL CheckDfError(errcode)

      !***** Create Data Fitting task *****
      errcode = dfdnewtask1d( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Perform cells search using STD method *****
      method = DF_METHOD_STD
      errcode = dfdsearchcells1d( task, method, nsite, site,             &
     &      sitehint, cell=cell )
      CALL CheckDfError(errcode)

      !***** Edit task parameters for linear spline construction *****
      errcode = dfdeditppspline1d( task, sorder, stype, bc_type,         &
     &      ic_type=ic_type, scoeff=scoeff, scoeffhint=scoeffhint )
      CALL CheckDfError(errcode)

      !***** Generate functions y = sin(2 * Pi * freq * x) *****
      freq = FFREQ
      errcode = dSinDataNotUniformGrid( y, x, freq, nx )
      CALL CheckDfError(errcode)

      !***** Construct linear spline using STD method *****
      comp_type = DF_PP_SPLINE
      method = DF_METHOD_STD
      errcode = dfdconstruct1d( task, comp_type, method )
      CALL CheckDfError(errcode)

      !***** Interpolate using user provided cell indices *****
      comp_type = DF_INTERP_USER_CELL
      method = DF_METHOD_PP
      ndorder = 1
      dorder(1) = 1
      errcode = dfdinterpolate1d( task, comp_type, method, nsite, site,  &
     &      sitehint, ndorder, dorder, r=r, rhint=rhint, cell=cell )
      CALL CheckDfError(errcode)

      !***** Check search results *****
      errcode = dFindCells( nx, x, nsite, site, ref_cell )
      DO i = 1, nsite
        ref_cell(i) = ref_cell(i) - 1
        IF ( ref_cell(i) /= cell(i) ) errnums = errnums+1
      END DO

      !***** Check interpolation results *****
      errcode = dCheckLinInterpRes( nx, x, ny, scoeff, nsite, site,              &
     &  ndorder, dorder, r, ref_r )
      IF (errcode < 0) errnums = errnums+1

      !***** Print results *****
      WRITE (*,901) "Number of break points :        ",nx
      WRITE (*,901) "Number of interpolation sites : ",nsite

      !***** Print given function *****
      WRITE (*,902) "  i   x(i)       y(i)"
      DO i = 1, nx
        WRITE (*,903) i, x(i), y(i)
      END DO

      !***** Print cell search and interpolation results *****
      WRITE (*,902) " Results of cell search and linear interpolation:"
      WRITE (*,902) " Site      Cell index      Interpolation result"
      WRITE (*,906) "                            Obtained"
      WRITE (*,904) "         Expected"
      DO i = 1, nsite
        WRITE (*,905) site(i), cell(i), ref_cell(i), r(i), ref_r(i)
      END DO

      !***** Delete Data Fitting task *****
      errcode = dfdeletetask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,902) "","Error: Computed interpolation results"
        WRITE (*,904) " are incorrect"
        STOP 1
      ELSE
        WRITE (*,902) "","Computed interpolation results"
        WRITE (*,904) " are correct"
      END IF
      STOP 0

  901 FORMAT (A,I0)
  902 FORMAT (/A)
  903 FORMAT (I3,F11.6,F11.6)
  904 FORMAT (99A)
  905 FORMAT (F11.6,I10,I10,F11.6,F11.6)
  906 FORMAT (A,$)

      END PROGRAM
