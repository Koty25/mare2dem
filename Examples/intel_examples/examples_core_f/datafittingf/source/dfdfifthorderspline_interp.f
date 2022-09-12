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
!    Interpolation, and computation of the 2nd derivative for the spline
!  of the fifth order Example Program Text
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
      INTEGER,PARAMETER :: N        = 4
      ! number of functions
      INTEGER,PARAMETER :: NNY      = 1
      !  number of sites for spline-based interpolation in one block
      INTEGER,PARAMETER :: NNSITE_BLOCK  = 10
      ! number of blocks
      INTEGER,PARAMETER :: NNBLOCKS  = 4
      ! size of array describing derivative orders to compute
      INTEGER,PARAMETER :: NNDORDER  = 3
      ! number of derivatives to compute
      INTEGER,PARAMETER :: NNDER     = 2
      ! spline order
      INTEGER,PARAMETER :: SPLINE_ORDER = 5
      ! total number of spline coefficients
      INTEGER,PARAMETER :: NNSCOEFF = (NNY*(N-1)*SPLINE_ORDER)
      ! total number of interpolation sites
      INTEGER,PARAMETER :: NNSITE   = (NNSITE_BLOCK*NNBLOCKS)
      ! left  limit of interpolation interval
      REAL(8),PARAMETER :: LLIM_X   = 0.0d0
      ! right limit of interpolation interval
      REAL(8),PARAMETER :: RLIM_X   = 1.5d0

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
      ! number of interpolaton sites in one block
      INTEGER :: nsite_bl
      ! number of blocks
      INTEGER :: nblocks
      ! number of interpolation sites
      INTEGER :: nsite
      ! additional info about interpolation sites
      INTEGER :: sitehint
      !  size of array describing derivative orders
      INTEGER :: ndorder
      ! interpolation results storage format
      INTEGER :: rhint
      ! left limit of the interpolation interval
      REAL(8) :: left
      ! right limit of the interpolation interval
      REAL(8) :: right
      ! array of break points
      REAL(8) :: x(N)
      ! function values
      REAL(8) :: y(NNY*N)
      ! array of spline coefficients
      REAL(8) :: scoeff(NNSCOEFF)
      ! array of interpolation sites
      REAL(8),TARGET :: site(NNSITE)
      ! spline evaluation and integration results
      REAL(8),TARGET :: r((NNDORDER+1)*NNSITE)
      ! type of calculations
      INTEGER :: type
      ! method that is used to perform calculations
      INTEGER :: method

      INTEGER :: i,j,nder,errnums
      INTEGER(4) :: errcode
      REAL(8) :: ref_r((NNDORDER+1) *NNSITE)
      REAL(8),POINTER :: site_ptr(:),r_ptr(:)
      INTEGER :: dorder(NNDORDER)

      errcode = 0
      errnums = 0

      !***** Initializing parameters for Data Fitting task *****

      sorder = SPLINE_ORDER
      stype  = DF_PP_NATURAL

      !***** Parameters describing interpolation interval *****
      left  = LLIM_X
      right = RLIM_X
      nx = N
      xhint = DF_QUASI_UNIFORM_PARTITION

      !***** Parameters describing function *****
      ny = NNY
      yhint = DF_NO_HINT

      !***** Parameters describing spline coefficients storage *****
      scoeff = (/0.0D0,  1.0D0,  0.0D0,  0.0D0,  0.0D0,                  &
     &           0.5D0, -1.0D0,  0.0D0,  0.0D0,  0.0D0,                  &
     &           0.0D0,  1.0D0,  0.0D0,  0.0D0,  0.0D0/)
      scoeffhint = DF_NO_HINT

      !***** Parameters describing interpolation sites *****

      nsite = NNSITE
      sitehint = DF_NON_UNIFORM_PARTITION

      !**** Parameter describing interpolation results storage *****
      rhint = DF_MATRIX_STORAGE_ROWS

      !**** Parameter describing array for derivative orders *****
      ndorder = NNDORDER
      dorder = (/1,0,1/)

      !***** Number of derivatives to compute *****
      nder = NNDER

      !***** Generate array of uniformly distributed break points *****
      errcode = dUniformData( x, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Generate interpolation sites ****
      errcode = dUniformRandData( site, left, right, nsite )
      CALL CheckDfError(errcode)

      !***** Create Data Fitting task *****
      errcode = dfdnewtask1d( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Register fifth order spline coefficients in Data Fitting task *****
      errcode = dfdeditppspline1d( task, sorder, stype, scoeff=scoeff,   &
     &                             scoeffhint=scoeffhint )
      CALL CheckDfError(errcode)

      !***** Interpolate, and compute 2nd order derivatives
      !      using STD method *****
      nblocks = NNBLOCKS
      nsite_bl = NNSITE_BLOCK

      site_ptr => site
      r_ptr => r

      type = DF_INTERP
      method = DF_METHOD_PP

      DO i = 1, nblocks
        errcode = dfdinterpolate1d( task, type, method,                  &
     &       nsite_bl, site_ptr, sitehint, ndorder,                      &
     &       dorder, r=r_ptr, rhint=rhint )
        CALL CheckDfError(errcode)

        site_ptr => site_ptr(nsite_bl+1:)
        r_ptr => r_ptr((nder*nsite_bl)+1:)
      END DO

      !***** Check results of interpolation *****
      errcode = dCheckFifthOrderInterpRes( nx, x, ny, scoeff,            &
     &    nsite, site, ndorder, dorder, r, ref_r )
      IF (errcode < 0) errnums = errnums+1

      !***** Print results *****
      WRITE (*,901) "Number of break points : ",nx

      !***** Print break points *****
      WRITE (*,902) " i    x(i)"
      DO j = 1, nx
        WRITE (*,903) " ",j," ",x(j)
      END DO

      !***** Print spline coefficients *****
      WRITE (*,904)                                                      &
     & "Coefficients are calculated for a polynomial of the form:"
      WRITE (*,905) "Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 +      &
     & Di*(x - x(i))^3 + Ei*(x - x(i))^4"
      WRITE (*,905) "    where x(i) <= x < x(i+1)"
      WRITE (*,902) "Spline coefficients:"
      WRITE (*,905)                                                      &
     & " i   Ai            Bi            Ci            Di            Ei"
      DO j = 0,nx-1-1
        WRITE (*,906) " ",j," ",scoeff(sorder*j+1),"   ",                &
     &  scoeff(sorder*j + 2),"   ",scoeff(sorder*j + 3),"   ",           &
     &  scoeff(sorder*j + 4),"   ",scoeff(sorder*j + 5)
      END DO

      !***** Print interpolation results ******
      WRITE (*,902) "Results of interpolation:"
      WRITE (*,915) "                    Function value"
      WRITE (*,905) "             Second derivative"
      WRITE (*,915) "    Sites         Obtained     Expected"
      WRITE (*,905) "      Obtained    Expected"
      DO i = 0,nsite-1
        WRITE (*,907,ADVANCE='NO') " ",site(i+1)
        DO j = 1, nder
          WRITE (*,908,ADVANCE='NO') "   ",r(i*nder + j),"  ",           &
     &                                 ref_r(i*nder + j)
        END DO
        WRITE (*,905) ""
      END DO

      !***** Delete Data Fitting task *****
      errcode = dfdeletetask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,909) "","Error: Computed interpolation results"
        WRITE (*,905) " are incorrect"
        STOP 1
      ELSE
        WRITE (*,909) "","Computed interpolation results"
        WRITE (*,905) " are correct"
      END IF
      STOP 0

901   FORMAT (A,I0)
902   FORMAT (/A)
903   FORMAT (A,I1,A,SP,F11.6)
904   FORMAT (/A)
905   FORMAT (99A)
915   FORMAT (A,$)
906   FORMAT (A,I1,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,F11.6,A,SP,     &
     &        F11.6)
907   FORMAT (A,SP,F11.6)
908   FORMAT (A,SP,F11.6,A,SP,F11.6)
909   FORMAT (//A)

      END PROGRAM
