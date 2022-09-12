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
!    Querying Data Fitting task parameters Example Program Text
!*******************************************************************************

      include 'mkl_df.f90'
      include "errcheck.inc"
      include "rescheck.inc"
      include "generatedata.inc"

      PROGRAM MKL_DF_TEST

      USE MKL_DF_TYPE
      USE MKL_DF

      ! number of breakpoints
      INTEGER,PARAMETER :: N       = 10
      ! number of datasets to interpolate
      INTEGER,PARAMETER :: NNY     = 2

      INTEGER(8),PARAMETER :: SIZEOF_REAL8 = 8
      INTEGER(8),PARAMETER :: INT32_MASK   = INT(Z'00000000FFFFFFFF', 8)

      ! Data Fitting task descriptor
      TYPE (DF_TASK) task

      ! number of break points
      INTEGER :: nx
      INTEGER :: nx_attr
      ! additional info about break points
      INTEGER :: xhint
      INTEGER :: xhint_attr
      ! number of functions
      INTEGER :: ny
      INTEGER :: ny_attr
      ! functions storage format
      INTEGER :: yhint
      INTEGER :: yhint_attr

      ! array of break points
      REAL(8) :: x(N)
      INTEGER :: x_attr
      ! function values
      REAL(8) :: y(N*NNY)
      INTEGER :: y_attr

      INTEGER :: idx

      ! parameters to be queried
      INTEGER :: nx_ret
      INTEGER(8) :: y_ret
      INTEGER(8) :: y1_ret

      INTEGER(8) :: loc_y
      INTEGER(8) :: loc_y1

      INTEGER :: errnums
      INTEGER(4) :: errcode

      !***** Initializing parameters for Data Fitting task *****

      errcode = 0

      !***** Parameters describing interpolation interval *****
      nx        = N
      xhint     = DF_UNIFORM_PARTITION

      x(1) = -1.0d0
      x(2) =  1.0d0

      !***** Parameters describing functions *****
      ny         = NNY
      yhint      = DF_MATRIX_STORAGE_ROWS

      !***** Create Data Fitting task *****
      errcode = dfdNewTask1D( task, nx, x, xhint, ny, y, yhint )
      CALL CheckDfError(errcode)

      !***** Query Data Fitting task parameters *****

      !***** Query value *****
      nx_attr = DF_NX
      errcode = dfiQueryVal( task, nx_attr, nx_ret )
      CALL CheckDfError(errcode)

      y_ret  = 0
      y1_ret = 0

      !***** Query pointer *****
      y_attr = DF_Y
      errcode = dfdQueryPtr( task, y_attr, y_ret )
      CALL CheckDfError(errcode)

      !***** Query pointer by index *****
      idx = 1
      errcode = dfdQueryIdxPtr( task, y_attr, idx, y1_ret )
      CALL CheckDfError(errcode)

      !***** Check queried parameters *****
      errnums = 0
      if ( nx /= nx_ret ) then
        errnums = errnums + 1
      end if

      loc_y  = loc(y)
      loc_y1 = loc_y + nx * SIZEOF_REAL8

      if (BIT_SIZE(loc(y)) == 32) then
      !***** Adjust pointer values on IA32 architecture *****
        loc_y = IAND(loc_y,INT32_MASK)
        y_ret = IAND(y_ret,INT32_MASK)
        loc_y1 = IAND(loc_y1,INT32_MASK)
        y1_ret = IAND(y1_ret,INT32_MASK)
      end if

      if ( loc_y /= y_ret ) then
        errnums = errnums + 1
      end if

      if ( loc_y1 /= y1_ret ) then
        errnums = errnums + 1
      end if

      WRITE (*,900,ADVANCE='NO') "                                  "
      WRITE (*,900) "          Expected          Obtained"
      WRITE (*,902) "Number of break points          : ", nx, nx_ret
      WRITE (*,903) "Address of function Y           : ", loc_y, y_ret
      WRITE (*,903) "Address of 1-st coordinate of Y : ", loc_y1, y1_ret

      !***** Delete Data Fitting task *****
      errcode = dfDeleteTask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,901) "Error: Not all requested parameters are correct"
        STOP 1
      ELSE
        WRITE (*,901) "All requested parameters are correct"
      END IF
      STOP 0

900   FORMAT (99A)
901   FORMAT (//A)
902   FORMAT (A,I18,I18)
903   FORMAT (A,Z18,Z18)
      END PROGRAM
