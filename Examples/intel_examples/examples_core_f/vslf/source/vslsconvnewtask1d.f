!===============================================================================
! Copyright 2003-2020 Intel Corporation.
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
!    vslsConvNewTask1D  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'

      PROGRAM MKL_VSL_TEST
      USE MKL_VSL_TYPE
      USE MKL_VSL
      IMPLICIT NONE
!
      TYPE(VSL_CONV_TASK) task
      INTEGER(KIND=4) status
      INTEGER mode
      INTEGER xshape,yshape,zshape
!
      PRINT *, 'EXAMPLE creating a new task'
!
      mode = VSL_CONV_MODE_AUTO
      xshape = 100
      yshape = 1000
      zshape = (xshape-1) + (yshape-1) + 1
      status = vslsconvnewtask1d(task,mode,xshape,yshape,zshape)
!
      IF (status .NE. VSL_STATUS_OK) THEN
         PRINT *, 'ERROR: bad status: ',status
         PRINT *, 'EXAMPLE FAILED'
         STOP 1
      ELSE
         PRINT *, 'EXAMPLE PASSED'
         STOP
      END IF
!
      END
