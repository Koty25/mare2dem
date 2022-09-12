!===============================================================================
! Copyright 2001-2020 Intel Corporation.
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

! Content:
! Example of 1-dimension linear convolution operation on double precision data
!*******************************************************************************/

      INCLUDE 'mkl_vsl.fi'

      PROGRAM vsldconvolution_1d_auto
      USE MKL_VSL_TYPE
      USE MKL_VSL
      IMPLICIT NONE
C
      INTEGER status,i
C
      TYPE(VSL_CONV_TASK) task
      INTEGER mode,xshape,yshape,zshape
      REAL*8 x(4),y(8),z(11),e(11)
      DATA x/1,2,3,4/
      DATA y/11,12,13,14,15,16,17,18/
      DATA z/0,0,0,0,0,0,0,0,0,0,0/
      DATA e/11,34,70,120,130,140,150,160,151,122,72/
      DATA xshape/4/, yshape/8/, zshape/11/

      mode = VSL_CONV_MODE_AUTO

C
C  Create task descriptor (create descriptor of problem)
C
      status = vsldconvnewtask1d(task,mode,xshape,yshape,zshape)
      IF (status .NE. VSL_STATUS_OK) THEN
          PRINT *, "ERROR: creation of job failed, exit with ", status
          STOP 1
      ENDIF
C
C  Execute task (Calculate 1 dimension convolution of two arrays)
C
      status = vsldconvexec1d(task,x,1,y,1,z,1)
      IF ( status .NE. VSL_STATUS_OK) THEN
          PRINT *,"ERROR: job status bad, exit with", status
          STOP 1
      ENDIF
C
C  Delete task object (delete descriptor of problem)
C
      status = vslconvdeletetask(task)
      IF (status .NE. VSL_STATUS_OK) THEN
          PRINT *,"ERROR: failed to delete task, exit with", status
          STOP 1
      ENDIF
C
C  Check results for correctness:
C
      DO i=1,zshape
         IF (abs(z(i)-e(i)) .GT. e(i)*1d-10) THEN
            PRINT *, "ERROR:"
            PRINT *, "      i :", i
            PRINT *, "    z(i):", z(i)
            PRINT *, "expected:", e(i)
            PRINT *, "EXAMPLE FAILED"
            STOP 1
         END IF
      END DO
C
      PRINT *, "EXAMPLE PASSED"
      STOP 0
      END
