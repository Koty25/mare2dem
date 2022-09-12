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
! Example of 2-dimension linear correlation operation on single precision data
!*******************************************************************************/
      INCLUDE 'mkl_vsl.fi'

      PROGRAM vslscorrelation_2d_fft
      USE MKL_VSL_TYPE
      USE MKL_VSL
      IMPLICIT NONE
C
      INTEGER status,i,j
C
      TYPE(VSL_CORR_TASK) task
      INTEGER mode,rank,xshape(2),yshape(2),zshape(2)
      INTEGER xstride(2),ystride(2),zstride(2)
      REAL*4 x(3,2),y(4,3),z(6,4),e(6,4)
      DATA x/1,1,1, 1,1,1/
      DATA y/1,1,1,1, 1,1,1,1, 1,1,1,1/
      DATA z/0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0/
      DATA e/1,2,3,3,2,1, 2,4,6,6,4,2, 2,4,6,6,4,2, 1,2,3,3,2,1/
      DATA xstride/1,3/, ystride/1,4/, zstride/1,6/
      DATA xshape/3,2/,  yshape/4,3/,  zshape/6,4/
      DATA rank/2/

      mode = VSL_CORR_MODE_FFT

C
C  Create task descriptor (create descriptor of problem)
C
      status = vslscorrnewtask(task,mode,rank,xshape,yshape,zshape)
      IF (status .NE. VSL_STATUS_OK) THEN
          PRINT *,"ERROR: creation of job failed, exit with ", status
          STOP 1
      ENDIF
C
C  Execute task (compute result of linear correlation)
C
      status = vslscorrexec(task,x,xstride,y,ystride,z,zstride)
      IF (status .NE. VSL_STATUS_OK) THEN
          PRINT *,"ERROR: job status bad, exit with", status
          STOP 1
      ENDIF
C
C  Delete task object (delete descriptor)
C
      status = vslcorrdeletetask(task)
      IF (status .NE. VSL_STATUS_OK) THEN
          PRINT *,"ERROR: failed to delete task, exit with", status
          STOP 1
      ENDIF
C
C  Check results for correctness:
C
      DO j=1,zshape(2)
        DO i=1,zshape(1)
            IF (abs(z(i,j)-e(i,j)) .GT. e(i,j)*1e-5) THEN
                PRINT *, "ERROR:"
                PRINT *, "        i :", i
                PRINT *, "        j :", j
                PRINT *, "    z(i,j):", z(i,j)
                PRINT *, "  expected:", e(i,j)
                PRINT *, "EXAMPLE FAILED"
                STOP 1
             END IF
        END DO
      END DO
C
      PRINT *, "EXAMPLE PASSED"
      STOP 0
      END
