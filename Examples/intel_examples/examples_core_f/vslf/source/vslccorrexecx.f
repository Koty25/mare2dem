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
!    vslcCorrExecX  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'

      PROGRAM MKL_VSL_TEST
      USE MKL_VSL_TYPE
      USE MKL_VSL
      IMPLICIT NONE
!
      TYPE(VSL_CORR_TASK) task
      INTEGER(KIND=4) status,i
      INTEGER mode,rank
      INTEGER xshape(1),yshape(1),zshape(1)
      INTEGER xstride(1),ystride(1),zstride(1)
      COMPLEX(KIND=4) x(100),y(1000),z(1099)
      LOGICAL ok
!
      DATA rank/1/
      DATA xshape/100/, yshape/1000/, zshape/1099/
      DATA xstride/1/, ystride/1/, zstride/1/
      DATA x/100*0D0/, y/1000*0D0/
!
      ok = .TRUE.
      PRINT *, 'EXAMPLE executing a correlation task'
!
      mode = VSL_CORR_MODE_AUTO
      status = vslccorrnewtaskx(task,mode,rank,xshape,yshape,zshape,    &
     &  x,xstride)
!
      status = vslccorrexecx(task,y,ystride,z,zstride)
!
      IF (status .NE. VSL_STATUS_OK) THEN
         PRINT *, 'ERROR: bad status: ',status
         PRINT *, 'EXAMPLE FAILED'
         ok = .FALSE.
      END IF
!
      DO i=1,1099
         IF (z(i) .NE. 0) THEN
            PRINT *, 'ERROR: wrong result: i=',i,', z(i)=',z(i)
            ok = .FALSE.
         END IF
      END DO
!
      IF (ok) THEN
         PRINT *, 'EXAMPLE PASSED'
         STOP
      ELSE
         PRINT *, 'EXAMPLE FAILED'
         STOP 1
      END IF
!
      END