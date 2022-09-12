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

      include 'mkl_vsl.f90'

      SUBROUTINE SCONF(init,                                            &
     &                 h,inc1h,x,inc1x,inc2x,y,inc1y,inc2y,             &
     &                 nh,nx,m,iy0,ny,                                  &
     &                 aux1,naux1,aux2,naux2)
      USE MKL_VSL_TYPE
      USE MKL_VSL
      IMPLICIT NONE
!
      INTEGER(KIND=4) init,inc2x,inc2y,m,iy0
      INTEGER inc1x,inc1y,inc1h,nh,nx,ny
      REAL(KIND=4) h(*),x(inc2x,*),y(inc2y,*)
      INTEGER(KIND=4) naux1,naux2
      REAL aux1(*),aux2(*)
!
      TYPE(VSL_CONV_TASK) task
      INTEGER(KIND=4) status,error,i
      INTEGER start(1)
      INTEGER mode
!
!     ignore aux1, aux2:
      IF (init .NE. 0) RETURN
!
      start(1) = iy0
      mode = VSL_CONV_MODE_FFT
      status = vslsconvnewtaskx1d(task,mode,nh,nx,ny,h,inc1h)
      status = vslconvsetstart(task,start)
!
!     task is implicitly committed at i==1
      DO i=1,m
      status = vslsconvexecx1d(task, x(1,i),inc1x, y(1,i),inc1y)
      END DO
!
      error  = vslconvdeletetask(task)
!
      IF (status .NE. VSL_STATUS_OK) THEN
      PRINT *, 'ERROR: sconf(): bad status=',status
      STOP 1
      END IF
!
      IF (error .NE. 0) THEN
      PRINT *, 'ERROR: sconf(): failed to destroy the task descriptor'
      STOP 1
      END IF
!
      RETURN
      END
