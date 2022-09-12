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

      SUBROUTINE SDCOR(h,inch,x,incx,y,incy,nh,nx,iy0,ny,id)
      USE MKL_VSL_TYPE
      USE MKL_VSL
      IMPLICIT NONE
!
      REAL(KIND=4) h(*),x(*),y(*)
      INTEGER(KIND=4) iy0,id
      INTEGER inch,incx,incy,nh,nx,ny
!
      TYPE(VSL_CORR_TASK) task
      INTEGER(KIND=4) status,error
      INTEGER start(1),decimation(1)
      INTEGER mode
!
      start(1) = iy0
      decimation(1) = id
      mode = VSL_CORR_MODE_DIRECT
      status = vslscorrnewtask1d(task,mode,nh,nx,ny)
      status = vslcorrsetstart(task,start)
      status = vslcorrsetdecimation(task,decimation)
      status = vslscorrexec1d(task,h,inch,x,incx,y,incy)
      error  = vslcorrdeletetask(task)
!
      IF (status .NE. VSL_STATUS_OK) THEN
      PRINT *, 'ERROR: sdcor(): bad status=',status
      STOP 1
      END IF
!
      IF (error .NE. 0) THEN
      PRINT *, 'ERROR: sdcor(): failed to destroy the task descriptor'
      STOP 1
      END IF
!
      RETURN
      END
