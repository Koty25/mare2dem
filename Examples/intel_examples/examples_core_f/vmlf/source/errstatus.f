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

!  Content:
!    vmlSetErrStatus/vmlGetErrStatus/vmlClearErrStatus  Example Program Text
!*******************************************************************************

      program MKL_VML_TEST

      include "mkl_vml.f90"

      integer errst
      integer(kind=4) terrst

      print *,"vmlSetErrStatus/vmlGetErrStatus/vmlClearErrStatus",      &
     &        " example program"
      print *,""
      print *,""

      errst=VMLGETERRSTATUS()

      print 10,"Initial value of vmlErrStatus: ",errst,"h"
      call PrintTextVmlErrStatus(errst)

      errst=VML_STATUS_BADMEM
      terrst=VMLSETERRSTATUS(errst)
      errst=VMLGETERRSTATUS()
      print 10,"Value of vmlErrStatus after using vmlSetErrStatus: ",   &
     &         errst,"h"
      call PrintTextVmlErrStatus(errst)

      terrst=VMLCLEARERRSTATUS()
      errst=VMLGETERRSTATUS()
      print 10,"Value of vmlErrStatus after using vmlClearErrStatus: ", &
     &         errst,"h"
      call PrintTextVmlErrStatus(errst)

10    format(A,Z16,A)

      end

      subroutine PrintTextVmlErrStatus(errst)

      include "mkl_vml.f90"

      integer errst

      if(errst.eq.VML_STATUS_OK       ) print *,"VML_STATUS_OK       "
      if(errst.eq.VML_STATUS_BADSIZE  ) print *,"VML_STATUS_BADSIZE  "
      if(errst.eq.VML_STATUS_BADMEM   ) print *,"VML_STATUS_BADMEM   "
      if(errst.eq.VML_STATUS_ERRDOM   ) print *,"VML_STATUS_ERRDOM   "
      if(errst.eq.VML_STATUS_SING     ) print *,"VML_STATUS_SING     "
      if(errst.eq.VML_STATUS_OVERFLOW ) print *,"VML_STATUS_OVERFLOW "
      if(errst.eq.VML_STATUS_UNDERFLOW) print *,"VML_STATUS_UNDERFLOW"

      end
