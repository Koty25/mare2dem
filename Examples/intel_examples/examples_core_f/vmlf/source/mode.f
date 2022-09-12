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
!    vmlSetMode/vmlGetMode  Example Program Text
!*******************************************************************************

      program MKL_VML_TEST

      include "mkl_vml.f90"

      integer mode
      integer tmode

      mode=0

      print *,"vmlSetMode/vmlGetMode example program"
      print *,""
      print *,""

      mode=VMLGETMODE()
      print 10,"Default value of vmlMode: ",mode,"h"
      call PrintTextVmlMode(mode)

      mode=IOR(VML_LA,IOR(VML_FLOAT_CONSISTENT,VML_ERRMODE_IGNORE))
      tmode=VMLSETMODE(mode)
      mode=VMLGETMODE()
      print 10,"Value of vmlMode after using vmlSetMode: ",mode,"h"
      call PrintTextVmlMode(mode)

10    format(A,Z4,A)

      end


      subroutine PrintTextVmlMode(input_mode)

      include "mkl_vml.f90"

      integer input_mode

      integer(kind=4) mode
      integer(kind=4) mask

      mode = INT(input_mode, kind=4)
      mask = VML_ACCURACY_MASK

      if(IAND(mode,mask).eq.1) print *,"VML_LA"
      if(IAND(mode,mask).eq.2) print *,"VML_HA"

      mask=VML_FPUMODE_MASK
      if(IAND(mode,mask).eq.0) print *,"VML_DEFAULT_PRECISION"
      if(IAND(mode,mask).eq.16) print *,"VML_FLOAT_CONSISTENT "
      if(IAND(mode,mask).eq.32) print *,"VML_DOUBLE_CONSISTENT"
      if(IAND(mode,mask).eq.48) print *,"VML_RESTORE"

      mask=VML_ERRMODE_MASK
      if(IAND(mode,VML_ERRMODE_IGNORE).eq.VML_ERRMODE_IGNORE)           &
     &  print *,"VML_ERRMODE_IGNORE"
      if(IAND(mode,VML_ERRMODE_NOERR).eq.VML_ERRMODE_NOERR)             &
     &  print *,"VML_ERRMODE_NOERR"
      if(IAND(mode,VML_ERRMODE_ERRNO).eq.VML_ERRMODE_ERRNO)             &
     &  print *,"VML_ERRMODE_ERRNO"
      if(IAND(mode,VML_ERRMODE_STDERR).eq.VML_ERRMODE_STDERR)           &
     &  print *,"VML_ERRMODE_STDERR"
      if(IAND(mode,VML_ERRMODE_EXCEPT).eq.VML_ERRMODE_EXCEPT)           &
     &  print *,"VML_ERRMODE_EXCEPT"
      if(IAND(mode,VML_ERRMODE_CALLBACK).eq.VML_ERRMODE_CALLBACK)       &
     &  print *,"VML_ERRMODE_CALLBACK"

      end
