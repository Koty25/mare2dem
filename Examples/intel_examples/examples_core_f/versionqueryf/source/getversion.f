!===============================================================================
! Copyright 1999-2020 Intel Corporation.
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
!      MKLGetVersionString example program to obtain a string that contains the
!      version information.
!*******************************************************************************

      program getversionstring

      include 'mkl_service.fi'
      character*198    buf

      call mkl_get_version_string(buf)
      write(*,'(/a)') buf

      end
