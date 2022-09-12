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
!    vslGetNumRegBrngs  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      TYPE (VSL_STREAM_STATE) :: stream
      integer(kind=4) brngsObt,brngsExp
      integer(kind=4) errcode

      brngsExp=16
      brngsObt=0

      brngsObt = vslgetnumregbrngs()

!     ***** Printing results *****
      print *,"Sample of vslGetNumRegBrngs"
      print *,"-------------------------------"
      print *,""
      if (brngsObt.NE.brngsExp) then
        print 13,"Error: returned value ", brngsObt," is incorrect      &
     &           (expected ",brngsExp," )!"
        stop 1
      else
        print 14,"Returned ", brngsObt," as expected."
      end if

13    format(A,I3,A,I3,A)
14    format(A,I3,A)

      end
