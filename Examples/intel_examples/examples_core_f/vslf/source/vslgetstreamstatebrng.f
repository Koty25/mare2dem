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
!    vslGetStreamStateBrng  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      TYPE (VSL_STREAM_STATE) :: stream
      integer seed,brngExp
      integer(kind=4) brngObt
      integer(kind=4) errcode

      brngExp=VSL_BRNG_WH+127
      brngObt=0

      seed=7777777

!     ***** Initialize stream *****
      errcode=vslnewstream  ( stream, brngExp, seed )
      call CheckVslError(errcode)

      brngObt = vslgetstreamstatebrng( stream )

!     ***** Printing results *****
      print *,"Sample of vslGetStreamStateBrng"
      print *,"-------------------------------"
      print *,""
      print *,"Parameters:"
      print 11,"    seed =   ",seed
      print 11,"    brng =   ",brngExp

      print *,""
      if (brngObt.NE.brngExp) then
        print 13,"Error: returned value ", brngObt," is incorrect       &
     &            (expected ",brngExp," )!"
        stop 1
      else
        print 14,"Returned ", brngObt," as expected."
      end if

!     ***** Deinitialize *****
      errcode=vsldeletestream( stream )
      call CheckVslError(errcode)

11    format(A,I10)
13    format(A,I10,A,I2,A)
14    format(A,I10,A)

      end
