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
!    vsUnpackM  Example Program Text
!*******************************************************************************

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      real(kind=4) sA(10)
      real(kind=4) sB1(10)
      real(kind=4) sB2(10)
      integer i
      integer ma(10)

      i=0

      do i=1,10
            sA(i)=i/1.0
            sB1(i)=0.0
            sB2(i)=0.0
            ma(i)=iand(i,1)
      end do

      call vsPackM(10,sA,ma,sB1)
      call VSUNPACKM(10,sB1,sB2,ma)

      print *,"vsUnpackM test/example program"
      print *,""
      print *,""

      print *,"           Before packing             After packing   ", &
     &        "       After Unpacking"
      print *,"======================================================", &
     &        "========================="
      do i=1,10
            print 10,sA(i),sB1(i),sB2(i)
      end do

10    format (F25.13,F25.13,F25.13)

      end
