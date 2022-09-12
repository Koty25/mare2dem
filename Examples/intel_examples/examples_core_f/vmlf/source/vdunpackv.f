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
!    vdUnpackV  Example Program Text
!*******************************************************************************

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      real(kind=8) dA(10)
      real(kind=8) dB1(10)
      real(kind=8) dB2(10)
      integer i
      integer ia(10)

      i=0

      do i=1,10
            dA(i)=i/1.0
            dB1(i)=0.0
            dB2(i)=0.0
            ia(i)=10-i
      end do

      call VDPACKV(10,dA,ia,dB1)
      call VDUNPACKV(10,dB1,dB2,ia)

      print *,"vdUnpackV test/example program"
      print *,""
      print *,""

      print *,"           Before packing             After packing   ", &
     &        "       After Unpacking"
      print *,"======================================================", &
     &        "========================="
      do i=1,10
            print 10,dA(i),dB1(i),dB2(i)
      end do

10    format (F25.13,F25.13,F25.13)

      end
