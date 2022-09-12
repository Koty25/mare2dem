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
!    vcPackV  Example Program Text
!*******************************************************************************

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      complex(kind=4) cA(10)
      complex(kind=4) cB1(10)
      integer ia(10)
      integer i

      i=0

      do i=1,10
            cA(i)=cmplx(i/1.0,i/1.0)
            cB1(i)=cmplx(0.0,0.0)
            ia(i)=10-i
      end do

      call VCPACKV(10,cA,ia,cB1)

      print *,"vcPackV test/example program"
      print *,""
      print *,""

      print *,"           Before packing             After packing"
      print *,"======================================================", &
     &        "========================="
      do i=1,10
            print 10,cA(i),cB1(i)
      end do

10    format(E15.7,E15.7,E15.7,E15.7)

      end
