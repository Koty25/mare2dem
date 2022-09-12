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
!    vzUnpackM  Example Program Text
!*******************************************************************************

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      complex(kind=8) zA(10)
      complex(kind=8) zB1(10)
      complex(kind=8) zB2(10)
      integer i
      integer ma(10)

      i=0

      do i=1,10
            zA(i)=cmplx(i/1.0,i/1.0)
            zB1(i)=cmplx(0.0,0.0)
            zB2(i)=cmplx(0.0,0.0)
            ma(i)=iand(i,1)
      end do

      call vzPackM(10,zA,ma,zB1)
      call VZUNPACKM(10,zB1,zB2,ma)

      print *,"vzUnpackM test/example program"
      print *,""
      print *,""

      print *,"     Before packing          After packing   ",          &
     &        "       After Unpacking"
      print *,"======================================================", &
     &        "========================="
      do i=1,10
            print 10,zA(i),zB1(i),zB2(i)
      end do

10    format(E12.3,E12.3,E12.3,E12.3,E12.3,E12.3)

      end
