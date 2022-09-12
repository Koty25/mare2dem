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
!      D N R M 2  Example Program Text
!*******************************************************************************

      program  DNRM2_MAIN
*
      integer          n, incx
      integer          xmax
      parameter       (xmax=10)
      double precision x(xmax), res
      integer          i
*       Intrinsic Functions
      intrinsic        abs
*       External Subroutines
      external         DNRM2, PrintVectorD
      double precision DNRM2
*
*       Executable Statements
*
      print*
      print*,'   D N R M 2  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, n, incx
      if ( (1+(n-1)*abs(incx)).gt.xmax) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, (x(i),i=1,1+(n-1)*abs(incx))
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      call PrintVectorD(0,n,x,incx,'X ')
*
*       Call DNRM2 subroutine
      res = DNRM2(n,x,incx)
*
      print*
      print*, '     OUTPUT DATA'
      print 101, res

      stop
 100  format(7x,'N=',i2)
 101  format(10x,'DNRM2 = ',f8.3)
 999  stop 1
      end
