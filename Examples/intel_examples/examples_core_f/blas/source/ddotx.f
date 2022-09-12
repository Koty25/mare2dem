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
!      D D O T  Example Program Text
!*******************************************************************************

      program  DDOT_MAIN
*
      integer          n, incx, incy
      integer          xmax, ymax
      parameter       (xmax=20, ymax=20)
      double precision sx(xmax), sy(ymax), res
      integer          i
*       Intrinsic Functions
      intrinsic        abs
*       External Subroutines
      external         DDOT, PrintVectorD
      double precision DDOT
*
*       Executable Statements
*
      print*
      print*,'   D D O T  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, n
      read*, incx, incy
      if ( ((1+(n-1)*abs(incx)).gt.xmax).or.
     $     ((1+(n-1)*abs(incy)).gt.ymax) ) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, (sx(i),i=1,1+(n-1)*abs(incx))
      read*, (sy(i),i=1,1+(n-1)*abs(incy))
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      call PrintVectorD(0,n,sx,incx,'X ')
      call PrintVectorD(0,n,sy,incy,'Y ')
*
*       Call DDOT subroutine
      res = DDOT(n,sx,incx,sy,incy)
*
      print*
      print*, '     OUTPUT DATA'
      print 101, res

      stop
 100  format(7x,'N=',i2)
 101  format(7x,'DDOT = ',f8.3)
 999  stop 1
      end