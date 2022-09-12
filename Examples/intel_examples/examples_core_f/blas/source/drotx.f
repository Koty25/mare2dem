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
!      D R O T  Example Program Text
!*******************************************************************************

      program  DROT_MAIN
*
      integer           n, incx, incy
      double precision  c, s
      integer           xmax, ymax
      parameter        (xmax=10, ymax=10)
      double precision  x(xmax), y(ymax)
      integer           i
*       Intrinsic Functions
      intrinsic         abs
*       External Subroutines
      external          DROT, PrintVectorD
*
*       Executable Statements
*
      print*
      print*,'   D R O T  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, n
      read*, incx, incy
      read*, c, s
      if ( ((1+(n-1)*abs(incx)).gt.xmax).or.
     $     ((1+(n-1)*abs(incy)).gt.ymax) ) then
          print*, ' Insufficient memory for arrays'
          goto 999
      end if
      read*, (x(i),i=1,1+(n-1)*abs(incx))
      read*, (y(i),i=1,1+(n-1)*abs(incy))
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      print 101, c, s
      call PrintVectorD(0,n,x,incx,'X ')
      call PrintVectorD(0,n,y,incy,'Y ')
*
*       Call DROT subroutine
      call DROT(n,x,incx,y,incy,c,s)
*
      print*
      print*, '     OUTPUT DATA'
      call PrintVectorD(1,n,x,incx,'X ')
      call PrintVectorD(1,n,y,incy,'Y ')

      stop
 100  format(7x,'N=',i2)
 101  format(7x,'C=',f5.2,'  S=',f5.2)
 999  stop 1
      end
