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
!      C A X P Y  Example Program Text
!*******************************************************************************

      program  CAXPY_MAIN
*
      integer    n, incx, incy
      complex    alpha
      integer    xmax, ymax
      parameter  (xmax=10, ymax=10)
      complex    x(xmax), y(ymax)
      integer    i
*       Intrinsic Functions
      intrinsic  abs
*       External Subroutines
      external   CAXPY, PrintVectorC
*
*       Executable Statements
*
      print*
      print*,'   C A X P Y  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, n, incx, incy
      read*, alpha
      if ( ((1+(n-1)*abs(incx)).gt.xmax).or.
     $     ((1+(n-1)*abs(incy)).gt.ymax) ) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, (x(i),i=1,1+(n-1)*abs(incx))
      read*, (y(i),i=1,1+(n-1)*abs(incy))

*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      print 101, alpha
      call PrintVectorC(0,n,x,incx,'X ')
      call PrintVectorC(0,n,y,incy,'Y ')

*       Call CAXPY subroutine
      call CAXPY(n,alpha,x,incx,y,incy)

      print*
      print*, '     OUTPUT DATA'
      call PrintVectorC(1,n,y,incy,'Y ')
      stop
 100  format(7x,'N=',i2)
 101  format(7x,'ALPHA = (',f4.1,',',f4.1,' )')
 999  stop 1
      end
