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
!      C S W A P  Example Program Text
!*******************************************************************************

      program CSWAP_MAIN
*
      integer      n, incx, incy
      integer      xmax, ymax
      parameter    (xmax=10,ymax=10)
      complex      x(xmax), y(ymax)
      integer      i
*       External Subroutines
      external     CSWAP, PrintVectorC
*
*   Executable Statementcs
*
      print*
      print*,'   C S W A P  EXAMPLE PROGRAM'
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
      read*, (x(i),i=1,1+(n-1)*abs(incx))
      read*, (y(i),i=1,1+(n-1)*abs(incy))
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      call PrintVectorC(0,n,x,incx,'X ')
      call PrintVectorC(0,n,y,incy,'Y ')
*
*      Call CSWAP subroutine
      call CSWAP(n,x,incx,y,incy)
*
      print*
      print*, '     OUTPUT DATA'
      call PrintVectorC(1,n,x,incx,'X ')
      call PrintVectorC(1,n,y,incy,'Y ')

      stop
 100  format(7x,'N=',i2)
 999  stop 1
      end
