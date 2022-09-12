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
!      D S P R 2  Example Program Text
!*******************************************************************************

      program   DSPR2_MAIN
*
      character*1      uplo
      integer          n, incx, incy
      integer          apmax, xmax, ymax
      parameter        (apmax=15, xmax=5, ymax=5)
      double precision alpha
      double precision ap(apmax), x(xmax), y(ymax)
      integer          i
*       Intrinsic Functions
      intrinsic        abs
*       External Subroutines
      external         DSPR2, PrintVectorD
*
*       Executable Statements
*
      print*
      print*, '   D S P R 2  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, n
      read*, incx, incy
      read*, alpha
      read 100, uplo
      if ( ((1+(n-1)*abs(incx)).gt.xmax).or.
     $     ((1+(n-1)*abs(incy)).gt.ymax).or.
     $     (n*(n+1))/2.gt.apmax ) then
          print*, ' Insufficient memory for arrays'
          goto 999
      end if
      read*, (x(i),i=1, 1+(n-1)*abs(incx))
      read*, (y(i),i=1, 1+(n-1)*abs(incy))
      read*, (ap(i),i=1,(n*(n+1))/2)
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 101, n
      print 102, alpha
      print 103, uplo
      call PrintVectorD(0,n,x,incx,'X ')
      call PrintVectorD(0,n,y,incy,'Y ')
      call PrintVectorD(1,(n*(n+1))/2,ap,1,'AP')
*
*      Call DSPR2 subroutine
      call DSPR2( uplo, n, alpha, x, incx, y, incy, ap)
*
      print*
      print*, '     OUTPUT DATA'
      call PrintVectorD(1,(n*(n+1))/2,ap,1,'AP')

      stop
 100  format(a1)
 101  format(7x,'N=',i1)
 102  format(7x,'ALPHA=',f4.1)
 103  format(7x,'UPLO=',a1)
 999  stop 1
      end