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
!      D T P S V  Example Program Text
!*******************************************************************************

      program   DTPSV_MAIN
*
      character*1      uplo, trans, diag
      integer          n, incx
      integer          apmax, xmax
      parameter        (apmax=15, xmax=10)
      double precision ap(apmax), x(xmax)
      integer          i
*       Intrinsic Functions
      intrinsic        abs
*       External Subroutines
      external         DTPSV, PrintVectorD
*
*       Executable Statements
*
      print*
      print*, '   D T P S V  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, n
      read*, incx
      read 100, uplo, trans, diag
      if ( (1+(n-1)*abs(incx)).gt.xmax.or.
     $     (n*(n+1))/2.gt.apmax ) then
          print*, ' Insufficient memory for arrays'
          goto 999
      end if
      read*, (x(i),i=1, 1+(n-1)*abs(incx))
      read*, (ap(i),i=1,(n*(n+1))/2)
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 101, n
      print 102, uplo, trans, diag
      call PrintVectorD(0,n,x,incx,'X ')
      call PrintVectorD(1,(n*(n+1))/2,ap,1,'AP')
*
*      Call DTPSV subroutine
      call DTPSV(uplo, trans, diag, n, ap, x, incx)
*
      print*
      print*, '     OUTPUT DATA'
      call PrintVectorD(1,n,x,incx,'X ')

      stop
 100  format(3(a1,1x))
 101  format(7x,'N=',i1)
 102  format(7x,'UPLO=',a1,'  TRANS=',a1,'  DIAG=',a1)
 999  stop 1
      end
