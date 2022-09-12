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
!      I D A M I N  Example Program Text
!*******************************************************************************

      program  IDAMIN_MAIN
*
      integer          n, incx
      integer          xmax
      parameter       (xmax=10)
      double precision x(xmax)
      integer          i, index
*       Intrinsic Functions
      intrinsic        abs
*       External Subroutines
      external         IDAMIN, PrintVectorD
      integer          IDAMIN
*
*       Executable Statements
*
      print*
      print*,'   I D A M I N  EXAMPLE PROGRAM'
*       Read input data from input file
      read*
      read*, n, incx
      if ((1+(n-1)*abs(incx)).gt.xmax) then
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
*       Call IDAMIN subroutine
      index = IDAMIN(n,x,incx)
*
      print*
      print*, '     OUTPUT DATA'
      print 101, index

      stop
 100  format(7x,'N=',i2)
 101  format(10x,'IDAMIN = ',i3)
 999  stop 1
      end
