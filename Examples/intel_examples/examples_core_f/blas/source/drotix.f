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
!      D R O T I  Example Program Text
!*******************************************************************************

      program  DROTI_MAIN
*
      integer          n
      double precision c, s
      integer          xmax, ymax, imax
      parameter        (xmax=10, ymax=10, imax=10)
      double precision x(xmax), y(ymax)
      integer          indx(imax)
      integer          i, indmax
*       External Subroutines
      external         DROTI, MaxValue, PrintVectorD
*
*       Executable Statements
*
      print*
      print*,'   D R O T I  EXAMPLE PROGRAM'
*       Read input data from input file
      read*
      read*, n
      read*, c, s
      if (n.gt.imax.or.n.gt.xmax) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, (x(i),i=1,n)
      read*, (indx(i),i=1,n)
      indmax = MaxValue(n,indx)
      if (indmax.gt.ymax) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, (y(i),i=1,indmax)
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      print 101, c, s
      call PrintVectorD(1,n,x,1,'X ')
      print 102, (indx(i),i=1,n)
      call PrintVectorD(1,indmax,y,1,'Y ')
*
*       Call DROTI subroutine
      call DROTI(n,x,indx,y,c,s)
*
      print*
      print*, '     OUTPUT DATA'
      call PrintVectorD(1,n,x,1,'X ')
      call PrintVectorD(1,indmax,y,1,'Y ')

      stop
 100  format(7x,'N=',i2)
 101  format(7x,'C=',f4.1,'  S=',f4.1)
 102  format(7x,'VECTOR INDX'/10x,10(i2,1x))
 999  stop 1
      end
