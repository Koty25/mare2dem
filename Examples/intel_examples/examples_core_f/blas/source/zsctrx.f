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
!      Z S C T R  Example Program Text
!*******************************************************************************

      program   ZSCTR_MAIN
*
      integer             nmax, ymax
      parameter           (nmax=10, ymax=10)
      integer             n
      integer             indx(nmax)
      complex*16          x(nmax), y(ymax)

      integer             i, indmax

*       External Function
      external            ZSCTR, MaxValue, PrintVectorZ
*
*       Executable Statements ..
*
      print*
      print*,'   Z S C T R  EXAMPLE PROGRAM'
*       Read input data from input file
      read*
      read*, n
      if (n.gt.nmax) then
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
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      call PrintVectorZ(1,n,x,1,'X ')
      print 101, (indx(i),i=1,n)
*
*       Call ZSCTR function
      call ZSCTR(n,x,indx,y)
*
      print*
      print*, '     OUTPUT DATA'
      call PrintVectorZ(1,indmax,y,1,'Y ')

      stop
 100  format(7x,'N=',i2)
 101  format(7x,'VECTOR INDX'/10x,10(i2,1x))
 999  stop 1
      end
