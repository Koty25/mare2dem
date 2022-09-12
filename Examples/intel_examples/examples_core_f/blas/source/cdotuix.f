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
!      C D O T U I  Example Program Text
!*******************************************************************************

      program   CDOTUI_MAIN
*
      integer             nmax, ymax
      parameter           (nmax=10, ymax = 15)
      integer             n
      integer             indx(nmax)
      complex             x(nmax), y(ymax)

      integer             i, indmax
      complex             res

*       External Function
      external            CDOTUI, MaxValue, PrintVectorC
      complex             CDOTUI
*
*       Executable Statements
*
      print*
      print*,'   C D O T U I  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, n
      if (n.gt.nmax) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, (x(i),i=1,n)
      read*, (indx(i),i=1,n)
      indmax = MaxValue(n, indx)
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
      call PrintVectorC(1,n,x,1,'X ')
      print 101, (indx(i),i=1,n)
      call PrintVectorC(1,indmax,y,1,'Y ')
*
*       Call CDOTUI function
      res = CDOTUI(n,x,indx,y)
*
      print*
      print*, '     OUTPUT DATA'
      print 102, res

      stop
 100  format(7x,'N=',i2)
 101  format(7x,'VECTOR INDX'/10x,10(i2,1x))
 102  format(11x,'CDOTUI = (',f6.2,',',f6.2,' )')
 999  stop 1
      end