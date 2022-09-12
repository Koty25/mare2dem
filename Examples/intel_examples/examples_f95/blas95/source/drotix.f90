!===============================================================================
! Copyright 2005-2020 Intel Corporation.
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

      use f95_precision, only: wp => dp
      use blas95, only: roti

      implicit none

      integer :: n
      real(wp) :: c, s

      real(wp), allocatable :: x(:), y(:)
      integer, allocatable :: indx(:)
      integer :: i, indmax
!     External Subroutines
      integer MaxValue
      external MaxValue, PrintVectorD

!     Executable Statements

      print*
      print*,'   D R O T I  EXAMPLE PROGRAM'

!     Read input data from input file
      read*
      read*, n
      read*, c, s
      allocate(x(n))
      allocate(indx(n))
      read*, (x(i),i=1,n)
      read*, (indx(i),i=1,n)
      indmax = MaxValue(n,indx)
      allocate(y(indmax))
      read*, (y(i),i=1,indmax)

!     Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      print 101, c, s
      call PrintVectorD(1,n,x,1,'X ')
      print 102, (indx(i),i=1,n)
      call PrintVectorD(1,indmax,y,1,'Y ')

!     Call DROTI subroutine
      call ROTI(x, indx, y, c, s)

      print*
      print*, '     OUTPUT DATA'
      call PrintVectorD(1,n,x,1,'X ')
      call PrintVectorD(1,indmax,y,1,'Y ')

      deallocate(x)
      deallocate(indx)
      deallocate(y)

 100  format(7x,'N=',i2)
 101  format(7x,'C=',f4.1,'  S=',f4.1)
 102  format(7x,'VECTOR INDX'/10x,10(i2,1x))
      stop
      end
