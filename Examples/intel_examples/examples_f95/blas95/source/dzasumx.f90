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
!      D Z A S U M  Example Program Text
!*******************************************************************************

      program  DZASUM_MAIN

      use f95_precision, only: wp => dp
      use blas95, only: asum

      implicit none

      integer :: n, incx
      integer :: nx, nx1, nx2

      complex(wp), allocatable :: x(:)
      real(wp) :: sum
      integer  :: i
!     Intrinsic Functions
      intrinsic abs
!     External Subroutines
      external PrintVectorZ

!     Executable Statements

      print*
      print*,'   D Z A S U M  EXAMPLE PROGRAM'
!     Read input data from input file
      read*
      read*, n, incx

      nx = 1+(n-1)*abs(incx)

      allocate(x(nx))
      read*, (x(i),i=1,nx)

!     Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      call PrintVectorZ(0,n,x,incx,'X ')

!     Call DZASUM subroutine
      if (incx > 0) then
          nx1 = 1
          nx2 = nx
      else
          nx1 = nx
          nx2 = 1
      end if
      sum = ASUM(x(nx1:nx2:incx))

      print*
      print*, '     OUTPUT DATA'
      print 101, sum

      deallocate(x)

 100  format(7x,'N=',i2)
 101  format(10x,'DZASUM = ',f7.3)
      stop
      end
