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
!      D N R M 2  Example Program Text
!*******************************************************************************

      program  DNRM2_MAIN

      use f95_precision, only: wp => dp
      use blas95, only: nrm2

      implicit none

      integer :: n, incx
      integer :: nx, nx1, nx2
      real(wp), allocatable :: x(:)
      real(wp) :: res
      integer :: i
!     Intrinsic Functions
      intrinsic abs
!     External Subroutines
      external PrintVectorD

!     Executable Statements

      print*
      print*,'   D N R M 2  EXAMPLE PROGRAM'

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
      call PrintVectorD(0,n,x,incx,'X ')

!     Call DNRM2 subroutine
      if (incx > 0) then
          nx1 = 1
          nx2 = nx
      else
          nx1 = nx
          nx2 = 1
      end if
      res = NRM2(x(nx1:nx2:incx))

      print*
      print*, '     OUTPUT DATA'
      print 101, res

      deallocate(x)

 100  format(7x,'N=',i2)
 101  format(10x,'DNRM2 = ',f8.3)
      stop
      end
