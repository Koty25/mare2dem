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
!      D R O T M   Example Program Text
!*******************************************************************************

      program   DROTM_MAIN

      use f95_precision, only: wp => dp
      use blas95, only: rotm

      implicit none

      integer :: n, incx, incy
      integer :: nx, ny, nx1, nx2, ny1, ny2
      real(wp) :: param(5)
      real(wp), allocatable :: x(:), y(:)
      integer :: i
!     Intrinsic Function
      intrinsic abs
!     External Subroutines
      external PrintVectorD

!     Executable Statements

      print*
      print*,'   D R O T M  EXAMPLE PROGRAM'

!     Read input data from input file
      read*
      read*, n
      read*, incx, incy

      nx = 1+(n-1)*abs(incx)
      ny = 1+(n-1)*abs(incy)

      read*, (param(i),i=1,5)
      allocate(x(nx))
      allocate(y(ny))
      read*, (x(i),i=1,nx)
      read*, (y(i),i=1,ny)

!     Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      print 101, (param(i),i=1,5)
      call PrintVectorD(0, n, x, incx, 'X ')
      call PrintVectorD(0, n, y, incy, 'Y ')

!     Call DROTM subroutine
      if (incx > 0) then
          nx1 = 1
          nx2 = nx
      else
          nx1 = nx
          nx2 = 1
      end if
      if (incy > 0) then
          ny1 = 1
          ny2 = ny
      else
          ny1 = ny
          ny2 = 1
      end if
      call rotm(x(nx1:nx2:incx), y(ny1:ny2:incy), param)

      print*
      print*, '     OUTPUT DATA'
      call PrintVectorD(1, n, x, incx, 'X ')
      call PrintVectorD(1, n, y, incy, 'Y ')

      deallocate(x)
      deallocate(y)

 100  format(7x,'N=',i2)
 101  format(7x,'VECTOR PARAM',/9x,5(f6.3,1x))
      stop
      end
