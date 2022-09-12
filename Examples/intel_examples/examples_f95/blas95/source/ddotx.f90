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
!      D D O T  Example Program Text
!*******************************************************************************

      program  DDOT_MAIN
      use f95_precision, only: wp => dp
      use blas95, only: dot

      implicit none

      integer :: n, incx, incy
      integer :: nx, ny, nx1, nx2, ny1, ny2
      real(wp), allocatable :: sx(:), sy(:)
      real(wp) :: res
      integer :: i

!     Intrinsic Functions
      intrinsic abs
!     External Subroutines
      external PrintVectorD

!     Executable Statements

      print*
      print*,'   D D O T  EXAMPLE PROGRAM'

!     Read input data from input file
      read*
      read*, n
      read*, incx, incy

      nx = 1+(n-1)*abs(incx)
      ny = 1+(n-1)*abs(incy)

      allocate(sx(nx))
      allocate(sy(ny))
      read*, (sx(i),i=1,nx)
      read*, (sy(i),i=1,ny)

!     Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      call PrintVectorD(0,n,sx,incx,'X ')
      call PrintVectorD(0,n,sy,incy,'Y ')

!     Call DDOT subroutine
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
      res = DOT(sx(nx1:nx2:incx),sy(ny1:ny2:incy))

      print*
      print*, '     OUTPUT DATA'
      print 101, res

      deallocate(sx)
      deallocate(sy)

 100  format(7x,'N=',i2)
 101  format(7x,'DDOT = ',f8.3)
      stop
      end
