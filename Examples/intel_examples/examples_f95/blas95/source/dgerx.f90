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
!      D G E R  Example Program Text
!*******************************************************************************

      program   DGER_MAIN
      use f95_precision, only: wp => dp
      use blas95, only: ger

      implicit none

      integer :: m, n, lda, incx, incy
      integer :: nx, ny, nx1, nx2, ny1, ny2
      real(wp) :: alpha
      real(wp), allocatable :: a(:,:), x(:), y(:)
      integer :: i, j

!     Intrinsic Functions
      intrinsic abs
!     External Subroutines
      external PrintVectorD, PrintArrayD

!     Executable Statements

      print*
      print*, '   D G E R  EXAMPLE PROGRAM'

!     Read input data from input file
      read*
      read*, m, n
      read*, incx, incy
      read*, alpha
      nx = 1+(m-1)*abs(incx)
      ny = 1+(n-1)*abs(incy)
      allocate(x(nx))
      allocate(y(ny))
      read*, (x(i),i=1, nx)
      read*, (y(i),i=1, ny)
      lda = m
      allocate(a(m, n))
      read*, ((a(i,j),j=1,n),i=1,m)

!     Print input data
      print*
      print*, '     INPUT DATA'
      print 101, m, n
      print 102, alpha
      call PrintVectorD(0,m,x,incx,'X ')
      call PrintVectorD(0,n,y,incy,'Y ')
      call PrintArrayD(0,0,m,n,a,lda,'A')

!     Call DGER subroutine
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
      call GER(a, x(nx1:nx2:incx), y(ny1:ny2:incy), alpha)

      print*
      print*, '     OUTPUT DATA'
      call PrintArrayD(0,0,m,n,a,lda,'A')

      deallocate(x)
      deallocate(y)
      deallocate(a)

 101  format(7x,'M=',i1,'  N=',i1)
 102  format(7x,'ALPHA=',f4.1)
      stop
      end
