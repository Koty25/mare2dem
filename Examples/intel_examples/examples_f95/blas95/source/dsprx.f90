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
!      D S P R  Example Program Text
!*******************************************************************************

      program   DSPR_MAIN

      use f95_precision, only: wp => dp
      use blas95, only: spr

      implicit none

      character(len = 1) :: uplo
      integer :: n, incx
      integer :: nx, nap, nx1, nx2
      real(wp) :: alpha
      real(wp), allocatable :: ap(:), x(:)
      integer :: i

!     Intrinsic Functions
      intrinsic abs
!     External Subroutines
      external PrintVectorD

!     Executable Statements

      print*
      print*, '   D S P R  EXAMPLE PROGRAM'

!     Read input data from input file
      read*
      read*, n
      read*, incx
      read*, alpha
      read 100, uplo

      nx = 1+(n-1)*abs(incx)
      nap = (n*(n+1))/2

      allocate(x(nx))
      allocate(ap(nap))
      read*, (x(i),i=1, nx)
      read*, (ap(i),i=1,nap)

!     Print input data
      print*
      print*, '     INPUT DATA'
      print 101, n
      print 102, alpha
      print 103, uplo
      call PrintVectorD(0,n,x,incx,'X ')
      call PrintVectorD(1,nap,ap,1,'AP')

!     Call DSPR subroutine
      if (incx > 0) then
          nx1 = 1
          nx2 = nx
      else
          nx1 = nx
          nx2 = 1
      end if
      call SPR(ap, x(nx1:nx2:incx), uplo, alpha)

      print*
      print*, '     OUTPUT DATA'
      call PrintVectorD(1,nap,ap,1,'AP')

      deallocate(x)
      deallocate(ap)

 100  format(a1)
 101  format(7x,'N=',i1)
 102  format(7x,'ALPHA=',f4.1)
 103  format(7x,'UPLO=',a1)
      stop
      end
