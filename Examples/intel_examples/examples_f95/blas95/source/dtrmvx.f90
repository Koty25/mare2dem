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
!      D T R M V  Example Program Text
!*******************************************************************************

      program   DTRMV_MAIN
      use f95_precision, only: wp => dp
      use blas95, only: trmv

      implicit none

      character(len = 1) :: uplo, trans, diag
      integer :: n, lda, incx
      integer :: nx, nx1, nx2
      real(wp), allocatable :: a(:,:), x(:)
      integer :: i, j

!     Intrinsic Functions
      intrinsic abs
!     External Subroutines
      external PrintVectorD, PrintArrayD

!     Executable Statements

      print*
      print*, '   D T R M V  EXAMPLE PROGRAM'

!     Read input data from input file
      read*
      read*, n
      read*, incx
      read 100, uplo, trans, diag
      nx = 1+(n-1)*abs(incx)
      lda = n
      allocate(x(nx))
      allocate(a(n, n))
      read*, (x(i),i=1, nx)
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        read*, ((a(i,j),j=i,n),i=1,n)
      else
        read*, ((a(i,j),j=1,i),i=1,n)
      end if

!     Print input data
      print*
      print*, '     INPUT DATA'
      print 101, n
      print 102, uplo, trans, diag
      call PrintVectorD(0,n,x,incx,'X ')
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        call PrintArrayD(0,1,n,n,a,lda,'A')
      else
        call PrintArrayD(0,-1,n,n,a,lda,'A')
      end if

!     Call DTRMV subroutine
      if (incx > 0) then
          nx1 = 1
          nx2 = nx
      else
          nx1 = nx
          nx2 = 1
      end if
      call TRMV(a, x(nx1:nx2:incx), uplo, trans, diag)

      print*
      print*, '     OUTPUT DATA'
      call PrintVectorD(1,n,x,incx,'X ')

      deallocate(x)
      deallocate(a)

 100  format(3(a1,1x))
 101  format(7x,'N=',i1)
 102  format(7x,'UPLO=',a1,'  TRANS=',a1,'  DIAG=',a1)
      stop
      end
