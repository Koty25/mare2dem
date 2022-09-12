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
!      D G B M V  Example Program Text
!*******************************************************************************

      program   DGBMV_MAIN
      use f95_precision, only: wp => dp
      use blas95, only: gbmv

      implicit none

      character(len = 1) :: trans
      integer :: m, n, k, kl, ku, incx, incy
      integer :: nx, ny, nx1, nx2, ny1, ny2
      integer :: lda
      real(wp) :: alpha, beta
      real(wp), allocatable :: a(:,:), x(:), y(:)
      integer :: i, j, ku1, kl1

!     Intrinsic Functions
      intrinsic abs
!     External Subroutines
      external PrintVectorD, PrintBandArrayD

!     Executable Statements

      print*
      print*, '   D G B M V  EXAMPLE PROGRAM'

!     Read input data from input file
      read*
      read*, m, n, kl, ku
      read*, incx, incy
      read*, alpha, beta
      read 100, trans

      if ((trans.eq.'N').or.(trans.eq.'n')) then
        nx = 1+(n-1)*abs(incx)
        ny = 1+(m-1)*abs(incy)
      else
        nx = 1+(m-1)*abs(incx)
        ny = 1+(n-1)*abs(incy)
      end if
      allocate(x(nx))
      allocate(y(ny))
      read*,(x(i),i=1, nx)
      read*,(y(i),i=1, ny)
      lda = kl+ku+1
      allocate(a(lda, n))
      a = 0
      if (ku.ge.n) then
         ku1 = n-1
      else
         ku1 = ku
      end if
      if (kl.ge.m) then
         kl1 = m-1
      else
         kl1 = kl
      end if
      do i=1, ku1+1
         if ((ku1-i+m+1).ge.n) then
            k = n
         else
            k = ku1-i+m+1
         end if
         read*, (a(ku-ku1+i,j),j=ku1-i+2,k)
      end do
      do i=ku1+2, ku1+kl1+1
         if ((m+ku1-i+1).ge.n) then
            k = n
         else
            k = m+ku1-i+1
         end if
         read*, (a(ku-ku1+i,j),j=1,k)
      end do

!     Print input data
      print*
      print*, '     INPUT DATA'
      print 101, m, n
      print 102, alpha, beta
      print 103, trans
      if ((trans.eq.'N').or.(trans.eq.'n')) then
        call PrintVectorD(0,n,x,incx,'X ')
        call PrintVectorD(0,m,y,incy,'Y ')
      else
        call PrintVectorD(0,m,x,incx,'X ')
        call PrintVectorD(0,n,y,incy,'Y ')
      end if
      call PrintBandArrayD(0,kl,ku,m,n,a,lda,'A')

!     Call DGBMV subroutine
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
      call GBMV(a, x(nx1:nx2:incx), y(ny1:ny2:incy), kl, m, alpha, beta, trans)

      print*
      print*, '     OUTPUT DATA'
      if ((trans.eq.'N').or.(trans.eq.'n')) then
        call PrintVectorD(0,m,y,incy,'Y ')
      else
        call PrintVectorD(0,n,y,incy,'Y ')
      end if

      deallocate(x)
      deallocate(y)
      deallocate(a)

 100  format(a1)
 101  format(7x,'M=',i1,'  N=',i1)
 102  format(7x,'ALPHA=',f4.1,'  BETA=',f4.1)
 103  format(7x,'TRANS=',a1)
      stop
      end
