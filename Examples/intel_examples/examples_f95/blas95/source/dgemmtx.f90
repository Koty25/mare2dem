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
!      D G E M M T  Example Program Text
!*******************************************************************************

      program   DGEMMT_MAIN

      use f95_precision, only: wp => dp
      use blas95, only: gemmt

      implicit none

      character(len = 1) :: uplo, transa, transb
      integer :: n, k, lda, ldb, ldc
      real(wp) :: alpha, beta
      real(wp), allocatable :: a(:,:), b(:,:), c(:,:)
      integer :: i, j

!     External Subroutines
      external PrintArrayD

!     Executable Statements

      print*
      print*,'   D G E M M T  EXAMPLE PROGRAM'

!     Read input data
      read*
      read*, n, k
      read*, alpha, beta
      read 100, uplo, transa, transb

      if ((transa.eq.'N').or.(transa.eq.'n')) then
        lda = n
        allocate(a(n, k))
        read*, ((a(i,j),j=1,k),i=1,n)
      else
        lda = k
        allocate(a(k, n))
        read*, ((a(i,j),j=1,n),i=1,k)
      end if
      if ((transb.eq.'N').or.(transb.eq.'n')) then
        ldb = k
        allocate(b(k, n))
        read*, ((b(i,j),j=1,n),i=1,k)
      else
        ldb = n
        allocate(b(n, k))
        read*, ((b(i,j),j=1,k),i=1,n)
      end if
      ldc = n
      allocate(c(n, n))
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        read*, ((c(i,j),j=i,n),i=1,n)
      else
        read*, ((c(i,j),j=1,i),i=1,n)
      end if

!     Print input data
      print*
      print*, '     INPUT DATA'
      print 101, n, k
      print 102, alpha, beta
      print 103, uplo, transa, transb
      if ((transa.eq.'N').or.(transa.eq.'n')) then
        call PrintArrayD(0,0,n,k,a,lda,'A')
      else
        call PrintArrayD(0,0,k,n,a,lda,'A')
      end if
      if ((transb.eq.'N').or.(transb.eq.'n')) then
        call PrintArrayD(0,0,k,n,b,ldb,'B')
      else
        call PrintArrayD(0,0,n,k,b,ldb,'B')
      end if
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        call PrintArrayD(0,1,n,n,c,ldc,'C')
      else
        call PrintArrayD(0,-1,n,n,c,ldc,'C')
      end if

!     Call DGEMMT subroutine
      call GEMMT(a, b, c, uplo, transa, transb, alpha, beta)

      print*
      print*, '     OUTPUT DATA'
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        call PrintArrayD(1,1,n,n,c,ldc,'C')
      else
        call PrintArrayD(1,-1,n,n,c,ldc,'C')
      end if

      deallocate(a)
      deallocate(b)
      deallocate(c)

 100  format(3(a1,1x))
 101  format(7x,'N=',i1,'  K=',i1)
 102  format(7x,'ALPHA=',f3.1,'  BETA=',f3.1)
 103  format(7x,'UPLO=',a1, '  TRANSA=',a1, '  TRANSB=',a1)
      stop
      end
