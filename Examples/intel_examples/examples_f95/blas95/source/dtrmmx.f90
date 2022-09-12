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
!      D T R M M  Example Program Text
!*******************************************************************************

      program   DTRMM_MAIN
      use f95_precision, only: wp => dp
      use blas95, only: trmm

      implicit none

      character(len = 1) :: side, uplo, transa, diag
      integer :: m, n, lda, ldb
      real(wp) :: alpha
      real(wp), allocatable :: a(:,:), b(:,:)
      integer :: i, j

!     External Subroutines
      external PrintArrayD

!     Executable Statements

      print*
      print*,'   D T R M M  EXAMPLE PROGRAM'
!     Read input data from input file
      read*
      read*, m, n
      read*, alpha
      read 100, side, uplo, transa, diag
      if ((side.eq.'L').or.(side.eq.'l')) then
        lda = m
        allocate(a(m, m))
        if ((uplo.eq.'U').or.(uplo.eq.'u')) then
          read*, ((a(i,j),j=i,m),i=1,m)
        else
          read*, ((a(i,j),j=1,i),i=1,m)
        end if
      else
        lda = n
        allocate(a(n, n))
        if ((uplo.eq.'U').or.(uplo.eq.'u')) then
          read*, ((a(i,j),j=i,n),i=1,n)
        else
          read*, ((a(i,j),j=1,i),i=1,n)
        end if
      end if
      ldb = m
      allocate(b(m, n))
      read*, ((b(i,j),j=1,n),i=1,m)

!     Print input data
      print*
      print*, '     INPUT DATA'
      print 101, m, n
      print 102, alpha
      print 103, side, uplo, transa, diag
      if ((side.eq.'L').or.(side.eq.'l')) then
        if ((uplo.eq.'U').or.(uplo.eq.'u')) then
          call PrintArrayD(0, 1, m, m, a, lda, 'A')
        else
          call PrintArrayD(0, -1, m, m, a, lda, 'A')
        end if
      else
        if ((uplo.eq.'U').or.(uplo.eq.'u')) then
          call PrintArrayD(0, 1, n, n, a, lda, 'A')
        else
          call PrintArrayD(0, -1, n, n, a, lda, 'A')
        end if
      end if
      call PrintArrayD(0, 0, m, n, b, ldb, 'B')

!     Call DTRMM subroutine
      call TRMM(a, b, side, uplo, transa, diag, alpha)

!     Print output data
      print*
      print*, '     OUTPUT DATA'
      call PrintArrayD(1, 0, m, n, b, ldb, 'B')

      deallocate(a)
      deallocate(b)

 100  format(4(a1,1x))
 101  format(7x,'M=',i1,'  N=',i1)
 102  format(7x,'ALPHA=',f3.1)
 103  format(7x,'SIDE=',a1, '  UPLO=',a1,'  TRANSA=',a1,'  DIAG=',a1)
      stop
      end
