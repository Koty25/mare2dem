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
!      D G E M M _ B A T C H  Example Program Text
!*******************************************************************************

      program   DGEMM_BATCH_MAIN
      use f95_precision, only: wp => dp
      use blas95, only: dgemm_batch
      use, intrinsic :: ISO_C_BINDING
      implicit none
      TYPE arrayD
        real(wp),DIMENSION(:,:),ALLOCATABLE :: ARRAY
      ENDTYPE arrayD
      integer i, j, grpcount, status
      integer ig, is, matidx, matrixcount
      TYPE(arrayD), allocatable :: a(:), b(:), c(:)
      integer(KIND=C_SIZE_T), allocatable :: a_array(:), b_array(:), c_array(:)
      integer, allocatable :: grpsize(:)
      integer, allocatable :: m(:), n(:), k(:)
      real(wp), allocatable :: alpha(:), beta(:)
      character(len = 1), allocatable :: transa(:), transb(:)
!     External Subroutines
      external PrintArrayD
!
!      Executable Statements
!
      print*
      print*,'   G E M M _ B A T C H  EXAMPLE PROGRAM'
!
!      Read input data
      read*
      read*, grpcount

      allocate(grpsize(grpcount), STAT = status)
      if ( status .NE. 0 ) then
        print*, ' Insufficient memory'
        goto 999
      end if
      allocate(m(grpcount), n(grpcount), k(grpcount), STAT = status)
      if ( status .NE. 0 ) then
        print*, ' Insufficient memory'
        goto 999
      end if
      allocate(alpha(grpcount), beta(grpcount), STAT = status)
      if ( status .NE. 0 ) then
        print*, ' Insufficient memory'
        goto 999
      end if
      allocate(transa(grpcount), transb(grpcount), STAT = status)
      if ( status .NE. 0 ) then
        print*, ' Insufficient memory'
        goto 999
      end if

      read*, (grpsize(i), i=1,grpcount)
      do ig = 1, grpcount 
        read*, m(ig), n(ig), k(ig)
      end do
      do ig = 1, grpcount 
        read*, alpha(ig), beta(ig)
      end do
      do ig = 1, grpcount 
        read 100, transa(ig), transb(ig)
      end do

      matrixcount = 0
      do ig = 1, grpcount 
        matrixcount = matrixcount + grpsize(ig)
      end do

      allocate(a(matrixcount), b(matrixcount), c(matrixcount), STAT = status)
      if ( status .NE. 0 ) then
        print*, ' Insufficient memory'
        goto 999
      end if
      allocate(a_array(matrixcount), b_array(matrixcount),c_array(matrixcount), STAT = status)
      if ( status .NE. 0 ) then
        print*, ' Insufficient memory'
        goto 999
      end if
     
      matidx = 1
      do ig = 1, grpcount
        do is = 1, grpsize(ig)
          if ((transa(ig).eq.'N').or.(transa(ig).eq.'n')) then
            allocate(a(matidx)%ARRAY(m(ig),k(ig)), STAT = status)
            if ( status .NE. 0 ) then
              print*, ' Insufficient memory'
              goto 999
            end if
            read*, ((a(matidx)%ARRAY(i,j),j=1,k(ig)),i=1,m(ig))
          else
            allocate(a(matidx)%ARRAY(k(ig),m(ig)), STAT = status)
            if ( status .NE. 0 ) then
              print*, ' Insufficient memory'
              goto 999
            end if
            read*, ((a(matidx)%ARRAY(i,j),j=1,m(ig)),i=1,k(ig))
          end if
          if ((transb(ig).eq.'N').or.(transb(ig).eq.'n')) then
            allocate(b(matidx)%ARRAY(k(ig),n(ig)), STAT = status)
            if ( status .NE. 0 ) then
              print*, ' Insufficient memory'
              goto 999
            end if
            read*, ((b(matidx)%ARRAY(i,j),j=1,n(ig)),i=1,k(ig))
          else
            allocate(b(matidx)%ARRAY(n(ig),k(ig)), STAT = status)
            if ( status .NE. 0 ) then
              print*, ' Insufficient memory'
              goto 999
            end if
            read*, ((b(matidx)%ARRAY(i,j),j=1,k(ig)),i=1,n(ig))
          end if
          allocate(c(matidx)%ARRAY(m(ig),n(ig)), STAT = status)
          if ( status .NE. 0 ) then
            print*, ' Insufficient memory'
            goto 999
          end if
          read*, ((c(matidx)%ARRAY(i,j),j=1,n(ig)),i=1,m(ig))

          a_array(matidx) = LOC(a(matidx)%ARRAY(1,1))
          b_array(matidx) = LOC(b(matidx)%ARRAY(1,1))
          c_array(matidx) = LOC(c(matidx)%ARRAY(1,1))

          matidx = matidx + 1
        end do 
      end do
!
!       Print input data

      matidx = 1
      do ig = 1, grpcount
        print*
        print*, '     INPUT DATA'
        print 101, m(ig), n(ig), k(ig)
        print 102, alpha(ig), beta(ig)
        print 103, transa(ig), transb(ig)
        do is = 1, grpsize(ig)
          if ((transa(ig).eq.'N').or.(transa(ig).eq.'n')) then
            call PrintArrayD(0,0,m(ig),k(ig),a(matidx)%ARRAY,m(ig),'A')
          else
            call PrintArrayD(0,0,k(ig),m(ig),a(matidx)%ARRAY,k(ig),'A')
          end if
          if ((transb(ig).eq.'N').or.(transb(ig).eq.'n')) then
            call PrintArrayD(0,0,k(ig),n(ig),b(matidx)%ARRAY,k(ig),'B')
          else
            call PrintArrayD(0,0,n(ig),k(ig),b(matidx)%ARRAY,n(ig),'B')
          end if
          call PrintArrayD(0,0,m(ig),n(ig),c(matidx)%ARRAY,m(ig),'C')
          matidx = matidx + 1
        end do
      end do
!
!      Call DGEMM_BATCH subroutine
      call dgemm_batch(a_array,b_array,c_array,m,n,k,grpsize,transa,transb,alpha,beta)
      print*
      print*, '     OUTPUT DATA'

      matidx = 1
      do ig = 1, grpcount
        do is = 1, grpsize(ig)
          call PrintArrayD(1,0,m(ig),n(ig),c(matidx)%ARRAY,m(ig),'C')
          matidx = matidx + 1
        end do
      end do

!      Deallocate memory
      deallocate(grpsize)
      deallocate(m, n, k)
      deallocate(alpha, beta)
      deallocate(transa, transb)
      deallocate(a_array, b_array,c_array)
      do i = 1, matrixcount
        deallocate(a(i)%ARRAY)
        deallocate(b(i)%ARRAY)
        deallocate(c(i)%ARRAY)
      end do
      deallocate(a, b, c)

 100  format(2(a1,1x))
 101  format(7x,'M=',i1,'  N=',i1,'  K=',i1)
 102  format(7x,'ALPHA=',f5.2,'  BETA=',f5.2)
 103  format(7x,'TRANSA=',a1, '  TRANSB=',a1)
      stop
 999  stop 1
      end
