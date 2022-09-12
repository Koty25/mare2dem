!===============================================================================
! Copyright 1999-2020 Intel Corporation.
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
!      S T R S M  Example Program Text
!*******************************************************************************

      program   STRSM_MAIN
*
      character*1        side, uplo, transa, diag
      integer            m, n, lda, ldb
      real               alpha
      integer            rmaxa, cmaxa, rmaxb, cmaxb
      parameter          (rmaxa=5, cmaxa=5, rmaxb=5, cmaxb=5)
      parameter          (lda=rmaxa, ldb=rmaxb)
      real               a(rmaxa, cmaxa), b(rmaxb, cmaxb)
      integer            i, j
*      External Subroutines
      external           STRSM, PrintArrayS
*
*      Executable Statements
*
      print*
      print*,'   S T R S M  EXAMPLE PROGRAM'
*       Read input data from input file
      read*
      read*, m, n
      read*, alpha
      read 100, side, uplo, transa, diag
      if ((side.eq.'L').or.(side.eq.'l')) then
        if ( m.gt.rmaxa.or.m.gt.cmaxa.or.
     $       m.gt.rmaxb.or.n.gt.cmaxb ) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        if ((uplo.eq.'U').or.(uplo.eq.'u')) then
          read*, ((a(i,j),j=i,m),i=1,m)
        else
          read*, ((a(i,j),j=1,i),i=1,m)
        end if
      else
        if ( n.gt.rmaxa.or.n.gt.cmaxa.or.
     $       m.gt.rmaxb.or.n.gt.cmaxb ) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        if ((uplo.eq.'U').or.(uplo.eq.'u')) then
          read*, ((a(i,j),j=i,n),i=1,n)
        else
          read*, ((a(i,j),j=1,i),i=1,n)
        end if
      end if
      read*, ((b(i,j),j=1,n),i=1,m)
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 101, m, n
      print 102, alpha
      print 103, side, uplo, transa, diag
      if ((side.eq.'L').or.(side.eq.'l')) then
        if ((uplo.eq.'U').or.(uplo.eq.'u')) then
          call PrintArrayS(0, 1, m, m, a, lda, 'A')
        else
          call PrintArrayS(0, -1, m, m, a, lda, 'A')
        end if
      else
        if ((uplo.eq.'U').or.(uplo.eq.'u')) then
          call PrintArrayS(0, 1, n, n, a, lda, 'A')
        else
          call PrintArrayS(0, -1, n, n, a, lda, 'A')
        end if
      end if
      call PrintArrayS(0, 0, m, n, b, ldb, 'B')
*
*      Call STRSM subroutine
      call STRSM(side, uplo, transa, diag, m, n, alpha, a, lda,
     $           b, ldb)
*
*       Print output data
      print*
      print*, '     OUTPUT DATA'
      call PrintArrayS(1, 0, m, n, b, ldb, 'B')

      stop
 100  format(4(a1,1x))
 101  format(7x,'M=',i1,'  N=',i1)
 102  format(7x,'ALPHA=',f3.1)
 103  format(7x,'SIDE=',a1, '  UPLO=',a1,'  TRANSA=',a1,'  DIAG=',a1)
 999  stop 1
      end
