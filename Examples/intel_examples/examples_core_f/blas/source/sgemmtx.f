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
!      S G E M M T  Example Program Text
!*******************************************************************************

      program   SGEMMT_MAIN
*
      character*1     uplo, transa, transb
      integer         n, k
      integer         lda, ldb, ldc
      real            alpha, beta
      integer         rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc
      parameter       (rmaxa=3, cmaxa=4, rmaxb=5, cmaxb=5,
     $                 rmaxc=3, cmaxc=5)
      parameter       (lda=rmaxa, ldb=rmaxb, ldc=rmaxc)
      real            a(rmaxa,cmaxa), b(rmaxb,cmaxb), c(rmaxc,cmaxc)
*       External Subroutines
      external        SGEMMT, PrintArrayS
*
*      Executable Statements
*
      print*
      print*,'   S G E M M T  EXAMPLE PROGRAM'
*
*      Read input data
      read*
      read*, n, k
      read*, alpha, beta
      read 100, uplo, transa, transb
      if ((transa.eq.'N').or.(transa.eq.'n')) then
        if ( n.gt.rmaxa.or.k.gt.cmaxa ) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((a(i,j),j=1,k),i=1,n)
      else
        if ( k.gt.rmaxa.or.n.gt.cmaxa ) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((a(i,j),j=1,n),i=1,k)
      end if
      if ((transb.eq.'N').or.(transb.eq.'n')) then
        if ( k.gt.rmaxb.or.n.gt.cmaxb ) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((b(i,j),j=1,n),i=1,k)
      else
        if ( n.gt.rmaxb.or.k.gt.cmaxb ) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((b(i,j),j=1,k),i=1,n)
      end if
      if ( n.gt.rmaxc.or.n.gt.cmaxc ) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        read*, ((c(i,j),j=i,n),i=1,n)
      else
        read*, ((c(i,j),j=1,i),i=1,n)
      end if
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 101, n, k
      print 102, alpha, beta
      print 103, uplo, transa, transb
      if ((transa.eq.'N').or.(transa.eq.'n')) then
        call PrintArrayS(0,0,n,k,a,lda,'A')
      else
        call PrintArrayS(0,0,k,n,a,lda,'A')
      end if
      if ((transb.eq.'N').or.(transb.eq.'n')) then
        call PrintArrayS(0,0,k,n,b,ldb,'B')
      else
        call PrintArrayS(0,0,n,k,b,ldb,'B')
      end if
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        call PrintArrayS(0,1,n,n,c,ldc,'C')
      else
        call PrintArrayS(0,-1,n,n,c,ldc,'C')
      end if
*
*      Call SGEMMT subroutine
      call SGEMMT(uplo, transa, transb, n, k, alpha, a, lda, b, ldb,
     $            beta, c, ldc)
*
      print*
      print*, '     OUTPUT DATA'
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        call PrintArrayS(1,1,n,n,c,ldc,'C')
      else
        call PrintArrayS(1,-1,n,n,c,ldc,'C')
      end if

      stop
 100  format(3(a1,1x))
 101  format(7x,'N=',i1,'  K=',i1)
 102  format(7x,'ALPHA=',f3.1,'  BETA=',f3.1)
 103  format(7x,'UPLO=',a1, '  TRANSA=',a1, '  TRANSB=',a1)
 999  stop 1
      end
