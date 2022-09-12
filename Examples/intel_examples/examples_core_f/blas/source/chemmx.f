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
!      C H E M M  Example Program Text
!*******************************************************************************

      program   CHEMM_MAIN
*
      character*1     side, uplo
      integer         m, n
      integer         lda, ldb, ldc
      complex         alpha, beta
      integer         rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc
      parameter       (rmaxa=5, cmaxa=4, rmaxb=5, cmaxb=5,
     $                rmaxc=5, cmaxc=5)
      parameter       (lda=rmaxa, ldb=rmaxb, ldc=rmaxc)
      complex         a(rmaxa,cmaxa), b(rmaxb,cmaxb), c(rmaxc,cmaxc)
*       External Subroutines
      external        CHEMM, PrintArrayC
*
*      Executable Statements
*
      print*
      print*,'   C H E M M  EXAMPLE PROGRAM'
*
*      Read input data
      read*
      read*, m, n
      read*, alpha, beta
      read 100, side, uplo
      if ((side.eq.'L').or.(side.eq.'l')) then
        if ( (m.gt.rmaxa).or.(m.gt.cmaxa) ) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        if ((uplo.eq.'U').or.(side.eq.'u')) then
          read*, ((a(i,j),j=i,m),i=1,m)
        else
          read*, ((a(i,j),j=1,i),i=1,m)
        end if
      else
        if ( (n.gt.rmaxa).or.(n.gt.cmaxa) ) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        if ((uplo.eq.'U').or.(side.eq.'u')) then
          read*, ((a(i,j),j=i,n),i=1,n)
        else
          read*, ((a(i,j),j=1,i),i=1,n)
        end if
      end if
      if ( (m.gt.rmaxb).or.(n.gt.cmaxb).or.
     $     (m.gt.rmaxc).or.(n.gt.cmaxc) ) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, ((b(i,j),j=1,n),i=1,m)
      read*, ((c(i,j),j=1,n),i=1,m)
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 101, m, n
      print 102, alpha, beta
      print 103, side, uplo
      if ((side.eq.'L').or.(side.eq.'l')) then
        if ((uplo.eq.'U').or.(side.eq.'u')) then
          call PrintArrayC(0,1,m,m,a,lda,'A')
        else
          call PrintArrayC(0,-1,m,m,a,lda,'A')
        end if
      else
        if ((uplo.eq.'U').or.(side.eq.'u')) then
          call PrintArrayC(0,1,n,n,a,lda,'A')
        else
          call PrintArrayC(0,-1,n,n,a,lda,'A')
        end if
      end if
      call PrintArrayC(0,0,m,n,b,ldb,'B')
      call PrintArrayC(0,0,m,n,c,ldc,'C')
*
*      Call CHEMM subroutine
      call CHEMM(side, uplo, m, n, alpha, a, lda, b, ldb,
     $           beta, c, ldc)
*
      print*
      print*, '     OUTPUT DATA'
      call PrintArrayC(1,0,m,n,c,ldc,'C')

      stop
 100  format(2(a1,1x))
 101  format(7x,'M=',i1,'  N=',i1)
 102  format(7x,'ALPHA=(',f4.1,',',f4.1,' )',
     $        '  BETA=(',f4.1,',',f4.1,' )')
 103  format(7x,'SIDE=',a1, '  UPLO=',a1)
 999  stop 1
      end
