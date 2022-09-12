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
!      S G E M M _ C O M P U T E  Example Program Text
!*******************************************************************************

      program   SGEMM_COMPUTE_MAIN
      use ISO_C_BINDING
*
      character*1     transa, transb
      integer         m, n, k, destsize
      integer         lda, ldb, ldc
      real            alpha, beta
      integer         rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc
      parameter       (rmaxa=3, cmaxa=4, rmaxb=4, cmaxb=5,
     $                 rmaxc=3, cmaxc=5)
      parameter       (lda=rmaxa, ldb=rmaxb, ldc=rmaxc)
      real            a(rmaxa,cmaxa), b(rmaxb,cmaxb), c(rmaxc,cmaxc)
      integer*1, allocatable  :: buffer(:)
      integer         allocation_status 

*       External Subroutines
      external        SGEMM_PACK_GET_SIZE
      INTEGER         SGEMM_PACK_GET_SIZE
      external        SGEMM_COMPUTE, SGEMM_PACK
      external        PrintArrayS

*
*      Executable Statements
*
      print*
      print*,'   S G E M M _ C O M P U T E  EXAMPLE PROGRAM'
*
*      Read input data
      read*
      read*, m, n, k
      read*, alpha, beta
      read 100, transa, transb
      if ((transa.eq.'N').or.(transa.eq.'n')) then
        if (m.gt.rmaxa.or.k.gt.cmaxa) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((a(i,j),j=1,k),i=1,m)
      else
        if (k.gt.rmaxa.or.m.gt.cmaxa) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((a(i,j),j=1,m),i=1,k)
      end if
      if ((transb.eq.'N').or.(transb.eq.'n')) then
        if (k.gt.rmaxb.or.n.gt.cmaxb) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((b(i,j),j=1,n),i=1,k)
      else
        if (n.gt.rmaxb.or.k.gt.cmaxb) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((b(i,j),j=1,k),i=1,n)
      end if
      if (m.gt.rmaxc.or.n.gt.cmaxc) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, ((c(i,j),j=1,n),i=1,m)
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 101, m, n, k
      print 102, alpha, beta
      print 103, transa, transb
      if ((transa.eq.'N').or.(transa.eq.'n')) then
        call PrintArrayS(0,0,m,k,a,lda,'A')
      else
        call PrintArrayS(0,0,k,m,a,lda,'A')
      end if
      if ((transb.eq.'N').or.(transb.eq.'n')) then
        call PrintArrayS(0,0,k,n,b,ldb,'B')
      else
        call PrintArrayS(0,0,n,k,b,ldb,'B')
      end if
      call PrintArrayS(0,0,m,n,c,ldc,'C')
*
*     Call SGEMM_PACK_GET_SIZE to allocate buffer for A
      destsize = SGEMM_PACK_GET_SIZE('A',m,n,k)
      allocate(buffer(destsize), stat = allocation_status)

      if(allocation_status.gt.0) then
         write(*,*) 'Cannot allocate buffer '
         stop 1
      endif
*
*     Call SGEMM_PACK subroutine to perform scaling and packing
      call SGEMM_PACK('A',transa,m,n,k,alpha,a,lda,buffer)
*
*     Call SGEMM_COMPUTE subroutine
      call SGEMM_COMPUTE('P',transb,m,n,k,buffer,lda,b,ldb,beta,c,ldc)
*
*     Free allocated storage
      deallocate(buffer)
*
      print*
      print*, '     OUTPUT DATA'
      call PrintArrayS(1,0,m,n,c,ldc,'C')

      stop
 100  format(2(a1,1x))
 101  format(7x,'M=',i1,'  N=',i1,'  K=',i1)
 102  format(7x,'ALPHA=',f4.1,'  BETA=',f4.1)
 103  format(7x,'TRANSA=',a1, '  TRANSB=',a1)
 999  stop 1
      end
