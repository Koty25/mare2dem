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
!      G E M M _ S 8 U 8 S 3 2 C O M P U T E  Example Program Text
!*******************************************************************************

      program   GEMM_S8U8S32_COMPUTE_MAIN

*
      character*1     transa, transb, offsetc
      integer         m, n, k, destsize, sizeco
      integer         lda, ldb, ldc
      real            alpha, beta
      integer         irmaxa, icmaxa, irmaxb, icmaxb, irmaxc, icmaxc
      parameter       (irmaxa=3, icmaxa=4, irmaxb=4, icmaxb=5,
     $                 irmaxc=3, icmaxc=5, icmaxdim=5)
      parameter       (lda=irmaxa, ldb=irmaxb, ldc=irmaxc)
      integer*1       a(irmaxa,icmaxa), b(irmaxb,icmaxb), ao, bo
      integer*4       c(irmaxc,icmaxc)
      integer*4       co(icmaxdim)
      integer*1, allocatable  :: buffer(:)
      integer         allocation_status 

*       External Subroutines
      external        GEMM_S8U8S32_COMPUTE, GEMM_S8U8S32_PACK
      external        PrintArrayI8, PrintArrayI32, PrintVectorI32
      external        GEMM_S8U8S32_PACK_GET_SIZE
      INTEGER         GEMM_S8U8S32_PACK_GET_SIZE

*
*      Executable Statements
*
      print*
      print*,'G E M M _ S 8 U 8 S 3 2 _ C O M P U T E  EXAMPLE PROGRAM'
*
*      Read input data
      read*
      read*, m, n, k
      read*, alpha, beta
      read 100, transa, transb, offsetc
      read*, ao, bo

      destsize = 0
      sizeco = 0

      if ((offsetc.eq.'F').or.(offsetc.eq.'f')) then
        if (1.gt.icmaxdim) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        sizeco = 1
      else if ((offsetc.eq.'R').or.(offsetc.eq.'r')) then
            if (n.gt.icmaxdim) then
               print*, ' Insufficient memory for arrays'
               goto 999
            end if
            sizeco = n
      else if ((offsetc.eq.'C').or.(offsetc.eq.'c')) then
               if (m.gt.icmaxdim) then
                  print*, ' Insufficient memory for arrays'
                  goto 999
               end if
               sizeco = m
      else
         print*, 'Offset values not valid'
         goto 999
            end if

      read*, (co(i),i=1,sizeco)

      if ((transa.eq.'N').or.(transa.eq.'n')) then
        if (m.gt.irmaxa.or.k.gt.icmaxa) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((a(i,j),j=1,k),i=1,m)
      else
        if (k.gt.irmaxa.or.m.gt.icmaxa) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((a(i,j),j=1,m),i=1,k)
      end if

      if ((transb.eq.'N').or.(transb.eq.'n')) then
        if (k.gt.irmaxb.or.n.gt.icmaxb) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((b(i,j),j=1,n),i=1,k)
      else
        if (n.gt.irmaxb.or.k.gt.icmaxb) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*, ((b(i,j),j=1,k),i=1,n)
      end if
      if (m.gt.irmaxc.or.n.gt.icmaxc) then
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
      print 103, transa, transb, offsetc
      print 104, ao, bo

      if ((offsetc.eq.'F').or.(offsetc.eq.'f')) then
        call PrintVectorI32(1,co,"co")
      else
         if ((offsetc.eq.'R').or.(offsetc.eq.'r')) then
            call PrintVectorI32(n,co,"co")
         else
            if ((offsetc.eq.'C').or.(offsetc.eq.'c')) then
               call PrintVectorI32(m,co,"co")
            end if
         end if
      end if

      if ((transa.eq.'N').or.(transa.eq.'n')) then
        call PrintArrayI8(0,0,m,k,a,lda,'A')
      else
        call PrintArrayI8(0,0,k,m,a,lda,'A')
      end if

      if ((transb.eq.'N').or.(transb.eq.'n')) then
        call PrintArrayI8(0,0,k,n,b,ldb,'B')
      else
        call PrintArrayI8(0,0,n,k,b,ldb,'B')
      end if
      call PrintArrayI32(0,0,m,n,c,ldc,'C')
*
*     Call GEMM_S8U8S32_PACK_GET_SIZE function to allocate buffer for A

      destsize = GEMM_S8U8S32_PACK_GET_SIZE('A',m,n,k)
      allocate(buffer(destsize), stat = allocation_status)

      if (allocation_status.gt.0) then
         write(*,*) 'Cannot allocate buffer '
         stop 1
      endif
*
*     Call GEMM_S8U8S32_PACK subroutine to perform packing
      call GEMM_S8U8S32_PACK('A',transa,m,n,k,a,lda,buffer)
*
*     call GEMM_S8U8S32_COMPUTE subroutine
      call GEMM_S8U8S32_COMPUTE('P',transb,offsetc,m,n,k,alpha,
     +buffer,lda,ao,b,ldb,bo,beta,c,ldc,co)
*
*     Free allocated storage
      deallocate(buffer)
*
      print*
      print*, '     OUTPUT DATA'
      call PrintArrayI32(1,0,m,n,c,ldc,'C')

      stop
 100  format(3(a1,1x))
 101  format(7x,'M=',i1,'  N=',i1,'  K=',i1)
 102  format(7x,'ALPHA=',f4.1,'  BETA=',f4.1)
 103  format(7x,'TRANSA=',a1, '  TRANSB=',a1, ' OFFSETC=',a1)
 104  format(7x,'ao=',i1,'  bo=',i1)
 999  stop 1
      end
