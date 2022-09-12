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
!      C G E M M 3 M _ B A T C H  Example Program Text
!*******************************************************************************

      program   CGEMM3M_BATCH_MAIN
*
      use, intrinsic :: ISO_C_BINDING
      integer          maxgc, maxmat
      parameter        (maxgc=5, maxmat=20)
      integer          rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc
      parameter        (rmaxa=5, cmaxa=5, rmaxb=5, cmaxb=5,
     $                 rmaxc=5, cmaxc=5)
      integer          lda(maxgc), ldb(maxgc), ldc(maxgc)
      integer          grpcount, ig, is, matidx
      integer          m(maxgc), n(maxgc), k(maxgc)
      integer          grpsize(maxgc)
      complex          a1, b1
      complex          alpha(maxgc), beta(maxgc)
      complex          a(rmaxa,cmaxa,maxmat)
      complex          b(rmaxb,cmaxb,maxmat)
      complex          c(rmaxc,cmaxc,maxmat)
      character*1      transa(maxgc), transb(maxgc)
      integer(KIND=C_SIZE_T) a_array(maxmat), b_array(maxmat), 
     $                       c_array(maxmat)
 
*       External Subroutines
      external         CGEMM3M_BATCH, PrintArrayC
*
*      Executable Statements
*
      print*
      print*,'   C G E M M 3 M _ B A T C H  EXAMPLE PROGRAM'
*
*      Read input data
      read*
      read*, grpcount

      if (grpcount.gt.maxgc) then
        print*, ' group count is larger than the limit'
        goto 999
      end if

      read*, (grpsize(i), i=1,grpcount)

      do ig = 1, grpcount 
        read*, m(ig), n(ig), k(ig)
      end do

      do ig = 1, grpcount 
        read*, a1, b1
        alpha(ig) = a1
        beta(ig) = b1
      end do

      do ig = 1, grpcount 
        read 100, transa(ig), transb(ig)
      end do

      do ig = 1, grpcount 
        lda(ig) = rmaxa
        ldb(ig) = rmaxb
        ldc(ig) = rmaxc
      end do

      matidx = 1
      do ig = 1, grpcount
        do is = 1, grpsize(ig)
          if ((transa(ig).eq.'N').or.(transa(ig).eq.'n')) then
            if (m(ig).gt.rmaxa.or.k(ig).gt.cmaxa) then
              print*, ' Insufficient memory for arrays'
              goto 999
            end if
            read*, ((a(i,j,matidx),j=1,k(ig)),i=1,m(ig))
          else
            if (k(ig).gt.rmaxa.or.m(ig).gt.cmaxa) then
              print*, ' Insufficient memory for arrays'
              goto 999
            end if
            read*, ((a(i,j,matidx),j=1,m(ig)),i=1,k(ig))
          end if
          if ((transb(ig).eq.'N').or.(transb(ig).eq.'n')) then
            if (k(ig).gt.rmaxb.or.n(ig).gt.cmaxb) then
              print*, ' Insufficient memory for arrays'
              goto 999
            end if
            read*, ((b(i,j,matidx),j=1,n(ig)),i=1,k(ig))
          else
            if (n(ig).gt.rmaxb.or.k(ig).gt.cmaxb) then
              print*, ' Insufficient memory for arrays'
              goto 999
            end if
            read*, ((b(i,j,matidx),j=1,k(ig)),i=1,n(ig))
          end if
          if (m(ig).gt.rmaxc.or.n(ig).gt.cmaxc) then
              print*, ' Insufficient memory for arrays'
              goto 999
          end if
          read*, ((c(i,j,matidx),j=1,n(ig)),i=1,m(ig))

          a_array(matidx) = LOC(a(1,1,matidx))
          b_array(matidx) = LOC(b(1,1,matidx))
          c_array(matidx) = LOC(c(1,1,matidx))

          matidx = matidx + 1
        end do 
      end do
*
*       Print input data
      matidx = 1
      do ig = 1, grpcount
        print*
        print*, '     INPUT DATA'
        print 101, m(ig), n(ig), k(ig)
        print 102, alpha(ig), beta(ig)
        print 103, transa(ig), transb(ig)
        do is = 1, grpsize(ig)
          if ((transa(ig).eq.'N').or.(transa(ig).eq.'n')) then
            call PrintArrayC(0,0,m(ig),k(ig),a(:,:,matidx),lda(ig),'A')
          else
            call PrintArrayC(0,0,k(ig),m(ig),a(:,:,matidx),lda(ig),'A')
          end if
          if ((transb(ig).eq.'N').or.(transb(ig).eq.'n')) then
            call PrintArrayC(0,0,k(ig),n(ig),b(:,:,matidx),ldb(ig),'B')
          else
            call PrintArrayC(0,0,n(ig),k(ig),b(:,:,matidx),ldb(ig),'B')
          end if
          call PrintArrayC(0,0,m(ig),n(ig),c(:,:,matidx),ldc(ig),'C')
          matidx = matidx + 1
        end do
      end do
*
*      Call CGEMM3M_BATCH subroutine
      call CGEMM3M_BATCH(transa,transb,m,n,k,alpha,a_array,lda,b_array,
     $                 ldb,beta,c_array,ldc,grpcount,grpsize)
*
      print*
      print*, '     OUTPUT DATA'

      matidx = 1
      do ig = 1, grpcount
        do is = 1, grpsize(ig)
          call PrintArrayC(1,0,m(ig),n(ig),c(:,:,matidx),ldc(ig),'C')
          matidx = matidx + 1
        end do
      end do

      stop
 100  format(2(a1,1x))
 101  format(7x,'M=',i1,'  N=',i1,'  K=',i1)
 102  format(7x,'ALPHA=(',f5.2,',',f5.2,' )',
     $          '  BETA=(',f5.2,',',f5.2,' )')
 103  format(7x,'TRANSA=',a1, '  TRANSB=',a1)
 999  stop 1
      end
