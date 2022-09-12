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
!      S T R S M _ B A T C H  Example Program Text
!*******************************************************************************

      program   STRSM_BATCH_MAIN
*
      use, intrinsic :: ISO_C_BINDING
      integer        maxgc, maxmat
      parameter      (maxgc=5, maxmat=20)
      integer        rmaxa, cmaxa, rmaxb, cmaxb
      parameter      (rmaxa=5, cmaxa=5, rmaxb=5, cmaxb=5)
      integer        lda(maxgc), ldb(maxgc)
      integer        grpcount, ig, is, matidx
      integer        m(maxgc), n(maxgc)
      integer        grpsize(maxgc)
      real           alpha(maxgc)
      real           a(rmaxa,cmaxa,maxmat)
      real           b(rmaxb,cmaxb,maxmat)
      character*1    side(maxgc), uplo(maxgc)
      character*1    transa(maxgc), diag(maxgc)
      integer(KIND=C_SIZE_T) a_array(maxmat), b_array(maxmat)

*       External Subroutines
      external         STRSM_BATCH, PrintArrayS
*
*      Executable Statements
*
      print*
      print*,'   S T R S M _ B A T C H  EXAMPLE PROGRAM'
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
        read*, m(ig), n(ig)
      end do

      do ig = 1, grpcount
        read*, alpha(ig)
      end do

      do ig = 1, grpcount
        read*, side(ig), uplo(ig), transa(ig), diag(ig)
      end do

      do ig = 1, grpcount
        lda(ig) = rmaxa
        ldb(ig) = rmaxb
      end do

      matidx = 1
      do ig = 1, grpcount
        do is = 1, grpsize(ig)
          if ((side(ig).eq.'L').or.(side(ig).eq.'l')) then
            if (m(ig).gt.rmaxa.or.m(ig).gt.cmaxa.or.
     $          m(ig).gt.rmaxb.or.n(ig).gt.cmaxb ) then
              print*, ' Insufficient memory for arrays'
              goto 999
            end if
            if ((uplo(ig).eq.'U').or.(uplo(ig).eq.'u')) then
              read*, ((a(i,j,matidx),j=i,m(ig)),i=1,m(ig))
            else
              read*, ((a(i,j,matidx),j=1,i),i=1,m(ig))
            end if
          else
            if (n(ig).gt.rmaxa.or.n(ig).gt.cmaxa.or.
     $          m(ig).gt.rmaxb.or.n(ig).gt.cmaxb ) then
              print*, ' Insufficient memory for arrays'
              goto 999
            end if
            if ((uplo(ig).eq.'U').or.(uplo(ig).eq.'u')) then
              read*, ((a(i,j,matidx),j=i,n(ig)),i=1,n(ig))
            else
              read*, ((a(i,j,matidx),j=1,i),i=1,n(ig))
            end if
          end if
          read*, ((b(i,j,matidx),j=1,n(ig)),i=1,m(ig))

          a_array(matidx) = LOC(a(1,1,matidx))
          b_array(matidx) = LOC(b(1,1,matidx))

          matidx = matidx + 1
        end do
      end do
*
*       Print input data

      matidx = 1
      do ig = 1, grpcount
        print*
        print*, '     INPUT DATA'
        print 101, m(ig), n(ig)
        print 102, alpha(ig)
        print 103, side(ig), uplo(ig), transa(ig), diag(ig)
        do is = 1, grpsize(ig)
          if ((side(ig).eq.'L').or.(side(ig).eq.'l')) then
            if ((uplo(ig).eq.'U').or.(uplo(ig).eq.'u')) then
             call PrintArrayS(0,1,m(ig),m(ig),a(:,:,matidx),lda(ig),'A')
            else
             call PrintArrayS(0,-1,m(ig),m(ig),a(:,:,matidx),
     $                        lda(ig),'A')
            end if
          else
            if ((uplo(ig).eq.'U').or.(uplo(ig).eq.'u')) then
              call PrintArrayS(0,1,n(ig),n(ig),a(:,:,matidx),
     $                         lda(ig),'A')
            else
              call PrintArrayS(0,-1,n(ig),n(ig),a(:,:,matidx),
     $                         lda(ig),'A')
            end if
          end if
          call PrintArrayS(0,0,m(ig),n(ig),b(:,:,matidx),ldb(ig),'B')
          matidx = matidx + 1
        end do
      end do
*
*      Call STRSM_BATCH subroutine
      call STRSM_BATCH(side,uplo,transa,diag,m,n,alpha,a_array,lda,
     $                 b_array,ldb,grpcount,grpsize)
*
      print*
      print*, '     OUTPUT DATA'

      matidx = 1
      do ig = 1, grpcount
        do is = 1, grpsize(ig)
          call PrintArrayS(1,0,m(ig),n(ig),b(:,:,matidx),ldb(ig),'B')
          matidx = matidx + 1
        end do
      end do

      stop
 101  format(7x,'M=',i1,'  N=',i1)
 102  format(7x,'ALPHA=',f5.2)
 103  format(7x,'SIDE=',a1, '  UPLO=',a1, '  TRANSA=',a1, '  DIAG=',a1)
 999  stop 1
      end
