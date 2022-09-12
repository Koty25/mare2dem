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

!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Sparse format
!             converters FORTRAN example
!
!*******************************************************************************
!
! Example program for using Intel MKL Sparse format converters
! The following Sparse  Sparse format converters are used in the example:
!
!          MKL_CDNSCSR
!          MKL_CCSRCOO
!          MKL_CCSRBSR
!          MKL_CCSRDIA
!          MKL_CCSRCSC
!          MKL_CCSRSKY
!
!*******************************************************************************

       program sparseformats
       IMPLICIT NONE

c*******************************************************************************
c     Definition arrays for sparse matrix formats
c*******************************************************************************

      integer m,n,lda,ldAbsr,nzmax,nnz,mblk,nn,mn,idiag,ndiag
      parameter ( m=4, n=4,lda=4, nzmax=8,mblk=2,
     &            ldAbsr=4, nn=2, mn=16, idiag=3, ndiag=4)
      complex ADNS(m,n),   ADNS_standard(m,n)
      complex Absr(nzmax), Absr_standard(nzmax)
      complex Acsc(nzmax), Acsc_standard(nzmax)
      complex Acsr(nzmax), Acsr_standard(nzmax)
      complex Acoo(nzmax), Acoo_standard(nzmax)
      complex Adia(idiag*ndiag), Adia_standard(idiag*ndiag)
      complex Asky(6)
      complex Askyl_standard(6)
      integer pointers(M+1)
      integer pointersl_standard(M+1)
      integer AI(M+1), AI_standard(M+1)
      integer AJ(nzmax), AJ_standard(nzmax)
      integer AJL_standard(6)
      integer AI1(M+1), AI1_standard(M+1)
      integer AJ1(nzmax), AJ1_standard(nzmax)
      integer AJB(MBLK), AJB_standard(MBLK)
      integer AIB(MBLK+1), AIB_standard(MBLK+1)
      integer IR(nzmax), IR_standard(nzmax)
      integer JC(nzmax), JC_standard(nzmax)
      integer distance(idiag), distance_standard(idiag)
      data Absr_standard / (5.0, 1.0), (9.0, 1.0),
     &                     (8.0, 1.0), (2.0, 1.0), 
     &                     (3.0, 1.0), (1.0, 1.0),
     &                     (6.0, 1.0), (4.0, 1.0) /
      data Acsc_standard / (5.0, 1.0), (9.0, 1.0),
     &                     (8.0, 1.0), (2.0, 1.0), 
     &                     (3.0, 1.0), (1.0, 1.0),
     &                     (6.0,1.0), (4.0, 1.0) /
      data Acsr_standard / (5.0,1.0), (8.0, 1.0),
     &                     (9.0, 1.0), (2.0, 1.0), 
     &                     (3.0, 1.0), (6.0, 1.0),
     &                     (1.0, 1.0), (4.0, 1.0)/
      data Acoo_standard / (5.0, 1.0), (8.0, 1.0),
     &                     (9.0, 1.0), (2.0, 1.0), 
     &                     (3.0, 1.0), (6.0, 1.0),
     &                     (1.0, 1.0), (4.0, 1.0)/
      data Adia_standard / (0.0, 0.0), (9.0, 1.0),
     &                     (0.0, 0.0), (1.0, 1.0),
     &                     (5.0, 1.0), (2.0, 1.0),
     &                     (3.0, 1.0), (4.0, 1.0),
     &                     (8.0, 1.0), (0.0, 0.0),
     &                     (6.0, 1.0), (0.0, 0.0)/
      data Askyl_standard /(5.0, 1.0), (9.0, 1.0),
     &                     (2.0, 1.0), (3.0, 1.0), 
     &                     (1.0, 1.0), (4.0,1.0)/
      data AI_standard / 1, 3, 5, 7, 9/
      data AI1_standard / 1, 3, 5, 7, 9/
      data AJ_standard / 1, 2, 1, 2, 3, 4, 3, 4/
      data AJ1_standard / 1, 2, 1, 2, 3, 4, 3, 4/
      data AJB_standard / 1,2/
      data AIB_standard /1, 2, 3/
      data IR_standard /1, 1, 2, 2, 3, 3, 4, 4/
      data JC_standard /1, 2, 1, 2, 3, 4, 3, 4/
      data distance_standard /-1, 0, 1/
      data pointersl_standard /1, 2, 4, 5, 7/
      data AJL_standard /1, 1, 2, 3, 3, 4/


c*************************************************************************************************
c    Declaration of local variables :
c*************************************************************************************************
      integer ifail,i,j,ij,info,id,nd
      integer nr,ibase1,ibase2,locat
      real  rfail
      integer job(8)

      print *, 'EXAMPLE PROGRAM FOR CONVERTER FROM ONE'
      print *, '       SPARSE FORMAT TO OTHER'
      ifail=0
      rfail=0.0
      info = 0
      locat = 2
      ibase1 =1
      ibase2 =1
      job(2)=ibase1
      job(3)=ibase2
      job(4)=locat
      job(5)=nzmax
c***************************************************************************************************
c TASK 1    Obtain compressed sparse row matrix from dense matrix
c***************************************************************************************************
      do j=1,n
        do i=1,m
           ADNS(i,j)=cmplx(0.0, 0.0)
        enddo
      enddo
      ADNS(1,1)=cmplx(5.0, 1.0)
      ADNS(1,2)=cmplx(8.0, 1.0)
      ADNS(2,1)=cmplx(9.0, 1.0)
      ADNS(2,2)=cmplx(2.0, 1.0)
      ADNS(3,3)=cmplx(3.0, 1.0)
      ADNS(3,4)=cmplx(6.0, 1.0)
      ADNS(4,3)=cmplx(1.0, 1.0)
      ADNS(4,4)=cmplx(4.0, 1.0)
      do j=1,n
        do i=1,m
           ADNS_standard(i,j)=ADNS(i,j)
        enddo
      enddo
      do i=1,n+1
         AI(i)=AI_standard(i)
      enddo
      job(1)=0
      job(6)=1
      call mkl_cdnscsr(job,m,n,Adns,lda,Acsr,AJ,AI,info)
      if (info.ne.0) goto 101
      do i=1,n+1
         ifail=AI(i)-AI_standard(i)
         if (ifail.ne.0) goto 101
      enddo
      do i=1,nzmax
         ifail=AJ(i)-AJ_standard(i)
         if (ifail.ne.0) goto 101
      enddo
      do i=1,nzmax
         rfail=abs(Acsr(i)-Acsr_standard(i))
         if (rfail.ne.0.0) goto 101
      enddo
c***************************************************************************************************
c TASK 2    Obtain dense matrix from compressed compressed sparse row matrix
c***************************************************************************************************
      do j=1,n
        do i=1,m
           ADNS(i,j)=cmplx(0.0, 0.0)
        enddo
      enddo
      do i=1,n+1
         AI(i)=AI_standard(i)
      enddo
      do i=1,nzmax
         AJ(i)=AJ_standard(i)
      enddo
      do i=1,nzmax
         Acsr(i)=Acsr_standard(i)
      enddo
      job(1)=1
      call mkl_cdnscsr(job,m,n,Adns,lda,Acsr,AJ,AI,info)

      if (info.ne.0) goto 102
      do j=1,n
         do i=1,m
            rfail=abs(ADNS(i,j)-ADNS_standard(i,j))
            if (rfail.ne.0.0) goto 102
         enddo
      enddo
c***************************************************************************************************
c TASK 3    Obtain sparse coordinate matrix from compressed compressed sparse row matrix
c***************************************************************************************************
      job(1)=0
      job(6)=3
      call mkl_ccsrcoo (job,n, Acsr, AJ,AI,nnz,Acoo, ir,jc,info)

      if (info.ne.0) goto 103
      do i=1,nzmax
         ifail=IR(i)-IR_standard(i)
         if (ifail.ne.0) goto 103
      enddo
      do i=1,nzmax
         ifail=JC(i)-JC_standard(i)
         if (ifail.ne.0) goto 103
      enddo
      do i=1,nzmax
         rfail=abs(Acoo(i)-Acoo_standard(i))
         if (rfail.ne.0.0) goto 103
      enddo
c***************************************************************************************************
c TASK 4    Obtain compressed compressed sparse row matrix from sparse coordinate matrix
c***************************************************************************************************
      job(1)=1
      job(6)=2
      call mkl_ccsrcoo (job,n, Acsr, AJ,AI,nnz,Acoo, ir,jc,info)

      if (info.ne.0) goto 104
      do i=1,n+1
        ifail=AI(i)-AI_standard(i)
        if (ifail.ne.0) goto 104
      enddo
      do i=1,nzmax
         ifail=AJ(i)-AJ_standard(i)
         if (ifail.ne.0) goto 104
      enddo
      do i=1,nzmax
         rfail=abs(Acsr(i)-Acsr_standard(i))
         if (rfail.ne.0.0) goto 104
      enddo
c***************************************************************************************************
c TASK 5   Obtain block sparse row matrix from compressed sparse row matrix
c***************************************************************************************************
      do i=1,nzmax
         Absr(i)=0.0
      enddo
      job(1) = 0
      job(6) = 1

      call mkl_ccsrbsr (job,m,mblk,ldAbsr,Acsr,AJ,AI,Absr,AJB,AIB,info)

      nr = 1 + (m-1) / mblk
      if (info.ne.0) goto 105
      do i=1,nr+1
         ifail=AIB(i)-AIB_standard(i)
         if (ifail.NE.0) goto 105
      enddo
      do i=1,mblk
         ifail=AJB(i)-AJB_standard(i)
         if (ifail.NE.0) goto 105
      enddo
      do i=1,nzmax
         rfail=abs(Absr(i)-Absr_standard(i))
         if (rfail.NE.0) goto 105
      enddo
c***************************************************************************************************
c TASK 6   Obtain compressed sparse row matrix from block sparse row matrix
c***************************************************************************************************
      job(1) = 1
      job(6) = 3
      call mkl_ccsrbsr (job,nn,mblk,ldAbsr,Acsr,AJ,AI,Absr,AJB,AIB,info)

      if (info.ne.0) goto 106
      do i=1,n+1
         ifail=AI(i)-AI_standard(i)
         if (ifail.ne.0) goto 106
      enddo
      do i=1,nzmax
         ifail=AJ(i)-AJ_standard(i)
         if (ifail.ne.0) goto 106
      enddo
      do i=1,nzmax
         rfail=abs(Acsr(i)-Acsr_standard(i))
         if (rfail.ne.0.0) goto 106
      enddo
c***************************************************************************************************
c TASK 7    Obtain sparse diagonal matrix from compressed sparse row matrix
c***************************************************************************************************
      do i=1,idiag
         distance(i)=distance_standard(i)
      enddo
      job(1) = 0
      job(6) = 0
      id = idiag
      nd = ndiag
      call mkl_ccsrdia(job,n,Acsr,AJ,AI,Adia,nd,distance,id,Acsr,
     &                 AJ,AI,info)

      if (info.ne.0) goto 107
      do i=1,idiag
         ifail=distance(i)-distance_standard(i)
         if (ifail.NE.0) goto 107
      enddo
      do j=1,idiag
         do i=1,ndiag
            ij = i + ndiag*j
            rfail=abs(Adia(i)-Adia_standard(i))
            if (rfail.NE.0) goto 107
         enddo
      enddo
c***************************************************************************************************
c TASK 8    Obtain compressed sparse  row matrix from sparse diagonal matrix
c***************************************************************************************************
      job(1) = 1
      job(6) = 0
      call mkl_ccsrdia(job,n,Acsr,AJ,AI,Adia,nd,distance,id,Acsr,
     &                 AJ,AI,info)

      if (info.ne.0) goto 108
      do i=1,n+1
         ifail=AI(i)-AI_standard(i)
         if (ifail.ne.0) goto 108
      enddo
      do i=1,nzmax
         ifail=AJ(i)-AJ_standard(i)
         if (ifail.ne.0) goto 108
      enddo
      do i=1,nzmax
         rfail=abs(Acsr(i)-Acsr_standard(i))
         if (rfail.ne.0.0) goto 108
      enddo
c***************************************************************************************************
c TASK 9    Obtain compressed sparse column matrix from compressed sparse row matrix
c***************************************************************************************************
      job(1) = 0
      job(6) = 1
      call mkl_ccsrcsc(job,n,Acsr,AJ,AI,Acsc,AJ1,AI1,info)

      if (info.ne.0) goto 109
      do i=1,n+1
         ifail=AI1(i)-AI1_standard(i)
         if (ifail.ne.0) goto 109
      enddo
      do i=1,nzmax
         ifail=AJ1(i)-AJ1_standard(i)
         if (ifail.ne.0) goto 109
      enddo
      do i=1,nzmax
         rfail=abs(Acsc(i)-Acsc_standard(i))
         if (rfail.ne.0.0) goto 109
      enddo
c***************************************************************************************************
c TASK 10   Obtain compressed  sparse  row matrix from compressed sparse column matrix
c***************************************************************************************************
      job(1) = 1
      job(6) = 1
      call mkl_ccsrcsc(job,n,Acsr,AJ,AI,Acsc,AJ1,AI1,info)

      if (info.ne.0) goto 110
      do i=1,n+1
         ifail=AI(i)-AI_standard(i)
         if (ifail.ne.0) goto 110
      enddo
      do i=1,nzmax
         ifail=AJ(i)-AJ_standard(i)
         if (ifail.ne.0) goto 110
      enddo
      do i=1,nzmax
         rfail=abs(Acsr(i)-Acsr_standard(i))
         if (rfail.ne.0.0) goto 110
      enddo
c***************************************************************************************************
c TASK 11    Obtain sparse skyline matrix for lower triangle from compressed sparse row matrix
c***************************************************************************************************
      job(1) = 0
      job(4) = 0
      job(5) = 6
      job(6) = 0
      call mkl_ccsrsky(job,n,Acsr,AJ,AI,Asky,pointers,info)

      if (info.ne.0) goto 111
      do i=1,n+1
         ifail=pointers(i)-pointersl_standard(i)
         if (ifail.ne.0) goto 111
      enddo
      nnz = pointers(n+1)-pointers(1);
      do i=1,nnz
         rfail=abs(Asky(i)-Askyl_standard(i))
         if (rfail.ne.0.0) goto 111
      enddo
c***************************************************************************************************
c TASK 12   Obtain compressed sparse row matrix for lower triangle from sparse skyline matrix
c***************************************************************************************************
      job(1) = 1
      job(4) = 0
      job(6) = 0
      call mkl_ccsrsky(job,n,Acsr,AJ,AI,Asky,pointers,info)

      if (info.ne.0) goto 112
      do i=1,n+1
         ifail=AI(i)-pointersl_standard(i)
         if (ifail.ne.0) goto 112
      enddo
      nnz = pointers(n+1)-pointers(1);
      do i=1,nnz
         ifail=AJ(i)-AJL_standard(i)
         if (ifail.ne.0) goto 112
      enddo
      do i=1,nnz
         rfail=abs(Acsr(i)-Askyl_standard(i))
         if (rfail.ne.0.0) goto 112
      enddo
c****************************************************************************************************
      print *,'         ALL EXAMPLES PASSED'
      stop 0

c FAILURE message to print if something went wrong

 101  print *,('Example FAILED to convert from dns to csr...')
      stop 1
 102  print *,('Example FAILED to convert from csr to dns...')
      stop 1
 103  print *,('Example FAILED to convert from csr to coo...')
      stop 1
 104  print *,('Example FAILED to convert from coo to csr...')
      stop 1
 105  print *,('Example FAILED to convert from bsr to csr...')
      stop 1
 106  print *,('Example FAILED to convert from csr to bsr...')
      stop 1
 107  print *,('Example FAILED to convert from csr to dia...')
      stop 1
 108  print *,('Example FAILED to convert from dia to csr...')
      stop 1
 109  print *,('Example FAILED to convert from csr to csc...')
      stop 1
 110  print *,('Example FAILED to convert from csc to csr...')
      stop 1
 111  print *,('Example FAILED to convert from csr to sky...')
      stop 1
 112  print *,('Example FAILED to convert from sky to csr...')
      stop 1
      end program sparseformats
