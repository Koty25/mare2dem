!===============================================================================
! Copyright 2010-2020 Intel Corporation.
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
!      D G E M 2 V U  Example Program Text
!*******************************************************************************

      program   DGEM2VU_MAIN
*
      integer          m, n, lda, incx1, incy1, incx2, incy2
      integer          mmax, nmax, x1max, y1max, x2max, y2max
      parameter        (mmax=4, nmax=5)
      parameter        (x1max=5, y1max=5, x2max=5, y2max=5)
      parameter        (lda=mmax)
      double precision alpha, beta
      double precision a(mmax,nmax)
      double precision x1(x1max), y1(y1max), x2(x2max), y2(y2max)
      integer          i, j
*       Intrinsic Functions
      intrinsic        abs
*       External Subroutines
      external         DGEM2VU, PrintVectorD, PrintArrayD
*
*       Executable Statements
*
      print*
      print*, '   D G E M 2 V U  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, m, n
      read*, incx1, incy1, incx2, incy2
      read*, alpha, beta
      if ( ((1+(n-1)*abs(incx1)).gt.x1max).or.
     $     ((1+(m-1)*abs(incy1)).gt.y1max).or.
     $     ((1+(m-1)*abs(incx2)).gt.x2max).or.
     $     ((1+(n-1)*abs(incy2)).gt.y2max) ) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*,(x1(i),i=1, 1+(n-1)*abs(incx1))
      read*,(y1(i),i=1, 1+(m-1)*abs(incy1))
      read*,(x2(i),i=1, 1+(m-1)*abs(incx2))
      read*,(y2(i),i=1, 1+(n-1)*abs(incy2))
      if (m.gt.mmax.or.n.gt.nmax) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, ((a(i,j),j=1,n),i=1,m)
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 101, m, n
      print 102, alpha, beta
      call PrintVectorD(0,n,x1,incx1,'X1 ')
      call PrintVectorD(0,m,y1,incy1,'Y1 ')
      call PrintVectorD(0,m,x2,incx2,'X2 ')
      call PrintVectorD(0,n,y2,incy2,'Y2 ')
      call PrintArrayD(0,0,m,n,a,lda,'A')
*
*      Call DGEM2VU subroutine
      call DGEM2VU(m,n,alpha,a,lda,x1,incx1,x2,incx2,
     $             beta,y1,incy1,y2,incy2)
*
      print*
      print*, '     OUTPUT DATA'
      call PrintVectorD(0,m,y1,incy1,'Y1 ')
      call PrintVectorD(0,n,y2,incy2,'Y2 ')

      stop
 101  format(7x,'M=',i1,'  N=',i1)
 102  format(7x,'ALPHA=',f5.2,'  BETA=',f5.2)
 999  stop 1
      end
