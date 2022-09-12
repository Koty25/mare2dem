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
!      S S B M V  Example Program Text
!*******************************************************************************

      program   SSBMV_MAIN
*
      character*1      uplo
      integer          n, k, lda, incx, incy
      integer          mmax, nmax, xmax, ymax
      parameter        (mmax=5, nmax=6, xmax=5, ymax=5)
      parameter        (lda=mmax)
      real             alpha, beta
      real             a(mmax,nmax), x(xmax), y(ymax)
      integer          i, j, k1
*       Intrinsic Functions
      intrinsic        abs
*       External Subroutines
      external         SSBMV, PrintVectorS, PrintBandArrayS
*
*       Executable Statements
*
      print*
      print*, '   S S B M V  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, n, k
      read*, incx, incy
      read*, alpha, beta
      read 100, uplo
      if ( ((1+(n-1)*abs(incx)).gt.xmax).or.
     $     ((1+(n-1)*abs(incy)).gt.ymax) ) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, (x(i),i=1, 1+(n-1)*abs(incx))
      read*, (y(i),i=1, 1+(n-1)*abs(incy))
      if ( (k+1).gt.mmax.or.n.gt.nmax ) goto 999
      if (k.ge.n) then
         k1 = n-1
      else
         k1 = k
      end if
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        read*, ((a(k-k1+i,j),j=k1+2-i,n),i=1,k1+1)
      else
        do i = 1, k1+1
           read*, (a(i,j),j=1,n+1-i)
        end do
      end if
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 101, n, k
      print 102, alpha, beta
      print 103, uplo
      call PrintVectorS(0,n,x,incx,'X ')
      call PrintVectorS(0,n,y,incy,'Y ')
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        call PrintBandArrayS(1,0,k,n,n,a,lda,'A')
      else
        call PrintBandArrayS(1,k,0,n,n,a,lda,'A')
      end if
*
*      Call SSBMV subroutine
      call SSBMV( uplo, n, k, alpha, a, lda, x, incx,
     $            beta, y, incy )
*
      print*
      print*, '     OUTPUT DATA'
      call PrintVectorS(1,n,y,incy,'Y ')

      stop
 100  format(a1)
 101  format(7x,'N=',i1,'  K=',i1)
 102  format(7x,'ALPHA=',f4.1,'  BETA=',f4.1)
 103  format(7x,'UPLO=',a1)
 999  stop 1
      end
