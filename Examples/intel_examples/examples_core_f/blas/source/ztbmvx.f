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
!      Z T B M V  Example Program Text
!*******************************************************************************

      program   ZTBMV_MAIN
*
      character*1      uplo, trans, diag
      integer          n, k, lda, incx
      integer          mmax, nmax, xmax
      parameter        (mmax=5, nmax=6, xmax=5)
      parameter        (lda=mmax)
      complex*16       a(mmax, nmax), x(xmax)
      integer          i, j, k1
*       Intrinsic Functions
      intrinsic        abs
*       External Subroutines
      external         ZTBMV, PrintVectorZ, PrintBandArrayZ
*
*       Executable Statements
*
      print*
      print*, '   Z T B M V  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, n, k
      read*, incx
      read 100, uplo, trans, diag
      if ( (1+(n-1)*abs(incx)).gt.xmax ) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, (x(i),i=1, 1+(n-1)*abs(incx))
      print 101, n, k
      print 102, uplo, trans, diag
      if ( (k+1).gt.mmax.or.n.gt.nmax ) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      if (k.ge.n) then
         k1 = n-1
      else
         k1 = k
      end if
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        read*, ((a(k-k1+i,j),j=k1+2-i,n),i=1,k1+1)
      else
        do i=1, k1+1
           read*, (a(i,j),j=1,n+1-i)
        end do
      end if
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 101, n, k
      print 102, uplo, trans, diag
      call PrintVectorZ(0,n,x,incx,'X ')
      if ((uplo.eq.'U').or.(uplo.eq.'u')) then
        call PrintBandArrayZ(1,0,k,n,n,a,lda,'A')
      else
        call PrintBandArrayZ(1,k,0,n,n,a,lda,'A')
      end if
*
*      Call ZTBMV subroutine
      call ZTBMV( uplo, trans, diag, n, k, a, lda, x, incx )
*
      print*
      print*, '     OUTPUT DATA'
      call PrintVectorZ(1,n,x,incx,'X ')

      stop
 100  format(3(a1,1x))
 101  format(7x,'N=',i1,'  K=',i1)
 102  format(7x,'UPLO=',a1,'  TRANS=',a1,'  DIAG=',a1)
 999  stop 1
      end
