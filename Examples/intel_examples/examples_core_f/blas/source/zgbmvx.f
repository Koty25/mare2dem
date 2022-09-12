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
!      Z G B M V  Example Program Text
!*******************************************************************************

      program   ZGBMV_MAIN
*
      character*1      trans
      integer          m, n, kl, ku, lda, incx, incy
      integer          mmax, nmax, xmax, ymax
      parameter        (mmax=8, nmax=5, xmax=8, ymax=8)
      parameter        (lda=mmax)
      complex*16        alpha, beta
      complex*16        a(mmax,nmax), x(xmax), y(ymax)
      integer          i, j, ku1, kl1
*       Intrinsic Functions
      intrinsic        abs
*       External Subroutines
      external         ZGBMV, PrintVectorZ, PrintBandArrayZ
*
*       Executable Statements
*
      print*
      print*, '   Z G B M V  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, m, n, kl, ku
      read*, incx, incy
      read*, alpha, beta
      read 100, trans
      if ((trans.eq.'N').or.(trans.eq.'n')) then
        if ( ((1+(n-1)*abs(incx)).gt.xmax).or.
     $       ((1+(m-1)*abs(incy)).gt.ymax) ) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*,(x(i),i=1, 1+(n-1)*abs(incx))
        read*,(y(i),i=1, 1+(m-1)*abs(incy))
      else
        if ( ((1+(m-1)*abs(incx)).gt.xmax).or.
     $       ((1+(n-1)*abs(incy)).gt.ymax) ) then
          print*, ' Insufficient memory for arrays'
          goto 999
        end if
        read*,(x(i),i=1, 1+(m-1)*abs(incx))
        read*,(y(i),i=1, 1+(n-1)*abs(incy))
      end if
      if ( (kl+ku+1).gt.mmax.or.n.gt.nmax ) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      if (ku.ge.n) then
         ku1 = n-1
      else
         ku1 = ku
      end if
      if (kl.ge.m) then
         kl1 = m-1
      else
         kl1 = kl
      end if
      do i=1, ku1+1
         if ((ku1-i+m+1).ge.n) then
            k = n
         else
            k = ku1-i+m+1
         end if
         read*, (a(ku-ku1+i,j),j=ku1-i+2,k)
      end do
      do i=ku1+2, ku1+kl1+1
         if ((m+ku1-i+1).ge.n) then
            k = n
         else
            k = m+ku1-i+1
         end if
         read*, (a(ku-ku1+i,j),j=1,k)
      end do

*       Print input data
      print*
      print*, '     INPUT DATA'
      print 101, m, n
      print 102, alpha, beta
      print 103, trans
      if ((trans.eq.'N').or.(trans.eq.'n')) then
        call PrintVectorZ(0,n,x,incx,'X ')
        call PrintVectorZ(0,m,y,incy,'Y ')
      else
        call PrintVectorZ(0,m,x,incx,'X ')
        call PrintVectorZ(0,n,y,incy,'Y ')
      end if
      call PrintBandArrayZ(0,kl,ku,m,n,a,lda,'A')
*
*      Call ZGBMV subroutine
      call ZGBMV(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
*
      print*
      print*, '     OUTPUT DATA'
      if ((trans.eq.'N').or.(trans.eq.'n')) then
        call PrintVectorZ(0,m,y,incy,'Y ')
      else
        call PrintVectorZ(0,n,y,incy,'Y ')
      end if

      stop
 100  format(a1)
 101  format(7x,'M=',i1,'  N=',i1)
 102  format(7x,'ALPHA=(',f4.1,',',f4.1,')  BETA=(',f4.1,',',f4.1,')')
 103  format(7x,'TRANS=',a1)
 999  stop 1
      end
