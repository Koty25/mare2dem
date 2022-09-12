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

!  Content:
!
!*******************************************************************************

      integer function MaxValue(n, x)
      integer    n, x(*)
      integer    i, value

      value = x(1)
      do  i = 2, n
         if (x(i).gt.value) then
             value = x(i)
         end if
      end do
      MaxValue = value
      return
      end

      subroutine PrintVectorS(flag, n, x, incx, name)
      integer       flag, n, incx
      character*2   name
      real          x(*)

*       Intrinsic Functions
      intrinsic        abs

      if (flag.eq.0) then
        print 100, name, name, incx
      else
        print 101, name
      end if
      print 102, (x(i),i=1,1+(n-1)*abs(incx))

 100  format(7x,'VECTOR ',a2,'   INC',a1,'=',i2)
 101  format(7x,'VECTOR ',a2)
 102  format(9x,10(f6.2,2x))
      return
      end

      subroutine PrintVectorD(flag, n, x, incx, name)
      integer           flag, n, incx
      character*2       name
      double precision  x(*)

*       Intrinsic Functions
      intrinsic         abs

      if (flag.eq.0) then
        print 100, name, name, incx
      else
        print 101, name
      end if
      print 102, (x(i),i=1,1+(n-1)*abs(incx))

 100  format(7x,'VECTOR ',a2,'   INC',a1,'=',i2)
 101  format(7x,'VECTOR ',a2)
 102  format(9x,10(f8.3,2x))
      return
      end

      subroutine PrintVectorC(flag, n, x, incx, name)
      integer       flag, n, incx
      character*2   name
      complex       x(*)

*       Intrinsic Functions
      intrinsic        abs

      if (flag.eq.0) then
        print 100, name, name, incx
      else
        print 101, name
      end if
      print 102, (x(i),i=1,1+(n-1)*abs(incx))

 100  format(7x,'VECTOR ',a2,'   INC',a1,'=',i2)
 101  format(7x,'VECTOR ',a2)
 102  format(9x,10(f6.2,',',f6.2,3x))
      return
      end

      subroutine PrintVectorZ(flag, n, x, incx, name)
      integer       flag, n, incx
      character*2   name
      complex*16    x(*)

*       Intrinsic Functions
      intrinsic        abs

      if (flag.eq.0) then
        print 100, name, name, incx
      else
        print 101, name
      end if
      print 102, (x(i),i=1,1+(n-1)*abs(incx))

 100  format(7x,'VECTOR ',a2,'   INC',a1,'=',i2)
 101  format(7x,'VECTOR ',a2)
 102  format(7x,10(f6.2,',',f6.2,3x))
      return
      end

      subroutine PrintArrayS(flag1, flag2, m, n, a, lda, name)
      integer       flag1, flag2, m, n, lda
      character*1   name
      real          a(lda,*)
      integer       i, j

      if (flag1.eq.0) then
        print 100, name, name, lda
      else
        print 101, name
      end if
      if (flag2.eq.0) then
        do i=1,m
            print 110, (a(i,j),j=1,n)
        end do
      else if (flag2.eq.1) then
        do i=1, m
             goto (10,20,30,40,50) i
  10       print 110, (a(i,j),j=i,m)
           goto 60
  20       print 120, (a(i,j),j=i,m)
           goto 60
  30       print 130, (a(i,j),j=i,m)
           goto 60
  40       print 140, (a(i,j),j=i,m)
           goto 60
  50       print 150, (a(i,j),j=i,m)
  60       continue
        end do
      else if(flag2.eq.-1) then
          do i=1, m
             print 110, (a(i,j),j=1,i)
          end do
        end if

 100  format(7x,'ARRAY ',a1,'   LD',a1,'=',i1)
 101  format(7x,'ARRAY ',a1)
 110  format(9x,10(f6.2,2x))
 120  format(17x,10(f6.2,2x))
 130  format(25x,10(f6.2,2x))
 140  format(33x,10(f6.2,2x))
 150  format(41x,10(f6.2,2x))
      return
      end

      subroutine PrintArrayD(flag1, flag2, m, n, a, lda, name)
      integer           flag1, flag2, m, n, lda
      character*1       name
      double precision  a(lda,*)
      integer           i, j

      if (flag1.eq.0) then
        print 100, name, name, lda
      else
        print 101, name
      end if
      if (flag2.eq.0) then
        do i=1,m
            print 110, (a(i,j),j=1,n)
        end do
      else if (flag2.eq.1) then
        do i=1, m
           goto (10,20,30,40,50) i
  10       print 110, (a(i,j),j=i,m)
           goto 60
  20       print 120, (a(i,j),j=i,m)
           goto 60
  30       print 130, (a(i,j),j=i,m)
           goto 60
  40       print 140, (a(i,j),j=i,m)
           goto 60
  50       print 150, (a(i,j),j=i,m)
  60       continue
        end do
      else if(flag2.eq.-1) then
          do i=1, m
             print 110, (a(i,j),j=1,i)
          end do
        end if

 100  format(7x,'ARRAY ',a1,'   LD',a1,'=',i1)
 101  format(7x,'ARRAY ',a1)
 110  format(9x,10(f8.3,2x))
 120  format(19x,10(f8.3,2x))
 130  format(29x,10(f8.3,2x))
 140  format(39x,10(f8.3,2x))
 150  format(49x,10(f8.3,2x))
      return
      end

      subroutine PrintArrayC(flag1, flag2, m, n, a, lda, name)
      integer           flag1, flag2, m, n, lda
      character*1       name
      complex           a(lda,*)
      integer           i, j

      if (flag1.eq.0) then
        print 100, name, name, lda
      else
        print 101, name
      end if
      if (flag2.eq.0) then
        do i=1,m
            print 110, (a(i,j),j=1,n)
        end do
      else if (flag2.eq.1) then
        do i=1, m
           goto (10,20,30,40,50) i
  10       print 110, (a(i,j),j=i,m)
           goto 60
  20       print 120, (a(i,j),j=i,m)
           goto 60
  30       print 130, (a(i,j),j=i,m)
           goto 60
  40       print 140, (a(i,j),j=i,m)
           goto 60
  50       print 150, (a(i,j),j=i,m)
  60       continue
        end do
      else if(flag2.eq.-1) then
          do i=1, m
             print 110, (a(i,j),j=1,i)
          end do
      end if

 100  format(7x,'ARRAY ',a1,'   LD',a1,'=',i1)
 101  format(7x,'ARRAY ',a1)
 110  format(9x,10(f6.2,',',f6.2,3x))
 120  format(25x,10(f6.2,',',f6.2,3x))
 130  format(41x,10(f6.2,',',f6.2,3x))
 140  format(57x,10(f6.2,',',f6.2,3x))
 150  format(73x,10(f6.2,',',f6.2,3x))
      return
      end

      subroutine PrintArrayZ(flag1, flag2, m, n, a, lda, name)
      integer           flag1, flag2, m, n, lda
      character*1       name
      double complex    a(lda,*)
      integer           i, j

      if (flag1.eq.0) then
        print 100, name, name, lda
      else
        print 101, name
      end if
      if (flag2.eq.0) then
        do i=1,m
            print 110, (a(i,j),j=1,n)
        end do
      else if (flag2.eq.1) then
        do i=1, m
           goto (10,20,30,40,50) i
  10       print 110, (a(i,j),j=i,m)
           goto 60
  20       print 120, (a(i,j),j=i,m)
           goto 60
  30       print 130, (a(i,j),j=i,m)
           goto 60
  40       print 140, (a(i,j),j=i,m)
           goto 60
  50       print 150, (a(i,j),j=i,m)
  60       continue
        end do
      else if(flag2.eq.-1) then
          do i=1, m
             print 110, (a(i,j),j=1,i)
          end do
      end if

 100  format(7x,'ARRAY ',a1,'   LD',a1,'=',i1)
 101  format(7x,'ARRAY ',a1)
 110  format(9x,10(f6.2,',',f6.2,3x))
 120  format(25x,10(f6.2,',',f6.2,3x))
 130  format(41x,10(f6.2,',',f6.2,3x))
 140  format(57x,10(f6.2,',',f6.2,3x))
 150  format(73x,10(f6.2,',',f6.2,3x))
      return
      end

      subroutine PrintBandArrayS(flag1, kl, ku, m, n, a, lda, name)
      integer       flag1, kl, ku, m, n, lda
      character*1   name
      real          a(lda,*)
      integer       i, j, k, k1, ku1, kl1, Nrow

      if (flag1.eq.0) then
        print 100, name, lda, kl, ku
      else
        print 101, name, lda
      end if
      if (ku.ge.n) then
         ku1 = n-1
      else
         ku1 = ku
      end if
      print*, '       N row'
      Nrow = 1
      do i = 1, ku-n+1
          print 102, Nrow
          Nrow = Nrow + 1
      end do
      k = ku1+1
      do i=1, ku1+1
         if ((ku1-i+m+1).ge.n) then
            k1 = n
         else
            k1 = ku1-i+m+1
         end if
         goto (10,20,30,40,50) k
  10     print 110, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  20     print 120, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  30     print 130, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  40     print 140, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  50     print 150, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
  60     continue
         Nrow = Nrow + 1
         k = k-1
      end do
      if (kl.ge.m) then
         kl1 = m-1
      else
         kl1 = kl
      end if
      do i=ku1+2, ku1+kl1+1
         if ((m+ku1-i+1).ge.n) then
            k1 = n
         else
            k1 = m+ku1-i+1
         end if
         print 110, Nrow, (a(ku-ku1+i,j),j=1,k1)
         Nrow = Nrow + 1
      end do

 100  format(7x,'ARRAY ',a1,'   LDA=',i1,'  KL=',i1,'  KU=',i1)
 101  format(7x,'ARRAY ',a1,'   LDA=',i1)
 102  format(7x,i2)
 110  format(7x,i2, 2x,10(f6.2,2x))
 120  format(7x,i2,10x,10(f6.2,2x))
 130  format(7x,i2,18x,10(f6.2,2x))
 140  format(7x,i2,26x,10(f6.2,2x))
 150  format(7x,i2,34x,10(f6.2,2x))
      return
      end

      subroutine PrintBandArrayD(flag1, kl, ku, m, n, a, lda, name)
      integer          flag1, kl, ku, n, lda
      character*1      name
      double precision a(lda,*)
      integer          i, j, k, k1, ku1, kl1, Nrow

      if (flag1.eq.0) then
        print 100, name, lda, kl, ku
      else
        print 101, name, lda
      end if
      if (ku.ge.n) then
          ku1 = n-1
      else
          ku1 = ku
      end if
      print*, '       N row'
      Nrow = 1
      do i = 1, ku-n+1
          print 102, Nrow
          Nrow = Nrow + 1
      end do
      k = ku1+1
      do i=1, ku1+1
         if ((ku1-i+m+1).ge.n) then
            k1 = n
         else
            k1 = ku1-i+m+1
         end if
         goto (10,20,30,40,50) k
  10     print 110, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  20     print 120, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  30     print 130, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  40     print 140, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  50     print 150, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
  60     continue
         Nrow = Nrow + 1
         k = k-1
      end do
      if (kl.ge.m) then
         kl1 = m-1
      else
         kl1 = kl
      end if
      do i=ku1+2, ku1+kl1+1
         if ((m+ku1-i+1).ge.n) then
            k1 = n
         else
            k1 = m+ku1-i+1
         end if
         print 110, Nrow, (a(ku-ku1+i,j),j=1,k1)
         Nrow = Nrow + 1
      end do

 100  format(7x,'ARRAY ',a1,'   LDA=',i1,'  KL=',i1,'  KU=',i1)
 101  format(7x,'ARRAY ',a1,'   LDA=',i1)
 102  format(7x,i2)
 110  format(7x,i2, 2x,10(f8.3,2x))
 120  format(7x,i2,12x,10(f8.3,2x))
 130  format(7x,i2,22x,10(f8.3,2x))
 140  format(7x,i2,32x,10(f8.3,2x))
 150  format(7x,i2,42x,10(f8.3,2x))
      return
      end

      subroutine PrintBandArrayC(flag1, kl, ku, m, n, a, lda, name)
      integer       flag1, kl, ku, n, lda
      character*1   name
      complex       a(lda,*)
      integer       i, j, k, k1, ku1, kl1, Nrow

      if (flag1.eq.0) then
        print 100, name, lda, kl, ku
      else
        print 101, name, lda
      end if
      if (ku.ge.n) then
         ku1 = n-1
      else
         ku1 = ku
      end if
      print*, '       N row'
      Nrow = 1
      do i = 1, ku-n+1
          print 102, Nrow
          Nrow = Nrow + 1
      end do
      k = ku1+1
      do i=1, ku1+1
         if ((ku1-i+m+1).ge.n) then
            k1 = n
         else
            k1 = ku1-i+m+1
         end if
         goto (10,20,30,40,50) k
  10     print 110, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  20     print 120, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  30     print 130, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  40     print 140, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  50     print 150, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
  60     continue
         Nrow = Nrow + 1
         k = k-1
      end do
      if (kl.ge.m) then
          kl1 = m-1
      else
          kl1 = kl
      end if
      do i=ku1+2, ku1+kl1+1
         if ((m+ku1-i+1).ge.n) then
            k1 = n
         else
            k1 = m+ku1-i+1
         end if
         print 110, Nrow, (a(ku-ku1+i,j),j=1,k1)
         Nrow = Nrow + 1
      end do

 100  format(7x,'ARRAY ',a1,'   LDA=',i1,'  KL=',i1,'  KU=',i1)
 101  format(7x,'ARRAY ',a1,'   LDA=',i1)
 102  format(7x,i2)
 110  format(7x,i2, 2x,10(f6.2,',',f6.2,3x))
 120  format(7x,i2,18x,10(f6.2,',',f6.2,3x))
 130  format(7x,i2,34x,10(f6.2,',',f6.2,3x))
 140  format(7x,i2,50x,10(f6.2,',',f6.2,3x))
 150  format(7x,i2,66x,10(f6.2,',',f6.2,3x))
      return
      end

      subroutine PrintBandArrayZ(flag1, kl, ku, m, n, a, lda, name)
      integer       flag1, kl, ku, n, lda
      character*1   name
      complex*16    a(lda,*)
      integer       i, j, k, k1, ku1, kl1, Nrow

      if (flag1.eq.0) then
        print 100, name, lda, kl, ku
      else
        print 101, name, lda
      end if
      if (ku.ge.n) then
          ku1 = n-1
      else
          ku1 = ku
      end if
      print*, '       N row'
      Nrow = 1
      do i = 1, ku-n+1
          print 102, Nrow
          Nrow = Nrow + 1
      end do
      k = ku1+1
      do i=1, ku1+1
         if ((ku1-i+m+1).ge.n) then
            k1 = n
         else
            k1 = ku1-i+m+1
         end if
         goto (10,20,30,40,50) k
  10     print 110, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  20     print 120, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  30     print 130, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  40     print 140, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
         goto 60
  50     print 150, Nrow, (a(ku-ku1+i,j),j=ku1+2-i,k1)
  60     continue
         Nrow = Nrow + 1
         k = k-1
      end do
      if (kl.ge.m) then
          kl1 = m-1
      else
          kl1 = kl
      end if
      do i=ku1+2, ku1+kl1+1
         if ((m+ku1-i+1).ge.n) then
            k1 = n
         else
            k1 = m+ku1-i+1
         end if
         print 110, Nrow, (a(ku-ku1+i,j),j=1,k1)
         Nrow = Nrow + 1
      end do

 100  format(7x,'ARRAY ',a1,'   LDA=',i1,'  KL=',i1,'  KU=',i1)
 101  format(7x,'ARRAY ',a1,'   LDA=',i1)
 102  format(7x,i2)
 110  format(7x,i2, 2x,10(f6.2,',',f6.2,3x))
 120  format(7x,i2,18x,10(f6.2,',',f6.2,3x))
 130  format(7x,i2,34x,10(f6.2,',',f6.2,3x))
 140  format(7x,i2,50x,10(f6.2,',',f6.2,3x))
 150  format(7x,i2,66x,10(f6.2,',',f6.2,3x))
      return
      end
