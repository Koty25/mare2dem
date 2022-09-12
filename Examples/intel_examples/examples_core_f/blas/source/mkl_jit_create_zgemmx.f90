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
!      MKL_JIT_CREATE_ZGEMM  Example Program Text
!*******************************************************************************

include "mkl_blas.f90"

      program   MKL_JIT_CREATE_ZGEMM_MAIN
!

      use, intrinsic :: ISO_C_BINDING, only : c_ptr, c_int, c_intptr_t, c_funptr, c_double, c_f_procpointer
      use MKL_JIT_BLAS

      TYPE(C_PTR)              jitter
      character*1              transa, transb
      integer(KIND=c_intptr_t) status
      integer                  m, n, k
      integer                  lda, ldb, ldc
      complex(KIND=c_double)   alpha, beta
      integer                  rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc
      parameter                (rmaxa=4, cmaxa=4, rmaxb=5, cmaxb=5, rmaxc=4, cmaxc=5)
      parameter                (lda=rmaxa, ldb=rmaxb, ldc=rmaxc)
      complex(KIND=c_double)   a(rmaxa,cmaxa), b(rmaxb,cmaxb), c(rmaxc,cmaxc)
      TYPE(C_FUNPTR)           c_kernel

      PROCEDURE(zgemm_jit_kernel_t),POINTER :: f_kernel

!       External Subroutines
      external         PrintArrayZ

!
!      Executable Statements
!
      print*
      print*,'   MKL_JIT_CREATE_ZGEMM EXAMPLE PROGRAM'

!
!      Read input data
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
!
!       Print input data
      print*
      print*, '     INPUT DATA'
      print 101, m, n, k
      print 102, alpha, beta
      print 103, transa, transb
      if ((transa.eq.'N').or.(transa.eq.'n')) then
        call PrintArrayZ(0,0,m,k,a,lda,'A')
      else
        call PrintArrayZ(0,0,k,m,a,lda,'A')
      end if
      if ((transb.eq.'N').or.(transb.eq.'n')) then
        call PrintArrayZ(0,0,k,n,b,ldb,'B')
      else
        call PrintArrayZ(0,0,n,k,b,ldb,'B')
      end if
      call PrintArrayZ(0,0,m,n,c,ldc,'C')
!
!      Call MKL_JIT_CREATE_ZGEMM subroutine

      status = MKL_JIT_CREATE_ZGEMM(jitter,transa,transb,m,n,k,alpha,lda,ldb,beta,ldc)

      if ((status.eq.MKL_JIT_ERROR)) then
         print*, 'jit error'
         goto 999
      end if

      c_kernel = mkl_jit_get_zgemm_ptr(jitter)

      ! convert C function pointer to fortran procedure pointer
      CALL c_f_procpointer (c_kernel, f_kernel)

      call f_kernel(jitter, a, b, c)

      ! free memory
      status = mkl_jit_destroy(jitter)

!
      print*
      print*, '     OUTPUT DATA'
      call PrintArrayZ(1,0,m,n,c,ldc,'C')

      stop
 100  format(2(a1,1x))
 101  format(7x,'M=',i1,'  N=',i1,'  K=',i1)
 102  format(7x,'ALPHA=(',f5.2,',',f5.2,' )', &
               '  BETA=(',f5.2,',',f5.2,' )')
 103  format(7x,'TRANSA=',a1, '  TRANSB=',a1)
 999  stop 1
      end
