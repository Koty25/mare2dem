!===============================================================================
! Copyright 2004-2020 Intel Corporation.
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

*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Sparse BLAS
*             Fortran-77 example
*
********************************************************************************
C----------------------------------------------------------------------
C Example program to show the use of the Intel MKL Sparse BLAS routine
C for addition of two compressed sparse row format
C---------------------------------------------------------------------
      PROGRAM addition_test
      IMPLICIT NONE
C..  Include interfaces for Intel MKL service functions
      INCLUDE "mkl_service.fi"
C..  Description of all variables
      CHARACTER trans
      INTEGER n, sort, job, ierr, nzmax
      INTEGER ia(9)
      INTEGER ja(18)
      REAL a(18), cden(8, 8), answer(8, 8) 
      REAL  zero
      REAL beta, ddum, tolerance, normA, normres
      INTEGER i, j ,js, idum, nnz, sizeint
      INTEGER ic, jc
      REAL c
      POINTER (IC_PTR, IC(1)), (JC_PTR, JC(1)), (C_PTR, C(1)) 
      REAL  SNRM2
      INTEGER*4 alignment
#ifdef _IA32
      INTEGER*4 alloc_size
#else
      INTEGER*8 alloc_size
#endif  
c ..  
c .. Fill all arrays containing matrix data.
c ..
      DATA n /8/, sort/0/, sizeint/8/
      DATA beta/1.0/, zero/0.0/
      DATA ia /1,5,8,10,12,15,17,18,19/
      DATA ja
     1 /1,  3,    6,7,
     2    2,3,  5,
     3      3,        8,
     4        4,    7,
     5          5,6,7,
     6            6,  8,
     7              7,
     8                8/
      DATA a
     1 /7.0,     1.0,          2.0,7.0,
     2       -4.0,8.0,     2.0,
     3            1.0,                    5.0,
     4                 7.0,     9.0,
     5                      5.0,1.0,5.0,
     6                           -1.0,     5.0,
     7                                11.0,
     8                                     5.0/
C..
C.. Compute A+A^T for the given matrix A
C..
      trans='t'
      tolerance= 1.D-8
c ..
c ..   Compute the correct answer directly
c ..
      do i=1, n
	do j=1, n
	  answer(i,j)=zero
        enddo
      enddo	
      do i=1, n
	do j=ia(i), ia(i+1)-1
	  js=ja(j)
	  answer(i, js)=answer(i, js)+ a(j)
          answer(js, i)=answer(js, i)+ a(j)
        enddo
      enddo
      normA=snrm2(ia(n+1)-1, a, 1)
      tolerance=2.D0*sqrt(dble(n))*tolerance*normA 
c ..
c .. Allocate pointer array and call the routine
c .. to get the actual size of arrays jc and c  
c ..
      alloc_size = sizeint*(n+1)
      ic_ptr=mkl_malloc(alloc_size, alignment) 
      if(ic_ptr.eq.0) then
         write(*,*) 'Cannot allocate pointer array of the length ',n+1
         stop 1
      endif
      job=1 
      call mkl_scsradd(trans, job, sort, n, n, a,ja,ia,
     *         beta, a, ja,ia, ddum, idum, ic, idum, idum) 
      nnz=ic(n+1)-1
      alloc_size = sizeint*nnz
      jc_ptr=mkl_malloc(alloc_size, alignment) 
      if(jc_ptr.eq.0) then
	 write(*,*) 'Cannot allocate integer array of the length ',nnz
         stop 1
      endif
      alloc_size = 4*nnz
      c_ptr=mkl_malloc(alloc_size, alignment) 
      if(c_ptr.eq.0) then
         write(*,*) 'Cannot allocate value array of the length ', nnz
         stop 1
      endif 
      job=2
      call mkl_scsradd(trans, job, sort, n, n, a,ja, ia,
     *         beta, a, ja,ia, c, jc, ic, nzmax, ierr) 
c ..
c     Check the correctness of the result
c ..
      do i=1, n
        do j=1, n
	  cden(i,j)=zero
	enddo
      enddo
      do i=1, n
        do j=ic(i), ic(i+1)-1
	  js=jc(j)
	  cden(i, js)=answer(i, js)- c(j)
	enddo
      enddo
      normres=snrm2(n*n, cden, 1)
      call mkl_free(ic_ptr)
      call mkl_free(jc_ptr)
      call mkl_free(c_ptr)   
      if(normres.le.tolerance)then
        print  *,' FIRST TEST PASSED '
      else 
         print  *,' FIRST TEST FAILED  '
         stop 1
      endif
c ..
c If the maximum number of non-zeros in the result matrix is known in advance, 
c the result matrix can be computed in the following way
c ..
      alloc_size=sizeint*(n+1)
      ic_ptr=mkl_malloc(alloc_size,alignment) 
      if(ic_ptr.eq.0) then
         write(*,*) 'Cannot allocate pointer array of the length ',n+1
         stop 1
      endif
      nzmax=n*n
      alloc_size =sizeint*nzmax
      jc_ptr=mkl_malloc(alloc_size, alignment) 
      if(jc_ptr.eq.0) then
         write(*,*) 'Cannot allocate integer array of the length ',nnz
         stop 1
      endif
      alloc_size =4*nzmax
      c_ptr=mkl_malloc(alloc_size, alignment) 
      if(c_ptr.eq.0) then
         write(*,*) 'Cannot allocate value array of the length ', nnz
         stop 1
      endif 
      job=0
      call mkl_scsradd(trans, job, sort, n, n, a,ja, ia,
     *         beta, a, ja,ia, c, jc, ic, nzmax, ierr) 
      if(ierr.ne.0)then 
	   print  *,' SECOND TEST FAILED  '
           stop 1
      endif
c ..
c     Check the correctness of the result
c ..
      do i=1, n
        do j=1, n
	  cden(i,j)=zero
	enddo
      enddo	
      do i=1, n
        do j=ic(i), ic(i+1)-1
	   js=jc(j)
	   cden(i, js)=answer(i, js)- c(j)
        enddo
      enddo
      normres=snrm2(n*n, cden, 1)
      call mkl_free(ic_ptr)
      call mkl_free(jc_ptr)
      call mkl_free(c_ptr)   
      if(normres.le.tolerance)then
         print  *,' SECOND TEST PASSED '
         stop 0
      else 
	 print  *,' SECOND TEST FAILED  '
         stop 1
      endif 
      end
