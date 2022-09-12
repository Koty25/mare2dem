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
*             Fortran example
*
********************************************************************************
C----------------------------------------------------------------------
C Example program to show the use of the Intel MKL Sparse BLAS routines
C for multiplication of two compressed sparse row format matrices
C---------------------------------------------------------------------
      PROGRAM multilication_test
      IMPLICIT NONE
C..  Include interfaces for Intel MKL service functions
      INCLUDE "mkl_service.fi"
C..  Description of all variables
      CHARACTER trans
      INTEGER n, sort, job, ierr, nzmax
      INTEGER ia(9)
      INTEGER ja(20)
      COMPLEX*16 a(20), cden(8, 8), answer(8, 8) 
      COMPLEX*16 zero
      COMPLEX*16 ddum
      REAL*8 tolerance, normA, normres
      INTEGER ic, jc
      COMPLEX*16 c
      POINTER (IC_PTR, IC(1)), (JC_PTR, JC(1)), (C_PTR, C(1))
      INTEGER i, j ,js, idum, nnz, sizeint
      REAL*8 DZNRM2
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
      DATA zero/(0.D0, 0.D0)/
      DATA ia /1,5,8,10,12,13,16,18,21/

      DATA ja
     1        /1,          3,                 6,    7,
     2               2,    3,          5,                   
     3                     3,                             8,            
     4                          4,                  7,     
     5               2,                                         
     6                     3,                 6,          8,         
     7               2,                             7,      
     8                     3,                       7,    8/
      DATA a
     1     /(7.d0, 1.d0), (1.d0,1.d0), (2.d0,1.d0), (7.d0,1.d0),
     2      (-4.d0,0.d0), (8.d0,1.d0), (2.d0,1.d0),                   
     3      (1.d0,1.d0),  (5.d0,1.d0),            
     4      (7.d0,0.d0),  (9.d0,1.d0),     
     5      (-4d0,1.d0),       
     6      (7.d0,1.d0),  (3.d0,1.d0), (8.d0,0.d0),         
     7      (1.d0,1.d0),  (11.d0,1.d0),      
     8      (-3.d0,1.d0), (2.d0,1.d0), (5.d0,0.d0)/ 
C..
C.. Compute A^T*A for the given matrix A
C..
      trans='t'
      job=1
      tolerance=1.E-16
      normA=DZNRM2(ia(n+1)-1, a, 1)
      tolerance=2.D0*dsqrt(dble(n))*tolerance*normA 
c ..
c ..   Compute the result matrix with the help of mkl_zcsrmultd
c ..
      call mkl_zcsrmultd(trans, n, n, n, a,ja,ia,
     *          a, ja, ia, answer, n)
c ..  
c ..  Compute the same matrix as a compressed sparse row matrix and 
c ..  allocate all arrays for the result matrix in advance
c ..  
      alloc_size = sizeint*(n+1)
      ic_ptr=mkl_malloc(alloc_size, alignment) 
      if(ic_ptr.eq.0) then
         write(*,*) 'Cannot allocate pointer array of the length ',n+1
         stop 1
      endif
      nzmax=n*n
      alloc_size = sizeint*nzmax
      jc_ptr=mkl_malloc(alloc_size, alignment) 
      if(jc_ptr.eq.0) then
         write(*,*) 'Cannot allocate integer array of the length ',nnz
         stop 1
      endif
      alloc_size = 16*nzmax
      c_ptr=mkl_malloc(alloc_size, alignment) 
      if(c_ptr.eq.0) then
         write(*,*) 'Cannot allocate value array of the length ', nnz
         stop 1
      endif
      job=0
      call mkl_zcsrmultcsr(trans, job, sort, n, n, n, a,ja,ia,
     *        a, ja, ia, c, jc, ic, nzmax, ierr)
c ..
c ..  Check the result of computations
c ..
      if(ierr.ne.0)then 
         print  *,' FIRST TEST FAILED  ', ierr
         stop 1
      endif
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
      normres=DZNRM2(n*n, cden, 1)
      call mkl_free(jc_ptr)
      call mkl_free(c_ptr)    
c 
      if(normres.le.tolerance)then
         print  *,' FIRST TEST PASSED '
      else
         print  *,' FIRST TEST FAILED  '
         print  *,' INCORRECT RESULT ', normres, tolerance
         stop 1
      endif 

c ..
c .. In order to save memory for the result matrix the multiplication can be done
c .. in two steps. The first step is to compute the size of column indices array
c .. Allocate only pointer array and call the routine
c .. to get the actual size of arrays jc and c  
c ..
      do i=1, n+1
        ic(i)=0
      enddo
      job=1
      call mkl_zcsrmultcsr(trans, job, sort, n, n, n, a,ja,ia,
     *        a, ja,ia, ddum, idum, ic, idum, idum)
      nnz=ic(n+1)-1
      alloc_size = sizeint*nnz
      jc_ptr=mkl_malloc(alloc_size, alignment) 
      if(jc_ptr.eq.0) then
	 write(*,*) 'Cannot allocate integer array of the length ',nnz
         stop 1
      endif
      alloc_size = 16*nnz
      c_ptr=mkl_malloc(alloc_size, alignment) 
      if(c_ptr.eq.0) then
         write(*,*) 'Cannot allocate value array of the length ', nnz
         stop 1
      endif
c ..
c .. The next step is to compute the result matrix
c .. 
      job=2

      call mkl_zcsrmultcsr(trans, job, sort, n, n, n, a, ja, ia,
     *        a, ja, ia, c, jc, ic, idum, idum) 
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
      normres=DZNRM2(n*n, cden, 1)
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
      stop
      end
