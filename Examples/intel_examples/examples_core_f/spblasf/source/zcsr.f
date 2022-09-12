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

!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Sparse BLAS
!             Fortran-77 example
!
!*******************************************************************************
!
! Example program for using Intel MKL Sparse BLAS Level 2 and 3
! for matrices represented in the compressed sparse row storage scheme.
! The following Sparse  Blas routines are used in the example:
!          MKL_ZCSRSM  MKL_ZCSRSV  MKL_ZCSRMM  MKL_ZCSRMV
!          MKL_ZCSRGEMV    MKL_ZCSRSYMV  MKL_ZCSRTRSV.   
!
! Consider the matrix A 
!
!                 |   1,1     -1,1    0   -3,1   0   |
!                 |  -2,1      5,1    0    0     0   |
!   A    =        |   0        0      4,1  6,1   4,1 |,
!                 |  -4,1      0      2,1  7,1   0   |
!                 |   0        8,1    0    0    -5,1 |
!
!  
! decomposed as 
!
!                      A = L + D + U,      
!
!  where L is the strict  lower triangle of A, U is the strictly  upper triangle
!  of A, D is the main diagonal. Namely 
!   
!        |   0    0    0    0     0   |       |  0   -1,1  0   -3,1  0   | 
!        |  -2,1  0    0    0     0   |       |  0    0    0    0    0   |
!   L  = |   0    0    0    0     0   |,  U=  |  0    0    0    6,1  4,1 |
!        |  -4,1  0    2,1  0     0   |       |  0    0    0    0    0   |
!        |   0    8,1  0    0     0   |       |  0    0    0    0    0   |
!
!
!           |   1,1  0    0    0    0   |
!           |   0    5,1  0    0    0   |
!   D    =  |   0    0    4,1  0    0   |.
!           |   0    0    0    7,1  0   |
!           |   0    0    0    0   -5   |
! 
!  The matrix A is represented in the compressed sparse row storage scheme with the help of three 
!  arrays  (see Appendix 'Sparse Matrix Storage') as follows:
!
!         values = (1,1 -1,1 -3,1 -2,1 5,1 4,1 6,1 4,1 -4,1 2,1 7,1 8,1 -5,1) 
!         columns = (1 2 4 1 2 3 4 5 1 3 4 2 5) 
!         rowIndex = (1  4  6  9  12 14)
!
!  It should be noted that two variations of the compressed sparse row storage scheme are supported by 
!  Intel MKL Sparse Blas (see 'Sparse Storage Formats for Sparse Blas level 2-3') :
!
!        1. variation accepted in the NIST Sparse Blas
!        2. variation accepted for DSS/PARDISO, CXML and many other libraries.  
!
!  The representation of the matrix A  given above is the PARDISO's variation. Two integer arrays 
!  pointerB and pointerE instead of the array rowIndex are used in the NIST variation of variation 
!  of the compressed sparse row format. Thus the arrays values and columns are the same for the both 
!  variations. The arrays pointerB and pointerE for the matrix A are defined as follows:
!                          pointerB = (1 4  6  9 12)
!                          pointerE = (4 6  9 12 14) 
!  It's easy to see that 
!                    pointerB(i)= rowIndex(i) for i=1, ..5;
!                    pointerE(i)= rowIndex(i+1) for i=1, ..5.
!
!
!  The purpose of the given example is to show
!  
!             1. how to call routines having interfaces suitable for the NIST's variation of the 
!                compressed sparse row format
!             2. how to form the arrays pointerB and pointerE for the NIST's variation of the 
!                compressed sparse row format using the  array rowIndex 
!             3. how to use minors of the matrix A by redefining the arrays pointerB and pointerE
!                but the arrays values and columns are the same.  
!
!  In what follows the symbol ' means taking of transposed. 
!
!  The test performs the following operations : 
!    
!       1. The code computes (L+D)'*S = F using MKL_ZCSRMM where S is a known 5 by 2  
!          matrix and then the code solves the system (L+D)'*X = F with the help of MKL_ZCSRSM. 
!          It's evident that X should be equal to S.
!
!       2. The code computes (U+I)'*S = F using MKL_ZCSRMV where S is a vector
!          and then the code calls MKL_ZCSRTRSV solves the system (U+I)'*X = F with the single right 
!          hand side. It's evident that X should be equal to S.
!
!       3. The code computes D*S = F using MKL_ZCSRMV where S is a vector
!          and then the code solves the system D*X = F with the single right hand side. 
!          It's evident that X should be equal to S.
!
!       4. The next step is the computation (U-U') S = F using MKL_ZCSRMV where S is 
!          a vector. It is easy to see that U-U' is a skew-symmetric matrix.
!
!       5. The next step is the computation (L+D+L') S = F using MKL_ZCSRSYMV where S is 
!          a vector. The vector is computed two times. At first, the sparse representation
!          of the whole matrix A is used. Then the vector is computed with the help of
!          sparse representation of L+D. These two calls must give the same vector.
!
!       6. The next step is the computation A'* S = F using MKL_ZCSRGEMV where S is 
!          a vector.
!     
!       7. Let's T be the upper 3 by 3 minor of the matrix A. Namely, T is the following matrix  
!
!                        |   1,1     -1,1    0   |
!          T    =        |  -2,1      5,1    0   |.
!                        |   0        0      4,1 |
!          The test performs the matrix-vector multiply T*S=F with the same arrays values, columns 
!          and pointerB used before for the whole matrix A. It is enough to change two values of 
!	   array pointerE in order to use the minor under consideration. Then the test solves the system
!          T*X =F using MKL_ZCSRSV. The routine MKL_ZCSRMV is used for getting matrix-vector multiply.    
!
! The code given below uses only one sparse representation for the all operations.
!
!*******************************************************************************
!    Declaration of arrays for sparse representation of  the matrix A
!    and the lower triangle of A in the compressed sparse row format:
!
!******************************************************************************* 
          implicit none
          integer  m,  nnz, mnew, nnz1
          parameter( m = 5,  nnz=13, mnew=3, nnz1=9)
          complex*16  values(nnz), values1(nnz1)
          integer columns(nnz), rowIndex(m+1), columns1(nnz1),
     &       rowIndex1(m+1)
          integer pointerB(m) , pointerE(m)
!*******************************************************************************
!    Sparse representation of the matrix A
!*******************************************************************************
          data values/(1.0d0, 1.0d0), (-1.0d0, 1.0d0),
     &      (-3.0d0, 1.0d0), (-2.0d0, 1.0d0), (5.0d0, 1.0d0),
     &      (4.0d0, 1.0d0), (6.0d0, 1.0d0), (4.0d0, 1.0d0),
     &      (-4.0d0, 1.0d0), (2.0d0, 1.0d0), (7.0d0, 1.0d0),
     &      (8.0d0, 1.0d0), (-5.0d0, 1.0d0)/ 
          data rowIndex/1, 4,  6,  9,  12, 14/
          data columns/1, 2, 4, 1, 2, 3, 4, 5, 1, 3, 4, 2, 5/
!*******************************************************************************
!    Sparse representation of the lower triangle L+D
!*******************************************************************************
          data values1/(1.d0, 1.d0),(-2.d0, 1.d0),
     &       (5.d0, 1.d0), (4.d0, 1.d0),
     &     (-4.d0, 1.d0), (2.d0, 1.d0), (7.d0, 1.d0),
     &     (8.d0,1.d0), (-5.d0, 1.d0)/
          data columns1/1,  1, 2, 3, 1, 3, 4, 2, 5/
          data rowIndex1/1, 2, 4, 5, 8, 10/
 
!*******************************************************************************
!    Declaration of local variables : 
!*******************************************************************************
          integer n
          parameter (n=2)
          complex*16 rhs(m, n), sol(m, n), temp(m, n) 
          data sol/(1.0d0, 1.0d0), (1.0d0, 1.0d0), (1.0d0, 1.0d0),
     &      (1.0d0, 1.0d0), (1.0d0, 1.0d0), (5.0d0, 1.0d0),
     &      (4.0d0, 1.0d0), (3.0d0, 1.0d0), (2.0d0, 1.0d0),
     &      (1.0d0, 1.0d0)/
          complex*16 alpha, beta
          data alpha/(1.0d0, 0.0d0)/, beta/(0.0d0, 0.0d0)/
          integer i, j
          print*
          print*, ' EXAMPLE PROGRAM FOR COMPRESSED SPARSE ROW
     &       FORMAT ROUTINES ' 

!*******************************************************************************
!Task 1.    Obtain matrix-matrix multiply (L+D)' *sol --> rhs
!    and solve triangular system   (L+D)' *temp = rhs with multiple right hand sides
!    Array temp must be equal to the array sol 
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_ZCSRMM '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 101, m, n
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input matrix '
          print 104, ((sol(i,j),j=1,n),i=1,m)
         
           call mkl_zcsrmm('t', m, n, m, alpha, 'tln',
     &       values, columns, rowIndex, rowindex(2), sol, m,  beta,
     &       rhs,  m)
          print*
          print*, '     OUTPUT DATA FOR MKL_ZCSRMM '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 104, ((rhs(i,j),j=1,n),i=1,m)
          print 100
          print*, ' Solve triangular system with obtained '
          print*, ' right hand side  '
          call mkl_zcsrsm('t', m, n, alpha, 'tln',
     &         values,  columns, rowIndex, rowindex(2), rhs, m, temp,
     &         m)

          print*
          print*, '     OUTPUT DATA FOR MKL_ZCSRSM '
          print 104, ((temp(i,j),j=1,n),i=1,m)
          print 100
!*******************************************************************************
! Task 2.    Obtain matrix-vector multiply (U+I)' *sol --> rhs
!    and solve triangular system   (U+I)' *temp = rhs with single right hand sides.
!    Array temp must be equal to the array sol. 
!
!    Let us form the arrays pointerB and pointerE for the NIST's variation of the 
!    compressed sparse row format using the  array rowIndex. 
!
!*******************************************************************************
          do i=1, m
            pointerB(i)=rowIndex(i)
            pointerE(i)=rowIndex(i+1)
          enddo 
          print*
          print*, '     INPUT DATA FOR MKL_ZCSRMV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i,1),i=1,m)
         
           call mkl_zcsrmv('t', m, m, alpha, 'tuu',
     &           values, columns, pointerB, pointerE, sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_ZCSRMV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
          print*, ' Solve triangular system with obtained '
          print*, ' right hand side  '
          call mkl_zcsrtrsv('u', 't', 'u', m, 
     &           values, rowIndex, columns, rhs, temp)
          print*
          print*, '     OUTPUT DATA FOR MKL_ZCSRTRSV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 105, (temp(i,1),i=1,m)
          print 100
!*******************************************************************************
! Task 3.   Obtain matrix-vector multiply D *sol --> rhs
!    and solve triangular system   D *temp = rhs with single right hand side
!    Array temp must be equal to the array sol 
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_ZCSRMV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i,2),i=1,m)
         
           call mkl_zcsrmv('n', m, m, alpha, 'dun',
     &          values, columns, pointerB, pointerE,  sol(1,2), beta,
     &          rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_ZCSRMV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
          print*, ' Multiply by inverse matrix '
          print*, ' with the help of MKL_ZCSRSV '


          call mkl_zcsrsv('t', m, alpha, 'dun',
     &         values, columns, rowIndex, rowIndex(2),  rhs, temp)

          print*
          print*, '     OUTPUT DATA FOR MKL_ZCSRSV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 105, (temp(i,1),i=1,m)
          print 100

!*******************************************************************************
! Task 4.  Obtain matrix-vector multiply (U -U')*sol --> rhs
!    Array temp must be equal to the array sol 
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_ZCSRMV '
          print*, '     WITH SKEW-SYMMETRIC MATRIX '
          print 102, alpha, beta
          print 103, 'n'
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)
         
           call mkl_zcsrmv('n', m, m, alpha, 'au',
     &          values, columns, rowIndex, rowIndex(2), sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_ZCSRMV '
          print*, '     WITH SKEW-SYMMETRIC MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
!*******************************************************************************
! Task 5.    Obtain matrix-vector multiply (L+D+L')*sol --> rhs with the help of 
!    MKL_ZCSRSYMV
! NOTE: The routine mkl_zcsrsymv as well as the similar dense Level 2
! routine zsymv has a possibilty to extract the required triangle from the input
! sparse matrix and perform symmetric matrix-vector multiply with the help of
! the triangle. Let the arrays values, rowIndex, columns be the sparse
! representation of A
!
!                |   (1,1)   (-1,1)    0   (-3,1)   0   |
!                |  (-2,1)    (5,1)    0      0     0   |
!  A    =        |     0        0    (4,1)  (6,1) (4,1) |,
!                |  (-4,1)      0    (2,1)  (7,1)   0   |
!                |     0      (8,1)    0      0  (-5,1) |
!
! Let the arrays values1, rowIndex1, columns1 which are the following
!     data values1/(1.0, 1.0),(-2.0, 1.0), (5.0, 1.0), (4.0, 1.0),
!     &        (-4.0, 1.0), (2.0, 1.0), (7.0, 1.0),  (8.0,1.0), (-5.0, 1.0)/
!   data column1/1,  1, 2, 3, 1, 3, 4, 2, 5/
!   data rowIndex/1, 2, 4, 5, 8, 10/
!
! be the sparse representation of the lower triangle of A
!
!
!           |   (1,1)      0      0     0      0   |
!           |  (-2,1)    (5,1)    0     0      0   |
!   L +D  = |     0        0    (4,1)   0      0   |,
!           |  (-4,1)      0    (2,1) (7,1)    0   |
!           |     0      (8,1)    0     0    (-5,1)|
!
!  The feature described above means that  the following two calls must give the
!  same output vector
!  call mkl_zcsrsymv(uplo, m, values, rowIndex, columns, sol_vec, rhs_vec);
!  call mkl_zcsrsymv(uplo, m, values1, rowIndex1, columns1, sol_vec, temp);
!
!  The test checks whether these two calls give the same output vector.
!
!*******************************************************************************
          print*
          print*, '     SYMMETRIC MATRIX-VECTOR MULTIPLY ROUTINE '
          print*, '     INPUT DATA FOR MKL_ZCSRSYMV '
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)
           call mkl_zcsrsymv('l', m,  values, rowIndex, columns,
     &           sol, rhs)
           call mkl_zcsrsymv('l', m,  values1, rowIndex1, columns1,
     &           sol, temp)
          print*
          print*, '     SYMMETRIC MATRIX-VECTOR MULTIPLY '
          print*, '     MKL_ZCSRSYMV CALLED TWO TIMES '
          do i=1, m
              print 104, rhs(i,1), temp(i,1)
          enddo
          print 100
!*******************************************************************************
! Task 6.   Obtain matrix-vector multiply A'*sol --> rhs with the help of MKL_ZCSRGEMV
!   
!*******************************************************************************
          print*
          print*, '   MATRIX-VECTOR MULTIPLY ROUTINE FOR GENERAL
     &       MATRICES '
          print*, '     INPUT DATA FOR MKL_ZCSRGEMV '
          print* , ' MATRIX TRANSPOSED '
 
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)
         
           call mkl_zcsrgemv('t', m,  values, rowIndex, columns,
     &           sol, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_ZCSRGEMV '
          print*, '     MATRIX-VECTOR MULTIPLY ROUTINE FOR
     &       GENERAL MATRICES  ' 
          print 105, (rhs(i,1),i=1,m)
          print 100

!*******************************************************************************
! Task 7.  Obtain matrix-vector multiply T*sol --> rhs with the help of MKL_ZCSRMV
!    where S is  3 by 3 minor of the matrix A starting with A(1,1) 
!    Let's us redefine two elements of the array pointerE in order to identify 
!    the needed minor. More precisely
!            pointerE(1) --> pointerE(1)-1
!            pointerE(3) --> pointerE(3)-2
!   
!*******************************************************************************
          pointerE(1) = pointerE(1) - 1
          pointerE(3) = pointerE(3) - 2
          print*
          print*, '     INPUT DATA FOR MKL_ZCSRMV '
          print*, '     WITH A MINOR OF GENERAL MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,mnew)
           call mkl_zcsrmv('n', mnew, mnew, alpha, 'tln', values, 
     &       columns,   pointerB, pointerE,  sol, beta, rhs)
         
          print*
          print*, '     OUTPUT DATA FOR MKL_ZCSRGEMV '
          print*, '     WITH A MINOR OF GENERAL MATRIX  '
          print 105, (rhs(i,1),i=1,mnew)
          print 100

          print*, ' Multiply by inverse to a minor of the matrix '
          print*, ' with the help of MKL_ZCSRSV '


          call mkl_zcsrsv('n', mnew, alpha, 'tln', values, 
     &        columns, pointerB, pointerE,  rhs, temp)

          print*
          print*, '     OUTPUT DATA FOR MKL_ZCSRSV '
          print*, '     WITH A MINOR OF GENERAL MATRIX '
          print 105,   (temp(i,1),i=1,mnew)
          print 100
 100      format('------------------------------------------------')
 101      format(7x,'M=',i1,'  N=',i1)
 102      format(7x,'ALPHA= (',f4.1,', ', f4.1, '),',' BETA= (' 
     &      ,f4.1,', ', f4.1, ')' ) 
 103      format(7x,'TRANS=',a1)
 104      format('(',f7.1,', ', f7.1, ')',4x, '(',f7.1,', ', f7.1, ')')
 105      format('(',f4.1,', ', f4.1, ')') 
          stop 
          end
