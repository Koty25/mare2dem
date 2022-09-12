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
! for matrices represented in the coordinate storage scheme.
! The following Sparse  Blas routines are used in the example:
!          MKL_SCOOSM  MKL_SCOOSV  MKL_SCOOMM  MKL_SCOOMV
!          MKL_SCOOGEMV    MKL_SCOOSYMV  MKL_SCOOTRSV.   
!
! Consider the matrix A (see Appendix 'Sparse Storage Formats for Sparse Blas 
! level 2-3')
!
!                 |   1       -1     -3    0     0   |
!                 |  -2        5      0    0     0   |
!   A    =        |   0        0      4    6     4   |,
!                 |  -4        0      2    7     0   |
!                 |   0        8      0    0    -5   |
!
!  
! decomposed as 
!
!                      A = L + D + U,      
!
!  where L is the strict  lower triangle of A, U is the strictly  upper triangle
!  of A, D is the main diagonal. Namely 
!   
!        |   0    0   0    0     0   |       |  0   -1   -3    0   0   | 
!        |  -2    0   0    0     0   |       |  0    0    0    0   0   |
!   L  = |   0    0   0    0     0   |,  U=  |  0    0    0    6   4   |
!        |  -4    0   2    0     0   |       |  0    0    0    0   0   |
!        |   0    8   0    0     0   |       |  0    0    0    0   0   |
!
!
!           |   1  0  0   0   0   |
!           |   0  5  0   0   0   |
!   D    =  |   0  0  4   0   0   |.
!           |   0  0  0   7   0   |
!           |   0  0  0   0  -5   |
! 
!  The matrix A given above is represented in the coordinate storage scheme with the help of three 
!  arrays of length nnz=13 (see Appendix 'Sparse Storage Formats for Sparse Blas level 2-3'
!
!          values = (1 -1 -3 -2 5 4 6 4 -4 2 7 8 -5) 
!          rows=    (1  1  1  2 2 3 3 3  4 4 4 5 5 ) 
!          columns = (1 2  3  1 2 3 4 5  1 3 4 2 5) 
!
!  In what follows the symbol ' means taking of transposed. 
!
!  The test performs the following operations : 
!    
!       1. The code computes (L+D)'*S = F using MKL_SCOOMM where S is a known 5 by 2  
!          matrix and then the code solves the system (L+D)'*X = F with the help of MKL_SCOOSM. 
!          It's evident that X should be equal to S.
!
!       2. The code computes (U+I)'*S = F using MKL_SCOOMV where S is a vector
!          and then the code calls MKL_SCOOTRSV solves the system (U+I)'*X = F with the single right 
!          hand side. It's evident that X should be equal to S.
!
!       3. The code computes D*S = F using MKL_SCOOMV where S is a vector
!          and then the code solves the system D*X = F with the single right hand side. 
!          It's evident that X should be equal to S.
!
!       4. The next step is the computation (U-U') S = F using MKL_SCOOMV where S is 
!          a vector. It is easy to see that U-U' is a skew-symmetric matrix.
!
!       5. The next step is the computation (L+D+L') S = F using MKL_SCOOSYMV where S is 
!          a vector. It is easy to see that L+D+L' is a symmetric matrix.
!
!       6. The next step is the computation A'* S = F using MKL_SCOOGEMV where S is 
!          a vector.
!
! The code given below uses only one sparse representation for the all operations.
!*******************************************************************************
!     Definition arrays for sparse representation of  the matrix A in 
!     the coordinate format: 
!******************************************************************************* 
          integer  m, nnz
          parameter( m = 5,  nnz=13)
          real*4 values(nnz)
          integer  rows(nnz), columns(nnz)
          data values/1.0, -1.0, -3.0, -2.0, 5.0,
     &        4.0, 6 .0, 4.0, -4.0, 2.0, 7.0, 8.0, -5.0/ 
          data rows/1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5/
          data columns/1, 2, 3, 1, 2, 3, 4, 5, 1, 3, 4, 2, 5/
!*******************************************************************************
!    Declaration of local variables : 
!*******************************************************************************
          integer n
          parameter (n=2)
          real*4 rhs(m, n), sol(m, n), temp(m, n) 
          data sol/1.0, 1.0, 1.0, 1.0, 1.0, 
     &    5.0, 4.0, 3.0, 2.0, 1.0/
          real*4 alpha, beta
          data alpha/1.0/, beta/0.0/
          integer i, j
          print*
          print*, ' EXAMPLE PROGRAM FOR COORDINATE FORMAT ROUTINES ' 

!*******************************************************************************
!    Task 1.    Obtain matrix-matrix multiply (L+D)' *sol --> rhs
!    and solve triangular system   (L+D)' *temp = rhs with multiple right hand sides
!    Array temp must be equal to the array sol 
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_SCOOMM '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 101, m, n
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input matrix '
          print 104, ((sol(i,j),j=1,n),i=1,m)
         
           call mkl_scoomm('t', m, n, m, alpha, 'tln',
     &           values, rows, columns, nnz, sol, m,  beta, rhs,  m)
          print*
          print*, '     OUTPUT DATA FOR MKL_SCOOMM '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 104, ((rhs(i,j),j=1,n),i=1,m)
          print 100
          print*, ' Solve triangular system with obtained '
          print*, ' right hand side  '
          call mkl_scoosm('t', m, n, alpha, 'tln',
     &           values, rows, columns, nnz, rhs, m, temp,  m)

          print*
          print*, '     OUTPUT DATA FOR MKL_SCOOSM '
          print 104, ((temp(i,j),j=1,n),i=1,m)
          print 100
!*******************************************************************************
!    Task 2.    Obtain matrix-vector multiply (U+I)' *sol --> rhs
!    and solve triangular system   (U+I)' *temp = rhs with single right hand sides
!    Array temp must be equal to the array sol 
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_SCOOMV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i,1),i=1,m)
         
           call mkl_scoomv('t', m, m, alpha, 'tuu',
     &           values, rows, columns, nnz, sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_SCOOMV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
          print*, ' Solve triangular system with obtained '
          print*, ' right hand side  '
          call mkl_scootrsv('u', 't', 'u', m, 
     &           values, rows, columns, nnz, rhs, temp)
          print*
          print*, '     OUTPUT DATA FOR MKL_SCOOTRSV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 105, (temp(i,1),i=1,m)
          print 100
!*******************************************************************************
!    Task 3.  Obtain matrix-vector multiply D *sol --> rhs
!    and solve triangular system   D *temp = rhs with single right hand side
!    Array temp must be equal to the array sol 
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_SCOOMV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i,2),i=1,m)
         
           call mkl_scoomv('t', m, m, alpha, 'dun',
     &          values, rows, columns, nnz,  sol(1,2), beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_SCOOMV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
          print*, ' Multiply by inverse matrix '
          print*, ' with the help of MKL_SCOOSV '


          call mkl_scoosv('t', m,  alpha, 'dun',
     &         values, rows, columns, nnz,  rhs, temp)

          print*
          print*, '     OUTPUT DATA FOR MKL_SCOOSV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 105, (temp(i,1),i=1,m)
          print 100

!*******************************************************************************
!    Task 4.  Obtain matrix-vector multiply (U -U')*sol --> rhs
!    Array temp must be equal to the array sol 
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_SCOOMV '
          print*, '     WITH SKEW-SYMMETRIC MATRIX '
          print 102, alpha, beta
          print 103, 'n'
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)
         
           call mkl_scoomv('n', m, m, alpha, 'au',
     &          values, rows, columns, nnz, sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_SCOOMV '
          print*, '     WITH SKEW-SYMMETRIC MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
!*******************************************************************************
!    Task 5.  Obtain matrix-vector multiply (L+D+L')*sol --> rhs whith the help of 
!    MKL_SCOOSYMV
!   
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_SCOOSYMV '
          print*, '     WITH SYMMETRIC MATRIX  '
          print 102, alpha, beta
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)
         
           call mkl_scoosymv('l', m,  values, rows, columns, nnz,
     &           sol, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_SCOOSYMV '
          print*, '     WITH SYMMETRIC MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
!*******************************************************************************
!    Task 6. Obtain matrix-vector multiply A'*sol --> rhs whith the help 
!    of MKL_SCOOGEMV
!   
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_SCOOGEMV '
          print*, '     WITH GENERAL MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)
         
           call mkl_scoogemv('t', m,  values, rows, columns, nnz,
     &           sol, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_SCOOGEMV '
          print*, '     WITH GENERAL MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100

 100      format('------------------------------------------------')
 101      format(7x,'M=',i1,'  N=',i1)
 102      format(7x,'ALPHA= ',f4.1,' BETA= ', f4.1)
 103      format(7x,'TRANS=',a1)
 104      format(2(f7.1, 3x))
 105      format(f4.1)
          stop 
          end
