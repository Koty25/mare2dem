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
! for matrices represented in the skyline storage scheme.
! The following Sparse  Blas routines are used in the example:
!          MKL_CSKYSM  MKL_CSKYSV  MKL_CSKYMM  MKL_CSKYMV.
!
! Consider the matrix A 
!
!                 |   1,1       -1,1     -3,1    0     0    |
!                 |  -2,1        5,1      0      0     0    |
!   A    =        |   0          0        4,1    6,1   4,1  |,
!                 |  -4,1        0        2,1    7,1   0    |
!                 |   0          8,1      0      0    -5,1  |
!
!  
! decomposed as 
!
!                      A = L + D + U,      
!
!  where L is the strict  lower triangle of A, U is the strictly  upper triangle
!  of A, D is the main diagonal. Namely 
!   
!        |   0    0   0    0     0   |       |  0   -1,1 -3,1  0   0   | 
!        |  -2,1  0   0    0     0   |       |  0    0    0    0   0   |
!   L  = |   0    0   0    0     0   |,  U=  |  0    0    0    6,1 4,1 |
!        |  -4,1  0   2,1  0     0   |       |  0    0    0    0   0   |
!        |   0    8,1 0    0     0   |       |  0    0    0    0   0   |
!
!
!           |   1,1  0    0     0     0     |
!           |   0    5,1  0     0     0     |
!   D    =  |   0    0    4,1   0     0     |.
!           |   0    0    0     7,1   0     |
!           |   0    0    0     0    -5,1   |
!
!  
!
!  The arrays val_lower and pointers_lower represent the lower triangle of the 
!  matrix A given above (e.g. they represent the matrix L+D)
!
!
!    val_lower=( 1,1   -2,1  5,1  4,1  -4,1  0,0  2,1  7,1  8,1  0,0  0,0  -5,1)
!    pointers_lower= ( 1  2   4  5  9  13)
!
!  The arrays val_upper and pointers_upper represent the upper triangle U+D of the 
!  matrix A given above (e.g. they represent the matrix U+D)
!
!    val_upper=( 1,1   -1,1   5,1  -3,1  0,0  4,1   6,1  7,1  4,1  0,0  -5,1)
!    pointers_upper= ( 1  2   4  7  9  12)
!
!  The code given below uses only these two sparse representations for the all operations.
!
!  In what follows the symbol ' means taking of transposed. 
!
!  The test performs the following operations: 
!    
!       1. The code computes (L+D)*S = F using MKL_CSKYMM where S is a known 5 by 2  
!          matrix and then the code solves the system A*X = F with the help of MKL_CSKYSM. 
!          It's evident that X should be equal to S.
!
!       2. The code computes (U+I)*S = F using MKL_CSKYMV where S is a vector
!          and then the code solves the system A*X = F with the single right hand side. 
!          It's evident that X should be equal to S.
!
!       3. The code computes D*S = F using MKL_CSKYMV where S is a vector
!          and then the code solves the system D*X = F with the single right hand side. 
!          It's evident that X should be equal to S.
!
!       4. The next step is the computation (L+D+L') S = F using MKL_CSKYMV where S is 
!          a vector. It is easy to see that L+D+L' is a symmetric matrix.
!
!       5. The next step is the computation (U-U') S = F using MKL_CSKYMV where S is 
!          a vector. It is easy to see that U-U' is a skew-symmetric matrix.
!
!*******************************************************************************
!     Definition of lower triangle of the matrix A including the main diagonal: 
!******************************************************************************* 
          integer  m
          parameter(m=5)
          complex*8 val_lower(12)
          integer pointers_lower(6)
          data val_lower/(1.0, 1.0), (-2.0, 1.0), (5.0, 1.0),
     &      (4.0, 1.0), (-4.0, 1.0), (0.0, 1.0),  (2.0, 1.0),
     &      (7.0, 1.0),  (8.0, 1.0),  (0.0, 1.0),  (0.0, 1.0),
     &      (-5.0, 1.0)/
          data  pointers_lower/ 1,  2,   4,  5,  9,  13/   
!*******************************************************************************
!     Definition of upper triangle of the matrix A: 
!*******************************************************************************
          complex*8 val_upper(11)   
          integer pointers_upper(6)
          data val_upper/(1.0, 1.0), (-1.0, 1.0), (5.0, 1.0),
     &      (-3.0, 1.0), (0.0, 1.0), (4.0, 1.0),  
     &      (6.0, 1.0),  (7.0, 1.0), (4.0, 1.0),
     &      (0.0, 1.0),  (-5.0, 1.0)/
          data pointers_upper/ 1,   2, 4,  7,  9,  12/
!*******************************************************************************
!    Declaration of local variables : 
!*******************************************************************************
          integer n
          parameter (n=2)
          complex*8 rhs(m, n), sol(m, n), temp(m, n)   
          data sol/(1.0, 1.0), (1.0, 1.0), (1.0, 1.0),
     &      (1.0, 1.0), (1.0, 1.0), (5.0, 1.0), (4.0, 1.0),
     &      (3.0, 1.0), (2.0, 1.0), (1.0, 1.0)/
          complex*8 alpha, beta
          data alpha/(1.0, 0.0)/, beta/(0.0, 0.0)/
          integer i, j
          print*
          print*, ' EXAMPLE PROGRAM FOR SKYLINE FORMAT ROUTINES '

!*******************************************************************************
!    Task 1. Obtain matrix-matrix multiply (L+D) *sol --> rhs
!    and solve triangular system   (L+D) *temp = rhs with multiple right hand sides
!    Array temp must be equal to the array sol 
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_CSKYMM '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 101, m, n
          print 102, alpha, beta
          print 103, 'n'
          print*, ' Input matrix '
          print 104, ((sol(i,j),j=1,n),i=1,m)
         
           call mkl_cskymm('n', m, n, m, alpha, 'tln',
     &           val_lower, pointers_lower, sol, m,  beta, rhs,  m)
          print*
          print*, '     OUTPUT DATA FOR MKL_CSKYMM '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 104, ((rhs(i,j),j=1,n),i=1,m)
          print 100
          print*, ' Solve triangular system with obtained '
          print*, ' right hand side  '


           call mkl_cskysm('n', m, n, alpha, 'tln',
     &           val_lower, pointers_lower, rhs, m, temp,  m)

          print*
          print*, '     OUTPUT DATA FOR MKL_CSKYSM '
          print 104, ((temp(i,j),j=1,n),i=1,m)
          print 100
!*******************************************************************************
!    Task 2. Obtain matrix-vector multiply (U+I) *sol --> rhs
!    and solve triangular system   (U+I) *temp = rhs with single right hand sides
!    Array temp must be equal to the array sol 
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_CSKYMV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 102, alpha, beta
          print 103, 'n'
          print*, ' Input vector '
          print 105, (sol(i,1),i=1,m)
         
           call mkl_cskymv('n', m, m, alpha, 'tuu',
     &           val_upper, pointers_upper, sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_CSKYMV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
          print*, ' Solve triangular system with obtained '
          print*, ' right hand side  '


          call mkl_cskysv('n', m,  alpha, 'tuu',
     &           val_upper, pointers_upper, rhs, temp)

          print*
          print*, '     OUTPUT DATA FOR MKL_CSKYSV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 105, (temp(i,1),i=1,m)
          print 100
!*******************************************************************************
!    Task 3. Obtain matrix-vector multiply D *sol --> rhs
!    and solve triangular system   D *temp = rhs with single right hand side
!    Array temp must be equal to the array sol 
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_CSKYMV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i,2),i=1,m)
         
           call mkl_cskymv('t', m, m, alpha, 'dun',
     &           val_upper, pointers_upper, sol(1,2), beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_CSKYMV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
          print*, ' Multiply by inverse diagonal '
          print*, ' with the help of MKL_CSKYSV '


          call mkl_cskysv('c', m,  alpha, 'dun',
     &           val_upper, pointers_upper, rhs, temp)

          print*
          print*, '     OUTPUT DATA FOR MKL_CSKYSV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 105, (temp(i,1),i=1,m)
          print 100

!*******************************************************************************
!   Task 4. Obtain matrix-vector multiply (L+I+L**T)*sol --> rhs
!   
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_CSKYMV '
          print*, '     WITH SYMMETRIC MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)
         
           call mkl_cskymv('t', m, m, alpha, 'slu',
     &           val_lower, pointers_lower, sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_CSKYSV '
          print*, '     WITH SYMMETRIC MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
!*******************************************************************************
!    Task 5. Obtain matrix-vector multiply (D-D**T)*sol --> rhs
!    Array temp must be equal to the array sol 
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_CSKYMV '
          print*, '     WITH SKEW-SYMMETRIC MATRIX '
          print 102, alpha, beta
          print 103, 'n'
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)
         
           call mkl_cskymv('n', m, m, alpha, 'au',
     &           val_upper, pointers_upper, sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_CSKYMV '
          print*, '     WITH SKEW-SYMMETRIC MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100

 100      format('------------------------------------------------')
 101      format(7x,'M=',i1,'  N=',i1)
 102      format(7x,'ALPHA= ',f4.1,' BETA= ', f4.1)
 103      format(7x,'TRANS=',a1)
 104      format(2(f4.1, 3x))
 105      format(f4.1)
          stop 
          end
