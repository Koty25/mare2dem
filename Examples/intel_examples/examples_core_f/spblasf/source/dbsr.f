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

!
!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Sparse BLAS
!             Fortran-77 example
!
!*******************************************************************************
!
! Example program for using Intel MKL Sparse BLAS Level 2 and 3
! for matrices represented in the block compressed sparse row storage scheme.
! The following Sparse  Blas routines are used in the example:
!          MKL_DBSRSM  MKL_DBSRSV  MKL_DBSRMM  MKL_DBSRMV
!          MKL_DBSRGEMV    MKL_DBSRSYMV  MKL_DBSRTRSV.
!
! Consider the matrix A (see Appendix 'Sparse Storage Formats for Sparse Blas
! level 2-3')
!
!                 |   1   -1      0   -3   |
!                 |  -2    5      0    0   |   |Q  W |
!   A    =        |   0    0      4    6   | = |0  R |,
!                 |   0    0      2    7   |
!
! decomposed as
!
!                      A = L + D + U,
!
!  where L is the strict lower triangle of A, U is the strictly upper triangle
!  of A, D is the main diagonal. Namely
!
!        |   0    0   0    0     0   |       |  0   -1    0   -3   0   |
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
!  The matrix A is represented in the compressed block sparse row storage scheme
!  with the help of three arrays (see Appendix 'Sparse Matrix Storage') as follows:
!
!         values = (1 -2 -1  5  0  0 -3  0  4  2  6  7)
!         columns = (1 2 1)
!         rowIndex = (1  3  4)
!
!  It should be noted that two variations of the compressed sparse row storage scheme are supported by
!  Intel MKL Sparse Blas (see 'Sparse Storage Formats for Sparse Blas level 2-3') :
!
!        1. variation simular to the NIST Sparse Blas
!        2. variation simular to DSS/PARDISO, CXML and many other libraries.
!
!  The representation of the matrix A  given above is the 2-th variation. Two integer arrays
!  pointerB and pointerE instead of the array rowIndex are used in the NIST variation of variation
!  of the block compressed sparse row format. Thus the arrays values and columns are the same for
!  the both variations. The arrays pointerB and pointerE for the matrix A are defined as follows:
!                          pointerB = (1 3)
!                          pointerE = (3 4)
!  It's easy to see that
!                    pointerB(i)= rowIndex(i) for i=1,2;
!                    pointerE(i)= rowIndex(i+1) for i=1,2.
!
!
!  The purpose of the given example is to show
!
!             1. how to call routines having interfaces suitable for the NIST's variation of the
!                block compressed sparse row format
!             2. how to form the arrays pointerB and pointerE for the NIST's variation of the
!                block compressed sparse row format using the  array rowIndex
!
!  In what follows the symbol ' means taking of transposed.
!
!  The test performs the following operations :
!
!       1. The code computes (L+D)'*S = F using MKL_DBSRMM where S is a known 5 by 2
!          matrix and then the code solves the system (L+D)'*X = F with the help of MKL_DBSRSM.
!          It's evident that X should be equal to S.
!
!       2. The code computes (U+I)'*S = F using MKL_DBSRMV where S is a vector
!          and then the code calls MKL_DBSRTRSV solves the system (U+I)'*X = F with the single right
!          hand side. It's evident that X should be equal to S.
!
!       3. The next step is the computation (U-U') S = F using MKL_DBSRMV where S is
!          a vector. It is easy to see that U-U' is a skew-symmetric matrix.
!
!       4. The next step is the computation (L+D+L') S = F using MKL_DBSRSYMV where S is
!          a vector. It is easy to see that L+D+L' is a symmetric matrix.
!
!       5. The next step is the computation A'* S = F using MKL_DBSRGEMV where S is
!          a vector.
!
!*******************************************************************************
!     Definition arrays for sparse representation of  the matrix A in
!     the compressed sparse row format:
!*******************************************************************************
          integer  m,  nnz, nnzb, lb
          parameter( m = 2,  nnz = 12, nnzb = 3, lb = 2)
          real*8  values(nnz)
          integer columns(nnzb), rowIndex(m+1)
          integer pointerB(m) , pointerE(m)
          data values/  1.0d0, -2.0d0, -1.0d0, 5.0d0,
     &                  0.0d0,  0.0d0, -3.0d0, 0.0d0,
     &                  4.0d0,  2.0d0,  6.0d0, 7.0d0 /
          data rowIndex/1, 3, 4/
          data columns/1, 2, 2/
!*******************************************************************************
!    Declaration of local variables :
!*******************************************************************************
          integer n
          parameter (n = 2)
          real*8 rhs(m*lb, n), sol(m*lb, n), temp(m*lb, n)
          data sol/1.D0, 1.D0, 1.D0, 1.D0,
     &             4.D0, 3.D0, 2.D0, 1.D0/
          real*8 alpha, beta
          data alpha/1.d0/, beta/0.d0/
          integer i, j
          print*
          print*, ' EXAMPLE PROGRAM FOR                         '
          print*, ' BLOCK COMPRESSED SPARSE ROW FORMAT ROUTINES '

!*******************************************************************************
!Task 1.    Obtain matrix-matrix multiply (L+D)' *sol --> rhs
!    and solve triangular system   (L+D)' *temp = rhs with multiple right hand sides
!    Array temp must be equal to the array sol
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_DBSRMM '
          print*, '     WITH TRIANGULAR MATRIX    '
          print 101, m, n
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input matrix '
          print 104, ((sol(i,j),j=1,n),i=1,m*lb)

           call mkl_dbsrmm('t', m, n, m, lb, alpha, 'tlnf',
     &         values, columns, rowIndex, rowindex(2), sol, m*lb,  beta,
     &         rhs,  m*lb)
          print*
          print*, '     OUTPUT DATA FOR MKL_DBSRMM '
          print*, '     WITH TRIANGULAR MATRIX     '
          print 104, ((rhs(i,j),j=1,n),i=1,m*lb)
          print 100
          print*, ' Solve triangular system with obtained '
          print*, ' right hand side                       '
          call mkl_dbsrsm('t', m, n, lb, alpha, 'tlnf',
     &           values,  columns, rowIndex, rowindex(2), rhs, m*lb,
     &           temp,  m*lb)

          print*
          print*, '     OUTPUT DATA FOR MKL_DBSRSM '
          print 104, ((temp(i,j),j=1,n),i=1,m*lb)
          print 100
!*******************************************************************************
! Task 2.    Obtain matrix-vector multiply (U+D)' *sol --> rhs
!    and solve triangular system   (U+D)' *temp = rhs with single right hand sides.
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
          print*, '     INPUT DATA FOR MKL_DBSRMV '
          print*, '     WITH TRIANGULAR MATRIX    '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i,1),i=1,m*lb)

           call mkl_dbsrmv('t', m, m, lb, alpha, 'tunf',
     &           values, columns, pointerB, pointerE, sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_DBSRMV '
          print*, '     WITH TRIANGULAR MATRIX     '
          print 105, (rhs(i,1),i=1,m*lb)
          print 100
          print*, ' Solve triangular system with obtained '
          print*, ' right hand side  '
          call mkl_dbsrtrsv('u', 't', 'n', m, lb,
     &           values, rowIndex, columns, rhs, temp)
          print*
          print*, '     OUTPUT DATA FOR MKL_DBSRTRSV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 105, (temp(i,1),i=1,m*lb)
          print 100
!*******************************************************************************
! Task 3.  Obtain matrix-vector multiply (U -U')*sol --> rhs
!    Array temp must be equal to the array sol
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_DBSRMV  '
          print*, '     WITH SKEW-SYMMETRIC MATRIX '
          print 102, alpha, beta
          print 103, 'n'
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m*lb)

           call mkl_dbsrmv('n', m, m, lb, alpha, 'au f',
     &          values, columns, rowIndex, rowIndex(2), sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_DBSRMV  '
          print*, '     WITH SKEW-SYMMETRIC MATRIX  '
          print 105, (rhs(i,1),i=1,m*lb)
          print 100
!*******************************************************************************
! Task 4.    Obtain matrix-vector multiply (L+D+L')*sol --> rhs whith the help of
!    MKL_DBSRSYMV
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_DBSRSYMV '
          print*, '     WITH SYMMETRIC MATRIX       '
          print 102, alpha, beta
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m*lb)

           call mkl_dbsrsymv('l', m, lb, values, rowIndex, columns,
     &           sol, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_DBSRSYMV '
          print*, '     WITH SYMMETRIC MATRIX        '
          print 105, (rhs(i,1),i=1,m*lb)
          print 100
!*******************************************************************************
! Task 5.   Obtain matrix-vector multiply A'*sol --> rhs whith the help of MKL_DBSRGEMV
!
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_DBSRGEMV '
          print*, '     WITH GENERAL MATRIX         '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m*lb)

           call mkl_dbsrgemv('t', m, lb, values, rowIndex, columns,
     &           sol, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_DBSRGEMV '
          print*, '     WITH GENERAL MATRIX          '
          print 105, (rhs(i,1),i=1,m*lb)
          print 100

 100      format('------------------------------------------------')
 101      format(7x,'M=',i1,'  N=',i1)
 102      format(7x,'ALPHA= ',f4.1,' BETA= ', f4.1)
 103      format(7x,'TRANS=',a1)
 104      format(2(f7.1, 3x))
 105      format(f4.1)
          stop
          end
