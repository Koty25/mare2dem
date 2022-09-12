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
!          MKL_DCSCSM  MKL_DCSCSV  MKL_DCSCMM  MKL_DCSCMV.
!
! Consider the matrix A (see Appendix 'Sparse Storage Formats for Sparse Blas
! level 2-3')
!
!                 |   1       -1      0   -3     0   |
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
!  The matrix A is represented in the compressed sparse row storage scheme with the help of three
!  arrays  (see Appendix 'Sparse Matrix Storage') as follows:
!
!         values = (1 -1 -3 -2 5 4 6 4 -4 2 7 8 -5)
!         rows = (1 2 4 1 2 3 4 5 1 3 4 2 5)
!         colIndex = (1  4  6  9  12 14)
!
!  It should be noted that two variations of the compressed sparse row storage scheme are supported by
!  Intel MKL Sparse Blas (see 'Sparse Storage Formats for Sparse Blas level 2-3') :
!
!        1. variation accepted in the NIST Sparse Blas
!        2. variation is the Harwell-Boeing sparse matrix format.
!
!  The representation of the matrix A  given above is the PARDISO's variation. Two integer arrays
!  pointerB and pointerE instead of the array colIndex are used in the NIST variation of variation
!  of the compressed sparse row format. Thus the arrays values and rows are the same for the both
!  variations. The arrays pointerB and pointerE for the matrix A are defined as follows:
!                          pointerB = (1 4  6  9 12)
!                          pointerE = (4 6  9 12 14)
!  It's easy to see that
!                    pointerB(i)= colIndex(i) for i=1, ..5;
!                    pointerE(i)= colIndex(i+1) for i=1, ..5.
!
!
!  The purpose of the given example is to show
!
!             1. how to call routines having interfaces suitable for the NIST's variation of the
!                compressed sparse row format
!             2. how to form the arrays pointerB and pointerE for the NIST's variation of the
!                compressed sparse row format using the  array colIndex
!             3. how to use minors of the matrix A by redefining the arrays pointerB and pointerE
!                but the arrays values and rows are the same.
!
!  In what follows the symbol ' means taking of transposed.
!
!  The test performs the following operations :
!
!       1. The code computes (L+D)'*S = F using MKL_DCSCMM where S is a known 5 by 2
!          matrix and then the code solves the system (L+D)'*X = F with the help of MKL_DCSCSM.
!          It's evident that X should be equal to S.
!
!       2. The code computes (U+I)'*S = F using MKL_DCSCMV where S is a vector
!          and then the code calls MKL_DCSCTRSV solves the system (U+I)'*X = F with the single right
!          hand side. It's evident that X should be equal to S.
!
!       3. The code computes D*S = F using MKL_DCSCMV where S is a vector
!          and then the code solves the system D*X = F with the single right hand side.
!          It's evident that X should be equal to S.
!
!       4. The next step is the computation (U-U') S = F using MKL_DCSCMV where S is
!          a vector. It is easy to see that U-U' is a skew-symmetric matrix.
!
!       5. The next step is the computation (L+D+L') S = F using MKL_DCSCSYMV where S is
!          a vector. It is easy to see that L+D+L' is a symmetric matrix.
!
!       6. The next step is the computation A'* S = F using MKL_DCSCGEMV where S is
!          a vector.
!
!       7. Let's T be the upper 3 by 3 minor of the matrix A. Namely, T is the following matrix
!
!                        |   1       -1      0   |
!          T    =        |  -2        5      0   |.
!                        |   0        0      4   |
!          The test performs the matrix-vector multiply T*S=F with the same arrays values. rows
!          and pointerB used before for the whole matrix A. It is enough to change two values of
!	   array pointerE in order to use the minor under consideration. Then the test solves the system
!          T*X =F using MKL_DCSCSV. The routine MKL_DCSCMV is used for getting matrix-vector multiply.
!
! The code given below uses only one sparse representation for the all operations.
!
!*******************************************************************************
!     Definition arrays for sparse representation of  the matrix A in
!     the compressed sparse column format:
!*******************************************************************************
          integer  m,  nnz, mnew
          parameter( m = 5,  nnz=13, mnew=3)
          real*8  values(nnz)
          integer rows(nnz), colIndex(m+1)
          integer pointerB(m) , pointerE(m)
          data values/1.d0, -2.d0, -4.d0, -1.d0, 5.d0,
     &        8.d0, 4 .d0, 2.d0, -3.d0, 6.d0, 7.d0, 4.d0, -5.d0/
          data colIndex/1, 4,  7,  9,  12, 14/
          data rows/1, 2, 4,  1, 2, 5,  3, 4,  1, 3, 4,  3, 5/
!*******************************************************************************
!    Declaration of local variables :
!*******************************************************************************
          integer n
          parameter (n=2)
          real*8 rhs(m, n), sol(m, n), temp(m, n)
          data sol/1.D0, 1.D0, 1.D0, 1.D0, 1.D0,
     &    5.D0, 4.D0, 3.D0, 2.D0, 1.D0/
          real*8 alpha, beta
          data alpha/1.d0/, beta/0.d0/
          integer i, j
          print*
          print*, ' EXAMPLE PROGRAM FOR COMPRESSED SPARSE ROW
     &              FORMAT ROUTINES '

!*******************************************************************************
!Task 1.    Obtain matrix-matrix multiply (L+D)' *sol --> rhs
!    and solve triangular system   (L+D)' *temp = rhs with multiple right hand sides
!    Array temp must be equal to the array sol
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_DCSCMM '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 101, m, n
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input matrix '
          print 104, ((sol(i,j),j=1,n),i=1,m)

           call mkl_dcscmm('t', m, n, m, alpha, 'tln',
     &       values, rows, colIndex, colIndex(2), sol, m,  beta, rhs,
     &       m)
          print*
          print*, '     OUTPUT DATA FOR MKL_DCSCMM '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 104, ((rhs(i,j),j=1,n),i=1,m)
          print 100
          print*, ' Solve triangular system with obtained '
          print*, ' right hand side  '
          call mkl_dcscsm('t', m, n, alpha, 'tln',
     &           values,  rows, colIndex, colIndex(2), rhs, m, temp,  m)

          print*
          print*, '     OUTPUT DATA FOR MKL_DCSCSM '
          print 104, ((temp(i,j),j=1,n),i=1,m)
          print 100
!*******************************************************************************
! Task 2.    Obtain matrix-vector multiply (U+I)' *sol --> rhs
!    and solve triangular system   (U+I)' *temp = rhs with single right hand sides.
!    Array temp must be equal to the array sol.
!
!    Let us form the arrays pointerB and pointerE for the NIST's variation of the
!    compressed sparse row format using the  array colIndex.
!
!*******************************************************************************
          do i=1, m
            pointerB(i)=colIndex(i)
            pointerE(i)=colIndex(i+1)
          enddo
          print*
          print*, '     INPUT DATA FOR MKL_DCSCMV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i,1),i=1,m)

           call mkl_dcscmv('t', m, m, alpha, 'tuu',
     &           values, rows, pointerB, pointerE, sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_DCSCMV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
          print*, ' Solve triangular system with obtained '
          print*, ' right hand side  '
          call mkl_dcscsv('t', m, alpha, 'tuu',
     &         values, rows, colIndex, colIndex(2),  rhs, temp)
          print*
          print*, '     OUTPUT DATA FOR MKL_DCSCTRSV '
          print*, '     WITH TRIANGULAR MATRIX  '
          print 105, (temp(i,1),i=1,m)
          print 100
!*******************************************************************************
! Task 3.   Obtain matrix-vector multiply D *sol --> rhs
!    and solve triangular system   D *temp = rhs with single right hand side
!    Array temp must be equal to the array sol
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_DCSCMV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i,2),i=1,m)

           call mkl_dcscmv('n', m, m, alpha, 'dun',
     &          values, rows, pointerB, pointerE,  sol(1,2), beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_DCSCMV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
          print*, ' Multiply by inverse matrix '
          print*, ' with the help of MKL_DCSCSV '


          call mkl_dcscsv('t', m, alpha, 'dun',
     &         values, rows, colIndex, colIndex(2),  rhs, temp)

          print*
          print*, '     OUTPUT DATA FOR MKL_DCSCSV '
          print*, '     WITH DIAGONAL MATRIX  '
          print 105, (temp(i,1),i=1,m)
          print 100

!*******************************************************************************
! Task 4.  Obtain matrix-vector multiply (U -U')*sol --> rhs
!    Array temp must be equal to the array sol
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_DCSCMV '
          print*, '     WITH SKEW-SYMMETRIC MATRIX '
          print 102, alpha, beta
          print 103, 'n'
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)

           call mkl_dcscmv('n', m, m, alpha, 'au',
     &          values, rows, colIndex, colIndex(2), sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_DCSCMV '
          print*, '     WITH SKEW-SYMMETRIC MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
!*******************************************************************************
! Task 5.    Obtain matrix-vector multiply (L+D+L')*sol --> rhs whith the help of
!    MKL_DCSCMV
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_DCSCMV '
          print*, '     WITH SYMMETRIC MATRIX  '
          print 102, alpha, beta
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)

          call mkl_dcscmv('n', m, m, alpha, 'sln',
     &          values, rows, colIndex, colIndex(2), sol, beta, rhs)
          print*
          print*, '     OUTPUT DATA FOR MKL_DCSCMV '
          print*, '     WITH SYMMETRIC MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100
!*******************************************************************************
! Task 6.   Obtain matrix-vector multiply A'*sol --> rhs whith the help of MKL_DCSCMV
!
!*******************************************************************************
          print*
          print*, '     INPUT DATA FOR MKL_DCSCMV '
          print*, '     WITH GENERAL MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,m)
          call mkl_dcscmv('t', m, m, alpha, 'g',
     &          values, rows, colIndex, colIndex(2), sol, beta, rhs)

          print*
          print*, '     OUTPUT DATA FOR MKL_DCSCMV '
          print*, '     WITH GENERAL MATRIX  '
          print 105, (rhs(i,1),i=1,m)
          print 100

!*******************************************************************************
! Task 7.  Obtain matrix-vector multiply T*sol --> rhs whith the help of MKL_DCSCMV
!    where S is  3 by 3 minor of the matrix A starting with A(1,1)
!    Let's us redefine three elements of the array pointerE in order to identify
!    the needed minor. More precisely
!            pointerE(1) --> pointerE(1)-1
!            pointerE(2) --> pointerE(2)-1
!            pointerE(3) --> pointerE(3)-1
!
!*******************************************************************************
          do i=1, mnew
              pointerE(i) = pointerE(i) - 1
          enddo
          print*
          print*, '     INPUT DATA FOR MKL_DCSCMV '
          print*, '     WITH A MINOR OF GENERAL MATRIX  '
          print 102, alpha, beta
          print 103, 't'
          print*, ' Input vector '
          print 105, (sol(i, 1),i=1,mnew)
           call mkl_dcscmv('n', mnew, mnew, alpha, 'tln', values,
     &       rows,   pointerB, pointerE,  sol, beta, rhs)

          print*
          print*, '     OUTPUT DATA FOR MKL_DCSCGEMV '
          print*, '     WITH A MINOR OF GENERAL MATRIX  '
          print 105, (rhs(i,1),i=1,mnew)
          print 100

          print*, ' Multiply by inverse to a minor of the matrix '
          print*, ' with the help of MKL_DCSCSV '


          call mkl_dcscsv('n', mnew, alpha, 'tln', values,
     &        rows, pointerB, pointerE,  rhs, temp)

          print*
          print*, '     OUTPUT DATA FOR MKL_DCSCSV '
          print*, '     WITH A MINOR OF GENERAL MATRIX '
          print 105,   (temp(i,1),i=1,mnew)
          print 100
 100      format('------------------------------------------------')
 101      format(7x,'M=',i1,'  N=',i1)
 102      format(7x,'ALPHA= ',f4.1,' BETA= ', f4.1)
 103      format(7x,'TRANS=',a1)
 104      format(2(f7.1, 3x))
 105      format(f4.1)
          stop
          end
