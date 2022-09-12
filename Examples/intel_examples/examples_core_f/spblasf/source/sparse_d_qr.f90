!===============================================================================
! Copyright 2019-2020 Intel Corporation.
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
!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Sparse QR Fortran
!             example for mkl_sparse_d_qr routine
!
!*******************************************************************************
!
! Consider the sparse rectangular matrix A: (see 'Sparse Storage Formats for
! Sparse BLAS Level 2 and Level 3 in the Intel MKL Reference Manual')
!
!                 |   1       -1      0   -3     0   |
!                 |  -2        5      0    0     0   |
!   A    =        |   0        0      0    6     4   |,
!                 |  -4        0      2    7     0   |
!                 |   0        8      0    0    -5   |
!                 |   0        0      1    0     0   |
!                 |   2        0      0    0    -1   |
!
!  The matrix A is represented in a zero-based compressed sparse row (CSR) storage
!  scheme with three arrays (see 'Sparse Matrix Storage Schemes' in the
!  Intel MKL Reference Manual) as follows:
!
!         values  =  ( 1 -1 -3 -2  5  6  4 -4  2  7  8 -5  1  2 -1 )
!         columns =  ( 0  1  3  0  1  3  4  0  2  3  1  4  2  0  4 )
!         rowIndex = ( 0        3     5     7        10    12 13    15 )
!
!  The example solve system  Ax = b using mkl_sparse_d_qr routine
!
!*******************************************************************************
PROGRAM SPARSE_D_QR_EXAMPLE

     USE MKL_SPBLAS
     USE MKL_SPARSE_QR
     USE ISO_C_BINDING
     IMPLICIT NONE

     INTEGER info
     INTEGER i, j
 
     TYPE(SPARSE_MATRIX_T) csrA
     TYPE(MATRIX_DESCR) descrA

     INTEGER NROWS, NCOLS
     INTEGER NNZ

     INTEGER, ALLOCATABLE :: csrColInd(:), csrRowPtr(:)
     DOUBLE PRECISION, ALLOCATABLE :: csrVal(:)
     DOUBLE PRECISION, ALLOCATABLE :: x(:), r(:), b(:)

     DOUBLE PRECISION res, b_norm, diff_norm, alpha, beta

     NROWS = 7
     NCOLS = 5
     NNZ = 15

     descrA % TYPE = SPARSE_MATRIX_TYPE_GENERAL

     ALLOCATE(csrRowPtr(NROWS+1))
     ALLOCATE(csrColInd(NNZ))
     ALLOCATE(csrVal(NNZ))

     csrRowPtr = (/ 0, 3, 5, 7, 10, 12, 13, 15 /)
     csrColInd = (/ 0, 1, 3, 0, 1, 3, 4, 0, 2, 3, 1, 4, 2, 0, 4 /)
     csrVal    = (/ 1.0,-1.0,-3.0,-2.0,5.0,6.0,4.0,-4.0,2.0,7.0,8.0,-5.0,1.0,2.0,-1.0 /)

     ALLOCATE(x(NCOLS))
     ALLOCATE(r(NROWS))
     ALLOCATE(b(NROWS))

     b = (/ -1.0,8.0,4.0,2.0,11.0,3.0,1.0 /)

     print*,'EXAMPLE PROGRAM FOR MKL_SPARSE_D_QR'
     print*,'---------------------------------------------------'
     print*,'Input matrix A:'
     do i = 1, NROWS
         print*,'row #',i
         do j = csrRowPtr(i)+1, csrRowPtr(i+1)
             print*,csrColInd(j),csrVal(j)
         enddo
     enddo

!   Create CSR matrix
    info = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ZERO,NROWS,NCOLS,csrRowPtr(1),csrRowPtr(2),csrColInd,csrVal)

!   Solve Ax=b using Sparse QR decomposition
    info = MKL_SPARSE_D_QR( SPARSE_OPERATION_NON_TRANSPOSE, csrA, descrA, SPARSE_LAYOUT_COLUMN_MAJOR, 1, x, ncols, b, nrows )

    print*,'Output vector x:'
    do i = 1, NCOLS
        print*,i,x(i)
    enddo
    
    alpha = 1.0
    beta  = 0.0
    info  = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE, alpha, csrA, descrA, x, beta, r)

    b_norm = 0.0
    diff_norm = 0.0
    do i = 1, NROWS
        b_norm    = b_norm + b(i)*b(i)
        diff_norm = diff_norm + (b(i)-r(i))*(b(i)-r(i))
    enddo
    res = sqrt(diff_norm)/(sqrt(b_norm)+1.)
    print *,'Residual is',res

!   Release internal representation of CSR matrix
    info = MKL_SPARSE_DESTROY(csrA)

     print*,'---------------------------------------------------'

END PROGRAM SPARSE_D_QR_EXAMPLE
