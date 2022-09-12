!===============================================================================
! Copyright 2013-2020 Intel Corporation.
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

!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) IE SpBLAS Fortran
!             native example
!
!*******************************************************************************
!
! Consider the matrix A (see 'Sparse Storage Formats for Sparse BLAS Level 2
! and Level 3 in the Intel MKL Reference Manual')
!
!                 |   1       -1      0   -3     0   |
!                 |  -2        5      0    0     0   |
!   A    =        |   0        0      4    6     4   |,
!                 |  -4        0      2    7     0   |
!                 |   0        8      0    0    -5   |
!
!  The matrix A is represented in a zero-based compressed sparse row (CSR) storage
!  scheme with three arrays (see 'Sparse Matrix Storage Schemes' in the 
!  Intel MKL Reference Manual) as follows:
!
!         values  = ( 1 -1 -3 -2 5 4 6 4 -4 2 7 8 -5 )
!         columns = ( 0 1 3 0 1 2 3 4 0 2 3 1 4 )
!         rowIndex = ( 0  3  5  8  11 13 )
!
!  The test performs the following operations :
!
!       The code computes A*S = F using MKL_SPARSE_S_MV
!          where A is a general sparse matrix and S and F are vectors.
!
!*******************************************************************************
PROGRAM SPMV
!   *****************************************************************************
!   Declaration and initialization of parameters for sparse representation of
!   the matrix A in CSR format:
!   *****************************************************************************
    USE MKL_SPBLAS
    IMPLICIT NONE

    INTEGER M, N, NNZ, i, info
!   *****************************************************************************
!   Sparse representation of the matrix A
!   *****************************************************************************
    INTEGER, ALLOCATABLE :: csrColInd(:), csrRowPtr(:)
    REAL, ALLOCATABLE :: csrVal(:)
!   Matrix descriptor
    TYPE(MATRIX_DESCR) descrA     ! Sparse matrix descriptor
!   CSR matrix representation 
    TYPE(SPARSE_MATRIX_T) csrA    ! Structure with sparse matrix
!   *****************************************************************************
!   Declaration of local variables:
!   *****************************************************************************
    REAL, ALLOCATABLE :: x(:), y(:)
    REAL alpha, beta

    M = 5
    N = 5
    NNZ = 13
    ALLOCATE(csrColInd(NNZ))
    ALLOCATE(csrRowPtr(M+1))
    ALLOCATE(csrVal(NNZ))
    ALLOCATE(x(M))
    ALLOCATE(y(M))
    csrVal = (/ 1.0,-1.0,-3.0,-2.0,5.0,4.0,6.0,4.0,-4.0,2.0,7.0,8.0,-5.0 /)
    csrColInd = (/ 0,1,3,0,1,2,3,4,0,2,3,1,4 /)
    csrRowPtr = (/ 0, 3, 5, 8, 11, 13 /)
    x = (/ 1.0, 5.0, 1.0, 4.0, 1.0 /)
    y = (/ 0.0, 0.0, 0.0, 0.0, 0.0 /)
    alpha = 1.0
    beta  = 0.0

    print*,'EXAMPLE PROGRAM FOR MKL_SPARSE_S_MV'
    print*,'---------------------------------------------------'
    print*,''
    print*,'INPUT DATA FOR MKL_SPARSE_S_MV'
    print*,'WITH GENERAL SPARSE MATRIX'
    print*,'ALPHA =',alpha,'BETA =',beta
    print*,'SPARSE_OPERATION_NON_TRANSPOSE'
    print*,'Input vector'
    do i = 1, M
        print*,x(i)
    enddo

!   Create CSR matrix
    i = MKL_SPARSE_S_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ZERO,M,N,csrRowPtr,csrRowPtr(2),csrColInd,csrVal)

!   Create matrix descriptor
    descrA % TYPE = SPARSE_MATRIX_TYPE_GENERAL

!   Analyze sparse matrix; chose proper kernels and workload balancing strategy
    info = MKL_SPARSE_OPTIMIZE(csrA)

!   Compute y = alpha * A * x + beta * y
    info = MKL_SPARSE_S_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,x,beta,y)

!   Release internal representation of CSR matrix
    info = MKL_SPARSE_DESTROY(csrA)

    print*,''
    print*,'OUTPUT DATA FOR sparseDcsrmv'
    do i = 1, M
        print*,y(i)
    enddo

    print*,'---------------------------------------------------'

END PROGRAM SPMV
