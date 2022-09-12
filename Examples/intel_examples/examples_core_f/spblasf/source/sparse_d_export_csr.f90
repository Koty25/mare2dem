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
!      Computes A*A = B using MKL_SPARSE_SPMM,
!      then export matrix B using MKL_SPARSE_D_EXPORT_CSR routine.
!
!*******************************************************************************
PROGRAM EXPORT_CSR

    USE MKL_SPBLAS
    USE ISO_C_BINDING
    IMPLICIT NONE
!   *****************************************************************************
!   Sparse representation of the matrices A and B:
!   *****************************************************************************
    INTEGER, ALLOCATABLE :: csrColInd(:), csrRowPtr(:)
    DOUBLE PRECISION, ALLOCATABLE :: csrVal(:)
!   CSR matrix structure
    TYPE(SPARSE_MATRIX_T) csrA, csrB
!   Variables used for exporting sparse matrix
    INTEGER        :: nrows, ncols
    INTEGER(C_INT) :: indexing
    TYPE(C_PTR)    :: rows_start_c, rows_end_c, col_indx_c, values_c
    INTEGER         , POINTER :: rows_start_f(:), rows_end_f(:), col_indx_f(:)
    DOUBLE PRECISION, POINTER :: values_f(:)
!   *****************************************************************************
!   Declaration of local variables:
!   *****************************************************************************
    INTEGER M, N, NNZ, i, j, info
    M = 5
    NNZ = 13
    ALLOCATE(csrColInd(NNZ))
    ALLOCATE(csrRowPtr(M+1))
    ALLOCATE(csrVal(NNZ))
    csrVal = (/ 1.0,-1.0,-3.0,-2.0,5.0,4.0,6.0,4.0,-4.0,2.0,7.0,8.0,-5.0 /)
    csrColInd = (/ 0,1,3,0,1,2,3,4,0,2,3,1,4 /)
    csrRowPtr = (/ 0,    3,  5,    8,   11,  13 /)

    print*,'EXAMPLE PROGRAM FOR MKL_SPARSE_D_EXPORT_CSR'
    print*,'---------------------------------------------------'
    print*,'Input matrix A:'
    do i = 1, M
        print*,'row #',i
        do j = csrRowPtr(i)+1, csrRowPtr(i+1)
            print*,csrColInd(j),csrVal(j)
        enddo
    enddo

!   Create CSR matrix
    info = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ZERO,M,M,csrRowPtr(1),csrRowPtr(2),csrColInd,csrVal)

!   Compute B = A*A
    info = MKL_SPARSE_SPMM(SPARSE_OPERATION_NON_TRANSPOSE, csrA, csrA, csrB)

!   Export CSR matrix
    info = MKL_SPARSE_D_EXPORT_CSR(csrB, indexing, nrows, ncols, rows_start_c, rows_end_c, col_indx_c, values_c)

!   Converting C into Fortran pointers
    call C_F_POINTER(rows_start_c, rows_start_f, [nrows])
    call C_F_POINTER(rows_end_c  , rows_end_f  , [nrows])
    call C_F_POINTER(col_indx_c  , col_indx_f  , [rows_end_f(nrows)])
    call C_F_POINTER(values_c    , values_f    , [rows_end_f(nrows)])

!   Printing resulting matrix
    print*,'---------------------------------------------------'
    print*,'Output matrix B = A*A:'
    do i = 1, nrows
        print*,'row #',i
        do j = rows_start_f(i)+1, rows_end_f(i)
            print*,col_indx_f(j),values_f(j)
        enddo
    enddo

!   Release internal representation of CSR matrix
    info = MKL_SPARSE_DESTROY(csrA)
    info = MKL_SPARSE_DESTROY(csrB)

    print*,'---------------------------------------------------'

END PROGRAM EXPORT_CSR
