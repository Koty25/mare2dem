/*******************************************************************************
* Copyright 2013-2020 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/*
*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) IE Sparse BLAS C
*             example for mkl_sparse_d_mv
*
********************************************************************************
*
* Consider the matrix A (see 'Sparse Storage Formats for Sparse BLAS Level 2
* and Level 3 in the Intel MKL Reference Manual')
*
*                 |   1       -1      0   -3     0   |
*                 |  -2        5      0    0     0   |
*   A    =        |   0        0      4    6     4   |,
*                 |  -4        0      2    7     0   |
*                 |   0        8      0    0    -5   |
*
*  The matrix A is represented in a zero-based compressed sparse row (CSR) storage
*  scheme with three arrays (see 'Sparse Matrix Storage Schemes' in the
*  Intel MKL Reference Manual) as follows:
*
*         values  =  ( 1 -1 -3 -2  5  4  6  4 -4  2  7  8 -5 )
*         columns =  ( 0  1  3  0  1  2  3  4  0  2  3  1  4 )
*         rowIndex = ( 0        3     5        8       11    13 )
*
*  The test computes the following operations :
*
*       A*x = y using mkl_sparse_d_mv
*       where A is a general sparse matrix and x and y are vectors
*
********************************************************************************
*/
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "mkl_spblas.h"

int main() {
    //*******************************************************************************
    //     Declaration and initialization of parameters for sparse representation of
    //     the matrix A in the compressed sparse row format:
    //*******************************************************************************
#define M 5
#define N 5
#define NNZ 13
    //*******************************************************************************
    //    Sparse representation of the matrix A
    //*******************************************************************************
    double csrVal[NNZ]    = { 1.0, -1.0,     -3.0,
                             -2.0,  5.0,
                                         4.0, 6.0, 4.0,
                             -4.0,       2.0, 7.0,
                                    8.0,          -5.0 };
    MKL_INT    csrColInd[NNZ] = { 0,      1,        3,
                              0,      1,
                                           2,   3,   4,
                              0,           2,   3,
                                      1,             4 };
    MKL_INT    csrRowPtr[M+1] = { 0, 3, 5, 8, 11, 13 };
    // Descriptor of main sparse matrix properties
    struct matrix_descr descrA;
    // // Structure with sparse matrix stored in CSR format
    sparse_matrix_t       csrA;
    //*******************************************************************************
    //    Declaration of local variables:
    //*******************************************************************************
    double x[N]  = { 1.0, 5.0, 1.0, 4.0, 1.0};
    double y[N]  = { 0.0, 0.0, 0.0, 0.0, 0.0};
    double alpha = 1.0, beta = 0.0;
    MKL_INT    i;

    printf( "\n EXAMPLE PROGRAM FOR mkl_sparse_d_mv \n" );
    printf( "---------------------------------------------------\n" );
    printf( "\n" );
    printf( "   INPUT DATA FOR mkl_sparse_d_mv    \n" );
    printf( "   WITH GENERAL SPARSE MATRIX     \n" );
    printf( "   ALPHA = %4.1f  BETA = %4.1f    \n", alpha, beta );
    printf( "   SPARSE_OPERATION_NON_TRANSPOSE \n" );
    printf( "   Input vector                   \n" );
    for ( i = 0; i < N; i++ )
    {
        printf( "%7.1f\n", x[i] );
    };

    // Create handle with matrix stored in CSR format
    mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ZERO,
                                    N,  // number of rows
                                    M,  // number of cols
                                    csrRowPtr,
                                    csrRowPtr+1,
                                    csrColInd,
                                    csrVal );

    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    // Analyze sparse matrix; choose proper kernels and workload balancing strategy
    mkl_sparse_optimize ( csrA );

    // Compute y = alpha * A * x + beta * y
    mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE,
                      alpha,
                      csrA,
                      descrA,
                      x,
                      beta,
                      y );

    // Release matrix handle and deallocate matrix
    mkl_sparse_destroy ( csrA );

    printf( "                                   \n" );
    printf( "   OUTPUT DATA FOR mkl_sparse_d_mv \n" );

    // y should be equal { -16.0, 23.0, 32.0, 26.0, 35.0 }
    for ( i = 0; i < N; i++ )
    {
        printf( "%7.1f\n", y[i] );
    };

    printf( "---------------------------------------------------\n" );
    return 0;
}
