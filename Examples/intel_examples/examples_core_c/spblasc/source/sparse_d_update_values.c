/*******************************************************************************
* Copyright 2019-2020 Intel Corporation.
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
*             example for mkl_sparse_d_update_values() routine
*
********************************************************************************
*
* Consider the matrix A (see 'Sparse Storage Formats for Sparse BLAS Level 2
* and Level 3 in the Intel MKL Reference Manual')
*
*                 |   1    1   1   1   0    0   |
*                 |   1    1   1   1   0    0   |
*   A    =        |   0    0   1   1   1    1   |
*                 |   0    0   1   1   1    1   |
*                 |   1    1   0   0   1    1   |
*                 |   1    1   0   0   1    1   |
*
*  The matrix A is represented in a zero-based compressed sparse row (CSR) storage
*  scheme with three arrays (see 'Sparse Matrix Storage Schemes' in the
*  Intel MKL Reference Manual) as follows:
*
*         values  =  ( 1  ...  1 )
*         columns =  ( 0  1  1  2  0  2 )
*         rowIndex = ( 0  2  4  6 )
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
    //     the matrix A in the block CSR format:
    //*******************************************************************************
#define M 3
#define NNZ 24
#define NNZB 6
#define LB 2
    //*******************************************************************************
    //    Sparse representation of the matrix A
    //*******************************************************************************
    double     bsrVal[NNZ] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
    double new_bsrVal[NNZ] = { 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
                               2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
                               2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0 };
    MKL_INT    bsrColInd[NNZB] = { 0, 1, 1, 2, 0, 2 };
    MKL_INT    bsrRowPtr[M+1]  = { 0, 2, 4, 6 };
    // Descriptor of main sparse matrix properties
    struct matrix_descr descrA;
    // // Structure with sparse matrix stored in CSR format
    sparse_matrix_t       bsrA;
    //*******************************************************************************
    //    Declaration of local variables:
    //*******************************************************************************
    double x[M*LB]  = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
    double y[M*LB]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double z[M*LB]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double alpha = 1.0, beta = 0.0;
    MKL_INT i;

    printf( "\n EXAMPLE PROGRAM FOR BSR format routines from IE Sparse BLAS\n" );
    printf( "-------------------------------------------------------\n" );

    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    // Create handle with matrix stored in BSR format
    mkl_sparse_d_create_bsr ( &bsrA, SPARSE_INDEX_BASE_ZERO,
                                     SPARSE_LAYOUT_ROW_MAJOR,
                                     M,
                                     M,
                                     LB,
                                     bsrRowPtr,
                                     bsrRowPtr+1,
                                     bsrColInd,
                                     bsrVal );

    // Analyze sparse matrix; choose proper kernels and workload balancing strategy
    mkl_sparse_set_mv_hint ( bsrA, SPARSE_OPERATION_NON_TRANSPOSE, descrA, 100 );
    mkl_sparse_optimize ( bsrA );

    //  Task: matrix-vector multiplication
    printf( "                                  \n" );
    printf( "   Task:                          \n" );
    printf( "   mkl_sparse_d_mv                \n" );
    printf( "   WITH GENERAL SPARSE MATRIX     \n" );
    printf( "   ALPHA = %4.1f  BETA = %4.1f    \n", alpha, beta );
    printf( "   SPARSE_OPERATION_NON_TRANSPOSE \n" );
    printf( "   Input vector                   \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f\n", x[i] );
    }

    mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE,
                      alpha,
                      bsrA,
                      descrA,
                      x,
                      beta,
                      y );
    printf( "   OUTPUT DATA FOR mkl_sparse_d_mv \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f\n", y[i] );
    }

    // Update values in a handle
    mkl_sparse_d_update_values ( bsrA, 0, NULL, NULL, new_bsrVal );

    // Analyze sparse matrix; choose proper kernels and workload balancing strategy
    mkl_sparse_set_mv_hint ( bsrA, SPARSE_OPERATION_NON_TRANSPOSE, descrA, 100 );
    mkl_sparse_optimize ( bsrA );

    //  Task: matrix-vector multiplication
    printf( "                                       \n" );
    printf( "   Task:                               \n" );
    printf( "   mkl_sparse_d_mv with updated values \n" );
    printf( "   WITH GENERAL SPARSE MATRIX          \n" );
    printf( "   ALPHA = %4.1f  BETA = %4.1f         \n", alpha, beta );
    printf( "   SPARSE_OPERATION_NON_TRANSPOSE      \n" );
    mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE,
                      alpha,
                      bsrA,
                      descrA,
                      x,
                      beta,
                      z );
    printf( "   OUTPUT DATA FOR mkl_sparse_d_mv \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f\n", z[i] );
    }

    // Release matrix handle and deallocate matrix
    mkl_sparse_destroy ( bsrA );

    printf( "-------------------------------------------------------\n" );
    return 0;
}
