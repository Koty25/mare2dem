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
*             example for CSR format
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
********************************************************************************
*/
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "mkl_spblas.h"

int main() {
    //*******************************************************************************
    //     Declaration and initialization of parameters for sparse representation of
    //     the matrix A in the CSR format:
    //*******************************************************************************
#define M 5
#define NRHS 2
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
    double x_m[M*NRHS]  = { 1.0, 5.0, 3.0, 4.0, 2.0, 2.0, 10.0, 6.0, 8.0, 4.0};
    double y_m[M*NRHS]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double x_v[M]  = { 3.0, 2.0, 5.0, 4.0, 1.0};
    double y_v[M]  = { 0.0, 0.0, 0.0, 0.0, 0.0};
    double tmp_v[M]  = { 0.0, 0.0, 0.0, 0.0, 0.0};
    double alpha = 1.0, beta = 0.0;
    MKL_INT    i;

    printf( "\n EXAMPLE PROGRAM FOR CSR format routines from IE Sparse BLAS\n" );
    printf( "-------------------------------------------------------\n" );

    // Create handle with matrix stored in CSR format
    mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ZERO,
                                    M,  // number of rows
                                    M,  // number of cols
                                    csrRowPtr,
                                    csrRowPtr+1,
                                    csrColInd,
                                    csrVal );

    // Analyze sparse matrix; choose proper kernels and workload balancing strategy
    mkl_sparse_optimize ( csrA );

//  Task 1: Obtain matrix-matrix multiply (L+D)' *x_v --> y_v
//          and solve triangular system   (L+D)' *tmp_v = y_v
//          Array tmp_v must be equal to the array x_v
    printf( "                                  \n" );
    printf( "   Task 1:                        \n" );
    printf( "   INPUT DATA FOR mkl_sparse_d_mv \n" );
    printf( "   WITH TRIANGULAR SPARSE MATRIX  \n" );
    printf( "   ALPHA = %4.1f  BETA = %4.1f    \n", alpha, beta );
    printf( "   SPARSE_OPERATION_NON_TRANSPOSE     \n" );
    printf( "   Input vector                   \n" );
    for ( i = 0; i < M; i++ )
    {
        printf( "%7.1f\n", x_v[i] );
    }

    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    descrA.mode = SPARSE_FILL_MODE_LOWER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;

    mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE,
                      alpha,
                      csrA,
                      descrA,
                      x_v,
                      beta,
                      y_v );
    printf( "   OUTPUT DATA FOR mkl_sparse_d_mv \n" );
    for ( i = 0; i < M; i++ )
    {
        printf( "%7.1f\n", y_v[i] );
    }

    printf("   Solve triangular system   \n");
    printf("   with obtained             \n");
    printf("   right hand side           \n");

    mkl_sparse_d_trsv ( SPARSE_OPERATION_NON_TRANSPOSE,
                      alpha,
                      csrA,
                      descrA,
                      y_v,
                      tmp_v );
    printf( "   OUTPUT DATA FOR mkl_sparse_d_trsv \n" );
    for ( i = 0; i < M; i++ )
    {
        printf( "%7.1f\n", tmp_v[i] );
    }
    printf( "-------------------------------------------------------\n" );

//  Task 2: Obtain matrix-matrix multiply (U+I)' *x_v --> y_v
//          and solve triangular system   (U+I)' *tmp_v = y_v
//          Array tmp_v must be equal to the array x_v
    printf( "                                  \n" );
    printf( "   Task 2:                        \n" );
    printf( "   INPUT DATA FOR mkl_sparse_d_mv \n" );
    printf( "   WITH TRIANGULAR SPARSE MATRIX  \n" );
    printf( "   ALPHA = %4.1f  BETA = %4.1f    \n", alpha, beta );
    printf( "   SPARSE_OPERATION_TRANSPOSE     \n" );
    printf( "   Input vector                   \n" );
    // Release matrix handle and deallocate matrix
    for ( i = 0; i < M; i++ )
    {
        printf( "%7.1f\n", x_v[i] );
    }

    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_UNIT;

    mkl_sparse_d_mv ( SPARSE_OPERATION_TRANSPOSE,
                      alpha,
                      csrA,
                      descrA,
                      x_v,
                      beta,
                      y_v );
    printf( "   OUTPUT DATA FOR mkl_sparse_d_mv \n" );
    for ( i = 0; i < M; i++ )
    {
        printf( "%7.1f\n", y_v[i] );
    }

    printf("   Solve triangular system   \n");
    printf("   with obtained             \n");
    printf("   right hand side           \n");

    mkl_sparse_d_trsv ( SPARSE_OPERATION_TRANSPOSE,
                      alpha,
                      csrA,
                      descrA,
                      y_v,
                      tmp_v );
    printf( "   OUTPUT DATA FOR mkl_sparse_d_trsv \n" );
    for ( i = 0; i < M; i++ )
    {
        printf( "%7.1f\n", tmp_v[i] );
    }
    printf( "-------------------------------------------------------\n" );

//  Task 3: Obtain matrix-matrix multiply A' *x_m --> y_m
//          A - zero-based indexing,
//          x_m - column major ordering
    printf( "                                  \n" );
    printf( "   Task 3:                        \n" );
    printf( "   INPUT DATA FOR mkl_sparse_d_mm \n" );
    printf( "   WITH GENERAL SPARSE MATRIX     \n" );
    printf( "   COLUMN MAJOR ORDERING for RHS  \n" );
    printf( "   ALPHA = %4.1f  BETA = %4.1f    \n", alpha, beta );
    printf( "   SPARSE_OPERATION_TRANSPOSE     \n" );
    printf( "   Input vectors                  \n" );
    for ( i = 0; i < M; i++ )
    {
        printf( "%7.1f, %7.1f\n", x_m[i], x_m[M+i] );
    }

    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    mkl_sparse_d_mm ( SPARSE_OPERATION_TRANSPOSE,
                      alpha,
                      csrA,
                      descrA,
                      SPARSE_LAYOUT_COLUMN_MAJOR,
                      x_m,
                      NRHS,   // number of right hand sides
                      M,      // ldx
                      beta,
                      y_m,
                      M );    // ldy

    printf( "   OUTPUT DATA FOR mkl_sparse_d_mv \n" );
    for ( i = 0; i < M; i++ )
    {
        printf( "%7.1f, %7.1f\n", y_m[i], y_m[M+i] );
    }
    printf( "-------------------------------------------------------\n" );

//  Task 4: Obtain matrix-matrix multiply A*x_m --> y_m
//          A - zero-based indexing,
//          x_m - row major ordering
    printf( "                                  \n" );
    printf( "   Task 4:                        \n" );
    printf( "   INPUT DATA FOR mkl_sparse_d_mm \n" );
    printf( "   WITH GENERAL SPARSE MATRIX     \n" );
    printf( "   ROW MAJOR ORDERING for RHS     \n" );
    printf( "   ALPHA = %4.1f  BETA = %4.1f    \n", alpha, beta );
    printf( "   SPARSE_OPERATION_TRANSPOSE     \n" );
    printf( "   Input vectors                  \n" );
    for ( i = 0; i < M; i++ )
    {
        printf( "%7.1f, %7.1f\n", x_m[2*i], x_m[2*i+1] );
    }

    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    mkl_sparse_d_mm ( SPARSE_OPERATION_TRANSPOSE,
                      alpha,
                      csrA,
                      descrA,
                      SPARSE_LAYOUT_ROW_MAJOR,
                      x_m,
                      NRHS,   // number of right hand sides
                      NRHS,   // ldx
                      beta,
                      y_m,
                      NRHS ); // ldy

    printf( "   OUTPUT DATA FOR mkl_sparse_d_mv \n" );
    for ( i = 0; i < M; i++ )
    {
        printf( "%7.1f, %7.1f\n", y_m[2*i], y_m[2*i+1] );
    }

    // Release matrix handle and deallocate matrix
    mkl_sparse_destroy ( csrA );

    printf( "-------------------------------------------------------\n" );
    return 0;
}
