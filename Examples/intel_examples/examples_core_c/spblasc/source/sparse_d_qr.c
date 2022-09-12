/*******************************************************************************
* Copyright 2018-2020 Intel Corporation.
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
*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Sparse QR C example
*
********************************************************************************
*
* Consider the sparse rectangular matrix A: (see 'Sparse Storage Formats for
* Sparse BLAS Level 2 and Level 3 in the Intel MKL Reference Manual')
*
*                 |   1       -1      0   -3     0   |
*                 |  -2        5      0    0     0   |
*   A    =        |   0        0      0    6     4   |,
*                 |  -4        0      2    7     0   |
*                 |   0        8      0    0    -5   |
*                 |   0        0      1    0     0   |
*                 |   2        0      0    0    -1   |
*
*  The matrix A is represented in a zero-based compressed sparse row (CSR) storage
*  scheme with three arrays (see 'Sparse Matrix Storage Schemes' in the
*  Intel MKL Reference Manual) as follows:
*
*         values  =  ( 1 -1 -3 -2  5  6  4 -4  2  7  8 -5  1  2 -1 )
*         columns =  ( 0  1  3  0  1  3  4  0  2  3  1  4  2  0  4 )
*         rowIndex = ( 0        3     5     7        10    12 13    15 )
*
*  The example solve system  Ax = b using mkl_sparse_d_qr routine
*
********************************************************************************
*/
#include <stdio.h>
#include <math.h>
#include "mkl.h"

int main() {
/* To avoid constantly repeating the part of code that checks inbound SparseBLAS functions' status,
   use macro CALL_AND_CHECK_STATUS

   SPARSE_STATUS_SUCCESS           = 0,    the operation was successful
   SPARSE_STATUS_NOT_INITIALIZED   = 1,    empty handle or matrix arrays
   SPARSE_STATUS_ALLOC_FAILED      = 2,    internal error: memory allocation failed
   SPARSE_STATUS_INVALID_VALUE     = 3,    invalid input value
   SPARSE_STATUS_EXECUTION_FAILED  = 4,    e.g. 0-diagonal element for triangular solver, etc.
   SPARSE_STATUS_INTERNAL_ERROR    = 5,    internal error
   SPARSE_STATUS_NOT_SUPPORTED     = 6     e.g. operation for double precision doesn't support other types
*/

#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          status = function;                                \
          if(status != SPARSE_STATUS_SUCCESS) {             \
          printf(error_message); fflush(0);                 \
          printf("exit status is %d\n", status); fflush(0); \
          error = 1;                                        \
          goto memory_free;                                 \
          }                                                 \
} while(0)

/* Declaration and initializations of variables */
    sparse_status_t status = SPARSE_STATUS_SUCCESS;
    MKL_INT i, j, error = 0;

    sparse_matrix_t csrA = NULL;
    struct matrix_descr descrA;

    MKL_INT nrows  = 7;
    MKL_INT ncols  = 5;

    /* setting matrix A (rectangular sparse general matrix) */
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    MKL_INT ia[8]  = { 0, 3, 5, 7, 10, 12, 13, 15 };
    MKL_INT ja[15] =
    { 0,    1,          3,
      0,    1,
                        3,    4,
      0,          2,    3,
            1,                4,
                  2,
      0,                      4 };
    double a[15] =
    { 1.0, -1.0,       -3.0,
     -2.0,  5.0,
                        6.0,  4.0,
     -4.0,        2.0,  7.0,
            8.0,             -5.0,
                  1.0,
      2.0,                   -1.0 };

    double x[5] = {0};
    double r[7] = {0};
    double b[7] = { -1.0, 8.0, 4.0, 2.0, 11.0, 3.0, 1.0 };
    double res = 1.0, b_norm, diff_norm;

/* Printing usable data */
    printf( "\n\n_______Example program for MKL_SPARSE_D_QR_______\n\n" );
    printf( " SOLVE  Ax = b, where matrix are stored in CSR format\n"   );
    printf( "\n MATRIX A:\nrow# : (value, column) (value, column) \n"   );
    for( i = 0; i < nrows; i++ )
    {
        printf("row#%d:", i + 1); fflush(0);
        for( j = ia[i]; j < ia[i+1]; j++ )
        {
            printf(" (%5.0f, %6d)", a[j], ja[j] ); fflush(0);
        }
        printf( "\n" );
    }
    printf( "\n RHS b:\nrow# : value \n"   );
    for( i = 0; i < nrows; i++ )
    {
        printf("row#%d : %5.0f\n", i, b[i] ); fflush(0);
    }

/* Create handle for matrix A stored in CSR format */
    CALL_AND_CHECK_STATUS ( mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ZERO, nrows, ncols, ia, ia+1, ja, a ),
                           "Error after MKL_SPARSE_D_CREATE_CSR, csrA \n" );

/* Solve Ax=b using Sparse QR decomposition */
    CALL_AND_CHECK_STATUS ( mkl_sparse_d_qr( SPARSE_OPERATION_NON_TRANSPOSE, csrA, descrA, SPARSE_LAYOUT_COLUMN_MAJOR, 1, x, ncols, b, nrows ),
                           "Error after MKL_SPARSE_D_QR, csrA \n" );

/* Printing solution vector */
    printf( "\n x:\nrow# : value \n"   );
    for( i = 0; i < ncols; i++ )
    {
        printf("row#%d : %5.0f\n", i, x[i] ); fflush(0);
    }

/* Checking correctness of Sparse QR result and free data */
    CALL_AND_CHECK_STATUS ( mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1., csrA, descrA, x, 0., r ),
                           "Error after MKL_SPARSE_D_MV, csrA \n" );

    b_norm = 0.0; diff_norm = 0.0;
    for ( i = 0; i < nrows; i++ ) {
        b_norm    += b[i]*b[i];
        diff_norm += (b[i]-r[i])*(b[i]-r[i]);
    }
    res = sqrt(diff_norm)/(sqrt(b_norm)+1.);
    if ( res < 1.e07 ) { printf("Residual norm is      OK, res = %e\n", res); fflush(0); error = 0; }
    else               { printf("Residual norm is too big, res = %e\n", res); fflush(0); error = 1; }

memory_free:
    if( mkl_sparse_destroy(csrA) != SPARSE_STATUS_SUCCESS )
    { printf(" Error after MKL_SPARSE_DESTROY, csrA \n"); fflush(0); return 1; }

    return error;
}
