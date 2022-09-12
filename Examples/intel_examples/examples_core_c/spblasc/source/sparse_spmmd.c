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
*             example for mkl_sparse_d_spmmd
*
********************************************************************************
*
* Consider the matrix A
*
*                 |  10     11      0     0     0   |
*                 |   0      0     12    13     0   |
*   A    =        |  15      0      0     0    14   |,
*                 |   0     16     17     0     0   |
*                 |   0      0      0    18    19   |
*
* and diagonal matrix B
*
*                 |   5      0      0     0     0   |
*                 |   0      6      0     0     0   |
*   B    =        |   0      0      7     0     0   |.
*                 |   0      0      0     8     0   |
*                 |   0      0      0     0     9   |
*
*  Both matrices A and B are stored in a zero-based compressed sparse row (CSR) storage
*  scheme with three arrays (see 'Sparse Matrix Storage Schemes' in the
*  Intel MKL Developer Reference) as follows:
*
*           values_A = ( 10  11  12  13  15  14  16  17  18  19 )
*          columns_A = (  0   1   2   3   0   4   1   2   3   4 )
*         rowIndex_A = (  0       2       4       6       8      10 )
*
*           values_B = ( 5  6  7  8  9  )
*          columns_B = ( 0  1  2  3  4  )
*         rowIndex_B = ( 0  1  2  3  4  5 )
*
*  The example computes two scalar products :
*
*         < (A*B)*x ,       y > = left,   using MKL_SPARSE_D_SPMMD and CBLAS_DDOT.
*         <     B*x , (A^t)*y > = right,  using MKL_SPARSE_D_MV and CBLAS_DDOT.
*
*         These products should return the same value. To obtain matrix C,
*         use MKL_SPARSE_D_EXPORT_CSR and print the result.
*
******************************************************************************/

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "mkl.h"

int main() {

#define M 5
#define NNZ 10
#define ALIGN 128

/* To avoid constantly repeating the part of code hat checks inbound SparseBLAS functions' status,
   use macro CALL_AND_CHECK_STATUS */
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = 1;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)

/* Declaration of values */
    double  *values_A = NULL, *values_B = NULL;
    MKL_INT *columns_A = NULL, *columns_B = NULL;
    MKL_INT *rowIndex_A = NULL, *rowIndex_B = NULL;

    double  *rslt_mv = NULL, *rslt_mv_trans = NULL, *x = NULL, *y = NULL, *C = NULL;

    double   left, right, residual;
    MKL_INT  i, j, ii, status;

    struct matrix_descr    descr_type_gen;
    sparse_matrix_t        csrA = NULL, csrB = NULL;

/* Allocation of memory */
    values_A = (double *)mkl_malloc(sizeof(double) * NNZ, ALIGN);
    columns_A = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * NNZ, ALIGN);
    rowIndex_A = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * (M + 1), ALIGN);

    values_B = (double *)mkl_malloc(sizeof(double) * M, ALIGN);
    columns_B = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * M, ALIGN);
    rowIndex_B = (MKL_INT *)mkl_malloc(sizeof(MKL_INT) * (M + 1), ALIGN);

    C = (double *)mkl_malloc(sizeof(double) * pow(M, 2), ALIGN);
    x = (double *)mkl_malloc(sizeof(double) * M, ALIGN);
    y = (double *)mkl_malloc(sizeof(double) * M, ALIGN);
    rslt_mv = (double *)mkl_malloc(sizeof(double) * M, ALIGN);
    rslt_mv_trans = (double *)mkl_malloc(sizeof(double) * M, ALIGN);

/* Set values of the variables*/
    descr_type_gen.type = SPARSE_MATRIX_TYPE_GENERAL;
    ii = 0, status = 0;
 //Matrix A 
    for( i = 0; i < NNZ; i++ )
          values_A[i] = i + 10;
    for( i = 0; i < NNZ; i++ )
          columns_A[i] = i % 5;
    rowIndex_A[0] = 0;
    for( i = 1; i < M + 1; i++ )
          rowIndex_A[i] = rowIndex_A[i - 1] + 2;

 //Matrix B
    for( i = 0; i < M; i++ )
          values_B[i] = i + 5;
    for( i = 0; i < M; i++ )
          columns_B[i] = i % 5;
    for( i = 0; i < M + 1; i++ )
          rowIndex_B[i] = i;

 //Vectors x and y
    for( i = 0; i < M; i++ )
    {
          x[i] = 1.0; y[i] = 1.0;
    }
/* Printing usable data */
    printf( "\n\n_______________Example program for MKL_SPARSE_D_SPMMD_________________\n\n" );
    printf( " COMPUTE  A * B = C, where matrices A and B are stored in CSR format\n" );
    printf( "\n MATRIX A:\nrow# : (value, column) (value, column)\n" );
    for( i = 0; i < M; i++ )
    {
        printf("row#%d:", i + 1); fflush(0);
        for( j = rowIndex_A[i]; j < rowIndex_A[i+1]; j++ )
        {
            printf(" (%5.0f, %6d)", values_A[ii], columns_A[ii] ); fflush(0);
            ii++;
        }
        printf( "\n" );
    }
    ii = 0;
    printf( "\n MATRIX B:\nrow# : (value, column)\n" );
    for( i = 0; i < M; i++ )
    {
        printf("row#%d:", i + 1); fflush(0);
        for( j = rowIndex_B[i]; j < rowIndex_B[i+1]; j++ )
        {
            printf(" (%5.0f, %6d)", values_B[ii], columns_B[ii] ); fflush(0);
            ii++;
        }
        printf( "\n" );
    }
    printf( "\n Check the resultant matrix C, using two scalar products\n" );
    printf( " (values of these scalar products must match).\n" );

/* Prepare arrays which are related to matrices.
   Create handles with matrices A and B stored in CSR format */
    CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &csrA, SPARSE_INDEX_BASE_ZERO, M, M, rowIndex_A, rowIndex_A+1, columns_A, values_A ),
                          "Error after MKL_SPARSE_D_CREATE_CSR, csrA \n");
    CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &csrB, SPARSE_INDEX_BASE_ZERO, M, M, rowIndex_B, rowIndex_B+1, columns_B, values_B ),
                          "Error after MKL_SPARSE_D_CREATE_CSR, csrB \n");

/* Analytic Routines for MKL_SPARSE_D_MV.
   HINTS: provides estimate of number and type of upcoming matrix-vector operations
   OPTIMIZE: analyze sparse matrix; choose proper kernels and workload balancing strategy */
    CALL_AND_CHECK_STATUS(mkl_sparse_set_mv_hint( csrA, SPARSE_OPERATION_TRANSPOSE,     descr_type_gen, 1 ),
                          "Error after MKL_SPARSE_SET_MV_HINT, csrA \n");
    CALL_AND_CHECK_STATUS(mkl_sparse_set_mv_hint( csrB, SPARSE_OPERATION_NON_TRANSPOSE, descr_type_gen, 1 ),
                          "Error after MKL_SPARSE_SET_MV_HINT, csrB \n");

    CALL_AND_CHECK_STATUS(mkl_sparse_optimize( csrA ),
                          "Error after MKL_SPARSE_OPTIMIZE, csrA \n");
    CALL_AND_CHECK_STATUS(mkl_sparse_optimize( csrB ),
                          "Error after MKL_SPARSE_OPTIMIZE, csrB \n");

/* Execution Routines */
/* Step 1:
          Need to compute the following variables:
                       C = A * B
                 rslt_mv = C * x
                    left = <rslt_mv, y>              */

    CALL_AND_CHECK_STATUS(mkl_sparse_d_spmmd( SPARSE_OPERATION_NON_TRANSPOSE, csrA, csrB, SPARSE_LAYOUT_ROW_MAJOR, C, M ),
                          "Error after MKL_SPARSE_D_SPMMD \n");
    cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1.0, C, M, x, 1, 0.0, rslt_mv, 1);
    left = cblas_ddot( M, rslt_mv, 1, y, 1 );

/* Step 2:
          Need to compute the following variables:
           rslt_mv       =     B * x
           rslt_mv_trans = (A)^t * y
                   right = <rslt_mv, rslt_mv_trans>  */

    CALL_AND_CHECK_STATUS(mkl_sparse_d_mv( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrB, descr_type_gen, x, 0.0, rslt_mv ),
                          "Error after MKL_SPARSE_D_MV, csrB*x  \n");
    CALL_AND_CHECK_STATUS(mkl_sparse_d_mv( SPARSE_OPERATION_TRANSPOSE,     1.0, csrA, descr_type_gen, y, 0.0, rslt_mv_trans),
                          "Error after MKL_SPARSE_D_MV, csrA*y  \n");
    right = cblas_ddot( M, rslt_mv, 1, rslt_mv_trans, 1);

/* Step 3:
          Compare values obtained for left and right  */
    residual = fabs(left - right)/(fabs(left)+1);

    printf( "\n The difference between < C*x , y > and < B*x , (A^t)*y > = %g,\n", residual );
    printf( " which means that MKL_SPARSE_D_SPMMD arrived correct at a solution.\n" );

/* Printing OUTPUT DATA */
    printf( "\n RESULTANT DENSE MATRIX C:\n" );
    for( i = 0; i < M; i++ )
    {
        printf("|");fflush(0);
        for( j = 0; j < M; j++ )
        {
             printf("%5.1f ", C[i * M + j]);fflush(0);
        }
        printf("|\n");fflush(0);
    }
    printf( "_____________________________________________________________________  \n" );

/* Deallocate memory */
memory_free:
 //Deallocate arrays for which we allocate memory ourselves.
    mkl_free( C ); mkl_free(rslt_mv_trans); mkl_free(rslt_mv); mkl_free(x); mkl_free(y);

 //Release matrix handle and deallocate arrays for which we allocate memory ourselves.
    if( mkl_sparse_destroy( csrA ) != SPARSE_STATUS_SUCCESS)
    { printf(" Error after MKL_SPARSE_DESTROY, csrA \n");fflush(0); status = 1; }
    mkl_free(values_A); mkl_free(columns_A); mkl_free(rowIndex_A);

    if( mkl_sparse_destroy( csrB ) != SPARSE_STATUS_SUCCESS)
    { printf(" Error after MKL_SPARSE_DESTROY, csrB \n");fflush(0); status = 1; }
    mkl_free(values_B); mkl_free(columns_B); mkl_free(rowIndex_B);


    return status;
}
