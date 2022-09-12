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
*   In this example a system of linear equations (Ax=b) is solved by preconditioned CG
*   method with symmetric Gauss-Zeidel preconditioner:
*   Compute r_0 = b - Ax_0
*   z_0 = r_0 and p_0 = z_0
*   while not converged
*       {
*               alpha_j = (r_j , z_j )/(Ap_j , p_j )
*               x_{j+1} = x_j + beta_j*p_j
*               r_{j+1} = r_j - beta_j*A*p_j
*               w_{j+1} = B^{-1}*r_{j+1}
*               beta_j = (r_{j+1}, z_{j+1})/(r_j , z_j )
*               beta_{j+1} = z_{j+1} + beta_j*p_j
*       }
*
*                       A = -L+D-L^t; B = (D-L)*D^{-1}*(D-L^t).
*       Matrix vector multiplication and solution of triangular system implemented
*   using Sparse BLAS Inspector-Executor functionality
*/
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "mkl_spblas.h"

int main() {
    // ******************************************************************************
    //     Declaration and initialization of parameters for sparse representation of
    //     the matrix A in the compressed sparse row format:
    // ******************************************************************************
#define M 6
#define NNZ 13
    // ******************************************************************************
    //    Sparse representation of the matrix A in CSR format
    // ******************************************************************************
    double csrVal[NNZ]    = { 4.0, -1.0,       -1.0,
                                    4.0, -1.0,       -1.0,
                                          4.0,             -1.0,
                                                4.0, -1.0,
                                                      4.0, -1.0,
                                                            4.0 };
    MKL_INT csrColInd[NNZ] = { 0, 1, 3, 1, 2, 4, 2, 5, 3, 4, 4, 5, 5 };
    MKL_INT csrRowPtr[M+1] = { 0, 3, 6, 8, 10, 12, 13 };
    // Descriptor of main sparse matrix properties
    struct matrix_descr descrA_matrix, descrA_precond;
    // // Structure with sparse matrix stored in CSR format
    sparse_matrix_t       csrA;
    // ******************************************************************************
    //    Declaration of right hand side and initial solution vector:
    // ******************************************************************************
    double rhs[M]  = { 2.0, 1.0, 2.0, 2.0, 1.0, 2.0};
    double x[M] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // ******************************************************************************
    //    Declaration of local variables:
    // ******************************************************************************
    double r[M], w[M], p[M], tmp[M], tmp2[M];
    double one = 1.0, zero = 0.0, alpha, beta, norm_of_correction, initial_norm_of_correction, temp1, temp2;
    MKL_INT    i, k;
    MKL_INT expected_calls = 8;

    printf("   In this example preconditioned CG method with symmetric Gauss - Zeidel preconditioner\n");
    printf("   solved following system of sparse linear system:\n");
    printf("\n");
    printf("      | 4.0, -1.0,  0.0, -1.0, 0.0, 0.0 ||x_1|   |2.0|\n");
    printf("      |-1.0,  4.0, -1.0,  0.0,-1.0, 0.0 ||x_2|   |1.0|\n");
    printf("      | 0.0, -1.0,  4.0,  0.0, 0.0,-1.0 ||x_3| = |2.0|\n");
    printf("      |-1.0,  0.0,  0.0,  4.0,-1.0, 0.0 ||x_4|   |2.0|\n");
    printf("      | 0.0, -1.0,  0.0, -1.0, 4.0,-1.0 ||x_5|   |1.0|\n");
    printf("      | 0.0,  0.0, -1.0,  0.0,-1.0, 4.0 ||x_6|   |2.0|\n");
    printf("\n");
    printf("   Initial solution is equal to zero array\n");
    printf("   Process stopped after reducing norm of correction to 10^{-3} \n");
    printf("\n");
    printf("\n");

    // Create handle with matrix stored in CSR format
    mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ZERO, M, M, csrRowPtr, csrRowPtr+1, csrColInd, csrVal );

    descrA_matrix.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    descrA_matrix.mode = SPARSE_FILL_MODE_UPPER;
    descrA_matrix.diag = SPARSE_DIAG_NON_UNIT;

    mkl_sparse_set_mv_hint  ( csrA, SPARSE_OPERATION_NON_TRANSPOSE, descrA_matrix, expected_calls );

    // Create matrix descriptors for matrix and preconditioner
    descrA_precond.mode = SPARSE_FILL_MODE_UPPER;
    descrA_precond.diag = SPARSE_DIAG_NON_UNIT;

    descrA_precond.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    mkl_sparse_set_sv_hint  ( csrA, SPARSE_OPERATION_TRANSPOSE, descrA_precond, expected_calls);
    descrA_precond.type = SPARSE_MATRIX_TYPE_DIAGONAL;
    mkl_sparse_set_mv_hint  ( csrA, SPARSE_OPERATION_NON_TRANSPOSE, descrA_precond, expected_calls);
    descrA_precond.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    mkl_sparse_set_sv_hint  ( csrA, SPARSE_OPERATION_NON_TRANSPOSE, descrA_precond, expected_calls); 

    // Analyze sparse matrix; choose proper kernels and workload balancing strategy
    mkl_sparse_optimize ( csrA );

    // initial residual equal to RHS cause of zero initial vector
    for( i = 0; i < M; i++) r[i] = rhs[i];

    // Calculation B^{-1}r_0
    {
       mkl_sparse_d_trsv ( SPARSE_OPERATION_TRANSPOSE, 1.0, csrA, descrA_precond, r, tmp);
       descrA_precond.type = SPARSE_MATRIX_TYPE_DIAGONAL;
       mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descrA_precond, tmp, zero, tmp2 );
       descrA_precond.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
       mkl_sparse_d_trsv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descrA_precond, tmp2, w);
    }

    for( i = 0; i < M; i++) p[i] = w[i];

    // Calculate initial norm of correction
    initial_norm_of_correction = 0.0;
    for (i = 0; i < M; i++) initial_norm_of_correction += w[i] * w[i];
    initial_norm_of_correction = sqrt(initial_norm_of_correction);
    norm_of_correction = initial_norm_of_correction;

    // Start of main PCG algorithm
    k = 0;
    temp1 = 0.0;
    for( i = 0; i < M; i++) temp1 += r[i]*w[i];

    while (norm_of_correction / initial_norm_of_correction > 1.e-3 && k<1000)
    {
        // Calculate A*p
        mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descrA_matrix, p, zero, tmp );

        // Calculate alpha_k
        temp2 = 0.0;
        for( i = 0; i < M; i++) temp2 += p[i]*tmp[i];
        alpha = temp1/temp2;

        for( i = 0; i < M; i++) x[i] += alpha*p[i];
        // Calculate r_k = r_k - alpha*A*p_k
        mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, -alpha, csrA, descrA_matrix, p, one, r );

        // Calculate w_k = B^{-1}r_k
        {
           mkl_sparse_d_trsv ( SPARSE_OPERATION_TRANSPOSE, 1.0, csrA, descrA_precond, r, tmp);
           descrA_precond.type = SPARSE_MATRIX_TYPE_DIAGONAL;
           mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descrA_precond, tmp, zero, tmp2 );
           descrA_precond.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
           mkl_sparse_d_trsv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descrA_precond, tmp2, w);
        }

        // Calculate current norm of correction
        norm_of_correction = 0.0;
        for (i = 0; i < M; i++) norm_of_correction += w[i] * w[i];
        norm_of_correction = sqrt(norm_of_correction);
        printf("relative norm of residual on %d iteration = %4.5e\n", ++k, norm_of_correction / initial_norm_of_correction);
        if (norm_of_correction <= 1.e-3) break;

        // Calculate beta_k
        temp2 = 0.0;
        for( i = 0; i < M; i++) temp2 += r[i]*w[i];
        beta = temp2/temp1;
        temp1 = temp2;

        // Calculate p_k = w_k-beta*p_k
        for( i = 0; i < M; i++) p[i] = w[i] + beta*p[i];
    }

    printf("\n");
    printf("Preconditioned CG process successfully converge and following solution have been obtained\n");
    for( i = 0; i < M; i++) printf("x_%d = %4.5f\n",i+1,x[i]);

    // Release matrix handle and deallocate matrix
    mkl_sparse_destroy ( csrA );

    return 0;
}
