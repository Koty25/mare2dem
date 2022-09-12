/*******************************************************************************
* Copyright 2020 Intel Corporation.
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
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include "mkl.h"
#include "mkl_graph.h"

#define M 5
#define NNZ 10
#define ALIGN 4096

// To avoid repeating the code that checks graph functions' status, use macro
// CALL_AND_CHECK_STATUS.
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != MKL_GRAPH_STATUS_SUCCESS)          \
          {                                                 \
          printf(error_message); printf("\n"); fflush(0);   \
          error = 1;                                        \
          goto memory_free;                                 \
          }                                                 \
} while(0)

//******************************************************************************
// Task:    Graph matrix-vector multiply
//          y = A*x using the adjacency matrix
//          A, where zeros denote missing edges.
//
//           |  10  11   0   0   0   |
//           |   0   0  12  13   0   |
//   A     = |  15   0   0   0  14   |
//           |   0  16  17   0   0   |
//           |   0   0   0  18  19   |
//
// We represent A in CSR format with the arrays A_values, A_rows_start,
// A_col_indx as follows,
//
//   A_values = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
//   A_col_indx = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4},
//   A_rows_start = {0, 2, 4, 6, 8}.
//
//******************************************************************************

int main()
{
    int64_t i, j;
    int32_t error = 0;

    float   *A_values = NULL;
    int64_t *A_rows_start = NULL;
    int64_t *A_col_indx = NULL;
    int     A_indexing;
    int64_t nrows = M;
    int64_t ncols = M;
    float   *x_values = NULL;
    float   *y_values = NULL;

    mkl_graph_matrix_t Amat = NULL;
    mkl_graph_vector_t xvec = NULL;
    mkl_graph_vector_t yvec = NULL; 

    int64_t          ydim;
    float            *yout_values = NULL;
    float            *yref_values = NULL;
    mkl_graph_type_t yout_values_type;
    const float      float_tol = 1.0e-5;

    int64_t nrows_exp;
    float   *vals_exp;

    // The user allocates memory for A, x, and y.
    A_rows_start    = (int64_t*)mkl_malloc((M + 1)*sizeof(int64_t), ALIGN);
    A_col_indx      = (int64_t*)mkl_malloc(  NNZ  *sizeof(int64_t), ALIGN);
    A_values        = (float*  )mkl_malloc(  NNZ  *sizeof(float),   ALIGN);
    x_values        = (float*  )mkl_malloc(  M    *sizeof(float),   ALIGN);
    y_values        = (float*  )mkl_malloc(  M    *sizeof(float),   ALIGN);
    yref_values     = (float*  )mkl_malloc(  M    *sizeof(float),   ALIGN);

    if (A_rows_start == NULL || A_col_indx == NULL || A_values == NULL ||
        x_values == NULL || y_values == NULL || yref_values == NULL )
    {
        printf("Cannot allocate data \n");
        error = -1;
        goto memory_free;
    }

    // Populate CSR arrays for A.
    for( i = 0; i < NNZ; i++ )
        A_values[i] = i + 10;
    for( i = 0; i < NNZ; i++ )
        A_col_indx[i] = i % 5;

    A_rows_start[0] = 0;
    for( i = 1; i < M + 1; i++ )
        A_rows_start[i] = A_rows_start[i - 1] + 2;
    A_indexing = A_rows_start[0];

    // Populate values for the vector x.
    for( i = 0; i < M; i++ )
        x_values[i] = i * 1.0;

    // Compute values for the reference result y.
    for( i = 0; i < M; i++ ) {
        yref_values[i] = 0.0;
        printf("row %" PRId64 "\n", i);
        for( j = A_rows_start[i] - A_indexing;
             j < A_rows_start[i+1] - A_indexing; j++)
        {
            yref_values[i] += A_values[j] * x_values[A_col_indx[j]
                              - A_indexing];
            printf("(%" PRId64 ", %3.3f) ", A_col_indx[j], A_values[j]);
        }
        printf("x_values[%" PRId64 "] = %3.3f \n", i, x_values[i]);
    }


    // Create an mkl_graph_matrix representing A using A_rows_start, A_col_indx,
    // and A_values.
    // NOTE: The mkl_graph_matrix assumes ownership of these arrays for editing,
    // but does NOT assume responsibility for freeing memory.
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_create(&Amat),
                          "Error after calling mkl_graph_matrix_create for A");

    CALL_AND_CHECK_STATUS(mkl_graph_matrix_set_csr(Amat, nrows, ncols,
                                                   A_rows_start,
                                                   MKL_GRAPH_TYPE_INT64,
                                                   A_col_indx,
                                                   MKL_GRAPH_TYPE_INT64,
                                                   A_values,
                                                   MKL_GRAPH_TYPE_FP32),
                          "Error after calling mkl_graph_matrix_set_csr for A");

    // For applications with many matrix-vector multiplies, calling
    // mkl_graph_optimize_mxv before use can decrease overall runtime.
    CALL_AND_CHECK_STATUS(mkl_graph_optimize_mxv(NULL, MKL_GRAPH_SEMIRING_PLUS_TIMES_FP32,
                                                 Amat, NULL, NULL, 1),
                          "Error after calling mkl_graph_optimize_mxv");

    // Create mkl_graph_vectors representing x and y using x_values and
    // y_values.
    // NOTE: The mkl_graph_vector objects assume ownership of these arrays for
    // editing, but does NOT assume responsibility for freeing memory.
    CALL_AND_CHECK_STATUS(mkl_graph_vector_create(&xvec),
                          "Error after calling mkl_graph_vector_create for x");
    CALL_AND_CHECK_STATUS(mkl_graph_vector_set_dense(xvec, ncols, x_values,
                                                     MKL_GRAPH_TYPE_FP32),
                          "Error after calling mkl_graph_vector_set_dense "
                          "for x");

    CALL_AND_CHECK_STATUS(mkl_graph_vector_create(&yvec),
                          "Error after calling mkl_graph_vector_create for y");
    CALL_AND_CHECK_STATUS(mkl_graph_vector_set_dense(yvec, ncols, y_values,
                                                     MKL_GRAPH_TYPE_FP32),
                          "Error after calling mkl_graph_vector_set_dense for y");

    // Compute y = A * x with NULL passed for the mask and the descriptor.
    CALL_AND_CHECK_STATUS(mkl_graph_mxv(yvec, NULL, MKL_GRAPH_ACCUMULATOR_NONE,
                                        MKL_GRAPH_SEMIRING_PLUS_TIMES_FP32,
                                        Amat, xvec, NULL,
                                        MKL_GRAPH_REQUEST_COMPUTE_ALL,
                                        MKL_GRAPH_METHOD_AUTO),
                          "Error after calling mkl_graph_mxv to compute "
                          "y = A * x");

    // Export the values of y.
    // NOTE: The pointer yout_values can be safely used to access the data but
    // cannot be safely used for editing. It's possible but not guaranteed that
    // this pointer will be equal to the input y_values.
    CALL_AND_CHECK_STATUS(mkl_graph_vector_get_dense(yvec, &ydim,
                                                     (void**)&yout_values,
                                                     &yout_values_type),
                          "Error after calling mkl_graph_vector_get_dense "
                          "for y");

    // Export the values of A.
    // NOTE: Any unneeded arguments can be safely passed as NULL.
    // NOTE: The pointer vals_exp can safely be used to access the data but
    // cannot safely be edited. It's possible but not guaranteed that this
    // pointer will be equal to the input A_values.
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_get_csr(Amat, NULL, NULL, NULL, NULL,
                                                   NULL, NULL,
                                                   (void**)&vals_exp, NULL),
                          "Error after calling mkl_graph_matrix_get_csr for A");
    for( i = 0; i < NNZ; i++) {
        printf("vals_exp[%" PRId64 "] = %3.3f \n", i, vals_exp[i]);
    }

    // Export number of rows in A.
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_get_property(Amat,
                                                       MKL_GRAPH_PROPERTY_NROWS,
                                                        &nrows_exp),
                          "Error after calling mkl_graph_matrix_get_property "
                          "for A");
    printf("nrows_exp = %" PRId64 "\n", nrows_exp);


    printf("Validation:\n");
    printf("y: dim = %" PRId64 ", type = %d \n", ydim, yout_values_type);
    double error_l2 = 0.0, error_max = 0.0;
    for( i = 0; i < ydim; i++) {
        double tmp = fabs(1.0 * (yout_values[i] - yref_values[i]));
        error_l2 += tmp * tmp;
        if ( tmp > error_max )
            error_max = tmp;
        printf("yout_values[%" PRId64 "] = %e yref_values[%" PRId64 "] = %e\n",
               i, yout_values[i], i, yref_values[i]);
    }
    error_l2 = sqrt(error_l2 / ydim);
    printf("error_l2 = %e error_max = %e \n", error_l2, error_max);

    if (error_l2 < float_tol && error_max < float_tol)
        printf("Success, test PASSED!\n");
    else {
        printf("Validation FAILED: error_l2 = %e error_max = %e while float_tol = "
               "%e!\n", error_l2, error_max, float_tol);
        error = -1;
        goto memory_free;
    }

memory_free:
    // Cleanup workflow should destroy all mkl_graph objects and free all
    // user-allocated memory used as input to mkl_graph objects. Access pointers
    // like vals_exp and yout_values that are set by calls to
    // mkl_graph_<object>_get_<representation> do not need to be freed.

    // Calls to mkl_graph_<object>_destroy de-allocate all internal data.
    mkl_graph_matrix_destroy(&Amat);
    mkl_graph_vector_destroy(&xvec);
    mkl_graph_vector_destroy(&yvec);
    // Free user-allocated data. mkl_graph_<object>_destroy leaves
    // user-allocated memory intact.
    mkl_free(A_rows_start);
    mkl_free(A_col_indx);
    mkl_free(A_values);
    mkl_free(x_values);
    mkl_free(y_values);
    mkl_free(yref_values);

    return error;
}

