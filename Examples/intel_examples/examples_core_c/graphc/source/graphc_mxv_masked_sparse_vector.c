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

#define PRINT_REFERENCE_OUTPUT

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
// Task:    Graph masked matrix-vector multiply y<m> = A*x using the adjacency matrix
//          A, where zeros denote missing edges, sparse input vector x, mask is also sparse.
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
//   A_values     = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
//   A_col_indx   = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4},
//   A_rows_start = {0, 2, 4, 6, 8}.
//
//******************************************************************************

int main()
{
    int32_t error = 0;
    int64_t i, count;

    float   *A_values     = NULL;
    int64_t *A_rows_start = NULL;
    int64_t *A_col_indx   = NULL;

    int64_t  nrows = M;
    int64_t  ncols = M;

    mkl_graph_vector_t yvec = NULL;
    mkl_graph_vector_t mvec = NULL;
    mkl_graph_matrix_t Amat = NULL;
    mkl_graph_vector_t xvec = NULL;

    int64_t          ydim;
    mkl_graph_type_t yout_indices_type;
    mkl_graph_type_t yout_values_type;

    int64_t  out_nnz_ref     = 0;
    float   *out_values_ref  = NULL;
    int64_t *out_indices_ref = NULL;

    int64_t  out_nnz     = 0;
    float   *out_values  = NULL;
    int64_t *out_indices = NULL;

    int32_t *m_values_sparse  = NULL;
    int64_t *m_indices_sparse = NULL;
    int64_t  m_nnz_sparse     = 3;

    float   *x_values_sparse  = NULL;
    int64_t *x_indices_sparse = NULL;
    int64_t  x_nnz_sparse     = 3;

    // creating the mask (sparse vector)
    m_values_sparse  = (int32_t*)mkl_malloc(m_nnz_sparse * sizeof(int32_t), ALIGN);
    m_indices_sparse = (int64_t*)mkl_malloc(m_nnz_sparse * sizeof(int64_t), ALIGN);

    if (m_values_sparse == NULL || m_indices_sparse == NULL) {
        printf("Failed to allocate memory for sparse mask arrays\n");
        error = -1;
        goto memory_free;
    }

    m_indices_sparse[0] = 1;
    m_indices_sparse[1] = 2;
    m_indices_sparse[2] = 3;

    m_values_sparse[0] = 1;
    m_values_sparse[1] = 1;
    m_values_sparse[2] = 1;

    printf("Mask as a sparse vector consists of:\n");
    for( i = 0; i < m_nnz_sparse; i++ ) {
        printf("(%" PRId64 ", %3.3f)\n", m_indices_sparse[i], m_values_sparse[i]);
    }

    // creating input sparse vector
    x_values_sparse  = (float*  )mkl_malloc(x_nnz_sparse * sizeof(float)  , ALIGN);
    x_indices_sparse = (int64_t*)mkl_malloc(x_nnz_sparse * sizeof(int64_t), ALIGN);

    if (x_values_sparse == NULL || x_indices_sparse == NULL) {
        printf("Failed to allocate memory for sparse vector arrays \n");
        error = -1;
        goto memory_free;
    }

    x_indices_sparse[0] = 1;
    x_indices_sparse[1] = 2;
    x_indices_sparse[2] = 3;

    x_values_sparse[0] = 0.1;
    x_values_sparse[1] = 1.1;
    x_values_sparse[2] = 2.1;

    printf("Input sparse vector consists of:\n");
    for( i = 0; i < x_nnz_sparse; i++ ) {
        printf("(%" PRId64 ", %3.3f)\n", x_indices_sparse[i], x_values_sparse[i]);
    }

    // The user allocates memory for A, x, and y.
    A_rows_start = (int64_t*)mkl_malloc((M + 1)*sizeof(int64_t), ALIGN);
    A_col_indx   = (int64_t*)mkl_malloc(  NNZ  *sizeof(int64_t), ALIGN);
    A_values     = (float*  )mkl_malloc(  NNZ  *sizeof(float),   ALIGN);

    if (A_rows_start == NULL || A_col_indx == NULL || A_values == NULL)
    {
        printf("Cannot allocate data for the matrix A\n");
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

    // Compute values for the reference result y.
    // First, compute the number of nonzero elements in the output
    out_nnz_ref = 0;
    for (i = 0; i < m_nnz_sparse; i++) {
        int64_t  m_index    = m_indices_sparse[i];
        int64_t *shifted_ja = A_col_indx + A_rows_start[m_index];
        int64_t  row_nnz    = A_rows_start[m_index + 1] - A_rows_start[m_index];

        int64_t k = 0, j = 0;
        while (k < row_nnz && j < x_nnz_sparse) {
            int64_t colA = shifted_ja[k];
            float   indX = x_indices_sparse[j];
            if (colA > indX) {
                j++;
            } else if (colA < indX) {
                k++;
            } else { // equal
                out_nnz_ref++;
                break;
            }
        }
    }
    printf("out_nnz_ref = %" PRId64 "\n", out_nnz_ref);

    // Second, compute the indices and values of the nonzero elements in the output

    out_values_ref  = (float*  )mkl_malloc(out_nnz_ref * sizeof(float)  , ALIGN);
    out_indices_ref = (int64_t*)mkl_malloc(out_nnz_ref * sizeof(int64_t), ALIGN);

    if (out_values_ref == NULL || out_indices_ref == NULL) {
        printf("Failed to allocate memory for reference output sparse vector arrays\n");
        error = -1;
        goto memory_free;
    }

    for (i = 0; i < out_nnz_ref; i++){
        out_indices_ref[i] = -1;
        out_values_ref [i] = 0;
    }

    // computing indices and values of the nonzero elements in the output vector
    count = 0;
    for (i = 0; i < m_nnz_sparse; i++) {
        int64_t  m_index    = m_indices_sparse[i];
        int64_t *shifted_ja = A_col_indx + A_rows_start[m_index];
        float   *shifted_a  = A_values + A_rows_start[m_index];
        int64_t  row_nnz    = A_rows_start[m_index + 1] - A_rows_start[m_index];

        int64_t k = 0, j = 0, found = 0;

        while (k < row_nnz && j < x_nnz_sparse) {
            int64_t colA = shifted_ja[k];
            float   indX = x_indices_sparse[j];
            if (colA > indX) {
                j++;
            } else if (colA < indX) {
                k++;
            } else { // equal
                if (out_indices_ref[count] == -1)
                    out_indices_ref[count] = m_index;
                out_values_ref[count] += shifted_a[k] * x_values_sparse[j];
                found = 1;
                j++;
                k++;
            }
        }
        if (found == 1)
            count++;
    }

#ifdef PRINT_REFERENCE_OUTPUT
    printf("Reference output vector:\n");
    for( i = 0; i < out_nnz_ref; i++ ) {
        printf("(%" PRId64 ", %3.3f)\n", out_indices_ref[i], out_values_ref[i]);
    }
#endif

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

    // Create mkl_graph_vectors representing x, y and the mask.
    // NOTE: The mkl_graph_vector objects assume ownership of these arrays for
    // editing, but does NOT assume responsibility for freeing memory.
    CALL_AND_CHECK_STATUS(mkl_graph_vector_create(&mvec),
                          "Error after calling mkl_graph_vector_create for the mask");
    CALL_AND_CHECK_STATUS(mkl_graph_vector_set_sparse(mvec, nrows, m_nnz_sparse,
                                                      m_indices_sparse, MKL_GRAPH_TYPE_INT64,
                                                      m_values_sparse, MKL_GRAPH_TYPE_INT32),
                          "Error after calling mkl_graph_vector_set_sparse "
                          "for the mask");

    CALL_AND_CHECK_STATUS(mkl_graph_vector_create(&xvec),
                          "Error after calling mkl_graph_vector_create for x");
    CALL_AND_CHECK_STATUS(mkl_graph_vector_set_sparse(xvec, ncols, x_nnz_sparse,
                                                      x_indices_sparse, MKL_GRAPH_TYPE_INT64,
                                                      x_values_sparse, MKL_GRAPH_TYPE_FP32),
                          "Error after calling mkl_graph_vector_set_sparse "
                          "for x");

    CALL_AND_CHECK_STATUS(mkl_graph_vector_create(&yvec),
                          "Error after calling mkl_graph_vector_create for y");

    // Compute y<mask> = A * x.
    CALL_AND_CHECK_STATUS(mkl_graph_mxv(yvec, mvec, MKL_GRAPH_ACCUMULATOR_NONE,
                                        MKL_GRAPH_SEMIRING_PLUS_TIMES_FP32,
                                        Amat, xvec, NULL,
                                        MKL_GRAPH_REQUEST_COMPUTE_ALL,
                                        MKL_GRAPH_METHOD_AUTO),
                          "Error after calling mkl_graph_mxv to compute "
                          "y<mask> = A * x");

    // Export the values of y.
    // NOTE: The pointer yout_values can be safely used to access the data but
    // cannot be safely used for editing. It's possible but not guaranteed that
    // this pointer will be equal to the input y_values.
    CALL_AND_CHECK_STATUS(mkl_graph_vector_get_sparse(yvec, &ydim, &out_nnz,
                                                     (void**)&out_indices, &yout_indices_type,
                                                     (void**)&out_values,  &yout_values_type),
                          "Error after calling mkl_graph_vector_get_sparse "
                          "for y");

    // Validation: comparing out and out_ref
    if (out_nnz_ref != out_nnz) {
        printf("[ERROR] out_nnz_ref = %" PRId64 " != %" PRId64 " = out_nnz \n", out_nnz_ref, out_nnz);
        error = -1;
        goto memory_free;
    }

    if (out_values == NULL || out_indices == NULL) {
        printf("[ERROR] One of returned pointers is NULL\n");
        error = -1;
        goto memory_free;
    }

    for (i = 0; i < out_nnz; i++) {
        if (out_indices[i] != out_indices_ref[i]) {
            error = -1;
            printf("[ERROR] out_indices[%" PRId64 "] = %" PRId64 " != %" PRId64 " = out_indices_ref[%" PRId64 "] \n", i, out_indices[i], out_indices_ref[i], i);
        } else {
            printf("out_indices[%" PRId64 "] = %" PRId64 ", out_indices_ref[%" PRId64 "] = %" PRId64 "\n", i, out_indices[i], i, out_indices_ref[i]);
        }
    }

    double error_l2 = 0.0, error_max = 0.0;
    const float float_tol = 1.0e-5;
    for (i = 0; i < out_nnz; i++) {
        double tmp = fabs(out_values[i] - out_values_ref[i]);
        error_l2 += tmp * tmp;
        if ( tmp > error_max )
            error_max = tmp;
        printf("out_values[%" PRId64 "] = %3.3f, out_values_ref[%" PRId64 "] = %3.3f\n", i, out_values[i], i, out_values_ref[i]);
    }
    error_l2 = sqrt(error_l2 / out_nnz);

    if (error_l2 < float_tol && error_max < float_tol)
        printf("Success, validation PASSED!\n");
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
    mkl_graph_vector_destroy(&mvec);
    mkl_graph_vector_destroy(&yvec);
    // Free user-allocated data. mkl_graph_<object>_destroy leaves
    // user-allocated memory intact.
    mkl_free(A_rows_start);
    mkl_free(A_col_indx);
    mkl_free(A_values);
    mkl_free(m_indices_sparse);
    mkl_free(m_values_sparse);
    mkl_free(x_indices_sparse);
    mkl_free(x_values_sparse);
    mkl_free(out_indices_ref);
    mkl_free(out_values_ref);

    return error;
}

