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

#define M 6
#define NNZ 9
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
// Task:    Count triangles in the graph G represented by adjacency
//          matrix A using mkl_graph_mxm.
//
//           |   0   1   0   1   1   1   |
//           |   1   0   0   0   1   0   |
//   A     = |   0   0   0   1   0   1   |
//           |   1   0   1   0   1   0   |
//           |   1   1   0   1   0   1   |
//           |   1   0   1   0   1   0   |
//
// For triangle counting we use the strict lower triangular portion of A,
// denoted by L:
//
//           |   0   0   0   0   0   0   |
//           |   1   0   0   0   0   0   |
//   L     = |   0   0   0   0   0   0   |
//           |   1   0   1   0   0   0   |
//           |   1   1   0   1   0   0   |
//           |   1   0   1   0   1   0   |
//
// We represent L in CSR format with the arrays L_values, L_rows_start,
// L_col_indx as follows,
//
//   L_values = {1, 1, 1, 1, 1, 1, 1, 1, 1},
//   L_col_indx = {0, 0, 2, 0, 1, 3, 0, 2, 4},
//   L_rows_start = {0, 0, 1, 1, 3, 6, 9}.
//
//******************************************************************************

int main()
{
    int64_t i;
    int32_t error = 0;

    int32_t *L_values = NULL;
    int64_t *L_rows_start = NULL;
    int64_t *L_col_indx = NULL;
    int64_t nrows = M;
    int64_t ncols = M;

    int64_t  nnz_exp;
    int32_t *vals_exp;

    mkl_graph_matrix_t Lmat = NULL;
    mkl_graph_matrix_t Cmat = NULL;

    L_rows_start = (int64_t*)mkl_malloc((M + 1)*sizeof(int64_t), ALIGN);
    L_col_indx   = (int64_t*)mkl_malloc(  NNZ  *sizeof(int64_t), ALIGN);
    L_values     = (int32_t*)mkl_malloc(  NNZ  *sizeof(int32_t), ALIGN);

    if (L_rows_start == NULL || L_col_indx == NULL || L_values == NULL) {
        printf("Cannot allocate data \n");
        error = -1;
        goto memory_free;
    }

    // Construct CSR arrays for L = strict lower triangle of A.
    for( i = 0; i < NNZ; i++ )
        L_values[i] = 1;

    // row 0
    // row 1
    L_col_indx[0] = 0;
    // row 2
    // row 3
    L_col_indx[1] = 0;
    L_col_indx[2] = 2;
    // row 4
    L_col_indx[3] = 0;
    L_col_indx[4] = 1;
    L_col_indx[5] = 3;
    // row 5
    L_col_indx[6] = 0;
    L_col_indx[7] = 2;
    L_col_indx[8] = 4;

    L_rows_start[0] = 0;
    L_rows_start[1] = L_rows_start[0] + 0;
    L_rows_start[2] = L_rows_start[1] + 1;
    L_rows_start[3] = L_rows_start[2] + 0;
    L_rows_start[4] = L_rows_start[3] + 2;
    L_rows_start[5] = L_rows_start[4] + 3;
    L_rows_start[6] = L_rows_start[5] + 3;

    // Create an mkl_graph_matrix representing L using L_rows_start, L_col_indx,
    // and L_values.
    // NOTE: The mkl_graph_matrix assumes ownership of these arrays for editing,
    // but does NOT assume responsibility for freeing memory.
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_create(&Lmat),
                          "Error after calling mkl_graph_matrix_create for L");

    CALL_AND_CHECK_STATUS(mkl_graph_matrix_set_csr(Lmat, nrows, ncols,
                                                   L_rows_start,
                                                   MKL_GRAPH_TYPE_INT64,
                                                   L_col_indx,
                                                   MKL_GRAPH_TYPE_INT64,
                                                   L_values,
                                                   MKL_GRAPH_TYPE_INT32),
                          "Error after calling mkl_graph_matrix_set_csr for L");

    // For applications with many matrix-matrix multiplies, calling
    // mkl_graph_optimize_mxm before use can decrease overall runtime.
    // CALL_AND_CHECK_STATUS(mkl_graph_optimize_mxm(Lmat, MKL_GRAPH_SEMIRING_PLUS_TIMES_INT32,
    //                                              Lmat, Lmat, NULL, ncalls),
    //                       "Error after calling mkl_graph_optimize_mxm");

    // Create an mkl_graph_matrix to store the result.
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_create(&Cmat),
                          "Error after calling mkl_graph_matrix_create for C");

    // Compute C<L> = L * L using mkl_graph_mxm, with NULL passed as the
    // descriptor.
    CALL_AND_CHECK_STATUS(mkl_graph_mxm(Cmat, Lmat, MKL_GRAPH_ACCUMULATOR_NONE,
                                        MKL_GRAPH_SEMIRING_PLUS_TIMES_INT32,
                                        Lmat, Lmat, NULL,
                                        MKL_GRAPH_REQUEST_COMPUTE_ALL,
                                        MKL_GRAPH_METHOD_AUTO),
                          "Error after calling mkl_graph_mxm to compute "
                          "C<L> = L * L");

    // Export number of nonzeros from Cmat.
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_get_property(Cmat,
                                                        MKL_GRAPH_PROPERTY_NNZ,
                                                        &nnz_exp),
                          "Error after calling mkl_graph_matrix_get_property "
                          "for C");
    printf("nnz_exp = %" PRId64 " \n", nnz_exp);

    // Export values from Cmat.
    // NOTE: Any unneeded arguments can be safely passed as NULL.
    // NOTE: The pointer vals_exp can be safely used to access the data but
    // cannot be safely used for editing.
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_get_csr(Cmat, NULL, NULL, NULL, NULL,
                                                   NULL, NULL,
                                                   (void**)&vals_exp, NULL),
                          "Error after calling mkl_graph_matrix_get_csr for C");

    // Use the exported values array to count the triangles.
    int64_t ntriangles = 0;
    for( i = 0; i < nnz_exp; i++) {
        printf("vals_exp[%" PRId64 "] = %d \n", i, vals_exp[i]);
        ntriangles += vals_exp[i];
    }

    printf("Validation:\n");
    if (ntriangles == 3)
        printf("Success, test PASSED!\n");
    else {
        printf("Validation FAILED: ntriangles = %" PRId64 " while there should "
               "be 3!\n", ntriangles);
        error = -1;
        goto memory_free;
    }

memory_free:
    // Cleanup workflow should destroy all mkl_graph objects and free all
    // user-allocated memory used as input to mkl_graph objects. Access pointers
    // like vals_exp that are set by calls to
    // mkl_graph_<object>_get_<representation> do not need to be freed.

    // Calls to mkl_graph_matrix_destroy de-allocate all internal data.
    mkl_graph_matrix_destroy(&Lmat);
    mkl_graph_matrix_destroy(&Cmat);
    // Free user-allocated data. mkl_graph_matrix_destroy leaves user-allocated
    // memory intact.
    mkl_free(L_rows_start);
    mkl_free(L_col_indx);
    mkl_free(L_values);

    return error;
}

