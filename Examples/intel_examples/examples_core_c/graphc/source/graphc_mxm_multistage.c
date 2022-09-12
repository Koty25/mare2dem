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

// hard-coded the sum of the output values and the number of entries in the output
// (used for validation)
#define OUTSUM 170
#define OUTNNZ 5


//******************************************************************************
// Task:    Compute non-masked matrix-matrix product C = L * L
//          with a standard PLUS_TIMES semiring and allocate memory for the
//          output matrix on the side of the user.
//
//           |   0   0   0   0   0   0   |
//           |   1   0   0   0   0   0   |
//   L     = |   0   0   0   0   0   0   | (in a dense form)
//           |   2   0   3   0   0   0   |
//           |   4   5   0   6   0   0   |
//           |   7   0   8   0   9   0   |
//
// We represent L in CSR format with the arrays L_values, L_rows_start,
// L_col_indx as follows,
//
//   L_values     = {1, 2, 3, 4, 5, 6, 7, 8, 9};
//   L_col_indx   = {0, 0, 2, 0, 1, 3, 0, 2, 4};
//   L_rows_start = {0, 0, 1, 1, 3, 6, 9};
//
//******************************************************************************

/* To avoid constantly repeating the part of code that checks inbound graph functions' status,
   use macro CALL_AND_CHECK_STATUS */
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != MKL_GRAPH_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          error = 1;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)

int main()
{
    int64_t i;
    int32_t error = 0;

    // Arrays for the CSR representation of the input matrix
    int32_t *L_values     = NULL;
    int64_t *L_rows_start = NULL;
    int64_t *L_col_indx   = NULL;

    // Temporary output variables and accessors for exported data
    int64_t  outnnz = 0;
    int64_t  outsum = 0;

    // Arrays for the CSR representation of the output matrix
    int32_t *C_values     = NULL;
    int64_t *C_rows_start = NULL;
    int64_t *C_col_indx   = NULL;

    mkl_graph_matrix_t Lmat = NULL, Cmat = NULL;

    // Steps necessary for checking for the memory leaks
    int AllocatedBuffers = 0;
    MKL_INT64 occupied_memory, leaked_memory;
    mkl_free_buffers();
    occupied_memory = mkl_mem_stat(&AllocatedBuffers);

    L_rows_start = (int64_t*)mkl_malloc((M + 1)*sizeof(int64_t), ALIGN);
    L_col_indx   = (int64_t*)mkl_malloc(  NNZ  *sizeof(int64_t), ALIGN);
    L_values     = (int32_t*)mkl_malloc(  NNZ  *sizeof(int32_t), ALIGN);

    if (L_rows_start == NULL || L_col_indx == NULL || L_values == NULL) {
        printf("Cannot allocate data for L \n");
        error = -1;
        goto memory_free;
    }

    // CSR arrays for L = strict lower triangle of A
    for( i = 0; i < NNZ; i++ )
        L_values[i] = 1 + i;

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
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_create(&Lmat),"Error after calling mkl_graph_matrix_create for L");

    // Set CSR arrays for the input matrix in the graph matrix
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_set_csr(Lmat, M, M, L_rows_start, MKL_GRAPH_TYPE_INT64, L_col_indx, MKL_GRAPH_TYPE_INT64, L_values, MKL_GRAPH_TYPE_INT32),
                          "Error after calling mkl_graph_matrix_set_csr for L");

    // Create an mkl_graph_matrix for the output matrix
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_create(&Cmat),"Error after calling mkl_graph_matrix_create for C");

    // Now multistage workflow is ready to start

    // 1. First, allocate memory for the output matrix. In this example, the user
    // sets CSR format for the output.
    // NOTE: The user must decide the output format for the matrix in order to
    // allocate enough memory for rows_start/cols_start.
    C_rows_start = (int64_t*)mkl_malloc((M + 1)* sizeof(int64_t), ALIGN);
    if (C_rows_start == NULL) {
        printf("Cannot allocate data for C_rows_start\n");
        error = -1;
        goto memory_free;
    }

    // 2. Set rows_start array for the output matrix in the output graph matrix
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_set_csr(Cmat, M, M, C_rows_start, MKL_GRAPH_TYPE_INT64, NULL, MKL_GRAPH_TYPE_INT64, NULL, MKL_GRAPH_TYPE_INT32),
                          "Error after calling mkl_graph_matrix_set_csr for C with ia only");

    // 3. First stage of non-masked mxm: computing the number of entries in the output.
    // NOTE: In multistage execution, the memory for the output arrays should be allocated by the user.
    // This means that the user must have a memory buffer for rows_start/ (in case of CSR) set in C
    // prior to calling this stage.
    CALL_AND_CHECK_STATUS(mkl_graph_mxm(Cmat, NULL, MKL_GRAPH_ACCUMULATOR_NONE, MKL_GRAPH_SEMIRING_PLUS_TIMES_INT32, Lmat, Lmat, NULL,
                                        MKL_GRAPH_REQUEST_FILL_NNZ, MKL_GRAPH_METHOD_GUSTAVSON),
                          "Error after calling mkl_graph_mxm to compute C = L * L with request FILL_NNZ");

    // 4. Requesting the computed nnz in the output matrix
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_get_property(Cmat, MKL_GRAPH_PROPERTY_NNZ, &outnnz),
                          "Error after calling mkl_graph_matrix_get_property for C");
    printf("[INFO] outnnz computed in the first stage = %" PRId64 " \n", outnnz);

    // 5. Allocating arrays for column indices and values of the output matrix
    C_col_indx = (int64_t*)mkl_malloc(outnnz * sizeof(int64_t), ALIGN);
    C_values   = (int32_t*)mkl_malloc(outnnz * sizeof(int32_t), ALIGN);
    if (C_col_indx == NULL || C_values == NULL) {
        printf("Cannot allocate data for C_col_indx and/or C_values\n");
        error = -1;
        goto memory_free;
    }

    // 6. Setting allocated buffers for column indices and output valus in the output graph matrix
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_set_csr(Cmat, M, M, C_rows_start, MKL_GRAPH_TYPE_INT64, C_col_indx, MKL_GRAPH_TYPE_INT64, C_values, MKL_GRAPH_TYPE_INT32),
                          "Error after calling mkl_graph_matrix_set_csr for C with ja and values");

    // 7. Second stage: computing entries (column indices and values) of the output matrix
    // NOTE: Again, in multistage execution, the memory for the output arrays should be allocated by the user.
    // This means that the user must have memory buffers for col_indx and values (in case of CSR) set in C
    // prior to calling this stage
    CALL_AND_CHECK_STATUS(mkl_graph_mxm(Cmat, NULL, MKL_GRAPH_ACCUMULATOR_NONE, MKL_GRAPH_SEMIRING_PLUS_TIMES_INT32, Lmat, Lmat, NULL,
                                        MKL_GRAPH_REQUEST_FILL_ENTRIES, MKL_GRAPH_METHOD_GUSTAVSON),
                          "Error after calling mkl_graph_mxm to compute C = L * L with request FILL_ENTRIES");

    // Now the result of non-masked mxm is ready, the user can query its properties and export the data if needed
    CALL_AND_CHECK_STATUS(mkl_graph_matrix_get_property(Cmat, MKL_GRAPH_PROPERTY_NNZ, &outnnz),
                          "Error after calling mkl_graph_matrix_get_property for C");
    printf("After multistage computations, outnnz = %" PRId64 " \n", outnnz);

    // NOTE: The user could export data from the output matrix via  mkl_graph_matrix_get_<format> routine.
    // However, for multistage exeuction it is not needed since the data for the output matrix is allocated
    // on the side of the user

    for (i = 0; i < M + 1; i++)
        printf("rows_start[%" PRId64 "] = %" PRId64 " \n", i, C_rows_start[i]);
    for( i = 0; i < outnnz; i++) {
        printf("col_indx[%" PRId64 "] = %" PRId64 " values[%" PRId64 "] = %d \n", i, C_col_indx[i], i, C_values[i]);
        outsum += C_values[i];
    }

    /*
        Correct output for non-masked L * L with PLUS_TIMES_INT32 semiring:
        row 0:
        row 1:
        row 2:
        row 3:
        row 4: (0, 17.0)           (2, 18.0)
        row 5: (0, 36.0) (1, 45.0)           (3, 54.0)
    */


    printf("Validation (short, only the sum of the values and the number of entries in the output):\n");
    // Check that the sum of the output entries (number of triangles for the masked case) is correct
    // and that the number of nonzeros in the output is correct (incl. the zombie removal)
    if (outsum == OUTSUM && outnnz == OUTNNZ)
        printf("Success, test PASSED!\n");
    else {
        printf("Validation FAILED:\n");
        printf("Is outsum = %" PRId64 " equal to %d? Check!\n", outsum, OUTSUM);
        printf("Is outnnz = %" PRId64 " equal to %d? Check!\n", outnnz, OUTNNZ);
        error = -1;
        goto memory_free;
    }

    // Releasing the memory
memory_free:
    mkl_graph_matrix_destroy(&Cmat);
    mkl_graph_matrix_destroy(&Lmat);
    mkl_free(L_rows_start);
    mkl_free(L_col_indx);
    mkl_free(L_values);
    mkl_free(C_rows_start);
    mkl_free(C_col_indx);
    mkl_free(C_values);

    mkl_free_buffers();
    leaked_memory = mkl_mem_stat(&AllocatedBuffers) - occupied_memory;
    printf("\n[CHECK MEMORY] leaked_memory = %" PRId64 "\n", leaked_memory); fflush(0);

    if (leaked_memory != 0) {
        printf("[ERROR]: Memory has leaked, test will be marked as failed\n");
        error = -1;
    }

    return error;
}

