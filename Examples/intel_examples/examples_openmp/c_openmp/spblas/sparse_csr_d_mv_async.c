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

/*
*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Sparse BLAS C OpenMP
*             offload example for mkl_sparse_d_mv with async execution.
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
*         columns =  ( 1  2  4  1  2  3  4  5  1  3  4  2  5 )
*         rowIndex = ( 1        4     6        9       12    14 )
*
*  The test computes the following operations :
*
*       A*x = y using mkl_sparse_d_mv omp offload with async execution
*       where A is a general sparse matrix and x and y are vectors
*
********************************************************************************
*/
#include <assert.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "mkl.h"
#include "mkl_omp_offload.h"

#define FREE(p)       \
    do {              \
        if (p) {      \
            free(p);  \
            p = NULL; \
        }             \
    } while (0)

int validation_result(double *ref, double *z, int length, int flps_per_value);

int main()
{
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
    double values[NNZ]    = { 1.0, -1.0,     -3.0,
                             -2.0,  5.0,
                                         4.0, 6.0, 4.0,
                             -4.0,       2.0, 7.0,
                                    8.0,          -5.0 };
    MKL_INT colIndex[NNZ] = { 1,      2,        4,
                              1,      2,
                                           3,   4,   5,
                              1,           3,   4,
                                      2,             5 };
    MKL_INT rowStart[N + 1] = {1, 4, 6, 9, 12, 14};
    // Descriptor of main sparse matrix properties
    struct matrix_descr descrA;
    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    // // Structure with sparse matrix stored in CSR format
    sparse_matrix_t csrA;
    //*******************************************************************************
    //    Declaration of local variables:
    //*******************************************************************************
    double x_array[M] = {1.0, 5.0, 1.0, 4.0, 1.0};
    double alpha = 1.0, beta = 0.0;
    MKL_INT i;

    double *a   = (double *)malloc(sizeof(double) * NNZ);
    MKL_INT *ja = (MKL_INT *)malloc(sizeof(MKL_INT) * NNZ);
    MKL_INT *ia = (MKL_INT *)malloc(sizeof(MKL_INT) * (N + 1));
    double *x   = (double *)malloc(sizeof(double) * M);

    double *y = (double *)malloc(sizeof(double) * N);
    double *z = (double *)malloc(sizeof(double) * N);

    if (!a || !ja || !ia || !x || !y || !z) {
        FREE(ia);
        FREE(ja);
        FREE(a);
        FREE(x);
        FREE(y);
        FREE(z);
        return 1;
    }

    for (i = 0; i < NNZ; i++) {
        a[i]  = values[i];
        ja[i] = colIndex[i];
    }
    for (i = 0; i < N + 1; i++) {
        ia[i] = rowStart[i];
    }

    for (i   = 0; i < M; i++)
        x[i] = x_array[i];
    for (i = 0; i < N; i++) {
        y[i] = 0.0;
        z[i] = 0.0;
    }

    printf("\n EXAMPLE PROGRAM FOR mkl_sparse_d_mv omp_offload async\n");
    printf("---------------------------------------------------\n");
    printf("\n");
    printf("   INPUT DATA FOR mkl_sparse_d_mv omp offload async   \n");
    printf("   WITH GENERAL SPARSE MATRIX     \n");
    printf("   ALPHA = %4.1f  BETA = %4.1f    \n", alpha, beta);
    printf("   SPARSE_OPERATION_NON_TRANSPOSE \n");
    printf("   Input vector                   \n");
    for (i = 0; i < M; i++) {
        printf("%7.1f\n", x[i]);
    };

    // Create handle with matrix stored in CSR format
    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ONE,
                            N, // number of rows
                            M, // number of cols
                            ia, ia + 1, ja, a);

    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    // Analyze sparse matrix; choose proper kernels and workload balancing strategy
    mkl_sparse_optimize(csrA);

    // Compute y = alpha * A * x + beta * y
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, csrA, descrA, x, beta, y);

    // Release matrix handle and deallocate matrix
    mkl_sparse_destroy(csrA);

    printf("                                   \n");
    printf("   OUTPUT DATA FOR mkl_sparse_d_mv \n");

    // y should be equal { -16.0, 23.0, 32.0, 26.0, 35.0 }
    for (i = 0; i < N; i++) {
        printf("%7.1f\n", y[i]);
    };

    printf("---------------------------------------------------\n");
    fflush(stdout);

    const int devNum = 0;

    sparse_matrix_t csrA_gpu;

    sparse_status_t status_gpu1;
    sparse_status_t status_gpu2;
    sparse_status_t status_gpu3;

// call create_csr/mv/destroy via omp_offload.
#pragma omp target data map(to:ia[0:N+1],ja[0:NNZ],a[0:NNZ],x[0:M]) map(tofrom:z[0:N]) device(devNum)
    {
        printf("Create CSR matrix via omp_offload\n");

#pragma omp target variant dispatch device(devNum) use_device_ptr(ia, ja, a) nowait
        {
            status_gpu1 = mkl_sparse_d_create_csr(&csrA_gpu, SPARSE_INDEX_BASE_ONE, N, M, ia,
                                                  ia + 1, ja, a);
        }

        printf("Compute mkl_sparse_d_mv via omp_offload\n");

#pragma omp target variant dispatch device(devNum) use_device_ptr(x, z) nowait
        {
            status_gpu2 = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, csrA_gpu, descrA,
                                          x, beta, z);
        }

        printf("Destroy the CSR matrix via omp_offload\n");

#pragma omp target variant dispatch device(devNum) nowait
        {
            status_gpu3 = mkl_sparse_destroy(csrA_gpu);
        }
    }

#pragma omp taskwait

    int flps_per_val = 2 * ((NNZ / N) + 1);
    int status       = 0;

    if (status_gpu1 != SPARSE_STATUS_SUCCESS)
        goto failed;
    if (status_gpu2 != SPARSE_STATUS_SUCCESS)
        goto failed;
    if (status_gpu3 != SPARSE_STATUS_SUCCESS)
        goto failed;

    printf("   OUTPUT DATA FOR mkl_sparse_d_mv_omp_offload async execution.\n");
    // z should be equal { -16.0, 23.0, 32.0, 26.0, 35.0 }
    for (i = 0; i < N; i++) {
        printf("%7.1f\n", z[i]);
    };
    printf("---------------------------------------------------\n");

    status = validation_result(y, z, N, flps_per_val);

cleanup:
    FREE(ia);
    FREE(ja);
    FREE(a);
    FREE(x);
    FREE(y);
    FREE(z);

    {
        const int statusAll = status | status_gpu1 | status_gpu2 | status_gpu3;
        printf("Test %s\n", statusAll == 0 ? "PASSED" : "FAILED");
        return statusAll;
    }

failed:
    printf("\tERROR: status_gpu1 = %d, status_gpu2 = %d, status_gpu3 = %d\n", status_gpu1,
           status_gpu2, status_gpu3);
    goto cleanup;
}

int validation_result(double *ref, double *z, int length, int flps_per_value)
{
    const double err_thresh = 1.0 * flps_per_value * DBL_EPSILON;
    printf(" Check if error is below err_thresh %.3lg\n", err_thresh);

    double maxerr = 0.0;
    for (int i = 0; i < length; i++) {
        double err = fabs(ref[i] - z[i]);
        if (err > maxerr)
            maxerr = err;
        if (!(maxerr < err_thresh)) {
            printf(" z[%d]: expected %.7g, got %.7g, err %.3lg\n", i, ref[i], z[i], err);
            printf(" Verification FAILED\n");
            return 1;
        }
    }
    printf(" Verified, maximum error was %.3lg\n", maxerr);

    return 0;
}
