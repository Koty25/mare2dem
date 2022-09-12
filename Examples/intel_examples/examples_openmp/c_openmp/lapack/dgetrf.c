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
*  Content:
*      dgetrf OpenMP Offload Example
*******************************************************************************/
#include <stdio.h>
#include <omp.h>
#include "mkl.h"
#include "mkl_omp_offload.h"

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

/* Default GPU device */
int dnum = 0;

int main() {

    if (sizeof(MKL_INT) != 8) {
        printf("ILP-64 interface required\n");
        return 1;
    }

    /* Target region necessary for initialization of OpenMP Offloading, see
     * documentation for further information. */
    MKL_INT *temp = (MKL_INT *)mkl_malloc(sizeof(MKL_INT), 4096);
    #pragma omp target map(temp[0 : 1])
    {}
    mkl_free(temp);

    MKL_INT m = 3;
    MKL_INT n = 3;
    MKL_INT lda = 4;

    MKL_INT array_size = lda * n;
    MKL_INT ipiv_size = MIN(m, n);

    double *A = (double *)mkl_malloc(sizeof(double) * array_size, 4096);
    MKL_INT *ipiv = (MKL_INT *)mkl_calloc(1, sizeof(MKL_INT) * ipiv_size, 4096);
    MKL_INT *info = (MKL_INT *)mkl_calloc(1, sizeof(MKL_INT), 4096);

    if (!A || !ipiv || !info) {
        printf("could not allocate host memory\n");
        return 2;
    }

    double init[] = {1, 4, 3, 1, 3, 5, 1, -1, 3};
    for (MKL_INT row = 0; row < m; row++) {
        for (MKL_INT col = 0; col < n; col++) {
            A[row + col * lda] = init[row + col * row];
        }
    }

    /* Execute getrf on GPU via variant dispatch construct */
    #pragma omp target data map(A[0 : array_size], ipiv[0 : ipiv_size], info[0 : 1]) device(dnum)
    {
        #pragma omp target variant dispatch device(dnum) use_device_ptr(A, ipiv, info)
        dgetrf(&m, &n, A, &lda, ipiv, info);
    }
    if (*info) {
        printf("dgetrf offload failed with info = %lld\n", *info);
    }

    printf("Matrix A Input\n");
    for (MKL_INT row = 0; row < m; row++) {
        for (MKL_INT col = 0; col < n; col++) {
            printf("%5.2f ", init[row + col * m]);
        }
        printf("\n");
    }
    printf("\n");
    printf("Matrix A Output\n");
    for (MKL_INT row = 0; row < m; row++) {
        for (MKL_INT col = 0; col < n; col++) {
            printf("%5.2f ", A[row + col * lda]);
        }
        printf("\n");
    }

    mkl_free(A);
    mkl_free(ipiv);
    mkl_free(info);
    return 0;
}
