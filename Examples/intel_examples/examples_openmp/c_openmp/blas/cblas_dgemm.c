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
!  Content:
!      oneMKL CBLAS_DGEMM OpenMP offload Example Program Text 
!******************************************************************************/

#include <stdio.h>
#include <omp.h>
#include "mkl.h"
#include "mkl_omp_offload.h"
#include "common.h"

int dnum = 0;

int main() {

    double *a, *b, *c, *c_ref, alpha, beta;
    MKL_INT m, n, k, lda, ldb, ldc, i, j;
    MKL_INT sizea, sizeb, sizec;
    
    alpha = 1.0;
    beta = 1.0;

    m = 1449;
    n = 1120;
    k = 1083;

    lda = m;
    ldb = k;
    ldc = m;

    sizea = lda * k;
    sizeb = ldb * n;
    sizec = ldc * n;
    
    // allocate matrices
    a = (double *)mkl_malloc((lda * k) * sizeof(double), 64);
    b = (double *)mkl_malloc(ldb * n * sizeof(double), 64);
    c = (double *)mkl_malloc(ldc * n * sizeof(double), 64);
    c_ref = (double *)mkl_malloc(ldc * n * sizeof(double), 64);

    if ((a == NULL) || (b == NULL) || (c == NULL)) {
        printf("Cannot allocate matrices\n");
        return 1;
    }

    // initialize matrices
    init_double_array(lda * k, a, 1);
    init_double_array(ldb * n, b, 1);

#pragma omp target map(c[0:sizec], c_ref[0:sizec])
    {
        for (i = 0; i < sizec; i++) {
            c[i] = 42;
            c_ref[i] = 42;
        }
    }
    
    MKL_INT bound_m = (m > 10) ? 10 : m;
    MKL_INT bound_n = (n > 10) ? 10 : n;
    
    // run gemm on host, use standard oneMKL interface   
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a, lda, b, ldb, beta, c_ref, ldc);

#pragma omp target data map(to:a[0:sizea],b[0:sizeb]) map(tofrom:c[0:sizec]) device(dnum)
    {

        // run gemm on gpu, use standard oneMKL interface within a variant dispatch construct
#pragma omp target variant dispatch device(dnum) use_device_ptr(a, b, c)
        {
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        }
        
    }
    
    double real;
    int err = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            real = c[i + ldc * j] - c_ref[i + ldc * j];
            real = (real > 0) ? real : -real;
            if (real > 0.0001) {
#ifdef MKL_ILP64
                printf("c[%lld][%lld] != c_ref[%lld][%lld], computed value is %f, reference value is %f, difference is %f\n",
                       i, j, i, j, c[i + ldc * j], c_ref[i + ldc * j], real);
#else
                printf("c[%d][%d] != c_ref[%d][%d], computed value is %f, reference value is %f, difference is %f\n",
                       i, j, i, j, c[i + ldc * j], c_ref[i + ldc * j], real);
#endif
                mkl_free(a);
                mkl_free(b);
                mkl_free(c);
                mkl_free(c_ref);
                return 1;
            }
        }
    }
 
    printf("Upper left corner of the C matrix:\n");
    printf("C matrix:\n");
    for (i = 0; i < bound_m; i++) {
        for (j = 0; j < bound_n; j++) {
            printf("%f\t", c[i + ldc * j]);
        }
        printf("\n");
    }

    printf("Reference matrix:\n");
    for (i = 0; i < bound_m; i++) {
        for (j = 0; j < bound_n; j++) {
            printf("%f\t", c_ref[i + ldc * j]);
        }
        printf("\n");
    }
   
    mkl_free(a);
    mkl_free(b);
    mkl_free(c);
    mkl_free(c_ref);
    return 0;
}
