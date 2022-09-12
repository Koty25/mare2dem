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
!      oneMKL ZGEMM OpenMP offload Example Program Text 
!******************************************************************************/

#include <stdio.h>
#include <omp.h>
#include "mkl.h"
#include "mkl_omp_offload.h"
#include "common.h"

int dnum = 0;

int main() {

    MKL_Complex16 *a, *b, *c, *c_ref, alpha, beta;
    MKL_INT m, n, k, lda, ldb, ldc, i, j;
    MKL_INT sizea, sizeb, sizec;
    
    alpha.real = 1.0; alpha.imag = 0.0;
    beta.real = 1.0; beta.imag = 0.0;

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
    a = (MKL_Complex16 *)mkl_malloc((lda * k) * sizeof(MKL_Complex16), 64);
    b = (MKL_Complex16 *)mkl_malloc(ldb * n * sizeof(MKL_Complex16), 64);
    c = (MKL_Complex16 *)mkl_malloc(ldc * n * sizeof(MKL_Complex16), 64);
    c_ref = (MKL_Complex16 *)mkl_malloc(ldc * n * sizeof(MKL_Complex16), 64);

    if ((a == NULL) || (b == NULL) || (c == NULL)) {
        printf("Cannot allocate matrices\n");
        return 1;
    }

    // initialize matrices
    init_double_complex_array(lda * k, a, 1);
    init_double_complex_array(ldb * n, b, 1);

#pragma omp target map(c[0:sizec], c_ref[0:sizec])
    {
        for (i = 0; i < sizec; i++) {
            c[i].real = 0.42;
            c[i].imag = 0.42;
            c_ref[i].real = 0.42;
            c_ref[i].imag = 0.42;
        }
    }
    
    MKL_INT bound_m = (m > 10) ? 10 : m;
    MKL_INT bound_n = (n > 10) ? 10 : n;
    
    // run gemm on host, use standard oneMKL interface
    zgemm("N", "N", &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c_ref, &ldc);

#pragma omp target data map(to:a[0:sizea],b[0:sizeb]) map(tofrom:c[0:sizec]) device(dnum)
    {

        // run gemm on gpu, use standard oneMKL interface within a variant dispatch construct
#pragma omp target variant dispatch device(dnum) use_device_ptr(a, b, c)
        {
            zgemm("N", "N", &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
        }
        
    }
    
    double real, imag;
    int err = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            real = c[i + ldc * j].real - c_ref[i + ldc * j].real;
            real = (real > 0) ? real : -real;
            imag = c[i + ldc * j].imag - c_ref[i + ldc * j].imag;
            imag = (imag > 0) ? imag : -imag;
            
            if ((real + imag) > 0.0001) {
#ifdef MKL_ILP64
                printf("c[%lld][%lld] != c_ref[%lld][%lld], computed value is %lf + i%lf, reference value is %lf + i%lf, difference is %lf + i%lf\n",
                       i, j, i, j, c[i + ldc * j].real, c[i + ldc * j].imag, c_ref[i + ldc * j].real, c_ref[i + ldc * j].imag, real, imag);
#else
                printf("c[%d][%d] != c_ref[%d][%d], computed value is %lf + i%lf, reference value is %lf + i%lf, difference is %lf + i%lf\n",
                       i, j, i, j, c[i + ldc * j].real, c[i + ldc * j].imag, c_ref[i + ldc * j].real, c_ref[i + ldc * j].imag, real, imag);
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
            printf("%lf + i%lf\t", c[i + ldc * j].real, c[i + ldc * j].imag);
        }
        printf("\n");
    }

    printf("Reference matrix:\n");
    for (i = 0; i < bound_m; i++) {
        for (j = 0; j < bound_n; j++) {
            printf("%lf + i%lf\t", c_ref[i + ldc * j].real, c_ref[i + ldc * j].imag);
        }
        printf("\n");
    }
   
    mkl_free(a);
    mkl_free(b);
    mkl_free(c);
    mkl_free(c_ref);
    return 0;
}
