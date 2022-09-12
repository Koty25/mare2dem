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
!      oneMKL SGER OpenMP offload Example Program Text 
!******************************************************************************/

#include <stdio.h>
#include <omp.h>
#include "mkl.h"
#include "mkl_omp_offload.h"
#include "common.h"

int dnum = 0;

int main() {

    float *a, *x, *y, *a_ref, alpha;
    MKL_INT m, n, lda, incx, incy, i, j;
    MKL_INT sizea, sizex, sizey;
    
    alpha = 1.0;

    m = 1449;
    n = 1120;

    lda = m;
    incx = 1;
    incy = 1;

    sizea = lda * n;
    sizex = m;
    sizey = n;
    
    // allocate matrices
    a = (float *)mkl_malloc(sizea * sizeof(float), 64);
    x = (float *)mkl_malloc(sizex * sizeof(float), 64);
    y = (float *)mkl_malloc(sizey * sizeof(float), 64);
    a_ref = (float *)mkl_malloc(sizea * sizeof(float), 64);

    if ((a == NULL) || (x == NULL) || (y == NULL) || (a_ref == NULL)) {
        printf("Cannot allocate matrices\n");
        return 1;
    }

    // initialize matrices
    init_single_array(sizey, y, 1);
    init_single_array(sizex, x, 1);

#pragma omp target map(a[0:sizea], a_ref[0:sizea])
    {
        for (i = 0; i < sizea; i++) {
            a[i] = 42;
            a_ref[i] = 42;
        }
    }
    
    MKL_INT bound_m = (m > 10) ? 10 : m;
    MKL_INT bound_n = (n > 10) ? 10 : n;
    
    // run gemm on host, use standard oneMKL interface   
    sger(&m, &n, &alpha, x, &incx, y, &incy, a_ref, &lda);

#pragma omp target data map(to:y[0:sizey],x[0:sizex]) map(tofrom:a[0:sizea]) device(dnum)
    {

        // run gemm on gpu, use standard oneMKL interface within a variant dispatch construct
#pragma omp target variant dispatch device(dnum) use_device_ptr(a, x, y)
        {
            sger(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
        }
        
    }
    
    float real;
    int err = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            real = a[i + lda * j] - a_ref[i + lda * j];
            real = (real > 0) ? real : -real;
            if (real > 0.0001) {
#ifdef MKL_ILP64
                printf("a[%lld][%lld] != a_ref[%lld][%lld], computed value is %f, reference value is %f, difference is %f\n",
                       i, j, i, j, a[i + lda * j], a_ref[i + lda * j], real);
#else
                printf("a[%d][%d] != a_ref[%d][%d], computed value is %f, reference value is %f, difference is %f\n",
                       i, j, i, j, a[i + lda * j], a_ref[i + lda * j], real);
#endif
                mkl_free(a);
                mkl_free(x);
                mkl_free(y);
                mkl_free(a_ref);
                return 1;
            }
        }
    }
 
    printf("Upper left corner of the A matrix:\n");
    printf("A matrix:\n");
    for (i = 0; i < bound_m; i++) {
        for (j = 0; j < bound_n; j++) {
            printf("%f\t", a[i + lda * j]);
        }
        printf("\n");
    }

    printf("Reference matrix:\n");
    for (i = 0; i < bound_m; i++) {
        for (j = 0; j < bound_n; j++) {
            printf("%f\t", a_ref[i + lda * j]);
        }
        printf("\n");
    }
       
    mkl_free(a);
    mkl_free(x);
    mkl_free(y);
    mkl_free(a_ref);
    return 0;
}
