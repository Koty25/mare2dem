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
!      oneMKL CBLAS_DGEMV OpenMP offload Example Program Text 
!******************************************************************************/

#include <stdio.h>
#include <omp.h>
#include "mkl.h"
#include "mkl_omp_offload.h"
#include "common.h"

int dnum = 0;

int main() {

    double *a, *x, *y, *y_ref, alpha, beta;
    MKL_INT m, n, lda, incx, incy, i, j;
    MKL_INT sizea, sizex, sizey;
    
    alpha = 1.0;
    beta = 1.0;

    m = 1449;
    n = 1120;

    lda = m;
    incx = 1;
    incy = 1;

    sizea = lda * n;
    sizex = n;
    sizey = m;
    
    // allocate matrices
    a = (double *)mkl_malloc(sizea * sizeof(double), 64);
    x = (double *)mkl_malloc(sizex * sizeof(double), 64);
    y = (double *)mkl_malloc(sizey * sizeof(double), 64);
    y_ref = (double *)mkl_malloc(sizey * sizeof(double), 64);

    if ((a == NULL) || (x == NULL) || (y == NULL) || (y_ref == NULL)) {
        printf("Cannot allocate matrices\n");
        return 1;
    }

    // initialize matrices
    init_double_array(sizea, a, 1);
    init_double_array(sizex, x, 1);

#pragma omp target map(y[0:sizey], y_ref[0:sizey])
    {
        for (i = 0; i < sizey; i++) {
            y[i] = 42;
            y_ref[i] = 42;
        }
    }
    
    MKL_INT bound_m = (m > 10) ? 10 : m;
    
    // run gemm on host, use standard oneMKL interface   
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, alpha, a, lda, x, incx, beta, y_ref, incy);

#pragma omp target data map(to:a[0:sizea],x[0:sizex]) map(tofrom:y[0:sizey]) device(dnum)
    {

        // run gemm on gpu, use standard oneMKL interface within a variant dispatch construct
#pragma omp target variant dispatch device(dnum) use_device_ptr(a, x, y)
        {
            cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, alpha, a, lda, x, incx, beta, y, incy);
        }
        
    }
    
    double real;
    int err = 0;
    for (i = 0; i < m; i++) {
        real = y[i] - y_ref[i];
        real = (real > 0) ? real : -real;
        if (real > 0.0001) {
#ifdef MKL_ILP64
            printf("y[%lld] != y_ref[%lld], computed value is %lf, reference value is %lf, difference is %lf\n",
                   i, i, y[i], y_ref[i], real);
#else
            printf("y[%d] != y_ref[%d], computed value is %lf, reference value is %lf, difference is %lf\n",
                   i, i, y[i], y_ref[i], real);
#endif
            mkl_free(a);
            mkl_free(x);
            mkl_free(y);
            mkl_free(y_ref);
            return 1;
        }
    }
 
    printf("First elements of the output vector Y:\n");
    printf("Y vector:\n");
    for (i = 0; i < bound_m; i++) {
        printf("%lf\t", y[i]);
    }
    printf("\n");
    
    printf("Reference vector:\n");
    for (i = 0; i < bound_m; i++) {
        printf("%lf\t", y_ref[i]);
    }
    printf("\n");
       
    mkl_free(a);
    mkl_free(x);
    mkl_free(y);
    mkl_free(y_ref);
    return 0;
}
