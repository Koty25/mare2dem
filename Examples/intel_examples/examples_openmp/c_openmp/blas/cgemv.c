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
!      oneMKL CGEMV OpenMP offload Example Program Text 
!******************************************************************************/

#include <stdio.h>
#include <omp.h>
#include "mkl.h"
#include "mkl_omp_offload.h"
#include "common.h"

int dnum = 0;

int main() {

    MKL_Complex8 *a, *x, *y, *y_ref, alpha, beta;
    MKL_INT m, n, lda, incx, incy, i, j;
    MKL_INT sizea, sizex, sizey;

    alpha.real = 1.0; alpha.imag = 0.0;
    beta.real = 1.0; beta.imag = 0.0;

    m = 1449;
    n = 1120;

    lda = m;
    incx = 1;
    incy = 1;

    sizea = lda * n;
    sizex = n;
    sizey = m;
    
    // allocate matrices
    a = (MKL_Complex8 *)mkl_malloc(sizea * sizeof(MKL_Complex8), 64);
    x = (MKL_Complex8 *)mkl_malloc(sizex * sizeof(MKL_Complex8), 64);
    y = (MKL_Complex8 *)mkl_malloc(sizey * sizeof(MKL_Complex8), 64);
    y_ref = (MKL_Complex8 *)mkl_malloc(sizey * sizeof(MKL_Complex8), 64);

    if ((a == NULL) || (x == NULL) || (y == NULL) || (y_ref == NULL)) {
        printf("Cannot allocate matrices\n");
        return 1;
    }

    // initialize matrices
    init_single_complex_array(sizea, a, 1);
    init_single_complex_array(sizex, x, 1);

#pragma omp target map(y[0:sizey], y_ref[0:sizey])
    {
        for (i = 0; i < sizey; i++) {
            y[i].real = 42;
            y_ref[i].real = 42;
            y[i].imag = 42;
            y_ref[i].imag = 42;
        }
    }
    
    MKL_INT bound_m = (m > 10) ? 10 : m;
    
    // run gemm on host, use standard oneMKL interface   
    cgemv("N", &m, &n, &alpha, a, &lda, x, &incx, &beta, y_ref, &incy);

#pragma omp target data map(to:a[0:sizea],x[0:sizex]) map(tofrom:y[0:sizey]) device(dnum)
    {

        // run gemm on gpu, use standard oneMKL interface within a variant dispatch construct
#pragma omp target variant dispatch device(dnum) use_device_ptr(a, x, y)
        {
            cgemv("N", &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
        }
        
    }
    
    float real, imag;
    int err = 0;
    for (i = 0; i < m; i++) {
        real = y[i].real - y_ref[i].real;
        real = (real > 0) ? real : -real;
        imag = y[i].imag - y_ref[i].imag;
        imag = (imag > 0) ? imag : -imag;
            
        if ((real + imag) > 0.001) {
#ifdef MKL_ILP64
            printf("c[%lld][%lld] != c_ref[%lld][%lld], computed value is %lf + i%lf, reference value is %lf + i%lf, difference is %lf + i%lf\n",
                   i, j, i, j, y[i].real, y[i].imag, y_ref[i].real, y_ref[i].imag, real, imag);
#else
            printf("c[%d][%d] != c_ref[%d][%d], computed value is %lf + i%lf, reference value is %lf + i%lf, difference is %lf + i%lf\n",
                   i, j, i, j, y[i].real, y[i].imag, y_ref[i].real, y_ref[i].imag, real, imag);
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
        printf("%lf + i%lf\t", y[i].real, y[i].imag);
    }
    printf("\n");
    
    printf("Reference vector:\n");
    for (i = 0; i < bound_m; i++) {
        printf("%lf + i%lf\t", y_ref[i].real, y_ref[i].imag);
    }
    printf("\n");
       
    mkl_free(a);
    mkl_free(x);
    mkl_free(y);
    mkl_free(y_ref);
    return 0;
}
