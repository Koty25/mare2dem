/*******************************************************************************
* Copyright 2011-2020 Intel Corporation.
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
! Content:
!       Example of using fftw_plan_dft_r2c_1d function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init_r(double *x, int N, int H);
static int verify_c(fftw_complex *x, int N, int H);

int main(void)
{
    /* Size of 1D transform */
    int N = 1024;

    /* Arbitrary harmonic used to verify FFT */
    int H = -1;

    /* FFTW plan handle */
    fftw_plan r2c = 0;

    /* Pointers to input and output data */
    double *x = 0;
    fftw_complex *y = 0;

    /* Execution status */
    int status = 0;

    printf("Example dp_plan_dft_r2c_1d\n");
    printf("1D real-to-complex out-of-place transform\n");
    printf("Configuration parameters:\n");
    printf(" N = %d\n", N);
    printf(" H = %d\n", H);

    printf("Allocate data arrays: real x(%i), complex y(%i)\n", N, N/2+1);
    x  = fftw_malloc(sizeof(double)*N);
    y  = fftw_malloc(sizeof(fftw_complex)*(N/2+1));
    if (0 == x || y == 0) goto failed;

    printf("Create FFTW plan for 1D r2c out-of-place FFT\n");
    r2c = fftw_plan_dft_r2c_1d(N, x, y, FFTW_ESTIMATE);
    if (0 == r2c) goto failed;

    printf("Initialize input for r2c transform\n");
    init_r(x, N, H);

    printf("Compute r2c FFT\n");
    fftw_execute(r2c);

    printf("Verify the result of forward FFT\n");
    status = verify_c(y, N, H);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plan\n");
    fftw_destroy_plan(r2c);

    printf("Free data arrays\n");
    fftw_free(x);
    fftw_free(y);

    printf("TEST %s\n",0==status ? "PASSED" : "FAILED");
    return status;

 failed:
    printf(" ERROR\n");
    status = 1;
    goto cleanup;
}

/* Compute (K*L)%M accurately */
static double moda(int K, int L, int M)
{
    return (double)(((long long)K * L) % M);
}

/* Initialize array x(N) to produce unit peaks at y(H) and y(N-H) */
static void init_r(double *x, int N, int H)
{
    double TWOPI = 6.2831853071795864769, phase, factor;
    int n;

    factor = ((N-H)%N == 0) ? 1.0 : 2.0;
    for (n = 0; n < N; n++)
    {
        phase  = moda(n,H,N) / N;
        x[n] = factor * cos( TWOPI * phase ) / N;
    }
}

/* Verify that x has unit peak at H */
static int verify_c(fftw_complex *x, int N, int H)
{
    double err, errthr, maxerr;
    int n;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5 * log( (double)N ) / log(2.0) * DBL_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    maxerr = 0;
    for (n = 0; n < N/2+1; n++)
    {
        double re_exp = 0.0, im_exp = 0.0, re_got, im_got;

        if ((n-H)%N == 0 || (-n-H)%N == 0)
        {
            re_exp = 1;
        }

        re_got = x[n][0];
        im_got = x[n][1];
        err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr))
        {
            printf(" x[%i]: ",n);
            printf(" expected (%.17lg,%.17lg), ",re_exp,im_exp);
            printf(" got (%.17lg,%.17lg), ",re_got,im_got);
            printf(" err %.3lg\n", err);
            printf(" Verification FAILED\n");
            return 1;
        }
    }
    printf(" Verified,  maximum error was %.3lg\n", maxerr);
    return 0;
}
