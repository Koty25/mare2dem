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
!       Example of using fftwf_plan_dft_1d function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init(fftwf_complex *x, int N, int H);
static int verify(fftwf_complex *x, int N, int H);

int main(void)
{
    /* Size of 1D transform */
    int N = 64;

    /* Arbitrary harmonic used to verify FFT */
    int H = -N/2;

    /* FFTW plan handles */
    fftwf_plan forward_plan = 0, backward_plan = 0;

    /* Pointer to input/output data */
    fftwf_complex *x = 0;

    /* Execution status */
    int status = 0;

    printf("Example sp_plan_dft_1d\n");
    printf("Forward and backward 1D complex inplace transform\n");
    printf("Configuration parameters:\n");
    printf(" N = %d\n", N);
    printf(" H = %d\n", H);

    printf("Allocate array for input data\n");
    x  = fftwf_malloc(sizeof(fftwf_complex)*N);
    if (0 == x) goto failed;

    printf("Create FFTW plan for 1D double-precision forward transform\n");
    forward_plan = fftwf_plan_dft(1, &N, x, x, FFTW_FORWARD, FFTW_ESTIMATE);
    if (0 == forward_plan) goto failed;

    printf("Create FFTW plan for 1D double-precision backward transform\n");
    backward_plan = fftwf_plan_dft(1, &N, x, x, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (0 == backward_plan) goto failed;

    printf("Initialize input for forward transform\n");
    init(x, N, H);

    printf("Compute forward FFT\n");
    fftwf_execute(forward_plan);

    printf("Verify the result of forward FFT\n");
    status = verify(x, N, H);
    if (0 != status) goto failed;

    printf("Initialize input for backward transform\n");
    init(x, N, -H);

    printf("Compute backward transform using new-array function\n");
    fftwf_execute_dft(backward_plan, x, x);

    printf("Verify the result of backward FFT\n");
    status = verify(x, N, H);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plans\n");
    fftwf_destroy_plan(forward_plan);
    fftwf_destroy_plan(backward_plan);

    printf("Free data array\n");
    fftwf_free(x);

    printf("TEST %s\n",0==status ? "PASSED" : "FAILED");
    return status;

 failed:
    printf(" ERROR\n");
    status = 1;
    goto cleanup;
}

/* Compute (K*L)%M accurately */
static float moda(int K, int L, int M)
{
    return (float)(((long long)K * L) % M);
}

/* Initialize array x[N] with harmonic H */
static void init(fftwf_complex *x, int N, int H)
{
    float TWOPI = 6.2831853071795864769f, phase;
    int n;

    for (n = 0; n < N; n++)
    {
        phase  = moda(n,H,N) / N;
        x[n][0] = cosf( TWOPI * phase ) / N;
        x[n][1] = sinf( TWOPI * phase ) / N;
    }
}

/* Verify that x has unit peak at H */
static int verify(fftwf_complex *x, int N, int H)
{
    float err, errthr, maxerr;
    int n;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 5.0f * logf( (float)N ) / logf(2.0f) * FLT_EPSILON;
    printf(" Verify the result, errthr = %.3g\n", errthr);

    maxerr = 0;
    for (n = 0; n < N; n++)
    {
        float re_exp = 0.0, im_exp = 0.0, re_got, im_got;

        if ((n-H)%N==0)
        {
            re_exp = 1;
        }

        re_got = x[n][0];
        im_got = x[n][1];
        err  = fabsf(re_got - re_exp) + fabsf(im_got - im_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr))
        {
            printf(" x[%i]: ",n);
            printf(" expected (%.7g,%.7g), ",re_exp,im_exp);
            printf(" got (%.7g,%.7g), ",re_got,im_got);
            printf(" err %.3g\n", err);
            printf(" Verification FAILED\n");
            return 1;
        }
    }
    printf(" Verified, maximum error was %.3g\n", maxerr);
    return 0;
}
