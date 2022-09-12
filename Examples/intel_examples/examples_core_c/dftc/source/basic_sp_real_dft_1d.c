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
! A simple example of single-precision real-to-complex out-of-place 1D
! FFT using Intel(R) Math Kernel Library (Intel(R) MKL) DFTI
!
!****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mkl_service.h"
#include "mkl_dfti.h"

static void init_r(float *x, int N, int H);
static int verify_c(MKL_Complex8 *x, int N, int H);
static void init_c(MKL_Complex8 *x, int N, int H);
static int verify_r(float *x, int N, int H);

/* Define the format to printf MKL_LONG values */
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif

int main(void)
{
    /* Size of 1D transform */
    const int N = 6;

    /* Arbitrary harmonic used to verify FFT */
    int H = -1;

    /* Execution status */
    MKL_LONG status = 0;

    /* Pointers to input and output data */
    float *x_real = NULL;
    MKL_Complex8 *x_cmplx = NULL;

    DFTI_DESCRIPTOR_HANDLE hand = NULL;

    char version[DFTI_VERSION_LENGTH];

    DftiGetValue(0, DFTI_VERSION, version);
    printf("%s\n", version);

    printf("Example basic_sp_real_dft_1d\n");
    printf("Forward-Backward single-precision real-to-complex"
           " out-of-place 1D transform\n");
    printf("Configuration parameters:\n");
    printf(" DFTI_PRECISION                = DFTI_SINGLE\n");
    printf(" DFTI_FORWARD_DOMAIN           = DFTI_REAL\n");
    printf(" DFTI_DIMENSION                = 1\n");
    printf(" DFTI_LENGTHS                  = {%d}\n", N);
    printf(" DFTI_PLACEMENT                = DFTI_NOT_INPLACE\n");
    printf(" DFTI_CONJUGATE_EVEN_STORAGE   = DFTI_COMPLEX_COMPLEX\n");

    printf("Create DFTI descriptor\n");
    status = DftiCreateDescriptor(&hand, DFTI_SINGLE, DFTI_REAL,
                                  1, (MKL_LONG)N);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set configuration: out-of-place\n");
    status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set configuration: CCE storage\n");
    status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE,
                          DFTI_COMPLEX_COMPLEX);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Commit the  descriptor\n");
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Allocate data arrays\n");
    x_real  = (float *)mkl_malloc(N*sizeof(float), 64);
    x_cmplx = (MKL_Complex8*)mkl_malloc((N/2+1)*sizeof(MKL_Complex8), 64);
    if (x_real == NULL || x_cmplx == NULL) goto failed;

    printf("Initialize data for real-to-complex FFT\n");
    init_r(x_real, N, H);

    printf("Compute forward transform\n");
    status = DftiComputeForward(hand, x_real, x_cmplx);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Verify the result\n");
    status = verify_c(x_cmplx, N, H);
    if (status != 0) goto failed;

    printf("Initialize data for complex-to-real FFT\n");
    init_c(x_cmplx, N, H);

    printf("Compute backward transform\n");
    status = DftiComputeBackward(hand, x_cmplx, x_real);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Verify the result\n");
    status = verify_r(x_real, N, H);
    if (status != 0) goto failed;

 cleanup:

    printf("Free DFTI descriptor\n");
    DftiFreeDescriptor(&hand);

    printf("Free data arrays\n");
    mkl_free(x_real);
    mkl_free(x_cmplx);

    printf("TEST %s\n", (status == 0) ? "PASSED" : "FAILED");
    return status;

 failed:
    printf(" ERROR, status = "LI"\n", status);
    status = 1;
    goto cleanup;
}

/* Compute (K*L)%M accurately */
static float moda(int K, int L, int M)
{
    return (float)(((long long)K * L) % M);
}

/* Initialize array x(N) to produce unit peaks at y(H) and y(N-H) */
static void init_r(float *x, int N, int H)
{
    float TWOPI = 6.2831853071795864769f, phase, factor;
    int n;

    factor = (2*(N-H)%N == 0) ? 1.0f : 2.0f;
    for (n = 0; n < N; n++) {
        phase  = moda(n, H, N) / N;
        x[n] = factor * cosf(TWOPI * phase) / N;
    }
}

/* Verify that x has unit peak at H */
static int verify_c(MKL_Complex8 *x, int N, int H)
{
    float err, errthr, maxerr;
    int n;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5f * logf((float) N) / logf(2.0f) * FLT_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    maxerr = 0.0f;
    for (n = 0; n < N/2+1; n++) {
        float re_exp = 0.0f, im_exp = 0.0f, re_got, im_got;

        if ((n-H)%N == 0 || (-n-H)%N == 0) re_exp = 1.0f;

        re_got = x[n].real;
        im_got = x[n].imag;
        err  = fabsf(re_got - re_exp) + fabsf(im_got - im_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr)) {
            printf(" x[%i]: ", n);
            printf(" expected (%.7g,%.7g), ", re_exp, im_exp);
            printf(" got (%.7g,%.7g), ", re_got, im_got);
            printf(" err %.3lg\n", err);
            printf(" Verification FAILED\n");
            return 1;
        }
    }
    printf(" Verified,  maximum error was %.3lg\n", maxerr);
    return 0;
}

/* Initialize array x(N) to produce unit peak at y(H) */
static void init_c(MKL_Complex8 *x, int N, int H)
{
    float TWOPI = 6.2831853071795864769f, phase;
    int n;

    for (n = 0; n < N/2+1; n++) {
        phase  = moda(n, H, N) / N;
        x[n].real =  cosf(TWOPI * phase) / N;
        x[n].imag = -sinf(TWOPI * phase) / N;
    }
}

/* Verify that x has unit peak at H */
static int verify_r(float *x, int N, int H)
{
    float err, errthr, maxerr;
    int n;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5f * logf((float) N) / logf(2.0f) * FLT_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    maxerr = 0.0f;
    for (n = 0; n < N; n++) {
        float re_exp = 0.0f, re_got;

        if ((n-H)%N == 0) re_exp = 1.0f;

        re_got = x[n];
        err  = fabsf(re_got - re_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr)) {
            printf(" x[%i]: ", n);
            printf(" expected %.7g, ", re_exp);
            printf(" got %.7g, ", re_got);
            printf(" err %.3lg\n", err);
            printf(" Verification FAILED\n");
            return 1;
        }
    }
    printf(" Verified,  maximum error was %.3lg\n", maxerr);
    return 0;
}
