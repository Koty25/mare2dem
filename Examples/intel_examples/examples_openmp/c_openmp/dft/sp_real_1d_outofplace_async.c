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
! An example of asynchronous SINGLE-precision real-to-complex out-of-place 1D
! FFT on a (GPU) device using the OpenMP target (offload) interface of Intel(R)
! oneMKL DFTI
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>

#include "mkl_dfti_omp_offload.h"

static void init_r(float *x, MKL_LONG N, MKL_LONG H);
static int verify_c(MKL_Complex8 *x, MKL_LONG N, MKL_LONG H);
static void init_c(MKL_Complex8 *x, MKL_LONG N, MKL_LONG H);
static int verify_r(float *x, MKL_LONG N, MKL_LONG H);

// Define the format to printf MKL_LONG values
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif

int main(void)
{
    const int devNum = 0;

    // Size of 1D FFT
    const MKL_LONG N = 8;

    const MKL_LONG halfNplus1 = N/2 + 1;

    // Arbitrary harmonic used to verify FFT
    MKL_LONG H = -1;

    MKL_LONG status    = 0;
    MKL_LONG statusGPU = 0;

    // Pointers to input and output data
    float        *x_real     = NULL;
    float        *xGPU_real  = NULL;
    MKL_Complex8 *x_cmplx    = NULL;
    MKL_Complex8 *xGPU_cmplx = NULL;

    DFTI_DESCRIPTOR_HANDLE descHandle    = NULL;
    DFTI_DESCRIPTOR_HANDLE descHandleGPU = NULL;

    printf("DFTI_LENGTHS                  = {" LI "}\n", N);
    printf("DFTI_PLACEMENT                = DFTI_NOT_INPLACE\n");
    printf("DFTI_CONJUGATE_EVEN_STORAGE   = DFTI_COMPLEX_COMPLEX\n");

    printf("Create DFTI descriptor\n");
    status = DftiCreateDescriptor(&descHandle, DFTI_SINGLE, DFTI_REAL, 1, N);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Create GPU DFTI descriptor\n");
    statusGPU = DftiCreateDescriptor(&descHandleGPU, DFTI_SINGLE, DFTI_REAL, 1, N);
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Set configuration: out-of-place\n");
    status = DftiSetValue(descHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set GPU configuration: out-of-place\n");
    statusGPU = DftiSetValue(descHandleGPU, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Set configuration: CCE storage\n");
    status = DftiSetValue(descHandle,
                          DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set GPU configuration: CCE storage\n");
    statusGPU = DftiSetValue(descHandleGPU,
                             DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Commit descriptor\n");
    status = DftiCommitDescriptor(descHandle);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Commit GPU descriptor\n");
#pragma omp target variant dispatch device(devNum)
    {
        statusGPU = DftiCommitDescriptor(descHandleGPU);
    }
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Allocate data arrays\n");
    x_real  = (float*)mkl_malloc(N*sizeof(float), 64);
    x_cmplx = (MKL_Complex8*)mkl_malloc((halfNplus1)*sizeof(MKL_Complex8), 64);
    if (x_real == NULL || x_cmplx == NULL) goto failed;

    printf("Allocate GPU data arrays\n");
    xGPU_real  = (float*)mkl_malloc(N*sizeof(float), 64);
    xGPU_cmplx = (MKL_Complex8*)mkl_malloc((halfNplus1)*sizeof(MKL_Complex8), 64);
    if (xGPU_real == NULL || xGPU_cmplx == NULL) goto failed;

    printf("Initialize data for real-to-complex FFT\n");
    init_r(x_real, N, H);

    printf("Initialize GPU data for real-to-complex FFT\n");
    init_r(xGPU_real, N, H);

    printf("Compute forward FFT\n");
    status = DftiComputeForward(descHandle, x_real, x_cmplx);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Compute GPU forward FFT\n");
#pragma omp target data map(to:xGPU_real[0:N]) \
                        map(from:xGPU_cmplx[0:halfNplus1]) device(devNum)
    {
#pragma omp target variant dispatch use_device_ptr(xGPU_real, xGPU_cmplx) \
                                    device(devNum) nowait
        {
            statusGPU = DftiComputeForward(descHandleGPU, xGPU_real, xGPU_cmplx);
        }
#pragma omp taskwait
    }
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Verify the complex result\n");
    status = verify_c(x_cmplx, N, H);
    if (status != 0) goto failed;

    printf("Verify the GPU complex result\n");
    statusGPU = verify_c(xGPU_cmplx, N, H);
    if (statusGPU != 0) goto failed;

    printf("Initialize data for complex-to-real FFT\n");
    init_c(x_cmplx, N, H);

    printf("Initialize GPU data for complex-to-real FFT\n");
    init_c(xGPU_cmplx, N, H);

    printf("Compute backward FFT\n");
    status = DftiComputeBackward(descHandle, x_cmplx, x_real);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Compute GPU backward FFT\n");
#pragma omp target data map(to:xGPU_cmplx[0:halfNplus1]) \
                        map(from:xGPU_real[0:N]) device(devNum)
    {
#pragma omp target variant dispatch use_device_ptr(xGPU_real, xGPU_cmplx) \
                                    device(devNum) nowait
        {
            statusGPU = DftiComputeBackward(descHandleGPU, xGPU_cmplx, xGPU_real);
        }
#pragma omp taskwait
    }
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Verify the result\n");
    status = verify_r(x_real, N, H);
    if (status != 0) goto failed;

    printf("Verify the GPU result\n");
    status = verify_r(xGPU_real, N, H);
    if (status != 0) goto failed;

 cleanup:

    printf("Free DFTI descriptor\n");
    DftiFreeDescriptor(&descHandle);

    printf("Free GPU DFTI descriptor\n");
    DftiFreeDescriptor(&descHandleGPU);

    printf("Free data arrays\n");
    mkl_free(x_real);
    mkl_free(x_cmplx);

    printf("Free GPU data arrays\n");
    mkl_free(xGPU_real);
    mkl_free(xGPU_cmplx);

    {
        const MKL_LONG statusBoth = status | statusGPU;
        printf("TEST %s\n", statusBoth == 0 ? "PASSED" : "FAILED");
        return statusBoth;
    }

 failed:
    printf(" ERROR,"
           " status = "    LI
           " statusGPU = " LI "\n", status, statusGPU);
    goto cleanup;
}

// Compute (K*L)%M accurately
static float moda(MKL_LONG K, MKL_LONG L, MKL_LONG M)
{
    return (float)(((long long)K * L) % M);
}

const float TWOPI = 6.2831853071795864769f;

// Initialize array x to produce unit peaks at y(H) and y(N-H)
static void init_r(float* x, MKL_LONG N, MKL_LONG H)
{
    const float factor = (2 * (N - H) % N == 0) ? 1.0f : 2.0f;
    for (MKL_LONG n = 0; n < N; ++n) {
        float phase = moda(n, H, N) / N;
        x[n] = factor * cosf(TWOPI * phase) / N;
    }
}

// Verify that x has unit peak at H
static int verify_c(MKL_Complex8 *x, MKL_LONG N, MKL_LONG H)
{
    const float errthr = 2.5f * logf((float) N) / logf(2.0f) * FLT_EPSILON;
    printf(" Verifying the result, max err threshold = %.3lg\n", errthr);

    float maxerr = 0.0f;
    for (MKL_LONG n = 0; n < N/2+1; ++n) {
        float re_exp = 0.0f, im_exp = 0.0f, re_got, im_got;

        if ((n-H)%N == 0 || (-n-H)%N == 0) re_exp = 1.0f;

        re_got = x[n].real;
        im_got = x[n].imag;
        float err  = fabsf(re_got - re_exp) + fabsf(im_got - im_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr)) {
            printf(" x[" LI "]: ", n);
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

// Initialize array x to produce unit peak at y(H)
static void init_c(MKL_Complex8 *x, MKL_LONG N, MKL_LONG H)
{
    for (MKL_LONG n = 0; n < N/2+1; n++) {
        float phase  = moda(n, H, N) / N;
        x[n].real =  cosf(TWOPI * phase) / N;
        x[n].imag = -sinf(TWOPI * phase) / N;
    }
}

// Verify that x has unit peak at H
static int verify_r(float *x, MKL_LONG N, MKL_LONG H)
{
    const float errthr = 2.5f * logf((float) N) / logf(2.0f) * FLT_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    float maxerr = 0.0f;
    for (MKL_LONG n = 0; n < N; n++) {
        float re_exp = 0.0f, re_got;

        if ((n-H)%N == 0) re_exp = 1.0f;

        re_got = x[n];
        float err  = fabsf(re_got - re_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr)) {
            printf(" x[" LI "]: ", n);
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
