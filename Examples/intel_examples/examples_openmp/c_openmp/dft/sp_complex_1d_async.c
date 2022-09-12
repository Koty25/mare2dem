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
! A basic example of asynchronous SINGLE-precision complex-to-complex in-place
! 1D FFT on a (GPU) device using the OpenMP target (offload) interface of
! Intel(R) oneMKL DFTI
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>

#include "mkl_dfti_omp_offload.h"

// Define the format to printf MKL_LONG values
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif

static float moda(MKL_LONG K, MKL_LONG L, MKL_LONG M);
static void init(MKL_Complex8 *x, MKL_LONG N, MKL_LONG H);
static int verify(MKL_Complex8 *x, MKL_LONG N, MKL_LONG H);

int main(void)
{
    const int devNum = 0;

    // Size of 1D transform
    const MKL_LONG N1 = 64;
    const MKL_LONG N2 = 32;

    // Arbitrary point to construct a harmonic at to be used as
    // initial data on which to perform an FFT
    const MKL_LONG H1 = -N1/2;
    const MKL_LONG H2 = -N2/2;

    MKL_LONG status     = 0;
    MKL_LONG status1GPU = 0;
    MKL_LONG status2GPU = 0;
    DFTI_DESCRIPTOR_HANDLE descHandle     = NULL;
    DFTI_DESCRIPTOR_HANDLE descHandle1GPU = NULL;
    DFTI_DESCRIPTOR_HANDLE descHandle2GPU = NULL;
    MKL_Complex8 *x     = NULL;
    MKL_Complex8 *x1GPU = NULL;
    MKL_Complex8 *x2GPU = NULL;

    printf("Create DFTI descriptor\n");
    status = DftiCreateDescriptor(&descHandle, DFTI_SINGLE, DFTI_COMPLEX, 1, N1);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Create GPU DFTI descriptor 1\n");
    status1GPU = DftiCreateDescriptor(&descHandle1GPU, DFTI_SINGLE, DFTI_COMPLEX, 1, N1);
    if (status1GPU != DFTI_NO_ERROR) goto failed;

    printf("Create GPU DFTI descriptor 2\n");
    status2GPU = DftiCreateDescriptor(&descHandle2GPU, DFTI_SINGLE, DFTI_COMPLEX, 1, N2);
    if (status2GPU != DFTI_NO_ERROR) goto failed;

    printf("Commit DFTI descriptor\n");
    status = DftiCommitDescriptor(descHandle);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Commit GPU DFTI descriptor 1\n");
#pragma omp target variant dispatch device(devNum)
    {
        status1GPU = DftiCommitDescriptor(descHandle1GPU);
    }
    if (status1GPU != DFTI_NO_ERROR) goto failed;

    printf("Commit GPU DFTI descriptor 2\n");
#pragma omp target variant dispatch device(devNum)
    {
        status2GPU = DftiCommitDescriptor(descHandle2GPU);
    }
    if (status2GPU != DFTI_NO_ERROR) goto failed;

    printf("Allocate memory for input array\n");
    x = (MKL_Complex8 *)mkl_malloc(N1*sizeof(MKL_Complex8), 64);
    if (x == NULL) goto failed;

    printf("Allocate memory for GPU input array 1\n");
    x1GPU = (MKL_Complex8 *)mkl_malloc(N1*sizeof(MKL_Complex8), 64);
    if (x1GPU == NULL) goto failed;

    printf("Allocate memory for GPU input array 2\n");
    x2GPU = (MKL_Complex8 *)mkl_malloc(N2*sizeof(MKL_Complex8), 64);
    if (x2GPU == NULL) goto failed;

    printf("Initialize input for forward FFT\n");
    init(x, N1, H1);
    init(x1GPU, N1, H1);
    init(x2GPU, N2, H2);

    printf("Compute forward FFT in-place\n");
    status = DftiComputeForward(descHandle, x);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Compute GPU forward FFT 1 in-place\n");
#pragma omp target data map(tofrom:x1GPU[0:N1], x2GPU[0:N2]) device(devNum)
    {
#pragma omp target variant dispatch use_device_ptr(x1GPU) device(devNum) nowait
        {
            status1GPU = DftiComputeForward(descHandle1GPU, x1GPU);
        }
    printf("Compute GPU forward FFT 2 in-place\n");
#pragma omp target variant dispatch use_device_ptr(x2GPU) device(devNum) nowait
        {
            status2GPU = DftiComputeForward(descHandle2GPU, x2GPU);
        }
#pragma omp taskwait
    }
    if (status1GPU != DFTI_NO_ERROR) goto failed;
    if (status2GPU != DFTI_NO_ERROR) goto failed;

    printf("Verify the result of forward FFT\n");
    status = verify(x, N1, H1);
    if (status != 0) goto failed;

    printf("Verify the result of GPU forward FFT 1\n");
    status1GPU = verify(x1GPU, N1, H1);
    if (status1GPU != 0) goto failed;

    printf("Verify the result of GPU forward FFT 2\n");
    status2GPU = verify(x2GPU, N2, H2);
    if (status2GPU != 0) goto failed;

 cleanup:
    DftiFreeDescriptor(&descHandle);
    DftiFreeDescriptor(&descHandle1GPU);
    DftiFreeDescriptor(&descHandle2GPU);

    mkl_free(x);
    mkl_free(x1GPU);
    mkl_free(x2GPU);

    {
        const MKL_LONG statusAll = status | status1GPU | status2GPU;
        printf("TEST %s\n", statusAll == 0 ? "PASSED" : "FAILED");
        return statusAll;
    }

 failed:
    printf(" ERROR,"
           " status = "    LI
           " status1GPU = " LI
           " status2GPU = " LI "\n", status, status1GPU, status2GPU);
    goto cleanup;
}

// Compute (K*L)%M accurately
static float moda(MKL_LONG K, MKL_LONG L, MKL_LONG M)
{
    return (float)(((long long)K * L) % M);
}

const float TWOPI = 6.2831853071795864769f;

// Initialize array with harmonic H
static void init(MKL_Complex8 *x, MKL_LONG N, MKL_LONG H)
{
    for (MKL_LONG n = 0; n < N; ++n) {
        float phase  = moda(n, H, N) / N;
        x[n].real = cosf(TWOPI * phase) / N;
        x[n].imag = sinf(TWOPI * phase) / N;
    }
}

// Verify that x has unit peak at H
static int verify(MKL_Complex8 *x, MKL_LONG N, MKL_LONG H)
{
    const float errthr = 5.0f * logf((float) N) / logf(2.0f) * FLT_EPSILON;
    printf(" Verifying the result, max err threshold = %.3lg\n", errthr);

    float maxerr = 0.0f;
    for (MKL_LONG n = 0; n < N; ++n) {
        float re_exp = 0.0f, im_exp = 0.0f, re_got, im_got;

        if ((n-H)%N==0) re_exp = 1.0f;

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
    printf(" Verified, maximum error was %.3lg\n", maxerr);
    return 0;
}
