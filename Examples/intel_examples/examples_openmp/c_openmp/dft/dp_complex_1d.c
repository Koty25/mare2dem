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
! A basic example of DOUBLE-precision complex-to-complex in-place 1D FFT on a
! (GPU) device using the OpenMP target (offload) interface of Intel(R) oneMKL
! DFTI
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

static double moda(MKL_LONG K, MKL_LONG L, MKL_LONG M);
static void init(MKL_Complex16 *x, MKL_LONG N, MKL_LONG H);
static int verify(MKL_Complex16 *x, MKL_LONG N, MKL_LONG H);

int main(void)
{
    const int devNum = 0;

    // Size of 1D transform
    const MKL_LONG N = 64;

    // Arbitrary point to construct a harmonic at to be used as
    // initial data on which to perform an FFT
    const MKL_LONG H = -N/2;

    MKL_LONG status    = 0;
    MKL_LONG statusGPU = 0;
    DFTI_DESCRIPTOR_HANDLE descHandle    = NULL;
    DFTI_DESCRIPTOR_HANDLE descHandleGPU = NULL;
    MKL_Complex16 *x    = NULL;
    MKL_Complex16 *xGPU = NULL;

    printf("Create DFTI descriptor\n");
    status = DftiCreateDescriptor(&descHandle, DFTI_DOUBLE, DFTI_COMPLEX, 1, N);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Create GPU DFTI descriptor\n");
    statusGPU = DftiCreateDescriptor(&descHandleGPU, DFTI_DOUBLE, DFTI_COMPLEX, 1, N);
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Commit DFTI descriptor\n");
    status = DftiCommitDescriptor(descHandle);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Commit GPU DFTI descriptor\n");
#pragma omp target variant dispatch device(devNum)
    {
        statusGPU = DftiCommitDescriptor(descHandleGPU);
    }
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Allocate memory for input array\n");
    x = (MKL_Complex16 *)mkl_malloc(N*sizeof(MKL_Complex16), 64);
    if (x == NULL) goto failed;

    printf("Allocate memory for GPU input array\n");
    xGPU = (MKL_Complex16 *)mkl_malloc(N*sizeof(MKL_Complex16), 64);
    if (xGPU == NULL) goto failed;

    printf("Initialize input for forward FFT\n");
    init(x, N, H);
    init(xGPU, N, H);

    printf("Compute forward FFT in-place\n");
    status = DftiComputeForward(descHandle, x);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Compute GPU forward FFT in-place\n");
#pragma omp target data map(tofrom:xGPU[0:N]) device(devNum)
    {
#pragma omp target variant dispatch use_device_ptr(xGPU) device(devNum)
        {
            statusGPU = DftiComputeForward(descHandleGPU, xGPU);
        }
    }
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Verify the result of forward FFT\n");
    status = verify(x, N, H);
    if (status != 0) goto failed;

    printf("Verify the result of GPU forward FFT\n");
    statusGPU = verify(xGPU, N, H);
    if (statusGPU != 0) goto failed;

 cleanup:
    DftiFreeDescriptor(&descHandle);
    DftiFreeDescriptor(&descHandleGPU);

    mkl_free(x);
    mkl_free(xGPU);

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
static double moda(MKL_LONG K, MKL_LONG L, MKL_LONG M)
{
    return (double)(((long long)K * L) % M);
}

const double TWOPI = 6.2831853071795864769;

// Initialize array with harmonic H
static void init(MKL_Complex16 *x, MKL_LONG N, MKL_LONG H)
{
    for (MKL_LONG n = 0; n < N; ++n) {
        double phase  = moda(n, H, N) / N;
        x[n].real = cos(TWOPI * phase) / N;
        x[n].imag = sin(TWOPI * phase) / N;
    }
}

// Verify that x has unit peak at H
static int verify(MKL_Complex16 *x, MKL_LONG N, MKL_LONG H)
{
    const double errthr = 5.0 * log((double) N) / log(2.0) * DBL_EPSILON;
    printf(" Verifying the result, max err threshold = %.3lg\n", errthr);

    double maxerr = 0.0;
    for (MKL_LONG n = 0; n < N; ++n) {
        double re_exp = 0.0, im_exp = 0.0, re_got, im_got;

        if ((n-H)%N==0) re_exp = 1.0;

        re_got = x[n].real;
        im_got = x[n].imag;
        double err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr)) {
            printf(" x[" LI "]: ", n);
            printf(" expected (%.17lg,%.17lg), ", re_exp, im_exp);
            printf(" got (%.17lg,%.17lg), ", re_got, im_got);
            printf(" err %.3lg\n", err);
            printf(" Verification FAILED\n");
            return 1;
        }
    }
    printf(" Verified, maximum error was %.3lg\n", maxerr);
    return 0;
}
