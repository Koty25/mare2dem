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
! A basic example of DOUBLE-precision complex-to-complex in-place 2D FFT on a
! (GPU) device using the OpenMP target (offload) interface of Intel(R) oneMKL
! DFTI
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>

#include "mkl_dfti_omp_offload.h"

static void init(MKL_Complex16 *data, MKL_LONG N1, MKL_LONG N2, MKL_LONG H1, MKL_LONG H2);
static int verify(MKL_Complex16 *data, MKL_LONG N1, MKL_LONG N2, MKL_LONG H1, MKL_LONG H2);

// Define the format to printf MKL_LONG values
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif

int main(void)
{
    const int devNum = 0;

    // Sizes of 2D FFT
    MKL_LONG N1 = 11, N2 = 8;
    MKL_LONG N[2] = {N2, N1};   // N2 rows x N1 cols

    // Arbitrary harmonic used to verify FFT
    MKL_LONG H1 = -1, H2 = -2;

    MKL_LONG status    = 0;
    MKL_LONG statusGPU = 0;
    DFTI_DESCRIPTOR_HANDLE descHandle    = NULL;
    DFTI_DESCRIPTOR_HANDLE descHandleGPU = NULL;
    MKL_Complex16 *data    = NULL;
    MKL_Complex16 *dataGPU = NULL;

    printf("DFTI_LENGTHS = {" LI " (rows), " LI " (cols)} \n", N2, N1);

    printf("Create DFTI descriptor\n");
    status = DftiCreateDescriptor(&descHandle, DFTI_DOUBLE, DFTI_COMPLEX, 2, N);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Create GPU DFTI descriptor\n");
    statusGPU = DftiCreateDescriptor(&descHandleGPU, DFTI_DOUBLE, DFTI_COMPLEX, 2, N);
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

    printf("Allocate input array\n");
    data = (MKL_Complex16*)mkl_malloc(N2*N1*sizeof(MKL_Complex16), 64);
    if (data == NULL) goto failed;

    printf("Allocate GPU input array\n");
    dataGPU = (MKL_Complex16*)mkl_malloc(N2*N1*sizeof(MKL_Complex16), 64);
    if (dataGPU == NULL) goto failed;

    printf("Initialize input for forward FFT\n");
    init(data, N1, N2, H1, H2);

    printf("Initialize GPU input for forward FFT\n");
    init(dataGPU, N1, N2, H1, H2);

    printf("Compute forward FFT\n");
    status = DftiComputeForward(descHandle, data);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Compute GPU forward FFT in-place\n");
#pragma omp target data map(tofrom:dataGPU[0:N2*N1]) device(devNum)
    {
#pragma omp target variant dispatch use_device_ptr(dataGPU) device(devNum)
        {
            statusGPU = DftiComputeForward(descHandleGPU, dataGPU);
        }
    }
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Verify the result of forward FFT\n");
    status = verify(data, N1, N2, H1, H2);
    if (status != 0) goto failed;

    printf("Verify the result of GPU forward FFT\n");
    statusGPU = verify(dataGPU, N1, N2, H1, H2);
    if (statusGPU != 0) goto failed;

    printf("Initialize input for backward FFT\n");
    init(data, N1, N2, -H1, -H2);

    printf("Initialize GPU input for backward FFT\n");
    init(dataGPU, N1, N2, -H1, -H2);

    printf("Compute backward FFT\n");
    status = DftiComputeBackward(descHandle, data);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Compute GPU backward FFT in-place\n");
#pragma omp target data map(tofrom:dataGPU[0:N2*N1]) device(devNum)
    {
#pragma omp target variant dispatch use_device_ptr(dataGPU) device(devNum)
        {
            statusGPU = DftiComputeBackward(descHandleGPU, dataGPU);
        }
    }
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Verify the result of backward FFT\n");
    status = verify(data, N1, N2, H1, H2);
    if (status != 0) goto failed;

    printf("Verify the result of GPU backward FFT\n");
    statusGPU = verify(dataGPU, N1, N2, H1, H2);
    if (statusGPU != 0) goto failed;

 cleanup:
    DftiFreeDescriptor(&descHandle);
    DftiFreeDescriptor(&descHandleGPU);

    mkl_free(data);
    mkl_free(dataGPU);

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

double TWOPI = 6.2831853071795864769;

// Initialize array with harmonic (H1, H2)
static void init(MKL_Complex16 *data, MKL_LONG N1, MKL_LONG N2, MKL_LONG H1, MKL_LONG H2)
{
    // Generalized strides for row-major addressing of data
    const MKL_LONG S1 = 1, S2 = N1;

    for (MKL_LONG n2 = 0; n2 < N2; n2++) {
        for (MKL_LONG n1 = 0; n1 < N1; n1++) {
            double phase = TWOPI * (  moda(n1, H1, N1) / N1
                                   + moda(n2, H2, N2) / N2  );
            MKL_LONG index = n2*S2 + n1*S1;
            data[index].real = cos(phase) / (N2*N1);
            data[index].imag = sin(phase) / (N2*N1);
        }
    }
}

// Verify that data has unit peak at (H1, H2)
static int verify(MKL_Complex16 *data, MKL_LONG N1, MKL_LONG N2, MKL_LONG H1, MKL_LONG H2)
{
    // Generalized strides for row-major addressing of data
    const MKL_LONG S1 = 1, S2 = N1;

    const double errthr = 5.0 * log((double) N2*N1) / log(2.0) * DBL_EPSILON;
    printf(" Verify the result, errthr = %.3lg\n", errthr);

    int isWrong = 0;
    double maxerr = 0.0;
    for (MKL_LONG n2 = 0; n2 < N2; ++n2) {
        for (MKL_LONG n1 = 0; n1 < N1; ++n1) {
            double re_exp = 0.0, im_exp = 0.0, re_got, im_got;

            if ((n1-H1)%N1==0 && (n2-H2)%N2==0) re_exp = 1.0f;

            MKL_LONG index = n2*S2 + n1*S1;
            re_got = data[index].real;
            im_got = data[index].imag;
            double err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
            if (err > maxerr) maxerr = err;
            if (!(err < errthr)) {
                printf(" data[" LI "][" LI "]: ", n2, n1);
                printf(" expected (%.17lg,%.17lg), ", re_exp, im_exp);
                printf(" got (%.17lg,%.17lg), ", re_got, im_got);
                printf(" err %.3lg\n", err);
                isWrong = 1;
            }
        }
    }

    if (isWrong) printf(" Verification FAILED\n");
    else printf(" Verified, maximum error was %.3lg\n", maxerr);

    return isWrong;
}
