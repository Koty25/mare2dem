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
! An example of using DFTI_NUMBER_OF_TRANSFORMS configuration parameter to do
! multiple SINGLE-precision complex-to-complex in-place 2D FFTs in a single call
! to DftiComputeForward on a CPU or a (GPU) device using the OpenMP target
! (offload) interface of Intel(R) oneMKL DFTI
!********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>

#include "mkl_dfti_omp_offload.h"

static void init(MKL_Complex8 *data, MKL_LONG M,
                 MKL_LONG N1, MKL_LONG N2,
                 MKL_LONG H1, MKL_LONG H2);
static int verify(MKL_Complex8 *data, MKL_LONG M,
                  MKL_LONG N1, MKL_LONG N2,
                  MKL_LONG H1, MKL_LONG H2);

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
    MKL_LONG N[2] = {7, 13};

    // Distance between first elements for multiple transforms
    const MKL_LONG inputDistance = N[0]*N[1];

    // Number of transforms
    MKL_LONG M = 3;

    const MKL_LONG numPoints = M * N[0]*N[1];

    // Arbitrary harmonic used to verify FFT
    MKL_LONG H[2] = {1, -1};

    MKL_LONG status    = 0;
    MKL_LONG statusGPU = 0;
    MKL_Complex8 *data    = NULL;
    MKL_Complex8 *dataGPU = NULL;
    DFTI_DESCRIPTOR_HANDLE descHandle    = NULL;
    DFTI_DESCRIPTOR_HANDLE descHandleGPU = NULL;

    printf("DFTI_DIMENSION      = 2\n");
    printf("DFTI_LENGTHS        = { " LI ", " LI " }\n", N[0], N[1]);
    printf("DFTI_NUMBER_OF_TRANSFORMS = " LI "\n", M);

    printf("Create DFTI descriptor for single-precision 2D FFT\n");
    status = DftiCreateDescriptor(&descHandle, DFTI_SINGLE, DFTI_COMPLEX, 2, N);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Create GPU DFTI descriptor for single-precision 2D FFT\n");
    statusGPU = DftiCreateDescriptor(&descHandleGPU, DFTI_SINGLE, DFTI_COMPLEX, 2, N);
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Set configuration: DFTI_NUMBER_OF_TRANSFORMS\n");
    status = DftiSetValue(descHandle, DFTI_NUMBER_OF_TRANSFORMS, M);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set GPU configuration: DFTI_NUMBER_OF_TRANSFORMS\n");
    statusGPU = DftiSetValue(descHandleGPU, DFTI_NUMBER_OF_TRANSFORMS, M);
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Set configuration: DFTI_INPUT_DISTANCE = " LI "\n", inputDistance);
    status = DftiSetValue(descHandle, DFTI_INPUT_DISTANCE, inputDistance);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set GPU configuration: DFTI_INPUT_DISTANCE = " LI "\n", inputDistance);
    statusGPU = DftiSetValue(descHandleGPU, DFTI_INPUT_DISTANCE, inputDistance);
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    // Output distance is ignored for in-place complex transforms.
    //
    // status = DftiSetValue(descHandle, DFTI_OUTPUT_DISTANCE, inputDistance);
    // if (status != DFTI_NO_ERROR) goto failed;
    //
    // statusGPU = DftiSetValue(descHandleGPU, DFTI_OUTPUT_DISTANCE, inputDistance);
    // if (statusGPU != DFTI_NO_ERROR) goto failed;

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
    data = (MKL_Complex8*)mkl_malloc(numPoints * sizeof(MKL_Complex8), 64);
    if (data == NULL) goto failed;

    printf("Allocate GPU input array\n");
    dataGPU = (MKL_Complex8*)mkl_malloc(numPoints * sizeof(MKL_Complex8), 64);
    if (dataGPU == NULL) goto failed;

    printf("Initialize input for forward FFT\n");
    init(data, M, N[1], N[0], H[1], H[0]);

    printf("Initialize GPU input for forward FFT\n");
    init(dataGPU, M, N[1], N[0], H[1], H[0]);

    printf("Compute forward FFT in-place\n");
    status = DftiComputeForward(descHandle, data);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Compute GPU forward FFT in-place\n");
#pragma omp target data map(tofrom:dataGPU[0:numPoints]) device(devNum)
    {
#pragma omp target variant dispatch use_device_ptr(dataGPU) device(devNum)
        {
            statusGPU = DftiComputeForward(descHandleGPU, dataGPU);
        }
    }
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Verify the result of forward FFT\n");
    status = verify(data, M, N[1], N[0], H[1], H[0]);
    if (status != 0) goto failed;

    printf("Verify the result of GPU forward FFT\n");
    statusGPU = verify(dataGPU, M, N[1], N[0], H[1], H[0]);
    if (statusGPU != 0) goto failed;

    printf("Initialize input for backward FFT\n");
    init(data, M, N[1], N[0], -H[1], -H[0]);

    printf("Initialize GPU input for forward FFT\n");
    init(dataGPU, M, N[1], N[0], -H[1], -H[0]);

    printf("Compute forward FFT in-place\n");
    status = DftiComputeBackward(descHandle, data);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Compute GPU forward FFT in-place\n");
#pragma omp target data map(tofrom:dataGPU[0:numPoints]) device(devNum)
    {
#pragma omp target variant dispatch use_device_ptr(dataGPU) device(devNum)
        {
            statusGPU = DftiComputeBackward(descHandleGPU, dataGPU);
        }
    }
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Verify the result of backward FFT\n");
    status = verify(data, M, N[1], N[0], H[1], H[0]);
    if (status != 0) goto failed;

    printf("Verify the result of GPU backward FFT\n");
    statusGPU = verify(dataGPU, M, N[1], N[0], H[1], H[0]);
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
static float moda(MKL_LONG K, MKL_LONG L, MKL_LONG M)
{
    return (float)(((long long)K * L) % M);
}

// Initialize array with harmonic {H1, H2}
static void init(MKL_Complex8 *data, MKL_LONG M,
                 MKL_LONG N1, MKL_LONG N2,
                 MKL_LONG H1, MKL_LONG H2)
{
    float TWOPI = 6.2831853071795864769f, phase;
    MKL_LONG m, n1, n2, index;

    /* Generalized strides for row-major addressing of data */
    MKL_LONG SM = N2*N1, S1 = 1, S2 = N1;

    for (m = 0; m < M; m++) {
        for (n2 = 0; n2 < N2; n2++) {
            for (n1 = 0; n1 < N1; n1++) {
                phase =  moda(n1, H1, N1) / N1;
                phase += moda(n2, H2, N2) / N2;
                index = m*SM + n2*S2 + n1*S1;
                data[index].real = cosf(TWOPI * phase) / (N2*N1);
                data[index].imag = sinf(TWOPI * phase) / (N2*N1);
            }
        }
    }
}


// Verify that data(n1,n2,m) are unit peaks at H1,H2
static int verify(MKL_Complex8 *data, MKL_LONG M,
                  MKL_LONG N1, MKL_LONG N2,
                  MKL_LONG H1, MKL_LONG H2)
{
    float err, errthr, maxerr;
    MKL_LONG m, n1, n2, index;

    // Generalized strides for row-major addressing of data
    MKL_LONG SM = N2*N1, S1 = 1, S2 = N1;

    errthr = 5.0f * logf((float) N2*N1) / logf(2.0f) * FLT_EPSILON;
    printf(" Verify the result, errthr = %.3lg\n", errthr);

    maxerr = 0.0f;
    for (m = 0; m < M; m++) {
        for (n2 = 0; n2 < N2; n2++) {
            for (n1 = 0; n1 < N1; n1++) {
                float re_exp = 0.0f, im_exp = 0.0f, re_got, im_got;

                if ((n1-H1)%N1==0 && (n2-H2)%N2==0) {
                    re_exp = 1.0f;
                }

                index = m*SM + n2*S2 + n1*S1;
                re_got = data[index].real;
                im_got = data[index].imag;
                err  = fabsf(re_got - re_exp) + fabsf(im_got - im_exp);
                if (err > maxerr) maxerr = err;
                if (!(err < errthr)) {
                    printf(" data[" LI "][" LI "][" LI "]: ", m, n2, n1);
                    printf(" expected (%.7g,%.7g), ", re_exp, im_exp);
                    printf(" got (%.7g,%.7g), ", re_got, im_got);
                    printf(" err %.3lg\n", err);
                    printf(" Verification FAILED\n");
                    return 1;
                }
            }
        }
    }
    printf(" Verified, maximum error was %.3lg\n", maxerr);
    return 0;
}
