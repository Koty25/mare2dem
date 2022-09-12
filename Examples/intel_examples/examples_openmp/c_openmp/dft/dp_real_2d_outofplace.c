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
! An example of DOUBLE-precision real-to-complex out-of-place 2D FFT on a (GPU)
! device using the OpenMP target (offload) interface of Intel(R) oneMKL DFTI
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>

#include "mkl_dfti_omp_offload.h"

static void init_r(double *data, MKL_LONG N1, MKL_LONG N2, MKL_LONG H1, MKL_LONG H2);
static void init_c(MKL_Complex16 *data, MKL_LONG N1, MKL_LONG N2, MKL_LONG H1, MKL_LONG H2);
static int verify_c(MKL_Complex16 *data, MKL_LONG N1, MKL_LONG N2, MKL_LONG H1, MKL_LONG H2);
static int verify_r(double *data, MKL_LONG N1, MKL_LONG N2, MKL_LONG H1, MKL_LONG H2);

// Define the format to printf MKL_LONG values
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif

int main(void)
{
    const int devNum = 0;

    // Size of 2D FFT
    MKL_LONG N1 = 4, N2 = 5;
    MKL_LONG N[2] = {N2, N1};

    const MKL_LONG numRealPts = N2*N1;
    const MKL_LONG numCmplxPts = N2*(N1/2+1);

    // Arbitrary harmonic used to verify FFT
    MKL_LONG H1 = -1, H2 = 2;

    MKL_LONG status    = 0;
    MKL_LONG statusGPU = 0;

    // Pointers to input and output data
    double        *data_real     = NULL;
    double        *dataGPU_real  = NULL;
    MKL_Complex16 *data_cmplx    = NULL;
    MKL_Complex16 *dataGPU_cmplx = NULL;

    DFTI_DESCRIPTOR_HANDLE descHandle = NULL;
    DFTI_DESCRIPTOR_HANDLE descHandleGPU = NULL;

    // Strides describe data layout in real and conjugate-even domain
    MKL_LONG rs[3] = {0, N1, 1};
    MKL_LONG cs[3] = {0, N1/2+1, 1};

    printf("DFTI_LENGTHS                  = {" LI ", " LI "}\n", N2, N1);
    printf("DFTI_PLACEMENT                = DFTI_NOT_INPLACE\n");
    printf("DFTI_CONJUGATE_EVEN_STORAGE   = DFTI_COMPLEX_COMPLEX\n");

    printf("Create DFTI descriptor\n");
    status = DftiCreateDescriptor(&descHandle, DFTI_DOUBLE, DFTI_REAL, 2, N);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Create GPU DFTI descriptor\n");
    statusGPU = DftiCreateDescriptor(&descHandleGPU, DFTI_DOUBLE, DFTI_REAL, 2, N);
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

    printf("Set input strides = ");
    printf("{" LI ", " LI ", " LI "}\n", rs[0], rs[1], rs[2]);
    status = DftiSetValue(descHandle, DFTI_INPUT_STRIDES, rs);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set GPU input strides = ");
    printf("{" LI ", " LI ", " LI "}\n", rs[0], rs[1], rs[2]);
    statusGPU = DftiSetValue(descHandleGPU, DFTI_INPUT_STRIDES, rs);
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Set output strides = ");
    printf("{" LI ", " LI ", " LI "}\n", cs[0], cs[1], cs[2]);
    status = DftiSetValue(descHandle, DFTI_OUTPUT_STRIDES, cs);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set GPU output strides = ");
    printf("{" LI ", " LI ", " LI "}\n", cs[0], cs[1], cs[2]);
    statusGPU = DftiSetValue(descHandleGPU, DFTI_OUTPUT_STRIDES, cs);
    if (statusGPU != DFTI_NO_ERROR) goto failed;

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

    data_real = (double*)mkl_malloc(numRealPts*sizeof(double), 64);
    data_cmplx = (MKL_Complex16*)mkl_malloc(numCmplxPts*sizeof(MKL_Complex16), 64);
    if (data_real == NULL || data_cmplx == NULL) goto failed;

    printf("Allocate GPU data arrays\n");
    dataGPU_real = (double*)mkl_malloc(numRealPts*sizeof(double), 64);
    dataGPU_cmplx = (MKL_Complex16*)mkl_malloc(numCmplxPts*sizeof(MKL_Complex16), 64);
    if (dataGPU_real == NULL || dataGPU_cmplx == NULL) goto failed;

    printf("Initialize data for r2c FFT\n");
    init_r(data_real, N1, N2, H1, H2);

    printf("Initialize GPU data for r2c FFT\n");
    init_r(dataGPU_real, N1, N2, H1, H2);

    printf("Compute real-to-complex FFT\n");
    status = DftiComputeForward(descHandle, data_real, data_cmplx);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Compute GPU forward FFT\n");
#pragma omp target data map(to:dataGPU_real[0:numRealPts]) \
                        map(from:dataGPU_cmplx[0:numCmplxPts]) device(devNum)
    {
#pragma omp target variant dispatch use_device_ptr(dataGPU_real, dataGPU_cmplx) \
                                    device(devNum)
        {
            statusGPU = DftiComputeForward(descHandleGPU, dataGPU_real, dataGPU_cmplx);
        }
    }
    if (statusGPU != DFTI_NO_ERROR) goto failed;

    printf("Verify the result\n");
    status = verify_c(data_cmplx, N1, N2, H1, H2);
    if (status != 0) goto failed;

    printf("Verify the GPU result\n");
    statusGPU = verify_c(dataGPU_cmplx, N1, N2, H1, H2);
    if (statusGPU != 0) goto failed;

 cleanup:

    printf("Free DFTI descriptor\n");
    DftiFreeDescriptor(&descHandle);

    printf("Free GPU DFTI descriptor\n");
    DftiFreeDescriptor(&descHandleGPU);

    printf("Free data arrays\n");
    mkl_free(data_real);
    mkl_free(data_cmplx);

    printf("Free GPU data arrays\n");
    mkl_free(dataGPU_real);
    mkl_free(dataGPU_cmplx);

    {
        const MKL_LONG statusBoth = status | statusGPU;
        printf("TEST %s\n", statusBoth == 0 ? "PASSED" : "FAILED");
        return statusBoth;
    }

 failed:
    printf(" ERROR, status = " LI ", statusGPU = " LI "\n", status, statusGPU);
    goto cleanup;
}

// Compute (K*L)%M accurately
static double moda(MKL_LONG K, MKL_LONG L, MKL_LONG M)
{
    return (double)(((long long)K * L) % M);
}

// Initialize array data(N) to produce unit peaks at data(H) and data(N-H)
static void init_r(double *data, MKL_LONG N1, MKL_LONG N2, MKL_LONG H1, MKL_LONG H2)
{
    double TWOPI = 6.2831853071795864769, phase, factor;
    MKL_LONG n1, n2, index;

    // Generalized strides for row-major addressing of data
    MKL_LONG S1 = 1, S2 = N1;

    factor = (2*(N1-H1)%N1==0 && 2*(N2-H2)%N2==0) ? 1.0 : 2.0;
    for (n2 = 0; n2 < N2; n2++) {
        for (n1 = 0; n1 < N1; n1++) {
            phase  = moda(n1, H1, N1) / N1;
            phase += moda(n2, H2, N2) / N2;
            index = n2*S2 + n1*S1;
            data[index] = factor * cos(TWOPI * phase) / (N2*N1);
        }
    }
}

// Verify that x has unit peak at H
static int verify_c(MKL_Complex16 *data, MKL_LONG N1, MKL_LONG N2, MKL_LONG H1, MKL_LONG H2)
{
    double err, errthr, maxerr;
    MKL_LONG n1, n2, index;

    // Generalized strides for row-major addressing of data
    MKL_LONG S1 = 1, S2 = N1/2+1;

    errthr = 2.5 * log((double) N2*N1) / log(2.0) * DBL_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    maxerr = 0.0;
    for (n2 = 0; n2 < N2; n2++) {
        for (n1 = 0; n1 < N1/2+1; n1++) {
            double re_exp = 0.0, im_exp = 0.0, re_got, im_got;

            if ((( n1-H1)%N1==0 && ( n2-H2)%N2==0) ||
                ((-n1-H1)%N1==0 && (-n2-H2)%N2==0)
            ) {
                re_exp = 1.0;
            }

            index = n2*S2 + n1*S1;
            re_got = data[index].real;
            im_got = data[index].imag;
            err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
            if (err > maxerr) maxerr = err;
            if (!(err < errthr)) {
                printf(" data[" LI "][" LI "]: ", n2, n1);
                printf(" expected (%.17lg,%.17lg), ", re_exp, im_exp);
                printf(" got (%.17lg,%.17lg), ", re_got, im_got);
                printf(" err %.3lg\n", err);
                printf(" Verification FAILED\n");
                return 1;
            }
        }
    }
    printf(" Verified,  maximum error was %.3lg\n", maxerr);
    return 0;
}
