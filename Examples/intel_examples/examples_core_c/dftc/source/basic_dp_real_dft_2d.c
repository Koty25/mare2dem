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
! A simple example of double-precision real-to-complex out-of-place 2D
! FFT using Intel(R) Math Kernel Library (Intel(R) MKL) DFTI
!
!****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mkl_service.h"
#include "mkl_dfti.h"

static void init_r(double *data, int N1, int N2, int H1, int H2);
static void init_c(MKL_Complex16 *data, int N1, int N2, int H1, int H2);
static int verify_c(MKL_Complex16 *data, int N1, int N2, int H1, int H2);
static int verify_r(double *data, int N1, int N2, int H1, int H2);

/* Define the format to printf MKL_LONG values */
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif

int main(void)
{
    /* Size of 2D transform */
    int N1 = 4, N2 = 5;
    MKL_LONG N[2] = {N2, N1};

    /* Arbitrary harmonic used to verify FFT */
    int H1 = -1, H2 = 2;

    /* Execution status */
    MKL_LONG status = 0;

    /* Pointers to input and output data */
    double *data_real = NULL;
    MKL_Complex16 *data_cmplx = NULL;

    /* Strides describe data layout in real and conjugate-even domain */
    MKL_LONG rs[3] = {0, N1, 1};
    MKL_LONG cs[3] = {0, N1/2+1, 1};

    DFTI_DESCRIPTOR_HANDLE hand = NULL;

    char version[DFTI_VERSION_LENGTH];

    DftiGetValue(0, DFTI_VERSION, version);
    printf("%s\n", version);

    printf("Example basic_dp_real_dft_2d\n");
    printf("Forward-Backward double-precision 2D real out-of-place FFT\n");
    printf("Configuration parameters:\n");
    printf(" DFTI_PRECISION                = DFTI_DOUBLE\n");
    printf(" DFTI_FORWARD_DOMAIN           = DFTI_REAL\n");
    printf(" DFTI_DIMENSION                = 2\n");
    printf(" DFTI_LENGTHS                  = {%i, %i}\n", N2, N1);
    printf(" DFTI_PLACEMENT                = DFTI_NOT_INPLACE\n");
    printf(" DFTI_CONJUGATE_EVEN_STORAGE   = DFTI_COMPLEX_COMPLEX\n");

    printf("Create DFTI descriptor\n");
    status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_REAL, 2, N);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set configuration: out-of-place\n");
    status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set configuration: CCE storage\n");
    status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE,
                          DFTI_COMPLEX_COMPLEX);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set input  strides = ");
    printf("{"LI", "LI", "LI"}\n", rs[0], rs[1], rs[2]);
    status = DftiSetValue(hand, DFTI_INPUT_STRIDES, rs);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set output strides = ");
    printf("{"LI", "LI", "LI"}\n", cs[0], cs[1], cs[2]);
    status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cs);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Commit the descriptor\n");
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Allocate data arrays\n");
    data_real = (double*)mkl_malloc(N2*N1*sizeof(double), 64);
    if (data_real == NULL) goto failed;
    data_cmplx = (MKL_Complex16*)mkl_malloc(N2*(N1/2+1)*sizeof(MKL_Complex16),
                                            64);
    if (data_cmplx == NULL) goto failed;

    printf("Initialize data for r2c transform\n");
    init_r(data_real, N1, N2, H1, H2);

    printf("Compute real-to-complex transform\n");
    status = DftiComputeForward(hand, data_real, data_cmplx);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Verify the result\n");
    status = verify_c(data_cmplx, N1, N2, H1, H2);
    if (status != 0) goto failed;

    printf("Reconfigure DFTI descriptor for backward transform\n");

    printf("Set input  strides = ");
    printf("{"LI", "LI", "LI"}\n", cs[0], cs[1], cs[2]);
    status = DftiSetValue(hand, DFTI_INPUT_STRIDES, cs);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set output strides = ");
    printf("{"LI", "LI", "LI"}\n", rs[0], rs[1], rs[2]);
    status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, rs);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Commit the descriptor\n");
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Initialize data for c2r transform\n");
    init_c(data_cmplx, N1, N2, H1, H2);

    printf("Compute backward transform\n");
    status = DftiComputeBackward(hand, data_cmplx, data_real);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Verify the result of the backward transform\n");
    status = verify_r(data_real, N1, N2, H1, H2);
    if (status != 0) goto failed;

 cleanup:

    printf("Free DFTI descriptor\n");
    DftiFreeDescriptor(&hand);

    printf("Free data arrays\n");
    mkl_free(data_real);
    mkl_free(data_cmplx);

    printf("TEST %s\n", (status == 0) ? "PASSED" : "FAILED");
    return status;

 failed:
    printf(" ERROR, status = "LI"\n", status);
    status = 1;
    goto cleanup;
}

/* Compute (K*L)%M accurately */
static double moda(int K, int L, int M)
{
    return (double)(((long long)K * L) % M);
}

/* Initialize array data(N) to produce unit peaks at data(H) and data(N-H) */
static void init_r(double *data, int N1, int N2, int H1, int H2)
{
    double TWOPI = 6.2831853071795864769, phase, factor;
    int n1, n2, index;

    /* Generalized strides for row-major addressing of data */
    int S1 = 1, S2 = N1;

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

/* Verify that x has unit peak at H */
static int verify_c(MKL_Complex16 *data, int N1, int N2, int H1, int H2)
{
    double err, errthr, maxerr;
    int n1, n2, index;

    /* Generalized strides for row-major addressing of data */
    int S1 = 1, S2 = N1/2+1;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
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
                printf(" data[%i][%i]: ", n2, n1);
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

/* Initialize array data(N) to produce unit peak at data(H) */
static void init_c(MKL_Complex16 *data, int N1, int N2, int H1, int H2)
{
    double TWOPI = 6.2831853071795864769, phase;
    int n1, n2, index;

    /* Generalized strides for row-major addressing of data */
    int S1 = 1, S2 = N1/2+1;

    for (n2 = 0; n2 < N2; n2++) {
        for (n1 = 0; n1 < N1/2+1; n1++) {
            phase  = moda(n1, H1, N1) / N1;
            phase += moda(n2, H2, N2) / N2;
            index = n2*S2 + n1*S1;
            data[index].real =  cos(TWOPI * phase) / (N2*N1);
            data[index].imag = -sin(TWOPI * phase) / (N2*N1);
        }
    }
}

/* Verify that x has unit peak at H */
static int verify_r(double *data, int N1, int N2, int H1, int H2)
{
    double err, errthr, maxerr;
    int n1, n2, index;

    /* Generalized strides for row-major addressing of data */
    int S1 = 1, S2 = N1;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5 * log((double) N2*N1) / log(2.0) * DBL_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    maxerr = 0.0;
    for (n2 = 0; n2 < N2; n2++) {
        for (n1 = 0; n1 < N1; n1++) {
            double re_exp = 0.0, re_got;

            if ((n1-H1)%N1==0 && (n2-H2)%N2==0) re_exp = 1.0;

            index = n2*S2 + n1*S1;
            re_got = data[index];
            err  = fabs(re_got - re_exp);
            if (err > maxerr) maxerr = err;
            if (!(err < errthr)) {
                printf(" data[%i][%i]: ", n2, n1);
                printf(" expected %.17lg, ", re_exp);
                printf(" got %.17lg, ", re_got);
                printf(" err %.3lg\n", err);
                printf(" Verification FAILED\n");
                return 1;
            }
        }
    }
    printf(" Verified,  maximum error was %.3lg\n", maxerr);
    return 0;
}
