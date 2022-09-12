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
! An example of using Intel(R) Math Kernel Library (Intel(R) MKL) DFTI
! configuration parameter DFTI_COMPLEX_STORAGE.
! The parameter is used to select how complex data is laid out in memory
! (DFTI_FORWARD_DOMAIN = DFTI_COMPLEX is required).
!
! Values:
! DFTI_COMPLEX_COMPLEX (default) - use array of complex data
! DFTI_REAL_REAL                 - use two arrays of real data (split complex)
!
!****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mkl_service.h"
#include "mkl_dfti.h"

static void init(double *data_re, double *data_im,
                 MKL_LONG N1, MKL_LONG N2,
                 MKL_LONG H1, MKL_LONG H2);
static int verify(double *data_re, double *data_im,
                  MKL_LONG N1, MKL_LONG N2,
                  MKL_LONG H1, MKL_LONG H2);

/* Define the format to printf MKL_LONG values */
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif


int main(void)
{
    /* Sizes of 2D transform */
    MKL_LONG N[] = {7, 13};

    /* Arbitrary harmonic to verify FFT */
    MKL_LONG H[] = {1, -1};

    /* Execution status */
    MKL_LONG status = 0;

    /* Pointers to input/output data */
    double *data_re = NULL, *data_im = NULL;

    DFTI_DESCRIPTOR_HANDLE hand = NULL;

    char version[DFTI_VERSION_LENGTH];


    DftiGetValue(0, DFTI_VERSION, version);
    printf("%s\n", version);
    printf("Example config_complex_storage\n");
    printf("Forward and backward split complex in-place 2D transform\n");
    printf("Configuration parameters:\n");
    printf(" DFTI_PRECISION       = DFTI_DOUBLE\n");
    printf(" DFTI_FORWARD_DOMAIN  = DFTI_COMPLEX\n");
    printf(" DFTI_DIMENSION       = 2\n");
    printf(" DFTI_LENGTHS         = { "LI", "LI" }\n", N[0], N[1]);
    printf(" DFTI_COMPLEX_STORAGE = DFTI_REAL_REAL\n");

    printf("Create DFTI descriptor\n");
    status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_COMPLEX, 2, N);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set configuration: split-complex\n");
    status = DftiSetValue(hand, DFTI_COMPLEX_STORAGE, DFTI_REAL_REAL);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Commit the descriptor\n");
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Allocate arrays\n");
    data_re = (double*)mkl_malloc(N[0]*N[1]*sizeof(double), 64);
    data_im = (double*)mkl_malloc(N[0]*N[1]*sizeof(double), 64);
    if (data_re == NULL || data_im == NULL) goto failed;

    printf("Initialize input data\n");
    init(data_re, data_im, N[1], N[0], H[1], H[0]);

    printf("Compute forward transform\n");
    status = DftiComputeForward(hand, data_re, data_im);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Verify the result\n");
    status = verify(data_re, data_im, N[1], N[0], H[1], H[0]);
    if (status != 0) goto failed;

    printf("Initialize input for backward transform\n");
    init(data_re, data_im, N[1], N[0], -H[1], -H[0]);

    printf("Compute backward transform\n");
    status = DftiComputeBackward(hand, data_re, data_im);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Verify the result of backward transform\n");
    status = verify(data_re, data_im, N[1], N[0], H[1], H[0]);
    if (status != 0) goto failed;

 cleanup:

    printf("Free DFTI descriptor\n");
    DftiFreeDescriptor(&hand);

    printf("Free data arrays\n");
    mkl_free(data_re);
    mkl_free(data_im);

    printf("TEST %s\n", (status == 0) ? "PASSED" : "FAILED");
    return status;

 failed:
    printf(" ERROR, status = "LI"\n", status);
    status = 1;
    goto cleanup;
}

/* Compute (K*L)%M accurately */
static double moda(MKL_LONG K, MKL_LONG L, MKL_LONG M)
{
    return (double)(((long long)K * L) % M);
}

/* Initialize array data(N) to produce unit peaks at data(H) */
static void init(double *data_re, double *data_im,
                 MKL_LONG N1, MKL_LONG N2,
                 MKL_LONG H1, MKL_LONG H2)
{
    double TWOPI = 6.2831853071795864769, phase;
    MKL_LONG n1, n2, index;

    /* Generalized strides for row-major addressing of data */
    int S1 = 1, S2 = N1;

    for (n2 = 0; n2 < N2; n2++) {
        for (n1 = 0; n1 < N1; n1++) {
            phase  = moda(n1, H1, N1) / N1;
            phase += moda(n2, H2, N2) / N2;
            index = n2*S2 + n1*S1;
            data_re[index] = cos(TWOPI * phase) / (N2*N1);
            data_im[index] = sin(TWOPI * phase) / (N2*N1);
        }
    }
}

/* Verify that data has unit peak at H */
static int verify(double *data_re, double *data_im,
                  MKL_LONG N1, MKL_LONG N2,
                  MKL_LONG H1, MKL_LONG H2)
{
    double err, errthr, maxerr;
    MKL_LONG n1, n2, index;

    /* Generalized strides for row-major addressing of data */
    int S1 = 1, S2 = N1;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 5.0 * log((double) N2*N1) / log(2.0) * DBL_EPSILON;
    printf(" Verify the result, errthr = %.3lg\n", errthr);

    maxerr = 0.0;
    for (n2 = 0; n2 < N2; n2++) {
        for (n1 = 0; n1 < N1; n1++) {
            double re_exp = 0.0, im_exp = 0.0, re_got, im_got;

            if ((n1-H1)%N1==0 && (n2-H2)%N2==0) re_exp = 1.0;

            index = n2*S2 + n1*S1;
            re_got = data_re[index];
            im_got = data_im[index];
            err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
            if (err > maxerr) maxerr = err;
            if (!(err < errthr)) {
                printf(" data_re["LI"]["LI"]: ", n2, n1);
                printf(" expected (%.17lg,%.17lg), ", re_exp, im_exp);
                printf(" got (%.17lg,%.17lg), ", re_got, im_got);
                printf(" err %.3lg\n", err);
                printf(" Verification FAILED\n");
                return 1;
            }
        }
    }
    printf(" Verified, maximum error was %.3lg\n", maxerr);
    return 0;
}
