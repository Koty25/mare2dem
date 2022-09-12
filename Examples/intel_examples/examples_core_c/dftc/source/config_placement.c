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
! configuration parameter DFTI_PLACEMENT.
! The parameter defines if the result overwrites the input data or not.
!
! Values:
!   DFTI_INPLACE (default) - result overwrites input data
!   DFTI_NOT_INPLACE       - result is placed in a separate array
!
! Note: When storage data types of forward and backward domains are
!       the same, the configuration parameters for the layout of input are
!       also used for the layout of output (e.g. output strides are
!       ignored). Otherwise, both input and output layout shall be defined
!       (for example real transform with conjugate even storage
!       set to DFTI_COMPLEX_COMPLEX).
!
! This example computes an in-place and an out-of-place 2D real transform.
!
!****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mkl_service.h"
#include "mkl_dfti.h"

static void init_r(double *x, MKL_LONG *N, MKL_LONG *S, MKL_LONG *H);
static int verify_c(MKL_Complex16 *x, MKL_LONG *N, MKL_LONG *S, MKL_LONG *H);

/* Define the format to printf MKL_LONG values */
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif

int main(void)
{
    /* Sizes of 2D transform */
    MKL_LONG N[2] = {128, 256};

    /* Arbitrary harmonic used to verify FFT */
    MKL_LONG H[2] = {-1, 2};

    /* Execution status */
    MKL_LONG status = 0;

    /* Pointers to input and output data */
    double *r = NULL;
    MKL_Complex16 *c = NULL;

    /* Strides define data layout in forward and backward domains */
    MKL_LONG rs_inplace[3] = {0, (N[1]/2+1)*2, 1};
    MKL_LONG cs_inplace[3] = {0, (N[1]/2+1), 1};
    MKL_LONG rs_notinplace[3] = {0, N[1], 1};
    MKL_LONG cs_notinplace[3] = {0, (N[1]/2+1), 1};

    DFTI_DESCRIPTOR_HANDLE hand = NULL;

    char version[DFTI_VERSION_LENGTH];


    DftiGetValue(0, DFTI_VERSION, version);
    printf("%s\n", version);
    printf("Example config_placement\n");
    printf("In-place and out-of-place 2D real FFT\n");
    printf(" Configuration parameters:\n");
    printf(" DFTI_PRECISION      = DFTI_DOUBLE\n");
    printf(" DFTI_FORWARD_DOMAIN = DFTI_REAL\n");
    printf(" DFTI_DIMENSION      = 2\n");
    printf(" DFTI_LENGTHS        = { "LI", "LI" }\n", N[0], N[1]);

    printf("======= In-place 2D real FFT =======\n");

    printf("Allocate array for input/output data\n");
    r = mkl_malloc(sizeof(double) * N[0]*(N[1]/2+1)*2, 64);
    if (r == NULL) goto failed;


    printf("Create DFTI descriptor\n");
    status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_REAL, 2, N);
    if (status != DFTI_NO_ERROR) goto failed;

    /*
     * In-place is default configuration, no need to set it.
     *
     * status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_INPLACE);
     * if (status != DFTI_NO_ERROR) goto failed;
     */

    printf("Set CCE storage\n");
    status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE,
                          DFTI_COMPLEX_COMPLEX);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set input  strides = "LI", "LI", "LI"\n",
           rs_inplace[0], rs_inplace[1], rs_inplace[2]);
    status = DftiSetValue(hand, DFTI_INPUT_STRIDES, rs_inplace);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set output strides = "LI", "LI", "LI"\n",
           cs_inplace[0], cs_inplace[1], cs_inplace[2]);
    status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cs_inplace);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Commit descriptor\n");
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Inititalize input for 2D real in-place transform\n");
    init_r(r, N, rs_inplace, H);

    printf("Compute real-to-complex in-place transform\n");
    status = DftiComputeForward(hand, r);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Verify 2D real in-place transform\n");
    status = verify_c((MKL_Complex16*)r, N, cs_inplace, H);
    if (status != 0) goto failed;


    printf("======= Reconfigure for out-of-place 2D real FFT =======\n");

    printf("Allocate array for output data\n");
    c = mkl_malloc(sizeof(MKL_Complex16) * N[0]*(N[1]/2+1), 64);
    if (c == NULL) goto failed;

    printf("Set out-of-place configuration\n");
    status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set input  strides = "LI", "LI", "LI"\n",
           rs_notinplace[0], rs_notinplace[1], rs_notinplace[2]);
    status = DftiSetValue(hand, DFTI_INPUT_STRIDES, rs_notinplace);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set output strides = "LI", "LI", "LI"\n",
           cs_notinplace[0], cs_notinplace[1], cs_notinplace[2]);
    status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cs_notinplace);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Commit descriptor\n");
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Inititalize input for 2D real out-of-place transform\n");
    init_r(r, N, rs_notinplace, H);

    printf("Compute real-to-complex out-of-place transform\n");
    status = DftiComputeForward(hand, r, c);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Verify 2D real in-place transform\n");
    status = verify_c(c, N, cs_notinplace, H);
    if (status != 0) goto failed;

 cleanup:

    printf("Free DFTI descriptor\n");
    DftiFreeDescriptor(&hand);

    printf("Free data arrays\n");
    mkl_free(r);
    mkl_free(c);

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

/* Initialize array data(N) to produce unit peaks at data(H) and data(N-H) */
static void init_r(double *data, MKL_LONG *N, MKL_LONG *S, MKL_LONG *H)
{
    double TWOPI = 6.2831853071795864769, phase, factor;
    MKL_LONG n1, n2, N1, N2, S1, S2, H1, H2, index;

    N1 = N[0];  S1 = S[1];  H1 = H[0];
    N2 = N[1];  S2 = S[2];  H2 = H[1];

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

/* Verify that data has unit peak at H */
static int verify_c(MKL_Complex16 *data, MKL_LONG *N, MKL_LONG *S, MKL_LONG *H)
{
    double err, errthr, maxerr;
    MKL_LONG n1, n2, N1, N2, S1, S2, H1, H2, index;

    N1 = N[0];  S1 = S[1];  H1 = H[0];
    N2 = N[1];  S2 = S[2];  H2 = H[1];

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5 * log((double) N2*N1) / log(2.0) * DBL_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    maxerr = 0.0;
    for (n1 = 0; n1 < N1; n1++) {
        for (n2 = 0; n2 < N2/2+1; n2++) {
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
                printf(" data["LI"]["LI"]: ", n1, n2);
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
