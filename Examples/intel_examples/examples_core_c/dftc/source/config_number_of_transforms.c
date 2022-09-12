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
! An example of using DFTI_NUMBER_OF_TRANSFORMS configuration parameter.
! The parameter defines how many identical transforms are computed by one call
! of DftiComputeForward or DftiComputeBackward function.
!
! Values:
! A positive integer (default 1)
!
!****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mkl_service.h"
#include "mkl_dfti.h"

static void init(MKL_Complex16 *data, MKL_LONG M,
                 MKL_LONG N1, MKL_LONG N2, MKL_LONG N3,
                 MKL_LONG H1, MKL_LONG H2, MKL_LONG H3);
static int verify(MKL_Complex16 *data, MKL_LONG M,
                  MKL_LONG N1, MKL_LONG N2, MKL_LONG N3,
                  MKL_LONG H1, MKL_LONG H2, MKL_LONG H3);

/* Define the format to printf MKL_LONG values */
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif

int main(void)
{
    /* Sizes of 3D transform */
    MKL_LONG N[3] = {7, 13, 5};

    /* Number of transforms */
    MKL_LONG M = 4;

    /* Arbitrary harmonic used to verify FFT */
    MKL_LONG H[3] = {1, -1, -2};

    /* Execution status */
    MKL_LONG status = 0;

    /* Pointer to input/output data */
    MKL_Complex16 *data = NULL;

    DFTI_DESCRIPTOR_HANDLE hand = NULL;

    /* Distance between first elements for multiple transforms */
    MKL_LONG dist;

    char version[DFTI_VERSION_LENGTH];


    DftiGetValue(0, DFTI_VERSION, version);
    printf("%s\n", version);
    printf("Example config_number_of_transforms\n");
    printf("Multiple in-place 3D FFT\n");
    printf("Configuration parameters:\n");
    printf(" DFTI_PRECISION      = DFTI_DOUBLE\n");
    printf(" DFTI_FORWARD_DOMAIN = DFTI_COMPLEX\n");
    printf(" DFTI_DIMENSION      = 3\n");
    printf(" DFTI_LENGTHS        = { "LI", "LI", "LI" }\n", N[0], N[1], N[2]);
    printf(" DFTI_NUMBER_OF_TRANSFORMS = "LI"\n", M);

    printf("Create DFTI descriptor for double-precision 3D transform\n");
    status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_COMPLEX, 3, N);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set configuration: DFTI_NUMBER_OF_TRANSFORMS\n");
    status = DftiSetValue(hand, DFTI_NUMBER_OF_TRANSFORMS, M);
    if (status != DFTI_NO_ERROR) goto failed;

    dist = N[0]*N[1]*N[2];

    printf("Set configuration: DFTI_INPUT_DISTANCE = "LI"\n", dist);
    status = DftiSetValue(hand, DFTI_INPUT_DISTANCE, dist);
    if (status != DFTI_NO_ERROR) goto failed;

    /*
     * Output distance is ignored for in-place complex transforms.
     *
     * status = DftiSetValue(hand, DFTI_INPUT_DISTANCE, dist);
     * if (status != DFTI_NO_ERROR) goto failed;
     */

    printf("Commit descriptor\n");
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) goto failed;

    /* Allocate input array */
    int count = M * N[0]*N[1]*N[2];
    data = (MKL_Complex16*)mkl_malloc(count * sizeof(MKL_Complex16), 64);
    if (data == NULL) goto failed;

    printf("Initialize input for forward transform\n");
    init(data, M, N[2], N[1], N[0], H[2], H[1], H[0]);

    printf("Compute forward transform\n");
    status = DftiComputeForward(hand, data);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Verify the result\n");
    status = verify(data, M, N[2], N[1], N[0], H[2], H[1], H[0]);
    if (status != 0) goto failed;

    printf("Initialize input for backward transform\n");
    init(data, M, N[2], N[1], N[0], -H[2], -H[1], -H[0]);

    printf("Compute backward transform\n");
    status = DftiComputeBackward(hand, data);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Verify the result\n");
    status = verify(data, M, N[2], N[1], N[0], H[2], H[1], H[0]);
    if (status != 0) goto failed;

 cleanup:

    printf("Free DFTI descriptor\n");
    DftiFreeDescriptor(&hand);

    printf("Free data array\n");
    mkl_free(data);

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

/* Initialize array with harmonic {H1, H2, H3} */
static void init(MKL_Complex16 *data, MKL_LONG M,
                 MKL_LONG N1, MKL_LONG N2, MKL_LONG N3,
                 MKL_LONG H1, MKL_LONG H2, MKL_LONG H3)
{
    double TWOPI = 6.2831853071795864769, phase;
    MKL_LONG m, n1, n2, n3, index;

    /* Generalized strides for row-major addressing of data */
    MKL_LONG SM = N3*N2*N1, S1 = 1, S2 = N1, S3 = N2*N1;

    for (m = 0; m < M; m++) {
        for (n3 = 0; n3 < N3; n3++) {
            for (n2 = 0; n2 < N2; n2++) {
                for (n1 = 0; n1 < N1; n1++) {
                    phase =  moda(n1, H1, N1) / N1;
                    phase += moda(n2, H2, N2) / N2;
                    phase += moda(n3, H3, N3) / N3;
                    index = m*SM + n3*S3 + n2*S2 + n1*S1;
                    data[index].real = cos(TWOPI * phase) / (N3*N2*N1);
                    data[index].imag = sin(TWOPI * phase) / (N3*N2*N1);
                }
            }
        }
    }
}


/* Verify that data(n1,n2,n3,m) are unit peaks at H1,H2,H3 */
static int verify(MKL_Complex16 *data, MKL_LONG M,
                  MKL_LONG N1, MKL_LONG N2, MKL_LONG N3,
                  MKL_LONG H1, MKL_LONG H2, MKL_LONG H3)
{
    double err, errthr, maxerr;
    MKL_LONG m, n1, n2, n3, index;

    /* Generalized strides for row-major addressing of data */
    MKL_LONG SM = N3*N2*N1, S1 = 1, S2 = N1, S3 = N2*N1;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 5.0 * log((double) N3*N2*N1) / log(2.0) * DBL_EPSILON;
    printf(" Verify the result, errthr = %.3lg\n", errthr);

    maxerr = 0.0;
    for (m = 0; m < M; m++) {
        for (n3 = 0; n3 < N3; n3++) {
            for (n2 = 0; n2 < N2; n2++) {
                for (n1 = 0; n1 < N1; n1++) {
                    double re_exp = 0.0, im_exp = 0.0, re_got, im_got;

                    if ((n1-H1)%N1==0 && (n2-H2)%N2==0 && (n3-H3)%N3==0) {
                        re_exp = 1.0;
                    }

                    index = m*SM + n3*S3 + n2*S2 + n1*S1;
                    re_got = data[index].real;
                    im_got = data[index].imag;
                    err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
                    if (err > maxerr) maxerr = err;
                    if (!(err < errthr)) {
                        printf(" data["LI"]["LI"]["LI"]["LI"]: ", m, n3, n2, n1);
                        printf(" expected (%.17lg,%.17lg), ", re_exp, im_exp);
                        printf(" got (%.17lg,%.17lg), ", re_got, im_got);
                        printf(" err %.3lg\n", err);
                        printf(" Verification FAILED\n");
                        return 1;
                    }
                }
            }
        }
    }
    printf(" Verified, maximum error was %.3lg\n", maxerr);
    return 0;
}
