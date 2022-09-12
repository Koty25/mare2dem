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
! An example of using DFTI_NUMBER_OF_USER_THREADS configuration parameter.
! The parameter specifies how many user threads (OS threads or OpenMP threads)
! share the descriptor for computation of FFT.
!
! Values:
! Any positive integer (default 1)
!
!****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mkl_service.h"
#include "mkl_dfti.h"

static void init(MKL_Complex16 *data, int M,
                 int N1, int N2, int N3,
                 int H1, int H2, int H3);
static int verify(MKL_Complex16 *data, int M,
                  int N1, int N2, int N3,
                  int H1, int H2, int H3);

/* Define the format to printf MKL_LONG values */
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif


int main(void)
{
    /* Sizes of 3D transform */
    int N1 = 5, N2 = 5, N3 = 5;
    MKL_LONG N[3] = {N3, N2, N1};

    /* Number of transforms to compute */
    int M = 100;

    /* Number of user threads sharing the descriptor */
    int NUT = 4;

    /* Arbitrary harmonic used to verify FFT */
    int H1 = -2, H2 = -3, H3 = -4;

    /* Execution status */
    MKL_LONG status = 0;

    /* Pointer to input/output data */
    MKL_Complex16 *data = NULL;

    DFTI_DESCRIPTOR_HANDLE hand = NULL;

    char version[DFTI_VERSION_LENGTH];

    /* Local variables */
    int m;
    MKL_LONG thr_status;

    DftiGetValue(0, DFTI_VERSION, version);
    printf("%s\n", version);
    printf("Example config_number_of_user_threads\n");
    printf("Multiple 3D in-place FFT using shared descriptor\n");
    printf("Configuration parameters:\n");
    printf(" DFTI_PRECISION              = DFTI_DOUBLE\n");
    printf(" DFTI_FORWARD_DOMAIN         = DFTI_COMPLEX\n");
    printf(" DFTI_DIMENSION              = 3\n");
    printf(" DFTI_LENGTHS                = {%d, %d, %d}\n", N3, N2, N1);
    printf(" Number of transforms      M = %d\n", M);
    printf(" DFTI_NUMBER_OF_USER_THREADS = %d\n", NUT);

    printf("Create DFTI descriptor\n");
    status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_COMPLEX, 3, N);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Set configuration: number of user threads\n");
    status = DftiSetValue(hand, DFTI_NUMBER_OF_USER_THREADS, NUT);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Commit descriptor\n");
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Allocate input/output array\n");
    data = (MKL_Complex16*)mkl_malloc(M * N3*N2*N1 * sizeof(MKL_Complex16), 64);
    if (data == NULL) goto failed;

    printf("Initialize input\n");
    init(data, M, N1, N2, N3, H1, H2, H3);

    printf("Compute forward transform by parallel user threads\n");
#if defined(_OPENMP)
#pragma omp parallel for shared(hand, data) private(m, thr_status)
#endif
    for (m = 0; m < M; ++m) {

        /*
         * If the actual size of parallel team of threads sharing 'hand' is
         * greater than 'NUT', the number of user threads set in the
         * descriptor, then the performance may be negatively affected.
         */

        thr_status = DftiComputeForward(hand, data + m * N3*N2*N1);

        /* Update global status only if this thread fails */
        if (thr_status != DFTI_NO_ERROR) status = thr_status;
    }
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Verify the result\n");
    status = verify(data, M, N1, N2, N3, H1, H2, H3);
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
static double moda(int K, int L, int M)
{
    return (double)(((long long)K * L) % M);
}

/* Initialize array with harmonic {H1, H2, H3} */
static void init(MKL_Complex16 *data, int M,
                 int N1, int N2, int N3,
                 int H1, int H2, int H3)
{
    double TWOPI = 6.2831853071795864769, phase;
    int m, n1, n2, n3, index;

    /* Generalized strides for row-major addressing of data */
    int SM = N3*N2*N1, S3 = N2*N1, S2 = N1, S1 = 1;

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

/* Verify that data(m,n1,n2,n3) are unit peaks at H1,H2,H3 */
static int verify(MKL_Complex16 *data, int M,
                  int N1, int N2, int N3,
                  int H1, int H2, int H3)
{
    double err, errthr, maxerr;
    int m, n1, n2, n3, index;

    /* Generalized strides for row-major addressing of data */
    int SM = N3*N2*N1, S3 = N2*N1, S2 = N1, S1 = 1;

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
                        printf(" data[%i][%i][%i][%i]: ", m, n3, n2, n1);
                        printf(" expected (%.17lg,%.17lg), ", re_exp, im_exp);
                        printf(" got (%.17lg,%.17lg), ", re_got,im_got);
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
