/*******************************************************************************
* Copyright 2012-2020 Intel Corporation.
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
! An example of using DFTI_THREAD_LIMIT configuration parameter.
! The parameter specifies maximum number of OpenMP threads FFT can use.
!
! Values:
!   0 (default) = use number of threads specified by
!                 mkl_[domain_]set_num_threads()
!   Any positive integer N = use not more than N threads
!
!****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mkl_dfti.h"
#include "mkl_service.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

#define REAL double
typedef struct { REAL real,imag; } COMPLEX;

static void init(COMPLEX *data,
                 int N1, int N2, int N3,
                 int H1, int H2, int H3);
static int verify(COMPLEX *data,
                  int N1, int N2, int N3,
                  int H1, int H2, int H3);

/* Define the format to printf MKL_LONG values */
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif


int run_dft
(
    int tid,                    /* Id of this thread */
    int tlimit,                 /* Thread limit */
    int N1, int N2, int N3,     /* Sizes of 3D transform */
    int H1, int H2, int H3      /* Arbitrary harmonic used to verify FFT */
)
{
    /* Execution status */
    MKL_LONG status = 0;

    /* Pointer to input/output data */
    COMPLEX *data = NULL;

    DFTI_DESCRIPTOR_HANDLE hand = NULL;

    printf("Thread %i: 3D in-place FFT on %i threads\n", tid, tlimit);
    printf("Thread %i: Create DFTI descriptor for %ix%ix%i FFT\n", tid,
           N3, N2, N1);
    {
        MKL_LONG N[3] = {N3, N2, N1};
        status = DftiCreateDescriptor(&hand,
                                      sizeof(REAL)==sizeof(float)
                                      ? DFTI_SINGLE : DFTI_DOUBLE,
                                      DFTI_COMPLEX, 3, N);
        if (status != DFTI_NO_ERROR) goto failed;
    }

    printf("Thread %i: Set thread limit %i\n", tid, tlimit);
    status = DftiSetValue(hand, DFTI_THREAD_LIMIT, tlimit);
    if (status != DFTI_NO_ERROR) goto failed;

    /* If tlimit > 1 check if we linked with sequential Intel(R) Math Kernel Library (Intel(R) MKL) */
    if (tlimit > 1) {
        /* Get thread limit of uncommitted descriptor */
        MKL_LONG tl;

        status = DftiGetValue(hand, DFTI_THREAD_LIMIT, &tl);
        if (status != DFTI_NO_ERROR) goto failed;

        printf("Thread %i: uncommitted descriptor thread limit %i %s\n",
               tid, (int)tl, tl==1 ? "(sequential MKL)" : "");
    }

    printf("Thread %i: commit descriptor\n", tid);
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) goto failed;

    /* Get thread limit of committed descriptor */
    {
        MKL_LONG tl;
        status = DftiGetValue(hand, DFTI_THREAD_LIMIT, &tl);
        if (status != DFTI_NO_ERROR) goto failed;

        printf("Thread %i: committed descriptor thread limit %i\n", tid,
               (int)tl);
    }

    printf("Thread %i: Allocate input/output array\n", tid);
    data = (COMPLEX*)mkl_malloc(N3*N2*N1 * sizeof(COMPLEX), 64);
    if (data == NULL) goto failed;

    printf("Thread %i: Initialize input\n", tid);
    init(data, N1, N2, N3, H1, H2, H3);

    printf("Thread %i: Compute forward transform\n", tid);
    status = DftiComputeForward(hand, data);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Thread %i: Verify the result\n", tid);
    status = verify(data, N1, N2, N3, H1, H2, H3);
    if (status != 0) goto failed;

 cleanup:

    printf("Thread %i: Free DFTI descriptor\n", tid);
    DftiFreeDescriptor(&hand);

    printf("Thread %i: Free data array\n", tid);
    mkl_free(data);

    printf("Thread %i: Subtest %s\n", tid, (status == 0) ? "Passed" : "Failed");
    return status;

 failed:
    printf("Thread %i: ERROR, status = "LI"\n", tid, status);
    status = 1;
    goto cleanup;
}

int main()
{
    int failed = 0;
    char version[DFTI_VERSION_LENGTH];


#if defined(_OPENMP)
    int NUT = 2; /* Number of parallel user threads */

    /* Enable nested parallel OpenMP sections (maybe oversubscribed) */
    omp_set_nested(1);
    omp_set_dynamic(0);
#endif

    /* Enable threading of Intel MKL called from OpenMP parallel sections */
    MKL_Set_Dynamic(0);

    DftiGetValue(0, DFTI_VERSION, version);
    printf("%s\n", version);

    printf("Example config_thread_limit\n");

#if defined(_OPENMP)
    printf("Run parallel FFTs on %i parallel threads\n",NUT);
#pragma omp parallel num_threads(NUT)
#else
    printf("Run parallel FFT on a single thread\n");
#endif
    {
        /* Two threads running DFT on different number of threads */
        int err;
#if defined(_OPENMP)
        int me = omp_get_thread_num();
        int team = omp_get_num_threads();
#else
        int me = 0;
        int team = 1;
#endif
        if (me == 0)
            printf("Thread %i: parallel team is %i threads\n", me, team);

        if (me) {
            err = run_dft(me, 2, 100, 200, 300, -1, -2, -3);
        } else {
            err = run_dft(me, 3, 300, 100, 200, -1, -2, -3);
        }
        if (err) {
            failed = err;
        }
    }

    printf("TEST %s\n", failed ? "FAILED" : "PASSED");
    return failed;
}

/* Compute (K*L)%M accurately */
static double moda(int K, int L, int M)
{
    return (double)(((long long)K * L) % M);
}

/* Initialize array with harmonic {H1, H2, H3} */
static void init(COMPLEX *x,
                 int N1, int N2, int N3,
                 int H1, int H2, int H3)
{
    double TWOPI = 6.2831853071795864769, phase;
    int n1, n2, n3, index;

    /* Generalized strides for row-major addressing of x */
    int S3 = N2*N1, S2 = N1, S1 = 1;

    for (n3 = 0; n3 < N3; n3++) {
        for (n2 = 0; n2 < N2; n2++) {
            for (n1 = 0; n1 < N1; n1++) {
                phase =  moda(n1, H1, N1) / N1;
                phase += moda(n2, H2, N2) / N2;
                phase += moda(n3, H3, N3) / N3;
                index = n3*S3 + n2*S2 + n1*S1;
                x[index].real = cos(TWOPI * phase) / (N3*N2*N1);
                x[index].imag = sin(TWOPI * phase) / (N3*N2*N1);
            }
        }
    }
}

/* Verify that data(n1,n2,n3) are unit peaks at H1,H2,H3 */
static int verify(COMPLEX *data,
                  int N1, int N2, int N3,
                  int H1, int H2, int H3)
{
    double err, errthr, maxerr;
    int n1, n2, n3, index;

    /* Generalized strides for row-major addressing of data */
    int S3 = N2*N1, S2 = N1, S1 = 1;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 5.0 * log((double) N3*N2*N1) / log(2.0)
             * (sizeof(REAL) == sizeof(float) ? FLT_EPSILON : DBL_EPSILON);
    printf(" Verify the result, errthr = %.3lg\n", errthr);

    maxerr = 0.0;
    for (n3 = 0; n3 < N3; n3++) {
        for (n2 = 0; n2 < N2; n2++) {
            for (n1 = 0; n1 < N1; n1++) {
                double re_exp = 0.0, im_exp = 0.0, re_got, im_got;

                if ((n1-H1)%N1==0 && (n2-H2)%N2==0 && (n3-H3)%N3==0) {
                    re_exp = 1.0;
                }

                index = n3*S3 + n2*S2 + n1*S1;
                re_got = data[index].real;
                im_got = data[index].imag;
                err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
                if (err > maxerr) maxerr = err;
                if (!(err <= errthr)) {
                    printf(" data[%i][%i][%i]: ", n3, n2, n1);
                    printf(" expected (%.17lg,%.17lg), ", re_exp, im_exp);
                    printf(" got (%.17lg,%.17lg), ", re_got, im_got);
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
