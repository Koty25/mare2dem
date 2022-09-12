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
!       Example of using fftw_plan_dft_c2r_2d function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init_c(fftw_complex *x, int N1, int N2, int H1, int H2);
static int verify_r(double *x, int N1, int N2, int H1, int H2);

int main(void)
{
    /* Sizes of 2D transform */
    int N1 = 20;
    int N2 = 25;

    /* Arbitrary harmonic used to verify FFT */
    int H1 = N1/2;
    int H2 = -2;

    /* FFTW plan handles */
    fftw_plan c2r = 0;

    /* Pointers to input and output data */
    fftw_complex *x = 0;
    double *y = 0;

    /* Execution status */
    int status = 0;

    printf("Example dp_plan_dft_c2r_2d\n");
    printf("2D real-to-complex out-of-place transform\n");
    printf("Configuration parameters:\n");
    printf(" N = {%d, %d}\n", N1, N2);
    printf(" H = {%d, %d}\n", H1, H2);

    printf("Allocate data arrays: complex x(%i), real y(%i)\n", N1*(N2/2+1), N1*N2);
    x  = fftw_malloc(sizeof(fftw_complex)*N1*(N2/2+1));
    y  = fftw_malloc(sizeof(double)*N1*N2);
    if (0 == x || 0 == y) goto failed;

    printf("Create FFTW plan for 2D c2r out-of-place FFT\n");
    c2r = fftw_plan_dft_c2r_2d(N1, N2, x, y, FFTW_ESTIMATE);
    if (0 == c2r) goto failed;

    printf("Initialize input for c2r transform\n");
    init_c(x, N1, N2, H1, H2);

    printf("Compute c2r FFT\n");
    fftw_execute(c2r);

    printf("Verify the result of c2r FFT\n");
    status = verify_r(y, N1, N2, H1, H2);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plan\n");
    fftw_destroy_plan(c2r);

    printf("Free data array\n");
    fftw_free(x);
    fftw_free(y);

    printf("TEST %s\n",0==status ? "PASSED" : "FAILED");
    return status;

 failed:
    printf(" ERROR\n");
    status = 1;
    goto cleanup;
}

/* Compute (K*L)%M accurately */
static double moda(int K, int L, int M)
{
    return (double)(((long long)K * L) % M);
}

/* Initialize array x(N) to produce unit peak at x(H) */
static void init_c(fftw_complex *x, int N1, int N2, int H1, int H2)
{
    double TWOPI = 6.2831853071795864769, phase;
    int n1, n2, S1, S2, index;

    /* Generalized strides for row-major addressing of x */
    S2 = 1;
    S1 = N2/2+1;

    for (n1 = 0; n1 < N1; n1++)
    {
        for (n2 = 0; n2 < N2/2+1; n2++)
        {
            phase  = moda(n1,H1,N1) / N1;
            phase += moda(n2,H2,N2) / N2;
            index = n1*S1 + n2*S2;
            x[index][0] =  cos( TWOPI * phase ) / (N1*N2);
            x[index][1] = -sin( TWOPI * phase ) / (N1*N2);
        }
    }
}

/* Verify that x has unit peak at H */
static int verify_r(double *x, int N1, int N2, int H1, int H2)
{
    double err, errthr, maxerr;
    int n1, n2, S1, S2, index;

    /* Generalized strides for row-major addressing of x */
    S2 = 1;
    S1 = N2;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5 * log( (double)N1*N2 ) / log(2.0) * DBL_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    maxerr = 0;
    for (n1 = 0; n1 < N1; n1++)
    {
        for (n2 = 0; n2 < N2; n2++)
        {
            double re_exp = 0.0, re_got;

            if ((n1-H1)%N1==0 && (n2-H2)%N2==0)
            {
                re_exp = 1;
            }

            index = n1*S1 + n2*S2;
            re_got = x[index];
            err  = fabs(re_got - re_exp);
            if (err > maxerr) maxerr = err;
            if (!(err < errthr))
            {
                printf(" x[%i][%i]: ",n1,n2);
                printf(" expected %.17lg, ",re_exp);
                printf(" got %.17lg, ",re_got);
                printf(" err %.3lg\n", err);
                printf(" Verification FAILED\n");
                return 1;
            }
        }
    }
    printf(" Verified,  maximum error was %.3lg\n", maxerr);
    return 0;
}
