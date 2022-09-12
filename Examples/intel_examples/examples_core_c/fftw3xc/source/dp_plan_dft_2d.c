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
!       Example of using fftw_plan_dft_2d function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init(fftw_complex *x, int N1,int N2, int H1,int H2);
static int verify(fftw_complex *x, int N1,int N2, int H1,int H2);

int main(void)
{
    /* Sizes of 2D transform */
    int N[2] = {32, 32};

    /* Arbitrary harmonic used to verify FFT */
    int H[2] = {-1, -2};

    /* FFTW plan handles */
    fftw_plan forward_plan = 0, backward_plan = 0;

    /* Pointer to input/output data */
    fftw_complex *x = 0;

    /* Execution status */
    int status = 0;

    printf("Example dp_plan_dft_2d\n");
    printf("Forward and backward 2D complex inplace transform\n");
    printf("Configuration parameters:\n");
    printf(" N = {%d, %d}\n", N[0], N[1]);
    printf(" H = {%d, %d}\n", H[0], H[1]);

    printf("Allocate array for input data\n");
    x  = fftw_malloc(sizeof(fftw_complex)*N[0]*N[1]);
    if (0 == x) goto failed;

    printf("Create FFTW plan for 2D double-precision forward transform\n");
    forward_plan = fftw_plan_dft(2, N, x, x, FFTW_FORWARD, FFTW_ESTIMATE);
    if (0 == forward_plan) goto failed;

    printf("Create FFTW plan for 2D double-precision backward transform\n");
    backward_plan = fftw_plan_dft(2, N, x, x, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (0 == backward_plan) goto failed;

    printf("Initialize input for forward transform\n");
    init(x, N[0],N[1], H[0],H[1]);

    printf("Compute forward FFT\n");
    fftw_execute(forward_plan);

    printf("Verify the result of forward FFT\n");
    status = verify(x, N[0],N[1], H[0],H[1]);
    if (0 != status) goto failed;

    printf("Initialize input for backward transform\n");
    init(x, N[0],N[1], -H[0],-H[1]);

    printf("Compute backward transform\n");
    fftw_execute(backward_plan);

    printf("Verify the result of backward FFT\n");
    status = verify(x, N[0],N[1], H[0],H[1]);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plans\n");
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);

    printf("Free data array\n");
    fftw_free(x);

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

/* Initialize array x[N1][N2] with harmonic H */
static void init(fftw_complex *x, int N1,int N2, int H1,int H2)
{
    double TWOPI = 6.2831853071795864769, phase;
    int n1, n2, index;

    /* Generalized strides for row-major addressing of x */
    int S1 = N2, S2 = 1;

    for (n1 = 0; n1 < N1; n1++)
    {
        for (n2 = 0; n2 < N2; n2++)
        {
            phase  = moda(n1,H1,N1) / N1;
            phase += moda(n2,H2,N2) / N2;
            index = n1*S1 + n2*S2;
            x[index][0] = cos( TWOPI * phase ) / (N1*N2);
            x[index][1] = sin( TWOPI * phase ) / (N1*N2);
        }
    }
}

/* Verify that x has unit peak at H */
static int verify(fftw_complex *x, int N1,int N2, int H1,int H2)
{
    double err, errthr, maxerr;
    int n1, n2, index;

    /* Generalized strides for row-major addressing of x */
    int S1 = N2, S2 = 1;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 5.0 * log( (double)N1*N2 ) / log(2.0) * DBL_EPSILON;
    printf(" Verify the result, errthr = %.3lg\n", errthr);

    maxerr = 0;
    for (n1 = 0; n1 < N1; n1++)
    {
        for (n2 = 0; n2 < N2; n2++)
        {
            double re_exp = 0.0, im_exp = 0.0, re_got, im_got;

            if ((n1-H1)%N1==0 && (n2-H2)%N2==0)
            {
                re_exp = 1;
            }

            index = n1*S1 + n2*S2;
            re_got = x[index][0];
            im_got = x[index][1];
            err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
            if (err > maxerr) maxerr = err;
            if (!(err < errthr))
            {
                printf(" x[%i][%i]: ",n1,n2);
                printf(" expected (%.17lg,%.17lg), ",re_exp,im_exp);
                printf(" got (%.17lg,%.17lg), ",re_got,im_got);
                printf(" err %.3lg\n", err);
                printf(" Verification FAILED\n");
                return 1;
            }
        }
    }
    printf(" Verified, maximum error was %.3lg\n", maxerr);
    return 0;
}
