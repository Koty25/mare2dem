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
!       Example of using fftw_plan_dft_c2r function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init_c(fftw_complex *x, int *N, int *H);
static int verify_r(double *x, int *N, int *H);

int main(void)
{
    /* Sizes of 4D transform */
    int N[4] = {18, 6, 6, 4};

    /* Arbitrary harmonic used to verify FFT */
    int H[4] = {-2, -3, -4, -5};

    /* FFTW plan handle */
    fftw_plan c2r = 0;

    /* Pointer to input/output data */
    fftw_complex *x = 0;

    /* Execution status */
    int status = 0;


    printf("Example dp_plan_dft_c2r\n");
    printf("4D complex-to-real in-place transform\n");
    printf("Configuration parameters:\n");
    printf(" N = {%d, %d, %d, %d}\n", N[0], N[1], N[2], N[3]);
    printf(" H = {%d, %d, %d, %d}\n", H[0], H[1], H[2], H[3]);

    printf("Allocate complex data array x(%i)\n",N[0]*N[1]*N[2]*(N[3]/2+1));
    x  = fftw_malloc(sizeof(fftw_complex)*N[0]*N[1]*N[2]*(N[3]/2+1));
    if (0 == x) goto failed;

    printf("Create FFTW plan for 4D c2r inplace FFT\n");
    c2r = fftw_plan_dft_c2r(4, N, x, (double*)x, FFTW_ESTIMATE);
    if (0 == c2r) goto failed;

    printf("Initialize input for c2r transform\n");
    init_c(x, N, H);

    printf("Compute c2r FFT\n");
    fftw_execute(c2r);

    printf("Verify the result of c2r FFT\n");
    status = verify_r((double*)x, N, H);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plan\n");
    fftw_destroy_plan(c2r);

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

/* Initialize array x(N) to produce unit peaks at x(H) */
static void init_c(fftw_complex *x, int *N, int *H)
{
    double TWOPI = 6.2831853071795864769, phase;
    int n1, n2, n3, n4, S1, S2, S3, S4, index;

    /* Generalized strides for row-major addressing of x */
    S4 = 1;
    S3 = (N[3]/2+1);
    S2 = N[2]*S3;
    S1 = N[1]*S2;

    for (n1 = 0; n1 < N[0]; n1++)
    {
        for (n2 = 0; n2 < N[1]; n2++)
        {
            for (n3 = 0; n3 < N[2]; n3++)
            {
                for (n4 = 0; n4 < N[3]/2+1; n4++)
                {
                    phase  = moda(n1,H[0],N[0]) / N[0];
                    phase += moda(n2,H[1],N[1]) / N[1];
                    phase += moda(n3,H[2],N[2]) / N[2];
                    phase += moda(n4,H[3],N[3]) / N[3];
                    index = n1*S1 + n2*S2 + n3*S3 + n4*S4;
                    x[index][0] =  cos( TWOPI * phase ) / (N[0]*N[1]*N[2]*N[3]);
                    x[index][1] = -sin( TWOPI * phase ) / (N[0]*N[1]*N[2]*N[3]);
                }
            }
        }
    }
}

/* Verify that x has unit peak at H */
static int verify_r(double *x, int *N, int *H)
{
    double err, errthr, maxerr;
    int n1, n2, n3, n4, S1, S2, S3, S4, index;

    /* Generalized strides for row-major addressing of x */
    S4 = 1;
    S3 = 2*(N[3]/2+1);
    S2 = N[2]*S3;
    S1 = N[1]*S2;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5 * log( (double)N[0]*N[1]*N[2]*N[3] ) / log(2.0) * DBL_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    maxerr = 0;
    for (n1 = 0; n1 < N[0]; n1++)
    {
        for (n2 = 0; n2 < N[1]; n2++)
        {
            for (n3 = 0; n3 < N[2]; n3++)
            {
                for (n4 = 0; n4 < N[3]; n4++)
                {
                    double re_exp = 0.0, re_got;

                    if ((n1-H[0])%N[0]==0 &&
                        (n2-H[1])%N[1]==0 &&
                        (n3-H[2])%N[2]==0 &&
                        (n4-H[3])%N[3]==0)
                    {
                        re_exp = 1;
                    }

                    index = n1*S1 + n2*S2 + n3*S3 + n4*S4;
                    re_got = x[index];
                    err  = fabs(re_got - re_exp);
                    if (err > maxerr) maxerr = err;
                    if (!(err < errthr))
                    {
                        printf(" x[%i][%i][%i][%i]: ",n1,n2,n3,n4);
                        printf(" expected %.17lg, ",re_exp);
                        printf(" got %.17lg, ",re_got);
                        printf(" err %.3lg\n", err);
                        printf(" Verification FAILED\n");
                        return 1;
                    }
                }
            }
        }
    }
    printf(" Verified,  maximum error was %.3lg\n", maxerr);
    return 0;
}
