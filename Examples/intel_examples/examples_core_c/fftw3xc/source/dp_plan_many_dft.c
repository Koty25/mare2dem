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
!       Example of using fftw_plan_many_dft function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init(fftw_complex *x,
                 int *N, int M, int *EN, int stride, int dist, int *H);
static int verify(fftw_complex *x,
                 int *N, int M, int *EN, int stride, int dist, int *H);

int main(void)
{
    /*
     * In this example we perform M in-place 3D FFT on the data
     * contained in a larger array.  The sizes of the FFT are defined
     * by array N[], the embedding array has dimensions defined by
     * array EN[] and EM, not necessarily in this order.
     */

    /* Sizes of 3D transform and the number of them */
    int N[3] = { 20, 15, 10 };
    int M = 4;

    /* Sizes of embedding array, stride and distance, to be defined */
    int EN[3], EM, stride, dist;

    /* Arbitrary harmonic used to verify FFT */
    int H[3] = { -1, -2, -3 };

    /* FFTW plan handle */
    fftw_plan plan = 0;

    /* Pointer to input/output data */
    fftw_complex *x = 0;

    /* Execution status */
    int status = 0;

    printf("Example dp_plan_many_dft\n");
    printf("Forward multiple 3D complex in-place FFT\n");
    printf("Configuration parameters:\n");
    printf(" N  = {%d, %d, %d}\n", N[0], N[1], N[2]);
    printf(" M  = %d\n", M);
    printf(" H  = {%d, %d, %d}\n", H[0], H[1], H[2]);

    printf("Define data layout for PARALLEL transforms\n");
    /*
     * Leading dimension is N[2] (parallel transforms): embedding array
     * has dimensions (EM,EN[0],EN[1],EN[2]).
     */
    stride = 1;
    EN[2]  = N[2]*stride       +1; /* +1 is arbitrary padding */
    EN[1]  = N[1]              +2; /* +2 is arbitrary padding */
    EN[0]  = N[0]              +3; /* +3 is arbitrary padding */
    dist   = EN[0]*EN[1]*EN[2] +4; /* +4 is arbitrary padding */
    EM     = M;
    printf(" EM=%i, dist=%i, EN={%i,%i,%i}, stride=%i\n",
           EM,dist,EN[0],EN[1],EN[2],stride);

    printf("Allocate x(%i)\n", EM*dist );
    x  = fftw_malloc(sizeof(fftw_complex) * EM*dist );
    if (0 == x) goto failed;

    printf("Create FFTW plan for forward transform\n");
    plan = fftw_plan_many_dft(3, N, M,
                              x, EN, stride, dist,
                              x, EN, stride, dist,
                              FFTW_FORWARD, FFTW_ESTIMATE);
    if (0 == plan) goto failed;

    printf("Initialize input for forward transform\n");
    init(x, N, M, EN, stride, dist, H);

    printf("Compute forward FFT\n");
    fftw_execute(plan);

    printf("Verify the result of forward FFT\n");
    status = verify(x, N, M, EN, stride, dist, H);
    if (0 != status) goto failed;

    printf("Destroy FFTW plan\n");
    fftw_destroy_plan(plan);

    printf("Free data array\n");
    fftw_free(x);


    printf("\nDefine data layout for VECTOR transforms\n");
    /*
     * Leading dimension is M (vector transforms): embedding array
     * has dimensions (EN[0],EN[1],EN[2],EM).
     */
    dist   = 1;
    EM     = M*dist         +1; /* +1 is arbitrary padding */
    stride = EM;                /* must have stride==EM */
    EN[2]  = N[2]           +2; /* +2 is arbitrary padding */
    EN[1]  = N[1]           +3; /* +3 is arbitrary padding */
    EN[0]  = N[0];

    printf(" EN={%i,%i,%i}, stride=%i, EM=%i, dist=%i\n",
           EN[0],EN[1],EN[2],stride,EM,dist);

    printf("Allocate x(%i)\n", EN[0]*EN[1]*EN[2]*stride );
    x  = fftw_malloc(sizeof(fftw_complex) * EN[0]*EN[1]*EN[2]*stride );
    if (0 == x) goto failed;

    printf("Create FFTW plan for forward transform\n");
    plan = fftw_plan_many_dft(3, N, M,
                              x, EN, stride, dist,
                              x, EN, stride, dist,
                              FFTW_FORWARD, FFTW_ESTIMATE);
    if (0 == plan) goto failed;

    printf("Initialize input for forward transform\n");
    init(x, N, M, EN, EM, dist, H);

    printf("Compute forward FFT\n");
    fftw_execute(plan);

    printf("Verify the result of forward FFT\n");
    status = verify(x, N, M, EN, EM, dist, H);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plan\n");
    fftw_destroy_plan(plan);

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

/* Initialize arrays x with harmonic H */
static void init(fftw_complex *x,
                 int *N, int M, int *EN, int stride, int dist, int *H)
{
    double TWOPI = 6.2831853071795864769, phase;
    int n1, n2, n3, m, N1, N2, N3, S1, S2, S3, SM, H1, H2, H3, index;

    SM = dist;
    S3 = stride;   N3 = N[2]; H3 = H[2];
    S2 = EN[2]*S3; N2 = N[1]; H2 = H[1];
    S1 = EN[1]*S2; N1 = N[0]; H1 = H[0];

    for (m = 0; m < M; m++)
    {
        for (n1 = 0; n1 < N1; n1++)
        {
            for (n2 = 0; n2 < N2; n2++)
            {
                for (n3 = 0; n3 < N3; n3++)
                {
                    phase  = moda(n1,H1,N1) / N1;
                    phase += moda(n2,H2,N2) / N2;
                    phase += moda(n3,H3,N3) / N3;
                    index = n1*S1 + n2*S2 + n3*S3 + m*SM;
                    x[index][0] = cos( TWOPI * phase ) / (N1*N2*N3);
                    x[index][1] = sin( TWOPI * phase ) / (N1*N2*N3);
                }
            }
        }
    }
}

/* Verify that x has unit peak at H */
static int verify(fftw_complex *x,
                  int *N, int M, int *EN, int stride, int dist, int *H)
{
    double err, errthr, maxerr;
    int n1, n2, n3, m, N1, N2, N3, S1, S2, S3, SM, H1, H2, H3, index;

    SM = dist;
    S3 = stride;   N3 = N[2]; H3 = H[2];
    S2 = EN[2]*S3; N2 = N[1]; H2 = H[1];
    S1 = EN[1]*S2; N1 = N[0]; H1 = H[0];

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 5.0 * log( (double)N1*N2*N3 ) / log(2.0) * DBL_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    maxerr = 0;
    for (m = 0; m < M; m++)
    {
        for (n1 = 0; n1 < N1; n1++)
        {
            for (n2 = 0; n2 < N2; n2++)
            {
                for (n3 = 0; n3 < N3; n3++)
                {
                    double re_exp = 0.0, im_exp = 0.0, re_got, im_got;

                    if ((n1-H1)%N1==0 && (n2-H2)%N2==0 && (n3-H3)%N3==0)
                    {
                        re_exp = 1;
                    }

                    index = n1*S1 + n2*S2 + n3*S3 + m*SM;
                    re_got = x[index][0];
                    im_got = x[index][1];
                    err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
                    if (err > maxerr) maxerr = err;
                    if (!(err < errthr))
                    {
                        printf(" x(n1=%i,n2=%i,n3=%i,m=%i): ",n1,n2,n3,m);
                        printf(" expected (%.17lg,%.17lg), ",re_exp,im_exp);
                        printf(" got (%.17lg,%.17lg), ",re_got,im_got);
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
