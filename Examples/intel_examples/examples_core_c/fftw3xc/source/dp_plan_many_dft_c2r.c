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
!       Example of using fftw_plan_many_dft_c2r function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init_c(fftw_complex *x,
                   int *N, int M, int *EN, int stride, int dist, int *H);
static int verify_r(double *x,
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
    int N[3] = { 12, 24, 6 };
    int M = 16;

    /* Sizes of embedding real array, stride and distance, to be defined */
    int EN[3], EM, stride, dist;

    /* Stride and distance for the multiple transform, to be defined */

    /* Arbitrary harmonic used to verify FFT */
    int H[3] = { 2, 3, 5 };

    /* FFTW plan handle */
    fftw_plan c2r = 0;

    /* Pointer to input/output data */
    fftw_complex *x = 0;

    /* Execution status */
    int status = 0;

    printf("Example dp_plan_many_dft_c2r\n");
    printf("Multiple 3D c2r in-place FFTs\n");
    printf("Configuration parameters:\n");
    printf(" N  = {%d, %d, %d}\n", N[0], N[1], N[2]);
    printf(" M  = %d\n", M);
    printf(" H  = {%d, %d, %d}\n", H[0], H[1], H[2]);

    printf("Define data layout for complex domain\n");
    /*
     * Leading dimension is N[2], embedding array
     * has dimensions (EM,EN[0],EN[1],EN[2]).
     */
    stride = 1;
    EN[2]  = (N[2]/2+1)*stride  +1;   /* +1 is arbitrary padding */
    EN[1]  =  N[1]              +2;   /* +2 is arbitrary padding */
    EN[0]  =  N[0]              +3;   /* +3 is arbitrary padding */
    dist   =  EN[0]*EN[1]*EN[2] +4;   /* +4 is arbitrary padding */
    EM     = M;
    printf(" EM=%i, dist=%i, EN={%i,%i,%i}, stride=%i\n",
           EM,dist,EN[0],EN[1],EN[2],stride);

    printf("Allocate complex x(%i)\n", EM*dist );
    x  = fftw_malloc(sizeof(fftw_complex) * EM*dist );
    if (0 == x) goto failed;

    printf("Create FFTW plan for c2r transform\n");
    {
        int ENr[3];
        ENr[0] = EN[0];
        ENr[1] = EN[1];
        ENr[2] = EN[2] * 2;
        c2r = fftw_plan_many_dft_c2r(3, N, M,
                                     x, EN, stride, dist,
                                     (double*)x, ENr, stride, dist*2,
                                     FFTW_ESTIMATE);
        if (0 == c2r) goto failed;
    }

    printf("Initialize input for c2r transform\n");
    init_c(x, N, M, EN, stride, dist, H);

    printf("Compute c2r FFT\n");
    fftw_execute(c2r);

    printf("Verify the result of c2r FFT\n");
    status = verify_r((double*)x, N, M, EN, stride, dist*2, H);
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

/* Initialize arrays x with harmonic H */
static void init_c(fftw_complex *x,
                   int *N, int M, int *EN, int stride, int dist, int *H)
{
    double TWOPI = 6.2831853071795864769, phase;
    int n1, n2, n3, m, S1, S2, S3, index;

    S3 = stride;
    S2 = EN[2];
    S1 = EN[1]*S2;

    for (m = 0; m < M; m++)
    {
        for (n1 = 0; n1 < N[0]; n1++)
        {
            for (n2 = 0; n2 < N[1]; n2++)
            {
                for (n3 = 0; n3 < N[2]/2+1; n3++)
                {
                    phase  = moda(n1,H[0],N[0]) / N[0];
                    phase += moda(n2,H[1],N[1]) / N[1];
                    phase += moda(n3,H[2],N[2]) / N[2];
                    index = n1*S1 + n2*S2 + n3*S3 + m*dist;
                    x[index][0] =  cos( TWOPI * phase ) / (N[0]*N[1]*N[2]);
                    x[index][1] = -sin( TWOPI * phase ) / (N[0]*N[1]*N[2]);
                }
            }
        }
    }
}

/* Verify that x has unit peak at H */
static int verify_r(double *x,
                  int *N, int M, int *EN, int stride, int dist, int *H)
{
    double err, errthr, maxerr;
    int n1, n2, n3, m, S1, S2, S3, index;

    S3 = stride;
    S2 = EN[2]*2;
    S1 = EN[1]*S2;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5 * log( (double)N[0]*N[1]*N[2] ) / log(2.0) * DBL_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    maxerr = 0;
    for (m = 0; m < M; m++)
    {
        for (n1 = 0; n1 < N[0]; n1++)
        {
            for (n2 = 0; n2 < N[1]; n2++)
            {
                for (n3 = 0; n3 < N[2]; n3++)
                {
                    double re_exp = 0.0, re_got;

                    if ((n1-H[0])%N[0]==0 &&
                        (n2-H[1])%N[1]==0 &&
                        (n3-H[2])%N[2]==0)
                    {
                        re_exp = 1;
                    }

                    index = n1*S1 + n2*S2 + n3*S3 + m*dist;
                    re_got = x[index];
                    err  = fabs(re_got - re_exp);
                    if (err > maxerr) maxerr = err;
                    if (!(err < errthr))
                    {
                        printf(" x(n1=%i,n2=%i,n3=%i,m=%i): ",n1,n2,n3,m);
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
