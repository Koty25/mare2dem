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
!       Example of using fftwf_plan_many_dft_r2c function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init_r(float *x,
                   int *N, int M, int *EN, int stride, int dist, int *H);
static int verify_c(fftwf_complex *x,
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
    int N[3] = { 11, 22, 33 };
    int M = 8;

    /* Sizes of embedding real array, stride and distance, to be defined */
    int EN[3], EM, stride, dist;

    /* Stride and distance for the multiple transform, to be defined */

    /* Arbitrary harmonic used to verify FFT */
    int H[3] = { 2, 3, 5 };

    /* FFTW plan handle */
    fftwf_plan r2c = 0;

    /* Pointer to input/output data */
    float *x = 0;

    /* Execution status */
    int status = 0;

    printf("Example sp_plan_many_dft_r2c\n");
    printf("Multiple 3D r2c in-place FFTs\n");
    printf("Configuration parameters:\n");
    printf(" N  = {%d, %d, %d}\n", N[0], N[1], N[2]);
    printf(" M  = %d\n", M);
    printf(" H  = {%d, %d, %d}\n", H[0], H[1], H[2]);

    printf("Define data layout for real domain\n");
    /*
     * Leading dimension is N[2], embedding array
     * has dimensions (EM,EN[0],EN[1],EN[2]).
     */
    stride = 1;
    EN[2]  = 2*( (N[2]/2+1)*stride +1 ); /* +1 is arbitrary padding */
    EN[1]  =     N[1]              +2;   /* +2 is arbitrary padding */
    EN[0]  =     N[0]              +3;   /* +3 is arbitrary padding */
    dist   =     EN[0]*EN[1]*EN[2] +4;   /* +4 is arbitrary padding */
    EM     = M;
    printf(" EM=%i, dist=%i, EN={%i,%i,%i}, stride=%i\n",
           EM,dist,EN[0],EN[1],EN[2],stride);

    printf("Allocate real x(%i)\n", EM*dist );
    x  = fftwf_malloc(sizeof(float) * EM*dist );
    if (0 == x) goto failed;

    printf("Create FFTW plan for r2c transform\n");
    {
        int ENc[3];
        ENc[0] = EN[0];
        ENc[1] = EN[1];
        ENc[2] = EN[2] / 2;
        r2c = fftwf_plan_many_dft_r2c(3, N, M,
                                     x, EN, stride, dist,
                                     (fftwf_complex*)x, ENc, stride, dist/2,
                                     FFTW_ESTIMATE);
        if (0 == r2c) goto failed;
    }

    printf("Initialize input for r2c transform\n");
    init_r(x, N, M, EN, stride, dist, H);

    printf("Compute r2c FFT\n");
    fftwf_execute(r2c);

    printf("Verify the result of r2c FFT\n");
    status = verify_c((fftwf_complex*)x, N, M, EN, stride, dist/2, H);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plan\n");
    fftwf_destroy_plan(r2c);

    printf("Free data array\n");
    fftwf_free(x);

    printf("TEST %s\n",0==status ? "PASSED" : "FAILED");
    return status;

 failed:
    printf(" ERROR\n");
    status = 1;
    goto cleanup;
}

/* Compute (K*L)%M accurately */
static float moda(int K, int L, int M)
{
    return (float)(((long long)K * L) % M);
}

/* Initialize arrays x with harmonic H */
static void init_r(float *x,
                   int *N, int M, int *EN, int stride, int dist, int *H)
{
    float TWOPI = 6.2831853071795864769f, phase, factor;
    int n1, n2, n3, m, S1, S2, S3, index;

    S3 = stride;
    S2 = EN[2];
    S1 = EN[1]*S2;
    if ((N[0]-H[0])%N[0]==0 && (N[1]-H[1])%N[1]==0 && (N[2]-H[2])%N[2]==0)
    {
        factor = 1;
    }
    else
    {
        factor = 2;
    }

    for (m = 0; m < M; m++)
    {
        for (n1 = 0; n1 < N[0]; n1++)
        {
            for (n2 = 0; n2 < N[1]; n2++)
            {
                for (n3 = 0; n3 < N[2]; n3++)
                {
                    phase  = moda(n1,H[0],N[0]) / N[0];
                    phase += moda(n2,H[1],N[1]) / N[1];
                    phase += moda(n3,H[2],N[2]) / N[2];
                    index = n1*S1 + n2*S2 + n3*S3 + m*dist;
                    x[index] = factor * cosf( TWOPI * phase ) / (N[0]*N[1]*N[2]);
                }
            }
        }
    }
}

/* Verify that x has unit peak at H */
static int verify_c(fftwf_complex *x,
                  int *N, int M, int *EN, int stride, int dist, int *H)
{
    float err, errthr, maxerr;
    int n1, n2, n3, m, S1, S2, S3, index;

    S3 = stride;
    S2 = EN[2]/2;
    S1 = EN[1]*S2;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5f * logf( (float)N[0]*N[1]*N[2] ) / logf(2.0f) * FLT_EPSILON;
    printf(" Check if err is below errthr %.3g\n", errthr);

    maxerr = 0;
    for (m = 0; m < M; m++)
    {
        for (n1 = 0; n1 < N[0]; n1++)
        {
            for (n2 = 0; n2 < N[1]; n2++)
            {
                for (n3 = 0; n3 < N[2]/2+1; n3++)
                {
                    float re_exp = 0.0, im_exp = 0.0, re_got, im_got;

                    if ((n1-H[0])%N[0]==0 &&
                        (n2-H[1])%N[1]==0 &&
                        (n3-H[2])%N[2]==0)
                    {
                        re_exp = 1;
                    }
                    else if ((-n1-H[0])%N[0]==0 &&
                             (-n2-H[1])%N[1]==0 &&
                             (-n3-H[2])%N[2]==0)
                    {
                        re_exp = 1;
                    }

                    index = n1*S1 + n2*S2 + n3*S3 + m*dist;
                    re_got = x[index][0];
                    im_got = x[index][1];
                    err  = fabsf(re_got - re_exp) + fabsf(im_got - im_exp);
                    if (err > maxerr) maxerr = err;
                    if (!(err < errthr))
                    {
                        printf(" x(n1=%i,n2=%i,n3=%i,m=%i): ",n1,n2,n3,m);
                        printf(" expected (%.7g,%.7g), ",re_exp,im_exp);
                        printf(" got (%.7g,%.7g), ",re_got,im_got);
                        printf(" err %.3g\n", err);
                        printf(" Verification FAILED\n");
                        return 1;
                    }
                }
            }
        }
    }
    printf(" Verified,  maximum error was %.3g\n", maxerr);
    return 0;
}
