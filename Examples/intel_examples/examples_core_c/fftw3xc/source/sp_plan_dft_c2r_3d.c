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
!       Example of using fftwf_plan_dft_c2r_3d function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init_c(fftwf_complex *x, int *N, int *H);
static int verify_r(float *x, int *N, int *H);

int main(void)
{
    /* Sizes of 3D transform */
    int N[3] = {5, 21, 29};

    /* Arbitrary harmonic used to verify FFT */
    int H[3] = {-3, -11, 13};

    /* FFTW plan handle */
    fftwf_plan c2r = 0;

    /* Pointers to input and output data */
    fftwf_complex *x = 0;
    float *y = 0;

    /* Execution status */
    int status = 0;


    printf("Example sp_plan_dft_c2r_3d\n");
    printf("3D complex-to-real out-of-place transform\n");
    printf("Configuration parameters:\n");
    printf(" N = {%d, %d, %d}\n", N[0], N[1], N[2]);
    printf(" H = {%d, %d, %d}\n", H[0], H[1], H[2]);

    printf("Allocate complex data array x(%i)\n", N[0]*N[1]*(N[2]/2+1));
    x  = fftwf_malloc(sizeof(fftwf_complex)*N[0]*N[1]*(N[2]/2+1));
    printf("Allocate real data array y(%i)\n", N[0]*N[1]*N[2]);
    y  = fftwf_malloc(sizeof(float)*N[0]*N[1]*N[2]);
    if (0 == x || 0 == y) goto failed;

    printf("Create FFTW plan for 3D c2r FFT\n");
    c2r = fftwf_plan_dft_c2r_3d(N[0], N[1], N[2], x, y, FFTW_ESTIMATE);
    if (0 == c2r) goto failed;

    printf("Initialize input for c2r transform\n");
    init_c(x, N, H);

    printf("Compute c2r FFT\n");
    fftwf_execute(c2r);

    printf("Verify the result of c2r FFT\n");
    status = verify_r(y, N, H);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plan\n");
    fftwf_destroy_plan(c2r);

    printf("Free data arrays\n");
    fftwf_free(x);
    fftwf_free(y);

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

/* Initialize array x(N) to produce unit peak at x(H) */
static void init_c(fftwf_complex *x, int *N, int *H)
{
    float TWOPI = 6.2831853071795864769f, phase;
    int n1, n2, n3, S1, S2, S3, index;

    /* Generalized strides for row-major addressing of x */
    S3 = 1;
    S2 = N[2]/2+1;
    S1 = N[1]*S2;

    for (n1 = 0; n1 < N[0]; n1++)
    {
        for (n2 = 0; n2 < N[1]; n2++)
        {
            for (n3 = 0; n3 < N[2]/2+1; n3++)
            {
                phase  = moda(n1,H[0],N[0]) / N[0];
                phase += moda(n2,H[1],N[1]) / N[1];
                phase += moda(n3,H[2],N[2]) / N[2];
                index = n1*S1 + n2*S2 + n3*S3;
                x[index][0] =  cosf( TWOPI * phase ) / (N[0]*N[1]*N[2]);
                x[index][1] = -sinf( TWOPI * phase ) / (N[0]*N[1]*N[2]);
            }
        }
    }
}

/* Verify that x has unit peak at H */
static int verify_r(float *x, int *N, int *H)
{
    float err, errthr, maxerr;
    int n1, n2, n3, S1, S2, S3, index;

    /* Generalized strides for row-major addressing of x */
    S3 = 1;
    S2 = N[2];
    S1 = N[1]*S2;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5f * logf( (float)N[0]*N[1]*N[2] ) / logf(2.0f) * FLT_EPSILON;
    printf(" Check if err is below errthr %.3g\n", errthr);

    maxerr = 0;
    for (n1 = 0; n1 < N[0]; n1++)
    {
        for (n2 = 0; n2 < N[1]; n2++)
        {
            for (n3 = 0; n3 < N[2]; n3++)
            {
                float re_exp = 0.0, re_got;

                if ((n1-H[0])%N[0]==0 &&
                    (n2-H[1])%N[1]==0 &&
                    (n3-H[2])%N[2]==0)
                {
                    re_exp = 1;
                }

                index = n1*S1 + n2*S2 + n3*S3;
                re_got = x[index];
                err  = fabsf(re_got - re_exp);
                if (err > maxerr) maxerr = err;
                if (!(err < errthr))
                {
                    printf(" x[%i][%i][%i]: ",n1,n2,n3);
                    printf(" expected %.7g, ",re_exp);
                    printf(" got %.7g, ",re_got);
                    printf(" err %.3g\n", err);
                    printf(" Verification FAILED\n");
                    return 1;
                }
            }
        }
    }
    printf(" Verified,  maximum error was %.3g\n", maxerr);
    return 0;
}
