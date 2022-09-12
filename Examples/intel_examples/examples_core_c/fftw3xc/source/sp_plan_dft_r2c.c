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
!       Example of using fftwf_plan_dft_r2c function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init_r(float *x, int *N, int *H);
static int verify_c(fftwf_complex *x, int *N, int *H);

int main(void)
{
    /* Sizes of 4D transform */
    int N[4] = {8, 16, 16, 4};

    /* Arbitrary harmonic used to verify FFT */
    int H[4] = {-1, -2, -3, -4};

    /* FFTW plan handle */
    fftwf_plan r2c = 0;

    /* Pointer to input/output data */
    float *x = 0;

    /* Execution status */
    int status = 0;


    printf("Example sp_plan_dft_r2c\n");
    printf("4D real-to-complex inplace transform\n");
    printf("Configuration parameters:\n");
    printf(" N = {%d, %d, %d, %d}\n", N[0], N[1], N[2], N[3]);
    printf(" H = {%d, %d, %d, %d}\n", H[0], H[1], H[2], H[3]);

    printf("Allocate real data array x(%i)\n",N[0]*N[1]*N[2]*(N[3]/2+1)*2);
    x  = fftwf_malloc(sizeof(float)*N[0]*N[1]*N[2]*(N[3]/2+1)*2);
    if (0 == x) goto failed;

    printf("Create FFTW plan for 4D r2c inplace FFT\n");
    r2c = fftwf_plan_dft_r2c(4, N, x, (fftwf_complex*)x, FFTW_ESTIMATE);
    if (0 == r2c) goto failed;

    printf("Initialize input for r2c transform\n");
    init_r(x, N, H);

    printf("Compute r2c FFT\n");
    fftwf_execute(r2c);

    printf("Verify the result of r2c FFT\n");
    status = verify_c((fftwf_complex*)x, N, H);
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

/* Initialize array x(N) to produce unit peaks at x(H) and x(N-H) */
static void init_r(float *x, int *N, int *H)
{
    float TWOPI = 6.2831853071795864769f, phase, factor;
    int n1, n2, n3, n4, S1, S2, S3, S4, index;

    /* Generalized strides for row-major addressing of x */
    S4 = 1;
    S3 = (N[3]/2+1)*2;
    S2 = N[2]*S3;
    S1 = N[1]*S2;

    if ((N[0]-H[0]%N[0])==0 &&
        (N[1]-H[1]%N[1])==0 &&
        (N[2]-H[2]%N[2])==0 &&
        (N[3]-H[3]%N[3])==0)
    {
        factor = 1;
    }
    else
    {
        factor = 2;
    }

    for (n1 = 0; n1 < N[0]; n1++)
    {
        for (n2 = 0; n2 < N[1]; n2++)
        {
            for (n3 = 0; n3 < N[2]; n3++)
            {
                for (n4 = 0; n4 < N[3]; n4++)
                {
                    phase  = moda(n1,H[0],N[0]) / N[0];
                    phase += moda(n2,H[1],N[1]) / N[1];
                    phase += moda(n3,H[2],N[2]) / N[2];
                    phase += moda(n4,H[3],N[3]) / N[3];
                    index = n1*S1 + n2*S2 + n3*S3 + n4*S4;
                    x[index] = factor * cosf( TWOPI * phase ) / (N[0]*N[1]*N[2]*N[3]);
                }
            }
        }
    }
}

/* Verify that x has unit peak at H */
static int verify_c(fftwf_complex *x, int *N, int *H)
{
    float err, errthr, maxerr;
    int n1, n2, n3, n4, S1, S2, S3, S4, index;

    /* Generalized strides for row-major addressing of x */
    S4 = 1;
    S3 = N[3]/2+1;
    S2 = N[2]*S3;
    S1 = N[1]*S2;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5f * logf( (float)N[0]*N[1]*N[2]*N[3] ) / logf(2.0f) * FLT_EPSILON;
    printf(" Check if err is below errthr %.3g\n", errthr);

    maxerr = 0;
    for (n1 = 0; n1 < N[0]; n1++)
    {
        for (n2 = 0; n2 < N[1]; n2++)
        {
            for (n3 = 0; n3 < N[2]; n3++)
            {
                for (n4 = 0; n4 < N[3]/2+1; n4++)
                {
                    float re_exp = 0.0, im_exp = 0.0, re_got, im_got;

                    if ((n1-H[0])%N[0]==0 &&
                        (n2-H[1])%N[1]==0 &&
                        (n3-H[2])%N[2]==0 &&
                        (n4-H[3])%N[3]==0)
                    {
                        re_exp = 1;
                    }
                    else if ((-n1-H[0])%N[0]==0 &&
                             (-n2-H[1])%N[1]==0 &&
                             (-n3-H[2])%N[2]==0 &&
                             (-n4-H[3])%N[3]==0)
                    {
                        re_exp = 1;
                    }

                    index = n1*S1 + n2*S2 + n3*S3 + n4*S4;
                    re_got = x[index][0];
                    im_got = x[index][1];
                    err  = fabsf(re_got - re_exp) + fabsf(im_got - im_exp);
                    if (err > maxerr) maxerr = err;
                    if (!(err < errthr))
                    {
                        printf(" x[%i][%i][%i][%i]: ",n1,n2,n3,n4);
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
