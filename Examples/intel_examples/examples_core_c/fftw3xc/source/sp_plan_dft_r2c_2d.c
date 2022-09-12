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
!       Example of using fftwf_plan_dft_r2c_2d function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init_r(float *x, int N1, int N2, int H1, int H2);
static int verify_c(fftwf_complex *x, int N1, int N2, int H1, int H2);

int main(void)
{
    /* Sizes of 2D transform */
    int N1 = 10;
    int N2 = 15;

    /* Arbitrary harmonic used to verify FFT */
    int H1 = 1;
    int H2 = N2/2;

    /* FFTW plan handles */
    fftwf_plan r2c = 0;

    /* Pointer to input/output data */
    float *x = 0;

    /* Execution status */
    int status = 0;

    printf("Example sp_plan_dft_r2c_2d\n");
    printf("2D real-to-complex inplace transform\n");
    printf("Configuration parameters:\n");
    printf(" N = {%d, %d}\n", N1, N2);
    printf(" H = {%d, %d}\n", H1, H2);

    printf("Allocate real data array x(%i)\n",2*N1*(N2/2+1));
    x  = fftwf_malloc(sizeof(float)*2*N1*(N2/2+1));
    if (0 == x) goto failed;

    printf("Create FFTW plan for 2D r2c inplace FFT\n");
    r2c = fftwf_plan_dft_r2c_2d(N1, N2, x, (fftwf_complex*)x, FFTW_ESTIMATE);
    if (0 == r2c) goto failed;

    printf("Initialize input for r2c transform\n");
    init_r(x, N1, N2, H1, H2);

    printf("Compute r2c FFT\n");
    fftwf_execute(r2c);

    printf("Verify the result of r2c FFT\n");
    status = verify_c((fftwf_complex*)x, N1, N2, H1, H2);
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
static void init_r(float *x, int N1, int N2, int H1, int H2)
{
    float TWOPI = 6.2831853071795864769f, phase, factor;
    int n1, n2, S1, S2, index;

    /* Generalized strides for row-major addressing of x */
    S2 = 1;
    S1 = (N2/2+1)*2;

    factor = ((N1-H1%N1)==0 && (N2-H2%N2)==0) ? 1.0f : 2.0f;
    for (n1 = 0; n1 < N1; n1++)
    {
        for (n2 = 0; n2 < N2; n2++)
        {
            phase  = moda(n1,H1,N1) / N1;
            phase += moda(n2,H2,N2) / N2;
            index = n1*S1 + n2*S2;
            x[index] = factor * cosf( TWOPI * phase ) / (N1*N2);
        }
    }
}

/* Verify that x has unit peak at H */
static int verify_c(fftwf_complex *x, int N1, int N2, int H1, int H2)
{
    float err, errthr, maxerr;
    int n1, n2, S1, S2, index;

    /* Generalized strides for row-major addressing of x */
    S2 = 1;
    S1 = N2/2+1;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5f * logf( (float)N1*N2 ) / logf(2.0f) * FLT_EPSILON;
    printf(" Check if err is below errthr %.3g\n", errthr);

    maxerr = 0;
    for (n1 = 0; n1 < N1; n1++)
    {
        for (n2 = 0; n2 < N2/2+1; n2++)
        {
            float re_exp = 0.0, im_exp = 0.0, re_got, im_got;

            if ((( n1-H1)%N1==0 && ( n2-H2)%N2==0) ||
                ((-n1-H1)%N1==0 && (-n2-H2)%N2==0))
            {
                re_exp = 1;
            }

            index = n1*S1 + n2*S2;
            re_got = x[index][0];
            im_got = x[index][1];
            err  = fabsf(re_got - re_exp) + fabsf(im_got - im_exp);
            if (err > maxerr) maxerr = err;
            if (!(err < errthr))
            {
                printf(" x[%i][%i]: ",n1,n2);
                printf(" expected (%.7g,%.7g), ",re_exp,im_exp);
                printf(" got (%.7g,%.7g), ",re_got,im_got);
                printf(" err %.3g\n", err);
                printf(" Verification FAILED\n");
                return 1;
            }
        }
    }
    printf(" Verified,  maximum error was %.3g\n", maxerr);
    return 0;
}
