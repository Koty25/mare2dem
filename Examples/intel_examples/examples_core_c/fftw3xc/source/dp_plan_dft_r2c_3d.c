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
!       Example of using fftw_plan_dft_r2c_3d function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init_r(double *x, int *N, int *H);
static int verify_c(fftw_complex *x, int *N, int *H);

int main(void)
{
    /* Sizes of 3D transform */
    int N[3] = {7, 14, 21};

    /* Arbitrary harmonic used to verify FFT */
    int H[3] = {-2, -3, -4};

    /* FFTW plan handle */
    fftw_plan r2c = 0;

    /* Pointer to input/output data */
    double *x = 0;

    /* Execution status */
    int status = 0;


    printf("Example dp_plan_dft_r2c_3d\n");
    printf("3D real-to-complex inplace transform\n");
    printf("Configuration parameters:\n");
    printf(" N = {%d, %d, %d}\n", N[0], N[1], N[2]);
    printf(" H = {%d, %d, %d}\n", H[0], H[1], H[2]);

    printf("Allocate real data array x(%i)\n",2*N[0]*N[1]*(N[2]/2+1));
    x  = fftw_malloc(sizeof(double)*2*N[0]*N[1]*(N[2]/2+1));
    if (0 == x) goto failed;

    printf("Create FFTW plan for 3D r2c inplace FFT\n");
    r2c = fftw_plan_dft_r2c_3d(N[0], N[1], N[2], x, (fftw_complex*)x, FFTW_ESTIMATE);
    if (0 == r2c) goto failed;

    printf("Initialize input for r2c transform\n");
    init_r(x, N, H);

    printf("Compute r2c FFT\n");
    fftw_execute(r2c);

    printf("Verify the result of r2c FFT\n");
    status = verify_c((fftw_complex*)x, N, H);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plan\n");
    fftw_destroy_plan(r2c);

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

/* Initialize array x(N) to produce unit peaks at x(H) and x(N-H) */
static void init_r(double *x, int *N, int *H)
{
    double TWOPI = 6.2831853071795864769, phase, factor;
    int n1, n2, n3, S1, S2, S3, index;

    /* Generalized strides for row-major addressing of x */
    S3 = 1;
    S2 = (N[2]/2+1)*2;
    S1 = N[1]*S2;

    if ((N[0]-H[0]%N[0])==0 &&
        (N[1]-H[1]%N[1])==0 &&
        (N[2]-H[2]%N[2])==0)
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
                phase  = moda(n1,H[0],N[0]) / N[0];
                phase += moda(n2,H[1],N[1]) / N[1];
                phase += moda(n3,H[2],N[2]) / N[2];
                index = n1*S1 + n2*S2 + n3*S3;
                x[index] = factor * cos( TWOPI * phase ) / (N[0]*N[1]*N[2]);
            }
        }
    }
}

/* Verify that x has unit peak at H */
static int verify_c(fftw_complex *x, int *N, int *H)
{
    double err, errthr, maxerr;
    int n1, n2, n3, S1, S2, S3, index;

    /* Generalized strides for row-major addressing of x */
    S3 = 1;
    S2 = N[2]/2+1;
    S1 = N[1]*S2;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5 * log( (double)N[0]*N[1]*N[2] ) / log(2.0) * DBL_EPSILON;
    printf(" Check if err is below errthr %.3lg\n", errthr);

    maxerr = 0;
    for (n1 = 0; n1 < N[0]; n1++)
    {
        for (n2 = 0; n2 < N[1]; n2++)
        {
            for (n3 = 0; n3 < N[2]/2+1; n3++)
            {
                double re_exp = 0.0, im_exp = 0.0, re_got, im_got;

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

                index = n1*S1 + n2*S2 + n3*S3;
                re_got = x[index][0];
                im_got = x[index][1];
                err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
                if (err > maxerr) maxerr = err;
                if (!(err < errthr))
                {
                    printf(" x[%i][%i][%i]: ",n1,n2,n3);
                    printf(" expected (%.17lg,%.17lg), ",re_exp,im_exp);
                    printf(" got (%.17lg,%.17lg), ",re_got,im_got);
                    printf(" err %.3lg\n", err);
                    printf(" Verification FAILED\n");
                    return 1;
                }
            }
        }
    }
    printf(" Verified,  maximum error was %.3lg\n", maxerr);
    return 0;
}
