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
!       Example of using fftwf_plan_guru_dft_c2r function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init_r(float *x, int *N, int *rs, int M, int rd, int *H);
static void init_c(fftwf_complex *x, int *N, int *rs, int M, int rd, int *H);
static int verify_r(float *x, int *N, int *rs, int M, int rd, int *H);
static int verify_c(fftwf_complex *x, int *N, int *rs, int M, int rd, int *H);

int main(void)
{
    /*
     * In this example we show how to compute multiple
     * three-dimensional in-place c2r FFTs by one call of FFTW.
     */

    /* Sizes of 3D transform and the number of them.
     * N[2] is the last dimension, it will be halved in complex domain.
     */
    int N[3] = { 15, 31, 7};
    int M = 15;

    /* Arbitrary harmonic used to verify FFT */
    int H[3] = { -9, -8, -7 };

    /* Strides and distance describe data layout for real and complex domains */
    int rstride[3], rdist;
    int cstride[3], cdist;

    /* Size of input/output array, in complex elements */
    int nc;

    /* FFTW plan handles */
    fftwf_plan r2c = 0, c2r = 0;

    /* Pointer to input/output data */
    fftwf_complex *x = 0;

    /* Execution status */
    int status = 0;

    printf("Example sp_plan_guru_dft_c2r\n");
    printf("Multiple 3D in-place r2c/c2r FFT\n");
    printf("Configuration parameters:\n");
    printf(" N  = {%d, %d, %d}\n", N[0], N[1], N[2]);
    printf(" M  = %d\n", M);
    printf(" H  = {%d, %d, %d}\n", H[0], H[1], H[2]);

    printf("Define layout of data in complex domain\n");
    /*
     * Model declaration x[ M ][ N[0] ][ N[1] ][ N[2]/2+1 ].
     */
    cstride[2] = 1;
    cstride[1] = (N[2]/2+1) * cstride[2];
    cstride[0] = N[1] * cstride[1];
    cdist      = N[0] * cstride[0];
    nc         = M * cdist;
    rstride[2] = 1;
    rstride[1] = 2*cstride[1];
    rstride[0] = 2*cstride[0];
    rdist      = 2*cdist;

    printf(" rdist=%i, rstride={%i,%i,%i}\n",
           rdist, rstride[0],rstride[1],rstride[2]);
    printf(" cdist=%i, cstride={%i,%i,%i}\n",
           cdist, cstride[0],cstride[1],cstride[2]);

    printf("Allocate complex x(%i)\n", nc );
    x  = fftwf_malloc(sizeof(fftwf_complex) * nc );
    if (0 == x) goto failed;

    printf("Create FFTW plan for c2r transform\n");
    {
        fftwf_iodim dim[3], M_dim;
        dim[2].n = N[2]; dim[2].is = cstride[2]; dim[2].os = rstride[2];
        dim[1].n = N[1]; dim[1].is = cstride[1]; dim[1].os = rstride[1];
        dim[0].n = N[0]; dim[0].is = cstride[0]; dim[0].os = rstride[0];
        M_dim.n  = M;    M_dim.is  = cdist;      M_dim.os  = rdist;
        c2r = fftwf_plan_guru_dft_c2r(3, dim, 1, &M_dim, x, (float*)x, FFTW_ESTIMATE);
        if (0 == c2r) goto failed;
    }

    printf("Create FFTW plan for r2c transform\n");
    {
        fftwf_iodim dim[3], M_dim;
        dim[2].n = N[2]; dim[2].is = rstride[2]; dim[2].os = cstride[2];
        dim[1].n = N[1]; dim[1].is = rstride[1]; dim[1].os = cstride[1];
        dim[0].n = N[0]; dim[0].is = rstride[0]; dim[0].os = cstride[0];
        M_dim.n  = M;    M_dim.is  = rdist;      M_dim.os  = cdist;
        r2c = fftwf_plan_guru_dft_r2c(3, dim, 1, &M_dim, (float*)x, x, FFTW_ESTIMATE);
        if (0 == r2c) goto failed;
    }

    printf("Initialize input for c2r FFT\n");
    init_c(x, N, cstride, M, cdist, H);

    printf("Compute c2r FFT\n");
    fftwf_execute(c2r);

    printf("Verify the result of c2r FFT\n");
    status = verify_r((float*)x, N, rstride, M, rdist, H);
    if (0 != status) goto failed;

    printf("Initialize input for r2c FFT\n");
    init_r((float*)x, N, rstride, M, rdist, H);

    printf("Compute r2c FFT using new-array function\n");
    fftwf_execute_dft_r2c(r2c, (float*)x, x);

    printf("Verify the result of r2c FFT\n");
    status = verify_c(x, N, cstride, M, cdist, H);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plans\n");
    fftwf_destroy_plan(r2c);
    fftwf_destroy_plan(c2r);

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
static void init_r(float *x, int *N, int *rs, int M, int rd, int *H)
{
    float TWOPI = 6.2831853071795864769f, phase, factor;
    int n1, n2, n3, m, S1, S2, S3, index;

    S3 = rs[2];
    S2 = rs[1];
    S1 = rs[0];
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
                    index = n1*S1 + n2*S2 + n3*S3 + m*rd;
                    x[index] = factor * cosf( TWOPI * phase ) / (N[0]*N[1]*N[2]);
                }
            }
        }
    }
}

/* Initialize x to be inverse of unit peak at H[0],H[1],H[2] */
static void init_c(fftwf_complex *x, int *N, int *cs, int M, int cd, int *H)
{
    float TWOPI = 6.2831853071795864769f, phase;
    int n1, n2, n3, m, S1, S2, S3, index;

    S3 = cs[2];
    S2 = cs[1];
    S1 = cs[0];
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
                    index = n1*S1 + n2*S2 + n3*S3 + m*cd;
                    x[index][0] =  cosf( TWOPI * phase ) / (N[0]*N[1]*N[2]);
                    x[index][1] = -sinf( TWOPI * phase ) / (N[0]*N[1]*N[2]);
                }
            }
        }
    }
}

/* Verify that x has unit peak at H */
static int verify_c(fftwf_complex *x, int *N, int *cs, int M, int cd, int *H)
{
    float err, errthr, maxerr;
    int n1, n2, n3, m, S1, S2, S3, index;

    S3 = cs[2];
    S2 = cs[1];
    S1 = cs[0];

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

                    index = n1*S1 + n2*S2 + n3*S3 + m*cd;
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

/* Verify that x has unit peak at H */
static int verify_r(float *x, int *N, int *rs, int M, int rd, int *H)
{
    float err, errthr, maxerr;
    int n1, n2, n3, m, S1, S2, S3, index;

    S3 = rs[2];
    S2 = rs[1];
    S1 = rs[0];

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
                for (n3 = 0; n3 < N[2]; n3++)
                {
                    float re_exp = 0.0, re_got;

                    if ((n1-H[0])%N[0]==0 &&
                        (n2-H[1])%N[1]==0 &&
                        (n3-H[2])%N[2]==0)
                    {
                        re_exp = 1;
                    }

                    index = n1*S1 + n2*S2 + n3*S3 + m*rd;
                    re_got = x[index];
                    err  = fabsf(re_got - re_exp);
                    if (err > maxerr) maxerr = err;
                    if (!(err < errthr))
                    {
                        printf(" x(n1=%i,n2=%i,n3=%i,m=%i): ",n1,n2,n3,m);
                        printf(" expected %.7g, ",re_exp);
                        printf(" got %.7g, ",re_got);
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
