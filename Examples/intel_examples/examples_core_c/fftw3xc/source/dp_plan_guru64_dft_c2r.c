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
!       Example of using fftw_plan_guru64_dft_c2r function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

/* Define format to printf ptrdiff_t values used by guru64 interfaces */
#if defined(_WIN32)
#define LI "%Ii" /* on Windows */
#else
#define LI "%zi" /* Otherwise */
#endif

static void init_r(double *x, ptrdiff_t *N, ptrdiff_t *rs, ptrdiff_t M,
                   ptrdiff_t rd, ptrdiff_t *H);
static void init_c(fftw_complex *x, ptrdiff_t *N, ptrdiff_t *rs,ptrdiff_t M,
                   ptrdiff_t rd, ptrdiff_t *H);
static int verify_r(double *x,  ptrdiff_t *N, ptrdiff_t *rs, ptrdiff_t M,
                    ptrdiff_t rd, ptrdiff_t *H);
static int verify_c(fftw_complex *x,  ptrdiff_t *N, ptrdiff_t *rs, ptrdiff_t M,
                    ptrdiff_t rd, ptrdiff_t *H);

int main(void)
{
    /*
     * In this example we show how to compute multiple
     * three-dimensional out-of-place c2r FFTs by one call of FFTW.
     */

    /* Sizes of 3D transform and the number of them.
     * N[2] is the last dimension, it will be halved in complex domain.
     */
    ptrdiff_t N[3] = { 6, 12, 36 };
    ptrdiff_t M = 10;

    /* Arbitrary harmonic used to verify FFT */
    ptrdiff_t H[3] = { 1, 2, 3 };

    /* Strides and distance describe data layout for real and complex domains */
    ptrdiff_t rstride[3], rdist;
    ptrdiff_t cstride[3], cdist;

    /* Size of input and output arrays */
    ptrdiff_t nr, nc;

    /* FFTW plan handles */
    fftw_plan r2c = 0, c2r = 0;

    /* Pointers to input and output data */
    fftw_complex *x = 0;
    double *y = 0;

    /* Execution status */
    int status = 0;

    printf("Example dp_plan_guru64_dft_c2r\n");
    printf("Multiple 3D out-of-place c2r/r2c FFT\n");
    printf("Configuration parameters:\n");
    printf(" N  = {"LI", "LI", "LI"}\n", N[0], N[1], N[2]);
    printf(" M  = "LI"\n", M);
    printf(" H  = {"LI", "LI", "LI"}\n", H[0], H[1], H[2]);

    printf("Define layout of data\n");
    /*
     * Model declaration complex x[ M ][ N[0] ][ N[1] ][ N[2]/2+1 ].
     */
    cstride[2] = 1;
    cstride[1] = (N[2]/2+1) * cstride[2];
    cstride[0] = N[1] * cstride[1];
    cdist      = N[0] * cstride[0];
    nc         = M * cdist;

    /*
     * Model declaration real y[ M ][ N[0] ][ N[1] ][ N[2] ].
     */
    rstride[2] = 1;
    rstride[1] = N[2]*rstride[2];
    rstride[0] = N[1]*rstride[1];
    rdist      = N[0]*rstride[0];
    nr         = M * rdist;

    printf(" rdist="LI", rstride={"LI","LI","LI"}\n",
           rdist, rstride[0],rstride[1],rstride[2]);
    printf(" cdist="LI", cstride={"LI","LI","LI"}\n",
           cdist, cstride[0],cstride[1],cstride[2]);

    printf("Allocate complex x("LI")\n", nc );
    x  = fftw_malloc(sizeof(fftw_complex) * nc );
    printf("Allocate real y("LI")\n", nr );
    y  = fftw_malloc(sizeof(double) * nr );
    if (0 == x || 0 == y) goto failed;

    printf("Create FFTW plan for c2r transform\n");
    {
        fftw_iodim64 dim[3], M_dim;
        dim[2].n = N[2]; dim[2].is = cstride[2]; dim[2].os = rstride[2];
        dim[1].n = N[1]; dim[1].is = cstride[1]; dim[1].os = rstride[1];
        dim[0].n = N[0]; dim[0].is = cstride[0]; dim[0].os = rstride[0];
        M_dim.n  = M;    M_dim.is  = cdist;      M_dim.os  = rdist;
        c2r = fftw_plan_guru64_dft_c2r(3, dim, 1, &M_dim, x, y, FFTW_ESTIMATE);
        if (0 == c2r) goto failed;
    }

    printf("Create FFTW plan for r2c transform\n");
    {
        fftw_iodim64 dim[3], M_dim;
        dim[2].n = N[2]; dim[2].is = rstride[2]; dim[2].os = cstride[2];
        dim[1].n = N[1]; dim[1].is = rstride[1]; dim[1].os = cstride[1];
        dim[0].n = N[0]; dim[0].is = rstride[0]; dim[0].os = cstride[0];
        M_dim.n  = M;    M_dim.is  = rdist;      M_dim.os  = cdist;
        r2c = fftw_plan_guru64_dft_r2c(3, dim, 1, &M_dim, y, x, FFTW_ESTIMATE);
        if (0 == r2c) goto failed;
    }

    printf("Initialize input for c2r FFT\n");
    init_c(x, N, cstride, M, cdist, H);

    printf("Compute c2r FFT\n");
    fftw_execute(c2r);

    printf("Verify the result of c2r FFT\n");
    status = verify_r(y, N, rstride, M, rdist, H);
    if (0 != status) goto failed;

    printf("Initialize input for r2c FFT\n");
    init_r(y, N, rstride, M, rdist, H);

    printf("Compute r2c FFT using new-array function\n");
    fftw_execute_dft_r2c(r2c, y, x);

    printf("Verify the result of r2c FFT\n");
    status = verify_c(x, N, cstride, M, cdist, H);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plans\n");
    fftw_destroy_plan(r2c);
    fftw_destroy_plan(c2r);

    printf("Free data arrays\n");
    fftw_free(x);
    fftw_free(y);

    printf("TEST %s\n",0==status ? "PASSED" : "FAILED");
    return status;

 failed:
    printf(" ERROR\n");
    status = 1;
    goto cleanup;
}

/* Compute (K*L)%M accurately, ptrdiff_t is 32-bit on ia32 */
static double moda(ptrdiff_t K, ptrdiff_t L, ptrdiff_t M)
{
    return (double)(((long long)K * L) % M);
}

/* Initialize arrays x with harmonic H */
static void init_r(double *x, ptrdiff_t *N, ptrdiff_t *rs, ptrdiff_t M,
                   ptrdiff_t rd, ptrdiff_t *H)
{
    double TWOPI = 6.2831853071795864769, phase, factor;
    ptrdiff_t n1, n2, n3, m, S1, S2, S3, index;

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
                    x[index] = factor * cos( TWOPI * phase ) / (N[0]*N[1]*N[2]);
                }
            }
        }
    }
}

/* Initialize x to be inverse of unit peak at H[0],H[1],H[2] */
static void init_c(fftw_complex *x, ptrdiff_t *N, ptrdiff_t *cs, ptrdiff_t M,
                   ptrdiff_t cd, ptrdiff_t *H)
{
    double TWOPI = 6.2831853071795864769, phase;
    ptrdiff_t n1, n2, n3, m, S1, S2, S3, index;

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
                    x[index][0] =  cos( TWOPI * phase ) / (N[0]*N[1]*N[2]);
                    x[index][1] = -sin( TWOPI * phase ) / (N[0]*N[1]*N[2]);
                }
            }
        }
    }
}

/* Verify that x has unit peak at H */
static int verify_c(fftw_complex *x, ptrdiff_t *N, ptrdiff_t *cs, ptrdiff_t M,
                    ptrdiff_t cd, ptrdiff_t *H)
{
    double err, errthr, maxerr;
    ptrdiff_t n1, n2, n3, m, S1, S2, S3, index;

    S3 = cs[2];
    S2 = cs[1];
    S1 = cs[0];

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

                    index = n1*S1 + n2*S2 + n3*S3 + m*cd;
                    re_got = x[index][0];
                    im_got = x[index][1];
                    err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
                    if (err > maxerr) maxerr = err;
                    if (!(err < errthr))
                    {
                        printf(" x(n1="LI",n2="LI",n3="LI",m="LI"): ",n1,n2,n3,m);
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
    printf(" Verified,  maximum error was %.3lg\n", maxerr);
    return 0;
}

/* Verify that x has unit peak at H */
static int verify_r(double *x, ptrdiff_t *N, ptrdiff_t *rs, ptrdiff_t M, ptrdiff_t rd,
                    ptrdiff_t *H)
{
    double err, errthr, maxerr;
    ptrdiff_t n1, n2, n3, m, S1, S2, S3, index;

    S3 = rs[2];
    S2 = rs[1];
    S1 = rs[0];

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

                    index = n1*S1 + n2*S2 + n3*S3 + m*rd;
                    re_got = x[index];
                    err  = fabs(re_got - re_exp);
                    if (err > maxerr) maxerr = err;
                    if (!(err < errthr))
                    {
                        printf(" x(n1="LI",n2="LI",n3="LI",m="LI"): ",n1,n2,n3,m);
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
