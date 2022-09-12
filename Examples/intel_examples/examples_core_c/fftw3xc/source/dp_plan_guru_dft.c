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
!       Example of using fftw_plan_guru_dft function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"

static void init_i(fftw_complex *x, fftw_iodim *ndims, fftw_iodim *vdims,
                   int H1, int H2, int H3);
static int verify_o(fftw_complex *x, fftw_iodim *ndims, fftw_iodim *vdims,
                    int H1, int H2, int H3);
static int  isize2(int n, fftw_iodim *ndims, int m, fftw_iodim *vdims);
static int  osize2(int n, fftw_iodim *ndims, int m, fftw_iodim *vdims);
static void swapio2(int n, fftw_iodim *ndims, int m, fftw_iodim *vdims);
static void print_dims(int n,fftw_iodim *dims);

int main(void)
{
    /*
     * In this example we show how to compute multiple
     * three-dimensional FFTs by one call of FFTW, given that input
     * and output data may have different layouts.
     */

    /* Sizes of 3D transform and the number of them */
    int N1 = 10, N2 = 20, N3 = 30;
    int M = 4;

    /* Arbitrary harmonic used to verify FFT */
    int H1 = -1, H2 = -2, H3 = -3;

    /* FFTW plan handles */
    fftw_plan fwd = 0, bwd = 0;

    /* FFTW iodims used to define the data layout */
    fftw_iodim dim[3], M_dim;

    /* Pointers to input and output data */
    fftw_complex *x = 0, *y = 0;

    /* Execution status */
    int status = 0;


    printf("Example dp_plan_guru_dft\n");
    printf("Forward and backward multiple 3D complex out-of-place FFT\n");
    printf("Configuration parameters:\n");
    printf(" N  = {%d, %d, %d}\n", N1, N2, N3);
    printf(" M  = %d\n", M);
    printf(" H  = {%d, %d, %d}\n", H1, H2, H3);

    printf("Define iodims for forward transform\n");
    dim[0].n    = N1;
    dim[1].n    = N2;
    dim[2].n    = N3;
    M_dim.n     = M;

    /*
     * Input array x: element x(n1,n2,n3,m) is located at
     * x[n1*dim[0].is+n2*dim[1].is+n3*dim[2].is+m*M_dim.is].
     *
     * We define the input strides as if x was declared
     * x[M][N1][N2][N3], but other layouts are possible.
     */
    dim[2].is   = 1;
    dim[1].is   = N3;
    dim[0].is   = N2*N3;
    M_dim.is    = N1*N2*N3;

    /*
     * Output array y: element y(n1,n2,n3,m) is located at
     * y[n1*dim[0].os+n2*dim[1].os+n3*dim[2].os+m*M_dim.os].
     *
     * We define the output strides as if y was declared
     * y[N3][N2][N1][M], but other layouts are possible.
     */
    M_dim.os    = 1;
    dim[0].os   = M;
    dim[1].os   = N1*M;
    dim[2].os   = N2*N1*M;

    printf(" dim   = ");  print_dims(3,dim); printf("\n");
    printf(" M_dim = ");  print_dims(1,&M_dim); printf("\n");

    printf("Allocate x(%i)\n", isize2(3,dim,1,&M_dim) );
    x  = fftw_malloc(sizeof(fftw_complex) * isize2(3,dim,1,&M_dim) );
    if (0 == x) goto failed;

    printf("Allocate y(%i)\n", osize2(3,dim,1,&M_dim) );
    y  = fftw_malloc(sizeof(fftw_complex) * osize2(3,dim,1,&M_dim) );
    if (0 == y) goto failed;

    printf("Create FFTW plan for forward transform\n");
    fwd = fftw_plan_guru_dft(3, dim, 1, &M_dim, x, y, FFTW_FORWARD,
                             FFTW_ESTIMATE);
    if (0 == fwd) goto failed;

    printf("Create FFTW plan for backward transform\n");
    swapio2(3,dim,1,&M_dim);
    bwd = fftw_plan_guru_dft(3, dim, 1, &M_dim, y, x, FFTW_BACKWARD,
                             FFTW_ESTIMATE);
    if (0 == bwd) goto failed;
    swapio2(3,dim,1,&M_dim);

    printf("Initialize input for forward transform\n");
    init_i(x, dim, &M_dim, H1,H2,H3);

    printf("Compute forward FFT\n");
    fftw_execute(fwd);

    printf("Verify the result of forward FFT\n");
    status = verify_o(y, dim, &M_dim, H1,H2,H3);
    if (0 != status) goto failed;

    printf("Initialize input for backward transform\n");
    swapio2( 3, dim, 1, &M_dim );
    init_i(y, dim, &M_dim, -H1,-H2,-H3);
    swapio2( 3, dim, 1, &M_dim );

    printf("Compute backward transform\n");
    fftw_execute(bwd);

    printf("Verify the result of backward FFT\n");
    swapio2( 3, dim, 1, &M_dim );
    status = verify_o(x, dim, &M_dim, H1,H2,H3);
    if (0 != status) goto failed;
    swapio2( 3, dim, 1, &M_dim );

 cleanup:

    printf("Destroy FFTW plans\n");
    fftw_destroy_plan(fwd);
    fftw_destroy_plan(bwd);

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

/* Compute the input size, in elements */
static int  isize2(int n, fftw_iodim *ndims, int m, fftw_iodim *vdims)
{
    int i, res = 1;
    for (i = 0; i < n; ++i)
    {
        if (res < ndims[i].n*ndims[i].is) res = ndims[i].n*ndims[i].is;
    }
    for (i = 0; i < m; ++i)
    {
        if (res < vdims[i].n*vdims[i].is) res = vdims[i].n*vdims[i].is;
    }
    return res;
}

/* Compute the output size, in elements */
static int  osize2(int n, fftw_iodim *ndims, int m, fftw_iodim *vdims)
{
    int i, res = 1;
    for (i = 0; i < n; ++i)
    {
        if (res < ndims[i].n*ndims[i].os) res = ndims[i].n*ndims[i].os;
    }
    for (i = 0; i < m; ++i)
    {
        if (res < vdims[i].n*vdims[i].os) res = vdims[i].n*vdims[i].os;
    }
    return res;
}

/* Swap i<->o layouts */
static void swapio2(int n, fftw_iodim *ndims, int m, fftw_iodim *vdims)
{
    int i, t;
    for (i = 0; i < n; ++i)
    {
        t = ndims[i].is;
        ndims[i].is = ndims[i].os;
        ndims[i].os = t;
    }
    for (i = 0; i < m; ++i)
    {
        t = vdims[i].is;
        vdims[i].is = vdims[i].os;
        vdims[i].os = t;
    }
}

/* Print iodims */
static void print_dims(int n, fftw_iodim *dims)
{
    int i;
    printf("%i:%i:%i", dims[0].n, dims[0].is, dims[0].os);
    for (i=1; i<n; ++i)
    {
        printf("x%i:%i:%i", dims[i].n, dims[i].is, dims[i].os);
    }
}

/* Compute (K*L)%M accurately */
static double moda(int K, int L, int M)
{
    return (double)(((long long)K * L) % M);
}

/* Initialize arrays x with harmonic H, using *.is for indexing */
static void init_i(fftw_complex *x, fftw_iodim *ndims, fftw_iodim *vdims,
                   int H1,int H2,int H3)
{
    double TWOPI = 6.2831853071795864769, phase;
    int n1, n2, n3, m, N1, N2, N3, M, S1, S2, S3, SM, index;

    N1 = ndims[0].n; S1 = ndims[0].is;
    N2 = ndims[1].n; S2 = ndims[1].is;
    N3 = ndims[2].n; S3 = ndims[2].is;
    M  = vdims[0].n; SM = vdims[0].is;
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

/* Verify that x has unit peak at H, using *.os for indexing */
static int verify_o(fftw_complex *x, fftw_iodim *ndims, fftw_iodim *vdims,
                    int H1,int H2,int H3)
{
    double err, errthr, maxerr;
    int n1, n2, n3, m, N1, N2, N3, M, S1, S2, S3, SM, index;

    N1 = ndims[0].n; S1 = ndims[0].os;
    N2 = ndims[1].n; S2 = ndims[1].os;
    N3 = ndims[2].n; S3 = ndims[2].os;
    M  = vdims[0].n; SM = vdims[0].os;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 5.0 * log( (double)N1*N2*N3 ) / log(2.0) * DBL_EPSILON;
    printf(" Verify the result, errthr = %.3lg\n", errthr);

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
