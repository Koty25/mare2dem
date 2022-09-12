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
!       Example of using fftwf_plan_r2r_1d function.
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "fftw3.h"
#include "fftw3_mkl.h"

static void init_redft00(float *x, int N, int H);
static void init_redft10(float *x, int N, int H);
static void init_redft01(float *x, int N, int H);
static void init_redft11(float *x, int N, int H);
static void init_rodft00(float *x, int N, int H);
static void init_rodft10(float *x, int N, int H);
static void init_rodft01(float *x, int N, int H);
static void init_rodft11(float *x, int N, int H);
static int verify(float *x, int N, int H);

int main(void)
{
    /* Size of 1D r2r transform, N>1 */
    int N = 8;

    /* Arbitrary harmonic used to verify FFT */
    int H = 1;

    /* FFTW plan handle */
    fftwf_plan r2r = 0;

    /* Pointer to input/output data */
    float *x = 0;

    /* Execution status */
    int status = 0;

    printf("Example sp_plan_r2r_1d\n");
    printf("1D real-to-real in-place transform\n");
    printf("Configuration parameters:\n");
    printf(" N = %d\n", N);
    printf(" H = %d\n", H);

    printf("Allocate data array x(%i)\n", N+1); /* extra element for MKl_RODFT00 kind */
    x  = fftwf_malloc(sizeof(float)*(N+1));
    if (0 == x) goto failed;

    printf("============ REDFT00 (DCT1) ============\n");

    printf("Create FFTW_REDFT00 plan for N=%i\n",N);
    r2r = fftwf_plan_r2r_1d(N, x, x, FFTW_REDFT00, FFTW_ESTIMATE);
    if (0 == r2r) goto failed;

    printf("Initialize input for REDFT00 transform\n");
    init_redft00(x, N, H);

    printf("Compute FFTW_REDFT00\n");
    fftwf_execute(r2r);

    printf("Verify the result of REDFT00\n");
    status = verify(x, N, H);
    if (0 != status) goto failed;

    printf("Destroy FFTW plan\n");
    fftwf_destroy_plan(r2r);

    printf("============ REDFT10 (DCT2) ============\n");

    printf("Create FFTW_REDFT10 plan for N=%i\n",N);
    r2r = fftwf_plan_r2r_1d(N, x, x, FFTW_REDFT10, FFTW_ESTIMATE);
    if (0 == r2r) goto failed;

    printf("Initialize input for REDFT10 transform\n");
    init_redft10(x, N, H);

    printf("Compute FFTW_REDFT10\n");
    fftwf_execute(r2r);

    printf("Verify the result of REDFT10\n");
    status = verify(x, N, H);
    if (0 != status) goto failed;

    printf("Destroy FFTW plan\n");
    fftwf_destroy_plan(r2r);

    printf("============ REDFT01 (DCT3) ============\n");

    printf("Create FFTW_REDFT01 plan for N=%i\n",N);
    r2r = fftwf_plan_r2r_1d(N, x, x, FFTW_REDFT01, FFTW_ESTIMATE);
    if (0 == r2r) goto failed;

    printf("Initialize input for REDFT01 transform\n");
    init_redft01(x, N, H);

    printf("Compute FFTW_REDFT01\n");
    fftwf_execute(r2r);

    printf("Verify the result of REDFT01\n");
    status = verify(x, N, H);
    if (0 != status) goto failed;

    printf("Destroy FFTW plan\n");
    fftwf_destroy_plan(r2r);

    printf("============ REDFT11 (DCT4) ============\n");

    printf("Create FFTW_REDFT11 plan for N=%i\n",N);
    r2r = fftwf_plan_r2r_1d(N, x, x, FFTW_REDFT11, FFTW_ESTIMATE);
    if (0 == r2r) goto failed;

    printf("Initialize input for REDFT11 transform\n");
    init_redft11(x, N, H);

    printf("Compute FFTW_REDFT11\n");
    fftwf_execute(r2r);

    printf("Verify the result of REDFT11\n");
    status = verify(x, N, H);
    if (0 != status) goto failed;

    printf("Destroy FFTW plan\n");
    fftwf_destroy_plan(r2r);

    printf("============ RODFT00 (DST1) ============\n");

    printf("Create FFTW_RODFT00 plan for N=%i\n",N);
    r2r = fftwf_plan_r2r_1d(N, x, x, FFTW_RODFT00, FFTW_ESTIMATE);
    if (0 == r2r) goto failed;

    printf("Initialize input for RODFT00 transform\n");
    init_rodft00(x, N, H);

    printf("Compute FFTW_RODFT00\n");
    fftwf_execute(r2r);

    printf("Verify the result of RODFT00\n");
    status = verify(x, N, H);
    if (0 != status) goto failed;

    printf("Destroy FFTW plan\n");
    fftwf_destroy_plan(r2r);

    printf("============ RODFT10 (DST2) ============\n");

    printf("Create FFTW_RODFT10 plan for N=%i\n",N);
    r2r = fftwf_plan_r2r_1d(N, x, x, FFTW_RODFT10, FFTW_ESTIMATE);
    if (0 == r2r) goto failed;

    printf("Initialize input for RODFT10 transform\n");
    init_rodft10(x, N, H);

    printf("Compute FFTW_RODFT10\n");
    fftwf_execute(r2r);

    printf("Verify the result of RODFT10\n");
    status = verify(x, N, H);
    if (0 != status) goto failed;

    printf("Destroy FFTW plan\n");
    fftwf_destroy_plan(r2r);

    printf("============ RODFT01 (DST3) ============\n");

    printf("Create FFTW_RODFT01 plan for N=%i\n",N);
    r2r = fftwf_plan_r2r_1d(N, x, x, FFTW_RODFT01, FFTW_ESTIMATE);
    if (0 == r2r) goto failed;

    printf("Initialize input for RODFT01 transform\n");
    init_rodft01(x, N, H);

    printf("Compute FFTW_RODFT01\n");
    fftwf_execute(r2r);

    printf("Verify the result of RODFT01\n");
    status = verify(x, N, H);
    if (0 != status) goto failed;

    printf("Destroy FFTW plan\n");
    fftwf_destroy_plan(r2r);

    printf("============ RODFT11 (DST4) ============\n");

    printf("Create FFTW_RODFT11 plan for N=%i\n",N);
    r2r = fftwf_plan_r2r_1d(N, x, x, FFTW_RODFT11, FFTW_ESTIMATE);
    if (0 == r2r) goto failed;

    printf("Initialize input for RODFT11 transform\n");
    init_rodft11(x, N, H);

    printf("Compute FFTW_RODFT11\n");
    fftwf_execute(r2r);

    printf("Verify the result of RODFT11\n");
    status = verify(x, N, H);
    if (0 != status) goto failed;

    printf("Destroy FFTW plan\n");
    fftwf_destroy_plan(r2r);

    printf("============ MKL_RODFT00 (DST1) ============\n");

    printf("Create MKL_RODFT00 plan for N=%i\n",N);
    r2r = fftwf_plan_r2r_1d(N, x, x, (fftwf_r2r_kind)MKL_RODFT00, FFTW_ESTIMATE);
    if (0 == r2r) goto failed;

    printf("Initialize input for RODFT00 transform\n");
    init_rodft00(x+1, N, H);

    printf("Compute RODFT00 using new-array execute function for variety\n");
    fftwf_execute_r2r(r2r, x, x);

    printf("Verify the result of RODFT00\n");
    status = verify(x+1, N, H);
    if (0 != status) goto failed;

 cleanup:

    printf("Destroy FFTW plan\n");
    fftwf_destroy_plan(r2r);

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

/* Compute K%L reduced to 0..L-1 */
static int modulo(int K, int L)
{
    return (K % L + L) % L;
}

/* Initialize array x(N) to produce unit peak at H */
static void init_redft00(float *x, int N, int H)
{
    float PI = 3.14159265358979323846f, factor;
    int n, j, NL;

    j = modulo(H,N);
    NL = 2*(N-1); /* logical size of transform */
    factor = ( j==0 || j==N-1 ) ? 1.0f : 2.0f;
    for (n = 0; n < N; n++)
    {
        x[n] = factor * cosf( 2*PI * moda(n,j,NL) / NL ) / NL;
    }
}

/* Initialize array x(N) to produce unit peak at H */
static void init_redft10(float *x, int N, int H)
{
    float PI = 3.14159265358979323846f, factor;
    int n, j, NL;

    j = modulo(H,N);
    NL = 2*N; /* logical size of transform */
    factor = ( j==0 ) ? 1.0f : 2.0f;
    for (n = 0; n < N; n++)
    {
        x[n] = factor * cosf( PI * moda(2*n+1,j,2*NL) / NL ) / NL;
    }
}

/* Initialize array x(N) to produce unit peak at H */
static void init_redft01(float *x, int N, int H)
{
    float PI = 3.14159265358979323846f, factor;
    int n, j, NL;

    j = modulo(H,N);
    NL = 2*N; /* logical size of transform */
    factor = 2;
    for (n = 0; n < N; n++)
    {
        x[n] = factor * cosf( PI * moda(n,2*j+1,2*NL) / NL ) / NL;
    }
}

/* Initialize array x(N) to produce unit peak at H */
static void init_redft11(float *x, int N, int H)
{
    float PI = 3.14159265358979323846f, factor;
    int n, j, NL;

    j = modulo(H,N);
    NL = 2*N; /* logical size of transform */
    factor = 2;
    for (n = 0; n < N; n++)
    {
        x[n] = factor * cosf( PI * moda(2*n+1,2*j+1,4*NL) / (2*NL) ) / NL;
    }
}

/* Initialize array x(N) to produce unit peak at H */
static void init_rodft00(float *x, int N, int H)
{
    float PI = 3.14159265358979323846f, factor;
    int n, j, NL;

    j = modulo(H,N);
    NL = 2*(N+1); /* logical size of transform */
    factor = 2;
    for (n = 0; n < N; n++)
    {
        x[n] = factor * sinf( 2*PI * moda(n+1,j+1,NL) / NL ) / NL;
    }
}

static void init_rodft10(float *x, int N, int H)
{
    float PI = 3.14159265358979323846f, factor;
    int n, j, NL;

    j = modulo(H,N);
    NL = 2*N; /* logical size of transform */
    factor = (j==N-1) ? 1.0f : 2.0f;
    for (n = 0; n < N; n++)
    {
        x[n] = factor * sinf( PI * moda(2*n+1,j+1,2*NL) / NL ) / NL;
    }
}

static void init_rodft01(float *x, int N, int H)
{
    float PI = 3.14159265358979323846f, factor;
    int n, j, NL;

    j = modulo(H,N);
    NL = 2*N; /* logical size of transform */
    factor = 2;
    for (n = 0; n < N; n++)
    {
        x[n] = factor * sinf( PI * moda(n+1,2*j+1,2*NL) / NL ) / NL;
    }
}

static void init_rodft11(float *x, int N, int H)
{
    float PI = 3.14159265358979323846f, factor;
    int n, j, NL;

    j = modulo(H,N);
    NL = 2*N; /* logical size of transform */
    factor = 2;
    for (n = 0; n < N; n++)
    {
        x[n] = factor * sinf( PI * moda(2*n+1,2*j+1,4*NL) / (2*NL) ) / NL;
    }
}

/* Verify that x has unit peak at H */
static int verify(float *x, int N, int H)
{
    float err, errthr, maxerr;
    int n;

    /*
     * Note, this simple error bound doesn't take into account error of
     * input data
     */
    errthr = 2.5f * logf( (float)N ) / logf(2.0f) * FLT_EPSILON;
    printf(" Check if err is below errthr %.3g\n", errthr);

    maxerr = 0;
    for (n = 0; n < N; n++)
    {
        float re_exp, re_got;

        re_exp = ((n-H)%N == 0 ? 1 : 0);
        re_got = x[n];

        err  = fabsf(re_got - re_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr))
        {
            printf(" x[%i]: ",n);
            printf(" expected %.7g, ",re_exp);
            printf(" got %.7g, ",re_got);
            printf(" err %.3g\n", err);
            printf(" Verification FAILED\n");
            return 1;
        }
    }
    printf(" Verified,  maximum error was %.3g\n", maxerr);
    return 0;
}
