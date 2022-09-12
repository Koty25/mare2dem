/*******************************************************************************
* Copyright 2006-2020 Intel Corporation.
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
!       Intel(R) Math Kernel Library (Intel(R) MKL) DFTI interface example
!       program (C-interface)
!
!       Examples support function set
!
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mkl_dfti_examples.h"

/*
**   Initialize x_in and copy to expected x_exp
*/

int init_input_and_expected_vectors_c(void* x_init, void* res_exp, long n)
{
    float* x_in  = x_init;
    float* x_exp = res_exp;
    long   i;
    float  step, step0;

    /*
    **  Put input data into X(1),...,X(n)
    */
    step = -(float)MATH_PI;
    step0 = (2.0f*(float)MATH_PI)/(float)n;
    step = step + step0;
    for (i = 0; i < n; i++) {
        x_in[2*i]   = sinf(step)*sqrtf(3.0f)/2.0f;
        x_in[2*i+1] = sinf(step)/sqrtf(3.0f);
        step = step + step0;
    }

    /*
    **  Put expected data into X_EXP(1),...,X_EXP(n)
    */
    for (i = 0; i < n; i++) {
        x_exp[2*i]   = x_in[2*i];
        x_exp[2*i+1] = x_in[2*i+1];
    }
    return 0;
}

int init_input_and_expected_vectors_z(void* x_init, void* res_exp, long n)
{
    double* x_in  = x_init;
    double* x_exp = res_exp;
    long    i;
    double  step, step0;

    /*
    **  Put input data into X(1),...,X(n)
    */
    step = -MATH_PI;
    step0 = (2.0*MATH_PI)/(double)n;
    step = step + step0;
    for (i = 0; i < n; i++) {
        x_in[2*i]   = sin(step)*sqrt(3.0)/2.0;
        x_in[2*i+1] = sin(step)/sqrt(3.0);
        step = step + step0;
    }

    /*
    **  Put expected data into X_EXP(1),...,X_EXP(n)
    */
    for (i = 0; i < n; i++) {
        x_exp[2*i]   = x_in[2*i];
        x_exp[2*i+1] = x_in[2*i+1];
    }
    return 0;
}

int init_real_vectors_s(void* x, long n)
{
    float*  x_in = x;
    long    i;
    float   step;

    /*
    **  Put input data into X(1),...,X(n)
    */
    for (i = 0; i < n; i++) {
          step = (float)(i+1);
          x_in[i] = sinf(step)*sqrtf(3.0f)/2.0f;
    }
    return 0;
}

int init_real_vectors_d(void* x, long n)
{
    double* x_in = x;
    long    i;
    double  step;

    /*
    **  Put input data into X(1),...,X(n)
    */
    for (i = 0; i < n; i++) {
          step = (double)(i+1);
          x_in[i] = sin(step)*sqrt(3.0)/2.0;
    }
    return 0;
}

int zero_init_c(void* x, long n)
{
    float*  x_in = x;
    long    i;

    for (i = 0; i < n; i++) {
            x_in[2*i]   = 0.0f;
            x_in[2*i+1] = 0.0f;
    }
    return 0;
}

int zero_init_z(void* x, long n)
{
    double* x_in = x;
    long    i;

    for (i = 0; i < n; i++) {
            x_in[2*i]   = 0.0;
            x_in[2*i+1] = 0.0;
    }
    return 0;
}

int zero_init_s(void* x, long n)
{
    float*  x_in = x;
    long    i;

    for (i = 0; i < n; i++)
        x_in[i] = 0.0f;
    return 0;
}

int zero_init_d(void* x, long n)
{
    double*  x_in = x;
    long     i;

    for (i = 0; i < n; i++)
        x_in[i] = 0.0;
    return 0;
}

int init_multiple_columns_c(void* x, long n, long multiple,
                            long first_index, long step_x)
{
    float*  x_in = x;
    float   step;
    long    i, j;

    x_in += 2*first_index;

    for (i = 0; i < multiple; i++) {
        for (j = 0; j < n; j++) {
            step = (float)(i+1)*(float)(i+j+1);
            x_in[2*j*step_x]   = sinf(step)*sqrtf(3.0f)/2.0f;
            x_in[2*j*step_x+1] = sinf(step)/sqrtf(3.0f);
        }
        x_in += 2;
    }
    return 0;
}

int init_3d_columns_c(void* x, long m, long n, long k,
                      long* strides)
{
    float*  x_2d = x;
    float*  x_in;
    float   step;
    long    i, j, l;
    long    first_index;
    long    step_m, step_n, step_k;

    first_index = 2*strides[0];
    step_m      = 2*strides[1];
    step_n      = 2*strides[2];
    step_k      = 2*strides[3];

    x_2d += first_index;
    for (l = 0; l < k; l++) {
         x_in = x_2d;
         for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                step = (float)(l+1)*(float)(i+1)*(float)(i+j+1);
                x_in[j*step_m]   = sinf(step)*sqrtf(3.0f)/2.0f;
                x_in[j*step_m+1] = sinf(step)/sqrtf(3.0f);
            }
            x_in += step_n;
         }
         x_2d += step_k;
    }
    return 0;
}

int init_3d_columns_z(void* x, long m, long n, long k,
                       long* strides)
{
    double* x_2d = x;
    double* x_in;
    double  step;
    long    i, j, l;
    long    first_index;
    long    step_m, step_n, step_k;

    first_index = 2*strides[0];
    step_m      = 2*strides[1];
    step_n      = 2*strides[2];
    step_k      = 2*strides[3];

    x_2d += first_index;
    for (l = 0; l < k; l++) {
         x_in = x_2d;
         for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                step = (double)(l+1)*(double)(i+1)*(double)(i+j+1);
                x_in[j*step_m]   = sin(step)*sqrt(3.0)/2.0;
                x_in[j*step_m+1] = sin(step)/sqrt(3.0);
            }
            x_in += step_n;
         }
         x_2d += step_k;
    }
    return 0;
}

int init_multiple_columns_z(void* x, long n, long multiple,
                             long first_index, long step_x)
{
    double*  x_in = x;
    double   step;
    long     i, j;

    x_in += 2*first_index;
    for (i = 0; i < multiple; i++) {
        for (j = 0; j < n; j++) {
            step = (double)(i+j+1);
            x_in[2*j*step_x]   = sin(step)*sqrt(3.0)/2.0;
            x_in[2*j*step_x+1] = sin(step)/sqrt(3.0);
        }
        x_in += 2;
    }
    return 0;
}

void init_multiple_columns_step_z(void* x, long n, long multiple,
                                  long step_y, long step_x) {
    double*  x_in = x;
    double   step;
    long     i, j;

    for (i = 0; i < multiple; i++) {
        for (j = 0; j < n; j++) {
            step = (double)(i+j+1);
            x_in[2*j*step_x]   = sin(step)*sqrt(3.0)/2.0;
            x_in[2*j*step_x+1] = sin(step)/sqrt(3.0);
        }
        x_in += 2*step_y;
    }
    return;
}
void init_multiple_columns_step_c(void* x, long n, long multiple,
                             long step_y, long step_x) {
    float*  x_in = x;
    float   step;
    long    i, j;

    for (i = 0; i < multiple; i++) {
        for (j = 0; j < n; j++) {
            step = (float)(i+j+1);
            x_in[2*j*step_x]   = sinf(step)*sqrtf(3.0f)/2.0f;
            x_in[2*j*step_x+1] = sinf(step)/sqrtf(3.0f);
        }
        x_in += 2*step_y;
    }
    return;
}

/*
**   Check result
*/

float check_result_c(void* x, void* res_exp, long n)
{
    float* x_in = x;
    float* x_exp = res_exp;
    float  maxerr, d;
    long   i;

    maxerr = 0;
    for (i = 0; i < n; i++) {
        d = x_exp[2*i]   - x_in[2*i];   if (d < 0) d = -d; if (d > maxerr) maxerr = d;
        d = x_exp[2*i+1] - x_in[2*i+1]; if (d < 0) d = -d; if (d > maxerr) maxerr = d;
    }
    return maxerr;
}


double check_result_z(void* x, void* res_exp, long n)
{
    double* x_in = x;
    double* x_exp = res_exp;
    double  maxerr, d;
    long    i;

    maxerr = 0;
    for (i = 0; i < n; i++) {
        d = x_exp[2*i]   - x_in[2*i];   if (d < 0) d = -d; if (d > maxerr) maxerr = d;
        d = x_exp[2*i+1] - x_in[2*i+1]; if (d < 0) d = -d; if (d > maxerr) maxerr = d;
    }

    return maxerr;
}

float check_result_s(void* x, void* y, long n)
{
    float*  x_in = x;
    float*  x_exp = y;
    float   maxerr, d;
    long    i;

    maxerr = 0;
    for (i = 0; i < n; i++) {
        d = x_exp[i] - x_in[i];
        if (d < 0)
            d = -d;
        if (d > maxerr)
            maxerr = d;
    }
    return maxerr;
}

double check_result_d(void* x, void* y, long n)
{
    double* x_in = x;
    double* x_exp = y;
    double  maxerr, d;
    long    i;

    maxerr = 0;
    for (i = 0; i < n; i++) {
        d = x_exp[i] - x_in[i];
        if (d < 0)
            d = -d;
        if (d > maxerr)
            maxerr = d;
    }
    return maxerr;
}

TYPE_PRECISION check_result_multiple(void* x, void* y, int n, int multiple, int istride, int idist)
{
    TYPE_PRECISION* x_in = x;
    TYPE_PRECISION* x_exp = y;
    TYPE_PRECISION  maxerr, d;
    int i, j, id, is;

    maxerr = 0;
    id = 0;
    for (i=0; i<multiple; i++) {
        for (j=0; j<n; j++) {
            is = id+j*istride;
            d = x_exp[is] - x_in[is];
            if (d < 0)
                d = -d;
            if (d > maxerr)
                maxerr = d;
        }
        id = id + idist;
    }
    return maxerr;
}

TYPE_PRECISION check_result_2d_r(TYPE_PRECISION* x_in, TYPE_PRECISION* x_exp, int m, int n)
{
    TYPE_PRECISION maxerr = 0, d = 0;
    int i = 0, j = 0, nf = 0;

    nf = 2*(n/2+1);

    for (j = 0; j < m; j++)
    {
        for (i = 0; i < n; i++)
        {
            d = x_exp[i+j*nf] - x_in[i+j*nf];
            if (d < 0) d = -d;
            if (d > maxerr) maxerr = d;
        }
    }
    return maxerr;
}

TYPE_PRECISION check_result_2d_r_multiple(TYPE_PRECISION* x_in, TYPE_PRECISION* x_exp, int m, int n,
                                          int multiple, int istride, int idist)
{
    TYPE_PRECISION maxerr = 0, d = 0;
    int i, j, nf, im;

    nf = 2*(n/2+1);

    for (im = 0; im < multiple; im++)
    {
        for (j = 0; j < m; j++)
        {
            for (i = 0; i < n; i++)
            {
                d = x_exp[im*idist + (i+j*nf)*istride] - x_in[im*idist + (i+j*nf)*istride];
                if (d < 0) d = -d;
                if (d > maxerr) maxerr = d;
            }
        }
    }
    return maxerr;
}

TYPE_PRECISION check_result_3d_r(TYPE_PRECISION* x_in, TYPE_PRECISION* x_exp, int m, int n, int k)
{
    TYPE_PRECISION maxerr = 0, d = 0;
    int i, j, nf, s, nf1;

    nf = 2*(k/2+1);
    nf1 = 2*(k/2+1)*n;

    for (s = 0; s < m; s++)
    {
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < k; i++)
            {
                d = x_exp[i+j*nf+s*nf1] - x_in[i+j*nf+s*nf1];
                if (d < 0) d = -d;
                if (d > maxerr) {
                    maxerr = d;
                }
            }
        }
    }
    return maxerr;
}

TYPE_PRECISION check_result_3d_r_multiple(TYPE_PRECISION* x_in, TYPE_PRECISION* x_exp, int m, int n, int k,
                                          int multiple, int istride, int idist)
{
    TYPE_PRECISION maxerr = 0, d = 0;
    int i, j, nf, s, nf1, im;

    nf = 2*(k/2+1);
    nf1 = 2*(k/2+1)*n;

    for (im = 0; im < multiple; im++) {
        for (s = 0; s < m; s++) {
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < k; i++)
                {
                    d = x_exp[im*idist + (i+j*nf+s*nf1)*istride] - x_in[im*idist + (i+j*nf+s*nf1)*istride];
                    if (d < 0) d = -d;
                    if (d > maxerr) {
                        maxerr = d;
                    }
                }
            }
        }
    }
    return maxerr;
}

void scaling(void* xx, TYPE_PRECISION Scale, int n)
{
    int i;
    TYPE_PRECISION* x_in;

    x_in = (TYPE_PRECISION*)xx;

    for (i=0; i<2*n; i++) {

        x_in[i] = x_in[i]*Scale;
    }
    return;
}

void scaling_multiple(void* xx, TYPE_PRECISION Scale,
                      int n, int multiple, int istride, int idist)
{
    int i, j, id, is;
    TYPE_PRECISION* x_in;

    x_in = (TYPE_PRECISION*)xx;
    id = 0;
    for (i=0; i<multiple; i++) {
        for (j=0; j<n; j++) {
            is = id+j*2*istride;
            x_in[is] = x_in[is]*Scale;
            x_in[is+1] = x_in[is+1]*Scale;
        }
        id = id + 2*idist;
    }
    return;
}

void scaling_r_multiple(void* xx, TYPE_PRECISION Scale,
                        int n, int multiple, int istride, int idist)
{
    int i, j, id, is;
    TYPE_PRECISION* x_in;

    x_in = (TYPE_PRECISION*)xx;
    id = 0;
    for (i=0; i<multiple; i++) {
        for (j=0; j<n; j++) {
            is = id+j*istride;
            x_in[is] = x_in[is]*Scale;
        }
        id = id + idist;
    }
    return;
}

void scaling_multiple_2d(void* xx, TYPE_PRECISION Scale, int m, int n,
                         int multiple, int stride, int dist)
{

    int i, j, im, id, is;
    TYPE_PRECISION* x_in;

    x_in = (TYPE_PRECISION*)xx;
    id = 0;

    for (im=0; im<multiple; im++) {
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                is = id + (i*n+j)*2*stride;
                x_in[is] = x_in[is]*Scale;
                x_in[is+1] = x_in[is+1]*Scale;
            }
        }
        id = id + 2*dist;
    }
    return;
}

void scaling_2d_r(TYPE_PRECISION* x_in, TYPE_PRECISION Scale, int m, int n)
{

    int i, j, nf;

    nf = 2*(n/2+1);

    for (j = 0; j < m; j++) {
        for (i = 0; i < n; i++) {
            x_in[i+j*nf] *= Scale;
        }
    }
    return;
}

void scaling_2d_r_multiple(TYPE_PRECISION* x_in, TYPE_PRECISION Scale, int m, int n,
                           int multiple, int stride, int dist)
{
    int i, j, nf, im;

    nf = 2*(n/2+1);

    for (im = 0; im < multiple; im++)
    {
        for (j = 0; j < m; j++)
        {
            for (i = 0; i < n; i++)
            {
                x_in[im*dist + (i+j*nf)*stride] *= Scale;
            }
        }
    }
    return;
}

void scaling_3d_r(TYPE_PRECISION* x_in, TYPE_PRECISION Scale, int m, int n, int k)
{
    int i, j, nf, s, nf1;

    nf = 2*(k/2+1);
    nf1 = 2*(k/2+1)*n;

    for (s = 0; s < m; s++) {
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < k; i++)
            {
                x_in[i+j*nf+s*nf1] *= Scale;
            }
        }
    }
    return;
}

void scaling_3d_r_multiple(TYPE_PRECISION* x_in, TYPE_PRECISION Scale, int m, int n, int k,
                           int multiple, int stride, int dist)
{
    int i, j, nf, s, nf1, im;

    nf = 2*(k/2+1);
    nf1 = 2*(k/2+1)*n;

    for (im = 0; im < multiple; im++) {
        for (s = 0; s < m; s++) {
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < k; i++)
                {
                    x_in[im*dist + (i + j*nf + s*nf1)*stride] *= Scale;
                }
            }
        }
    }
    return;
}

void scaling_r(void* xx, TYPE_PRECISION Scale, int n)
{
    int i;
    TYPE_PRECISION* x_in;

    x_in = (TYPE_PRECISION*)xx;

    for (i=0; i<n; i++) {

        x_in[i] = x_in[i]*Scale;
    }
    return;
}

int init_3d_columns_s(void* x, long m, long n, long k, long* strides)
{
    float*  x_2d = x;
    float*  x_in;
    float   step;
    long    i, j, l;
    long    first_index;
    long    step_m, step_n, step_k;

    first_index = strides[0];
    step_m      = strides[1];
    step_n      = strides[2];
    step_k      = strides[3];

    x_2d += first_index;
    for (l = 0; l < k; l++) {
         x_in = x_2d;
         for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                step = (float)(l+1)*(float)(i+1)*(float)(i+j+1);
                x_in[j*step_m] = sinf(step)*sqrtf(3.0f)/2.0f;
            }
            x_in += step_n;
         }
         x_2d += step_k;
    }
    return 0;
}

int init_3d_columns_d(void* x, long m, long n, long k, long* strides)
{
    double*  x_2d = x;
    double*  x_in;
    double   step;
    long     i, j, l;
    long     first_index;
    long     step_m, step_n, step_k;

    first_index = strides[0];
    step_m      = strides[1];
    step_n      = strides[2];
    step_k      = strides[3];

    x_2d += first_index;
    for (l = 0; l < k; l++) {
         x_in = x_2d;
         for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                step = (double)(l+1)*(double)(i+1)*(double)(i+j+1);
                x_in[j*step_m] = sin(step)*sqrt(3.0)/2.0;
            }
            x_in += step_n;
        }
        x_2d += step_k;
    }
    return 0;
}

void die_unless(int expr, const char *file, int line, const char *expr_str)
{
    if (!expr)
    {
        fprintf(stderr,"Expectation '%s' failed at %s:%i\n", expr_str, file, line);
        fprintf(stderr,"TEST FAILED\n");
        exit(1);
    }
}
