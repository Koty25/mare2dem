/*******************************************************************************
* Copyright 2003-2020 Intel Corporation.
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

#include <stdlib.h>

/*
//  Generate double precision data set using VSL Uniform RNG
//  with left bound a and right bound b
*/
int dGenerateUniformData(double *x, MKL_INT n, double a, double b)
{
    VSLStreamStatePtr stream;
    MKL_INT brng = VSL_BRNG_WH;
    MKL_INT seed = 777;
    int errcode;

    MKL_INT method = VSL_RNG_METHOD_UNIFORM_STD;

    /***** Initialize *****/
    errcode = vslNewStream( &stream, brng,  seed );
    CheckVslError(errcode);

    /***** Call RNG *****/
    errcode = vdRngUniform( method, stream, n, x, a, b );
    CheckVslError(errcode);

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError(errcode);

    return errcode;
}

/*
//  Generate single precision data set using VSL Uniform RNG
//  with left bound a and right bound b
*/
int sGenerateUniformData(float* x, MKL_INT n, float a, float b)
{
    VSLStreamStatePtr stream;
    MKL_INT brng = VSL_BRNG_WH;
    MKL_INT seed = 777;
    int errcode;

    MKL_INT method = VSL_RNG_METHOD_UNIFORM_STD;

    /***** Initialize *****/
    errcode = vslNewStream( &stream, brng,  seed );
    CheckVslError(errcode);

    /***** Call RNG *****/
    errcode = vsRngUniform( method, stream, n, x, a, b );
    CheckVslError(errcode);

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError(errcode);

    return errcode;
}


/*
//  Generate double precision data set using VSL Gaussian RNG
//  with mean a and standard deviation sigma
*/
int dGenerateGaussianData(double* x, MKL_INT dim, MKL_INT n, double a, double sigma)
{
    VSLStreamStatePtr stream;
    MKL_INT brng = VSL_BRNG_WH;
    MKL_INT seed = 777;
    int errcode;

    MKL_INT method = VSL_RNG_METHOD_GAUSSIAN_ICDF;

    /***** Initialize *****/
    errcode = vslNewStream( &stream, brng,  seed );
    CheckVslError(errcode);

    /***** Call RNG *****/
    errcode = vdRngGaussian( method, stream, n*dim, x, a, sigma );
    CheckVslError(errcode);

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError(errcode);

    return errcode;
}

/*
//  Generate single precision data set using VSL Gaussian RNG
//  with mean a and standard deviation sigma
*/
int sGenerateGaussianData(float* x, MKL_INT dim, MKL_INT n, float a, float sigma)
{
    VSLStreamStatePtr stream;
    MKL_INT brng = VSL_BRNG_WH;
    MKL_INT seed = 777;
    int errcode;

    MKL_INT method = VSL_RNG_METHOD_GAUSSIAN_ICDF;

    /***** Initialize *****/
    errcode = vslNewStream( &stream, brng,  seed );
    CheckVslError(errcode);

    /***** Call RNG *****/
    errcode = vsRngGaussian( method, stream, n*dim, x, a, sigma );
    CheckVslError(errcode);

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError(errcode);

    return errcode;
}

/*
//  Generate double precision data set using VSL GaussianMV RNG
//  with mean a and covariation matrix C
*/
int dGenerateGaussianMVData(double* x, MKL_INT dim, MKL_INT n, double* a, double* C)
{
    VSLStreamStatePtr stream;
    MKL_INT brng = VSL_BRNG_WH;
    MKL_INT seed = 777;
    int errcode = 0;

    MKL_INT method = VSL_RNG_METHOD_GAUSSIANMV_ICDF;

    /***** Following variables are used in Cholesky factorization routine *****/
    char    uplo = 'U';
    MKL_INT num = dim;
    double* T;
    MKL_INT lda = dim;
    MKL_INT info;

    int i;

    T = (double*)malloc((dim * dim) * sizeof(double));
    if (! T) {
        /* Memory allocation failed */
        return -1;
    }

    for(i = 0; i < dim*dim; i++)
    {
        T[i] = C[i];
    }

    /***** Intel(R) MKL Choelsky factorization routine call *****/
    dpotrf( &uplo, &num, T, &lda, &info );

    /***** Initialize *****/
    errcode = vslNewStream( &stream, brng,  seed );
    CheckVslError(errcode);

    /***** Call RNG *****/
    errcode = vdRngGaussianMV( method, stream, n, x, dim,
                               VSL_MATRIX_STORAGE_FULL, a, T );
    CheckVslError(errcode);

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError(errcode);

    free(T);
    return errcode;
}

/*
//  Generate single precision data set using VSL GaussianMV RNG
//  with mean a and covariation matrix C
*/
int sGenerateGaussianMVData(float* x, MKL_INT dim, MKL_INT n, float* a, float* C)
{
    VSLStreamStatePtr stream;
    MKL_INT brng = VSL_BRNG_WH;
    MKL_INT seed = 777;
    int errcode = 0;

    MKL_INT method = VSL_RNG_METHOD_GAUSSIANMV_ICDF;

    /***** Following variables are used in Cholesky factorization routine *****/
    char    uplo = 'U';
    MKL_INT num = dim;
    float* T;
    MKL_INT lda = dim;
    MKL_INT info;

    int i;

    T = (float*)malloc((dim * dim) * sizeof(float));
    if (! T) {
        /* Memory allocation failed */
        return -1;
    }

    for(i = 0; i < dim*dim; i++)
    {
        T[i] = C[i];
    }

    /***** Intel(R) MKL Choelsky factorization routine call *****/
    spotrf( &uplo, &num, (float*)T, &lda, &info );

    /***** Initialize *****/
    errcode = vslNewStream( &stream, brng,  seed );
    CheckVslError(errcode);

    /***** Call RNG *****/
    errcode = vsRngGaussianMV( method, stream, n, x, dim, VSL_MATRIX_STORAGE_FULL, a, T );
    CheckVslError(errcode);

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError(errcode);

    free(T);
    return errcode;
}

/*
// Make Cholesky factorization of covariance matrix C and create new
// stream
*/
VSLStreamStatePtr dInitGaussianMVDataGenerator(MKL_INT dim,
                                  double* C, double* T)
{
    VSLStreamStatePtr stream;
    MKL_INT brng = VSL_BRNG_WH;
    MKL_INT seed = 777;
    int errcode = 0;

    /***** Following variables are used in Cholesky factorization routine *****/
    char    uplo = 'U';
    MKL_INT num = dim;
    MKL_INT lda = dim;
    MKL_INT info;

    int i;

    for(i = 0; i < dim*dim; i++)
    {
        T[i] = C[i];
    }

    /***** Intel(R) MKL Choelsky factorization routine call *****/
    dpotrf( &uplo, &num, T, &lda, &info );

    /***** Initialize *****/
    errcode = vslNewStream( &stream, brng,  seed );
    CheckVslError(errcode);

    return stream;
}

/*
// Make Cholesky factorization of covariance matrix C and create new
// stream
*/
VSLStreamStatePtr sInitGaussianMVDataGenerator(MKL_INT dim,
                                  float* C, float* T)
{
    VSLStreamStatePtr stream;
    MKL_INT brng = VSL_BRNG_WH;
    MKL_INT seed = 777;
    int errcode = 0;

    /***** Following variables are used in Cholesky factorization routine *****/
    char    uplo = 'U';
    MKL_INT num = dim;
    MKL_INT lda = dim;
    MKL_INT info;

    int i;

    for(i = 0; i < dim*dim; i++)
    {
        T[i] = C[i];
    }

    /***** Intel(R) MKL Choelsky factorization routine call *****/
    spotrf( &uplo, &num, T, &lda, &info );

    /***** Initialize *****/
    errcode = vslNewStream( &stream, brng,  seed );
    CheckVslError(errcode);

    return stream;
}

/*
//  Generate double precision data block using VSL GaussianMV RNG
*/
int dGenerateGaussianMVDataBlock(VSLStreamStatePtr stream, double* x,
                                 MKL_INT dim, MKL_INT n, double* a, double* T)
{
    int errcode = 0;
    MKL_INT method = VSL_RNG_METHOD_GAUSSIANMV_ICDF;

    /***** Call RNG *****/
    errcode = vdRngGaussianMV( method, stream, n, x, dim,
                               VSL_MATRIX_STORAGE_FULL, a, T );
    CheckVslError(errcode);

    return errcode;
}

/*
//  Generate single precision data block using VSL GaussianMV RNG
*/
int sGenerateGaussianMVDataBlock(VSLStreamStatePtr stream, float* x,
                                  MKL_INT dim, MKL_INT n, float* a, float* T)
{
    int errcode = 0;
    MKL_INT method = VSL_RNG_METHOD_GAUSSIANMV_ICDF;

    /***** Call RNG *****/
    errcode = vsRngGaussianMV( method, stream, n, x, dim,
                               VSL_MATRIX_STORAGE_FULL, a, T );
    CheckVslError(errcode);

    return errcode;
}

/*
//  Generates contaminated double precision data set using
//  GaussianMV RNG with mean a and covariance C for the majotity of values.
//  The outliers are generated with ratio 'ratio' using GaussianMV RNG with
//  mean m and covariance coeff * C.
*/
int dGenerateContaminatedDataset(double* x, MKL_INT dim, MKL_INT n, double* a,
                              double* C, int ratio, double m, double coeff)
{
    double *buffer;
    double *r1, *r2, *mean, *T;
    int i, j, mc, errcode = VSL_ERROR_OK;

    /***** Compute number of outliers *****/
    mc = (n * ratio) / 100;

    /***** Allocate memory *****/
    buffer = (double*)malloc((dim*(n-mc) + dim*mc + dim + dim*dim)*sizeof(double));
    if (! buffer) {
        /* Memory allocation failed */
        printf("Error : Memory allocation failed\n");
        return -1;
    }
    r1   = buffer;
    r2   = r1 + dim * (n-mc);
    mean = r2 + dim * mc;
    T    = mean + dim;



    for(i = 0; i < dim * dim; i++)
    {
        T[i] = C[i];
    }

    /* Generate dataset from major distribution */
    if ( n-mc ) {
        errcode = dGenerateGaussianMVData( r1, dim, n-mc, a, T );
        CheckVslError( errcode );
    }

    /* Generate outliers */
    if ( mc ) {
        for(i = 0; i < dim; i++)
        {
            mean[i] = m;
        }

        for(i = 0; i < dim * dim; i++)
        {
            T[i] = coeff * C[i];
        }

        errcode = dGenerateGaussianMVData( r2, dim, mc, mean, T );
        CheckVslError( errcode );
    }

    /* Form final dataset */
    for(i = 0; i < n-mc; i++)
    {
        for(j = 0; j < dim; j++)
        {
            x[j * n + i] = r1[i * dim + j];
        }
    }


    for(i = 0; i < mc; i++)
    {
        for(j = 0; j < dim; j++)
        {
            x[j * n + (i + n-mc)] = r2[i * dim + j];
        }
    }

    free(buffer);

    return errcode;
}

/*
//  Generate contaminated single precision data set using
//  GaussianMV RNG with mean a and covariance C for the majotity of values.
//  The outliers are generated with ratio 'ratio' using GaussianMV RNG with
//  mean m and covariance coeff*C.
*/
int sGenerateContaminatedDataset(float* x, MKL_INT dim, MKL_INT n, float* a,
                                 float* C, int ratio, float m, float coeff)
{
    float *buffer;
    float *r1, *r2, *mean, *T;
    int i, j, mc, errcode = VSL_ERROR_OK;

    /***** Compute number of outliers *****/
    mc = (n * ratio) / 100;

    /***** Allocate memory *****/
    buffer = (float*)malloc((dim*(n-mc) + dim*mc + dim + dim*dim) *
                            sizeof(float));
    if (! buffer) {
        /* Memory allocation failed */
        printf("Error : Memory allocation failed\n");
        return -1;
    }
    r1   = buffer;
    r2   = r1 + dim * (n-mc);
    mean = r2 + dim * mc;
    T    = mean + dim;

    for(i = 0; i < dim * dim; i++)
    {
        T[i] = C[i];
    }

    /* Generate dataset from major distribution */
    if ( n-mc ) {
        errcode = sGenerateGaussianMVData( r1, dim, n - mc, a, T );
        CheckVslError( errcode );
    }

    /* Generate outliers */
    if ( mc ) {
        for(i = 0; i < dim; i++)
        {
            mean[i] = m;
        }

        for(i = 0; i < dim * dim; i++)
        {
            T[i] = coeff * C[i];
        }

        errcode = sGenerateGaussianMVData( r2, dim, mc, mean, T );
        CheckVslError( errcode );
    }

    /* Form final dataset */
    for(i = 0; i < n-mc; i++)
    {
        for(j = 0; j < dim; j++)
        {
            x[j * n + i] = r1[i * dim + j];
        }
    }


    for(i = 0; i < mc; i++)
    {
        for(j = 0; j < dim; j++)
        {
            x[j * n + (i + n - mc)] = r2[i * dim + j];
        }
    }

    free(buffer);

    return errcode;
}

/*

*/
int dGenerateMissingValuesInput( double *x, MKL_INT dim, MKL_INT n, int epsilon,
                                 int *patt, int npatt, int* miss_vals,
                                 double *a, double *C,
                                 MKL_INT *n_miss_vals )
{
    double av_weight, threshold, eps;
    long long iNAN;
    double    dNAN;
    double floor;
    MKL_INT nmv = 0;
    int i, j, ipatt, errcode;
    double *xt;
    double *tmp_rand;

    /***** Initialize local variables *****/
    eps  = (double)epsilon / 100.0;
    iNAN = VSL_SS_DNAN;
    dNAN = *((double *)(&iNAN));

    /* Calculate average patterns weight */
    av_weight = 0.0;
    for(i = 0; i < npatt; i++)
    {
        for(j = 0; j < dim; j++)
        {
            av_weight += (double)patt[i*dim + j];
        }
    }
    av_weight /= (npatt * dim);

    threshold = eps / av_weight;

    /***** Allocate memory *****/
    xt       = (double*)malloc( n * dim * sizeof(double) );
    tmp_rand = (double*)malloc( n * sizeof(double) );
    if ( (!tmp_rand) || (!xt) )
    {
        /* Memory allocation failed */
        if(xt!=NULL) free(xt);
        if(tmp_rand!=NULL) free(tmp_rand);
        return -1;
    }

    /* Generating random numbers from multivariate normal distribution */
    errcode = dGenerateGaussianMVData( xt, dim, n, a, C );
    CheckVslError(errcode);
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < dim; j++)
        {
            x[j*n + i] = xt[i*dim + j];
        }
    }

    /* Integration of missing values into matrix of observations */

    /* Generating random numbers from uniform distribution to form
       missing values */
    errcode = dGenerateUniformData( tmp_rand, n, 0.0, (double)npatt );
    CheckVslError(errcode);

    nmv = 0;
    for(i = 0; i < n; i++)
    {
        if (tmp_rand[i] / npatt < threshold)
        {
            /* Here if there are missing values in i-th observation */
            vdFloor(1, &tmp_rand[(i + 1) % n], &floor);
            ipatt = (int)floor;
            miss_vals[i] = ipatt;
            for(j = 0; j < dim; j++)
            {
                if(patt[ipatt*dim + j])
                {
                    x[i + j*n] = dNAN;
                    nmv++;
                }
            }
        }
    }

    *n_miss_vals = nmv;

    free( tmp_rand );
    free( xt );

    return errcode;
}

/*

*/
int sGenerateMissingValuesInput( float *x, MKL_INT dim, MKL_INT n, int epsilon,
                                 int *patt, int npatt, int* miss_vals,
                                 float *a, float *C,
                                 MKL_INT *n_miss_vals )
{
    float av_weight, threshold, eps;
    long  iNAN;
    float sNAN;
    float floor;
    MKL_INT nmv = 0;
    int i, j, ipatt, errcode;
    float *xt;
    float *tmp_rand;

    /***** Initialize local variables *****/
    eps  = (float)epsilon / 100.0;
    iNAN = VSL_SS_SNAN;
    sNAN = *((float *)(&iNAN));

    /* Calculate average patterns weight */
    av_weight = 0.0;
    for(i = 0; i < npatt; i++)
    {
        for(j = 0; j < dim; j++)
        {
            av_weight += (float)patt[i*dim + j];
        }
    }
    av_weight /= (npatt * dim);

    threshold = eps / av_weight;

    /***** Allocate memory *****/
    xt       = (float*)malloc( n * dim * sizeof(float) );
    tmp_rand = (float*)malloc( n * sizeof(float) );
    if ( (!tmp_rand) || (!xt) )
    {
        /* Memory allocation failed */
        if(xt!=NULL) free(xt);
        if(tmp_rand!=NULL) free(tmp_rand);
        return -1;
    }

    /* Generating random numbers from multivariate normal distribution */
    errcode = sGenerateGaussianMVData( xt, dim, n, a, C );
    CheckVslError(errcode);
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < dim; j++)
        {
            x[j*n + i] = xt[i*dim + j];
        }
    }

    /* Integration of missing values into matrix of observations */

    /* Generating random numbers from uniform distribution to form
       missing values */
    errcode = sGenerateUniformData( tmp_rand, n, 0.0, (float)npatt );
    CheckVslError(errcode);

    nmv = 0;
    for(i = 0; i < n; i++)
    {
        if (tmp_rand[i] / npatt < threshold)
        {
            /* Here if there are missing values in i-th observation */
            vsFloor(1, &tmp_rand[(i + 1) % n], &floor);
            ipatt = (int)floor;
            miss_vals[i] = ipatt;
            for(j = 0; j < dim; j++)
            {
                if(patt[ipatt*dim + j])
                {
                    x[i + j*n] = sNAN;
                    nmv++;
                }
            }
        }
    }

    *n_miss_vals = nmv;

    free( tmp_rand );
    free( xt );

    return errcode;
}


/*
//  Generate double precision data set with outliers using
//  GaussianMV RNG with mean a and covariance C for the majotity of values.
//  The outliers are generated with ratio epsilon using GaussianMV RNG with
//  mean a_outl and covariance C.
*/
int dGenerateOutliersInput( double *x, MKL_INT dim, MKL_INT n, double epsilon,
                            double *a, double *a_outl, double *C,
                            double *outlIndex, MKL_INT *cntoutl )
{
    MKL_INT i, k;
    double *x_outl;
    double *tmp_rand;
    MKL_INT co = 0;
    int errcode = 0;

    x_outl   = (double*)malloc( n * dim * sizeof(double) );
    tmp_rand = (double*)malloc( n * sizeof(double) );
    if ( (!x_outl) || (!tmp_rand) )
    {
        if(x_outl!=NULL) free(x_outl);
        if(tmp_rand!=NULL) free(tmp_rand);
        return -1;
    }

    /* Generating random numbers from multivariate normal distribution */
    errcode = dGenerateGaussianMVData( x, dim, n, a, C );
    CheckVslError(errcode);

    /* Generating random numbers from multivariate normal distribution
       to form outliers */
    errcode = dGenerateGaussianMVData( x_outl, dim, n, a_outl, C );
    CheckVslError(errcode);

    /* Generating random numbers from uniform distribution to form outliers */
    errcode = dGenerateUniformData( tmp_rand, n, 0.0, 1.0 );
    CheckVslError(errcode);

    /* Integration of outliers into matrix of observations */
    co = 0;
    for(i = 0; i < n; i++)
    {
        if( tmp_rand[i] < epsilon )
        {
            co++;
            outlIndex[i] = 0.0;
            for(k = 0; k < dim; k++)
            {
                x[i * dim + k] = x_outl[i * dim + k];
            }
        }
        else
        {
            outlIndex[i] = 1.0;
        }
    }

    *cntoutl = co;

    free( tmp_rand );
    free( x_outl );

    return errcode;
}


/*
//  Generate single precision data set with outliers using
//  GaussianMV RNG with mean a and covariance C for the majotity of values.
//  The outliers are generated with ratio epsilon using GaussianMV RNG with
//  mean a_outl and covariance C.
*/
int sGenerateOutliersInput( float *x, MKL_INT dim, MKL_INT n, float epsilon,
                            float *a, float *a_outl, float *C,
                            float *outlIndex, MKL_INT *cntoutl )
{
    MKL_INT i, k;
    float *x_outl;
    float *tmp_rand;
    MKL_INT co = 0;
    int errcode = 0;

    x_outl   = (float*)malloc( n * dim * sizeof(float) );
    tmp_rand = (float*)malloc( n * sizeof(float) );
    if ( (!x_outl) || (!tmp_rand) )
    {
        if(x_outl!=NULL) free(x_outl);
        if(tmp_rand!=NULL) free(tmp_rand);
        return -1;
    }

    /* Generating random numbers from multivariate normal distribution */
    errcode = sGenerateGaussianMVData( x, dim, n, a, C );
    CheckVslError(errcode);

    /* Generating random numbers from multivariate normal distribution
       to form outliers */
    errcode = sGenerateGaussianMVData( x_outl, dim, n, a_outl, C );
    CheckVslError(errcode);

    /* Generating random numbers from uniform distribution to form outliers */
    errcode = sGenerateUniformData( tmp_rand, n, 0.0, 1.0 );
    CheckVslError(errcode);

    /* Integration of outliers into matrix of observations */
    co = 0;
    for(i = 0; i < n; i++)
    {
        if( tmp_rand[i] < epsilon )
        {
            co++;
            outlIndex[i] = 0.0;
            for(k = 0; k < dim; k++)
            {
                x[i * dim + k] = x_outl[i * dim + k];
            }
        }
        else
        {
            outlIndex[i] = 1.0;
        }
    }

    *cntoutl = co;

    free( tmp_rand );
    free( x_outl );

    return errcode;
}
