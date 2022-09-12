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

/* file : statchars.inc */
/*
// Define a macro for absolute value computation
*/
#define ABS(x)    (((x) >= 0) ? (x) : (-(x)))

/*
// Computes P-values for double precision mean estimate
*/
void dComputePvalsMean(MKL_INT dim, MKL_INT n, double* mean,
                       double* a, double* C, double* pvals)
{
    double s, sqrt;
    int i;

    for(i = 0; i < dim; i++)
    {
        s = C[i * dim + i] / (double)n;
        vdSqrt( 1, &s, &sqrt );
        s = ABS(mean[i] - a[i]) / sqrt;
        vdCdfNorm( 1, &s, &pvals[i] );
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }
}

/*
// Computes P-values for single precision mean estimate
*/
void sComputePvalsMean(MKL_INT dim, MKL_INT n, float* mean,
                       float* a, float* C, float* pvals)
{
    float s, sqrt;
    int i;

    for(i = 0; i < dim; i++)
    {
        s = C[i * dim + i] / (float)n;
        vsSqrt( 1, &s, &sqrt );
        s = ABS(mean[i] - a[i]) / sqrt;
        vsCdfNorm( 1, &s, &pvals[i] );
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }
}

/*
// Computes P-value for single precision variance estimate.
*/
double dPvalVariance(MKL_INT n, double var, double D)
{
    double s, pval, sqrt;

    s = 2.0 / (double)(n - 1);
    vdSqrt( 1, &s, &sqrt );
    s = ABS((var - D) / (D * sqrt));
    vdCdfNorm( 1, &s, &pval );
    pval = 2.0 * (1.0 - pval);

    return pval;
}

/*
// Computes P-value for single precision variance estimate
*/
float sPvalVariance(MKL_INT n, float var, float D)
{
    float s, pval, sqrt;

    s = 2.0 / (float)(n - 1);
    vsSqrt( 1, &s, &sqrt );
    s = ABS((var - D) / (D * sqrt));
    vsCdfNorm( 1, &s, &pval );
    pval = 2.0 * (1.0 - pval);

    return pval;
}

/*
//  Computes P-value for double precision variance estimates
*/
void dComputePvalsVariance(MKL_INT dim, MKL_INT n, double* cov,
                           MKL_INT cov_storage, double* C, double* pvals)
{
    int i, p_index, c_index;

    for(i = 0; i < dim; i++)
    {
        p_index = i * dim + i;
        if (cov_storage == VSL_SS_MATRIX_STORAGE_FULL)
            c_index = i * dim + i;
        else
            if (cov_storage == VSL_SS_MATRIX_STORAGE_L_PACKED)
                c_index = (i + 1) * (i + 2) / 2 - 1;
            else // cov_storage == VSL_SS_MATRIX_STORAGE_U_PACKED
                c_index = (n * (n + 1) - (n - i) * (n - i + 1)) / 2;
        pvals[p_index] = dPvalVariance( n, cov[c_index], C[p_index] );
    }
}

/*
//  Computes P-value for single precision variance estimates
*/
void sComputePvalsVariance(MKL_INT dim, MKL_INT n, float* cov,
                           MKL_INT cov_storage, float* C, float* pvals)
{
    int i, p_index, c_index;

    for(i = 0; i < dim; i++)
    {
        p_index = i * dim + i;
        if (cov_storage == VSL_SS_MATRIX_STORAGE_FULL)
            c_index = i * dim + i;
        else
            if (cov_storage == VSL_SS_MATRIX_STORAGE_L_PACKED)
                c_index = (i + 1) * (i + 2) / 2 - 1;
            else // cov_storage == VSL_SS_MATRIX_STORAGE_U_PACKED
                c_index = (n * (n + 1) - (n - i) * (n - i + 1)) / 2;
        pvals[p_index] = sPvalVariance( n, cov[c_index], C[p_index] );
    }
}

/*
//  Computes P-value for double precision covariance estimate
*/
double dPvalCovariance(MKL_INT n, double cov, double c_ii, double c_ij,
                       double c_jj)
{
    double s, pval, sqrt;

    s = (double)(n - 1) / ((c_ii * c_jj) + (c_ij * c_ij)) ;
    vdSqrt(1, &s, &sqrt);
    s = ABS(cov - c_ij) * sqrt;
    vdCdfNorm(1, &s, &pval);
    pval = 2.0 * (1.0 - pval);

    return pval;
}

/*
//  Computes P-value for single precision covariance estimate
*/
float sPvalCovariance(MKL_INT n, float cov, float c_ii, float c_ij,
                      float c_jj)
{
    float s, pval, sqrt;

    s = (float)(n - 1) / ((c_ii * c_jj) + (c_ij * c_ij)) ;
    vsSqrt(1, &s, &sqrt);
    s = ABS(cov - c_ij) * sqrt;
    vsCdfNorm(1, &s, &pval);
    pval = 2.0 * (1.0 - pval);

    return pval;
}

/*
//  Compute p-values for double precision covariance matrix estimate
*/
void dComputePvalsCovariance(MKL_INT dim, MKL_INT n, double* cov,
                             MKL_INT cov_storage, double* C, double* pvals)
{
    int i, j, c_index;

    for(i = 0; i < dim; i++)
    {
        if (cov_storage == VSL_SS_MATRIX_STORAGE_FULL)
            c_index = i * dim;
        else
            if (cov_storage == VSL_SS_MATRIX_STORAGE_L_PACKED)
                c_index = i * (i + 1) / 2;
            else // cov_storage == VSL_SS_MATRIX_STORAGE_U_PACKED
                c_index = (n * (n + 1) - (n - i) * (n - i + 1)) / 2;
        for(j = 0; j < i; j++)
        {
            pvals[i * dim + j] = pvals[j * dim + i] =
                dPvalCovariance(n, cov[c_index++],
                C[i * dim + i], C[i * dim + j], C[j * dim + j]);
        }
    }
}

/*
//  Compute p-values for single precision covariance matrix estimate
*/
void sComputePvalsCovariance(MKL_INT dim, MKL_INT n, float* cov,
                             MKL_INT cov_storage, float* C, float* pvals)
{
    int i, j, c_index;

    for(i = 0; i < dim; i++)
    {
        if (cov_storage == VSL_SS_MATRIX_STORAGE_FULL)
            c_index = i * dim;
        else
            if (cov_storage == VSL_SS_MATRIX_STORAGE_L_PACKED)
                c_index = i * (i + 1) / 2;
            else // cov_storage == VSL_SS_MATRIX_STORAGE_U_PACKED
                c_index = (n * (n + 1) - (n - i) * (n - i + 1)) / 2;
        for(j = 0; j < i; j++)
        {
            pvals[i * dim + j] = pvals[j * dim + i] =
                sPvalCovariance(n, cov[c_index++],
                C[i * dim + i], C[i * dim + j], C[j * dim + j]);
        }
    }
}

/* Compute double precision raw moment estimate */
double dComputeRawMoment(int order, double mean, double sigma);

/*
// Compute p-values for double precision raw moments estimates
*/
void dComputePvalsRawMoments(MKL_INT dim, MKL_INT n, double* raw, int order,
                             double* a, double* C, double* pvals)
{
    int i;
    double sigma, Es, Es2, D, sqrtD;
    double *d;

    d = (double*)malloc( dim * sizeof(double) );
    if(d==NULL) exit(1);

    for(i = 0; i < dim; i++)
    {
        vdSqrt( 1, &C[i * dim + i], &sigma );

        Es = dComputeRawMoment( order, a[i], sigma );
        Es2 = dComputeRawMoment( 2 * order, a[i], sigma );
        D = (Es2 - Es * Es) / (double)n;
        vdSqrt( 1, &D, &sqrtD );
        d[i] = ABS(raw[i] - Es) / sqrtD;
    }

    vdCdfNorm( dim, d, pvals );

    for(i = 0; i < dim; i++)
    {
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }

    free( d );
}

/* Compute single precision raw moment estimate */
float sComputeRawMoment(int order, float mean, float sigma);

/*
// Compute p-values for single precision raw moments estimates
*/
void sComputePvalsRawMoments(MKL_INT dim, MKL_INT n, float* raw, int order,
                             float* a, float* C, float* pvals)
{
    int i;
    float sigma, Es, Es2, D, sqrtD;
    float *d;

    d = (float*)malloc( dim * sizeof(float) );
    if(d==NULL) exit(1);

    for(i = 0; i < dim; i++)
    {
        vsSqrt( 1, &C[i * dim + i], &sigma );

        Es = sComputeRawMoment( order, a[i], sigma );
        Es2 = sComputeRawMoment( 2 * order, a[i], sigma );
        D = (Es2 - Es * Es) / (float)n;
        vsSqrt( 1, &D, &sqrtD );
        d[i] = ABS(raw[i] - Es) / sqrtD;
    }

    vsCdfNorm( dim, d, pvals );

    for(i = 0; i < dim; i++)
    {
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }

    free( d );
}

/* Compute odd (repeated) factorial */
int oddFactorial(int n);

/* Compute theoretical central moment estrimate */
double dComputeCentralMoment(int order, double D);

/*
// Compute p-values for double precision central moments estimates
*/
void dComputePvalsCentralMoments(MKL_INT dim, MKL_INT n, double* central,
                                 int order, double* a, double* C,
                                 double* pvals)
{
    int i;
    double D, sqrtD;
    double *d;
    double mu_2, mu_2k, mu_kp1, mu_km1, mu_k;

    d = (double*)malloc( dim * sizeof(double) );
    if(d==NULL) exit(1);

    for (i = 0; i < dim; i++)
    {
        mu_2 = C[i * dim + i];
        mu_2k  = dComputeCentralMoment( 2*order,   mu_2 );
        mu_kp1 = dComputeCentralMoment( order + 1, mu_2 );
        mu_km1 = dComputeCentralMoment( order - 1, mu_2 );
        mu_k   = dComputeCentralMoment( order,     mu_2 );
        D = mu_2k - 2.0 * mu_kp1 * mu_km1 - mu_k * mu_k +
            order * order * mu_2 * mu_km1;
        D /= (double)n;
        vdSqrt( 1, &D, &sqrtD );
        d[i] = ABS(central[i] - mu_k) / sqrtD;
    }

    vdCdfNorm( dim, d, pvals );

    for(i = 0; i < dim; i++)
    {
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }

    free( d );
}

/* Compute theoretical central moment estrimate */
float sComputeCentralMoment(int order, float D);

/*
// Compute p-values for single precision central moments estimates
*/
void sComputePvalsCentralMoments(MKL_INT dim, MKL_INT n, float* central,
                                 int order, float* a, float* C,
                                 float* pvals)
{
    int i;
    float D, sqrtD;
    float *d;
    float mu_2, mu_2k, mu_kp1, mu_km1, mu_k;

    d = (float*)malloc( dim * sizeof(float) );
    if(d==NULL) exit(1);

    for (i = 0; i < dim; i++)
    {
        mu_2 = C[i * dim + i];
        mu_2k  = sComputeCentralMoment( 2*order,   mu_2 );
        mu_kp1 = sComputeCentralMoment( order + 1, mu_2 );
        mu_km1 = sComputeCentralMoment( order - 1, mu_2 );
        mu_k   = sComputeCentralMoment( order,     mu_2 );
        D = mu_2k - 2.0 * mu_kp1 * mu_km1 - mu_k * mu_k +
            order * order * mu_2 * mu_km1;
        D /= (float)n;
        vsSqrt( 1, &D, &sqrtD );
        d[i] = ABS(central[i] - mu_k) / sqrtD;
    }

    vsCdfNorm( dim, d, pvals );

    for(i = 0; i < dim; i++)
    {
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }

    free( d );
}


/*
//  Compute p-values for double precision kurtosis estimates
*/
void dComputePvalsKurtosis(MKL_INT dim, MKL_INT n, double* kurt,
                           double* C, double* pvals)
{
    int i;
    double sigma2, Es, D, sqrtD;
    double mu2, mu4, mu6, mu8, mu2_2, mu4_2;
    double *d;

    d = (double*)malloc( dim * sizeof(double) );
    if(d==NULL) exit(1);

    for (i = 0; i < dim; i++)
    {
        sigma2 = C[i * dim + i];
        Es = dComputeCentralMoment( 4, sigma2 );
        Es /= (sigma2 * sigma2);
        Es -= 3.0;

        mu2 = dComputeCentralMoment( 2, sigma2 );
        mu4 = dComputeCentralMoment( 4, sigma2 );
        mu6 = dComputeCentralMoment( 6, sigma2 );
        mu8 = dComputeCentralMoment( 8, sigma2 );
        mu2_2 = mu2 * mu2;
        mu4_2 = mu4 * mu4;

        D = (mu2_2 * mu8 - 4.0 * mu2 * mu4 * mu6 + 4.0 * mu4 * mu4_2 -
             mu2_2 * mu4_2) / ((double)n * mu2_2 * mu2_2 * mu2_2);

        vdSqrt( 1, &D, &sqrtD );

        d[i] = ABS(kurt[i] - Es) / sqrtD;
    }

    vdCdfNorm( dim, d, pvals );

    for(i = 0; i < dim; i++)
    {
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }

    free( d );
}

/*
//  Compute p-values for single precision kurtosis estimates
*/
void sComputePvalsKurtosis(MKL_INT dim, MKL_INT n, float* kurt,
                           float* C, float* pvals)
{
    int i;
    float sigma2, Es, D, sqrtD;
    float mu2, mu4, mu6, mu8, mu2_2, mu4_2;
    float *d;

    d = (float*)malloc( dim * sizeof(float) );
    if(d==NULL) exit(1);

    for (i = 0; i < dim; i++)
    {
        sigma2 = C[i * dim + i];
        Es = sComputeCentralMoment( 4, sigma2 );
        Es /= (sigma2 * sigma2);
        Es -= 3.0;

        mu2 = sComputeCentralMoment( 2, sigma2 );
        mu4 = sComputeCentralMoment( 4, sigma2 );
        mu6 = sComputeCentralMoment( 6, sigma2 );
        mu8 = sComputeCentralMoment( 8, sigma2 );
        mu2_2 = mu2 * mu2;
        mu4_2 = mu4 * mu4;

        D = (mu2_2 * mu8 - 4.0 * mu2 * mu4 * mu6 + 4.0 * mu4 * mu4_2 -
             mu2_2 * mu4_2) / ((float)n * mu2_2 * mu2_2 * mu2_2);

        vsSqrt( 1, &D, &sqrtD );

        d[i] = ABS(kurt[i] - Es) / sqrtD;
    }

    vsCdfNorm( dim, d, pvals );

    for(i = 0; i < dim; i++)
    {
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }

    free( d );
}

/*
//  Compute p-values for double precision skewness estimates
*/
void dComputePvalsSkewness(MKL_INT dim, MKL_INT n, double* skew,
                           double* C, double* pvals)
{
    int i;
    double sigma2, D, sqrtD;
    double mu2, mu4, mu6, mu2_2, mu2_3, mu4_2;
    double *d;

    d = (double*)malloc( dim * sizeof(double) );
    if(d==NULL) exit(1);

    for (i = 0; i < dim; i++)
    {
        sigma2 = C[i * dim + i];
        mu2 = dComputeCentralMoment( 2, sigma2 );
        mu4 = dComputeCentralMoment( 4, sigma2 );
        mu6 = dComputeCentralMoment( 6, sigma2 );
        mu2_2 = mu2 * mu2;
        mu2_3 = mu2 * mu2_2;
        mu4_2 = mu4 * mu4;

        D = (mu2_2 *mu6 - 6.0 * mu2_3 * mu4 + 9.0 * mu2_2 * mu2_3) /
            ((double)n * mu2_2 * mu2_3);

        vdSqrt( 1, &D, &sqrtD );

        d[i] = ABS(skew[i]) / sqrtD;
    }

    vdCdfNorm( dim, d, pvals );

    for(i = 0; i < dim; i++)
    {
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }

    free( d );
}

/*
//  Compute p-values for single precision skewness estimates
*/
void sComputePvalsSkewness(MKL_INT dim, MKL_INT n, float* skew,
                           float* C, float* pvals)
{
    int i;
    float sigma2, D, sqrtD;
    float mu2, mu4, mu6, mu2_2, mu2_3, mu4_2;
    float *d;

    d = (float*)malloc( dim * sizeof(float) );
    if(d==NULL) exit(1);

    for (i = 0; i < dim; i++)
    {
        sigma2 = C[i * dim + i];
        mu2 = sComputeCentralMoment( 2, sigma2 );
        mu4 = sComputeCentralMoment( 4, sigma2 );
        mu6 = sComputeCentralMoment( 6, sigma2 );
        mu2_2 = mu2 * mu2;
        mu2_3 = mu2 * mu2_2;
        mu4_2 = mu4 * mu4;

        D = (mu2_2 *mu6 - 6.0 * mu2_3 * mu4 + 9.0 * mu2_2 * mu2_3) /
            ((float)n * mu2_2 * mu2_3);

        vsSqrt( 1, &D, &sqrtD );

        d[i] = ABS(skew[i]) / sqrtD;
    }

    vsCdfNorm( dim, d, pvals );

    for(i = 0; i < dim; i++)
    {
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }

    free( d );
}

/*
//  Compute p-values for double precision variation estimates
*/
void dComputePvalsVariation(MKL_INT dim, MKL_INT n, double* variation,
                            double* a, double* C, double* pvals)
{
    int i;
    double s, sqrt, sigma, V;
    double *d;

    d = (double*)malloc( dim * sizeof(double) );
    if(d==NULL) exit(1);

    for(i = 0; i < dim; i++)
    {
        s = 2.0 / (double)(n - 1);
        vdSqrt( 1, &s, &sqrt );
        vdSqrt( 1, &C[i * dim + i], &sigma );
        V = sigma / a[i];
        d[i] = ABS((variation[i] - V) / (V * sqrt));
    }

    vdCdfNorm( dim, d, pvals );

    for(i = 0; i < dim; i++)
    {
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }

    free( d );
}

/*
//  Compute p-values for single precision variation estimates
*/
void sComputePvalsVariation(MKL_INT dim, MKL_INT n, float* variation,
                            float* a, float* C, float* pvals)
{
    int i;
    float s, sqrt, sigma, V;
    float *d;

    d = (float*)malloc( dim * sizeof(float) );
    if(d==NULL) exit(1);

    for(i = 0; i < dim; i++)
    {
        s = 2.0 / (float)(n - 1);
        vsSqrt( 1, &s, &sqrt );
        vsSqrt( 1, &C[i * dim + i], &sigma );
        V = sigma / a[i];
        d[i] = ABS((variation[i] - V) / (V * sqrt));
    }

    vsCdfNorm( dim, d, pvals );

    for(i = 0; i < dim; i++)
    {
        pvals[i] = 2.0 * (1.0 - pvals[i]);
    }

    free( d );
}

/*
// Implementation of the routine for computation of double precision
// raw moments estimates
*/
double dComputeRawMoment(int order, double mean, double sigma)
{
    double mu=0.0, tmp;
    double mean2, mean4, mean6, mean8;
    double sigma2, sigma4, sigma6, sigma8;

    switch (order)
    {
    case 1:
        mu = mean;
        break;
    case 2:
        vdSqr( 1, &mean, &mean2 );
        vdSqr( 1, &sigma, &sigma2 );
        mu = mean2 + sigma2;
        break;
    case 3:
        vdSqr( 1, &mean, &mean2 );
        vdSqr( 1, &sigma, &sigma2 );
        mu = mean * (mean2 + 3.0 * sigma2);
        break;
    case 4:
        vdSqr( 1, &mean, &mean2 );
        vdSqr( 1, &sigma, &sigma2 );
        vdPowx( 1, &mean, 4.0, &mean4 );
        vdPowx( 1, &sigma, 4.0, &sigma4 );
        mu = mean4 + 6.0 * mean2 * sigma2 + 3.0 * sigma4;
        break;
    case 6:
        vdSqr( 1, &mean, &mean2 );
        vdSqr( 1, &sigma, &sigma2 );
        vdPowx( 1, &mean, 6.0, &mean6 );
        vdPowx( 1, &sigma, 6.0, &sigma6 );
        mu = mean6 + 15.0 * mean2 * sigma2 * (mean2 + 3.0 * sigma2) +
             15.0 * sigma6;
        break;
    case 8:
        vdSqr( 1, &mean, &mean2 );
        vdSqr( 1, &sigma, &sigma2 );
        vdPowx( 1, &mean, 4.0, &mean4 );
        vdPowx( 1, &sigma, 4.0, &sigma4 );
        vdPowx( 1, &mean, 8.0, &mean8 );
        vdPowx( 1, &sigma, 8.0, &sigma8 );
        tmp = mean2 * sigma2;
        mu = mean8 + 14.0 * tmp * (2.0 * mean4 + 15.0 * tmp + 30.0 * sigma4) +
             105.0 * sigma8;
        break;
    default:
        break;
    }
    return mu;
}

/*
// Implementation of the routine for computation of single precision
// raw moments estimates
*/
float sComputeRawMoment(int order, float mean, float sigma)
{
    float mu=0.0, tmp;
    float mean2, mean4, mean6, mean8;
    float sigma2, sigma4, sigma6, sigma8;

    switch (order)
    {
    case 1:
        mu = mean;
        break;
    case 2:
        vsSqr( 1, &mean, &mean2 );
        vsSqr( 1, &sigma, &sigma2 );
        mu = mean2 + sigma2;
        break;
    case 3:
        vsSqr( 1, &mean, &mean2 );
        vsSqr( 1, &sigma, &sigma2 );
        mu = mean * (mean2 + 3.0 * sigma2);
        break;
    case 4:
        vsSqr( 1, &mean, &mean2 );
        vsSqr( 1, &sigma, &sigma2 );
        vsPowx( 1, &mean, 4.0, &mean4 );
        vsPowx( 1, &sigma, 4.0, &sigma4 );
        mu = mean4 + 6.0 * mean2 * sigma2 + 3.0 * sigma4;
        break;
    case 6:
        vsSqr( 1, &mean, &mean2 );
        vsSqr( 1, &sigma, &sigma2 );
        vsPowx( 1, &mean, 6.0, &mean6 );
        vsPowx( 1, &sigma, 6.0, &sigma6 );
        mu = mean6 + 15.0 * mean2 * sigma2 * (mean2 + 3.0 * sigma2) +
             15.0 * sigma6;
        break;
    case 8:
        vsSqr( 1, &mean, &mean2 );
        vsSqr( 1, &sigma, &sigma2 );
        vsPowx( 1, &mean, 4.0, &mean4 );
        vsPowx( 1, &sigma, 4.0, &sigma4 );
        vsPowx( 1, &mean, 8.0, &mean8 );
        vsPowx( 1, &sigma, 8.0, &sigma8 );
        tmp = mean2 * sigma2;
        mu = mean8 + 14.0 * tmp * (2.0 * mean4 + 15.0 * tmp + 30.0 * sigma4) +
             105.0 * sigma8;
        break;
    default:
        break;
    }
    return mu;
}

/*
// Implementation of the routine for odd (repeated) factorial computation
// starts here
*/
int oddFactorial(int n)
{
    if (n == 1)
        return 1;
    else
        return (n * oddFactorial(n - 2));
}

/*
// Implementation of the routine for computation of double precision
// central moments estimates
*/
double dComputeCentralMoment(int order, double D)
{
    int factorial;
    double tmp;
    double mu = 0.0;

    if (order % 2 == 0)
    {
        vdPowx( 1, &D, (double)order / 2.0, &tmp );
        factorial = oddFactorial( order - 1 );
        mu = (double)factorial * tmp;
    }

    return mu;
}

/*
// Implementation of the routine for computation of single precision
// central moments estimates
*/
float sComputeCentralMoment(int order, float D)
{
    int factorial;
    float tmp;
    float mu = 0.0;

    if (order % 2 == 0)
    {
        vsPowx( 1, &D, (float)order / 2.0, &tmp );
        factorial = oddFactorial( order - 1 );
        mu = (float)factorial * tmp;
    }

    return mu;
}

/*
// Implementation of double precision sweep operatior
//
// Reference: J. L. Schafer, Analysis of Incomplete Multivariate Data, 1999
*/
void dSweep(int k, int dim, double* g)
{
    int i, j;
    double h;

    g[k*dim + k] = -1.0 / g[k*dim + k];

    h = g[k*dim + k];
    for(i = 0; i < k; i++)
    {
        g[i*dim + k] *= h;
        g[i*dim + k] = -g[i*dim + k];
        g[k*dim + i] = g[i*dim + k];
    }

    for(i = k + 1; i < dim; i++)
    {
        g[i*dim + k] *= h;
        g[i*dim + k] = -g[i*dim + k];
        g[k*dim + i] = g[i*dim + k];
    }

    for(i = 0; i < dim; i++)
    {
        if (i != k)
        {
            for(j = 0; j < i; j++)
            {
                if (j != k)
                {
                    g[i*dim + j] += (g[k*dim + i] * g[j*dim + k] / h);
                    g[j*dim + i] = g[i*dim + j];
                }
            }
            g[i*dim + i] += (g[k*dim + i] * g[i*dim + k] / h);
        }
    }
}

/*
// Implementation of single precision sweep operatior
//
// Reference: J. L. Schafer, Analysis of Incomplete Multivariate Data, 1999
*/
void sSweep(int k, int dim, float* g)
{
    int i, j;
    float h;

    g[k*dim + k] = -1.0 / g[k*dim + k];

    h = g[k*dim + k];
    for(i = 0; i < k; i++)
    {
        g[i*dim + k] *= h;
        g[i*dim + k] = -g[i*dim + k];
        g[k*dim + i] = g[i*dim + k];
    }

    for(i = k + 1; i < dim; i++)
    {
        g[i*dim + k] *= h;
        g[i*dim + k] = -g[i*dim + k];
        g[k*dim + i] = g[i*dim + k];
    }

    for(i = 0; i < dim; i++)
    {
        if (i != k)
        {
            for(j = 0; j < i; j++)
            {
                if (j != k)
                {
                    g[i*dim + j] += (g[k*dim + i] * g[j*dim + k] / h);
                    g[j*dim + i] = g[i*dim + j];
                }
            }
            g[i*dim + i] += (g[k*dim + i] * g[i*dim + k] / h);
        }
    }
}
