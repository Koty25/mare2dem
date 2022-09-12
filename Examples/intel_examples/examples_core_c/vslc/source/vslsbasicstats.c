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

/*
!  Content:
!    Calculation of basic statistics Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "statchars.inc"

#define DIM     4        /* Task dimension */
#define N       1000     /* Number of observations */

#define P_THRESHOLD     0.01

float C[DIM][DIM] = {
    { 1.0, 0.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0, 0.0 },
    { 0.0, 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 0.0, 1.0 }
};
float a[DIM] = { 5.0, 5.0, 5.0, 5.0 };

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim;
    MKL_INT n;
    MKL_INT x_storage;
    MKL_INT cov_storage;
    MKL_INT cor_storage;
    float x[N][DIM];
    float cov[DIM][DIM], cor[DIM][DIM];
    float mean[DIM];
    float min_estimate[DIM], max_estimate[DIM];
    float raw2[DIM], raw3[DIM], raw4[DIM];
    float cen2[DIM], cen3[DIM], cen4[DIM];
    float skewness[DIM], kurtosis[DIM], variation[DIM];
    int i, j;
    int errcode;
    unsigned MKL_INT64 estimate = 0;
    int errnums = 0;

    float pval_mean[DIM];
    float pval_cov[DIM][DIM];
    float pval_raw2[DIM], pval_raw3[DIM], pval_raw4[DIM];
    float pval_cen2[DIM], pval_cen3[DIM], pval_cen4[DIM];
    float pval_kurt[DIM], pval_skew[DIM], pval_var[DIM];

    /***** Initializing parameters for Summary Statistics task *****/
    dim         = DIM;
    n           = N;
    x_storage   = VSL_SS_MATRIX_STORAGE_COLS;
    cov_storage = VSL_SS_MATRIX_STORAGE_FULL;
    cor_storage = VSL_SS_MATRIX_STORAGE_FULL;

    for(i = 0; i < dim; i++)
    {
        min_estimate[i] = x[0][i];
        max_estimate[i] = x[0][i];
    }

    /***** Generate data set using VSL Gaussian RNG *****/
    errcode = sGenerateGaussianMVData( (float*)x, dim, n, a, (float*)C );
    CheckVslError(errcode);

    /***** Create Summary Statistics task *****/
    errcode = vslsSSNewTask( &task, &dim, &n, &x_storage, (float*)x, 0, 0 );
    CheckVslError(errcode);

    /***** Edit task parameters for min and max computation *****/
    errcode = vslsSSEditTask( task, VSL_SS_ED_MIN, min_estimate );
    CheckVslError(errcode);

    errcode = vslsSSEditTask( task, VSL_SS_ED_MAX, max_estimate );
    CheckVslError(errcode);

    /***** Edit task parameters for computatin of mean estimate and 2nd, 3rd
           and 4th raw and central moments estimates *****/
    errcode = vslsSSEditMoments( task, mean, raw2, raw3, raw4,
                                 cen2, cen3, cen4 );
    CheckVslError(errcode);

    /***** Edit task parameters for kurtosis, skewness and variation
           computation *****/
    errcode = vslsSSEditTask( task, VSL_SS_ED_KURTOSIS, kurtosis );
    CheckVslError(errcode);

    errcode = vslsSSEditTask( task, VSL_SS_ED_SKEWNESS, skewness );
    CheckVslError(errcode);

    errcode = vslsSSEditTask( task, VSL_SS_ED_VARIATION, variation );
    CheckVslError(errcode);

    /***** Initialization of the task parameters using FULL_STORAGE
           for covariance/correlation matrices *****/
    errcode = vslsSSEditCovCor( task, mean, (float*)cov, &cov_storage ,
                                (float*)cor, &cor_storage );
    CheckVslError(errcode);

    /***** Minimum and maximum are included in the list of estimates
           to compute *****/
    estimate |= VSL_SS_MIN | VSL_SS_MAX;

    /***** Mean and 2nd, 3rd and 4th raw and central moments are included
           in the list of estimates to compute *****/
    estimate |= VSL_SS_MEAN |
        VSL_SS_2R_MOM | VSL_SS_3R_MOM | VSL_SS_4R_MOM |
        VSL_SS_2C_MOM | VSL_SS_3C_MOM |
        VSL_SS_4C_MOM;

    /***** Kurtosis, skewness and variation are included in the list
           of estimates to compute *****/
    estimate |= VSL_SS_KURTOSIS | VSL_SS_SKEWNESS | VSL_SS_VARIATION;

    /***** Covariance and correlation matrices are included in the list
           of estimates to compute *****/
    estimate |= VSL_SS_COV | VSL_SS_COR;

    /***** Compute the estimates using FAST method *****/
    errcode = vslsSSCompute( task, estimate, VSL_SS_METHOD_FAST );
    CheckVslError(errcode);

    /***** Testing stat characteristics of the computed estimates *****/

    /* Comparison of observations with min and max estimates */
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < n; j++)
        {
            if(x[j][i] < min_estimate[i]) errnums++;
            if(x[j][i] > max_estimate[i]) errnums++;
        }
    }

    /* Compute p-values for mean estimates */
    sComputePvalsMean( dim, n, mean, a, (float*)C, pval_mean );
    /* Compute p-values for variance estimates */
    sComputePvalsVariance( dim, n, (float*)cov, cov_storage, (float*)C,
                           (float*)pval_cov );
    /* Compute p-values for covariance estimates */
    sComputePvalsCovariance( dim, n, (float*)cov, cov_storage, (float*)C,
                             (float*)pval_cov );
    /* Compute p-values for raw moments estimates */
    sComputePvalsRawMoments( dim, n, raw2, 2, a, (float*)C, pval_raw2 );
    sComputePvalsRawMoments( dim, n, raw3, 3, a, (float*)C, pval_raw3 );
    sComputePvalsRawMoments( dim, n, raw4, 4, a, (float*)C, pval_raw4 );
    /* Compute p-values for central moments estimates */
    sComputePvalsCentralMoments( dim, n, cen2, 2, a, (float*)C,
                                 pval_cen2 );
    sComputePvalsCentralMoments( dim, n, cen3, 3, a, (float*)C,
                                 pval_cen3 );
    sComputePvalsCentralMoments( dim, n, cen4, 4, a, (float*)C,
                                 pval_cen4 );
    /* Compute p-values for kurtosis, skewness and variation estimates */
    sComputePvalsVariation( dim, n, variation, a, (float*)C, pval_var );
    sComputePvalsSkewness( dim, n, skewness, (float*)C, pval_skew );
    sComputePvalsKurtosis( dim, n, kurtosis, (float*)C, pval_kurt );

    /***** Checking the validity of p-values for all estimates *****/
    for(i = 0; i < dim; i++)
    {
        if (pval_mean[i] < P_THRESHOLD) errnums++;
        if (pval_raw2[i] < P_THRESHOLD) errnums++;
        if (pval_raw3[i] < P_THRESHOLD) errnums++;
        if (pval_raw4[i] < P_THRESHOLD) errnums++;
        if (pval_cen2[i] < P_THRESHOLD) errnums++;
        if (pval_cen3[i] < P_THRESHOLD) errnums++;
        if (pval_cen4[i] < P_THRESHOLD) errnums++;
        if (pval_kurt[i] < P_THRESHOLD) errnums++;
        if (pval_skew[i] < P_THRESHOLD) errnums++;
        if (pval_var[i] < P_THRESHOLD) errnums++;
        if (pval_cov[i][i] < P_THRESHOLD) errnums++;
        for(j = 0; j < i; j++)
        {
            if (pval_cov[i][j] < P_THRESHOLD) errnums++;
        }
    }

    /***** Printing results *****/
    printf("Task dimension : %d\n", (int)dim);
    printf("Number of observations : %d\n\n", (int)n);

    /***** Printing computed minimum, maximum, mean and moments estimates *****/
    printf("               Min        Max        Mean       2nd_raw   ");
    printf("  3rd_raw      4th_raw      2nd_cen    3rd_cen    4th_cen\n");
    for(i = 0; i < dim; i++)
    {
        printf("Variable #%i:  %+lf  %+lf  %+lf  %+lf  %+lf  %+lf", i + 1,
            min_estimate[i], max_estimate[i], mean[i],
            raw2[i], raw3[i], raw4[i]);
        printf("  %+lf  %+lf  %+lf\n", cen2[i], cen3[i], cen4[i]);
    }

    /***** Printing computed kurtosis, skewness and variation estimates *****/
    printf("\n               Kurtosis   Skewness   Variation\n");
    for(i = 0; i < dim; i++)
    {
        printf("Variable #%i:  %+lf  %+lf  %+lf\n", i + 1,
            kurtosis[i], skewness[i], variation[i]);
    }

    /***** Printing computed covariance and correlation matrices *****/
    printf("\n Computed covariance matrix              ");
    printf("   Computed correlation matrix\n");
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+9lf ", cov[i][j]);
        }

        printf("   ");

        for(j = 0; j < dim; j++)
        {
            printf("%+9lf ", cor[i][j]);
        }
        printf("\n");
    }

    /***** Printing p-values for mean and moments estimates *****/
    printf("\n\nP-values of the computed estimates\n\n");
    printf("               Mean       2nd_raw  ");
    printf("  3rd_raw    4th_raw    2nd_cen    3rd_cen    4th_cen\n");
    for(i = 0; i < dim; i++)
    {
        printf("Variable #%i:  %+lf  %+lf  %+lf  %+lf  %+lf  %+lf  %+lf\n",
            i + 1, pval_mean[i], pval_raw2[i], pval_raw3[i], pval_raw4[i],
            pval_cen2[i], pval_cen3[i], pval_cen4[i]);
    }

    /***** Printing p-values for kurtosis, skewness and variation estimates *****/
    printf("\n               Kurtosis   Skewness   Variation\n");
    for(i = 0; i < dim; i++)
    {
        printf("Variable #%i:  %+lf  %+lf  %+lf\n", i + 1,
            pval_kurt[i], pval_skew[i], pval_var[i]);
    }

    /***** Printing p-values for covariance matrix estimate *****/
    printf("\n Covariance matrix\n");
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+9lf ", pval_cov[i][j]);
        }
        printf("\n");
    }

    /***** Printing summary of the test *****/
    if (errnums == 0)
    {
        printf("\n\nAll the computed estimates agree with theory.\n");
    }
    else
    {
        printf("\n\nError: At least one of the computed estimates");
        printf(" disagrees with theory.\n");
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return errcode;
}
