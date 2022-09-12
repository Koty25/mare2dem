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
!    Computation of robust covariance matrix and mean Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "statchars.inc"

#define DIM         5       /* Task dimension */
#define N           5000   /* Number of observations */

#define P_THRESHOLD 0.001

#define RATIO       2       /* Ratio of outliers in the dataset */
#define M           100.0   /* Mean of the outliers */
#define COEFF       1.0     /* Coefficient to compute covarince of outliers */

/***** Robust method parameters *****/
#define BD_POINT    0.4
#define ARP         0.001
#define ACCURACY    0.001
#define ITER_NUM    5

/***** Parameters for major distribution *****/
double C[DIM][DIM] = {
    { 1.0, 0.1, 0.1, 0.1, 0.1 },
    { 0.1, 2.0, 0.1, 0.1, 0.1 },
    { 0.1, 0.1, 1.0, 0.1, 0.1 },
    { 0.1, 0.1, 0.1, 2.0, 0.1 },
    { 0.1, 0.1, 0.1, 0.1, 1.0 }
};

double a[DIM] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim;
    MKL_INT n;
    MKL_INT x_storage;
    MKL_INT cov_storage;
    MKL_INT rcov_storage;
    double x[DIM][N];
    double mean[DIM], rmean[DIM];
    double cov[DIM*(DIM+1)/2], rcov[DIM*(DIM+1)/2];
    double pval_c[DIM][DIM], pval_r[DIM][DIM];
    int i, j, k, errcode;
    int errnums = 0;

    MKL_INT robust_params_n;
    double robust_method_params[VSL_SS_TBS_PARAMS_N];

    /***** Initializing parameters for Summary Statistics task *****/
    dim              = DIM;
    n                = N;
    x_storage        = VSL_SS_MATRIX_STORAGE_ROWS;
    cov_storage      = VSL_SS_MATRIX_STORAGE_L_PACKED;
    rcov_storage     = VSL_SS_MATRIX_STORAGE_L_PACKED;

    /***** Generate dataset *****/
    errcode = dGenerateContaminatedDataset( (double*)x, dim, n, a,
                                            (double*)C, RATIO, M, COEFF );
    CheckVslError(errcode);

    /***** Create Summary Statistics task *****/
    errcode = vsldSSNewTask( &task, &dim, &n, &x_storage, (double*)x, 0, 0 );
    CheckVslError(errcode);

    errcode = vsldSSEditCovCor( task, mean, cov, &cov_storage, 0, 0 );
    CheckVslError(errcode);

    /***** Initialization of the task parameters
           for robust covariance estimator *****/
    robust_params_n         = VSL_SS_TBS_PARAMS_N;
    robust_method_params[0] = BD_POINT;
    robust_method_params[1] = ARP;
    robust_method_params[2] = ACCURACY;
    robust_method_params[3] = ITER_NUM;

    errcode =  vsldSSEditRobustCovariance( task, &rcov_storage,
                                           &robust_params_n,
                                           robust_method_params,
                                           rmean, rcov );
    CheckVslError(errcode);

    /***** Compute covariance matrix using FAST method  *****/
    errcode = vsldSSCompute( task, VSL_SS_COV |
                             VSL_SS_ROBUST_COV,
                             VSL_SS_METHOD_FAST | VSL_SS_METHOD_TBS );
    CheckVslError(errcode);

    /***** Printing results *****/
    printf("Task dimension : %d\n", (int)dim);
    printf("Number of observations : %d\n\n", (int)n);

    /***** Print original mean and covariance matrix *****/
    printf("Original covariance matrix\n");
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%lf ", C[i][j]);
        }
        printf("\n");
    }

    printf("\nOriginal vector of means\n");
    for(i = 0; i < dim; i++)
    {
        printf("%lf, ", a[i]);
    }
    printf("\n\n");

    /***** Print classical mean and covariance matrix estimate *****/
    printf("Classical covariance estimate\n");
    k = 0;
    for ( i = 0; i < dim; i++ )
    {
        for(j = 0; j <= i; j++)
        {
            printf("%lf ", cov[k++]);
        }
        printf("\n");
    }

    printf("\nClassical mean estimate\n");
    for (i = 0; i < dim; i++)
    {
        printf("%lf, ", mean[i]);
    }
    printf("\n\n");

    /***** Print robust mean and covariance matrix estimate *****/
    printf("Robust covariance estimate:\n");
    k = 0;
    for (i = 0; i < dim; i++)
    {
        for(j = 0; j <= i; j++)
        {
            printf("%lf ", rcov[k++]);
        }
        printf("\n");
    }

    printf("\nRobust mean estimate:\n");
    for(i = 0; i < dim; i++)
    {
        printf("%lf, ", rmean[i]);
    }
    printf("\n");

    /***** Testing stat characteristics of classic and robust
           covariance matrices *****/
    dComputePvalsVariance( dim, n,  cov, cov_storage, (double*)C,
                           (double*)pval_c );
    dComputePvalsVariance( dim, n, rcov, rcov_storage, (double*)C,
                           (double*)pval_r );
    dComputePvalsCovariance( dim, n,  cov, cov_storage, (double*)C,
                            (double*)pval_c );
    dComputePvalsCovariance( dim, n, rcov, rcov_storage, (double*)C,
                             (double*)pval_r );

    for(i = 0; i < dim; i++)
    {
        for(j = 0; j <= i; j++)
        {
            if (pval_r[i][j] < P_THRESHOLD) errnums++;
        }
    }

    printf("\n\nP-values of the computed classic covariance matrix\n");

    for(i = 0; i < dim; i++)
    {
        for(j = 0; j <= i; j++)
        {
            printf("%9lf ", pval_c[i][j]);
        }
        printf("\n");
    }

    printf("\n\nP-values of the computed robust covariance matrix\n");

    for(i = 0; i < dim; i++)
    {
        for(j = 0; j <= i; j++)
        {
            printf("%9lf ", pval_r[i][j]);
        }
        printf("\n");
    }

    /***** Printing summary of the test *****/
    if (errnums == 0)
    {
        printf("\n\nRobust covariance estimate agrees with theory\n");
    }
    else
    {
        printf("\n\nError: Robust covariance estimate");
        printf(" disagrees with theory\n");
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return 0;
}
