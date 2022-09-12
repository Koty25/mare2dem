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
!    Calculation of group/pooled covariance matrices Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "statchars.inc"

#define DIM     3      /* Task dimension */
#define N       10000  /* Number of observations */
#define G       2      /* Number of groups */
#define GN      2      /* Number of group covariance matrices */

#define P_THRESHOLD    0.005

double C[DIM][DIM] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
};

double m[DIM] = { 0.0, 0.0, 0.0 };

/***** 1st and 2nd group matrix to be returned *****/
MKL_INT group_matrix_indices[G] = { 1, 1 };

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim;
    MKL_INT n;
    MKL_INT x_storage;
    MKL_INT cov_storage;
    MKL_INT pld_cov_storage;
    MKL_INT grp_cov_storage;
    double x[DIM][N];
    double mean[DIM], pld_mean[DIM], grp_mean[DIM*GN];
    double cov[DIM*DIM], pld_cov[DIM*DIM], grp_cov[DIM*DIM*GN];
    double a = 0.0, sigma = 1.0;
    MKL_INT group_indices[N];
    int i, j, errcode;
    int errnums = 0;

    double pval_pld_mean[DIM], pval_grp_mean[DIM*GN];
    double pval_pld_cov[DIM*DIM], pval_grp_cov[DIM*DIM*GN];

    /***** Initializing parameters for Summary Statistics task *****/
    dim              = DIM;
    n                = N;
    x_storage        = VSL_SS_MATRIX_STORAGE_ROWS;
    cov_storage      = VSL_SS_MATRIX_STORAGE_FULL;
    pld_cov_storage  = VSL_SS_MATRIX_STORAGE_FULL;
    grp_cov_storage  = VSL_SS_MATRIX_STORAGE_FULL;

    /***** Generate data set using VSL Gaussian RNG with mean a = 0 and
           stdev = 1 *****/
    errcode = dGenerateGaussianData( (double*)x, dim, n, a, sigma );
    CheckVslError(errcode);

    /***** Create Summary Statistics task *****/
    errcode = vsldSSNewTask( &task, &dim, &n, &x_storage, (double*)x, 0, 0 );
    CheckVslError(errcode);

    /***** Initialization of the task parameters for pooled covariance
           estimator *****/
    errcode = vsliSSEditTask( task, VSL_SS_ED_POOLED_COV_STORAGE,
        &pld_cov_storage );
    CheckVslError(errcode);

    /***** Initialization of the task parameters for group covariance
           estimator *****/
    errcode = vsliSSEditTask( task, VSL_SS_ED_GROUP_COV_STORAGE,
        &grp_cov_storage );
    CheckVslError(errcode);

    /***** Initialization of the task parameters using FULL_STORAGE
        for covariance/correlation matrices *****/
    errcode = vsldSSEditCovCor( task, mean, cov, &cov_storage,
                                NULL, NULL );
    CheckVslError(errcode);

    /***** Dividing elements into odd and even *****/
    for(i = 0; i < n; i++)
    {
        group_indices[i] = i % 2;
    }

    /***** Initialization of the task parameters for pooled and
           group covariance estimators *****/
    errcode = vsldSSEditPooledCovariance( task, group_indices,
        pld_mean, pld_cov, group_matrix_indices, grp_mean, grp_cov );
    CheckVslError(errcode);

    /***** Compute covariance matrices using FAST method  *****/
    errcode = vsldSSCompute( task,
        VSL_SS_COV |
        VSL_SS_POOLED_COV | VSL_SS_GROUP_COV,
        VSL_SS_METHOD_1PASS );
    CheckVslError(errcode);

    /***** Testing stat characteristics of mean and covariance matrices *****/
    /* Compute p-values for group mean estimates */
    dComputePvalsMean( dim, n, grp_mean, m, (double*)C, pval_grp_mean );
    dComputePvalsMean( dim, n, &grp_mean[dim], m, (double*)C,
                       &pval_grp_mean[dim] );
    /* Compute p-values for group variance estimates */
    dComputePvalsVariance( dim, n, grp_cov, grp_cov_storage,
                           (double*)C, pval_grp_cov );
    dComputePvalsVariance( dim, n, &grp_cov[dim * dim], grp_cov_storage,
                           (double*)C, &pval_grp_cov[dim * dim] );
    /* Compute p-values for group covariance estimates */
    dComputePvalsCovariance( dim, n, grp_cov, grp_cov_storage,
                             (double*)C, pval_grp_cov );
    dComputePvalsCovariance( dim, n, &grp_cov[dim * dim], grp_cov_storage,
                             (double*)C, &pval_grp_cov[dim * dim] );
    /* Compute p-values for pooled mean estimates */
    dComputePvalsMean( dim, n, pld_mean, m, (double*)C, pval_pld_mean );
    /* Compute p-values for pooled variance estimates */
    dComputePvalsVariance( dim, n, pld_cov, pld_cov_storage, (double*)C,
                           pval_pld_cov );
    /* Compute p-values for pooled covariance estimates */
    dComputePvalsCovariance( dim, n, pld_cov, pld_cov_storage, (double*)C,
                             pval_pld_cov );

    for(i = 0; i < dim; i++)
    {
        if (pval_grp_mean[i] < P_THRESHOLD) errnums++;
        if (pval_grp_mean[i + dim] < P_THRESHOLD) errnums++;
        if (pval_pld_mean[i] < P_THRESHOLD) errnums++;
        for(j = 0; j <= i; j++)
        {
            if (pval_grp_cov[i * dim + j] < P_THRESHOLD) errnums++;
            if (pval_grp_cov[i * dim + j + dim * dim] < P_THRESHOLD) errnums++;
            if (pval_pld_cov[i * dim + j] < P_THRESHOLD) errnums++;
        }
    }

    /***** Printing results *****/
    printf("Task dimension : %d\n", (int)dim);
    printf("Number of observations : %d\n", (int)n);

    /***** Print exact covariance matrix and mean *****/
    printf("\nExact covariance matrix:\n");
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+lf ", C[i][j]);
        }
        printf("\n");
    }

    printf("\nExact means:\n");
    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", m[i]);
    }

    /***** Print computed covariance matrix and mean estimates *****/
    printf("\n\nComputed covariance matrix:\n");
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            printf("%+lf ", cov[i * dim + j]);
        }
        printf("\n");
    }

    printf("\nComputed means:\n");
    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", mean[i]);
    }

    /***** Print group covariance matrices and mean estimates *****/
    printf("\n\nGroup covariance matrices:\n");
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+lf ", grp_cov[i * dim + j]);
        }

        printf("     ");

        for(j = 0; j < dim; j++)
        {
            printf("%+lf ", grp_cov[i * dim + j + dim * dim]);
        }
        printf("\n");
    }

    printf("\nGroup means:\n");
    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", grp_mean[i]);
    }

    printf("     ");

    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", grp_mean[i + dim]);
    }

    /***** Print pooled covariance matrix and mean estimates *****/
    printf("\n\nPooled covariance matrix:\n");
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+lf ", pld_cov[i * dim + j]);
        }
        printf("\n");
    }

    printf("\nPooled means:\n");
    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", pld_mean[i]);
    }

    /***** Print P-values of the group covariance matrices *****/
    printf("\n\n\nP-values of the computed group covariance matrices:\n");
    for (i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+lf ", pval_grp_cov[i * dim + j]);
        }

        printf("     ");

        for (j = 0; j < dim; j++)
        {
            printf("%+lf ", pval_grp_cov[i * dim + j + dim * dim]);
        }
        printf("\n");
    }

    printf("\nP-values of the computed group mean:\n");
    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", pval_grp_mean[i]);
    }

    printf("     ");

    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", pval_grp_mean[i + dim]);
    }

    /***** Printing P-values of the pooled covariance matrix *****/
    printf("\n\nP-values of the computed pooled covariance matrix:\n");

    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%9lf ", pval_pld_cov[i * dim + j]);
        }
        printf("\n");
    }

    printf("\nP-values of the computed pooled mean:\n");
    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", pval_pld_mean[i]);
    }

    /***** Printing summary of the test *****/
    if (errnums == 0) {
        printf("\n\nPooled and group covariance matrices");
        printf(" and mean estimates agree with theory\n");
    }
    else {
        printf("\n\nError: Pooled and group covariance matrices");
        printf(" and mean estimates disagree with theory\n");
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return 0;
}
