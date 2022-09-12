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
!    Calculation of group/pooled means Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "statchars.inc"

#define DIM     3      /* Task dimension */
#define N       10000  /* Number of observations */
#define G       2      /* Number of groups */
#define GN      2      /* Number of group means */

#define P_THRESHOLD    0.005

float C[DIM][DIM] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
};

float m[DIM] = { 0.0, 0.0, 0.0 };

/***** 1st and 2nd group means to be returned *****/
MKL_INT group_mean_indices[G] = { 1, 1 };

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim=DIM, n=N, x_storage;
    float x[DIM*N], pld_mean[DIM], grp_mean[DIM*GN];
    MKL_INT group_indices[N];

    float a = 0.0, sigma = 1.0;

    int i, errcode, ret_value;
    int errnums = 0;

    float pval_pld_mean[DIM], pval_grp_mean[DIM*GN];

    /***** Generate data set using VSL Gaussian RNG with mean a = 0 and
           stdev = 1 *****/
    errcode = sGenerateGaussianData( x, dim, n, a, sigma );
    CheckVslError(errcode);

    /***** Initializing parameters for Summary Statistics task *****/
    dim              = DIM;
    n                = N;
    x_storage        = VSL_SS_MATRIX_STORAGE_ROWS;

    /***** Create Summary Statistics task *****/
    errcode = vslsSSNewTask( &task, &dim, &n, &x_storage, x, 0, 0 );
    CheckVslError(errcode);

    /***** Dividing elements into odd and even *****/
    for(i = 0; i < n; i++)
    {
        group_indices[i] = i % 2;
    }

    /***** Initialization of the task parameters for pooled and
           group mean estimators *****/
    errcode = vslsSSEditPooledCovariance( task, group_indices,
        pld_mean, 0, group_mean_indices, grp_mean, 0 );
    CheckVslError(errcode);

    /***** Compute group and pooled mean using 1PASS method  *****/
    errcode = vslsSSCompute( task,VSL_SS_POOLED_MEAN | VSL_SS_GROUP_MEAN,
                                  VSL_SS_METHOD_1PASS );
    CheckVslError(errcode);

    /***** Testing stat characteristics of means *****/
    /* Compute p-values for group mean estimates */
    sComputePvalsMean( dim, n, grp_mean, m, (float*)C, pval_grp_mean );
    sComputePvalsMean( dim, n, &grp_mean[dim], m, (float*)C,&pval_grp_mean[dim] );
    /* Compute p-values for pooled mean estimates */
    sComputePvalsMean( dim, n, pld_mean, m, (float*)C, pval_pld_mean );

    for(i = 0; i < dim; i++)
    {
        if (pval_grp_mean[i] < P_THRESHOLD) errnums++;
        if (pval_grp_mean[i + dim] < P_THRESHOLD) errnums++;
        if (pval_pld_mean[i] < P_THRESHOLD) errnums++;
    }

    /***** Printing results *****/
    printf("Task dimension : %d\n", (int)dim);
    printf("Number of observations : %d\n", (int)n);

    /***** Print exact mean *****/
    printf("\nExact means:\n");
    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", m[i]);
    }

    /***** Print group mean estimates *****/
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

    /***** Print pooled mean estimates *****/
    printf("\nPooled means:\n");
    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", pld_mean[i]);
    }

    /***** Print P-values of the group means *****/
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

    /***** Printing P-values of the pooled mean *****/
    printf("\nP-values of the computed pooled mean:\n");
    for(i = 0; i < dim; i++)
    {
        printf("%+lf ", pval_pld_mean[i]);
    }

    /***** Printing summary of the test *****/
    if (errnums == 0) {
        printf("\n\nPooled and group mean estimates ");
        printf("agree with theory\n");
        ret_value = 0;
    }
    else {
        printf("\n\nPooled and group mean estimates ");
        printf("disagree with theory\n");
        ret_value = 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return ret_value;
}
