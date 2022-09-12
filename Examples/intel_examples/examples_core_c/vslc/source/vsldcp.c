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
!    Calculation of cross-product matrix Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "statchars.inc"

#define DIM     4        /* Task dimension */
#define N       1000     /* Number of observations */

#define P_THRESHOLD     0.01

double C[DIM][DIM] = {
    { 1.0, 0.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0, 0.0 },
    { 0.0, 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 0.0, 1.0 }
};
double a[DIM] = { 5.0, 5.0, 5.0, 5.0 };

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim;
    MKL_INT n;
    MKL_INT x_storage;
    MKL_INT cov_storage, cp_storage;

    double x[N][DIM];
    double mean[DIM], cp[DIM][DIM], cov[DIM][DIM];
    int i, j, errcode;
    unsigned MKL_INT64 estimate = 0;
    int errnums = 0;

    double pval_cov[DIM][DIM];

    /***** Initializing parameters for Summary Statistics task *****/
    dim         = DIM;
    n           = N;
    x_storage   = VSL_SS_MATRIX_STORAGE_COLS;
    cov_storage = VSL_SS_MATRIX_STORAGE_FULL;
    cp_storage  = VSL_SS_MATRIX_STORAGE_FULL;

    /***** Generate data set using VSL GaussianMV RNG *****/
    errcode = dGenerateGaussianMVData( (double*)x, dim, n, a, (double*)C );
    CheckVslError(errcode);

    /***** Create Summary Statistics task *****/
    errcode = vsldSSNewTask( &task, &dim, &n, &x_storage, (double*)x, 0, 0 );
    CheckVslError(errcode);

    /***** Initialization of the task parameters related to
           cross-product matrix computation *****/
    vsldSSEditCP(task, mean, 0, (double*)cp, &cp_storage);
    CheckVslError(errcode);

    /***** Cross-product matrix is included in the list of estimates
           to compute *****/
    estimate = VSL_SS_CP;

    /***** Compute the estimates using FAST method *****/
    errcode = vsldSSCompute( task, estimate, VSL_SS_METHOD_FAST );
    CheckVslError(errcode);


    /***** Edit task parameters for computation of covariance *****/
    errcode = vsldSSEditTask( task, VSL_SS_ED_COV, (double*)cov );
    CheckVslError(errcode);
    errcode = vsliSSEditTask( task, VSL_SS_ED_COV_STORAGE, &cov_storage );
    CheckVslError(errcode);

    /***** Convert cross-product matrix into correlation matrix *****/
    errcode = vsldSSCompute( task, VSL_SS_COV, VSL_SS_METHOD_CP_TO_COVCOR );
    CheckVslError(errcode);

    /***** Testing stat characteristics of the computed estimates *****/
    /* Compute p-values for covariance computed from cross-product */
    dComputePvalsVariance( dim, n, (double*)cov, cov_storage, (double*)C,
                           (double*)pval_cov );
    dComputePvalsCovariance( dim, n, (double*)cov, cov_storage, (double*)C,
                             (double*)pval_cov );

    /***** Checking the validity of p-values for all estimates *****/
    for(i = 0; i < dim; i++)
    {
        if (pval_cov[i][i] < P_THRESHOLD) errnums++;
        for(j = 0; j < i; j++)
        {
            if (pval_cov[i][j] < P_THRESHOLD) errnums++;
        }
    }

    /***** Printing results *****/
    printf("Task dimension : %d\n", (int)dim);
    printf("Number of observations : %d\n\n", (int)n);

    /***** Printing computed cross-product matrix *****/
    printf("\n Computed cross-product matrix\n");
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+9lf ", cp[i][j]);
        }
        printf("\n");
    }

    /***** Printing p-values for covariance matrix estimate *****/
    printf("\n P-values for covariance obtained from cross-product matrix\n");
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
        printf("\n\nAll the computed estimates agree with theory\n");
    }
    else
    {
        printf("\n\nError: At least one of the computed estimates");
        printf(" disagrees with theory\n");
        return 1;
    }

    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return errcode;
}
