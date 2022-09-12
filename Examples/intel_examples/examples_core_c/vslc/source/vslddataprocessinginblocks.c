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
!    Data processing in chunks  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "statchars.inc"

#define DIM             5           /* Task dimension */
#define N               100000      /* Number of observations */
#define NBLOCKS         100         /* Number of data portions */
#define BLOCK_SIZE      (N / NBLOCKS)

#define P_THRESHOLD     0.05

/* Exact covariance matrix */
double C[DIM][DIM] = {
    { 1.0,  0.05, 0.05, 0.05, 0.05 },
    { 0.05,  1.0, 0.05, 0.05, 0.05 },
    { 0.05, 0.05,  1.0, 0.05, 0.05 },
    { 0.05, 0.05, 0.05,  1.0, 0.05 },
    { 0.05, 0.05, 0.05, 0.05,  1.0 }
};

/* Exact means vector */
double a[DIM] = { -1.0, 1.0, 2.0, -3.0, 4.0 };

/* Array of accumulated weights */
double W[2] = { 0.0, 0.0 };

MKL_INT indices[DIM] = { 0, 1, 1, 0, 1 };

int main(void)
{
    VSLSSTaskPtr task;
    VSLStreamStatePtr stream;
    MKL_INT dim;
    MKL_INT n;
    MKL_INT nblocks, block_size;
    MKL_INT x_storage;
    MKL_INT cov_storage;
    double x[BLOCK_SIZE*DIM];
    double cov[DIM][DIM];
    double weights[N];
    double mean[DIM];
    double r2[DIM], c2[DIM];
    double max_est[DIM], min_est[DIM];
    double T[DIM*DIM];
    int i, j, k, errcode;
    unsigned MKL_INT64 estimate = 0;
    int errnums = 0;

    double pval_mean[DIM];
    double pval_cov[DIM][DIM];
    double pval_r2[DIM], pval_c2[DIM];

    /***** Initializing parameters for Summary Statistics task *****/
    dim         = DIM;
    n           = N;
    nblocks     = NBLOCKS;
    block_size  = BLOCK_SIZE;
    x_storage   = VSL_SS_MATRIX_STORAGE_COLS;
    cov_storage = VSL_SS_MATRIX_STORAGE_FULL;

    for(i = 0; i < n; i++)
    {
        weights[i] = 1.0;
    }

    /***** Generate data set using VSL GaussianMV RNG *****/
    stream = dInitGaussianMVDataGenerator( dim, (double*)C, T );

    /***** Create Summary Statistics task *****/
    errcode = vsldSSNewTask( &task, &dim, &block_size,  &x_storage,
                             x, weights, indices );
    CheckVslError(errcode);

    /***** Register array of weights in the task *****/
    errcode = vsldSSEditTask( task, VSL_SS_ED_ACCUM_WEIGHT, W );
    CheckVslError(errcode);

    /***** Edit task parameters for computing of mean estimate and
           2nd raw and central moments estimates *****/
    errcode = vsldSSEditMoments( task, mean, r2, 0, 0, c2, 0, 0 );
    CheckVslError(errcode);

    /***** Initialization of the task parameters using FULL_STORAGE
           for covariance matrix computation *****/
    errcode = vsldSSEditCovCor( task, mean, (double*)cov, &cov_storage, 0, 0 );
    CheckVslError(errcode);

    /***** Edit task parameters for min and max computation *****/
    errcode = vsldSSEditTask( task, VSL_SS_ED_MAX, max_est);
    CheckVslError(errcode);

    errcode = vsldSSEditTask( task, VSL_SS_ED_MIN, min_est );
    CheckVslError(errcode);

    /***** Minimum and maximum are included in the list of estimates
           to compute *****/
    estimate |= VSL_SS_MIN | VSL_SS_MAX;

    /***** Mean and 2nd raw and central moments are included in the list
           of estimates to compute *****/
    estimate |= VSL_SS_MEAN | VSL_SS_2R_MOM | VSL_SS_2C_MOM;

    /***** Covariance matrix is included in the list of estimates
           to compute *****/
    estimate |= VSL_SS_COV;

    for(i = 0; i < nblocks; i++)
    {
        /***** Generate new portion of data using VSL GaussianMV RNG *****/
        errcode = dGenerateGaussianMVDataBlock( stream, x, dim, block_size,
                                                a, T );
        CheckVslError(errcode);

        if (i == 0)
        {
            for(j = 0; j < dim; j++)
            {
                min_est[j] = x[j];
                max_est[j] = x[j];
            }
        }

        /***** Compute the estimates using FAST method *****/
        errcode = vsldSSCompute( task, estimate, VSL_SS_METHOD_FAST );
        CheckVslError(errcode);

        /* Comparison of observations with min and max estimates */
        for(k = 0; k < dim; k++)
        {
            if (indices[k] != 0)
            {
                for(j = 0; j < block_size; j++)
                {
                    if(x[j * dim + k] < min_est[k]) errnums++;
                    if(x[j * dim + k] > max_est[k]) errnums++;
                }
            }
        }
    }

    /***** Testing stat characteristics of the computed estimates *****/
    /* Compute p-values for mean estimates */
    dComputePvalsMean( dim, n, mean, a, (double*)C, pval_mean );
    /* Compute p-values for variance estimates */
    dComputePvalsVariance( dim, n, (double*)cov, cov_storage,
                           (double*)C, (double*)pval_cov );
    /* Compute p-values for covariance estimates */
    dComputePvalsCovariance( dim, n, (double*)cov, cov_storage,
                             (double*)C, (double*)pval_cov );
    /* Compute p-values for raw moments estimates */
    dComputePvalsRawMoments( dim, n, r2, 2, a, (double*)C,
                             pval_r2 );
    /* Compute p-values for central moments estimates */
    dComputePvalsCentralMoments( dim, n, c2, 2, a, (double*)C,
                                 pval_c2 );

    /***** Checking the validity of p-values for all estimates *****/
    for(i = 0; i < dim; i++)
    {
        if (indices[i] != 0)
        {
            if (pval_mean[i] < P_THRESHOLD) errnums++;
            if (pval_r2[i]   < P_THRESHOLD) errnums++;
            if (pval_c2[i]   < P_THRESHOLD) errnums++;

            for(j = 0; j <= i; j++)
            {
                if (indices[j] != 0)
                {
                    if (pval_cov[i][j] < P_THRESHOLD) errnums++;
                }
            }
        }
    }

    /***** Printing results *****/
    printf("Task dimension : %d\n", (int)dim);
    printf("Number of observations : %d\n", (int)n);
    printf("Number of blocks : %d\n\n", (int)nblocks);

    /***** Printing computed minimum, maximum, mean and moments estimates *****/
    printf("               Min        Max        Mean       2nd_raw  ");
    printf("  2nd_cen\n");
    for(i = 0; i < dim; i++)
    {
        if (indices[i] != 0)
        {
            printf("Variable #%i:  %+lf  %+lf  %+lf  %+lf  %+lf\n", i + 1,
                min_est[i], max_est[i], mean[i],
                r2[i], c2[i]);
        }
    }

    /***** Printing computed covariance matrix *****/
    printf("\n Computed covariance matrix\n");
    for(i = 0; i < dim; i++)
    {
        if (indices[i])
        {
            printf("Variable #%i:  ", i + 1);
            for(j = 0; j < dim; j++)
            {
                if(indices[j])
                {
                    printf("%+9lf ", cov[i][j]);
                }
            }
            printf("\n");
        }
    }

    /***** Printing p-values for mean and moments estimates *****/
    printf("\n\nP-values of the computed estimates\n\n");
    printf("               Mean       2nd_raw    2nd_cen\n");
    for(i = 0; i < dim; i++)
    {
        if (indices[i] != 0)
        {
            printf("Variable #%i:  %+lf  %+lf  %+lf\n",
                i + 1, pval_mean[i], pval_r2[i], pval_c2[i]);
        }
    }

    /***** Printing p-values for covariance matrix estimate *****/
    printf("\n Covariance matrix\n");
    for(i = 0; i < dim; i++)
    {
        if (indices[i])
        {
            printf("Variable #%i:  ", i + 1);
            for(j = 0; j < dim; j++)
            {
                if(indices[j])
                {
                    printf("%+9lf ", pval_cov[i][j]);
                }
            }
            printf("\n");
        }
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

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError(errcode);

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return 0;
}
