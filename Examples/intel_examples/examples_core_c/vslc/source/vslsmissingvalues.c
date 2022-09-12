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
!  Support of missing values, Multiple Imputation algorithm Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"

#define DIM         4      /* Task dimension                   */
#define N           1000   /* Number of observations           */
#define EPSILON     20     /* Missing values percentage        */
#define NPATT       9      /* Number of different patterns     */
#define M           5      /* Number of sets of imputed values */

#define EM_ITER_NUM 100.0  /* Number of iterations of EM algorithm */
#define DA_ITER_NUM 10     /* Number of iterations of DA algorithm */
#define EM_ACCURACY 0.001  /* Accuracy of EM algorithm             */

/* Number of simulated missing values */
#define MI_SIMVALS_N    20000

/* Initial number of est */
#define MI_INIT_EST_N   (DIM + DIM * (DIM + 1) / 2)

/* Number of est */
#define MI_EST_N        M * DA_ITER_NUM * DIM * (DIM + 3) / 2

/* Exact covariance matrix */
float C[DIM][DIM] = {
    {   1.0, 0.125, 0.125, 0.125 },
    { 0.125,   1.0, 0.125, 0.125 },
    { 0.125, 0.125,   1.0, 0.125 },
    { 0.125, 0.125,  0.125,  1.0 }
};

/* Exact mean */
float a[DIM] = { 0.0, 0.0, 0.0, 0.0 };

/* Matrix of patterns */
int patt[NPATT][DIM] = {
    { 1, 0, 0, 0 },
    { 0, 1, 0, 0 },
    { 0, 0, 1, 0 },
    { 0, 0, 0, 1 },
    { 1, 0, 1, 0 },
    { 0, 0, 1, 1 },
    { 1, 1, 0, 1 },
    { 1, 0, 1, 1 },
    { 1, 1, 1, 1 }
};

float W[2] = {0.0, 0.0};

int main()
{
    VSLSSTaskPtr task;
    MKL_INT n;
    MKL_INT dim;
    MKL_INT m;
    MKL_INT n_miss_vals;
    MKL_INT index;
    MKL_INT x_storage;
    MKL_INT q_storage;
    int npatt;
    int eps;
    float x[DIM][N];
    float mean[DIM*M], raw2[DIM*M], variance[DIM*M];
    float qmean[DIM], qvariance[DIM];
    float B[DIM];               /* "between-imputation" variance vector */
    float U[DIM];               /* "within-imputation" variance vector */
    float T[DIM];               /* "total" variance vector */
    float sqrtT[DIM];
    float cint_left[DIM], cint_right[DIM];
    int miss_vals[N];
    int ipatt;
    int i, j, k, errcode;
    int errnums = 0;

    MKL_INT nparams;
    MKL_INT n_simvals, n_est, n_init_est;
    float params[VSL_SS_MI_PARAMS_SIZE];
    float simvals[MI_SIMVALS_N], est[MI_EST_N], init_est[MI_INIT_EST_N];

    /***** Initializing parameters for Summary Statistics task *****/
    n          = N;
    dim        = DIM;
    m          = M;
    x_storage  = VSL_SS_MATRIX_STORAGE_ROWS;
    q_storage  = VSL_SS_MATRIX_STORAGE_COLS;
    npatt      = NPATT;
    eps        = EPSILON;
    nparams    = VSL_SS_MI_PARAMS_SIZE;
    n_simvals  = MI_SIMVALS_N;
    n_est      = MI_EST_N;
    n_init_est = MI_INIT_EST_N;

    for(i = 0; i < dim; i++)
    {
        qmean[i] = 0.0;
    }

    for(i = 0; i < dim*m; i++)
    {
        mean[i]     = 0.0;
        raw2[i]     = 0.0;
        variance[i] = 0.0;
    }

    /***** Define covariance matrix and mean of input data *****/
    k = 0;
    for(i = 0; i < dim; i++)
    {
        init_est[i] = a[i];
        for(j = i; j < dim; j++)
        {
            init_est[dim + k] = C[i][j];
            k++;
        }
    }

    /***** Generate dataset with missing values *****/
    for(i = 0; i < n; i++)
    {
        miss_vals[i] = -1;
    }

    errcode = sGenerateMissingValuesInput( (float*)x, dim, n, eps,
                                           (int*)patt, npatt, miss_vals,
                                           a, (float*)C, &n_miss_vals );
    CheckVslError( errcode );

    /***** Create Summary Statistics task *****/
    errcode = vslsSSNewTask( &task, &dim, &n, &x_storage, (float*)x, 0, 0 );
    CheckVslError( errcode );

    /***** Initialization of the task parameters for missing values
           estimator *****/
    n_simvals = n_miss_vals * m;

    params[0] = EM_ITER_NUM;
    params[1] = DA_ITER_NUM;
    params[2] = EM_ACCURACY;
    params[3] = (float)m;
    params[4] = (float)n_miss_vals;

    errcode = vslsSSEditMissingValues( task, &nparams, params,
                                       &n_init_est, init_est, 0, 0,
                                       &n_simvals, simvals, &n_est, est );
    CheckVslError( errcode );

    /***** Compute missing values est using MULTIPLE IMPUTATION
           method *****/
    errcode = vslsSSCompute( task, VSL_SS_MISSING_VALS,
                             VSL_SS_METHOD_MI );
    CheckVslError( errcode );

    /***** Compute mean and variance for all imputations *****/
    index = 0;
    for(k = 0; k < m; k++)
    {
        /* Recover of the input matrix */
        for(j = 0; j < n; j++)
        {
            ipatt = miss_vals[j];
            if (ipatt >= 0)
            {
                /* Here if there are missing values in j-th observation */
                for(i = 0; i < dim; i++)
                {
                    if (patt[ipatt][i])
                    {
                        x[i][j] = simvals[index];
                        index++;
                    }
                }
            }
        }

        W[0] = 0.0;
        W[1] = 0.0;

        /* Initializing of the task parameters for mean, 2nd raw moment and
           variance estimators computation */
        errcode = vslsSSEditMoments( task, mean + k * dim, raw2 + k * dim,
                                     0, 0, variance + k * dim, 0, 0 );
        CheckVslError( errcode );

        /* Compute mean and variance using fast method */
        errcode = vslsSSCompute( task, VSL_SS_MEAN | VSL_SS_2C_MOM,
                                 VSL_SS_METHOD_FAST );
        CheckVslError( errcode );
    } //for(k = 0; k < m; k++)

    /***** Testing stat characteristics of computed estimates *****/

    /* Register matrix of means estimates as observations matrix in
       Summary Statistics task */
    errcode = vslsSSEditTask( task, VSL_SS_ED_OBSERV, mean );
    CheckVslError( errcode );

    /* Set proper number of observations */
    errcode = vsliSSEditTask( task, VSL_SS_ED_OBSERV_N, &m );
    CheckVslError( errcode );

    /* Register proper observation storage format in the task */
    errcode = vsliSSEditTask( task, VSL_SS_ED_OBSERV_STORAGE,
                              &q_storage );
    CheckVslError( errcode );

    /* Initializing of the task parameters for mean, 2nd raw moment and
       variance estimators computation */
    errcode = vslsSSEditMoments( task, qmean, raw2, 0, 0, B, 0, 0);
    CheckVslError( errcode );

    W[0] = 0.0;
    W[1] = 0.0;

    /***** Compute mean, 2nd raw moments and variance est for
           matrix of means est for all datasets *****/
    errcode = vslsSSCompute( task, VSL_SS_MEAN | VSL_SS_2C_MOM,
                             VSL_SS_METHOD_FAST );
    CheckVslError( errcode );

    /* Calculate average of variance esimates */
    for(j = 0; j < dim; j++)
    {
        qvariance[j] = 0;
        for(i = 0; i < m; i++)
        {
            qvariance[j] += variance[i*dim + j];
        }
        qvariance[j] /= (float)m;
    }

    /* Calculate "within imputation" variance */
    for(j = 0; j < dim; j++)
    {
        U[j] = qvariance[j] / n;
    }

    /* Calculate "total" variance */
    for(j = 0; j < dim; j++)
    {
        T[j] = U[j] + (1.0 + 1.0 / (float)m) * B[j];
    }

    /* Compute borders of 95% confidence interval for mean estimates */
    vsSqrt( dim, T, sqrtT );
    for(j = 0; j < dim; j++)
    {
        cint_right[j] = qmean[j] + 2*sqrtT[j];
        cint_left[j]  = qmean[j] - 2*sqrtT[j];
    }

    for(j = 0; j < dim; j++)
    {
        if( (cint_left[j] > a[j]) || (cint_right[j] < a[j]) )
        {
            errnums++;
        }
    }

    /***** Printing results *****/
    printf("Task dimension :         %d\n", (int)dim);
    printf("Number of observations : %d\n\n", (int)n);
    printf("Number of initial estimates :        %d\n", (int)n_init_est);
    printf("Number of simulated missing values : %d\n", (int)n_simvals);
    printf("Number of estimates :                %d\n", (int)n_est);
    printf("Number of missing values :           %d\n", (int)n_miss_vals);

    /***** Print exact means and variances *****/
    printf("\n Exact mean:\n        ");
    for(i = 0; i < dim; i++)
    {
        printf(" %lf ", a[i]);
    }

    printf("\n Exact variance:\n        ");
    for(i = 0; i < dim; i++)
    {
        printf(" %lf ", C[i][i]);
    }

    /***** Print computes means and variances for all sets
           of imputed values *****/
    printf("\n\nSet:     Mean:\n");
    for(i = 0; i < m; i++)
    {
        printf("   %d    ", i + 1);
        for(j = 0; j < dim; j++)
        {
            printf("%+6f ", mean[i*dim + j]);
        }
        printf("\n");
    }

    printf("\nAverage ");
    for(j = 0; j < dim; j++)
    {
        printf("%+6f ", qmean[j]);
    }

    printf("\n\n\nSet:     Variance:\n");
    for(i = 0; i < m; i++)
    {
        printf("   %d    ", i + 1);
        for(j = 0; j < dim; j++)
        {
            printf(" %lf ", variance[i*dim + j]);
        }
        printf("\n");
    }

    printf("\nAverage ");
    for(j = 0; j < dim; j++)
    {
        printf(" %lf ", qvariance[j]);
    }
    printf("\n");

    /***** Print between-imputation, within-imputation and total
           variances *****/
    printf("\nBetween-imputation variance:\n");
    printf("\n        ");
    for(j = 0; j < dim; j++)
    {
        printf("%-6f ", B[j]);
    }
    printf("\n");

    printf("\nWithin-imputation variance:\n");
    printf("\n        ");
    for(j = 0; j < dim; j++)
    {
        printf("%-6f ", U[j]);
    }
    printf("\n");

    printf("\nTotal variance:\n");
    printf("\n        ");
    for(j = 0; j < dim; j++)
    {
        printf("%#6f ", T[j]);
    }
    printf("\n");

    /***** Print borders of 95% confidence intervals for mean estimates *****/
    printf("\n95%% confidence interval:\n");
    printf("\n right  ");
    for(j = 0; j < dim; j++)
    {
        printf("%+6f ", cint_right[j]);
    }
    printf("\n left   ");
    for(j = 0; j < dim; j++)
    {
        printf("%+6f ", cint_left[j]);
    }
    printf("\n");

    /***** Printing summary of the test *****/
    if (errnums == 0)
    {
        printf("\n\nComputed missing values estimates");
        printf(" agree with theory\n");
    }
    else
    {
        printf("\n\nError: Computed missing values estimates");
        printf(" disagree with theory\n");
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return 0;
}
