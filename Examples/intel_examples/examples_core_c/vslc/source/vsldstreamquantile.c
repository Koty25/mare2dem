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
!    Calculation of quantiles, streaming data case  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "statchars.inc"

#define DIM     3      /* Task dimension */
#define N       1000   /* Number of observations */
#define M       9      /* Number of deciles */

#define EPS     0.01   /* Accuracy of quantile computation */

/* 1st component of the vector is processed */
MKL_INT indices[DIM] = { 1, 0, 0 };

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim;
    MKL_INT n;
    MKL_INT m;
    MKL_INT params_n;
    MKL_INT x_storage;
    MKL_INT o_storage;
    double x[DIM][N], xT[N][DIM];
    double o_stat[DIM][N];
    double o_quant[M], eq[DIM][M], sq[DIM][M];
    double min_abs_eq, min_abs_sq;
    double params;
    double a = 0.0, sigma = 1.0;
    int l1, l2;
    int i, j, k, errcode;
    int errnums = 0;

    /***** Initializing parameters for Summary Statistics task *****/
    dim           = DIM;
    n             = N;
    m             = M;
    x_storage     = VSL_SS_MATRIX_STORAGE_ROWS;
    o_storage     = VSL_SS_MATRIX_STORAGE_ROWS;
    params        = EPS;
    params_n      = VSL_SS_SQUANTS_ZW_PARAMS_N;

    /***** Generate transposed data set using VSL Gaussian RNG
    with mean a = 0 and stdev = 1 *****/
    errcode = dGenerateGaussianData( (double*)xT, dim, n, a, sigma );
    CheckVslError(errcode);

    for( j = 0; j < dim; j++ )
    {
        for( i = 0; i < n; i++ )
        {
            x[j][i] = xT[i][j];
        }
    }

    /***** Compute of deciles *****/
    for ( i = 0; i < m; i++ ) o_quant[i] = (double)(i + 1) / (double)(m + 1);

    /***** Create Summary Statistics task *****/
    errcode = vsldSSNewTask( &task, &dim, &n, &x_storage, (double*)x, 0,
                             indices );
    CheckVslError(errcode);

    /***** Edit task parameters for streaming quantile computation *****/
    errcode = vsldSSEditStreamQuantiles( task, &m, o_quant, (double*)sq,
                                         &params_n, &params );
    CheckVslError(errcode);

    /***** Compute streaming quantiles with accuracy EPS *****/
    errcode = vsldSSCompute( task, VSL_SS_STREAM_QUANTS,
                             VSL_SS_METHOD_SQUANTS_ZW );
    CheckVslError(errcode);

    /***** Edit task parameters for standard deciles and order statistics
           computation *****/
    errcode = vsldSSEditQuantiles( task, &m, o_quant, (double*)eq,
                                   (double*)o_stat, &o_storage );
    CheckVslError(errcode);

    /***** Compute exact quantiles and order statistica using FAST method *****/
    errcode = vsldSSCompute( task, VSL_SS_QUANTS|VSL_SS_ORDER_STATS,
                             VSL_SS_METHOD_FAST );
    CheckVslError(errcode);

    /***** Check the correctness of computed streaming quantiles *****/
    for( i = 0; i < dim; i++ )
    {
        if ( indices[i] )
        {
            for ( j = 0; j < m; j++ )
            {
                min_abs_eq = ABS( o_stat[i][0] - eq[i][j] );
                min_abs_sq = ABS( o_stat[i][0] - sq[i][j] );
                l1 = 0;
                l2 = 0;
                for ( k = 1; k < n; k++ )
                {
                    if ( ABS( o_stat[i][k] - eq[i][j] ) < min_abs_eq )
                    {
                        min_abs_eq = ABS( o_stat[i][k] - eq[i][j] );
                        l1 = k;
                    }
                    if ( ABS( o_stat[i][k] - sq[i][j] ) < min_abs_sq )
                    {
                        min_abs_sq = ABS( o_stat[i][k] - sq[i][j] );
                        l2 = k;
                    }
                }

                if ( ABS( l1 - l2 ) >= params * n ) errnums++;
            }
        }
    }

    /***** Printing results *****/
    printf("Task dimension         : %d\n", (int)dim);
    printf("Number of observations : %d\n\n", (int)n);

    /***** Printing part of the initial matrix of observations *****/
    printf("\n1st 4 and last 4 observations in source matrix\n");
    for ( i = 0; i < dim; i++ )
    {
        for ( j = 0; j < 4; j++ )
        {
            printf("%+.3f ", x[i][j] );
        }
        printf("     ...      ");
        for ( j = n - 5; j < n; j++ )
        {
            printf("%+.3f ", x[i][j] );
        }
        printf("\n");
    }

    /***** Printing computed exact quantiles *****/
    printf("\nStandard deciles of the observations:\n");
    printf("                 ");
    for ( j = 0; j < m; j++ ) printf("D%i     ", j + 1 );
    printf("\n");
    k = 0;
    for ( i = 0; i < dim; i++ )
    {
        if (indices[i])
        {
            printf("Variable #%d   ", i + 1);
            for ( j = 0; j < m; j++ ) printf("%+.3f ", eq[k][j] );
            printf("\n");
            k++;
        }
    }
    printf("\n");

    /***** Printing computed streaming quantiles *****/
    printf("\nStreaming deciles of the observations:\n");
    printf("                 ");
    for ( j = 0; j < m; j++ ) printf("D%i     ", j + 1 );
    printf("\n");
    k = 0;
    for ( i = 0; i < dim; i++ )
    {
        if (indices[i])
        {
            printf("Variable #%d   ", i + 1);
            for ( j = 0; j < m; j++ ) printf("%+.3f ", sq[k][j] );
            printf("\n");
            k++;
        }
    }

    /***** Printing summary of the test *****/
    if (errnums == 0) {
        printf("\n\nComputed streaming quantiles agree with theory.\n");
    }
    else {
        printf("\n\nError: Computed streaming quantiles disagree with");
        printf("  theory.\n");
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return 0;
}
