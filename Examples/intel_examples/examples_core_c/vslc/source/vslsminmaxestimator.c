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
!    Calculation of min/max estimates  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"

#define DIM     3      /* Task dimension */
#define N       1000   /* Number of observations */

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim;
    MKL_INT n;
    MKL_INT storage;
    float x[N][DIM];
    float min_est[DIM], max_est[DIM];
    float a = 0.0, sigma = 1.0;
    int i, j, errcode;
    int errnums = 0;

    /***** Initializing parameters for Summary Statistics task *****/
    dim     = DIM;
    n       = N;
    storage = VSL_SS_MATRIX_STORAGE_COLS;

    /***** Generate data set using VSL Gaussian RNG
    with mean a = 0 and stdev = 1 *****/
    errcode = sGenerateGaussianData( (float*)x, dim, n, a, sigma );
    CheckVslError(errcode);

    /***** Set initial values of the estimates *****/
    for(i = 0; i < dim; i++)
    {
        min_est[i] = max_est[i] = x[0][i];
    }

    /***** Create Summary Statistics task *****/
    errcode = vslsSSNewTask( &task, &dim, &n, &storage, (float*)x, 0, 0 );
    CheckVslError(errcode);

    /***** Edit task parameters for min and max computation *****/
    errcode = vslsSSEditTask( task, VSL_SS_ED_MIN, min_est );
    CheckVslError(errcode);

    errcode = vslsSSEditTask( task, VSL_SS_ED_MAX, max_est );
    CheckVslError(errcode);

    /***** Compute min and max estimates using FAST method *****/
    errcode = vslsSSCompute( task, VSL_SS_MIN|VSL_SS_MAX, VSL_SS_METHOD_FAST );
    CheckVslError(errcode);

    /***** Comparison of observations with min and max estimates *****/
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < n; j++)
        {
            if(x[j][i] < min_est[i]) errnums++;
            if(x[j][i] > max_est[i]) errnums++;
        }
    }

    /***** Printing results *****/
    printf("Task dimension : %d\n", (int)dim);
    printf("Number of observations : %d\n\n", (int)n);

    printf("               Min        Max\n");
    for(i = 0; i < dim; i++)
    {
        printf("Variable #%i:  %+lf  %+lf\n", i + 1,
            min_est[i], max_est[i]);
    }

    if (errnums == 0) {
        printf("\n\nAll observations are within ranges for all dimensions\n");
    }
    else {
        printf("\nError: There are %i observations beyond the ranges\n", errnums);
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return 0;
}
