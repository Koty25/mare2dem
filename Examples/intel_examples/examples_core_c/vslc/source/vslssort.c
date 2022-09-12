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
!    Sorting data array  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "mkl_vsl.h"
#include "errcheck.inc"
#include "generatedata.inc"

#define P       3    /* Number of observations */
#define N       10   /* Task dimension */
#define SEED    1    /* Initial value for stream initialization */

int main()
{
    VSLSSTaskPtr task;
    MKL_INT n;
    MKL_INT p;
    MKL_INT x_storage;
    MKL_INT sorted_x_storage;
    float x[P][N];                    /* Matrix of observations */
    float y[P][N];                    /* Output  matrix for sorted data */
    float lBound = 0.0, rBound = 10;  /* Bounds for uniform generator */
    int errnums = 0;
    int i, j, errcode;
    VSLStreamStatePtr stream;

    /***** Initializing parameters for Summary Statistics task *****/
    n = N;
    p = P;
    x_storage        = VSL_SS_MATRIX_STORAGE_ROWS;
    sorted_x_storage = VSL_SS_MATRIX_STORAGE_ROWS;

    /***** Generate data set using VSL Uniform RNG *****/
    /***** Initialize *****/
    errcode = vslNewStream(&stream, VSL_BRNG_MCG31, SEED);
    CheckVslError(errcode);

    /***** Call RNG *****/
    errcode = vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, N * P, (float *) x, lBound, rBound);
    CheckVslError(errcode);

    /***** Create Summary Statistics task *****/
    errcode = vslsSSNewTask(&task, &p, &n, &x_storage, (float *) x, 0, 0);
    CheckVslError(errcode);

    /***** Edit task parameters for sorting *****/
    errcode = vslsSSEditTask(task, VSL_SS_ED_SORTED_OBSERV, (float *) y);
    CheckVslError(errcode);

    /***** Edit task parameters for sorting *****/
    errcode = vsliSSEditTask(task, VSL_SS_ED_SORTED_OBSERV_STORAGE, &sorted_x_storage);
    CheckVslError(errcode);

    /***** Sort data using radix method *****/
    errcode = vslsSSCompute(task, VSL_SS_SORTED_OBSERV, VSL_SS_METHOD_RADIX);
    CheckVslError(errcode);

    /***** Check the correctness of sorting *****/
    for( j = 0; j < p; j++ )
    {
        for ( i = 0; i < n - 1; i++ )
        {
            if ( y[j][i] > y[j][i + 1] )
            {
                errnums++;
            }
        }
    }

    /***** Printing results *****/
    printf("Task dimension : %d\n", (int)p);
    printf("Number of observations : %d\n\n", (int)n);

    /***** Printing of the initial matrix *****/
    printf("\n Initial matrix: \n");
    for ( j = 0; j < p; j++ )
    {
        for ( i = 0; i < n; i++ )
        {
            printf("%+.3f   ", x[j][i]);
        }
        printf("\n");
    }
    printf("\n");

    /***** Printing of the sorted matrix *****/
    printf("\n Sorted matrix: \n");
    for ( j = 0; j < p; j++ )
    {
        for ( i = 0; i < n; i++ )
        {
            printf("%+.3f   ", y[j][i]);
        }
        printf("\n");
    }

    /***** Printing summary of the test *****/
    if (errnums == 0)
    {
        printf("\n\n Sorting is correct \n\n");
    }
    else
    {
        printf("\n\n Error: Sorting is incorrect \n\n");
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask(&task);
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return 0;
}
