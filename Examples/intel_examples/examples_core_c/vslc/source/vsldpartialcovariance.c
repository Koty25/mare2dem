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
!    Calculation of partial covariance matrix  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "statchars.inc"

#define DIM      4            /* Task dimension */
#define PART_DIM (DIM / 2)    /* Partial covariance dimension */

#define EPSILON  1e-8

double cov[DIM][DIM] = {
    {  1.0,  0.1, 0.15,  0.1 },
    {  0.1,  2.0,  0.1,  0.1 },
    { 0.15,  0.1,  1.0,  0.1 },
    {  0.1,  0.1,  0.1,  1.0 }
};

MKL_INT pcov_index[DIM] = { 1, 1, -1, -1 };

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim;
    MKL_INT pdim;
    MKL_INT cov_storage;
    MKL_INT pcov_storage;
    double cp_cov[DIM][DIM];
    double pcov[PART_DIM][PART_DIM];
    double th_pcov[PART_DIM][PART_DIM];
    int i, i1, j, j1, errcode;
    int errnums = 0;

    /***** Initializing parameters for Summary Statistics task *****/
    dim          = DIM;
    pdim         = PART_DIM;
    cov_storage  = VSL_SS_MATRIX_STORAGE_FULL;
    pcov_storage = VSL_SS_MATRIX_STORAGE_FULL;

    for(i = 0; i < PART_DIM; i++)
    {
        for(j = 0; j < PART_DIM; j++)
        {
            pcov[i][j] = 0.0;
        }
    }

    /***** Create Summary Statistics task *****/
    errcode = vsldSSNewTask( &task, &dim, 0, 0, 0, 0, 0 );
    CheckVslError(errcode);

    /***** Edit task parameters for partial covariance matrix computation *****/
    errcode = vsldSSEditPartialCovCor( task, (MKL_INT*)pcov_index,
                                       (double*)cov, &cov_storage, 0, 0,
                                       (double*)pcov, &pcov_storage, 0,
                                       &pcov_storage );
    CheckVslError(errcode);

    /***** Compute partial covariance matrix using FAST method *****/
    errcode = vsldSSCompute( task, VSL_SS_PARTIAL_COV,
                             VSL_SS_METHOD_FAST );
    CheckVslError(errcode);

    /***** Printing results *****/
    printf("Task dimension : %d\n\n", (int)dim);

    /* Print input covariance matrix */
    printf(" Covariance matrix\n");
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+lf ", cov[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    /* Print computed partial covariance matrix estimate */
    printf(" Computed partial covariance matrix\n");
    for(i = 0; i < pdim; i++)
    {
        for(j = 0; j < pdim; j++)
        {
            printf("%+lf ", pcov[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    /***** Testing stat characteristics of partial covariance matrix *****/
    /* Compute theoretical partial covariance estimate using sweep operator */
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            cp_cov[i][j] = cov[i][j];
        }
    }

    for(i = 0; i < dim; i++)
    {
        if (pcov_index[i] == -1)
            dSweep( i, dim, (double*)cp_cov );
    }

    i1 = 0;
    j1 = 0;
    for(i = 0; i < dim; i++)
    {
        if (pcov_index[i] == 1)
        {
            j1 = 0;
            for(j = 0; j < dim; j++)
            {
                if (pcov_index[j] == 1)
                {
                    th_pcov[i1][j1] = cp_cov[i][j];
                    j1++;
                }
            }
            i1++;
        }
    }

    /* Print theoretical partial covariance estimate */
    printf(" Theoretical partial covariance matrix\n");
    for(i = 0; i < pdim; i++)
    {
        for(j = 0; j < pdim; j++)
        {
            printf("%+lf ", th_pcov[i][j]);
        }
        printf("\n");
    }

    /* Check the correctness of computed partial covariance matrix */
    for(i = 0; i < pdim; i++)
    {
        for(j = 0; j < pdim; j++)
        {
            if(ABS(pcov[i][j] - th_pcov[i][j]) > EPSILON) errnums++;
        }
    }

    /***** Printing summary of the test *****/
    if (errnums == 0)
    {
        printf("\n\nComputed partial covariance matrix estimate");
        printf(" agrees with theory\n");
    }
    else
    {
        printf("\n\nError: Computed partial covariance matrix estimate");
        printf(" disagrees with theory\n");
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return 0;
}
