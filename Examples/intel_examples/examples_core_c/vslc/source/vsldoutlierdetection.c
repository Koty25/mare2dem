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
!    Outliers detection  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "statchars.inc"

#define DIM         10             /* Task dimension */
#define N           10000          /* Number of observations */
#define EPSILON     2              /* Percent of outliers in the observation */
#define OUTL_MEAN   30.0

/* BACON algorithm parameters */
#define NPARAMS     3
#define ALPHA       0.001 / N
#define BETA        0.005
#define INIT_METHOD VSL_SS_METHOD_BACON_MAHALANOBIS_INIT

int main()
{
    VSLSSTaskPtr task;
    MKL_INT n;
    MKL_INT dim;
    MKL_INT x_storage;
    double x[N*DIM];
    double weights[N];
    double C[DIM*DIM];
    double a[DIM];
    double a_outl[DIM];
    double ro, eps;
    MKL_INT cntoutl;
    int ndetoutl, ncorrdetoutl, nwrongdetoutl, nlostoutl;
    int i, j, errcode;
    int errnums = 0;

    MKL_INT nparams;
    double BACONparams[NPARAMS];
    double BACONweights[N];

    /***** Initializing parameters for Summary Statistics task *****/
    n           = N;
    dim         = DIM;
    nparams     = NPARAMS;
    x_storage   = VSL_SS_MATRIX_STORAGE_COLS;
    ro          = 1.0 / ( (double)dim * 30.0 );
    eps         = (double)EPSILON / 100.0;

    /***** Generate data set *****/
    /* Definition of covariance matrix and means for input data */
    for (i = 0; i < dim; i++)
    {
        a[i]  = 0.0;
        a_outl[i] = OUTL_MEAN;

        for(j = 0; j < i; j++)
        {
            C[i * dim + j] = ro;
            C[j * dim + i] = ro;
        }

        C[i * dim + i] = 1.0;
    }

    /* Generate matrix of observations containing outliers */
    errcode = dGenerateOutliersInput( x, dim, n, eps, a, a_outl,
                                      C, weights, &cntoutl );
    CheckVslError( errcode );

    /***** Create Summary Statistics task *****/
    errcode = vsldSSNewTask( &task, &dim, &n,  &x_storage, x, 0, 0 );
    CheckVslError( errcode );

    /***** Initializing Bacon algorithm parameters *****/
    BACONparams[0] = INIT_METHOD;
    BACONparams[1] = ALPHA;
    BACONparams[2] = BETA;

    for(i = 0; i < n; i++)
        BACONweights[i] = 0.0;

    /***** Initialization of the task parameters for Bacon algorithm *****/
    errcode = vsldSSEditOutliersDetection( task, &nparams, BACONparams,
                                           BACONweights );
    CheckVslError( errcode );

    /***** Detect outliers using BACON method *****/
    errcode = vsldSSCompute( task, VSL_SS_OUTLIERS,
                             VSL_SS_METHOD_BACON );
    CheckVslError( errcode );

    /***** Compute stat characteristics of the test *****/
    ndetoutl      = 0;
    ncorrdetoutl  = 0;
    nwrongdetoutl = 0;
    nlostoutl     = 0;

    for(i = 0; i < n; i++)
    {
        if (BACONweights[i] == 0.0)
        {
            ndetoutl++;
            if (weights[i] == 0.0)
                ncorrdetoutl++;
            else
                nwrongdetoutl++;
        }
        else // BACONweights[i] != 0.0
        {
            if (weights[i] == 0.0)
                nlostoutl++;
        }
    }
    errnums = nwrongdetoutl + nlostoutl;

    /***** Printing results *****/
    printf("Task dimension         : %d\n", (int)dim);
    printf("Number of observations : %d\n\n", (int)n);
    printf("Tail of Chi2 distribution alpha = %.8f\n", ALPHA);
    printf("Stopping criteria beta          = %.8f\n", BETA);

    /***** Printing summary results of BACON algorithm *****/
    printf("\n\nResults of BACON algorithm:\n");
    printf("\tnumber of outliers in dataset : %d\n", (int)cntoutl);
    printf("\tnumber of detected outliers   : %d\n",ndetoutl);

    printf("\nCheck of the output of BACON algorithm:\n");
    printf("\tRatio of correctly detected  outliers  = %3.2f%%\n",
          ( cntoutl > 0 ) ? 100.0 * (double)ncorrdetoutl / cntoutl : 0 );
    printf("\tRatio of incorrectly detected outliers = %3.2f%%\n",
          ( ndetoutl > 0 ) ? 100.0 * (double)nwrongdetoutl / ndetoutl : 0 );
    printf("\tRatio of missed outliers               = %3.2f%%\n",
          ( cntoutl > 0 ) ? 100.0 * (double)nlostoutl / cntoutl : 0 );

    /***** Printing summary of the test *****/
    if (errnums == 0)
    {
        printf("\nAll outliers are correctly detected.\n");
    }
    else
    {
        printf("\nError : Some of outliers are incorrectly detected");
        printf(" or missed.\n");
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError( errcode );

    MKL_Free_Buffers();

    return 0;
}
