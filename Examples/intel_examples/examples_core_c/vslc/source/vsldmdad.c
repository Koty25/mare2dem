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
!    Calculation of median absolute deviation  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include <math.h>
#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"

#define DIM     3       /* Task dimension */
#define N       1000    /* Number of observations */

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim;
    MKL_INT n;
    MKL_INT x_storage;

    double x[DIM][N];       /* matrix of observations */
    double mdad[DIM];
    double a = 0.0, sigma = 1.0;
    int i, errcode;
    int errnums = 0;
    double dn=(double)N;

    double tD, tQ, tD2, s, sig, sD, deltaD;

    /***** Initializing parameters for Summary Statistics task *****/
    dim              = DIM;
    n                = N;
    x_storage        = VSL_SS_MATRIX_STORAGE_ROWS;

    /***** Generate transposed data set using VSL Gaussian RNG
    with mean a = 0 and stdev = 1 *****/
    errcode = dGenerateGaussianData( (double*)x, dim, n, a, sigma );
    CheckVslError(errcode);

    /***** Create Summary Statistics task *****/
    errcode = vsldSSNewTask( &task, &dim, &n, &x_storage, (double*)x, 0, 0 );
    CheckVslError(errcode);

    /***** Provide array of median absolute deviation *****/
    errcode = vsldSSEditTask( task, VSL_SS_ED_MDAD, mdad );
    CheckVslError(errcode);

    /***** Compute median absolute deviation using FAST method *****/
    errcode = vsldSSCompute( task, VSL_SS_MDAD, VSL_SS_METHOD_FAST );
    CheckVslError(errcode);

    /***** Check the correctness of computed median absolute deviations *****/
    /***** Testing relies on property that for Gaussian distribution
           standard deviation estimate ~= 1.4826 * mdad ******/
    tD=sigma*sigma;
    tQ=720.0*sigma*sigma*sigma*sigma;
    tD2=tD*tD;
    s=((tQ-tD2)/dn)-(2.0*(tQ-2*tD2)/(dn*dn))+((tQ-3*tD2)/(dn*dn*dn));

    for ( i = 0; i < dim; i++ )
    {
        sig = 1.4826 * mdad[i];
        sD = sig * sig;
        deltaD = fabs((tD-sD) / sqrt(s));
        if  ( deltaD > 3.0 ) errnums++;
    }

    /***** Printing results *****/
    printf("Task dimension : %d\n", (int)dim);
    printf("Number of observations : %d\n\n", (int)n);

    /***** Printing computed median absolute deviations *****/
    printf("\nMedian absolute deviation for all variables\n");
    for ( i = 0; i < dim; i++ ) printf("MdAD%i  ", i + 1 );
    printf("\n");
    for ( i = 0; i < dim; i++ )
    {
        printf("%.3f  ", mdad[i] );
    }
    printf("\n");

    /***** Printing summary of the test *****/
    if (errnums == 0) {
        printf("\nComputed median absolute deviations");
        printf(" agree with theory.\n");
    }
    else {
        printf("\nError: Computed median absolute deviations");
        printf(" disagree with theory\n");
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return 0;
}
