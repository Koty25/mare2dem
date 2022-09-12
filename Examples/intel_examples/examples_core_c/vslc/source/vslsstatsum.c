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
!    Calculation of raw/central sums Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "statchars.inc"

#define DIM     4        /* Task dimension */
#define N       1000     /* Number of observations */

#define P_THRESHOLD     0.01

float C[DIM][DIM] = {
    { 1.0, 0.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0, 0.0 },
    { 0.0, 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 0.0, 1.0 }
};
float a[DIM] = { 5.0, 5.0, 5.0, 5.0 };

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim;
    MKL_INT n;
    MKL_INT x_storage;
    float x[N][DIM];

    float sum[DIM], mean[DIM];
    float rs2[DIM], rs3[DIM], rs4[DIM];
    float rm2[DIM], rm3[DIM], rm4[DIM];

    float cs2[DIM], cs3[DIM], cs4[DIM];
    float cm2[DIM], cm3[DIM], cm4[DIM];

    int i, errcode;
    unsigned MKL_INT64 estimate = 0;
    int errnums = 0;

    float pval_mean[DIM];
    float pval_rm2[DIM], pval_rm3[DIM], pval_rm4[DIM];
    float pval_cm2[DIM], pval_cm3[DIM], pval_cm4[DIM];

    /***** Initializing parameters for Summary Statistics task *****/
    dim         = DIM;
    n           = N;
    x_storage   = VSL_SS_MATRIX_STORAGE_COLS;

    /***** Generate data set using VSL GaussianMV RNG *****/
    errcode = sGenerateGaussianMVData( (float*)x, dim, n, a, (float*)C );
    CheckVslError(errcode);

    /***** Create Summary Statistics task *****/
    errcode = vslsSSNewTask( &task, &dim, &n, &x_storage, (float*)x, 0, 0 );
    CheckVslError(errcode);

    /***** Edit task parameters for mean to be used in sum
          computations *****/
    errcode = vslsSSEditTask( task, VSL_SS_ED_MEAN, mean );

    /***** Edit task parameters for computation of sum, 2nd, 3rd
           and 4th raw and central sums *****/
    errcode = vslsSSEditSums( task, sum, rs2, rs3, rs4, cs2, cs3, cs4 );
    CheckVslError(errcode);

    /***** Sum, 2nd, 3rd and 4th raw and central sums are included
           in the list of estimates to compute *****/
    estimate |= VSL_SS_SUM |
        VSL_SS_2R_SUM | VSL_SS_3R_SUM | VSL_SS_4R_SUM |
        VSL_SS_2C_SUM | VSL_SS_3C_SUM | VSL_SS_4C_SUM;

    /***** Compute the sums using 1PASS method *****/
    errcode = vslsSSCompute( task, estimate, VSL_SS_METHOD_1PASS );
    CheckVslError(errcode);

    /***** Edit task parameters for computation of mean, 2nd, 3rd
           and 4th raw and central moments *****/
    errcode = vslsSSEditMoments( task, mean, rm2, rm3, rm4, cm2, cm3, cm4 );
    CheckVslError(errcode);

    /***** Convert sums into moments *****/
    estimate = VSL_SS_MEAN |
        VSL_SS_2R_MOM | VSL_SS_3R_MOM | VSL_SS_4R_MOM |
        VSL_SS_2C_MOM | VSL_SS_3C_MOM | VSL_SS_4C_MOM;

    errcode = vslsSSCompute( task, estimate, VSL_SS_METHOD_SUM_TO_MOM );
    CheckVslError(errcode);

    /***** Testing stat characteristics of the computed estimates *****/

    /* Compute p-values for sum estimates */
    sComputePvalsMean( dim, n, mean, a, (float*)C, pval_mean );

    /* Compute p-values for raw moments estimates */
    sComputePvalsRawMoments( dim, n, rm2, 2, a, (float*)C, pval_rm2 );
    sComputePvalsRawMoments( dim, n, rm3, 3, a, (float*)C, pval_rm3 );
    sComputePvalsRawMoments( dim, n, rm4, 4, a, (float*)C, pval_rm4 );
    /* Compute p-values for central moments estimates */
    sComputePvalsCentralMoments( dim, n, cm2, 2, a, (float*)C,pval_cm2 );
    sComputePvalsCentralMoments( dim, n, cm3, 3, a, (float*)C,pval_cm3 );
    sComputePvalsCentralMoments( dim, n, cm4, 4, a, (float*)C,pval_cm4 );

    /***** Checking the validity of p-values for all estimates *****/
    for(i = 0; i < dim; i++)
    {
        if (pval_mean[i] < P_THRESHOLD) errnums++;
        if (pval_rm2 [i] < P_THRESHOLD) errnums++;
        if (pval_rm3 [i] < P_THRESHOLD) errnums++;
        if (pval_rm4 [i] < P_THRESHOLD) errnums++;
        if (pval_cm2 [i] < P_THRESHOLD) errnums++;
        if (pval_cm3 [i] < P_THRESHOLD) errnums++;
        if (pval_cm4 [i] < P_THRESHOLD) errnums++;
    }

    /***** Printing results *****/
    printf("Task dimension : %d\n", (int)dim);
    printf("Number of observations : %d\n\n", (int)n);

    /***** Printing computed sum estimates *****/
    printf("                   Sum            2RS            3RS            ");
    printf("4RS            2CS         3CS          4RS\n");
    for(i = 0; i < dim; i++)
    {
        printf("Variable #%i:  %+f  %+f  %+f  %+f", i + 1,
                 sum[i],  rs2[i], rs3[i], rs4[i]);
        printf("  %+f  %+f  %+f\n", cs2[i], cs3[i], cs4[i]);
    }

    /***** Printing p-values for moments obtained from sums *****/
    printf("\n\nP-values of the moments computed from sums\n\n");
    printf("                 MEAN       2RM  ");
    printf("      3RM        4RM        2CM        3CM        4CM\n");
    for(i = 0; i < dim; i++)
    {
        printf("Variable #%i:  %+f  %+f  %+f  %+f  %+f  %+f  %+f\n",
            i + 1, pval_mean[i], pval_rm2[i], pval_rm3[i], pval_rm4[i],
            pval_cm2[i], pval_cm3[i], pval_cm4[i]);
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
