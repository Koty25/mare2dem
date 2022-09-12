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
!    Parametrization of correlation matrix  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"

#define DIM   3          /* Task dimension */
#define LWORK 25*DIM

#define TEST_THRESHOLD -1.0E-6

/* Distorted correlation matrix */
double cor[DIM][DIM] = {
    {  1.0,  0.95,  0.7 },
    { 0.95,   1.0, 0.29 },
    {  0.7,  0.29,  1.0 }
};

int main()
{
    VSLSSTaskPtr task;
    MKL_INT dim;
    MKL_INT cor_storage;
    MKL_INT pcor_storage;
    double p_cor[DIM][DIM];
    double copy_cor[DIM][DIM];
    int i, j, errcode;
    int errnums = 0;

    /***** Following variables are used in routine which finds eigenvalues
           of simmetric matrix *****/
    double eigenvals[DIM], work[LWORK];
    MKL_INT lwork, info;
    char jobz, uplo;

    /***** Initializing parameters for Summary Statistics task *****/
    dim          = DIM;
    cor_storage  = VSL_SS_MATRIX_STORAGE_FULL;
    pcor_storage = VSL_SS_MATRIX_STORAGE_FULL;
    errcode      = 0;

    /***** Create Summary Statistics task *****/
    errcode = vsldSSNewTask( &task, &dim, 0, 0, 0, 0, 0 );
    CheckVslError(errcode);

    /***** Edit task parameters for parameterization of correlation *****/
    errcode = vsldSSEditCorParameterization( task,
              (double*)cor, &cor_storage, (double*)p_cor, &pcor_storage );
    CheckVslError(errcode);

    /***** Parametrize correlation *****/
    errcode = vsldSSCompute( task, VSL_SS_PARAMTR_COR,
                             VSL_SS_METHOD_SD );
    CheckVslError(errcode);

    /***** Compute eigenvalues of distorted correlation matrix *****/
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            copy_cor[i][j] = cor[i][j];
        }
    }

    lwork = LWORK;
    jobz = 'N';
    uplo = 'U';
    dsyev( &jobz, &uplo, &dim, (double*)copy_cor, &dim,
           eigenvals, work, &lwork, &info );
    CheckVslError(info);

    /***** Printing results *****/
    printf("Task dimension : %d\n", (int)dim);

    /***** Print distorted correlation matrix and it's eigenvalues *****/
    printf("\nDistorted correlation matrix\n");
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+1.5f ", cor[i][j]);
        }
        printf("\n");
    }

    printf("\nEigenvalues of the distorted correlation matrix\n");
    for(j = 0; j < dim; j++)
    {
        printf("%+1.5f ", eigenvals[j]);
    }
    printf("\n");

    /***** Compute eigevalue of parametrized correlation matrix *****/
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            copy_cor[i][j] = p_cor[i][j];
        }
    }

    lwork = LWORK;
    jobz = 'N';
    uplo = 'U';
    dsyev( &jobz, &uplo, &dim, (double*)copy_cor, &dim,
           eigenvals, work, &lwork, &info );
    CheckVslError(info);

    /***** Print parametrized correlation matrix and it's eigenvalues *****/
    printf("\nParameterized correlation matrix\n");
    for (i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            printf("%+1.5f ", p_cor[i][j]);
        }
        printf("\n");
    }

    printf("\nEigenvalues of the parameterized correlation matrix\n");
    for(j = 0; j < dim; j++)
    {
        printf("%+1.5f ", eigenvals[j]);
    }
    printf("\n\n");

    for(i = 0; i < dim; i++)
    {
        if (eigenvals[i] < TEST_THRESHOLD) errnums++;
    }

    if (errnums == 0)
    {
        printf("\nAll eigenvalues of parametrized correlation are in the");
        printf(" expected range\n");
    }
    else
    {
        printf("\nError: Parameterized correlation matrix has %i", errnums);
        printf(" negative eigenvalues\n");
        return 1;
    }

    /***** Delete Summary Statistics task *****/
    errcode = vslSSDeleteTask( &task );
    CheckVslError(errcode);

    MKL_Free_Buffers();

    return 0;
}
