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
!    vdRngGaussianMV  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include <math.h>

#include "mkl.h"
#include "errcheck.inc"
#include "statcheck.inc"

#define BRNG    VSL_BRNG_MCG31 // VSL basic generator to be used
#define SEED    7777777        // Initial value for stream initialization
#define NDIM    3              // Number of dimensions
#define N       10000          // Number of NDIM-dimensional vectors to generate
#define NN      5              // Number of NDIM-dimensional vectors to print

/*****************************************************
 Variance-covariance matrix for test
 (should be symmetric,positive-definite)
*****************************************************/
/* This is full storage for dpotrf subroutine */
static double C[NDIM][NDIM]={ 16.0, 8.0,  4.0,
                              8.0, 13.0, 17.0,
                              4.0, 17.0, 62.0  };

/*****************************************************
 Variance-covariance matrix for test
 (should be symmetric,positive-definite)
*****************************************************/
static double a[NDIM]={ 3.0, 5.0, 2.0 };
static double r[N][NDIM];

int main()
{
    VSLStreamStatePtr stream;

    double  CC[NDIM][NDIM]; // Copy of C
    int    i, j, errcode, status;

    /* Following variables are used in Cholesky factorization subroutine */
    char    uplo;
    MKL_INT n;
    double* T;
    MKL_INT lda;
    MKL_INT info;

    /* Sample characteristics */
    double dbMean[NDIM], dbVar[NDIM];
    double dbCovXY, dbCovXZ, dbCovYZ;

    double dbS[NDIM], dbS2[NDIM];

    double S[NDIM], D2[NDIM], Q[NDIM];
    double DeltaM[NDIM], DeltaD[NDIM];

    /*****************************************************
     Printing variance-covariance matrix and mean vector
     Also making the copy CC from C
    *****************************************************/
    printf("Variance-covariance matrix C:\n");
    for(i=0;i<NDIM;i++) {
        for(j=0;j<NDIM;j++) {
            CC[i][j] = C[i][j];
            printf("%4.1f ",C[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("Mean vector a:\n");
    for(i=0;i<NDIM;i++) {
        printf("%4.1f ",a[i]);
    }
    printf("\n\n");

    /*****************************************************
     Cholesky factorization -
     transform source symmetric matrix
     to lower triangular one
    *****************************************************/
    printf("VSL_MATRIX_STORAGE_FULL\n");
    printf("-----------------------\n\n");

    /* Full maxtrix storage (SPOTRF subroutine) */
    uplo = 'U'; // Due to FORTRAN stores matrices by columns while C by rows, ->
                // -> 'U' defines lower triangular matrix
    n=NDIM;
    lda=NDIM;
    T=(double*)C;

    /* Intel(R) MKL Choelsky factorization routine call */
    dpotrf( &uplo, &n, T, &lda, &info );

    /*****************************************************
     Stream initialization
    *****************************************************/
    errcode = vslNewStream( &stream, BRNG, SEED );
    CheckVslError( errcode );

    /*****************************************************
     Generating random numbers
     from multivariate normal distribution
    *****************************************************/
    errcode = vdRngGaussianMV( VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2, stream, N, (double *)r, NDIM, VSL_MATRIX_STORAGE_FULL, a, T );
    CheckVslError( errcode );

    /*****************************************************
     Printing random numbers
    *****************************************************/
    printf("Results (first %d of %d):\n", NN, N);
    printf("---------------------------\n");
    for(i=0;i<NN;i++) {
        printf("r[%d]=(",i);
        for(j=0;j<NDIM;j++) {
            printf("%5.2f ",r[i][j]);
        }
        printf(")\n");
    }
    printf("\n");

    /*****************************************************
     Testing sample characteristics
    *****************************************************/
    dCalculateGaussianMVSampleCharacteristics(N, NDIM, (double *)r,
            dbS, dbS2, dbMean, dbVar, &dbCovXY, &dbCovXZ, &dbCovYZ);

    /* Printing results */
    printf("Sample characteristics:\n");
    printf("-----------------------\n");
    printf("       Sample           Theory\n");
    printf("Mean :(%4.1f,%4.1f,%4.1f) (%4.1f,%4.1f,%4.1f)\n",dbMean[0],dbMean[1],dbMean[2],a[0],a[1],a[2]);
    printf("Var. :(%4.1f,%4.1f,%4.1f) (%4.1f,%4.1f,%4.1f)\n",dbVar[0],dbVar[1],dbVar[2],CC[0][0],CC[1][1],CC[2][2]);
    printf("CovXY: %4.1f               %4.1f\n",dbCovXY,CC[0][1]);
    printf("CovXZ: %4.1f               %4.1f\n",dbCovXZ,CC[0][2]);
    printf("CovYZ: %4.1f               %4.1f\n\n\n",dbCovYZ,CC[1][2]);

    /* Checking results */
    status = dGaussianMVCheckResults(NDIM, N, a, (double *)CC, dbMean, dbVar, S, D2, Q, DeltaM, DeltaD);

    if(status)
    {
        printf ("Error: sample moments\n");
        printf ("disagree with theory\n");
        printf ("     DeltaM: %7.3f %7.3f %7.3f \n", DeltaM[0], DeltaM[1], DeltaM[2]);
        printf ("     DeltaD: %7.3f %7.3f %7.3f \n", DeltaD[0], DeltaD[1], DeltaD[2]);
        printf ("     ( at least one of the Deltas > 3.0)\n");
        return 1;
    }
    else
    {
        printf ("Sample moments\n");
        printf ("agree with theory\n");
        printf ("     DeltaM: %7.3f %7.3f %7.3f \n", DeltaM[0], DeltaM[1], DeltaM[2]);
        printf ("     DeltaD: %7.3f %7.3f %7.3f \n", DeltaD[0], DeltaD[1], DeltaD[2]);
        printf ("     ( All Deltas < 3.0)\n");
    }

    /*****************************************************
     Stream finalization
    *****************************************************/
    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );

    return 0;
}
