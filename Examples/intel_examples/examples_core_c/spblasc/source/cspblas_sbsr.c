/*******************************************************************************
* Copyright 2005-2020 Intel Corporation.
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
!
!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Sparse BLAS C example
!
!*******************************************************************************
!
! Example program for using Intel MKL Sparse BLAS Level 2 and 3
! for matrices represented in the block compressed sparse row storage scheme.
! The following Sparse  Blas routines are used in the example:
!          MKL_SBSRSM  MKL_SBSRSV  MKL_SBSRMM  MKL_SBSRMV
!          MKL_SBSRGEMV    MKL_SBSRSYMV  MKL_SBSRTRSV.
!
! Consider the matrix A (see Appendix 'Sparse Storage Formats for Sparse Blas
! level 2-3')
!
!                 |   1   -1      0   -3   |
!                 |  -2    5      0    0   |   |Q  W |
!   A    =        |   0    0      4    6   | = |0  R |,
!                 |   0    0      2    7   |
!
! where
!           |  1  -1 |         |  0  -3 |        |  4   6 |
!     Q =   | -2   5 |,  W =   |  0   0 |, R =   |  2   7 |,
!
! decomposed as
!
!                      A = L + D + U,
!
!  where L is the strict lower triangle of A, U is the strictly upper triangle
!  of A, D is the main diagonal. Namely
!
!        |   0    0   0    0     0   |       |  0   -1    0   -3   0   |
!        |  -2    0   0    0     0   |       |  0    0    0    0   0   |
!   L  = |   0    0   0    0     0   |,  U=  |  0    0    0    6   4   |
!        |  -4    0   2    0     0   |       |  0    0    0    0   0   |
!        |   0    8   0    0     0   |       |  0    0    0    0   0   |
!
!
!           |   1  0  0   0   0   |
!           |   0  5  0   0   0   |
!   D    =  |   0  0  4   0   0   |.
!           |   0  0  0   7   0   |
!           |   0  0  0   0  -5   |
!
!
!
!  The matrix A is represented in the compressed block sparse row storage scheme
!  with the help of three arrays (see Appendix 'Sparse Matrix Storage') as follows:
!
!         values = (1 -2 -1  5  0  0 -3  0  4  2  6  7)
!         columns = (1 2 1)
!         rowIndex = (1  3  4)
!
!  It should be noted that two variations of the compressed sparse row storage scheme are supported by
!  Intel MKL Sparse Blas (see 'Sparse Storage Formats for Sparse Blas level 2-3') :
!
!        1. variation similar to the NIST Sparse Blas
!        2. variation similar to DSS/PARDISO, CXML and many other libraries.
!
!  The representation of the matrix A  given above is the 2-th variation. Two integer arrays
!  pointerB and pointerE instead of the array rowIndex are used in the NIST variation
!  of the block compressed sparse row format. Thus the arrays values and columns are the same for
!  the both variations. The arrays pointerB and pointerE for the matrix A are defined as follows:
!                          pointerB = (1 3)
!                          pointerE = (3 4)
!  It's easy to see that
!                    pointerB[i]= rowIndex[i] for i=0,1;
!                    pointerE[i]= rowIndex[i+1] for i=0,1.
!
!
!  The purpose of the given example is to show
!
!             1. how to call routines having interfaces suitable for the NIST's variation of the
!                block compressed sparse row format
!             2. how to form the arrays pointerB and pointerE for the NIST's variation of the
!                block compressed sparse row format using the  array rowIndex
!
!  In what follows the symbol ' means taking of transposed.
!
!  The test performs the following operations(Here all routine call
!  for zero-based modification of data) :
!
!       1. The code computes (L+D)'*S = F using MKL_SBSRMM where S is a known 5 by 2
!          matrix and then the code solves the system (L+D)'*X = F with the help of MKL_SBSRSM.
!          It's evident that X should be equal to S.
!
!       2. The code computes (U+D)'*S = F using MKL_SBSRMV where S is a vector
!          and then the code calls MKL_SBSRTRSV solves the system (U+I)'*X = F with the single right
!          hand side. It's evident that X should be equal to S.
!
!       3. The next step is the computation (U-U') S = F using MKL_SBSRMV where S is
!          a vector. It is easy to see that U-U' is a skew-symmetric matrix.
!
!       4. The next step is the computation (L+D+L') S = F using MKL_SBSRSYMV where S is
!          a vector. It is easy to see that L+D+L' is a symmetric matrix.
!
!       5. The next step is the computation A'* S = F using MKL_SBSRGEMV where S is
!          a vector.
!
!*******************************************************************************
*/
#include <stdio.h>
#include "mkl_types.h"
#include "mkl_spblas.h"

int main() {
//*******************************************************************************
//     Definition arrays for sparse representation of  the matrix A in
//     the compressed sparse row format:
//*******************************************************************************
#define M    2
#define NNZ  12
#define NNZB 3
#define LB   2
    MKL_INT m = M, lb = LB;
    float   values[NNZ]   = {1.0, -1.0, -2.0, 5.0, 0.0, -3.0, 0.0, 0.0, 4.0, 6.0, 2.0, 7.0};
    MKL_INT columns[NNZB] = {0, 1, 1};
    MKL_INT rowIndex[M+1] = {0, 2, 3};
    MKL_INT pointerB[M] = {0, 2}, pointerE[M] = {2, 3};
//*******************************************************************************
//    Declaration of local variables :
//*******************************************************************************
#define N 2
    MKL_INT     n = N;
    float       sol[M*LB][N]    = {{1.0, 4.0}, {1.0, 3.0}, {1.0, 2.0}, {1.0, 1.0}};
    float       rhs[M*LB][N]    = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    float       temp[M*LB][N]   = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    float       sol_vec[M*LB]   = {1.0, 1.0, 1.0, 1.0};
    float       rhs_vec[M*LB]   = {0.0, 0.0, 0.0, 0.0};
    float       temp_vec[M*LB]  = {0.0, 0.0, 0.0, 0.0};
    float       alpha = 1.0, beta = 0.0;
    MKL_INT     i, j;
    char        transa, uplo, nonunit;
    char        matdescra[6];

    printf("\n EXAMPLE PROGRAM FOR BLOCK COMPRESSED\n");
    printf("\n SPARSE ROW FORMAT ROUTINES          \n");
//*******************************************************************************
//Task 1.    Obtain matrix-matrix multiply (L+D)' *sol --> rhs
//    and solve triangular system   (L+D)'
//    *temp = rhs with multiple right hand sides
//    Array temp must be equal to the array sol
//*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_SBSRMM \n");
    printf("   WITH TRIANGULAR MATRIX    \n");
    printf("     M = %1.1i   N = %1.1i\n", m, n);
    printf("     ALPHA = %4.1f  BETA = %4.1f \n", alpha, beta);
    printf("     TRANS = '%c' \n", 'N');
    printf("   Input matrix              \n");
    for (i = 0; i < m * lb; i++) {
        for (j = 0; j < n; j++) {
            printf("%7.1f", sol[i][j]);
        };
        printf("\n");
    };

    transa = 't';
    matdescra[0] = 't';
    matdescra[1] = 'l';
    matdescra[2] = 'n';
    matdescra[3] = 'c';

    mkl_sbsrmm(&transa, &m, &n, &m, &lb, &alpha, matdescra, values, columns, rowIndex, &(rowIndex[1]), &(sol[0][0]), &n,  &beta, &(rhs[0][0]), &n);
    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_SBSRMM\n");
    printf("   WITH TRIANGULAR MATRIX    \n");
    for (i = 0; i < m * lb; i++) {
        for (j = 0; j < n; j++) {
            printf("%7.1f", rhs[i][j]);
        };
        printf("\n");
    };
    printf("-----------------------------------------------\n");
    printf("   Solve triangular system   \n");
    printf("   with obtained             \n");
    printf("   right hand side           \n");
    mkl_sbsrsm(&transa, &m, &n, &lb, &alpha, matdescra, values, columns, rowIndex, &(rowIndex[1]), &(rhs[0][0]), &n, &(temp[0][0]), &n);

    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_SBSRSM\n");
    for (i = 0; i < m * lb; i++) {
        for (j = 0; j < n; j++) {
            printf("%7.1f", temp[i][j]);
        };
        printf("\n");
    };
    printf("-----------------------------------------------\n");

//*******************************************************************************
// Task 2.    Obtain matrix-vector multiply (U+D)' *sol --> rhs
//    and solve triangular system   (U+D)' *temp = rhs with single right hand sides.
//    Array temp must be equal to the array sol.
//
//    Let us form the arrays pointerB and pointerE for the NIST's variation of the
//    compressed sparse row format using the  array rowIndex.
//
//*******************************************************************************
    for (j = 0; j < m; j++) {
        pointerB[j] = rowIndex[j];
        pointerE[j] = rowIndex[j+1];
    };
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_SBSRMV \n");
    printf("   WITH TRIANGULAR MATRIX    \n");
    printf("     ALPHA = %4.1f  BETA = %4.1f \n", alpha, beta);
    printf("     TRANS = '%c' \n", 'T');
    printf("   Input vector              \n");
    for (i = 0; i < m * lb; i++) {
        printf("%7.1f\n", sol_vec[i]);
    };

    transa = 't';
    matdescra[0] = 't';
    matdescra[1] = 'u';
    matdescra[2] = 'n';
    matdescra[3] = 'c';


    mkl_sbsrmv(&transa, &m, &m, &lb, &alpha, matdescra, values, columns, pointerB, pointerE, sol_vec, &beta, rhs_vec);
    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_SBSRMV\n");
    printf("   WITH TRIANGULAR MATRIX    \n");
    for (i = 0; i < m * lb; i++) {
        printf("%7.1f\n", rhs_vec[i]);
    };
    printf("-----------------------------------------------\n");
    printf("   Solve triangular system   \n");
    printf("   with obtained             \n");
    printf("   right hand side           \n");

    uplo        = 'u';
    nonunit     = 'n';
    mkl_cspblas_sbsrtrsv(&uplo, &transa, &nonunit, &m, &lb, values, rowIndex, columns, rhs_vec, temp_vec);

    printf("                             \n");
    printf("   OUTPUT DATA FOR           \n");
    printf("   MKL_CSPBLAS_SBSRTRSV      \n");
    printf("   WITH TRIANGULAR MATRIX    \n");

    for (i = 0; i < m*lb; i++) {
        printf("%7.1f\n", temp_vec[i]);
    };
    printf("-----------------------------------------------\n");
//*******************************************************************************
// Task 3.  Obtain matrix-vector multiply (U -U')*sol --> rhs
//    Array temp must be equal to the array sol
//*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_SBSRMV \n");
    printf("   WITH SKEW-SYMMETRIC MATRIX\n");
    printf("     ALPHA = %4.1f  BETA = %4.1f \n", alpha, beta);
    printf("     TRANS = '%c' \n", 'N');
    printf("   Input vector              \n");
    for (i = 0; i < m*lb; i++) {
        printf("%7.1f\n", sol_vec[i]);
    };

    transa = 'n';
    matdescra[0] = 'a';
    matdescra[1] = 'u';
    matdescra[2] = 'n';
    matdescra[3] = 'c';

    mkl_sbsrmv(&transa, &m, &m, &lb, &alpha, matdescra, values, columns, rowIndex, &(rowIndex[1]), sol_vec, &beta, rhs_vec);
    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_SCSRMV\n");
    printf("   WITH SKEW-SYMMETRIC MATRIX\n");
    for (i = 0; i < m * lb; i++) {
        printf("%7.1f\n", rhs_vec[i]);
    };
    printf("-----------------------------------------------\n");
//*******************************************************************************
// Task 4.    Obtain matrix-vector multiply (L+D+L')*sol --> rhs with the help of
//    MKL_SBSRSYMV
//*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR            \n");
    printf("   MKL_CSPBLAS_SBSRSYMV      \n");
    printf("   WITH SYMMETRIC MATRIX     \n");
    printf("     ALPHA = %4.1f  BETA = %4.1f \n", alpha, beta);
    printf("   Input vector              \n");
    for (i = 0; i < m*lb; i++) {
        printf("%7.1f\n", sol_vec[i]);
    };
    uplo = 'l';

    mkl_cspblas_sbsrsymv(&uplo, &m, &lb, values, rowIndex, columns, sol_vec, rhs_vec);
    printf("                             \n");
    printf("   OUTPUT DATA FOR           \n");
    printf("   MKL_CSPBLAS_SBSRGEMV      \n");
    printf("   WITH GENERAL MATRIX       \n");
    for (i = 0; i < m*lb; i++) {
        printf("%7.1f\n", rhs_vec[i]);
    };
    printf("-----------------------------------------------\n");

//*******************************************************************************
// Task 5.   Obtain matrix-vector multiply A'*sol --> rhs with the help of MKL_SBSRGEMV
//
//*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR            \n");
    printf("   MKL_CSPBLAS_SBSRGEMV      \n");
    printf("   WITH GENERAL MATRIX       \n");
    printf("     ALPHA = %4.1f  BETA = %4.1f \n", alpha, beta);
    printf("     TRANS = '%c' \n", 'T');
    printf("   Input vector              \n");
    for (i = 0; i < m * lb; i++) {
        printf("%7.1f\n", sol_vec[i]);
    };

    transa = 't';

    mkl_cspblas_sbsrgemv(&transa, &m, &lb, values, rowIndex, columns, sol_vec, rhs_vec);
    printf("                             \n");
    printf("   OUTPUT DATA FOR           \n");
    printf("   MKL_CSPBLAS_SBSRGEMV      \n");
    printf("   WITH GENERAL MATRIX       \n");
    for (i = 0; i < m * lb; i++) {
        printf("%7.1f\n", rhs_vec[i]);
    };
    printf("-----------------------------------------------\n");

    return 0;
}
