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
!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Sparse BLAS C example
!
!*******************************************************************************
!
! Example program for using Intel MKL Sparse BLAS Level 2 and 3
! for matrices represented in the coordinate storage scheme.
! The following Sparse  Blas routines are used in the example:
!          MKL_ZCOOSM  MKL_ZCOOSV  MKL_ZCOOMM  MKL_ZCOOMV
!          MKL_CSPBLAS_ZCOOGEMV    MKL_CSPBLAS_ZCOOSYMV
!                  MKL_CSPBLAS_ZCOOTRSV.
!
! Consider the matrix A (see Appendix 'Sparse Storage Formats for Sparse Blas
! level 2-3')
!
!                 |   1,1     -1,1   -3,1  0     0   |
!                 |  -2,1      5,1    0    0     0   |
!   A    =        |   0        0      4,1  6,1   4,1 |,
!                 |  -4,1      0      2,1  7,1   0   |
!                 |   0        8,1    0    0    -5,1 |
!
!
! decomposed as
!
!                      A = L + D + U,
!
!  where L is the strict  lower triangle of A, U is the strictly  upper triangle
!  of A, D is the main diagonal. Namely
!
!        |   0    0   0    0     0   |       |  0   -1,1 -3,1  0   0   |
!        |  -2,1  0   0    0     0   |       |  0    0    0    0   0   |
!   L  = |   0    0   0    0     0   |,  U=  |  0    0    0    6,1 4,1 |
!        |  -4,1  0   2,1  0     0   |       |  0    0    0    0   0   |
!        |   0    8   0    0     0   |       |  0    0    0    0   0   |
!
!
!           |   1,1 0   0   0   0   |
!           |   0   5,1 0   0   0   |
!   D    =  |   0   0   4,1 0   0   |.
!           |   0   0   0   7,1 0   |
!           |   0   0   0   0  -5,1 |
!
!  The matrix A given above is represented in the coordinate storage scheme with the help of three
!  arrays of length nnz=13 (see Appendix 'Sparse Storage Formats for Sparse Blas level 2-3'
!
!          values  = (1,1 -1,1 -3,1 -2,1 5,1 4,1 6,1 4,1 -4,1 2,1 7,1 8,1 -5,1)
!          rows    = (0  0  0  1 1 2 2 2  3 3 3 4  4)
!          columns = (0  1  2  0 1 2 3 4  0 2 3 1  4)
!
!  In what follows the symbol ' means taking of transposed.
!
!  The test performs the following operations :
!
!       1. The code computes (L+D)'*S = F using MKL_ZCOOMM where S is a known 5 by 2
!          matrix and then the code solves the system (L+D)'*X = F with the help of MKL_ZCOOSM.
!          It's evident that X should be equal to S.
!
!       2. The code computes (U+I)'*S = F using MKL_ZCOOMV where S is a vector
!          and then the code calls MKL_ZCOOTRSV solves the system (U+I)'*X = F with the single right
!          hand side. It's evident that X should be equal to S.
!
!       3. The code computes D*S = F using MKL_ZCOOMV where S is a vector
!          and then the code solves the system D*X = F with the single right hand side.
!          It's evident that X should be equal to S.
!
!       4. The next step is the computation (U-U') S = F using MKL_ZCOOMV where S is
!          a vector. It is easy to see that U-U' is a skew-symmetric matrix.
!
!       5. The next step is the computation (L+D+L') S = F using MKL_CSPBLAS_ZCOOSYMV
!                       where S is a vector. It is easy to see that L+D+L' is a symmetric matrix.
!
!       6. The next step is the computation A'* S = F using MKL_CSPBLAS_ZCOOGEMV
!                       where S is a vector.
!
! The code given below uses only one sparse representation for the all operations.
!*******************************************************************************
*/
#include <stdio.h>
#include "mkl_types.h"
#include "mkl_spblas.h"

int main () {
//*******************************************************************************
//     Definition arrays for sparse representation of  the matrix A in
//     the coordinate format:
//*******************************************************************************
#define M 5
#define NNZ 13
    MKL_INT         m = M, nnz = NNZ;
    MKL_Complex16   values[NNZ]   = { {1.0, 1.0}, {-1.0, 1.0}, {-3.0, 1.0}, {-2.0, 1.0}, {5.0, 1.0}, {4.0, 1.0}, {6.0, 1.0}, {4.0, 1.0}, {-4.0, 1.0}, {2.0, 1.0}, {7.0, 1.0}, {8.0, 1.0}, {-5.0, 1.0} };
    MKL_INT         columns[NNZ]  = {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4};
    MKL_INT         rows[NNZ]     = {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4};
//*******************************************************************************
//    Declaration of local variables :
//*******************************************************************************
#define N 2
    MKL_INT         n = N;
    MKL_Complex16   sol[M][N]       = { { {1.0, 1.0}, {5.0, 1.0} }, { {1.0, 1.0}, {4.0, 1.0} }, { {1.0, 1.0}, {3.0, 1.0} }, { {1.0, 1.0}, {2.0, 1.0} }, { {1.0, 1.0}, {1.0, 1.0} } };
    MKL_Complex16   rhs[M][N]       = { { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} } };
    MKL_Complex16   temp[M][N]      = { { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} } };
    MKL_Complex16   sol_vec[M]      = { {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0} };
    MKL_Complex16   rhs_vec[M]      = { {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0} };
    MKL_Complex16   temp_vec[M]     = { {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0} };
    MKL_Complex16   alpha = {1.0, 0.0}, beta = {0.0, 0.0};
    MKL_INT         i, j;
    char            transa, uplo, nonunit;
    char            matdescra[6];

    printf("\n EXAMPLE PROGRAM FOR COORDINATE FORMAT ROUTINES \n");

//*******************************************************************************
//    Task 1.    Obtain matrix-matrix multiply (L+D)' *sol --> rhs
//    and solve triangular system   (L+D)' *temp = rhs with multiple right hand sides
//    Array temp must be equal to the array sol
//*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_ZCOOMM \n");
    printf("   WITH TRIANGULAR MATRIX    \n");
    printf("     M = %1.1i   N = %1.1i\n", m, n);
    printf("     ALPHA = {%4.1f, %4.1f}  BETA = {%4.1f, %4.1f} \n", alpha.real, alpha.imag, beta.real, beta.imag);
    printf("     TRANS = '%c' \n", 'T');
    printf("   Input matrix              \n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("{%7.1f, %7.1f}   ", sol[i][j].real, sol[i][j].imag);
        };
        printf("\n");
    };

    transa = 't';
    matdescra[0] = 't';
    matdescra[1] = 'l';
    matdescra[2] = 'n';
    matdescra[3] = 'c';


    mkl_zcoomm(&transa, &m, &n, &m, &alpha, matdescra, values, rows, columns, &nnz, &(sol[0][0]), &n,  &beta, &(rhs[0][0]),  &n);

    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_ZCOOMM\n");
    printf("   WITH TRIANGULAR MATRIX    \n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("{%7.1f, %7.1f}   ", rhs[i][j].real, rhs[i][j].imag);
        };
        printf("\n");
    };
    printf("-----------------------------------------------\n");
    printf("   Solve triangular system   \n");
    printf("   with obtained             \n");
    printf("   right hand side           \n");
    mkl_zcoosm(&transa, &m, &n, &alpha, matdescra, values, rows, columns, &nnz, &(rhs[0][0]), &n, &(temp[0][0]), &n);

    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_ZCOOSM\n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("{%7.1f, %7.1f}   ", temp[i][j].real, temp[i][j].imag);
        };
        printf("\n");
    };
    printf("-----------------------------------------------\n");

//*******************************************************************************
//    Task 2.    Obtain matrix-vector multiply (U+I)' *sol --> rhs
//    and solve triangular system   (U+I)' *temp = rhs with single right hand sides
//    Array temp must be equal to the array sol
//*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_ZCOOMV \n");
    printf("   WITH TRIANGULAR MATRIX    \n");
    printf("     ALPHA = {%4.1f, %4.1f}  BETA = {%4.1f, %4.1f} \n", alpha.real, alpha.imag, beta.real, beta.imag);
    printf("     TRANS = '%c' \n", 'T');
    printf("   Input vector              \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", sol_vec[i].real, sol_vec[i].imag);
    };

    transa = 't';
    matdescra[0] = 't';
    matdescra[1] = 'u';
    matdescra[2] = 'u';
    matdescra[3] = 'c';

    mkl_zcoomv(&transa, &m, &m, &alpha, matdescra, values, rows, columns, &nnz, sol_vec, &beta, rhs_vec);

    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_ZCOOMV\n");
    printf("   WITH TRIANGULAR MATRIX    \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", rhs_vec[i].real, rhs_vec[i].imag);
    };
    printf("-----------------------------------------------\n");
    printf("   Solve triangular system   \n");
    printf("   with obtained             \n");
    printf("   right hand side           \n");

    uplo            = 'u';
    nonunit         = 'u';
    mkl_cspblas_zcootrsv(&uplo, &transa, &nonunit, &m, values, rows, columns, &nnz, rhs_vec, temp_vec);

    printf("                             \n");
    printf("   OUTPUT DATA FOR           \n");
    printf("   MKL_CSPBLAS_ZCOOTRSV      \n");
    printf("   WITH TRIANGULAR MATRIX    \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", temp_vec[i].real, temp_vec[i].imag);
    };
    printf("-----------------------------------------------\n");

//*******************************************************************************
//    Task 3.  Obtain matrix-vector multiply D *sol --> rhs
//    and solve triangular system   D *temp = rhs with single right hand side
//    Array temp must be equal to the array sol
//*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_ZCOOMV \n");
    printf("   WITH DIAGONAL MATRIX    \n");
    printf("     M = %1.1i   N = %1.1i\n", m, n);
    printf("     TRANS = '%c' \n", 'T');
    printf("   Input vector              \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", sol_vec[i].real, sol_vec[i].imag);
    };

    transa = 'n';
    matdescra[0] = 'd';
    matdescra[1] = 'u';
    matdescra[2] = 'n';
    matdescra[3] = 'c';

    mkl_zcoomv(&transa, &m, &m, &alpha, matdescra, values, rows, columns, &nnz, sol_vec, &beta, rhs_vec);

    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_ZCOOMV\n");
    printf("   WITH DIAGONAL MATRIX      \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", rhs_vec[i].real, rhs_vec[i].imag);
    };
    printf("-----------------------------------------------\n");
    printf("   Multiply by inverse      \n");
    printf("   matrix with the help     \n");
    printf("   of MKL_ZCSRSV            \n");

    mkl_zcoosv(&transa, &m, &alpha, matdescra, values, rows, columns, &nnz,  rhs_vec, temp_vec);

    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_ZCOOSV\n");
    printf("   WITH DIAGONAL MATRIX      \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", temp_vec[i].real, temp_vec[i].imag);
    };
    printf("-----------------------------------------------\n");


//*******************************************************************************
//    Task 4.  Obtain matrix-vector multiply (U -U')*sol --> rhs
//    Array temp must be equal to the array sol
//*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_ZCOOMV \n");
    printf("   WITH SKEW-SYMMETRIC MATRIX\n");
    printf("     ALPHA = {%4.1f, %4.1f}  BETA = {%4.1f, %4.1f} \n", alpha.real, alpha.imag, beta.real, beta.imag);
    printf("     TRANS = '%c' \n", 'N');
    printf("   Input vector              \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", sol_vec[i].real, sol_vec[i].imag);
    };

    transa = 'n';
    matdescra[0] = 'a';
    matdescra[1] = 'u';
    matdescra[2] = 'n';
    matdescra[3] = 'c';

    mkl_zcoomv(&transa, &m, &m, &alpha, matdescra, values, rows, columns, &nnz, sol_vec, &beta, rhs_vec);

    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_ZCOOMV\n");
    printf("   WITH SKEW-SYMMETRIC MATRIX\n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", rhs_vec[i].real, rhs_vec[i].imag);
    };
    printf("-----------------------------------------------\n");
//*******************************************************************************
//    Task 5.  Obtain matrix-vector multiply (L+D+L')*sol --> rhs with the help of
//    MKL_CSPBLAS_ZCOOSYMV
//
//*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR            \n");
    printf("   MKL_CSPBLAS_ZCOOSYMV      \n");
    printf("   WITH SYMMETRIC MATRIX     \n");
    printf("     ALPHA = {%4.1f, %4.1f}  BETA = {%4.1f, %4.1f} \n", alpha.real, alpha.imag, beta.real, beta.imag);
    printf("   Input vector              \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", sol_vec[i].real, sol_vec[i].imag);
    };
    uplo = 'l';

    mkl_cspblas_zcoosymv(&uplo, &m, values, rows, columns, &nnz, sol_vec, rhs_vec);

    printf("                             \n");
    printf("   OUTPUT DATA FOR           \n");
    printf("   MKL_CSPBLAS_ZCOOSYMV      \n");
    printf("   WITH SYMMETRIC MATRIX     \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", rhs_vec[i].real, rhs_vec[i].imag);
    };
    printf("-----------------------------------------------\n");

//*******************************************************************************
//    Task 6. Obtain matrix-vector multiply A'*sol --> rhs with the help
//    of MKL_CSPBLAS_ZCOOGEMV
//
//*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR            \n");
    printf("   MKL_CSPBLAS_ZCOOGEMV      \n");
    printf("   WITH GENERAL MATRIX       \n");
    printf("     ALPHA = {%4.1f, %4.1f}  BETA = {%4.1f, %4.1f} \n", alpha.real, alpha.imag, beta.real, beta.imag);
    printf("     TRANS = '%c' \n", 'T');
    printf("   Input vector              \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", sol_vec[i].real, sol_vec[i].imag);
    };

    transa = 't';

    mkl_cspblas_zcoogemv(&transa, &m, values, rows, columns, &nnz, sol_vec, rhs_vec);

    printf("                             \n");
    printf("   OUTPUT DATA FOR           \n");
    printf("   MKL_CSPBLAS_ZCOOGEMV      \n");
    printf("   WITH GENERAL MATRIX       \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", rhs_vec[i].real, rhs_vec[i].imag);
    };
    printf("-----------------------------------------------\n");

    return 0;
}
