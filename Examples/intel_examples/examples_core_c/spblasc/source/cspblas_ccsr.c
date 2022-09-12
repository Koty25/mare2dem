/*******************************************************************************
* Copyright 2006-2020 Intel Corporation.
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
! for matrices represented in the compressed sparse row storage scheme.
! The following Sparse  Blas routines are used in the example:
!          MKL_CCSRSM  MKL_CCSRSV  MKL_CCSRMM  MKL_CCSRMV
!          MKL_CSPBLAS_CCSRGEMV    MKL_CSPBLAS_CCSRSYMV
!                  MKL_CSPBLAS_CCSRTRSV.
!
! Consider the matrix A (see Appendix 'Sparse Storage Formats for Sparse Blas
! level 2-3')
!
!                 |   1,1     -1,1    0   -3,1   0   |
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
!  where L is the strict lower triangle of A, U is the strictly  upper triangle
!  of A, D is the main diagonal. Namely
!
!        |   0    0   0    0     0   |       |  0   -1,1  0   -3,1 0   |
!        |  -2,1  0   0    0     0   |       |  0    0    0    0   0   |
!   L  = |   0    0   0    0     0   |,  U=  |  0    0    0    6,1 4,1 |
!        |  -4,1  0   2,1  0     0   |       |  0    0    0    0   0   |
!        |   0    8,1 0    0     0   |       |  0    0    0    0   0   |
!
!
!           |   1,1  0    0   0   0   |
!           |   0    5,1  0   0   0   |
!   D    =  |   0    0    4,1 0   0   |.
!           |   0    0    0   7,1 0   |
!           |   0    0    0   0  -5,1 |
!
!  The matrix A is represented in the zero-based compressed sparse row storage
!  scheme with the help of three arrays  (see Appendix 'Sparse Matrix Storage')
!  as follows:
!
!         values = (1,1 -1,1 -3,1 -2,1 5,1 4,1 6,1 4,1 -4,1 2,1 7,1 8,1 -5,1)
!         columns = (0 1 3 0 1 2 3 4 0 2 3 1 4)
!         rowIndex = (0  3  5  8  11 13)
!
!  It should be noted that two variations of the compressed sparse row storage
!  scheme are supported by Intel MKL Sparse Blas (see 'Sparse Storage Formats
!  for Sparse Blas level 2-3') :
!
!        1. variation accepted in the NIST Sparse Blas (zero-based modification)
!        2. variation accepted for many other libraries (zero-based documentation).
!
!  The representation of the matrix A  given above is the 2's variation. Two
!  integer arrays pointerB and pointerE instead of the array rowIndex are used in
!  the NIST variation of variation of the compressed sparse row format. Thus the
!  arrays values and columns are the same for the both variations. The arrays
!  pointerB and pointerE for the matrix A are defined as follows:
!                          pointerB = (0 3  5  8 11)
!                          pointerE = (3 5  8 11 13)
!  It's easy to see that
!                    pointerB[i]= rowIndex[i] for i=0, ..4;
!                    pointerE[i]= rowIndex[i+1] for i=0, ..4.
!
!
!  The purpose of the given example is to show
!
!             1. how to form the arrays pointerB and pointerE for the NIST's
!                                variation of the  compressed sparse row format using
!                the  array rowIndex
!             2. how to use minors of the matrix A by redefining the arrays
!                                pointerB and pointerE but the arrays values and columns are
!                                the same.
!
!  In what follows the symbol ' means taking of transposed.
!
!  The test performs the following operations :
!
!       1. The code computes (L+D)'*S = F using MKL_CCSRMM where S is a known
!                  5 by 2 matrix and then the code solves the system (L+D)'*X = F with
!                  the help of MKL_CCSRSM. It's evident that X should be equal to S.
!
!       2. The code computes (U+I)'*S = F using MKL_CCSRMV where S is a vector
!          and then the code calls MKL_CSPBLAS_CCSRTRSV solves the system
!                  (U+I)'*X = F with the single right hand side. It's evident
!                       that X should be equal to S.
!
!       3. The code computes D*S = F using MKL_CCSRMV where S is a vector
!          and then the code solves the system D*X = F with the single
!                  right hand side. It's evident that X should be equal to S.
!
!       4. The next step is the computation (U-U') S = F using MKL_CCSRMV where
!          S is a vector. It is easy to see that U-U' is a skew-symmetric matrix.
!
!       5. The next step is the computation (L+D+L') S = F using
!                  MKL_CSPBLAS_CCSRSYMV where S is a vector.
!          The vector is computed two times. At first, the sparse representation
!          of the whole matrix A is used. Then the vector is computed with the help of
!          sparse representation of L+D. These two calls must give the same vector.
!
!       6. The next step is the computation A'* S = F using MKL_CSPBLAS_CCSRGEMV
!                  where S is a vector.
!
!       7. Let's T be the upper 3 by 3 minor of the matrix A. Namely, T is
!                  the following matrix
!
!                        |   1       -1      0   |
!          T    =        |  -2        5      0   |.
!                        |   0        0      4   |
!          The test performs the matrix-vector multiply T*S=F with the same
!                  arrays values, columns and pointerB used before for the whole
!                  matrix A. It is enough to change two values of array pointerE
!                  in order to use the minor under consideration. Then the test
!                  solves the system T*X =F using MKL_CCSRSV. The routine MKL_CCSRMV
!                  is used for getting matrix-vector multiply.
!
! The code given below uses only one sparse representation for the all operations.
!
!*******************************************************************************
*/
#include <stdio.h>
#include "mkl_types.h"
#include "mkl_spblas.h"

int main() {
//*******************************************************************************
//     Declaration and initialization of parameters for sparse representation of
//     the matrix A in  the compressed sparse row format:
//*******************************************************************************
#define M 5
#define NNZ 13
#define NNZ1 9
#define MNEW 3
    MKL_INT         m = M, mnew = MNEW;
    //*******************************************************************************
    //    Sparse representation of the matrix A
    //*******************************************************************************
    MKL_Complex8    values[NNZ]   = { {1.0, 1.0}, {-1.0, 1.0}, {-3.0, 1.0}, {-2.0, 1.0}, {5.0, 1.0}, {4.0, 1.0}, {6.0, 1.0}, {4.0, 1.0}, {-4.0, 1.0}, {2.0, 1.0}, {7.0, 1.0}, {8.0, 1.0}, {-5.0, 1.0} };
    MKL_INT         columns[NNZ]  = {0, 1, 3, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4};
    MKL_INT         rowIndex[M+1] = {0, 3,  5,  8,  11, 13};
    //*******************************************************************************
    //    Sparse representation of the lower triangle L+D
    //*******************************************************************************

    MKL_Complex8   values1[NNZ1]  = { {1.0, 1.0}, {-2.0, 1.0}, {5.0,1.0}, {4.0,1.0}, {-4.0,1.0}, {2.0,1.0}, {7.0,1.0}, {8.0,1.0}, {-5.0, 1.0} };
    MKL_INT        columns1[NNZ1] = {0,  0, 1, 2, 0, 2, 3, 1, 4};
    MKL_INT        rowIndex1[M+1] = {0, 1, 3, 4, 7, 9};
    //*******************************************************************************
    //    Declaration of local variables :
    //*******************************************************************************
    #define N 2
    MKL_INT         n = N;
    MKL_Complex8    sol[M][N]       = { { {1.0, 1.0}, {5.0, 1.0} }, { {1.0, 1.0}, {4.0, 1.0} }, { {1.0, 1.0}, {3.0, 1.0} }, { {1.0, 1.0}, {2.0, 1.0} }, { {1.0, 1.0}, {1.0, 1.0} } };
    MKL_Complex8    rhs[M][N]       = { { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} } };
    MKL_Complex8    temp[M][N]      = { { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} }, { {0.0, 1.0}, {0.0, 1.0} } };
    MKL_Complex8    sol_vec[M]      = { {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0} };
    MKL_Complex8    rhs_vec[M]      = { {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0} };
    MKL_Complex8    temp_vec[M]     = { {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0} };
    MKL_Complex8    alpha = {1.0, 0.0}, beta = {0.0, 0.0};
    MKL_INT         i, j;
    MKL_INT         pointerB[M] , pointerE[M];
    char            transa, uplo, nonunit;
    char            matdescra[6];

    printf("\n EXAMPLE PROGRAM FOR COMPRESSED SPARSE ROW FORMAT ROUTINES \n");
    //*******************************************************************************
    //Task 1.    Obtain matrix-matrix multiply (L+D)' *sol --> rhs
    //    and solve triangular system   (L+D)'
    //        *temp = rhs with multiple right hand sides
    //    Array temp must be equal to the array sol
    //*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_CCSRMM \n");
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

    mkl_ccsrmm(&transa, &m, &n, &m, &alpha, matdescra, values, columns, rowIndex, &(rowIndex[1]), &(sol[0][0]), &n,  &beta, &(rhs[0][0]), &n);

    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_CCSRMM\n");
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
    mkl_ccsrsm(&transa, &m, &n, &alpha, matdescra, values, columns, rowIndex, &(rowIndex[1]), &(rhs[0][0]), &n, &(temp[0][0]), &n);

    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_CCSRSM\n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("{%7.1f, %7.1f}   ", temp[i][j].real, temp[i][j].imag);
        };
        printf("\n");
    };
    printf("-----------------------------------------------\n");

    //*******************************************************************************
    // Task 2.    Obtain matrix-vector multiply (U+I)' *sol --> rhs
    //    and solve triangular system   (U+I)' *temp = rhs with single
    //    right hand sides. Array temp must be equal to the array sol.
    //
    //    Let us form the arrays pointerB and pointerE for the NIST's
    //        variation of the compressed sparse row format using the array
    //    rowIndex.
    //
    //*******************************************************************************
    for (i = 0; i < m; i++) {
        pointerB[i] = rowIndex[i];
        pointerE[i] = rowIndex[i+1];
    };
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_CCSRMV \n");
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

    mkl_ccsrmv(&transa, &m, &m, &alpha, matdescra, values, columns, pointerB, pointerE, sol_vec, &beta, rhs_vec);

    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_CCSRMV\n");
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
    mkl_cspblas_ccsrtrsv(&uplo, &transa, &nonunit, &m, values, rowIndex, columns, rhs_vec, temp_vec);

    printf("                             \n");
    printf("   OUTPUT DATA FOR           \n");
    printf("   MKL_CSPBLAS_CCSRTRSV      \n");
    printf("   WITH TRIANGULAR MATRIX    \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", temp_vec[i].real, temp_vec[i].imag);
    };
    printf("-----------------------------------------------\n");
    //*******************************************************************************
    // Task 3.   Obtain matrix-vector multiply D *sol --> rhs
    //    and solve triangular system   D *temp = rhs with single right hand side
    //    Array temp must be equal to the array sol
    //*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_CCSRMV \n");
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

    mkl_ccsrmv(&transa, &m, &m, &alpha, matdescra, values, columns, pointerB, pointerE, sol_vec, &beta, rhs_vec);
    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_CCSRMV\n");
    printf("   WITH DIAGONAL MATRIX      \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", rhs_vec[i].real, rhs_vec[i].imag);
    };
    printf("-----------------------------------------------\n");
    printf("   Multiply by inverse      \n");
    printf("   matrix with the help     \n");
    printf("   of MKL_CCSRSV            \n");


    mkl_ccsrsv(&transa, &m, &alpha, matdescra, values, columns, pointerB, pointerE, rhs_vec, temp_vec);
    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_CCSRSV\n");
    printf("   WITH DIAGONAL MATRIX      \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", temp_vec[i].real, temp_vec[i].imag);
    };
    printf("-----------------------------------------------\n");

    //*******************************************************************************
    // Task 4.  Obtain matrix-vector multiply (U -U')*sol --> rhs
    //    Array temp must be equal to the array sol
    //*******************************************************************************
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_CCSRMV \n");
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
    mkl_ccsrmv(&transa, &m, &m, &alpha, matdescra, values, columns, rowIndex, &(rowIndex[1]), sol_vec, &beta, rhs_vec);
    printf("                             \n");
    printf("   OUTPUT DATA FOR MKL_CCSRMV\n");
    printf("   WITH SKEW-SYMMETRIC MATRIX\n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", rhs_vec[i].real, rhs_vec[i].imag);
    };
    printf("-----------------------------------------------\n");
    //*******************************************************************************
    // Task 5.    Obtain matrix-vector multiply (L+D+L')*sol --> rhs with
    //                              the help of mkl_cspblas_ccsrsymv
    // NOTE: The routine mkl_cspblas_ccsrsymv as well as the similar dense Level 2
    // routine csymv has a possibilty to extract the required triangle from the input
    // sparse matrix and perform symmetric matrix-vector multiply with the help of
    // the triangle. Let the arrays values, rowIndex, columns be the sparse
    // representation of A
    //
    //
    //                |   (1,1)   (-1,1)    0   (-3,1)   0   |
    //                |  (-2,1)    (5,1)    0      0     0   |
    //  A    =        |     0        0    (4,1)  (6,1) (4,1) |,
    //                |  (-4,1)      0    (2,1)  (7,1)   0   |
    //                |     0      (8,1)    0      0  (-5,1) |
    //
    // Let the arrays values1, rowIndex1, columns1 which are the following
    //
    //   values1[NNZ1]={1.0, 1.0, -2.0,1.0, 5.0,1.0,  4.0,1.0,   -4.0,1.0,
    //           2.0,1.0,  7.0,1.0,   8.0,1.0,  -5.0, 1.0};
    //   columns1[9]={0,  0, 1, 2, 0, 2, 3, 1, 4};
    //   rowIndex1[6]={0, 1, 3, 4, 7, 9};
    // be the sparse representation of the lower triangle of A
    //
    //           |   (1,1)      0      0     0      0   |
    //           |  (-2,1)    (5,1)    0     0      0   |
    //   L +D  = |     0        0    (4,1)   0      0   |,
    //           |  (-4,1)      0    (2,1) (7,1)    0   |
    //           |     0      (8,1)    0     0    (-5,1)|
    //
    //  The feature described above means that  the following two calls must give the
    //  same output vector
    //  mkl_cspblas_ccsrsymv(&uplo, &m, values, rowIndex, columns, sol_vec, rhs_vec);
    //  mkl_cspblas_ccsrsymv(&uplo, &m, values1, rowIndex1, columns1, sol_vec, temp_vec);
    //
    //  The next test checks whether these two calls give the same output vector.
    //*******************************************************************************
    printf("                             \n");
    printf("   SYMMETRIC MATRIX-VECTOR MULTIPLY ROUTINE\n");
    printf("   INPUT DATA FOR            \n");
    printf("   MKL_CSPBLAS_CCSRSYMV      \n");
    printf("   Input vector              \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", sol_vec[i].real, sol_vec[i].imag);
    };
    uplo = 'l';

    mkl_cspblas_ccsrsymv(&uplo, &m, values, rowIndex, columns, sol_vec, rhs_vec);
    mkl_cspblas_ccsrsymv(&uplo, &m, values1, rowIndex1, columns1, sol_vec, temp_vec);

    printf("                             \n");
    printf(" OUTPUT VECTOR FOR SYMMETRIC MATRIX-VECTOR MULTIPLY ROUTINE\n");
    printf("   MKL_CSPBLAS_CCSRSYMV  CALLED TWO TIMES    \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}  {%7.1f, %7.1f} \n", rhs_vec[i].real, rhs_vec[i].imag, temp_vec[i].real, temp_vec[i].imag);
    };
    printf("-----------------------------------------------\n");
    //*******************************************************************************
    // Task 6.   Obtain matrix-vector multiply A'*sol --> rhs with
    //                              the help of MKL_CSPBLAS_CCSRGEMV
    //*******************************************************************************
    printf("                             \n");
    printf("   MATRIX-VECTOR MULTIPLY ROUTINE FOR GENERAL MATRICES\n");
    printf("   INPUT DATA FOR            \n");
    printf("   MKL_CSPBLAS_CCSRGEMV      \n");
    printf("     TRANS = '%c' \n", 'T');
    printf("   Input vector              \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", sol_vec[i].real, sol_vec[i].imag);
    };

    transa = 't';

    mkl_cspblas_ccsrgemv(&transa, &m, values, rowIndex, columns, sol_vec, rhs_vec);
    printf("                             \n");
    printf(" MATRIX-VECTOR MULTIPLY ROUTINE FOR GENERAL MATRICES\n");
    printf("   OUTPUT DATA FOR           \n");
    printf("   MKL_CSPBLAS_CCSRGEMV      \n");
    for (i = 0; i < m; i++) {
        printf("{%7.1f, %7.1f}\n", rhs_vec[i].real, rhs_vec[i].imag);
    };
    printf("-----------------------------------------------\n");

    //*******************************************************************************
    // Task 7.  Obtain matrix-vector multiply T*sol --> rhs with the help of
    //        MKL_CCSRMV where S is  3 by 3 minor of the matrix A starting with A(1,1)
    //    Let's us redefine two elements of the array pointerE in order to identify
    //    the needed minor. More precisely
    //            pointerE(1) --> pointerE(1)-1
    //            pointerE(3) --> pointerE(3)-2
    //
    //*******************************************************************************
    pointerE[0] = pointerE[0] - 1;
    pointerE[2] = pointerE[2] - 2;
    printf("                             \n");
    printf("   INPUT DATA FOR MKL_CCSRMV \n");
    printf("   WITH A MINOR OF GENERAL   \n");
    printf("   MATRIX                    \n");
    printf("     ALPHA = {%4.1f, %4.1f}  BETA = {%4.1f, %4.1f} \n", alpha.real, alpha.imag, beta.real, beta.imag);
    printf("     TRANS = '%c' \n", 'T');
    printf("   Input vector              \n");
    for (i = 0; i < mnew; i++) {
        printf("{%7.1f, %7.1f}\n", sol_vec[i].real, sol_vec[i].imag);
    };

    transa = 'n';
    matdescra[0] = 't';
    matdescra[1] = 'l';
    matdescra[2] = 'n';
    matdescra[3] = 'c';
    mkl_ccsrmv(&transa, &mnew, &mnew, &alpha, matdescra, values, columns, pointerB, pointerE, sol_vec, &beta, rhs_vec);

    printf("                                       \n");
    printf("   OUTPUT DATA FOR MKL_CCSRMV          \n");
    printf("   WITH A MINOR OF GENERAL MATRIX      \n");
    for (i = 0; i < mnew; i++) {
        printf("{%7.1f, %7.1f}\n", rhs_vec[i].real, rhs_vec[i].imag);
    };
    printf("-----------------------------------------------\n");

    printf("    Multiply by inverse to a    \n");
    printf("    minor of the matrix with    \n");
    printf("    the help of MKL_CCSRSV      \n");

    mkl_ccsrsv(&transa, &mnew, &alpha, matdescra, values, columns, pointerB, pointerE, rhs_vec, temp_vec);

    printf("                                 \n");
    printf("   OUTPUT DATA FOR MKL_CCSRSV    \n");
    printf("   WITH A MINOR OF GENERAL MATRIX\n");
    for (i = 0; i < mnew; i++) {
        printf("{%7.1f, %7.1f}\n", temp_vec[i].real, temp_vec[i].imag);
    };
    printf("-----------------------------------------------\n");
    return 0;
}
