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
*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Sparse BLAS C example
*
********************************************************************************
*
* Example program for using Intel MKL Sparse BLAS Level 2 and 3
* for matrices represented in the diagonal storage scheme.
* The following Sparse  Blas routines are used in the example: MKL_SDIAMM
*
* Consider the matrix A (see Appendix 'Sparse Storage Formats for Sparse Blas
* level 2-3')
*
*                 |   1       -1     -3    0     0   |
*                 |  -2        5      0    0     0   |
*   A    =        |   0        0      4    6     4   |,
*                 |  -4        0      2    7     0   |
*                 |   0        8      0    0    -5   |
*
*  The matrix A given above is represented in the diagonal storage scheme with the help of two
*  arrays (see Appendix 'Sparse Storage Formats for Sparse Blas level 2-3'
*
*              distance = (-3  -1  0   1   2  )
*                         | *   *  1  -1  -3  |
*                         | *  -2  5   0   0  |
*                values = | *   0  4   6   4  |
*                         |-4   2  7   0   *  |
*                         | 8   0 -5   *   *  |
*
*                 |   1   5   |
*                 |   1   4   |
*   S    =        |   1   3   |
*                 |   1   2   |
*                 |   1   1   |
*
*
*                 |  -3  - 8  |
*                 |   3   10  |
*   Sol  =        |  14   28  |
*                 |   5    0  |
*                 |   3   27  |
*
*
*
*  The test performs the following operations :
*
*       1. The next step is the computation A* S = F using MKL_SDIAMM where S is
*          a matrix 2*5. S is stored in C rows-oriented style.
*
* The code given below uses only one sparse representation for the all operations.*/

/*******************************************************************************
*     Definition arrays for sparse representation of  the matrix A in
*     the diagonal format:
*******************************************************************************/
#include <stdio.h>
#include "mkl_spblas.h"
#include <mkl_trans.h>
#include "mkl_types.h"

int main()
{
#define M 5 /* Number of rows in sparse matrix A */
#define N 2 /* Number of columns in matrix S */
#define NDIAG 5 /* Number of non-zero diagonal in matrix A */

    /* Matrix A defined by 2 arrays - values and distance */
    float  values[NDIAG*M] = {
        0.0,  0.0,  0.0, -4.0,  8.0,
        0.0, -2.0,  0.0,  2.0,  0.0,
        1.0,  5.0,  4.0,  7.0, -5.0,
        -1.0,  0.0,  6.0,  0.0,  0.0,
        -3.0,  0.0,  4.0,  0.0,  0.0
    };
    MKL_INT distance[NDIAG] = {-3, -1, 0,  1, 2};
    /* Matrix S. We should place matrix S column by column */
    float s[N*M] = {1.0, 5.0,
        1.0, 4.0,
        1.0, 3.0,
        1.0, 2.0,
        1.0, 1.0};

    /* Matrix H. We should place matrix H column by column */
    float h[N*M] = {- 3.0, - 8.0,
        3.0,  10.0,
        14.0,  28.0,
        5.0,   0.0,
        3.0,  27.0};

    /* Matrix Sol. We should place matrix Sol column by column */
    float sol[N*M] = {  0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0};

/*******************************************************************************
*    Declaration of local variables :
*******************************************************************************/
    MKL_INT i, j;
    MKL_INT n, m, ndiag;
    float alpha;    /* alpha */
    float beta;     /* beta */
    char *Ntrans;   /* Char[1] described transposed/non-transposed case */
    char *descra;   /* Char[6] described matrix properties */
    n = N;
    m = M;
    ndiag  = NDIAG;
    alpha  = 1.0;
    beta   = 0.0;
    Ntrans = "N";           /* Non-transposed case      */
    descra = "g_____";      /* General type of matrix A */
/*******************************************************************************
*    Task 1. Obtain matrix-matrix multiply  A*S
*    Array h must be equal to the array sol
*******************************************************************************/
    printf ("\n EXAMPLE PROGRAM FOR DIAGONAL FORMAT ROUTINES \n");
    printf ("\n     INPUT DATA FOR MKL_SDIAMM                \n");
    printf (  " Number of A rows     = %i\n", (int)m);
    printf (  " Number of S columns  = %i\n", (int)n);
    printf (  " alpha = %f\n", alpha);
    printf (  " beta  = %f\n", beta);
    printf (  " Input matrix                                 \n");

    /* Please pay attention: We store column by column to memory */
    printf (  " Input matrix                                 \n");
    for ( i = 0; i < M; i++ )
    {
        printf("| ");
        for ( j = 0; j < N; j++ )
        {
            printf(" %8.2f ", s[i * N + j]);
        }
        printf(" |\n");
    }
    /* Transpose matrix s */
    mkl_simatcopy(
                 'R', /* Row-major ordering */
                 'T', /* Matrix will be transposed*/
                 m,   /* Number of rows in source matrix */
                 n,   /* Number of columns in source matrix */
                 1.0, /* don't scale matrix */
                 s,   /* matrix which will be transposed*/
                 n,   /* source lda */
                 m    /* destination lda */
                 );
    /* Transpose matrix sol */
    mkl_simatcopy(
                 'R', /* Row-major ordering */
                 'T', /* Matrix will be transposed*/
                 m,   /* Number of rows in source matrix */
                 n,   /* Number of columns in source matrix */
                 1.0, /* don't scale matrix */
                 sol, /* matrix which will be transposed*/
                 n,   /* source lda */
                 m    /* destination lda */
                 );
    /* Compute sparse matrix-matrix multiplication */
    mkl_sdiamm(
              Ntrans, /* = 'N' we don't transpose matrix */
              &m,     /* = 5 Number of rows of sparse matrix */
              &n,     /* = 2 Number of columns in dense matrix C */
              &m,     /* = 5 Number of columns in sparse matrix */
              &alpha, /* = 1.0 */
              descra, /* = "G_____" */
              values, /* Array with non-zeros elements of sparse matrix A */
              &m,     /* Distance between diagonal in array values */
              distance, /* Array with distances between i diagonal and main diagonal */
              &ndiag, /* = 5 Number of diagonal with non-zeros elements */
              s,      /* Dense array S */
              &m,     /* Distance between columns of matrix S */
              &beta,  /* = 0.0 */
              sol,    /* Output array */
              &m      /* Distance between columns of matrix S */
              );
    /* Transpose matrix sol back */
    mkl_simatcopy(
                 'C', /* Column-major ordering */
                 'T', /* Matrix will be transposed*/
                 m,   /* Number of rows in source matrix */
                 n,   /* Number of columns in source matrix */
                 1.0, /* don't scale matrix */
                 sol, /* matrix which will be transposed*/
                 m,   /* source lda */
                 n    /* destination lda */
                 );
    printf ("\n     OUTPUT DATA FOR MKL_SDIAMM                \n");

    /* Please pay attention: We get column by column in memory */
    for ( i = 0; i < M; i++ )
    {
        printf("| ");
        for ( j = 0; j < N; j++ )
        {
            printf(" %8.2f ", sol[i * N + j]);
        }
        printf(" |\n");
    }
    printf ("\n     EXPECTED SOLUTION FOR MKL_SDIAMM           \n");
    for ( i = 0; i < M; i++ )
    {
        printf("| ");
        for ( j = 0; j < N; j++ )
        {
            printf(" %8.2f ", h[i * N + j]);
        }
        printf(" |\n");
    }

    return 0;
}
