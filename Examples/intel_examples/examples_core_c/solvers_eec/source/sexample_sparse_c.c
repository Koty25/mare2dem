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
!   Content: Example for Intel(R) Math Kernel Library (Intel(R) MKL) Extended
!            Eigensolvers (sparse format, single precision)
!
!*******************************************************************************
!
! The following routines are used in the example:
!          SGEMM  SFEAST_SCSREV  SFEAST_SCSRGV  FEASTINIT
!
! Consider the matrix A
!                 |  5   2   1   1   0   0   0   0   0   0   0  |
!                 |  2   6   3   1   1   0   0   0   0   0   0  |
!                 |  1   3   6   3   1   1   0   0   0   0   0  |
!                 |  1   1   3   6   3   1   1   0   0   0   0  |
!                 |  0   1   1   3   6   3   1   1   0   0   0  |
!    A    =       |  0   0   1   1   3   6   3   1   1   0   0  |,
!                 |  0   0   0   1   1   3   6   3   1   1   0  |
!                 |  0   0   0   0   1   1   3   6   3   1   1  |
!                 |  0   0   0   0   0   1   1   3   6   3   1  |
!                 |  0   0   0   0   0   0   1   1   3   6   2  |
!                 |  0   0   0   0   0   0   0   1   1   2   5  |
!
! stored as sparse matrix.
! B is a unit matrix:
!                 |  1   0   0   0   0   0   0   0   0   0   0  |
!                 |  0   1   0   0   0   0   0   0   0   0   0  |
!                 |  0   0   1   0   0   0   0   0   0   0   0  |
!                 |  0   0   0   1   0   0   0   0   0   0   0  |
!                 |  0   0   0   0   1   0   0   0   0   0   0  |
!    B    =       |  0   0   0   0   0   1   0   0   0   0   0  |.
!                 |  0   0   0   0   0   0   1   0   0   0   0  |
!                 |  0   0   0   0   0   0   0   1   0   0   0  |
!                 |  0   0   0   0   0   0   0   0   1   0   0  |
!                 |  0   0   0   0   0   0   0   0   0   1   0  |
!                 |  0   0   0   0   0   0   0   0   0   0   1  |
!
!  In what follows the symbol ' represents a transpose operation.
!
!  The test performs the following operations:
!
!  Step 1. Calls  FEASTINIT  to define the default values for the input
!          FEAST parameters.
!
!  Step 2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using
!          SFEAST_SCSREV.
!
!  Step 3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!           are the expected eigenvalues  and E(i) are eigenvalues computed
!           by SFEAST_SCSREV().
!
!  Step 4. The code computes the maximum absolute value of elements
!          of the matrix Y = X' *X - I, where X is the matrix of eigenvectors
!          computed by SFEAST_SCSREV.
!          SGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!
!  Step 5. The  code solves  the generalized eigenvalue problem Ax=eBx using
!          SFEAST_SCSRGV.
!
!  Step 6. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!          are the expected eigenvalues  and E(i) are eigenvalues computed
!          by SFEAST_SCSRGV().
!
!  Step 7. The code computes the maximum absolute value of the elements of
!          the matrix  Y = X' * X - I, where X is the matrix of eigenvectors
!          computed by SFEAST_SCSRGV.
!          SGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!
!*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"

#define max(a, b) (a) < (b) ? (b): (a)

int main()
{
    char          UPLO = 'F'; /* Type of matrix: (F=full matrix, L/U - lower/upper triangular part of matrix) */
    /* Matrix A of size N in CSR format. We use size of matrix N and 3 arrays to store matrix in CSR format */
    const MKL_INT N = 11;
    MKL_INT       rows[12] = { 1, 5, 10, 16, 23, 30, 37, 44, 51, 57, 62, 66 };
    MKL_INT       cols[65] = {    1,   2,   3,   4,
                                  1,   2,   3,   4,   5,
                                  1,   2,   3,   4,   5,   6,
                                  1,   2,   3,   4,   5,   6,   7,
                                       2,   3,   4,   5,   6,   7,   8,
                                            3,   4,   5,   6,   7,   8,   9,
                                                 4,   5,   6,   7,   8,   9,  10,
                                                      5,   6,   7,   8,   9,  10,  11,
                                                           6,   7,   8,   9,  10,  11,
                                                                7,   8,   9,  10,  11,
                                                                     8,   9,  10,  11
                            };

    float         val[65] = {   5.0, 2.0, 1.0, 1.0,
                                2.0, 6.0, 3.0, 1.0, 1.0,
                                1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                     1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                          1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                               1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                                    1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                                         1.0, 1.0, 3.0, 6.0, 3.0, 1.0,
                                                              1.0, 1.0, 3.0, 6.0, 2.0,
                                                                   1.0, 1.0, 2.0, 5.0
                            };

    /* Matrix B of size N in CSR format. We use size of matrix N and 3 arrays to store matrix in CSR format */
    MKL_INT       rowsb[12] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
    MKL_INT       colsb[11] = { 1,
                                    2,
                                        3,
                                            4,
                                                5,
                                                    6,
                                                        7,
                                                            8,
                                                                9,
                                                                    10,
                                                                        11
                              };
    float         valb[11] = {  1.0,
                                    1.0,
                                        1.0,
                                            1.0,
                                                1.0,
                                                    1.0,
                                                        1.0,
                                                            1.0,
                                                                1.0,
                                                                    1.0,
                                                                        1.0
                            };

    /* Declaration of FEAST variables */
    MKL_INT      fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */
    float        Emin, Emax;    /* Lower/upper bound of search interval [Emin,Emax] */

    float        epsout;        /* Relative error on the trace */
    MKL_INT      loop;          /* Number of refinement loop */
    MKL_INT      L = 9;
    MKL_INT      M0;            /* Initial guess for subspace dimension to be used */
    MKL_INT      M;             /* Total number of eigenvalues found in the interval */

    float        E[11];         /* Eigenvalues */
    float        X[121];        /* Eigenvectors */
    float        res[11];       /* Residual */

    /* Declaration of local variables */
    MKL_INT      info;          /* Errors */
    float        Eig[11];       /* Eig - array for storing exact eigenvalues */
    float        R[11];         /* R = |E-Eig| */
    float        Y[121];        /* Y=(X')*X-I */

    char         SGEMMC = 'T';  /* Character for GEMM routine, transposed case */
    char         SGEMMN = 'N';  /* Character for GEMM routine, non-transposed case */
    float        one = 1.0;     /* alpha parameter for GEMM */
    float        zero = 0.0;    /* beta  parameter for GEMM */
    MKL_INT      ldx = 11;      /* Leading dimension for source arrays in GEMM */
    MKL_INT      ldy = 11;      /* Leading dimension for destination array in GEMM */

    MKL_INT      i, j;
    float        trace, smax, eigabs;

    /* Exact eigenvalues in range (3.0, 7.0) */
    Eig[0] = 3.1715728752538100;
    Eig[1] = 4.0000000000000000;
    Eig[2] = 4.0000000000000000;
    Eig[3] = 4.1292484841890931;
    Eig[4] = 4.4066499006731521;
    Eig[5] = 6.0000000000000000;
    for (i=6; i<N; i++)
    {
        Eig[i] = 0.0;
    }

    printf("\n FEAST SFEAST_SCSREV AND SFEAST_SCSRGV EXAMPLE\n");
    /* Initialize matrix X */
    for (i=0; i<N*N; i++)
    {
        X[i] = zero;
    }

    printf("Sparse matrix size %i\n", (int)N);

    /* Search interval [Emin,Emax] */
    Emin = 3.0;
    Emax = 7.0;
    printf("Search interval [ %.7e, %.7e  ]  \n", Emin, Emax);

    M0   = L;
    M    = L;
    loop = 0;
    info = 0;
    epsout = 0.0;

    /* Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters */
    feastinit(
        fpm /* OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
        );

    fpm[0] =  1; /* Extended Eigensolver routines print runtime status to the screen. */
    fpm[1] = 24; /* The number of contour points Ne = 24.*/
    fpm[6] =  4; /* Error trace single precision stopping criteria 10e-fmp[6] */

    /* Step 2. Solve the standard Ax = ex eigenvalue problem. */
    printf("Testing sfeast_scsrev routine:\n");
    sfeast_scsrev(
        &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
        &N,      /* IN: Size of the problem */
        val,     /* IN: CSR matrix A, values of non-zero elements */
        rows,    /* IN: CSR matrix A, index of the first non-zero element in row */
        cols,    /* IN: CSR matrix A, columns indices for each non-zero element */
        fpm,     /* IN/OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
        &epsout, /* OUT: Relative error of on the trace */
        &loop,   /* OUT: Contains the number of refinement loop executed */
        &Emin,   /* IN: Lower bound of search interval */
        &Emax,   /* IN: Upper bound of search interval */
        &M0,     /* IN: The initial guess for subspace dimension to be used. */
        E,       /* OUT: The first M entries of Eigenvalues */
        X,       /* IN/OUT: The first M entries of Eigenvectors */
        &M,      /* OUT: The total number of eigenvalues found in the interval */
        res,     /* OUT: The first M components contain the relative residual vector */
        &info    /* OUT: Error code */
        );
    printf("FEAST OUTPUT INFO %d \n",info);
    if ( info != 0 )
    {
        printf("Routine sfeast_scsrev returns code of ERROR: %i", (int)info);
        return 1;
    }

    /* Step 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
    * are the expected eigenvalues and E(i) are eigenvalues computed by SFEAST_SCSREV(). */
    printf("Number of eigenvalues found %d \n", M);
    printf("Computed      |    Expected    \n");
    printf("Eigenvalues   |    Eigenvalues \n");
    eigabs = 0.0;
    for (i=0; i<M; i++)
    {
        R[i] = fabs(E[i] - Eig[i]);
        eigabs = max(eigabs, R[i]);
        printf("%.7e  %.7e \n", E[i], Eig[i]);
    }
    printf("Max value of | computed eigenvalue(i) - expected eigenvalues(i) | %.7e \n", eigabs);

    /* Step 4. The code computes the maximum absolute value of elements
     * of the matrix Y = X' *X - I, where X is the matrix of eigenvectors
     *  computed by SFEAST_SCSREV.
     *
     * Call BLAS to compute Y = X' * X  */
    sgemm(
        &SGEMMC, /* IN: 'T', transposed case*/
        &SGEMMN, /* IN: 'N', non-transposed case*/
        &M,      /* IN: Number of rows in matrix Y */
        &M,      /* IN: Number of columns in matrix X */
        &N,      /* IN: Number of columns in matrix Y */
        &one,    /* IN: alpha = 1.0 */
        X,       /* IN: Source #1 for GEMM, will be transposed */
        &ldx,    /* IN: Leading dimension of Source 1 */
        X,       /* IN: Source #2 for GEMM */
        &ldx,    /* IN: Leading dimension of Source 2 */
        &zero,   /* IN: beta = 0.0 */
        Y,       /* OUT: Destination */
        &ldy     /* IN: Leading dimension of Destination */
        );

    /* Compute Y = Y - I */
    for (i=0; i<M; i++)
    {
        Y[i*M + i] -= 1.0;
    }

    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("#Search interval [Emin,Emax] %.7e %.7e\n", Emin, Emax);
    printf("#mode found/subspace %d %d \n", M, M0);
    printf("#iterations %d \n", loop);
    trace = 0.0;
    for (i=0; i<M; i++)
    {
        trace += E[i];
    }
    printf("TRACE %.7e \n", trace);
    printf("Relative error on the Trace %.7e \n", epsout );
    printf("Index/Eigenvalues/Residuals\n");
    for (i=0; i<M; i++)
    {
        printf("%d  %.7e %.7e \n", i, E[i], res[i]);
    }
    smax = 0.0;
    for (i=0; i<M; i++)
    {
        for (j=0; j<M; j++)
        {
            smax = max(smax, fabs(Y[i*M + j]));
        }
    }
    printf( "Max(X' * X - I) = %.7e \n", smax);

    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    /*!!!!!!!!!!!!!!! GENERALIZED EIGENVALUE PROBLEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

    /* Reset initial parameters for new problem */
    M0 = L;
    for (i=0; i<N; i++)
    {
        E[i] = 0.0;
    }
    for (i=0; i<N*N; i++)
    {
        X[i] = zero;
    }

    /* Step 5. Solve the generalized eigenvalue problem Ax=eBx by SFEAST_SCSRGV */
    printf("Testing sfeast_scsrgv  \n");
    sfeast_scsrgv(
        &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
        &N,      /* IN: Size of the problem */
        val,     /* IN: CSR matrix A, values of non-zero elements */
        rows,    /* IN: CSR matrix A, index of the first non-zero element in row */
        cols,    /* IN: CSR matrix A, columns indices for each non-zero element */
        valb,    /* IN: CSR matrix B, values of non-zero elements */
        rowsb,   /* IN: CSR matrix B, index of the first non-zero element in row */
        colsb,   /* IN: CSR matrix B, columns indices for each non-zero element */
        fpm,     /* IN: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
        &epsout, /* OUT: Relative error of on the trace */
        &loop,   /* OUT: Contains the number of refinement loop executed */
        &Emin,   /* IN: Lower bound of search interval */
        &Emax,   /* IN: Upper bound of search interval */
        &M0,     /* IN/OUT: The initial guess for subspace dimension to be used. */
        E,       /* OUT: The first M entries of Eigenvalues */
        X,       /* OUT: The first M entries of Eigenvectors */
        &M,      /* OUT: The total number of eigenvalues found in the interval */
        res,     /* OUT: The first M components contain the relative residual vector */
        &info    /* OUT: Error code */
        );

    printf("FEAST OUTPUT INFO %d \n" ,info);
    if ( info != 0 )
    {
        printf("Routine sfeast_scsrgv return error: %i", (int)info);
        return 1;
    }
    /* Step 6. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
     *  are the expected eigenvalues  and E(i) are eigenvalues computed
     *  by  SFEAST_SCSRGV */
    printf("Number of eigenvalues found %d \n", M);
    printf("Computed      |    Expected    \n");
    printf("Eigenvalues   |    Eigenvalues \n");
    eigabs = 0.0;
    for (i=0; i<M; i++)
    {
        R[i] = fabs(E[i] - Eig[i]);
        eigabs = max(eigabs, R[i]);
        printf("%.7e  %.7e \n", E[i], Eig[i]);
    }
    printf("Max value of | computed eigenvalue(i) - expected eigenvalues(i) | %.7e \n", eigabs);

    /* Step 7. The code computes the maximum absolute value of the elements of
     * the matrix  Y = X' * X - I, where X is the matrix of eigenvectors
     * computed by SFEAST_SCSRGV.
     *
     * Call BLAS to compute X' * X */
    sgemm(
        &SGEMMC, /* IN: 'T', transposed case*/
        &SGEMMN, /* IN: 'N', non-transposed case*/
        &M,      /* IN: Number of rows in matrix Y */
        &M,      /* IN: Number of columns in matrix X */
        &N,      /* IN: Number of columns in matrix Y */
        &one,    /* IN: alpha = 1.0 */
        X,       /* IN: Source #1 for GEMM, will be transposed */
        &ldx,    /* IN: Leading dimension of Source 1 */
        X,       /* IN: Source #2 for GEMM */
        &ldx,    /* IN: Leading dimension of Source 2 */
        &zero,   /* IN: beta = 0.0 */
        Y,       /* OUT: Destination */
        &ldy     /* IN: Leading dimension of Destination */
        );

    /* Compute Y = Y - I */
    for (i=0; i<M; i++)
    {
        Y[i*M + i] -= 1.0;
    }

    /* Check the orthogonality of X' * X */
    smax = 0.0;
    for (i=0; i<M; i++)
    {
        for (j=0; j<M; j++)
        {
            smax = max(smax, fabs(Y[i*M + j]));
        }
    }

    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("# Search interval [Emin,Emax] %.7e %.7e\n", Emin, Emax);
    printf("# mode found/subspace %d %d \n",M,M0);
    printf("# iterations %d \n",loop);
    printf("Max(X' * X - I) = %.7e \n", smax);
    trace = 0.0;
    for (i=0; i<M; i++)
    {
        trace += E[i];
    }
    printf("TRACE %.7e \n", trace);
    printf("Relative error on the Trace %.7e \n", epsout);
    printf("Index/Eigenvalues/Residuals\n");
    for (i=0; i<M; i++)
    {
        printf("%d  %.7e %.7e \n", i, E[i], res[i]);
    }

    return 0;
}
