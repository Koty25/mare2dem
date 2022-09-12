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
!            Eigensolvers (dense format, single complex precision)
!
!*******************************************************************************
!
! The following routines are used in the example:
!          CGEMM  CFEAST_HEEV  CFEAST_HEGV  FEASTINIT.
!
! Consider the matrix A:
!
!                 |  10   1+2i  0    0    0    0    0    0    0    0    |
!                 |  1-2i  9   2+3i  0    0    0    0    0    0    0    |
!                 |  0    2-3i  8   3+4i  0    0    0    0    0    0    |
!                 |  0    0    3-4i  7   4+5i  0    0    0    0    0    |
!                 |  0    0    0    4-5i  6   5+6i  0    0    0    0    |
!    A    =       |  0    0    0    0    5-6i  5   6+7i  0    0    0    |,
!                 |  0    0    0    0    0    6-7i  4   7+8i  0    0    |
!                 |  0    0    0    0    0    0    7-8i  3   8+9i  0    |
!                 |  0    0    0    0    0    0    0    8-9i  2   9+10i |
!                 |  0    0    0    0    0    0    0    0    9-10i 1    |
!
! stored as dense matrix.
! B is a unit matrix:
!
!                 |  1   0   0   0   0   0   0   0   0   0  |
!                 |  0   1   0   0   0   0   0   0   0   0  |
!                 |  0   0   1   0   0   0   0   0   0   0  |
!                 |  0   0   0   1   0   0   0   0   0   0  |
!                 |  0   0   0   0   1   0   0   0   0   0  |
!    B    =       |  0   0   0   0   0   1   0   0   0   0  |.
!                 |  0   0   0   0   0   0   1   0   0   0  |
!                 |  0   0   0   0   0   0   0   1   0   0  |
!                 |  0   0   0   0   0   0   0   0   1   0  |
!                 |  0   0   0   0   0   0   0   0   0   1  |
!
!  In what follows the symbol ' represents a conjugate transpose operation.
!
!  The test performs the following operations :
!
!  Step 1. Calls  FEASTINIT  to define the default values for the input
!          FEAST parameters.
!
!  Step 2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using
!          CFEAST_HEEV.
!
!  Step 3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!          are the expected eigenvalues  and E(i) are eigenvalues computed
!          by CFEAST_HEEV().
!
!  Step 4. The code computes the maximum absolute value of elements
!          of the matrix Y = X' *X - I, where X is the matrix of eigenvectors
!          computed by CFEAST_HEEV.
!          CGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!
!  Step 5. The  code solves  the generalized eigenvalue problem Ax=eBx using
!          CFEAST_HEGV.
!
!  Step 6. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!          are the expected eigenvalues  and E(i) are eigenvalues computed
!          by CFEAST_HEGV().
!
!  Step 7. The code computes the maximum absolute value of the elements of
!          the matrix  Y = X' * X - I, where X is the matrix of eigenvectors
!          computed by CFEAST_HEGV.
!          CGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!
!*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"

#define max(a, b) (a) < (b) ? (b): (a)

int main()
{
    char          UPLO = 'F'; /* Type of matrix: (F=full matrix, L/U - lower/upper triangular part of the matrix) */
    /* Matrix A in dense format, size N by N, column-style Fortran storage for Feast*/
    const MKL_INT      N = 10;
    MKL_Complex8       A[100] = {
          {10.0, 0.0}, {1.0,-2.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0,  0.0}, {0.0,  0.0}, 
          { 1.0, 2.0}, {9.0, 0.0}, {2.0,-3.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0,  0.0}, {0.0,  0.0},
          { 0.0, 0.0}, {2.0, 3.0}, {8.0, 0.0}, {3.0,-4.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0,  0.0}, {0.0,  0.0},
          { 0.0, 0.0}, {0.0, 0.0}, {3.0, 4.0}, {7.0, 0.0}, {4.0,-5.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0,  0.0}, {0.0,  0.0},
          { 0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {4.0, 5.0}, {6.0, 0.0}, {5.0,-6.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0,  0.0}, {0.0,  0.0},
          { 0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {5.0, 6.0}, {5.0, 0.0}, {6.0,-7.0}, {0.0, 0.0}, {0.0,  0.0}, {0.0,  0.0}, 
          { 0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {6.0, 7.0}, {4.0, 0.0}, {7.0,-8.0}, {0.0,  0.0}, {0.0,  0.0}, 
          { 0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {7.0, 8.0}, {3.0, 0.0}, {8.0, -9.0}, {0.0,  0.0}, 
          { 0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {8.0, 9.0}, {2.0,  0.0}, {9.0,-10.0},
          { 0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {9.0, 10.0}, {1.0,  0.0}   
    }; 

    /* Matrix B in dense format, size N by N*/
    MKL_Complex8  B[100];

    /* Declaration of FEAST variables */
    MKL_INT       fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */
    float         Emin, Emax;    /* Lower/upper bound of search interval [Emin,Emax] */

    float         epsout;        /* Relative error on the trace */
    MKL_INT       loop;          /* Number of refinement loop */
    const MKL_INT L = 8;
    MKL_INT       M0;            /* Initial guess for subspace dimension to be used */
    MKL_INT       M;             /* Total number of eigenvalues found in the interval */

    float         E[10];         /* Eigenvalues */
    MKL_Complex8  X[100];        /* Eigenvectors */
    float         res[10];       /* Residual */

    /* Declaration of local variables */
    MKL_INT       info;          /* Errors */
    float         Eig[10];       /* Eig - array for storing exact eigenvalues */
    float         R[10];         /* R = |E-Eig| */
    MKL_Complex8  Y[100];        /* Y=(X')*X-I */
    
    char          CGEMMC='C';    /* Character for GEMM routine, conjugate transposed case */
    char          CGEMMN='N';    /* Character for GEMM routine, non-transposed case */
    MKL_Complex8  one  = {1.0, 0.0}; /* alpha parameter for GEMM */
    MKL_Complex8  zero = {0.0, 0.0}; /* beta  parameter for GEMM */
    MKL_INT       ldx = 10;      /* Leading dimension for source arrays in GEMM */
    MKL_INT       ldy = 10;      /* Leading dimension for destination array in GEMM */

    MKL_INT       i, j;
    float         trace, smax, eigabs;

    /* Exact eigenvalues in range (2.0, 12.0) */
    Eig[0] = 2.231051;
    Eig[1] = 6.058517;
    Eig[2] = 9.109751;
    Eig[3] = 11.703148;

    for (i=4; i<N; i++)
    {
        Eig[i] = 0.0;
    }

    printf("\n     FEAST CFEAST_HEEV AND CFEAST_HEGV EXAMPLE\n"); 
    /* Initialize matrices B and X */
    for (i=0; i<N*N; i++)
    {
        B[i] = zero;
        X[i] = zero;
    }

    for (i=0; i<N; i++)
    {
        B[i+i*N] = one;
    }

    printf(" Dense matrix size %i\n", (int) N);

    /* Search interval [Emin,Emax] */
    Emin = 2.0;
    Emax = 12.0;
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

    /* Step 2. Solve the standard Ax = ex eigenvalue problem. */
    printf(" Testing cfeast_heev routine:\n");
    cfeast_heev(
        &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
        &N,      /* IN: Size of the problem */
        A,       /* IN: dense matrix A */
        &N,      /* IN: The first dimension of the matrix A */
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
        printf("Routine cfeast_heev returns code of ERROR: %i", (int)info);
        return 1;
    }

    /* Step 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
    * are the expected eigenvalues and E(i) are eigenvalues computed by CFEAST_HEEV(). */
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
    printf("Max value of | computed eigenvalue(i) - expected eigenvalue(i) | %.7e \n", eigabs);

    /* Step 4. The code computes the maximum absolute value of elements
     * of the matrix Y = X' *X - I, where X is the matrix of eigenvectors
     * computed by CFEAST_HEEV.
     *
     * Call BLAS to compute Y = X' * X  */

    cgemm(
        &CGEMMC, /* IN: 'C', conjugate transposed case*/
        &CGEMMN, /* IN: 'N', non-transposed case*/
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
        Y[i*M + i].real -= 1.0;
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
            smax=max(smax, sqrt(Y[i*M + j].imag*Y[i*M + j].imag+Y[i*M + j].real*Y[i*M + j].real) );
        }
    }
    printf("Max(X' * X - I) = %.7e \n", smax);

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

    /* Step 5. Solve the generalized eigenvalue problem Ax=eBx by CFEAST_HEGV */
    printf("Testing cfeast_hegv  \n");
    cfeast_hegv(
        &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
        &N,      /* IN: Size of the problem */
        A,       /* IN: dense matrix A */
        &N,      /* IN: The first dimension of the matrix A */
        B,       /* IN: dense matrix B */
        &N,      /* IN: The first dimension of the matrix B */
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
        printf("Routine cfeast_hegv return error: %i", (int)info);
        return 1;
    }
    /* Step 6. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
     *  are the expected eigenvalues  and E(i) are eigenvalues computed
     *  by  CFEAST_HEGV */
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
    printf("Max value of | computed eigenvalue(i) - expected eigenvalue(i) | %.7e \n", eigabs);

    /* Step 7. The code computes the maximum absolute value of the elements of
     * the matrix  Y = X' * X - I, where X is the matrix of eigenvectors
     * computed by CFEAST_HEGV.
     *
     * Call BLAS to compute X' * X */

    cgemm(
        &CGEMMC, /* IN: 'C', conjugate transposed case*/
        &CGEMMN, /* IN: 'N', non-transposed case*/
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
        Y[i*M + i].real -= 1.0;
    }

    /* Check the orthogonality of X' * X */
    smax = 0.0;
    for (i=0; i<M; i++)
    {
        for (j=0; j<M; j++)
        {
            smax = max(smax, sqrt(Y[i*M + j].imag * Y[i*M + j].imag + Y[i*M + j].real * Y[i*M + j].real));
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
