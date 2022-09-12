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
!   Content: Example for k Max/Min eigenvalue problem based on
!            Intel(R) Math Kernel Library (Intel(R) MKL) Extended
!            Eigensolver (CSR sparse format, float precision)
!
!*******************************************************************************
!
! The following routines are used in the example:
!          MKL_SPARSE_S_EV
!
! Consider the 4x4 matrix A
!
!                 |  6   2   0   0   |
!                 |  2   3   0   0   |
!     A   =       |  0   0   2  -1   |
!                 |  0   0  -1   2   |
!
! stored as sparse matrix.
!
!
!  The test calls mkl_sparse_s_ev routine to find several largest singular
!  values and corresponding right-singular vectors. Orthogonality of singular
!  vectors is tested using SGEMM routine
!
!*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"
#include "mkl_solvers_ee.h"
#define max(a, b) (a) < (b) ? (b): (a)

int main()
{
    /* Matrix A of size N in CSR format */
    MKL_INT     N = 4;                                   /* number of rows in matrix A */
    MKL_INT ia[5] = {1,3,5,7,9};                         /* ia array from CSR format */
    MKL_INT ja[8] = {1,2,1,2,3,4,3,4};                   /* ja array from CSR format */
    float    a[8] = {6.0,2.0,2.0,3.0,2.0,-1.0,-1.0,2.0}; /* val array from CSR format */

    float   Eig[4] = {1.0, 2.0, 3.0, 7.0}; /* Exact eigenvalues */

    /* mkl_sparse_s_ev input parameters */
    char         which = 'S'; /* Which eigenvalues to calculate. ('L' - largest (algebraic) eigenvalues, 'S' - smallest (algebraic) eigenvalues) */
    MKL_INT      pm[128];     /* This array is used to pass various parameters to Extended Eigensolver Extensions routines. */
    MKL_INT      k0  = 3;     /* Desired number of max/min eigenvalues */

    /* mkl_sparse_s_ev output parameters */
    MKL_INT      k;      /* Number of eigenvalues found (might be less than k). */
    float        E[4];   /* Eigenvalues */
    float        X[4*4]; /* Eigenvectors */
    float        res[4]; /* Residual */

    /* Local variables */
    MKL_INT      info;               /* Errors */
    MKL_INT      compute_vectors = 1;/* Flag to compute eigenvectors */
    MKL_INT      tol = 5;            /* Tolerance */
    float        Y[4*4];             /* Y=(X')*X-I */
    MKL_INT      i, j;
    float        smax;

    /* Input variables for DGEMM */
    char         DGEMMC = 'T';       /* Character for GEMM routine, transposed case */
    char         DGEMMN = 'N';       /* Character for GEMM routine, non-transposed case */
    float        one  = 1.0;         /* alpha parameter for GEMM */
    float        zero = 0.0;         /* beta  parameter for GEMM */
    MKL_INT      ldx  = N;           /* Leading dimension for source arrays in GEMM */
    MKL_INT      ldy;                /* Leading dimension for destination array in GEMM */

    /* Sparse BLAS IE variables */
    sparse_matrix_t A = NULL; /* Handle containing sparse matrix in internal data structure */
    struct matrix_descr descr; /* Structure specifying sparse matrix properties */

    /* Create handle for matrix A stored in CSR format */
    descr.type = SPARSE_MATRIX_TYPE_GENERAL; /* Full matrix is stored */
    mkl_sparse_s_create_csr ( &A, SPARSE_INDEX_BASE_ONE, N, N, ia, ia+1, ja, a );

    /* Step 2. Call mkl_sparse_ee_init to define default input values */
    mkl_sparse_ee_init(pm);

    pm[1] = tol; /* Set tolerance */
    pm[6] = compute_vectors;

    /* Step 3. Solve the standard Ax = ex eigenvalue problem. */
    info = mkl_sparse_s_ev(&which, pm, A, descr, k0, &k, E, X, res);

    printf("mkl_sparse_s_ev output info %d \n",info);
    if ( info != 0 )
    {
        printf("Routine mkl_sparse_s_ev returns code of ERROR: %i", (int)info);
        return 1;
    }

    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("#mode found/subspace %d %d \n", k, k0);
    printf("Index/Exact Eigenvalues/Estimated Eigenvalues/Residuals\n");
    for (i=0; i<k; i++)
    {
        printf("   %d  %.15e %.15e %.15e \n",i, Eig[i], E[i], res[i]);
    }

    if ( compute_vectors )
    {
        /* The code computes the maximum absolute value of elements
         * of the matrix Y = X' *X - I, where X is the matrix of eigenvectors
         *  computed by mkl_sparse_s_ev.
         *
         * Call BLAS to compute Y = X' * X  */

        ldy = k;
        sgemm(
            &DGEMMC, /* IN: 'T', transposed case*/
            &DGEMMN, /* IN: 'N', non-transposed case*/
            &k,      /* IN: Number of rows in matrix Y */
            &k,      /* IN: Number of columns in matrix Y */
            &N,      /* IN: Number of rows in matrix X */
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
        for (i=0; i<k; i++)
        {
            Y[i*k+i] -= 1.0;
        }

        smax = 0.0;
        for (i=0; i<k; i++)
        {
            for (j=0; j<k; j++)
            {
                smax = max(smax, fabs(Y[i*k+j]));
            }
        }
        printf( "Max(X' * X - I) = %.15e \n", smax);
    }
    mkl_sparse_destroy(A);
    return 0;
}
