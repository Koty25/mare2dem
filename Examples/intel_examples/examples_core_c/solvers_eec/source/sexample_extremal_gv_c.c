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
!   Content: Example for k Max/Min generalized eigenvalue problem based on
!            Intel(R) Math Kernel Library (Intel(R) MKL)
!            Ax = l*Bx, A=A*, B=B*>0
!            Extended Eigensolver (CSR sparse format, float precision)
!
!*******************************************************************************
!
! The following routines are used in the example:
!          MKL_SPARSE_S_GV
!
! Consider the 4x4 matrices A, B
!
!                 |  9   0   0  |
!                 |  0   9   0  |
!     A   =       |  0   0   27  |
!
!                 |  5  -4   0  |
!                 | -4   5   0  |
!     B   =       |  0   0   1  |
!
! stored as sparse matrix.
!
!
!  The test calls mkl_sparse_s_gv routine to find several largest singular
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
    /* Matrices A and B of size N in CSR format */
    MKL_INT     N = 3;              /* number of rows/columns in matrix A,B */
    MKL_INT ia[4] = {1,2,3,4};      /* ia array from CSR format */
    MKL_INT ja[3] = {1,2,3};        /* ja array from CSR format */
    float    a[3] = {9.0,9.0,27.0}; /* val array from CSR format */

    MKL_INT ib[4] = {1,3,5,6};               /* ib array from CSR format */
    MKL_INT jb[5] = {1,2,1,2,3};             /* jb array from CSR format */
    float    b[5] = {5.0,-4.0,-4.0,5.0,1.0}; /* valb array from CSR format */

    float   Eig[3] = {1.0, 9.0, 27.0}; /* Exact eigenvalues */

    /* mkl_sparse_s_gv input parameters */
    char         which = 'S'; /* Which eigenvalues to calculate. ('L' - largest (algebraic) eigenvalues, 'S' - smallest (algebraic) eigenvalues) */
    MKL_INT      pm[128];     /* This array is used to pass various parameters to Extended Eigensolver Extensions routines. */
    MKL_INT      k0  = 2;     /* Desired number of max/min eigenvalues */

    /* mkl_sparse_s_gv output parameters */
    MKL_INT      k;      /* Number of eigenvalues found (might be less than k). */
    float        E[3];   /* Eigenvalues */
    float        X[3*3]; /* Eigenvectors */
    float        res[3]; /* Residual */

    /* Local variables */
    MKL_INT      info;               /* Errors */
    MKL_INT      compute_vectors = 0;/* Flag to compute eigenvectors */
    MKL_INT      tol = 5;            /* Tolerance */
    MKL_INT      i;

    /* Sparse BLAS IE variables */
    sparse_matrix_t A = NULL, B = NULL; /* Handle containing sparse matrix in internal data structure */
    struct matrix_descr descrA, descrB; /* Structure specifying sparse matrix properties */

    /* Create handle for matrix A stored in CSR format */
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL; /* Full matrix is stored */
    mkl_sparse_s_create_csr ( &A, SPARSE_INDEX_BASE_ONE, N, N, ia, ia+1, ja, a );

    descrB.type = SPARSE_MATRIX_TYPE_GENERAL; /* Full matrix is stored */
    mkl_sparse_s_create_csr ( &B, SPARSE_INDEX_BASE_ONE, N, N, ib, ib+1, jb, b );

    /* Step 2. Call mkl_sparse_ee_init to define default input values */
    mkl_sparse_ee_init(pm);

    pm[1] = tol; /* Set tolerance */
    pm[6] = compute_vectors;

    /* Step 3. Solve the standard Ax = ex eigenvalue problem. */
    info = mkl_sparse_s_gv(&which, pm, A, descrA, B, descrB, k0, &k, E, X, res);

    printf("mkl_sparse_s_gv output info %d \n",info);
    if ( info != 0 )
    {
        printf("Routine mkl_sparse_s_gv returns code of ERROR: %i", (int)info);
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
    mkl_sparse_destroy(A);
    mkl_sparse_destroy(B);
    return 0;
}
