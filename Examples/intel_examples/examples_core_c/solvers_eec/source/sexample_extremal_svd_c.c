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
!   Content: Example for k Max/Min singular value problem based on
!            Intel(R) Math Kernel Library (Intel(R) MKL)
!            Extended Eigensolver (CSR sparse format, single precision)
!
!*******************************************************************************
!
! The following routines are used in the example:
!          MKL_SPARSE_S_SVD
!
! Consider the 4x5 matrix A
!
!                 |  1   0   0   0   2  |
!                 |  0   0   3   0   0  |
!     A   =       |  0   0   0   1   0  |
!                 |  0   2   0   0   0  |
!
! stored as sparse matrix.
!
!
!  The test calls mkl_sparse_s_svd routine to find several largest singular 
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
    MKL_INT N = 4;               /* number of rows in matrix A */
    MKL_INT M = 5;               /* number of columns in matrix A */

    MKL_INT ia[5] = {1,
            3,
                4,
                5,
                    6};            /* ia array from CSR format */
    MKL_INT ja[5] = {1,
            5,
                3,
                4,
                2};            /* ja array from CSR format */
    float   a[5] = { 1.0, 
                    2.0, 
                3.0,
                 1.0,
            2.0};             /* val array from CSR format */
            
    float   Eig[4] = {1.0, 2.0, 2.2360679, 3.0 }; /* Exact singular values */
    
    /* Input variables for svd() */
    char    whichE = 'S';           /* Input: Which singular values to calculate. ('L' - largest (algebraic) singular values, 'S' - smallest (algebraic) singular values) */
    char    whichV = 'L';           /* Input: Which singular vectors to calculate. ('R' - right singular vectors, 'L' - left singular vectors */    
    MKL_INT pm[128];                /* Input: Array used to pass various parameters to Extended Eigensolver Extensions routines. */    
    MKL_INT k0     = 3;             /* Input: Number of largest(smallest) singular values to calculate */
    MKL_INT k;                      /* Output: Total number of sv calculated */
    float   E[5], XL[4*4], XR[5*5]; /* Output: E - k entries of singular values, XL - k corresponding left-singular vectors, XR - k corresponding right-singular vectors */
    float   res[5];                 /* Output: First k components contain the relative residual vector */
    
    /* Input variables for SGEMM */
    char        XGEMMC = 'T';  /* Character for GEMM routine, transposed case */
    char        XGEMMN = 'N';  /* Character for GEMM routine, non-transposed case */
    float       one = 1.0;     /* alpha parameter for GEMM */
    float       zero = 0.0;    /* beta  parameter for GEMM */
    
    /* Sparse BLAS IE variables */
    sparse_matrix_t A = NULL; /* Handle containing sparse matrix in internal data structure */
    struct matrix_descr descr; /* Structure specifying sparse matrix properties */
    
    /* Local variables */
    MKL_INT     i, j, info;            
    MKL_INT     compute_vectors = 1; /* Flag to compute singular vectors */
    MKL_INT     tol = 5;             /* Tolerance */    
    float       smax;
    float       Y[3*3];             /* Y=(XL')*XL-I */    

    /* Create handle for matrix A stored in CSR format */
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    mkl_sparse_s_create_csr ( &A, SPARSE_INDEX_BASE_ONE, N, M, ia, ia+1, ja, a );    

    /* Call mkl_sparse_ee_init to define default input values */
    mkl_sparse_ee_init(pm);
    pm[1] = tol; 
    pm[6] = compute_vectors;
    pm[7] = 1; /* Use absolute stopping critertia. Iteration is stopped if norm(AATx-λx) < 10-pm[1] */
    
    /* Call mkl_sparse_s_svd to obtain k0 singular values and singular vectors */    
    info = mkl_sparse_s_svd(&whichE, &whichV, pm, A, descr, k0, &k, E, XL, XR, res);    
    
    printf("mkl_sparse_s_svd output info %d \n",info);
    if ( info != 0 )
    {
        printf("Routine mkl_sparse_s_svd returns code of ERROR: %i", (int)info);
        return 1;
    }

    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("Number of singular values found %d \n", k);
    printf("Calculated SV  | Expected |  Residuals \n");        
    for (i=0; i<k; i++){
        printf("   %11.4f |  %5.4f  | %e\n", E[i], Eig[i], res[i]);fflush(0);
    }
    
    /* Call BLAS to compute Y = XL' * XL to check orthogonality of singular vectors */
    sgemm(
        &XGEMMC,    /* IN: 'T', transposed case*/
        &XGEMMN,    /* IN: 'N', non-transposed case*/
        &k,         /* IN: Number of rows in matrix Y */
        &k,         /* IN: Number of columns in matrix XL */
        &N,         /* IN: Number of columns in matrix Y */
        &one,       /* IN: alpha = 1.0 */
        XL,         /* IN: Source #1 for GEMM, will be transposed */
        &N,         /* IN: Leading dimension of Source 1 */
        XL,         /* IN: Source #2 for GEMM */
        &N,         /* IN: Leading dimension of Source 2 */
        &zero,      /* IN: beta = 0.0 */
        Y,          /* OUT: Destination */
        &k  /* IN: Leading dimension of Destination */
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
    printf( "Max(XL' * XL - I) = %.15e \n", smax);
    mkl_sparse_destroy(A);
    return 0;
}
